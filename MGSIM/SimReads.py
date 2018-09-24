from __future__ import print_function
# import
## batteries
import os
import sys
import re
import time
import subprocess
from glob import glob
from shutil import rmtree
from functools import partial
from multiprocessing import Pool
from distutils.spawn import find_executable
## 3rd party
import pandas as pd
from Bio import SeqIO
## application
from MGSIM import Utils


def tidy_taxon_names(x):
    x = re.sub(r'[()\/:;, ]+', '_', x)
    return x

def load_genome_table(in_file):
    """Loading genome table
    """
    df = pd.read_csv(in_file, sep='\t')
    
    ## check headers
    diff = set(['Taxon','Fasta']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # getting genome sizes
    df['Genome_size'] = [_genome_size(x) for i,x in df.iterrows()]

    # tidy taxon names
    df['Taxon'] = df['Taxon'].apply(tidy_taxon_names)
    
    return df

def _genome_size(x):
    """Get the total bp for the genome sequence of all genomes
        
    Parameters
    ----------
    """
    bp = 0
    for record in SeqIO.parse(x['Fasta'], 'fasta'):
        bp += len(record.seq)
    return bp

def load_abund_table(in_file):
    """Loading abundance table
    """
    df = pd.read_csv(in_file, sep='\t')
    
    ## check headers
    diff = set(['Community','Taxon','Perc_rel_abund']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # tidy taxon names
    df['Taxon'] = df['Taxon'].apply(tidy_taxon_names)

    return df

def sample_taxon_list(genome_table, abund_table):
    """Creating [sample,taxon] list of lists
    """
    # joining tables
    df = abund_table.merge(genome_table, on=['Taxon'])
    # convert to a list of lists
    sample_taxon = []
    cols = ['Community', 'Taxon', 'Genome_size', 'Fasta', 'Perc_rel_abund']
    for i,x in df[cols].iterrows():
        sample_taxon.append(x.tolist())
    return sample_taxon

def sim_illumina(sample_taxon, output_dir, seq_depth, art_params,
                 temp_dir, nproc=1, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    x : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    """
    # check that simulator exists
    exe = 'art_illumina'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))
    
    # calculating fold per genome
    for x in sample_taxon:
        genome_size = float(x[2])
        perc_rel_abund = float(x[4])
        try:
            _ = art_params['--paired']
            paired = 2
        except KeyError:
            try:
                art_params['--mflen']
                paired = 2
            except KeyError:
                paired = 1
        read_length  = art_params['--len'] * paired
        fold = perc_rel_abund / 100.0 * seq_depth * read_length / genome_size
        x.append(fold)

    # directories
    ## temp dir
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    ## output dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    # simulate per sample
    sys.stderr.write('Simulating reads...\n')
    func = partial(sim_art,
                   art_params=art_params,
                   temp_dir=temp_dir,
                   debug=debug)
    if debug is True:
        res = map(func, sample_taxon)
    else:
        p = Pool(nproc)
        res = p.map(func, sample_taxon)
    res = list(res)

    # combining all reads by sample
    sys.stderr.write('Combining simulated reads by sample...\n')    
    comms = list(set([x[0] for x in sample_taxon]))
    func = partial(combine_reads_by_sample,
                   temp_dir=temp_dir,
                   file_prefix='illumina',
                   output_dir=output_dir,
                   nproc=nproc,
                   debug=debug)
    if debug is True:
        res = map(func, comms)
    else:
        p = Pool(nproc)
        res = p.map(func, comms)
    res = list(res)

    # removing temp dir
    rmtree(temp_dir)
    
    # status
    for sample_list in res:
        for file_name in sample_list:
            if file_name is not None:
                print('File written: {}'.format(file_name))

def sim_art(x, art_params, temp_dir, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    x : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    """
    community = x[0]
    taxon = x[1]
    fasta = x[3]
    perc_rel_abund = x[4]
    fold = x[5]
    
    # output
    ## temporary directories
    out_dir = os.path.join(temp_dir, str(community), taxon)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    ## temporary files
    output_prefix = os.path.join(out_dir, 'illumina')
    
    # art command
    art_params = ' '.join(['{} {}'.format(k,v) for k,v in art_params.items()])
    cmd = 'art_illumina {art_params} -f {fold} -i {input} -o {output_prefix}'
    cmd = cmd.format(art_params=art_params,
                     fold=fold,
                     input=fasta,
                     output_prefix=output_prefix)

    ## system call
    if debug is True:
        sys.stderr.write('CMD: ' + cmd + '\n')
    try:
        res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
    res = res.stdout.decode()
    if debug is True:
        sys.stderr.write(res + '\n')
    
    return res
   
def combine_reads_by_sample(sample, temp_dir, file_prefix, output_dir, nproc=1, debug=False):
    # find files
    ## read1
    p = os.path.join(temp_dir, str(sample), '*', file_prefix + '1.fq')
    R1_files = glob(p)
    if len(R1_files) == 0:
        p = os.path.join(temp_dir, str(sample), '*', file_prefix + '.fq')
        R1_files = glob(p)
        R2_files = None
    else:
        ## read2
        p = os.path.join(temp_dir, str(sample), '*', file_prefix + '2.fq')
        R2_files = glob(p)

    # writing files
    output_dir = os.path.join(output_dir, str(sample))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    ## read1
    R1_files = _combine_reads(R1_files, output_dir, 'R1.fq')
    ## read2
    if R2_files is not None:
        R2_files = _combine_reads(R2_files, output_dir, 'R2.fq')

    return [R1_files, R2_files]
        
def _combine_reads(read_files, output_dir, output_file):
    """
    Parameters
    ----------
    """
    output_file = os.path.join(output_dir, output_file)
    with open(output_file, 'w') as outFH:
        for in_file in read_files:
            taxon = os.path.split(os.path.split(in_file)[0])[1]
            for i,record in enumerate(SeqIO.parse(in_file, 'fastq')):
                # renaming fastq read
                name =  '{}__SEQ{}'.format(taxon, i)
                record.id = name
                record.description = name
                SeqIO.write(record, outFH, 'fastq')

    return output_file
