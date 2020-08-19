from __future__ import print_function
# import
## batteries
import os
import sys
import re
import time
import logging
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

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions
def tidy_taxon_names(x):
    """Remove special characters from taxon names
    """
    x = re.sub(r'[()\/:;, ]+', '_', x)
    return x

def load_genome_table(in_file, nproc=1):
    """Loading genome table
    Parameters
    ----------
    in_file : str
        input file path
    """
    nproc = int(nproc)
    df = pd.read_csv(in_file, sep='\t')
    
    ## check headers
    diff = set(['Taxon','Fasta']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # getting genome sizes
    if nproc > 1:
        p = Pool(nproc)
        df['Genome_size'] = p.map(_genome_size, [x for i,x in df.iterrows()])
    else:
        df['Genome_size'] = [_genome_size(x) for i,x in df.iterrows()]

    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)
    
    return df

def _genome_size(x):
    """Get the total bp for the genome sequence of all genomes
        
    Parameters
    ----------
    x : pd.Series
       Series that includes 'Fasta' in the index
    """
    bp = 0
    for record in SeqIO.parse(x['Fasta'], 'fasta'):
        bp += len(record.seq)
    return bp

def load_abund_table(in_file):
    """Loading abundance table
    Parameters
    ----------
    in_file : str
        input file path
    """
    df = pd.read_csv(in_file, sep='\t')
    
    ## check headers
    diff = set(['Community','Taxon','Perc_rel_abund']) - set(df.columns.values)
    if len(diff) > 0:
        diff = ','.join(diff)
        raise ValueError('Cannot find table columns: {}'.format(diff))

    # tidy taxon names
    df['Taxon'] = df['Taxon'].astype(str).apply(tidy_taxon_names)

    return df

def sample_taxon_list(genome_table, abund_table):
    """Creating [sample,taxon] list of lists
    Parameters
    ----------
    genome_table : pd.DataFrame
    abund_table : pd.DataFrame
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
                 temp_dir, nproc=1, rndSeed=None, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    sample_taxon : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    output_dir : str
        Output director for all read files
    seq_depth : int
        Sequencing depth per sample
    art_params : dict
        Parameters provided to art_illumina
    temp_dir : str
        Temporary file directory
    nproc : int
        Number of parallel processes
    debug : bool
        Debug mode
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
    logging.info('Simulating reads...')
    func = partial(sim_art,
                   art_params=art_params,
                   temp_dir=temp_dir,
                   rndSeed=rndSeed,
                   debug=debug)
    if debug is True:
        fq_files = map(func, sample_taxon)
    else:
        p = Pool(nproc)
        fq_files = p.map(func, sample_taxon)
    fq_files = list(fq_files)

    # combining all reads by sample
    logging.info('Combining simulated reads by sample...')    
    comms = list(set([x[0] for x in sample_taxon]))
    func = partial(combine_reads_by_sample,
                   fq_files=fq_files,
                   temp_dir=temp_dir,
                   file_prefix='illumina',
                   output_dir=output_dir,
                   debug=debug)
    if debug is True:
        res = map(func, comms)
    else:
        p = Pool(nproc)
        res = p.map(func, comms)
    res = list(res)

    # removing temp dir
    logging.info('Removing temp directory...')
    rmtree(temp_dir)
    
    # status
    for sample_list in res:
        for file_name in sample_list:
            if file_name is not None:
                logging.info('File written: {}'.format(file_name))

def sim_art(x, art_params, temp_dir, rndSeed=None, debug=False):
    """Simulate illumina reads
    Parameters
    ----------
    x : list
        [Community,Taxon,Genome_size,Fasta,Perc_rel_abund,Fold]
    """
    community = str(x[0])
    taxon = str(x[1])
    fasta = str(x[3])
    perc_rel_abund = float(x[4])
    fold = float(x[5])
    
    # output
    ## temporary directories
    out_dir = os.path.join(temp_dir, str(community), taxon)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    ## temporary files
    output_prefix = os.path.join(out_dir, 'illumina')
    
    # art command
    art_params = ' '.join(['{} {}'.format(k,v) for k,v in art_params.items()])
    cmd = 'art_illumina {art_params} --noALN -f {fold} -i {input} -o {output_prefix}'
    cmd = cmd.format(art_params=art_params,
                     fold=fold,
                     rndSeed=rndSeed,
                     input=fasta,
                     output_prefix=output_prefix)
    if rndSeed is not None:
        cmd += ' --rndSeed {}'.format(rndSeed)

    ## system call
    if debug is True:
        sys.stderr.write('CMD: ' + cmd + '\n')
    try:
        res = subprocess.run(cmd, check=True, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        raise e
    if debug is True:
        sys.stderr.write(res.stderr.decode() + '\n')
        sys.stderr.write(res.stdout.decode() + '\n')

    # check that files have been created
    R0_file = output_prefix + '.fq'
    R1_file = output_prefix + '1.fq'
    R2_file = output_prefix + '2.fq'
    if os.path.isfile(R1_file) and os.path.isfile(R2_file):
        return [community, R1_file, R2_file]
    elif os.path.isfile(R0_file):
        return [community, R0_file]
    else:
        msg = 'Cannot find art_illumina output files!'
        raise ValueError(msg)
   
def combine_reads_by_sample(sample, fq_files, temp_dir, file_prefix, output_dir, debug=False):
    """ Concat all sample-taxon read files into per-sample read files
    Parameters
    ----------
    sample : str
        Sample ID
    temp_dir : str
        Temporary directory path
    file_prefix : str
        Output file prefix
    output_dir : str
        Output directory path
    debug : bool
        Debug mode
    """
    # sorting fastq files by read pair
    sample = str(sample)
    R1_files = [x[1] for x in fq_files if x[0] == sample]
    try:
        R2_files = [x[2] for x in fq_files if x[0] == sample]
    except IndexError:
        R2_files = None

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
    """Combine fastq read files into 1 read file.
    Parameters
    ----------
    read_files : list
        All read files to combine
    output_dir : str
        Output directory path
    output_file : str
        Output file path
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
            # delete temporary file
            os.remove(in_file)

    return output_file
