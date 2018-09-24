from __future__ import print_function
# import
## batteries
import sys, os
import time
from functools import partial
from multiprocessing import Pool
from distutils.spawn import find_executable
## 3rd party
import pandas as pd
from Bio import SeqIO
## application
from MGSIM import Utils

    
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

    return df

def _genome_size(x):
    """Get the total bp for the genome sequence of all genomes
        
    Parameters
    ----------
    """
    bp = 0
    for record in SeqIO.parse(x['Fasta'], "fasta"):
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
    
    return df

def sample_taxon_list(genome_table, abund_table):
    """Creating [sample,taxon] list of lists
    """
    # joining tables
    df = abund_table.merge(genome_table, on=['Taxon'])
    # convert to a list of lists
    sample_taxon = []
    cols = ['Taxon', 'Genome_size', 'Perc_rel_abund']
    for i,x in df[cols].iterrows():
        sample_taxon.append(x.tolist())
    return sample_taxon

def sim_illumina(sample_taxon, seq_depth, art_params,
                 temp_dir, nproc=1, debug=False):
    """Simulate illumina reads
    """
    # check that simulator exists
    exe = 'art_illumina'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))

    # calculating fold per genome
    for x in sample_taxon:
        perc_rel_abund = x[2]
        genome_size = x[1]
        fold = perc_rel_abund / 100 * seq_depth / genome_size
        x.append(fold)
        
    # simulate per sample
    func = partial(sim_art,
                   art_params=art_params,
                   temp_dir=temp_dir)
    if debug is True:
        read_fasta = map(func, sample_taxon)
    else:
        p = Pool(nproc)
        read_fasta = p.map(func, sample_taxon)
    read_fasta = list(read_fasta)

def sim_art(x, art_params, temp_dir):
    """Simulate illumina reads
    Parameters
    ----------
    x : list
        [Taxon,Genome_size,Perc_rel_abund,Fold]
    """
    cmd = 'art_illumina {art_params} -f {fold} -i {input} -o {output_prefix}'
    print(x)
    print(cmd)
    
