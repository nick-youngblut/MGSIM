from __future__ import print_function
# import
## batteries
import os
import sys
import re
import time
import uuid
import logging
import subprocess
from glob import glob
from pprint import pprint
from shutil import rmtree
from random import shuffle
from functools import partial
from multiprocessing import Pool
from distutils.spawn import find_executable
## 3rd party
import numpy as np
import pandas as pd
import pyfastx
from Bio import SeqIO
#from pyfaidx import Fasta
from scipy.stats import truncnorm
## application
from MGSIM import Utils

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# functions
def bin_barcodes(barcodes, binsize=1000):
    """Binning barcodes into chunks
    Parameters
    ----------
    barcodes : iterable
        Iterable of barcodes
    binsize : int
        Size of bin for grouping barcodes
    
    Returns
    -------
    yields list of barcode (1 bin)
    """
    binsize = int(float(binsize))
    bins = np.digitize(np.arange(0,barcodes.shape[0]),
                       np.arange(0,barcodes.shape[0],binsize))
    return [barcodes[bins == x] for x in np.unique(bins)]

def sim_barcode_reads(barcodes, genome_table, abund_table, tmp_dir, args,
                      rndSeed=None, debug=False):
    """Simulate reads for each fragment associated with each barcode
    Parameters
    ----------
    barcodes : iterable
        Iterable of barcodes [str]
    genome_table : pd.dataframe
        Table that maps genome fasta to each taxon
    abund_table : pd.dataframe
        Table listing relative abundances of each taxon
    args: dict
        Command line args 
    rndSeed: int
        --rndSeed for ART
    debug: bool
        Debug mode
    """
    # taxon abundance table of community of interest
    genome_table = genome_table.merge(abund_table,on=['Taxon'])
    ## filtering out zeros
    genome_table = genome_table[genome_table['Perc_rel_abund'] > 0].reset_index()
    
    # simulate fragment reads on a per-barcode level
    frag_tsv_files = []
    read_fq_files = []
    for barcode in barcodes:
        if debug:
            logging.info('Processing barcode: {}'.format(barcode))
        # number of fragments
        n_frags = np.random.normal(loc=float(args['--frag-bc-mean']),
                                   scale=float(args['--frag-bc-sd']))
        n_frags = int(round(n_frags, 0))
        
        # selecting reference genomes (with replacement) for generating random fragments 
        refs = select_refs(n_frags, genome_table)
        
        # sim fragment sizes
        refs = sim_frags(refs, n_frags,
                         frag_size_loc=float(args['--frag-size-mean']),
                         frag_size_sd=float(args['--frag-size-sd']),
                         frag_size_min=int(args['--frag-size-min']),
                         frag_size_max=int(args['--frag-size-max']))

        # parsing fragments
        fasta_file = os.path.join(tmp_dir, barcode + '.fq')
        tsv_file = os.path.join(tmp_dir, barcode + '.tsv')
        with open(fasta_file, 'w') as outFr, open(tsv_file, 'w') as outFt:
            for taxon in refs['Taxon'].unique():
                parse_frags(refs.loc[refs['Taxon'] == taxon],
                            barcode=barcode, outFr=outFr, outFt=outFt)
        
        # simulating reads
        art_params = {'--paired' : args['--art-paired'],
                      '--len' : args['--art-len'],
                      '--mflen' : args['--art-mflen'],
                      '--sdev' : args['--art-sdev'],
                      '--seqSys' : args['--art-seqSys']}
        if rndSeed is not None:
            art_params['--rndSeed'] = rndSeed
        fasta_files = sim_illumina(fasta_file, str(barcode),
                                   seq_depth=float(args['--seq-depth']),
                                   total_barcodes=int(float(args['--barcode-total'])),
                                   art_params=art_params,
                                   tmp_dir=tmp_dir,
                                   debug=False)        
        read_fq_files.append(fasta_files)
        frag_tsv_files.append(tsv_file)

    return [read_fq_files, frag_tsv_files]

def format_params(d):
    dd = {}
    for k,v in d.items():
        if v == True:
            dd[k] = ''
        elif v == False:
            continue
        else:
            dd[k] = v
    return dd

def sim_illumina(frag_fasta, barcode, seq_depth, total_barcodes,
                 art_params, tmp_dir, debug=False):
    """Simulating illumina reads for each fragment from each barcode
    Parameters
    ----------
    frag_fasta : str
        Fasta file of simulated gDNA fragments for a particular barcode
    barcode : str
        Barcode ID
    seq_depth : int
        Total sequencing depth
    total_barcodes : int
        Total number of barcodes
    art_params : dict
        art_illumina parameters
    tmp_dir : str
        temporary directory path
    debug : bool
        debug mode 
    """    
    output_prefix = os.path.splitext(frag_fasta)[0] + '_'

    # formatting params
    art_params = format_params(art_params)
    
    # calculate fold coverage to simulate
    fold = calc_fold(frag_fasta, seq_depth, total_barcodes, art_params, tmp_dir=tmp_dir)
    
    # art_illumina command
    art_params = ' '.join(['{} {}'.format(k,v) for k,v in art_params.items()])
    cmd = 'art_illumina {art_params} --noALN --id {barcode}'
    cmd += ' -f {fold} -i {input} -o {output_prefix}'
    cmd = cmd.format(art_params=art_params,
                     barcode=barcode + '-',
                     fold=fold,
                     input=frag_fasta,
                     output_prefix=output_prefix)

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
        return [R1_file, R2_file]
    elif os.path.isfile(R0_file):
        return [R0_file]
    else:
        msg = 'Cannot find art_illumina output files!'
        raise ValueError(msg)    

def get_total_seq_len(fasta_file):
    """Simple function that uses pyfastx to quickly read in a fasta,
    and then the sum of sequence lengths is returned
    """
    x = [len(seq) for h,seq in pyfastx.Fasta(fasta_file, build_index=False)]
    return sum(x)

def calc_fold(frag_fasta, seq_depth, total_barcodes, art_params, tmp_dir):
    """Calculate fold coverage to simulate per barcode
    Parameters
    ----------
    frag_fasta : str
        Fasta file of simulated gDNA fragments for a particular barcode
    seq_depth : int
        Total sequencing depth
    total_barcodes : int
        Total number of barcodes
    art_params : dict
        art_illumina parameters
    tmp_dir : str
        temp_dir
    """
    # fold coverage
    ## total length of fragments
    total_frag_len = get_total_seq_len(frag_fasta)
    ## length of reads
    try:
        read_len = int(float(art_params['--len']))
    except KeyError:
        raise KeyError('Cannot find --len param for art_illumina')
    try:
        _ = art_params['--paired']
        is_paired = True
    except KeyError:
        is_paired = False
    read_len = read_len * 2 if is_paired else read_len
    ## fold calc
    seqs_per_barcode = seq_depth / float(total_barcodes)
    fold  = seqs_per_barcode * read_len / total_frag_len
    # return
    return fold

def combine_frag_tsv(tsv_files, output_dir, debug=False):
    """ Concat all genome fragment tsv files
    Parameters
    ----------
    tsv_files : iterable
        Iterable of fragment tsv files
    output_dir : str
        Final output directory for read files
    debug : bool
        Debug mode
    """
    # writing files
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    out_file = os.path.join(output_dir, 'fragments.tsv')
    ## header
    with open(out_file, 'w') as outF:
        outF.write('\t'.join(['Barcode', 'Frag_ID',
                              'Genome', 'Contig',
                              'Frag_start', 'Frag_end']) + '\n')
    ## body
    for F in list(set(tsv_files)):
        with open(F, 'r') as inF, open(out_file, 'a') as outF:
            for line in inF:
                outF.write(line)
        # deleting temp file
        try:
            os.unlink(F)
        except OSError:
            logging.warning('Could not remove temp file: {}'.format(F))
    logging.info('File written: {}'.format(out_file))

def combine_reads(fq_files, output_dir, name_fmt='{readID} BX:Z:{barcodeID}',
                  seq_depth=None, debug=False):
    """ Concat all read files
    Parameters
    ----------
    fq_files : iterable
        Iterable of fastq read files
    output_dir : str
        Final output directory for read files
    name_fmt : str
        How to format the read name. "{readID} BX:Z:{barcodeID}" will generate names such as "@ST-J00101:121:HYCGGBBXX:5:1101:30533:1191 BX:Z:A58B91C07D86"
    debug : bool
        Debug mode
    """    
    # split by read pairs
    R1_files = [x[0] for x in fq_files]
    try:
        R2_files = [x[1] for x in fq_files]
    except IndexError:
        R2_files = None
        
    # writing files
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    ## read1
    R1_files = _combine_reads(R1_files, output_dir, 'R1.fq', name_fmt, seq_depth)
    ## read2
    if R2_files is not None:
        R2_files = _combine_reads(R2_files, output_dir, 'R2.fq', name_fmt, seq_depth)

def fix_nucs(seq):
    ambig = {'A', 'T', 'G', 'C'}
    seq = ['N' if x not in ambig else x for x in seq.rstrip()]
    return ''.join(seq) + '\n'
        
def _combine_reads(read_files, out_dir, out_file, name_fmt, seq_depth=None):
    """Combining temporary read files.
    eg., @ST-J00101:121:HYCGGBBXX:5:1101:30533:1191 BX:Z:A58B91C07D86
    """
    if seq_depth is not None:
        cnt = 0
    out_file = os.path.join(out_dir, out_file)
    with open(out_file, 'w') as outF:
        for i,F in enumerate(read_files):
            if seq_depth is not None and cnt > seq_depth:
                break
            with open(F, 'r') as inF:
                seq_len = 0
                qual_len = 0
                for ii,line in enumerate(inF):
                    if ii % 4 == 0:   # read header
                        if seq_depth is not None:
                            cnt += 1
                            if cnt > seq_depth:
                                break
                        line = line.split('-')
                        X = name_fmt.format(readID=line[0], barcodeID=line[1])
                        outF.write(X + '\n')
                        seq_len = 0
                        qual_len = 0
                        continue
                    elif ii % 4 == 1:  # read sequence
                        seq_len = len(line.rstrip())
                        line = fix_nucs(line)
                    elif ii % 4 == 3:  # read quality
                        qual_len = len(line.rstrip())
                        if seq_len != qual_len:
                            msg = 'Line {} in file {}: read-len != qual-len'
                            raise ValueError(msg.format(ii, F))
                    outF.write(line)
            # deleting temp file
            try:
                os.unlink(F)
            except OSError:
                logging.warning('Could not remove temp file: {}'.format(F))
                        
    logging.info('File written: {}'.format(out_file))
        
def select_refs(n_frags, genome_table):
    """Selecting genomes (with replacement)
    -- Columns in genome_table --
    "Fasta" = genome file
    "Taxon" = taxon to select
    "Perc_rel_abund" = weight
    """
    probs = genome_table.Perc_rel_abund / np.sum(genome_table.Perc_rel_abund)
    refs = np.random.choice(np.arange(genome_table.shape[0]),
                            size=n_frags,
                            replace=True,
                            p=probs)
    return genome_table.loc[refs]
    
def sim_frags(refs, n_frags, frag_size_loc, frag_size_sd,
              frag_size_min, frag_size_max):
    """Simulating fragment sizes
    Parameters
    ----------
    refs : pd.dataframe
        Dataframe of reference taxa
    n_frags : int
        Number of fragments
    frag_size_loc : float
        Mean fragment size
    frag_size_sd : float
        Stdev fragment size
    frag_size_min : float
        Min fragmen size (trun-norm distribution)
    frag_size_max : float
        Max fragmen size (trun-norm distribution)
    """
    a = (frag_size_min - frag_size_loc) / frag_size_sd  # fraction of range smaller than loc
    b = (frag_size_max - frag_size_loc) / frag_size_sd  # fraction of range greater than loc
    refs.loc[:,'Frag_size'] = truncnorm.rvs(a, b,
                                            size=n_frags,
                                            loc=frag_size_loc,
                                            scale=frag_size_sd)
    return refs

def read_fasta(fasta_file):
    """Fast reading into memory of fasta file
    return: dict
    """
    seqs = {h:seq for h,seq in pyfastx.Fasta(fasta_file, build_index=False)}
    return seqs
    
def parse_frags(refs, barcode, outFr, outFt):
    """Parsing fragment from a genome and writing them to a file.
    Giving each simulated fragment a UUID. 

    Parameters
    ----------
    refs : pd.dataframe
        Dataframe of reference taxa
    barcode : str
        Barcode ID
    outFr : writable file handle 
        Read sequences written to the file
    outFt : writable file handle 
        Frag sequence info written to the file
    """
    regex = re.compile(r'[ ;:,(){|]+')
    # fasta file of genome
    taxon = refs['Taxon'].unique()[0]
    taxon = regex.sub('_', taxon)
    fasta_file = refs['Fasta'].unique()[0]
    # loading fasta
    f = read_fasta(fasta_file)
    contig_ids = list(f.keys())
    for x in contig_ids:
        assert('|' not in x)
        
    contigs = []
    for idx,row in refs.iterrows():
        frag_size = int(row['Frag_size'])
        shuffle(contig_ids)
        for contig_id in contig_ids:
            contig_len = len(f[contig_id])
            if contig_len >= frag_size:
                frag_uuid = str(uuid.uuid4()).replace('-', '')
                frag_max_end = contig_len - frag_size
                frag_max_end = 1 if frag_max_end < 1 else frag_max_end
                frag_start = np.random.randint(0, frag_max_end)
                frag_end = frag_start + frag_size
                ## writing sequence
                outFr.write('>{}\n{}\n'.format(frag_uuid,
                                              f[contig_id][frag_start:frag_end]))
                contigs.append(contig_id)
                # writing tsv of positions
                outFt.write('\t'.join([str(barcode),
                                       str(frag_uuid),
                                       str(taxon),
                                       str(contig_id),
                                       str(frag_start),
                                       str(frag_end)]) + '\n')
                break
                
    # assert that frag was created
    if len(contigs) < refs.shape[0]:
        msg = 'For taxon "{}", cannot find contig with length >= fragment size ({})'
        raise ValueError(taxon, msg.format(frag_size))

def barcodes(n):
    """Return an array (len=n) of barcode IDs
    # barcode naming: 
      The 12-bp long i5-barcode is composed of 5bp "A-part"
      connected to a 6bp "C-part" with a constant 1bp A-handle.
      The 13-bp long i7-barcode is composed of 6bp "D-part" 
      connected to a 6bp "B-part" with a constant 1bp C-handle.
      In total we have 96 different A, B, C, D parts = 85 million combinations.
    Parameters
    ----------
    n : int
        number of barcodes
    """
    func = np.vectorize(_barcode_ids)
    return func(np.arange(n))

def _barcode_ids(start=1, end=96):
    x = _barcode_id('A') + _barcode_id('B') + _barcode_id('C') + _barcode_id('D')
    return x

def _barcode_id(Char, start=1, end=96):
    x = Char + '{0:02d}'.format(np.random.randint(start,end,dtype='int'))
    return x
