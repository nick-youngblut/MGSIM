#!/usr/bin/env python

"""
ht_reads: simulating haplotagging reads

Usage:
  ht_reads [options] <genome_table> <abund_table> <output_dir>
  ht_reads -h | --help
  ht_reads --version

Options:
  <genome_table>          Taxon genome info.
  <abund_table>           Taxon abundance info.
  <output_dir>            Output directory.
  --barcode-total=<bt>    Number of barcodes 
                          [Default: 1e3]
  --barcode-chunks=<bc>   Chunk size for parallel processing of barcodes
                          [Default: 1e2]
  --seq-depth=<d>         Number of (paired) Illumina reads per sample.
                          [Default: 1e5]
  --frag-size-mean=<fsm>  Mean fragment size of input gDNA (bp)
                          [Default: 10000]
  --frag-size-sd=<fss>    Fragment size standard deviation of input gDNA (bp)
                          [Default: 1000]
  --frag-size-min=<fsa>   Min fragment size of input gDNA (bp)
                          [Default: 8000]
  --frag-size-max=<fsa>   Min fragment size of input gDNA (bp)
                          [Default: 20000]
  --frag-bc-mean=<fbm>    Mean fragment-to-barcode count
                          [Default: 5]
  --frag-bc-sd=<fbs>      Standard deviation of fragment-to-barcode count
                          [Default: 1]
  --art-paired            art_illumina "paired"  parameter.
  --art-len=<al>          art_illumina "len" parameter.
                          [Default: 150] 
  --art-mflen=<am>        art_illumina "mflen" parameter.
                          Use 0 to turn off
                          [Default: 200]
  --art-sdev=<ad>         art_illumina "sdev" parameter.
                          [Default: 10]
  --art-seqSys=<as>       art_iilumina "seqSys" parameter.
                          [Default: HS25]
  --rndSeed=<rs>          Random Seed for Art. If None, then randomly set.
                          [Default: None]
  --tmp-dir=<td>          Temporary directory. 
                          [Default: .sim_reads]
  --read-name=<bft>       Read name format
                          [Default: {readID} BX:Z:{barcodeID}]
  -n=<n>                  Number of cpus. 
                          [Default: 1]
  --debug                 Debug mode (no subprocesses verbose output).
  -h --help               Show this screen.
  --version               Show version.

Description:
  Simulating reads for each taxon in each synthetic community

  abund_table
  -----------
  * tab-delimited
  * must contain 3 columns
    * "Community" = community ID (ie., sample ID)
    * "Taxon" = taxon name
    * "Perc_rel_abund" = percent relative abundance of the taxon
  * NOTE: For ht_reads, only the first community will be used (no multi-samples)

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  Output
  ------
  * tsv file of simulated genome fragments
    * These fragments were used as reference for read simulation
  * Read1 fastq [& Read2 fastq]
    * Default read header naming: <fragment UUID> BX:Z:<barcode_ID>
"""

# import
## batteries
from docopt import docopt
import sys,os
import re
import time
import logging
from pprint import pprint
from shutil import rmtree
from functools import partial
import multiprocessing as mp
from distutils.spawn import find_executable
## 3rd party
import numpy as np
## application
from MGSIM import SimReads
from MGSIM import SimHtReads
## logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def _flatten(list_of_lists):
    flattened_list = []
    for x in list_of_lists:
        for y in x:
            flattened_list.append(y)
    return flattened_list

def flatten(list_of_lists):
    flattened_list1 = []
    flattened_list2 = []
    for x in list_of_lists:
        for i,y in enumerate(x):
            if i % 2 == 0:
                flattened_list1.append(y)
            else:
                flattened_list2.append(y)
    return [_flatten(flattened_list1),
            _flatten(flattened_list2)]

def main(args):
    """ Main interface
    """
    # check executables
    exe = 'art_illumina'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))
    
    # load tables
    ## table mapping taxa to genomes
    logging.info('Loading genome table...')
    genome_table = SimReads.load_genome_table(args['<genome_table>'], nproc=args['-n'])
    ## table of taxon relative abundances
    logging.info('Loading taxon abundance table...')
    abund_table = SimReads.load_abund_table(args['<abund_table>'])

    # Sim Params
    ## random seed
    if args['--rndSeed'] is None or args['--rndSeed'] == 'None':
        rndSeed = None
    else:
        rndSeed = int(args['--rndSeed'])

    # simulating per community
    comms = abund_table.Community.unique()
    for comm in comms:
        logging.info('Simulating reads for community: {}'.format(comm))
        comm_tbl = abund_table.loc[abund_table.Community == comm,]
        sim_per_community(args, comm_id = comm,
                          abund_table = comm_tbl,
                          genome_table = genome_table,
                          rndSeed = rndSeed)

def sim_per_community(args, comm_id, abund_table, genome_table, rndSeed):
    # I/O
    outdir = os.path.join(args['<output_dir>'], str(comm_id))
    ## temp dir
    tmp_dir = os.path.join(args['--tmp-dir'], str(comm_id))
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    
    ## batching by barcodes
    ### creating barcode ID array
    n_barcodes = int(float(args['--barcode-total']))
    logging.info('Creating {} barcodes...'.format(n_barcodes))
    barcodes = SimHtReads.barcodes(n_barcodes)
    ### simulating barcodes
    func = partial(SimHtReads.sim_barcode_reads,
                   genome_table=genome_table,
                   abund_table=abund_table,
                   tmp_dir=tmp_dir,
                   args=args,
                   rndSeed=rndSeed,
                   debug=args['--debug'])
    barcodes = SimHtReads.bin_barcodes(barcodes, args['--barcode-chunks'])
    msg = 'Processing barcodes in {} chuncks of {} barcodes...'
    logging.info(msg.format(len(barcodes), args['--barcode-chunks']))
    if args['--debug'] or args['-n'] < 2:
        files = map(func, barcodes)
    else:
        Pool = mp.Pool(args['-n']) 
        files = Pool.map(func, barcodes)
        Pool.close()
    # flatten batched output
    (fq_files,tsv_files) = flatten(files)
    # combining all frag tsv
    logging.info('Combining all fragment info tables (n={})'.format(len(tsv_files)))
    SimHtReads.combine_frag_tsv(tsv_files, outdir)
    # combining all reads
    logging.info('Combining all read fastq files (n={})'.format(len(fq_files)))    
    SimHtReads.combine_reads(fq_files, outdir, name_fmt=args['--read-name'],
                             seq_depth=int(float(args['--seq-depth'])))
    # removing temp directory    
    if args['--debug'] is False:
        logging.info('Removing temporary directory: {}'.format(args['--tmp-dir']))
        rmtree(args['--tmp-dir'], ignore_errors=True)
        if os.path.isdir(args['--tmp-dir']):
            msg = 'Could not remove directory: {}'
            logging.warning(msg.format(args['--tmp-dir']))
            
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    args['-n'] = int(args['-n'])
    os.environ['OMP_NUM_THREADS'] = '1'
    main(args)
   
