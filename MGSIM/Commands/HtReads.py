#!/usr/bin/env python

"""
htreads: simulating haplotagging reads

Usage:
  htreads [options] <genome_table> <abund_table> <output_dir>
  htreads -h | --help
  htreads --version

Options:
  <genome_table>          Taxon genome info.
  <abund_table>           Taxon abundance info.
  <output_dir>            Output directory.
  --barcode-total=<bt>    Number of barcodes 
                          [Default: 1e2]
  --barcode-chunks=<bc>   Chunk size for parallel processing of barcodes
                          [Default: 1e1]
  --seq-depth=<d>         Number of (paired) Illumina reads per sample.
                          [Default: 1e5]
  --frag-size-mean=<fsm>  Mean fragment size of input gDNA (bp)
                          [Default: 10000]
  --frag-size-sp=<fss>    Fragment size spread of input gDNA (bp)
                          [Default: 10000]
  --frag-size-min=<fsa>   Min fragment size of input gDNA (bp)
                          [Default: 8000]
  --frag-size-max=<fsa>   Min fragment size of input gDNA (bp)
                          [Default: 20000]
  --frag-bc-mean=<fbm>    Mean fragment-to-barcode count
                          [Default: 5]
  --frag-bc-sd=<fbs>      Stdev fragment-to-barcode count
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
  --tmp-dir=<td>          Temporary directory. 
                          [Default: .sim_reads]
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
  * NOTE: For HtReads, only the first community will be used (no multi-samples)

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
    * read header naming: <fragment UUID> BX:Z:<barcode_ID>
"""

# import
## batteries
from docopt import docopt
import sys,os
import re
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
    """Notes on haplotagging read simulation
    # algorithm
    * for each barcode batch (eg., 1000 barcodes); (in parallel)
      * for each barcode-ID
        * create temporary directory
          * all temp files in this directory
        * create reference for simulation
          * select genomes (with replacement)
             * weighted by taxon abundances (DNA pool)
             * see https://stackoverflow.com/questions/10803135/weighted-choice-short-and-simple
          * for each genome selection (possible duplicate genomes):
            * random selection of genome region
              * package: pyfasta             
              * fragment length based on tuncnorm distribution
                 * https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.truncnorm.html#scipy.stats.truncnorm
            * output:
              * multi-fasta file of fragments 
                 * ><genome>__start<start>-end<end>
              * table: <barcode> <genome> <start> <end>
        * simulate reads on reference
          * number_of_reads = seq_depth / barcodes
          * art_illumina on multi-fasta file
            * `--id` => read prefix
            * `--out` => barcode
        * combine read files (by batch)
      * combine batch read files
    """
    # check executables
    exe = 'art_illumina'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))
    ## temp dir
    if not os.path.isdir(args['--tmp-dir']):
        os.makedirs(args['--tmp-dir'])
    
    # load tables
    ## table mapping taxa to genomes
    genome_table = SimReads.load_genome_table(args['<genome_table>'])
    ## table of taxon relative abundances
    abund_table = SimReads.load_abund_table(args['<abund_table>'])
    x = abund_table.Community[0]
    abund_table = abund_table.loc[abund_table.Community == x,]

    ## batching by barcodes
    ### creating barcode ID array
    n_barcodes = int(float(args['--barcode-total']))
    barcodes = SimHtReads.barcodes(n_barcodes)
    ### simulating barcodes
    func = partial(SimHtReads.sim_barcode_reads,
                   genome_table=genome_table,
                   abund_table=abund_table,
                   args=args)
    barcodes = SimHtReads.bin_barcodes(barcodes,args['--barcode-chunks'])
    if args['--debug']:
        files = map(func, barcodes)
    else:
        Pool = mp.Pool(int(args['-n']))
        files = Pool.map(func, barcodes)
    # flatten batched output
    (fq_files,tsv_files) = flatten(files)
    # combining all frag tsv
    SimHtReads.combine_frag_tsv(tsv_files, args['<output_dir>'])
    # combining all reads
    SimHtReads.combine_reads(fq_files, args['<output_dir>'])
    # removing temp directory
    if args['--debug'] is False:
        rmtree(args['--tmp-dir'])

        
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
