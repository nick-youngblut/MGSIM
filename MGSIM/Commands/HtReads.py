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
  --art-paired            art_iilumina "paired  parameter.
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
  > NOTE: For HtReads, only the first community will be used (no multi-samples)

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  Output
  ------
  * A set of read files for each sample
    * directory structure: OUTPUT_DIR/COMMUNITY/read_files
    * read sequences are named by the taxon they originate from
"""

test = """
"""

# import
## batteries
from docopt import docopt
import sys,os
import re
from shutil import rmtree
from functools import partial
import multiprocessing as mp
from distutils.spawn import find_executable
## 3rd party
import numpy as np
## application
from MGSIM import SimReads
from MGSIM import SimHtReads

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
    # checking params
    if args['--art-paired'] is True:
        args['--art-paired'] = ''
    else:
        args.pop('--art-paired', None)
    if int(args['--art-mflen']) <= 0:
        args.pop('--art-mflen', None)

    ## temp dir
    if not os.path.isdir(args['--tmp-dir']):
        os.makedirs(args['--tmp-dir'])
    ## check that simulator exists
    exe = 'art_illumina'
    if find_executable(exe) == '':
        raise IOError('Cannot find executable: {}'.format(exe))
    
    # load tables
    ## table mapping taxa to genomes
    genome_table = SimReads.load_genome_table(args['<genome_table>'])
    ## table of taxon relative abundances
    abund_table = SimReads.load_abund_table(args['<abund_table>'])
    x = abund_table.Community[0]
    abund_table = abund_table.loc[abund_table.Community == x,]

    ## batching by barcodes
    n_barcodes = int(float(args['--barcode-total']))
    barcodes = np.arange(n_barcodes)   # numeric array of barcodes
    Pool = mp.Pool(int(args['-n']))
    func = partial(SimHtReads.sim_barcode_reads,
                   genome_table=genome_table,
                   abund_table=abund_table,
                   args=args)
    barcodes = SimHtReads.bin_barcodes(barcodes,args['--barcode-chunks'])
    if args['--debug']:
        ret = map(func, barcodes)
    else:
        ret = Pool.map(func, barcodes)
    ret = list(ret)
    # combine all reads
    SimHtReads.combine_reads(args['--tmp-dir'], '', args['<output_dir>'], n_barcodes)
    # remove temp directory
    if args['--debug'] is None:
        rmtree(args['--tmp-dir'])
        

    # TODO: direct temp file following and deletion
        
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
