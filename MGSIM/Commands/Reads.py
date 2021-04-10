#!/usr/bin/env python

"""
reads: simulating reads

Usage:
  reads [options] <genome_table> <abund_table> <output_dir>
  reads -h | --help
  reads --version

Options:
  <abund_table>       Taxon abundance info.  
  <genome_table>      Taxon genome info.
  <output_dir>        Output directory.
  --sr-seq-depth=<d>  Number of (paired) Illumina reads per sample.
                      [Default: 1e5]
  --art-paired        art_iilumina --paired parameter.
  --art-len=<al>      art_illumina --len parameter.
                      [Default: 150] 
  --art-mflen=<am>    art_illumina --mflen parameter.
                      Use 0 to turn off
                      [Default: 200]
  --art-sdev=<ad>     art_illumina --sdev parameter.
                      [Default: 10]
  --art-seqSys=<as>   art_illumina --seqSys parameter.
                      [Default: HS25]
  --pb-seq-depth=<d>  Number of PacBio reads per sample.
                      [Default: 0]
  --pb-min-len=<ml>   Minimum read length for PacBio reads.
                      [Default: 500]  
  --np-seq-depth=<d>  Number of Nanopore reads per sample.
                      [Default: 0]
  --np-min-len=<ml>   Minimum read length for Nanopore reads.
                      [Default: 500]  
  --np-error=<ml>     NanoSim-H error profile
                      [Default: 'ecoli_R9_1D']  
  --tmp-dir=<td>      Temporary directory
                      [Default: .sim_reads]
  --rndSeed=<rs>      Random Seed for Art. If None, then randomly set.
                      [Default: None]
  -n=<n>              Number of cpus. 
                      [Default: 1]
  --debug             Debug mode (no subprocesses; verbose output)
  -h --help           Show this screen.
  --version           Show version.

Description:
  Simulating reads for each taxon in each synthetic community.
  Read types (simulator): 
    * Illumina (ART)
    * PacBio (SimLord)
    * Nanopore (NanoSim-H)

  abund_table
  -----------
  * tab-delimited
  * must contain 3 columns
    * "Community" = community ID (ie., sample ID)
    * "Taxon" = taxon name
    * "Perc_rel_abund" = percent relative abundance of the taxon

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

# import
## batteries
from docopt import docopt
import sys,os
import re
import logging
from functools import partial
import multiprocessing as mp
import logging
## application
from MGSIM import SimReads
## logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def main(args):
    # load tables
    genome_table = SimReads.load_genome_table(args['<genome_table>'])
    abund_table = SimReads.load_abund_table(args['<abund_table>'])
    # create pairwise
    sample_taxon = SimReads.sample_taxon_list(genome_table,
                                              abund_table)
    # simulating reads
    ## art params
    art_params = {'--paired' : args['--art-paired'],
                  '--len' : int(args['--art-len']),
                  '--mflen' : float(args['--art-mflen']),
                  '--sdev' : float(args['--art-sdev']),
                  '--seqSys' : args['--art-seqSys']}
    if art_params['--paired'] is True:
        art_params['--paired'] = ''
    else:
        art_params.pop('--paired', None)
    if int(art_params['--mflen']) <= 0:
        art_params.pop('--mflen', None)
    ## simlord params
    sl_params = {'--min-readlength' : args['--pb-min-len']}
    ## nanosim-h params
    ns_params = {'--min-len' : args['--np-min-len'],
                 '--profile' : args['--np-error']}
    ### random seed
    if args['--rndSeed'] is None or args['--rndSeed'] == 'None':
        rndSeed = None
    else:
        rndSeed = int(args['--rndSeed'])
    ## read simulate
    ### illumina
    if float(args['--sr-seq-depth']) > 0:
        SimReads.sim_illumina(sample_taxon,
                              output_dir=args['<output_dir>'],
                              seq_depth=float(args['--sr-seq-depth']),
                              art_params=art_params,
                              temp_dir=args['--tmp-dir'],
                              nproc=int(float(args['-n'])),
                              rndSeed=rndSeed,
                              debug=args['--debug'])
    ### pacbio
    if float(args['--pb-seq-depth']) > 0:
        SimReads.sim_pacbio(sample_taxon,
                            output_dir=args['<output_dir>'],
                            seq_depth=float(args['--pb-seq-depth']),
                            sl_params=sl_params,
                            temp_dir=args['--tmp-dir'],
                            nproc=int(float(args['-n'])),
                            debug=args['--debug'])
    ### nanopore
    if float(args['--np-seq-depth']) > 0:
        SimReads.sim_nanopore(sample_taxon,
                              output_dir=args['<output_dir>'],
                              seq_depth=float(args['--np-seq-depth']),
                              ns_params=ns_params,
                              temp_dir=args['--tmp-dir'],
                              nproc=int(float(args['-n'])),
                              debug=args['--debug'])
    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
