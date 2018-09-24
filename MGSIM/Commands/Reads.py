#!/usr/bin/env python

"""
reads: simulating reads

Usage:
  reads [options] <genome_table> <abund_table>
  reads -h | --help
  reads --version

Options:
  <abund_table>       Taxon abundance info.  
  <genome_table>      Taxon genome info.
  --sr-seq-depth=<d>  Number of (paired) Illumina reads per sample.
                      [Default: 1e6]
  --art-params=<ap>   art parameters
                      [Default: '-ss HS25 -p -l 150 -m 200 -s 10']
  --tmp-dir=<td>      Temporary directory
                      [Default: '.sim_reads']
  -n=<n>              Number of cpus. 
                      [Default: 1]
  --debug             Debug mode (no multiprocessing).
  -h --help           Show this screen.
  --version           Show version.

Description:

"""

# import
## batteries
from docopt import docopt
import sys,os
import re
from functools import partial
import multiprocessing as mp
## application
from MGSIM import SimReads

def main(args):
    # load genome_table
    # load genome_abundance info
    # create pairwise [sample,taxon,genome_fasta] list of lists
    # (parallel) call art_illumina
    ## convert abund to genome fold
    ### fold = abund * seq_depth / genome_length
    ## art_illumina -ss HS25 -p -l 150 -f 30 -m 200 -s 10 -f {fold} -i {input} -o {params.out_prefix}

    # load tables
    genome_table = SimReads.load_genome_table(args['<genome_table>'])
    abund_table = SimReads.load_abund_table(args['<abund_table>'])
    # create pairwise
    sample_taxon = SimReads.sample_taxon_list(genome_table,
                                              abund_table)
    # simulating reads
    SimReads.sim_illumina(sample_taxon,
                          seq_depth=float(args['--sr-seq-depth']),
                          art_params=args['--art-params'],
                          temp_dir=args['--tmp-dir'],
                          nproc=int(args['-n']),
                          debug=args['--debug'])
    
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
