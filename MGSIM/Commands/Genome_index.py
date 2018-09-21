#!/usr/bin/env python

"""
genome_index: index genomes for in-silico PCR; required for 
               amplicon-fragment simulation

Usage:
  genome_index [options] <genomeList>
  genome_index -h | --help
  genome_index --version

Options:
  <genomeList>    A file listing: taxonName<tab>genomeSeqFileName
  --fp=<fp>       Full path to genomeSeqFiles (if not in genomeList file).
  --K_value=<kv>  Kvalue for indexing. [Default: 9]
  --np=<np>       Number of genomes to process in parallel. [Default: 1]
  --twobit=<tb>   faToTwoBit executable. [Default: faToTwoBit]
  --fmt=<fm>      UniFastaFormat.py executable. [Default: UniFastaFormat.py]
  --idx=<id>      mfe_index_db.py executable. [Default: mfe_index_db.py]
  --quiet         Limit stderr output. 
  -h --help       Show this screen.
  --version       Show version.
  --debug         Debug mode (no parallel processes)

Description:
  The genomeList file (tab-delim; 'taxon_name<tab>file_name')
  should list all genome sequence files. Format:
   * fasta file format
   * 1 genome per file; 
   * genomes can be multi-chromosome/scaffold

  This script is not needed if shotgun-fragments are to be simulated
  (instead of amplicon-fragments).

"""

# import
## batteries
from docopt import docopt
## 3rd party
## application
from MGSIM import IndexGenomes
    

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    IndexGenomes.main(args)
   
