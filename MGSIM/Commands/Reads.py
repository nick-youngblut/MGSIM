#!/usr/bin/env python

"""
reads: simulating reads

Usage:
  reads [options] <abund_table>
  reads -h | --help
  reads --version

Options:
  <abund_table>   Taxon abundance info.  
                  ('-' if from STDIN)
  -n=<n>          Number of cpus. 
                  [Default: 1]
  --debug         Debug mode (no multiprocessing).
  -h --help       Show this screen.
  --version       Show version.

Description:

"""

# import
## batteries
from docopt import docopt
import sys,os
import re
from functools import partial
import multiprocessing as mp
## 3rd party
from Bio import SeqIO



def main(args):
    pass    

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
