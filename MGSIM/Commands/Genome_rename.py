#!/usr/bin/env python

"""
genome_rename: formatting genome sequences in a multi-fasta file for MGSIM

Usage:
  genome_rename [options] <genome_fasta>...
  genome_rename -h | --help
  genome_rename --version

Options:
  <genome_fasta>  Fasta file(s) containing genome sequences.
                  ('-' if from STDIN)
  --prefix=<p>    Output file prefix. 
                  It must differ from the input file paths.
                  The path will be created if it doesn't exist.
                  [Default: .]
  -a=<a>          Number of ambiguous nucleotides allowed in a genome. 
                  [Default: 0]
  -n=<n>          Number of cpus. 
                  [Default: 1]
  --debug         Debug mode (no multiprocessing).
  -h --help       Show this screen.
  --version       Show version.

Description:
  Reformating the sequence names so that they work with MGSIM.
  The sequence names are reformated to conform to naming rules
  for the downstream analyses.
 
  Output
  ------
  Edited fasta written to STDOUT.

  WARNINGS!
  ---------
  There are no checks that the file is actually in fasta format.

"""

# import
## batteries
from docopt import docopt
import sys,os
import re
import logging
from functools import partial
import multiprocessing as mp
## 3rd party
from Bio import SeqIO
## logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def seq_rename(inFile, ambig_cutoff=0, prefix='.'):
    (inFileDir, inFileName) = os.path.split(inFile)
    inFileDir = os.path.abspath(inFileDir)
    prefix = os.path.abspath(prefix)
    outFile = os.path.join(prefix, inFileName)

    if inFileDir == prefix:
        msg = 'ERROR: input and output directories are the same'
        sys.stderr.write(msg + '\n')
        return None

    # regexes
    re0 = re.compile(r'.+ complete genome. ')
    re1 = re.compile(r'\W')
    re2 = re.compile(r'^_*(.*?)_*$')
    re3 = re.compile(r'_*complete_genome')
    re4 = re.compile(r'(.{78}).+')
    ambig_chars = re.compile(r'[RYSWKMBVDH]')

    # iterating through sequence
    ambig_cnt = 0
    with open(outFile, 'w') as outFH:
        for i,record in enumerate(SeqIO.parse(inFile, 'fasta')):
            name = record.name
            name = re0.sub('', name)
            name = re1.sub('_', name)
            name = re2.sub(r'>\1', name)
            name = re3.sub('', name)
            name = re4.sub(r'\1', name)
            name = name.lstrip('>') + '__seq{}'.format(i)
            record.id = name
            record.description = name
            ambig_cnt += len(ambig_chars.findall(str(record.seq.upper())))
            SeqIO.write(record, outFH, 'fasta')

    if ambig_cnt > int(ambig_cutoff):
        msg = 'No. of ambig chars ({}) is > than ambig-cutoff ({})'
        raise ValueError(msg.format(ambig_cnt, ambig_cutoff))
        
    msg = 'File written: {}\n'
    sys.stderr.write(msg.format(outFile))


def main(args):
    # input
    if args['<genome_fasta>'][0] == '-':
        args['<genome_fasta>'] = [x.rstrip() for x in sys.stdin]

    # output
    args['--prefix'] = os.path.abspath(args['--prefix'])
    if not os.path.isdir(args['--prefix']):
        os.makedirs(args['--prefix'])
        
    # rename
    if args['--debug']:
        for f in args['<genome_fasta>']:
            seq_rename(f, ambig_cutoff=args['-a'], prefix=args['--prefix'])
    else:
        p = mp.Pool(int(args['-n']))
        seq_rename_p = partial(seq_rename, ambig_cutoff=args['-a'], prefix=args['--prefix'])
        p.map(seq_rename_p, args['<genome_fasta>'])
    

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
