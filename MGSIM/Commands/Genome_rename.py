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
from functools import partial
import multiprocessing as mp


def seq_rename(inFile, prefix='.'):
    (inFileDir, inFileName) = os.path.split(inFile)
    inFileDir = os.path.abspath(inFileDir)
    prefix = os.path.abspath(prefix)

    if inFileDir == prefix:
        msg = 'ERROR: input and output directories are the same'
        sys.stderr.write(msg + '\n')
        return None

    inFH = open(inFile, 'rb')             
    outFile = os.path.join(prefix, inFileName)
    outFH = open(outFile, 'wb')

    # regexes
    re0 = re.compile(r'.+ complete genome. ')  # added for tutorial
    re1 = re.compile(r'\W')
    re2 = re.compile(r'^_*(.*?)_*$')
    re3 = re.compile(r'_*complete_genome')
    re4 = re.compile(r'(.{78}).+')

    # iterating through file
    seq_names = dict()
    for line in inFH:
        line = line.rstrip()
        if line.startswith('>'):  # seq name
            line = re0.sub('', line)
            line = re1.sub('_', line)
            line = re2.sub(r'>\1', line)
            line = re3.sub('', line)
            line = re4.sub(r'\1', line)
            try:
                _ = seq_names[line]
                line = '_'.join(line, str(len(seq_names.keys())))
            except KeyError:
                pass
            seq_names[line] = 1                
        outFH.write(line + '\n')

    # end
    inFH.close()
    outFH.close()
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
            seq_rename(f, prefix=args['--prefix'])
    else:
        p = mp.Pool(int(args['-n']))
        seq_rename_p = partial(seq_rename, prefix=args['--prefix'])
        p.map(seq_rename_p, args['<genome_fasta>'])
    

def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
