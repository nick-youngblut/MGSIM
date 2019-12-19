#!/usr/bin/env python

# import
## batteries
import os
import sys
## 3rd party
from docopt import docopt
#import Utils
## application
from MGSIM.Commands import Communities
from MGSIM.Commands import Genome_download
from MGSIM.Commands import Genome_rename
from MGSIM.Commands import Reads
from MGSIM.Commands import HtReads

def main(args=None):
    """Main entry point for application
    """
    if args is None:
        args = sys.argv[1:]
    
    docs = """
MGSIM: simulate metagenomes for multiple synthetic communities

Usage:
  MGSIM <command> [<args>...]
  MGSIM -l | --list
  MGSIM -h | --help
  MGSIM --version

Options:
  -l --list     List subcommands.
  -h --help     Show this screen.
  --version     Show version.

Commands:
  Use the `list` option.
Description:
  Simulate metagenomes for multiple synthetic communities.
  See the sub-command documentation for more information on features.
    """
    # arg parse
    args = docopt(docs,
                  version='0.1',
                  options_first=True)

    # dict of all subcommands
    cmds = {'communities' : Communities,
            'genome_download' : Genome_download,
            'genome_rename' : Genome_rename,
            'reads' : Reads,
            'ht_reads' : HtReads}
    
    # list subcommands
    if args['--list']:
        cmd_list = '\n'.join(sorted(cmds.keys(), key=str.lower))
        print('#-- Commands --#')
        print(cmd_list)
        exit()

    # running subcommand    
    try:
        func = cmds[args['<command>']]
    except KeyError:
        msg = 'ERROR: command "{}" does not exist'
        print(msg.format(args['<command>']))
        exit()        
    func.opt_parse(args['<args>'])

    
if __name__ == '__main__':
    main()
