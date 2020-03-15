#!/usr/bin/env python

"""
genome_download: downloading genomes

Usage:
  genome_download [options] <accession_table>
  genome_download -h | --help
  genome_download --version

Options:
  <accessin_table>  Taxon-accession table (see Description).
                    Use '-' if from STDIN.
  -d=<d>            Output directory. [Default: .]
  -e=<e>            Email to use for NCBI queries. [Default: blank@gmail.com]
  -a=<a>            Number of ambiguous nucleotides allowed in a genome. [Default: 0]
  -n=<n>            Number of cpus. [Default: 1]
  -t=<t>            Number of tries to download genomes. [Default: 10]
  -r                Rename genome sequences based on taxon name?
  --debug           Debug mode (no multiprocessing).
  -h --help         Show this screen.
  --version         Show version.

Description:
  Taxon-accession table
  ---------------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Accession" = NCBI accession used for downloading 
      * Possible accessions:
         * ncbi nucleotide db
         * ncbi assembly db
         * ftp url to genome (direct download)
  * other columns are allowed

  Output
  ------
  * Genome fasta files written to the specified output directory
  * A table mapping taxa to the download genome fasta file is written to STDOUT
"""

# import
import sys,os
import logging
## batteries
from docopt import docopt
from MGSIM import Genome_Download
## logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


# opt parse
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    Genome_Download.main(args)
   
