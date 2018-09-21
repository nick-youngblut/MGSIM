#!/usr/bin/env python

"""
fragments: simulate genome fragments that would be found in a isopycnic gradient

Usage:
  fragments [options] <genomeList>
  fragments -h | --help
  fragments --version

Options:
  <genomeList>  A file listing: taxonName<tab>genomeSeqFileName
  --fp=<fp>     Full path to genomeSeqFiles (if not in genomeList file).
  --rtr=<rtr>   Read template length range (min,max).
                How big can the read template be (the part of the DNA fragment 
                that is actually sequenced)?
                [Default: 200,1200]
  --rtl=<rtl>   Read template length distribution (see Description).
                [Default: uniform,250,250]
  --nf=<nf>     Number of fragments to simulate per genome.
                If value ends in 'X', then the number of reads needed
                to obtain that coverage for each genome will be used
                [Default: 10000]
  --fld=<ld>    Fragment length distribution (see Description).
                [Default: skewed-normal,9000,2500,-5]
  --flr=<m>     Fragment length range (min,max). ('None' = no bounds)
                [Default: 4000,None]
  --fr=<fr>     Fasta of forward & reverse primers (if simulating amplicons).
  --np=<np>     Number of genomes to process in parallel.
                [Default: 1]
  --tbl         Write fragments as a table instead of a pickled python object.
  --MFE=<mfe>   MFEprimer executable. [Default: MFE_primer.py]
  --debug       Debug mode (turn off parallel processing)
  --version     Show version.
  -h --help     Show this screen.

Description:
  Simulate the genomic fragments that would be found in an isopycnic gradient.
  The location and G+C of each simulated fragment is written to a table.

  Genome sequence file names in the <genomeList> file should either include 
  the full path to the file or the path needs to be provided by the '--fp' flag.

  '--rtr' and '--rtl':
    'read template' refers to the template nucleotide molecule from
    which the read originated. Template size is needed to constrain the location
    and lengths of the simulated genome fragments that contain these templates
    (G+C calculated
    from these fragments).
    Templates could either be an amplicon (e.g. 16S rRNA sequencing)
    or a genomic fragment (e.g. shotgun metagenomics). '--rtl' constrains
    the size of this template. If primers are provided, the in-silico generated
    amplicons are filtered to those that fit in the specified read template
    range. '--rtl' is used to determine the template lengths of shotgun
    metagenomic reads (not contrained by provided primers).

    ** Example for 16S rRNA amplicons **
    Assumning amplicon lengths of 500-650 bp, use: '--rtr 500,650'.
    This will filter out all in-silico amplicons outside of this range.

  Selection of scaffold/chromosome
  --------------------------------
  For genoems with multiple scaffolds/chromosomes, the scaffold is chosen
  randomly if simulating shotgun fragments and selected randomly from the list
  of amplicons if simulating amplicon fragments.

  Distributions
  -------------
  * distribution functions from modules: np.random or skew_normal
  normal:
    Parameters: 
     * location (mean)
     * scale (standard deviation)

  uniform:
    Parameters: 
     * low
     * high
    
  skewed-normal:
    Parameters: 
     * location (mean)
     * scale (stardard deviation)
     * shape (skew)
    Example: --fld  skewed-normal,11000,1000,-1000

  truncated-normal:
    Parameters: 
     * location
     * scale
     * low
     * high

  Notes
  -----
  The simulated fragment size is constrained by the genome sequence
  template size. This may be why you get fragments that don't fall
  into the provided fragment length distribution.

  Output
  ------
  If --tbl: tab-delim file written to STDOUT, else: a pickled version of
  the table object

"""

# import
## batteries
from docopt import docopt
import sys,os
import re
## application libraries
from MGSIM import Fragments


def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
       
    sp = re.compile(' *, *')    
    args = {k:sp.split(str(v)) if sp.search(str(v)) 
            else v for k,v in args.items()}
    Fragments.main(args)
        
    

        
