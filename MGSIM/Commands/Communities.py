#!/usr/bin/env python
"""
communities: simulate taxon abundances in synthetic communities

Usage:
  communities [options] <genomeList> <prefix>
  communities -h | --help
  communities --version

Options:
  <genomeList>        GenomeList table (see Description)
                      Use '-' if from STDIN.
  <prefix>            Output file prefix
  --n-comm=<nc>       Number of communities to simulate.
                      [default: 1]
  --richness=<r>      The number of taxa in each library.
                      Values of 0-1 will be interpreted as a fraction of the 
                      total taxon pool.
                      [default: 1]
  --abund-dist=<a>    The statistical distribution used for selecting relative
                      abundances.
                      (see numpy.random for a list of distributions).
                      [default: lognormal]
  --abund-dist-p=<p>  Abundance distribution parameters.
                      (see numpy.random for distribution params).
                      [default: mean:10,sigma:2]
  --shared-perc=<sp>  The percent of taxa shared in each community.
                      Percent set by the community with the smallest richness.
                      Example: if smallest community is 10 taxa,
                               a shared percent of 20% = 2 shared taxa.
                      The total taxon pool must be large enough to accommodate
                      all un-shared taxa.
                      [default: 100]
  --perm-perc=<pp>    How much to vary the rank-abundances between communities. 
                      Percentage = percent of taxa to permute.
                      [default: 0]
  --config=<c>        Config file for setting community-specific parameters
                      (& global params).
                      Community-specific parameters can include:
                      (n_comm, richness, abund_dist)
  --debug             Debug mode
  --version           Show version.
  -h --help           Show this screen.

Description:
  Simulating the alpha and and beta diversities of >=1 community.

  genomeList
  ----------
    A text file listing one taxon name per line. 
    Anything after a <tab> will be ignored. 

  Output
  ------
    A tab-delimited table of taxon abundances for each library is written to
    STDOUT.
"""

# import
## batteries
from docopt import docopt
import sys,os
## application
from MGSIM.SimComms import SimComms


# functions
def main(uargs):
    # init    
    SC = SimComms(taxon_list = uargs['<genomeList>'],
                  perm_perc = uargs['--perm-perc'],
                  shared_perc = uargs['--shared-perc'],
                  richness = uargs['--richness'],
                  abund_dist = uargs['--abund-dist'],
                  abund_dist_params = uargs['--abund-dist-p'],
                  n_comm = uargs['--n-comm'],
                  config = uargs['--config'])
    
    # making communities
    for comm_id in SC.keys():        
        SC.make_comm(comm_id)

    # permuting based on permute perc
    for comm_id, comm in SC.items():
        if comm.richness > 1:
            SC.permute(comm, uargs['--perm-perc'])

    # writing out abundance table
    out_file = uargs['<prefix>'] + '_abund.txt'
    SC.write_comm_table(out_file)
            
    # creating abundance weighted by genome length
    SC.weighted_abundances()

    # writing out weighted abundance tables
    out_file = uargs['<prefix>'] + '_wAbund.txt'
    SC.write_comm_table(out_file)

        
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
