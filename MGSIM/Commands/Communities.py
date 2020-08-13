#!/usr/bin/env python
"""
communities: simulate taxon abundances in synthetic communities

Usage:
  communities [options] <genome_table> <prefix>
  communities -h | --help
  communities --version

Options:
  <genome_table>      Genome table (see Description)
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
  --beta-div=<m>...   beta-diversity measures to calculate among communities.
                      See scipy.spatial.distance for options.
                      [default: braycurtis jaccard]
  --config=<c>        Config file for setting community-specific parameters
                      (& global params).
                      Community-specific parameters can include:
                      (n_comm, richness, abund_dist)
  --rnd-seed=<rs>     Seed for reprodicible simulations.
                      [default: None]
  --debug             Debug mode
  --version           Show version.
  -h --help           Show this screen.

Description:
  Simulating taxon abundances (& the relative contribution to the DNA pool)
  for each synthetic community. 

  Community alpha- and beta-diveristy can be manipulated.

  genome_table
  ------------
  * tab-delimited
  * must contain 2 columns
    * "Taxon" = taxon name
    * "Fasta" = genome fasta file path
  * other columns are allowed

  Output
  ------
  * "PREFIX_abund.txt" = taxon relative abundances for each community
    * this is the relative number of genome copies for each taxon
  * "PREFIX_wAbund.txt" = taxon relative abundances weighted by genome size
    * this is the fraction of the DNA pool for each taxon
  * "PREFIX_beta-div.txt" = beta diversity among communities
    * see the beta-div parameter for selecting beta-diversity measures
"""

# import
## batteries
from docopt import docopt
import sys,os
import logging
## application
from MGSIM.SimComms import SimComms
## logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# functions
def main(uargs):
    # init
    SC = SimComms(taxon_list = uargs['<genome_table>'],
                  perm_perc = uargs['--perm-perc'],
                  shared_perc = uargs['--shared-perc'],
                  richness = uargs['--richness'],
                  abund_dist = uargs['--abund-dist'],
                  abund_dist_params = uargs['--abund-dist-p'],
                  n_comm = uargs['--n-comm'],
                  config = uargs['--config'],
                  rnd_seed = uargs['--rnd-seed'])
    
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

    # beta-diversity among communities
    out_file = uargs['<prefix>'] + '_beta-div.txt'
    SC.beta_diversity(measures=uargs['--beta-div'], outfile=out_file)
    
        
def opt_parse(args=None):
    if args is None:        
        args = docopt(__doc__, version='0.1')
    else:
        args = docopt(__doc__, version='0.1', argv=args)
    main(args)
   
