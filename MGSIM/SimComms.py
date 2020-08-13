from __future__ import print_function
# import
## batteries
import sys, os
import numpy as np
import pandas as pd
try:
    from StringIO import StringIO
except ModuleNotFoundError:
    from io import StringIO
import random
import re
import logging
from functools import partial
from itertools import chain,combinations
from operator import itemgetter
from collections import OrderedDict,defaultdict
## 3rd party
import scipy.stats as stats
from scipy.spatial import distance
from configobj import ConfigObj, flatten_errors
from Bio import SeqIO


# utility functions
def str2dict(s):
    """Parsing string (format: 'item:value,item:value')
    to create a dict object.    
    """    
    if hasattr(s, 'split'):
        l = re.split('[:,]', s)
        try:
            return {k.lower():float(v) for k,v in zip(l[0::2],l[1::2])}
        except TypeError:
            msg = 'distribution parameter values must be ints or floats.'
            raise TypeError(msg)
    else:
        return s
                
def random_insert_seq(l, seq):
    """Insert seq items at random locations in list.

    Paramters
    ---------
    l : list
        target list object
    seq : iteralbe
        items to insert in the list

    Returns
    -------
    list : target list with `seq` values randomly inserted
    """
    insert_locs = random.sample(range(len(l) + len(seq)), len(seq))
    inserts = dict(zip(insert_locs, seq))
    inputs = iter(l)
    return [inserts[pos] if pos in inserts else next(inputs)
            for pos in range(len(l) + len(seq))]

def power_neg(*args, **kwargs):
    return 1 - np.random.power(*args, **kwargs)

    
    
class _Comm(object):
    """Parent class for other classes in the module.
    """
    def __init__(self):
        pass
    
    # property/setter
    @property
    def abund_dist(self):
        return self._abund_dist
    @abund_dist.setter
    def abund_dist(self, x):
        self._abund_dist = str(x)
    @property
    def richness(self):
        return self._richness
    @richness.setter
    def richness(self, x):
        x = float(x)
        if x <= 1:
            # x = fraction of total taxa pool
            # setting x as that fraction number of taxa
            x = len(self.taxon_pool) * x
        self._richness = int(round(x,0))
        if self._richness < 1:
            self._richness = 1


class SimComms(_Comm):
    """Class for simulating taxon count data of communities.
    """
    def __init__(self, taxon_list, perm_perc, shared_perc,
                 richness, abund_dist, abund_dist_params,
                 n_comm, config=None, rnd_seed=None,
                 *args, **kwargs):
        """
        Parameters
        ----------
        See gradientComms
        """
        _Comm.__init__(self, *args, **kwargs)

        self._load_taxon_list(taxon_list)        
        self.perm_perc = perm_perc
        self.shared_perc = shared_perc
        self.richness = richness
        self.abund_dist = abund_dist
        self.abund_dist_params = str2dict(abund_dist_params)
        self.config = config
        self.n_comm = n_comm
        self.rnd_seed = rnd_seed

        if not self.rnd_seed is None and self.rnd_seed != 'None':
            random.seed(int(self.rnd_seed))
        
        # loading config; setting community parameters
        if config is not None:
            self.comm_params = self._load_config()
        else:
            self.comm_params = dict()
        self._set_comm_params()

        # lowering comm richness if taxon pool is limiting
        ## otherwise, shared_perc option throws and error
        if self.shared_perc < 100:
            self._lower_richness()

        # shared taxa
        self._set_shared_taxa()
                            
    def _get_configspec(self, strIO=True):
        """Return configspec set for instance.
        Parameters
        ----------
        strIO : bool
            return configspec as a StringIO instance
        Returns
        -------
        configspec object 
        """
        configspec = """
        [__many__]
            richness = float(0,inf, default=None)        
            abund_dist = string(default='exponential,1,0.5')
            start = float(default=None)
            end = float(default=None)
            loc = float(default=None)
            scale = float(default=None)
            sigma = float(default=None)
        """
        
        if strIO == True:
            return StringIO(configspec)
        else:
            return configspec
 
    def _load_config(self):    
        assert hasattr(self, 'config'), "No config attribute found."

        configspec = self._get_configspec()
        return ConfigObj(self.config, configspec=configspec)        
        
    def _set_comm_params(self):
        """Setting community-specific params including applying global params.
        """
        # adding to comm params if not enough set by config
        n_config_comms = len(self.comm_params.keys())
        n_diff = self.n_comm - n_config_comms
        
        for i in range(n_diff):
            self.comm_params[str(n_config_comms + i + 1)] = dict()
        
        for k,v in self.comm_params.items():
            # checking for params
            if ('richness' not in v.keys() or
                v['richness'] is None):
                v['richness'] = self.richness
            if ('abund_dist' not in v.keys() or
                v['abund_dist'] is None):
                v['abund_dist'] = self.abund_dist
            if ('abund_dist_p' not in v.keys() or
                v['abund_dist_p'] is None):
                v['abund_dist_p'] = self.abund_dist_params
            v['abund_dist_p'] = str2dict(v['abund_dist_p'])
                
    def _set_shared_taxa(self):
        """A list of taxa shared among all communities.
        The taxon list (pool) is reduced to just unshared taxa.
        """
        self.shared_taxa = self._drawFromTaxonPool(self.n_shared)
        
    def _load_taxon_list(self, fileName):
        """Loading taxon list file. Taxa order is randomly shuffled.
        
        Parameters
        ----------
        fileName : str
             name of taxon file
        """
        # load file
        if fileName == '-':
            inF = sys.stdin
        else:
            inF = open(fileName)
        df = pd.read_csv(inF, sep='\t')
        ## check headers
        diff = set(['Taxon','Fasta']) - set(df.columns.values)
        if len(diff) > 0:
            diff = ','.join(diff)
            raise ValueError('Cannot find table columns:{}'.format(diff))

        # setting taxon_pool
        self.taxon_pool = df['Taxon']
        random.shuffle(self.taxon_pool)

        # setting genome fasta dict {Taxon : Fasta}
        self.genome_fasta = {}
        for i,x in df.iterrows():
            self.genome_fasta[x['Taxon']] = x['Fasta']
        
    def _drawFromTaxonPool(self, n):
        """Draw from taxon pool, returning n-taxa;
        those taxa are removed from the pool.

        Parameters
        ----------
        n : int
            number of taxa to draw

        Returns
        -------
        list : [taxon_name1, taxon_nameN, ...]
        """
        assert n <= len(self.taxon_pool), \
            'Cannot draw {} taxa from taxon pool'.format(n)
        taxa = self.taxon_pool[:n]
        self.taxon_pool = self.taxon_pool[n:]
        return taxa


    def _lower_richness(self):
        """Lowering the richness of each community if the number of
        unique taxa (un-shared) + shared taxa is greater than the taxon
        pool from which to draw taxa.
        """
        rich = []
        for k,v in self.comm_params.items():
            try:
                rich.append(v['richness'])
            except KeyError:
                msg = 'Cannot find "richness" attribute for Comm {}'
                raise KeyError(msg.format(k))
        n_unique = np.sum([x - self.n_shared for x in rich])
        n_taxa_pool = len(self.taxon_pool)
        n_comm = len(rich)
        n_less = 0
        if n_unique + self.n_shared > n_taxa_pool:
            n_less = n_unique + self.n_shared - n_taxa_pool
            n_less = np.ceil(n_less/ n_comm) 
            n_less = int(n_less) + 1
        else:
            return 0
        for k,v in self.comm_params.items():
            new_rich = v['richness'] - n_less
            msg = 'WARNING: lowering richness ({} -> {}) for Community {}\n' + \
                  '  because the taxon pool is not large enough for the\n' + \
                  '  amount of un-shared taxa (set by --shared_perc)\n'
            sys.stderr.write(msg.format(v['richness'], new_rich, k))
            v['richness'] = new_rich
        self._n_shared = None
               
    def make_comm(self, comm_id):
        """Make a Comm object.
        
        Parameters
        ----------
        comm_id : str
            Community name from comm_params attrib
        """
        # assertions
        comm_id = str(comm_id)
        if not hasattr(self, 'comms'):
            self.comms = OrderedDict() #dict()
            
        try:
            self.comm_params[comm_id]
        except KeyError:
            raise KeyError('Cannot find community ID "{}" in '\
                           'community params\n'.format(comm_id))
        
        # init comm objects
        self.comms[comm_id] = Comm(comm_id, self)
        
    def write_comm_table(self, out_file=None, Long=True):
        """Joining comm objects into 1 dataframe and printing.
        Writing table to STDOUT.

        Parameters
        ----------
        out_file : str
            Output file name
        Long : bool
            Write table in long format
        """
        df =  pd.concat([x.taxa for x in self.values()],
                        axis=1, sort=False)

        write_index = True
        df.columns = self.keys()
        if Long == True:
            write_index = False
            # melting
            val_vars = list(df.columns)
            df['taxon'] = df.index
            df = pd.melt(df, id_vars=['taxon'], value_vars=val_vars)
            # ordering columns
            df.columns = ['taxon_name', 'library', 'rel_abund_perc']            
            # sorting 
            df.sort_values(by=['library','rel_abund_perc'], 
                           axis=0, ascending=[True,False], inplace=True)
            # converting any NAs to zeros
            df.fillna(0, inplace=True)
            # tidy taxon names
            func = lambda x: re.sub(r'[()\/:;, ]+', '_', x)
            df['taxon_name'] = df['taxon_name'].astype(str).apply(func)
            # getting rank by community (grouping by community)
            df['rank'] = df.groupby(['library'])['rel_abund_perc']\
                           .rank(method='first',ascending=False).astype('int')
            df = df[['library','taxon_name','rel_abund_perc','rank']]
            # renaming
            df.rename(columns={'library' : 'Community',
                       'taxon_name' : 'Taxon',
                       'rel_abund_perc' : 'Perc_rel_abund',
                       'rank' : 'Rank'}, inplace=True)
                        
        # writing dataframe
        if out_file is None:
            df.to_csv(sys.stdout, sep='\t', na_rep=0,
                      float_format='%.9f', index=write_index)
        else:
            df.to_csv(out_file, sep='\t', na_rep=0,
                      float_format='%.9f', index=write_index)            

    def weighted_abundances(self):
        """Calculate weighted abundance: rel_abund * genome_size
        
        Parameters
        ----------
        """
        self._set_genome_sizes()

        self.w_abund_tbl = {}
        for comm_id,comm in self.items():
            w_abunds = {}
            for x in comm.taxa.index:
                w_abunds[x] = comm.taxa[x] * self.genome_sizes[x]
            sum_abunds = sum(w_abunds.values())
            w_abunds = {k : v / sum_abunds * 100.0 for k,v in w_abunds.items()}
            self.comms[comm_id].taxa = pd.Series(w_abunds)
            
    def _set_genome_sizes(self):
        """Get the total bp for the genome sequence of all genomes
        
        Parameters
        ----------
        """
        self.genome_sizes = {}
        for taxon,fasta in self.genome_fasta.items():
            self.genome_sizes[taxon] = 0
            for record in SeqIO.parse(fasta, "fasta"):
                self.genome_sizes[taxon] += len(record.seq)        

    def beta_diversity(self, measures=['braycurtis'], outfile=None):
        """Calculate weighted abundance: rel_abund * genome_size
        
        Parameters
        ----------
        measure : list
          which scipy.spatial.distance function to use
        outfile : str
          if not None, then file written to that path
        """
        # taxon abudance table
        df =  pd.concat([x.taxa for x in self.values()],
                        axis=1, sort=False)
        df.fillna(0, inplace=True)
        # distance 
        D = defaultdict(dict)
        psbl_msr = {'braycurtis' : distance.braycurtis,
                    'jaccard' : distance.jaccard,
                    'euclidean' : distance.euclidean}
        for c1,c2 in combinations(df.columns, 2):
            for m,func in psbl_msr.items():
                if m in measures:
                    if m == 'jaccard':
                        D[m][(c1,c2)] = func(df[c1].astype(bool),
                                             df[c2].astype(bool))
                    else:
                        D[m][(c1,c2)] = func(df[c1], df[c2])
        if outfile is not None:
            with open(outfile, 'w') as outF:
                h = ['commX', 'commY', 'measure', 'distance']
                outF.write('\t'.join(h) + '\n')
                for m,v in D.items():                    
                    for c,d in v.items():
                        d = str(round(d,6))
                        c = [str(x) for x in list(c)]
                        outF.write('\t'.join(c + [m, d]) + '\n')
        return D                          
                
    @staticmethod
    def permute(comm, perm_perc):
        """Permute a certain percentage of the taxa abundances.
        Permuting just the indices of the series objects.
        In-place edit of comm table

        Parameters
        ----------
        comm : comm table object
        perm_perc : float
            percent of taxa to permute
        """
        # assertions
        perm_perc = float(perm_perc)
        assert (perm_perc >= 0 and perm_perc <= 100),\
            'perm_perc is not in range [0,100]'
        assert hasattr(comm, 'taxa'), \
            'No "taxa" attribute for comm {}'.format(comm.comm_id)

        # variables
        n_perm = int(round(perm_perc / 100 * comm.n_taxa,0))
        
        # permuting index of comm
        perm_idx = random.sample(range(comm.n_taxa), n_perm)
        perm_ig = itemgetter(perm_idx)
        n_perm_idx = set(range(comm.n_taxa)) - set(perm_idx)
                
        if len(n_perm_idx) > 0:
            n_perm_ig = itemgetter(*n_perm_idx)        
            # altering pandas series of taxa & abundances
            comm.taxa.index = random_insert_seq(n_perm_ig(comm.taxa.index),
                                                perm_ig(comm.taxa.index))
        else:
            # altering pandas series of taxa & abundances
            comm.taxa.index = random_insert_seq([],
                                                perm_ig(comm.taxa.index))

    # dict functions
    def items(self):
        return self.comms.items()
    def keys(self):
        try:
            return self.comms.keys()
        except AttributeError:
            return np.sort(list(self.comm_params.keys()))
    def values(self):
        return self.comms.values()
    
    @property
    def n_comm(self):        
        return self._n_comm
    @n_comm.setter
    def n_comm(self, x):
        try:
            self._n_comm = int(x)
        except ValueError:
            raise ValueError('n_comm must be an integer')

    @property
    def perm_perc(self):
        return self._perm_perc
    @perm_perc.setter
    def perm_perc(self, x):
        x = float(x)
        assert (x >= 0 and x <= 100), 'shared_perc must be in range 0-100'
        self._perm_perc = x        

    @property
    def shared_perc(self):
        return self._shared_perc
    @shared_perc.setter
    def shared_perc(self, x):
        x = float(x)
        assert (x >= 0 and x <= 100), 'shared_perc must be in range 0-100'
        self._shared_perc = x

    @property
    def min_richness(self):
        """The minimum richness of any community as defined by comm_params.
        """
        if not hasattr(self, '_min_richness'):
            setattr(self, '_min_richness', None)
        if self._min_richness is None:
            richness_vals = []
            for k,v in self.comm_params.items():
                try:
                    richness_vals.append(int(v['richness']))
                except KeyError:
                    raise KeyError('Cannot find richness attribute for '\
                                   'comm_id "{}"'.format(k))
            self._min_richness = min(richness_vals)
        return self._min_richness

    @property
    def n_shared(self):
        """The number of taxa that should be shared;
        defined by shared_perc * min richness of any community.
        """
        if not hasattr(self, '_n_shared'):
            setattr(self, '_n_shared', None)
        if self._n_shared is None:
            self._n_shared = self.min_richness * (self.shared_perc / 100.0)
            self._n_shared = int(round(self._n_shared,0))
        return self._n_shared

    @property
    def n_taxa_remaining(self):
        """The number of taxa that remain in taxa pool.
        """
        if not hasattr(self, '_n_taxa_remaining'):
            setattr(self, '_n_taxa_remaining', None)
        return len(self.taxon_pool)
                    
    
class Comm(_Comm):
    """Community class"""

    def __init__(self, comm_id, Comms, *args, **kwargs): 
        """
        Parameters
        ----------
        comm_id : str
            community ID
        Comms : community object
        """
        _Comm.__init__(self, *args, **kwargs)
        self.comm_id = comm_id                                                     
        self.params = Comms.comm_params[comm_id]
        self.n_shared = Comms.n_shared
        self.taxon_pool = Comms.taxon_pool
        self.richness = self.params['richness']

        # assertions
        if self.richness > self.n_shared + Comms.n_taxa_remaining:
            sys.exit('ERROR: Comm_ID {}\n'\
            '  Community richness is set too high! It is > taxon pool.\n'\
            '  There is not enough taxa to make the desired communities.\n' \
            '  You must reduce richness or increase perc_shared.\n'\
            '  NOTE: shared_perc is based on the community with the min. ' \
            'richness.\n'.format(comm_id))
        
        # selecting additional taxa beyond those shared by all comms
        ## unique taxa inserted rand in list while keeping shared taxa r-abund
        n_unique = self.richness - Comms.n_shared
        assert n_unique >= 0, 'ERROR: Comm_ID {}: the number ' \
            'of unique taxa is < 0'.format(comm_id)
        self.taxa = random_insert_seq(Comms.shared_taxa,
                                      Comms._drawFromTaxonPool(n_unique))

        
        # drawing relative abundances from the user-defined distribution
        abund_dist = self._get_abund_dist(self.params['abund_dist'],
                                          self.params['abund_dist_p'])        
        rel_abunds = abund_dist(size=self.n_taxa)
        rel_abunds = np.sort(rel_abunds / sum(rel_abunds) * 100)[::-1]
    
        # making a series for the taxa
        self.taxa = pd.Series(rel_abunds, index=self.taxa)        

    def __repr__(self):
        return self.taxa.__repr__()
        
    def _get_abund_dist(self, dist, params):
        try:
            distFunc = getattr(np.random, dist)
        except AttributeError:
            msg = 'Distribution "{}" is not supported'.format(dist)

        if dist == 'power':
            distFunc = power_neg
                        
        try:
            return partial(distFunc, **params)
        except TypeError:
            param_str = [':'.join([str(k),str(v)]) for k,v in params.items()]
            param_str = ','.join(param_str)
            msg = 'Params "{}" do not work with function "{}"'\
                  .format(param_str, dist)
            raise TypeError(msg)

    @property
    def n_taxa(self):
        return len(self.taxa)
