"""Classes and subclasses for genome fragment simulation"""

# import
## batteries
import os
import sys
import random
import functools

## 3rd party
import numpy as np
import skew_normal
import scipy.stats as stats
from rtnorm import rtnorm



def truncated_normal(location, scale, low, high, size=1):
    """Creating a function to randomly sample from a truncated normal distribution.
    Function returned from stats.logistic.ppf

    Parameters
    ----------
    loc : float
        location parameter of distribution
    scale : float
        standard deviation parameters of distribution
    low, high : float
        min|max value that can be sampled from full distribution
    size : int
        number of random samples to return

    Returns
    -------
    function : distribution function
    """
    nrm = stats.logistic.cdf(high)-stats.logistic.cdf(low)
    yr = np.random.rand(size)*(nrm)+stats.logistic.cdf(low)
    return stats.logistic.ppf(yr, loc=location)

    
class SimFrags(object):
    """Class for genome fragment simulation."""

    def __init__(self, fld, flr, rtl=None):
        """
        Paramters
        ---------
        fld : function
            fragment length distribution (type,moments...)
        flr : tuple
            fragment length range (min,max)
        rtl : tuple
            read template length distribution (type,moments..)
        """        
        self.flr = [float('inf') if
                    x.lower() == 'none' or
                    x.lower == 'inf' or
                    x is None else
                    int(x) for x in flr]
        if self.flr[0] == float('inf'):
            self.flr[0] = 1
        
        # setting fld distribution function
        self.fld = self._set_user_dist(fld)
        msg = 'fragment length distribution (fld) should not be None'
        assert self.fld is not None, msg
        
        # setting rtl distribution function
        if rtl is not None:
            self.rtl = self._set_user_dist(rtl)

        
    def simFrag(self, genome):
        """Simulation of a single fragment.

        Parameters
        ----------
        genome : genome object

        Returns
        -------
        list : [scaffold,fragStart,fragLen,fragSequence]
        """        

        tryCnt = 0
        while 1:
            # get read template position
            if genome.primerFile is not None:  # amplicon
                readTempPos = self._get_randPos_amp(genome)
            else:   # shotgun
                readTempPos = self._get_randPos_shotgun(genome)
            #msg = 'readTempPos should be: "scaf,start,end"'
            #assert len(readTempPos) == 3, msg
            

            # if no read template position (no amplicon): return None
            if readTempPos[0] is None:
                return ["No_amplicons", 'NA', 'NA', '']
            else:
                scafID = readTempPos[0]

            # read template sequence size
            if genome is not None:
                tempSize = genome.get_seq_len(scafID)
            else:
                tempSize = np.inf
                
            # simulating fragment around read template
            fragPos = self._simFragPos(readTempPos, tempSize)                        
            
            # parsing fragment from genome index
            fragSeq = genome.get_seq(*fragPos)

            # getting fragment info
            fragGC = genome.calcGC(fragSeq)
            fragSeqLen = len(fragSeq)            
            
            # checking to see if fragment acutally in range (accounting for edge cases)
            if tryCnt >= 1000:
                msg = 'Exceeded {} tries to find frag'\
                      ' length in min-max range.' \
                      ' You may need to adjust --flr.' 
                raise ValueError(msg.format(tryCnt)) 
            elif (fragSeqLen >= self.minFragSize
                and fragSeqLen <= self.maxFragSize):
                break
            else:
                tryCnt += 1
                continue
                
        return [fragPos[0], fragPos[1], fragSeqLen, fragGC]

        
    def _simFragPos(self, readTempPos, tempSize):
        """Simulate the position of the fragment containing the read template.
        Location selected using a uniform distribution.
        Fragment length selected using user-provide distribution
        (set during initilization).
        Multiple iterations will be performed to try and get 
        a fragment in range (minFragSize, maxFragSize).

        Parameters
        ----------
        readTempPos : list
            [scaffoldID,readTempStart,readTempEnd]                       

        Returns
        -------
        list : [scafID,fragStart,fragEnd]
        """
        # assertions
        if readTempPos[0] is None:
            return [None] * 3

        # read template size 
        ## (just size of read; not full template sequence)
        readTempSize = readTempPos[2] - readTempPos[1] + 1
                
        # nfrag size
        tryCnt = 0
        while 1:
            if tryCnt >= 1000:
                msg = 'Exceeded {} tries to find frag'\
                      ' length in min-max range.' \
                      ' You may need to adjust --flr.' 
                raise ValueError(msg.format(tryCnt))

            # draw from fragment size distribution
            fragSize = int(self.fld(size=1))
            if fragSize <= readTempSize:
                tryCnt += 1
                continue

            # frag start (position upstream from read template)        
            randPos = np.random.randint(0, fragSize - readTempSize)
            fragStart = readTempPos[1] - randPos
            fragEnd = fragStart + fragSize - 1
            # constraining start + end to within template fragment size
            if fragStart < 0:
                fragStart = 0
            if fragEnd + 1 > tempSize:
                fragEnd = fragStart + tempSize - 1

            fragLen = fragEnd - fragStart + 1

            if fragSize > readTempSize and \
               fragLen >= self.minFragSize and \
               fragLen <= self.maxFragSize:
                break
            else:
                tryCnt += 1
                continue
                
        # frag start should be < frag end
        if fragStart >= fragEnd:
            raise TypeError('fragStart >= fragEnd')
        
        return readTempPos[0], fragStart, fragEnd
                

    def _set_user_dist(self, userDist):
        """Setting user defined distribution. Using numpy distribution functions.
        
        Parameters
        ----------
        userDist : tuple
            User defined distribution with moment info.
            Example: ('normal',10,1)
        Returns
        -------
        function : partial numpy distribution function with moment values provided
        """
        userDist[0] = str(userDist[0]).lower()
        userDist[1:] = [int(x) for x in userDist[1:]]
        
        if userDist[0] == 'normal':
            assert len(userDist) >= 3, ('loc and scale must be provided'
                                        'for "normal" distribution')
            return functools.partial(np.random.normal,
                                     loc=userDist[1],
                                     scale=userDist[2])
        elif userDist[0] == 'uniform':
            assert len(userDist) >= 3, ('low and high must be provided'
                                        'for "uniform" distribution')
            return functools.partial(np.random.uniform,
                                     low=userDist[1],
                                     high=userDist[2])
        elif userDist[0] == 'skewed-normal':
            assert len(userDist) >= 4, ('loc, scale, and shape must be provided'
                                        'for "skewed-normal" distribution')
            return functools.partial(skew_normal.rnd_skewnormal,
                                     location=userDist[1],
                                     scale=userDist[2],
                                     shape=userDist[3])
        elif userDist[0] == 'truncated-normal':
            assert len(userDist) >= 5, ('loc, scale, low, and high must be provided'
                                       'for "truncated-normal" distribution')
            return functools.partial(rtnorm,
                                     mu=userDist[1],
                                     sigma=userDist[2],
                                     a=userDist[3],
                                     b=userDist[4])
        else:
            msg = 'Distribution "{}" is not supported.'
            raise ValueError(msg.format(userDist[0]))
                        
        
    def _get_randPos_amp(self, genome):
        """Getting the genomic position of a randomly selected amplicon.
        The amplicon is selected randomly from the list of amplicons.
        None values are returned if no amplicons exist.
        
        Parameters
        ----------
        genome : genome object

        Returns
        -------
        list : [read_scaffold_ID,start,end] 
        """
        nAmps = genome.nAmplicons        
        
        if nAmps <= 0:
            return [None] * 3
        elif nAmps == 1:
            i = 0
        else:
            i = np.random.randint(0, nAmps-1)
        row = genome.MFEprimerRes.iloc[i]
        row = row[['HitID','BindingStart','BindingStop']]        
        return [x for x in row]
        

    def _get_randPos_shotgun(self, genome):
        """Randomly selecting start-end position of read template.
        The scaffold/chromosome is randomly selected.
        start selected from a uniform distribution.
        end seleded based on read template length distribution.
        
        Parameters
        ----------
        genome : genomeobject

        Returns
        -------
        list : [scafName,start,end]
               scaffold and start-end of read template
        """
        assert hasattr(self, 'rtl'), '"rtl" attribute required'
        
        # randomly selecting a genome scaffold
        scafName = random.choice(genome.fastaIdx.keys())
        scaf = genome.fastaIdx[scafName]

        # start
        start = np.random.randint(scaf.start, scaf.stop)
        #pymcDist.rdiscrete_uniform(scaf.start, scaf.stop)
        
        # end        
        mpSize = self.rtl(size=1)
        assert mpSize > 0, 'Read template size must be > 0'

        # scaffoldName, readTempStart, readTempEnd
        return scafName, start, start + int(mpSize) -1

        
    @property
    def minFragSize(self):
        return self.flr[0]

    @property
    def maxFragSize(self):
        return self.flr[1]
                                    
