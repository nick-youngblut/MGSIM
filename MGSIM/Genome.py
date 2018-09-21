"""Class for simulating fragments from genomes"""

# import
## batteries
import sys, os
import subprocess
from StringIO import StringIO
import itertools

## 3rd party
import pyfasta
#import pydna
from Bio import SeqIO
import pandas as pd
from intervaltree import Interval, IntervalTree
import numpy as np

## appliation
import Utils
from SimFrags import SimFrags


            
class Genome(object):
    """subclass of genomes: 1 genome entry"""

    def __init__(self, inFile, taxonName, primerFile=None):
        """
        Parameters
        ----------
        inFile : str
            name of the genome sequence file
        taxonName : str
            name of taxon
        primerFile : str
            file name of primers in fasta format
        """
        self.fileName = inFile
        self.taxonName = taxonName
        self.primerFile = primerFile
        self._nAmplicons = None
        self._length = None
        
        # checking that the genome files exist & not empty
        Utils.checkExists(self.fileName)
        Utils.checkEmpty(self.fileName)
        
        # create fasta index
        self.fastaIdx = pyfasta.Fasta(self.get_fileName())
        
        
    def callMFEprimer(self, rtr, MFEprimerExe='./MFEprimer.py'):
        """Cal MFEprimer to get amplicons.
        in-place edit of MFEprimerRes object

        Parameters
        ----------
        rtr : tuple
            read template size range (min,max)
        MFEprimerExe : str
            path of MFEprimer.py
        """        
        # Input check
        try:
            (minTemplateLen,maxTemplateLen) = rtr[:3]
        except IndexError:
            raise IndexError('"rtr" must contain (min,max)')

        # system call
        cmd = '{} -i {} -d {} --tab --size_start {} ' \
              '--size_stop {} --ppc 70 --tm_start 50 --tm_stop 70'
        cmd = cmd.format(MFEprimerExe,
                         self.primerFile,
                         self.get_fileName(),
                         minTemplateLen,
                         maxTemplateLen)
        ret = subprocess.check_output(cmd, shell=True)
        
        # load results as dataframe
        self.MFEprimerRes = pd.read_csv(StringIO(ret), sep='\t')
        
        # editing values
        try:
            self.MFEprimerRes.BindingStart.astype(int)
        except AttributeError:
            msg = 'MFEprimerRes attribute error. Was genome indexing successful?'
            raise AttributeError(msg)
        self.MFEprimerRes.BindingStop.astype(int)        

        # checking values
        self._MFEprimerRes_chk()


    def _MFEprimerRes_chk(self):
        """Assertions on output table from MFEprimer."""
        if self.MFEprimerRes is None:
            return None
        
        # binding start/end < 0?
        for myCol in ['BindingStart', 'BindingStop']:            
            x = self.MFEprimerRes.loc[self.MFEprimerRes[myCol] < 0].shape
            if x[0] > 0:
                msg = myCol + " < 1"
                raise TypeError, msg
        
        
    def filterOverlaps(self, overlapPercCutoff=70):
        """Filtering out amplicons that substantially overlap.
        The amplicon with the highest PPC with be kept.
        The MFEprimerRes attribute must be set.
        in-place edit of MFEprimerRes object (table object filtered of overlaps)
        
        Parmeters
        ---------
        overlapPercCutoff : float
            percent of overlap to consider 'substantially' overlapping
        """
        if self.MFEprimerRes is None:  
            msg = 'genome object does not have MFEprimerRes attribute.' + \
                  ' Run MFEprimer() first'
            raise AttributeError, msg
            
        # making interval tree
        tree = IntervalTree()
        
        # loading intervals
        for count,row in self.MFEprimerRes.iterrows():
            # sanity check for + strand
            if row['BindingStart'] > row['BindingStop']:
                raise TypeError('MFEprimer binding start-stop is not + strand')
            tree.addi(row['BindingStart'], row['BindingStop'],
                      [count, row['PPC'], row['Size']])
                        
        # finding all that substantially overlap; keeping one with > PPC
        tree2 = tree.copy()
        for iv1 in tree.iter():
            # skipping if already removed from tree2
            if not iv1 in tree2:
                continue
                
            overlaps = tree.search(iv1.begin, iv1.end)

            # skipping those that poorly overlap 
            lowOverlap = set()
            for iv2 in overlaps:
                if iv1.overlaps(iv2):
                    percOverlaps = self._calcPercOverlap(iv1, iv2)
                    if percOverlaps[0] < overlapPercCutoff:
                        lowOverlap.add(iv2)
            overlaps = overlaps - lowOverlap  # just list of substantially overlapping

            # skipping those that have been already removed
            prevRm = set([x for x in overlaps if x not in tree2])
            overlaps = overlaps - prevRm
            
            # removing all substantially overlapping intervals with lower PPC
            if len(overlaps) > 1:
                overlaps = sorted(overlaps, 
                                  key=lambda x: x.data[1], 
                                  reverse=True)
                for o in overlaps[1:]:
                    if o in tree2:
                        tree2.remove(o)
            else:
                pass
 
        # selecting columns
        iv_idx = [x.data[0] for x in tree2.iter()]
        self.MFEprimerRes = self.MFEprimerRes.iloc[iv_idx]                
                      
        
    @staticmethod
    def _calcPercOverlap(iv1, iv2):
        """Calculate the overlap between intervals (iv1, iv2).
        
        Parmeters
        ---------
        iv1, iv2 : interval

        Returns
        -------
        tuple : (percent overlap relative to iv1, relative to iv2)
        """
        if not iv1.overlaps(iv2):
            return [0.0, 0.0]
        elif iv1 == iv2:
            return [100.0, 100.0]
        else:
            overlapSize = len(set(range(iv1.begin,iv1.end +1)) &
                              set(range(iv2.begin,iv2.end +2)))
            overlapSize = float(overlapSize)
            iv1PercOverlap = (iv1.end - iv1.begin + 1) / overlapSize * 100.0
            iv2PercOverlap = (iv2.end - iv2.begin + 1) / overlapSize * 100.0
            return [iv1PercOverlap, iv2PercOverlap]


    @staticmethod
    def calcGC(seq):
        """Calculating GC content of a sequence string. Return as percent G+C.
        Only unambiguous nucleotide codes (ATGC) will be counted.
        
        Parameters
        ----------
        seq : str
            DNA sequence

        Returns
        -------
        float : %GC
        """
        seq = str(seq).upper()
        aCount = seq.count('A')
        tCount = seq.count('T')        
        cCount = seq.count('C')
        gCount = seq.count('G')

        try:
            b = float(aCount + tCount + cCount + gCount)
            return float(cCount + gCount) / b * 100.0
        except ZeroDivisionError:
            return 'NA'
            
                    
    # getters/setters/iters
    def get_fileName(self, rmPath=False):
        """Get fileName, without path if rmPath=True."""
        if rmPath is True:
            return os.path.split(self.fileName)[1]
        else:
            return self.fileName    
                                
    def get_seq(self, scaffold, start, end):
        """Getting sequence from genome. 0-indexing.
        start-end will return at section of the sequence.
        
        Parameters
        ----------
        scaffold : str
            scaffold id
        start, end : int
            sequence start-end positions

        Returns
        -------
        str -- sequence
        """
        try:
            return str(self.fastaIdx[scaffold][start:end+1])
        except AttributeError:
            raise AttributeError('No fastaIdx attribute for genome object')
        except KeyError:
            msg = 'ScaffoldID "{}" not found in genome fasta index'
            raise KeyError(msg.format(scaffold))


    def get_seq_len(self, seqID):
        """Get length for sequence 'seqID'.

        Parameters
        ----------
        seqID : str
            ID of sequence

        Returns
        -------
        int : sequence length
        """
        if not hasattr(self, '_fastaIdx_lens'):
            self._fastaIdx_lens = {seqID:len(seq) for seqID,seq 
                                   in self.fastaIdx.items()}
        
        try:
            return self._fastaIdx_lens[seqID]
        except KeyError:
            msg = 'Cannot find "{}" in sequences'
            raise KeyError, msg.format(seqID)
        

    def iter_seqRecs(self):
        """Iterator of sequence records.

        Returns
        -------
        iterator : each sequence record
        """
        with open(self.fileName) as inF:
            for rec in SeqIO.parse(inF, 'fasta'):
                yield rec


    @property
    def MFEprimerRes(self):
        """Get the results of MFEprimer. Returns a dataframe"""
        try:
            return self._MFEprimerRes
        except AttributeError:
            return None

    @MFEprimerRes.setter
    def MFEprimerRes(self, value):
        self._MFEprimerRes = value


    @property
    def primerFile(self):                            
        try:
            return self._primerFile
        except AttributeError:
            return None

    @primerFile.setter
    def primerFile(self, value):
        if (value is not None) and (not os.path.exists(value)):
            raise IOError, 'Primer file not found!'
        self._primerFile = value
            
    @property
    def length(self):
        if self._length is None:
            self._length = sum([len(x) for x in self.iter_seqRecs()])
        return self._length

    @property
    def nAmplicons(self):
        if self._nAmplicons is None:
            try:
                self._nAmplicons = self.MFEprimerRes.shape[0]
            except (ValueError, AttributeError):
                self._nAmplicons = None
        return self._nAmplicons
