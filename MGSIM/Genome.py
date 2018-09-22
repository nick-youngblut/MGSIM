"""Class for simulating fragments from genomes"""

# import
## batteries
import sys, os
import subprocess
from StringIO import StringIO
import itertools

## 3rd party
from Bio import SeqIO
import pandas as pd
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
        #self.fastaIdx = pyfasta.Fasta(self.get_fileName())
        
        
                                
    # getters/setters/iters
    def get_fileName(self, rmPath=False):
        """Get fileName, without path if rmPath=True."""
        if rmPath is True:
            return os.path.split(self.fileName)[1]
        else:
            return self.fileName    
                                

