"""Fragment KDE classes"""

# import
## batteries
import sys,os
import math
import logging
import functools
import cPickle as pickle
from multiprocessing import Pool
## 3rd party
import numpy as np
import pandas as pd
import scipy.stats as stats
## Application
from Genome import Genome
from SimFrags import SimFrags
import Utils

# logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)



def main(args):
    """
    Parmeters
    ---------
    args : dict
        See ``fragments`` subcommand
    """
    # list of genome files
    genomeList =  Utils.parseGenomeList(args['<genomeList>'], 
                                        filePath=args['--fp'])
        
    # analyzing each genome (in parallel)    
    pfunc = functools.partial(by_genome, args=args)
    
    # difussion calc in parallel
    pool = Pool(int(args['--np']))
    if args['--debug']:
        fragList = map(pfunc, genomeList)
    else:
        fragList = pool.map(pfunc, genomeList)

    # writing out table
    if args['--tbl']:
        write_fragList(fragList)
    else:
        dill.dump(fragList, sys.stdout)


def load_frags_table(inFH, sep='\t'):
    """Load frag info table as a dict of dicts of 2d lists.
    
    Parameters
    ----------
    inFH : file handle
    sep : str
        value delimiter
    
    Returns
    -------
    frags : dict
        {taxon_name : {scaffold : [fragStart, fragEnd, GC]}}
    """    
    header_vals = set(['taxon_name','scaffoldID','fragStart',
                       'fragLength','fragGC'])
    
    d = dict()
    lineNum = 0
    for line in inFH.readlines():
        lineNum += 1
        line = line.rstrip().split(sep)

        #header
        if lineNum == 1:            
            if not (header_vals == set(line) or header_vals < set(line)):
                msg = 'The fragGC table does not have all'\
                      ' required columns:\n\t{}'\
                      .format(','.join(header_vals))
                raise IOError(msg)
            header_idx = {line[i]:i for i in xrange(len(line))}
        # body            
        else:
            taxon_name = line[header_idx['taxon_name']]
            try:
                type(d[taxon_name])
            except KeyError:
                d[taxon_name] = dict()
                d[taxon_name]['fragLength'] = []
                d[taxon_name]['fragGC'] = []

            # fragment length
            fragLength = line[header_idx['fragLength']]
            try:
                fragLength = int(fragLength)
            except ValueError:
                continue

            # fragment GC content
            fragGC = line[header_idx['fragGC']]
            try:
                fragGC = float(fragGC)
            except ValueError:
                continue

            # adding to dict
            d[taxon_name]['fragLength'].append(fragLength)
            d[taxon_name]['fragGC'].append(fragGC)

    return d

            
def load_frags_pickle(inFH):
    """Load frag GC info assuming a pickled python object
    produced by SIPSim fragGC.
    
    Parameters
    ----------
    inFH : file handle

    Returns
    -------
    d : dict
        {taxon_name : {info : [values]}}
    """
    fojb =  pickle.load(inFH)

    d = dict()
    for x in fojb:
        taxon_name = x[0]
        if d.has_key(taxon_name):
            msg =  'WARNING: {} was found multiple times\n'
            sys.stderr.write(msg.format(taxon_name))
        d[taxon_name] = dict()
        d[taxon_name]['fragLength'] = []
        d[taxon_name]['fragGC'] = []
            
        for scaf,v in x[1].items():            
            for z in v:
                # fragStart, fragLength, fragGC
                d[taxon_name]['fragLength'].append(z[1])
                d[taxon_name]['fragGC'].append(z[2])              
    return d


def load_frags(fileName):
    """Load fragment data (pickled) table.
    
    Parameters
    ----------
    fileName : str
        name of the fragment data table

    Returns
    -------
    d : dict{dict} 
        {taxon_name:{key:value}}
    """
    try:
        inFH = open(fileName, 'r')
    except IOError:
        inFH = sys.stdin

    try:
        frag_data = load_frags_pickle(inFH)
    except (pickle.UnpicklingError, EOFError):
        try:
            inFH.seek(0)
        except IOError:
            msg = ('Illegal seek; either you piped in a non-pickled table or'
                   'your file name is incorrect')
            raise IOError(msg)
        frag_data = load_frags_table(inFH)            

    inFH.close()
    
    return frag_data
     
def by_genome(x, args):
    """All processing conducted per genome.

    Parameters
    ----------
    x : list
        [inFile,taxonName]
        inFile -- genome sequence file name
        taxonName -- taxon name of genome
    args : dict
       user-provided args 

    Returns
    -------
    l2d -- list of lists
        for each fragment: [taxonName,scaf,start,end,GC]
    """
    taxonName,inFile = x
    # status
    sys.stderr.write('Processing: "{}"\n'.format(taxonName))

    # making genome object
    assert '--fr' in args, '"--fr" must be provided in args'
    genome = Genome(inFile, taxonName, args['--fr'])
    
    # MFEprimer.py executable
    MFEprimerExe = args['--MFE']
    
    # sequenced read template location: amplicons
    if genome.primerFile is not None:
        # in-silico PCR
        assert '--rtr' in args, '"--rtr" must be in args'
        genome.callMFEprimer(rtr=args['--rtr'], MFEprimerExe=MFEprimerExe)
    
        # filtering overlapping in-silico amplicons
        genome.filterOverlaps()
                
    # simulating fragments    
    simFO = SimFrags(fld=args['--fld'], flr=args['--flr'], rtl=args['--rtl'])
    nFragsMade = 0
    fragList = dict()
    ## if no amplicons
    if genome.nAmplicons == 0:
        pass
    ## if using coverage
    elif args['--nf'].endswith('X') or args['--nf'].endswith('x'):
        coverage = float(args['--nf'].rstrip('xX'))
        fragLenCov = genome.length * coverage
        fragLenTotal = 0
        while 1:
            (scaf,fragStart,fragLen,fragGC) = simFO.simFrag(genome)
            try:
                type(fragList[scaf])
            except KeyError:
                fragList[scaf] = []
                                
            if fragStart == "NA":
                break
            elif fragLenTotal > fragLenCov:
                break
            fragLenTotal += fragLen 

            nFragsMade += 1
            fragList[scaf].append([fragStart, fragLen, fragGC])            
    ## if using fixed number of fragments
    else:            
        for i in xrange(int(args['--nf'])):
            (scaf,fragStart,fragLen,fragGC) = simFO.simFrag(genome)

            try:
                type(fragList[scaf])
            except KeyError:
                fragList[scaf] = []

            if fragStart == "NA":
                break

            nFragsMade += 1
            fragList[scaf].append([fragStart, fragLen, fragGC])
                
    # status
    sys.stderr.write('  Genome name: {}\n'.format(genome.taxonName))                
    sys.stderr.write('  Genome length (bp): {}\n'.format(genome.length))
    if args['--nf']:
        msg = '  Number of amplicons: {}\n'
        sys.stderr.write(msg.format(genome.nAmplicons))
    msg = '  Number of fragments simulated: {}\n'
    sys.stderr.write(msg.format(nFragsMade))
                
    return [genome.taxonName, fragList]


def write_fragList(fragList):
    """Write out fragList as a tab-delim table."""
    print '\t'.join(['taxon_name','scaffoldID','fragStart',
                     'fragLength','fragGC'])            
    for x in fragList:
        taxon_name = x[0]
        for scaf,v in x[1].items():
            for y in v:                
                print '\t'.join([taxon_name, scaf] + [str(i) for i in y])


        

