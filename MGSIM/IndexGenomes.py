#!/usr/bin/env python

# import
## batteries
import sys,os
import platform
from multiprocessing import Pool
## 3rd party
import functools
## application libraries
from Genome import Genome
import Utils

    
def index_genome(x, faToTwoBitExe, fmtExe, idxExe, K_value, quiet=False):
    """index a genome with MFEprimer indexing scripts. 
    This function is just making system calls.
    
    Parameters
    ----------
    x : tuple
        (taxonName, genomeFile)
        taxonName -- taxon name of genome fasta
        genomeFile -- file name of genome fasta
    faToTwoBitExe : str
        path of faToTwoBit executable
    fmtExe : str
        path of UniFastaFormat.py executable
    idxExe : str
        path of mfe_index_db.py executable
    K_value : int
        k value used for indexing
    quiet : bool
        quiet all messages
    """
    taxonName,genomeFile = x
    # status
    if not quiet:
        sys.stderr.write('Indexing: "{}"\n'.format(taxonName))
        
    # begin indexing
    cmd = '{} -i {}'.format(fmtExe, genomeFile)
    Utils.sys_call(cmd)    
    
    # faToTwoBit
    cmd = '{} {} {}'.format(faToTwoBitExe,
                            genomeFile + '.unifasta',
                            genomeFile + '.2bit')
    Utils.sys_call(cmd)
    
    # index
    cmd = '{} -f {} -k {}'.format(idxExe,
                                  genomeFile + '.unifasta',
                                  K_value)
    Utils.sys_call(cmd)

    # cleanup
    os.remove(genomeFile + '.unifasta')


def main(args):
    """
    Parameters
    ----------
    args : dict
        See ``genome_index`` subcommand
    """    
    # loading genome list
    genomeList = Utils.parseGenomeList(args['<genomeList>'], filePath=args['--fp'])
                         
    # setting function for parallel calling
    pfunc = functools.partial(index_genome,                              
                              faToTwoBitExe=args['--twobit'],
                              fmtExe=args['--fmt'],
                              idxExe=args['--idx'],
                              K_value=args['--K_value'],
                              quiet=args['--quiet'])                                         


    # indexing genomes in parallel
    pool = Pool(int(args['--np']))
    if args['--debug']:
        KDE_BD = map(pfunc, genomeList)
    else:
        KDE_BD = pool.map(pfunc, genomeList)

    # status
    sys.stderr.write('#-- All genomes indexed --#\n')

        
