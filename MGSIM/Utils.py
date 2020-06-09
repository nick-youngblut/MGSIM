"""Utility scripts for application"""

# import
## batteries
import os,sys
import re
import time
import logging
import platform
import subprocess
from pprint import pprint
from itertools import chain
from functools import partial
import random
import glob

## 3rd party
import numpy as np
import pandas as pd

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

def get_os():
    """ Get operating system; only works for unix-like machines
    """
    OS = platform.uname()[0]
    if OS == 'Linux':
        OS = 'linux'
    elif OS == 'Darwin':
        OS = 'mac'
    else:
        sys.stderr.write('OS: "{}" not supported\n'.format(OS))
    return OS

def is_file(fileName):
    """ Does file exist? with custom output message
    """
    if os.path.isfile(fileName) is False:
        raise IOError('"{}" does not exist'.format(fileName))
        
def sys_call(cmd, quiet=False):
    """System call of command.
    
    Parameters
    ----------
    cmd : str
        The command to run as a system call
    quiet : bool
        Suppress the system command output

    Returns
    -------
    output : str
       system call output
    err : str
       system call error
    """
    try:
        if quiet:
            DEVNULL = open(os.devnull, 'w')
            proc = subprocess.Popen([cmd], shell=True, stdout=DEVNULL)
        else:
            proc = subprocess.Popen([cmd], shell=True)
    except subprocess.CalledProcessError:
        pass # handle errors in the called executable
    except OSError:
        raise OSError('No executable for command: "{}"\n'.format(cmd))
    output, err = proc.communicate()
    return output, err

def checkExists(f):
    """ Check that the file `f` exists."""
    if not os.path.isfile(f):
        msg = '"{}" not found. Did you provide the full PATH?'
        raise IOError(msg.format(f))

def checkEmpty(f):
    """ Check that the file `f` is not empty"""
    if os.stat(f).st_size == 0:
        msg = '"{}" is empty!'
        raise IOError(msg.format(f))

def parseGenomeList(inFile, filePath=None, check_exists=True):
    """Parsing the genome list file.

    Parameters
    ----------
    inFile : str
        file name of genome list file
    filePath : str
        The absolute path to genome sequence files. a
    check_exists : bool
        Check if genome sequence files exist
    
    Returns:
    list : [[taxonName, genomeFile], ...]
    """
    # parse file as list
    genomeList = []
    with open(inFile, 'rb') as inF:
        for line in inF:
            row = line.rstrip().split('\t')
            
            if row[0] == '' or row[1] == '':
                raise IOError("Necessary row value is empty!")

            if len(row) < 2:
                raise ValueError('Need format: "taxonName<tab>fileName";'
                                 'for row: "{}"'.format(row))
            else:
                (taxonName,fileName) = row[:2]
                
            # path to genome file
            if filePath is not None:
                fileName = os.path.join(filePath, fileName)

            # checking for file existence
            if check_exists:
                checkExists(fileName)
                
            #genomeList[fileName] = taxonName
            genomeList.append((taxonName,fileName))
                
    return genomeList


def describe_builtin(obj):
    """ Describe a builtin function if obj.__doc__
    available.

    Parameters
    ----------
    obj : python object
    
    Returns
    -------
    iterator : builtin args
    """
    #wi('+Built-in Function: %s' % obj.__name__)
    # Built-in functions cannot be inspected by
    # inspect.getargspec. We have to try and parse
    # the __doc__ attribute of the function.
    docstr = obj.__doc__
    args = ''
    
    if docstr:
        items = docstr.split('\n')
        if items:
            func_descr = items[0]
            s = func_descr.replace(obj.__name__,'')
            idx1 = s.find('(')
            idx2 = s.find(')',idx1)
            if idx1 != -1 and idx2 != -1 and (idx2>idx1+1):
                args = s[idx1+1:idx2]
                #wi('\t-Method Arguments:', args)
                for arg in args:
                    yield arg
                
    if args=='':
        yield None


def parseKeyValueString(x):
    """Parse a string in format: 'key1:value1,key2:value2,keyN:valueN'.
    Values assumed to be numeric.
    
    Parameters
    ----------
    x : string
        Required format: 'key:value,key:value,...' or 'key=value,key=value,...'

    Returns
    -------
    dict : {key:value, ...}
    """
    if x is None or x == 'None':
        return {}
    x = x.replace(' ','')
    l = re.split('[=:,]', x)
    return {k.lower():float(v) for k,v in zip(l[::2],l[1::2])}


