#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import pytest
## 3rd party
import numpy as np
import pandas as pd
## package
from MGSIM import Utils
from MGSIM.Commands import Reads as Reads_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('MGSIM', 'reads', '-h')
    assert ret.success

def test_main(script_runner):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(data_dir, 'temp_read_files')
    output_prefix = os.path.join(data_dir, 'TEST')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    assert ret.success

def test_main_multiProc(script_runner):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(data_dir, 'temp_read_files')
    output_prefix = os.path.join(data_dir, 'TEST')
    ret = script_runner.run('MGSIM', 'reads',
                            '-n', '2',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    assert ret.success

def test_main_unpaired(script_runner):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(data_dir, 'temp_read_files')
    output_prefix = os.path.join(data_dir, 'TEST')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-mflen', '0',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    assert ret.success

