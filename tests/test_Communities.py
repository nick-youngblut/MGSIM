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
from MGSIM.Commands import Communities as Communities_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_help(script_runner):
    ret = script_runner.run('MGSIM', 'communities', '-h')
    assert ret.success

def test_main(script_runner):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(data_dir, 'comm')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '2', '--rnd-seed', '34847')
    assert ret.success

def test_main_ncomm1(script_runner):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(data_dir, 'comm')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '1', '--richness', '1')
    assert ret.success

