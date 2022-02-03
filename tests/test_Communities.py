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

def chk_suc(ret):
    # checking cmd return
    if not ret.success:
        ret.print()
    assert ret.success

# tests
def test_help(script_runner):
    ret = script_runner.run('MGSIM', 'communities', '-h')
    chk_suc(ret)

def test_main(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(str(tmp_path), 'comm')    
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '2', '--rnd-seed', '34847')
    chk_suc(ret)

def test_main_ncomm1(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(str(tmp_path), 'comm')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '1', '--richness', '1')
    chk_suc(ret)

def test_main_shared(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list_n12.txt')
    output_prefix = os.path.join(str(tmp_path), 'comm-n12')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '2',
                            '--richness', '1',
                            '--shared-perc', '0.8')
    chk_suc(ret)
    
def test_main_ncomm_bias1(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list_n12.txt')
    output_prefix = os.path.join(str(tmp_path), 'comm-n12')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '2',
                            '--richness', '1',
                            '--group-bias', '0.5')
    chk_suc(ret)

def test_main_ncomm_bias2(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list_n12.txt')
    output_prefix = os.path.join(str(tmp_path), 'comm-n12')
    ret = script_runner.run('MGSIM', 'communities',
                            genomeList, output_prefix,
                            '--n-comm', '2',
                            '--richness', '1',
                            '--group-bias', '-0.7')
    chk_suc(ret)
