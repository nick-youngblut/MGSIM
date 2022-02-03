#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import pytest
## 3rd party
import pandas as pd
## package
from MGSIM import Utils
from MGSIM.Commands import Genome_rename as Genome_rename_CMD

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
    ret = script_runner.run('MGSIM', 'genome_rename', '-h')
    chk_suc(ret)
    
def test_main(script_runner, tmp_path):
    fasta_files = ['Escherichia_coli_K-12_MG1655.fna',
                   'Clostridium_perfringens_ATCC_13124.fna']
    fasta_files = [os.path.join(data_dir, x) for x in fasta_files]
    prefix = os.path.join(str(tmp_path), 'renamed')
    ret = script_runner.run('MGSIM', 'genome_rename',
                            '--debug', '--prefix', prefix,
                            *fasta_files)
    chk_suc(ret)

def test_main_ambig(script_runner, tmp_path):
    fasta_files = ['Escherichia_coli_K-12_MG1655_ambig.fna']
    fasta_files = [os.path.join(data_dir, x) for x in fasta_files]
    prefix = os.path.join(str(tmp_path), 'renamed')
    ret = script_runner.run('MGSIM', 'genome_rename',
                            '--debug', '--prefix', prefix,
                            *fasta_files)
    assert not ret.success

