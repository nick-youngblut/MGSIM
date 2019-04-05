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
from MGSIM.Commands import Genome_download as Genome_download_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    acc_tbl = os.path.join(data_dir, 'genome_download.txt')
    args = ['-d', data_dir, acc_tbl, '--debug']
    Genome_download_CMD.opt_parse(args)

def test_cmd_rename():
    acc_tbl = os.path.join(data_dir, 'genome_download.txt')
    args = ['-d', data_dir, acc_tbl, '--debug', '-r']
    Genome_download_CMD.opt_parse(args)

def test_cmd_badAcc():
    acc_tbl = os.path.join(data_dir, 'genome_download_badAcc.txt')
    args = ['-d', data_dir, acc_tbl, '--debug', '-r']
    Genome_download_CMD.opt_parse(args)

def test_cmd_dupTaxa():
    acc_tbl = os.path.join(data_dir, 'genome_download_dupTaxa.txt')
    args = ['-d', data_dir, acc_tbl, '--debug', '-r']
    with pytest.raises(ValueError):
        Genome_download_CMD.opt_parse(args)
