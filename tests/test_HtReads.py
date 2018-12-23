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
from MGSIM import SimHtReads
from MGSIM.Commands import HtReads as HtReads_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_barcode_gen():
    """Testing the generation of barcodes
    """
    n_barcodes = 1000
    barcodes = SimHtReads.barcodes(n_barcodes)
    assert len(barcodes) == n_barcodes
    assert type(barcodes) is np.ndarray


def test_main():
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(data_dir, 'temp_read_files')
    output_prefix = os.path.join(data_dir, 'TEST')
    args = ['--debug',
            '--art-paired',
            '--tmp-dir', temp_dir,
            '--barcode-total', 20,
            '--barcode-chunks', 2,
            '--seq-depth', 1e3, 
            genome_table, abund_table, output_prefix]
    HtReads_CMD.opt_parse(args)


def test_main_multi():
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(data_dir, 'temp_read_files')
    output_prefix = os.path.join(data_dir, 'TEST')
    args = ['-n', 4,
            '--art-paired',
            '--tmp-dir', temp_dir,
            '--seq-depth', 1e4, 
            genome_table, abund_table, output_prefix]
    HtReads_CMD.opt_parse(args)
