#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import uuid
import pytest
import subprocess
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
tmp_dir = os.path.join(test_dir, 'tmp')
out_dir = os.path.join(test_dir, 'output')

def chk_suc(ret):
    # checking cmd return
    if not ret.success:
        ret.print()
    assert ret.success

def validate_fastq(base_dir):
    # validating fastq
    ## Read1
    output_R = os.path.join(base_dir, '1', 'R1.fq')
    if not os.path.isfile(output_R):
        raise IOError('Cannot find file: {}'.format(output_R))
    ret = subprocess.run(['fqtools', 'validate', output_R])
    assert ret.returncode == 0
    ## Read2
    output_R = os.path.join(base_dir, '1', 'R2.fq')
    if not os.path.isfile(output_R):
        raise IOError('Cannot find file: {}'.format(output_R))
    ret = subprocess.run(['fqtools', 'validate', output_R])
    assert ret.returncode == 0

# tests
def test_barcode_gen():
    """
    Testing the generation of barcodes
    """
    n_barcodes = 1000
    barcodes = SimHtReads.barcodes(n_barcodes)
    assert len(barcodes) == n_barcodes
    assert type(barcodes) is np.ndarray

def test_help(script_runner):
    ret = script_runner.run('MGSIM', 'ht_reads', '-h')
    chk_suc(ret)

def test_main(script_runner, tmp_path):
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), str(uuid.uuid4()))
    output_prefix = os.path.join(str(tmp_path), 'TEST-main')
    
    ret = script_runner.run('MGSIM', 'ht_reads', '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--barcode-total', '20',
                            '--barcode-chunks', '2',
                            '--seq-depth', '1e3',
                            '--rndSeed', '8294',
                            genome_table, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)

def test_main_multi(script_runner, tmp_path):
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), str(uuid.uuid4()))
    output_prefix = os.path.join(str(tmp_path), 'TEST-main-multi')
    
    ret = script_runner.run('MGSIM', 'ht_reads',
                            '--art-paired', '-n', '2',
                            '--tmp-dir', temp_dir,
                            '--seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genome_table, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)

def test_main_zeros(script_runner, tmp_path):
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund_zeros.txt')
    temp_dir = os.path.join(str(tmp_path), str(uuid.uuid4()))
    output_prefix = os.path.join(str(tmp_path), 'TEST-main-zeros')
    
    ret = script_runner.run('MGSIM', 'ht_reads', '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--barcode-total', '20',
                            '--barcode-chunks', '2',
                            '--seq-depth', '1e3',
                            '--rndSeed', '8294',
                            genome_table, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)
    
def test_main_prefix(script_runner, tmp_path):
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), str(uuid.uuid4()))
    output_prefix = os.path.join(str(tmp_path), 'TEST-main-prefix')
    
    ret = script_runner.run('MGSIM', 'ht_reads', '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--barcode-total', '20',
                            '--barcode-chunks', '2',
                            '--seq-depth', '1e3',
                            '--read-name', '{readID}_BX:Z:{barcodeID}',
                            '--rndSeed', '8294',
                            genome_table, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)

def test_main_large_frag_sd(script_runner, tmp_path):
    genome_table = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), str(uuid.uuid4()))
    output_prefix = os.path.join(str(tmp_path), 'TEST-main-largeFragSize')
    
    ret = script_runner.run('MGSIM', 'ht_reads', '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--frag-size-mean', '10000',
                            '--frag-size-sd', '9000',
                            '--barcode-total', '20',
                            '--barcode-chunks', '2',
                            '--seq-depth', '5e3',
                            genome_table, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)     
