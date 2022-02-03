#!/usr/bin/env python
# import
from __future__ import print_function
## batteries
import os
import sys
import pytest
import subprocess
## 3rd party
import numpy as np
import pandas as pd
## package
from MGSIM import Utils
from MGSIM.Commands import Reads as Reads_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

def validate_fastq(out_dir, paired=True, seq_type='illumina', suffix='.fq'):
    # validating fastq
    output_R = os.path.join(out_dir, seq_type, '1', 'R1' + suffix)
    ret = subprocess.run(['fqtools', 'validate', output_R])
    assert ret.returncode == 0
    if paired is True:
        output_R = os.path.join(out_dir, seq_type, '1', 'R2' + suffix)
        ret = subprocess.run(['fqtools', 'validate', output_R])
        assert ret.returncode == 0

def chk_suc(ret):
    # checking cmd return
    if not ret.success:
        ret.print()
    assert ret.success
    
        
# tests
def test_help(script_runner):
    ret = script_runner.run('MGSIM', 'reads', '-h')
    chk_suc(ret)
    
def test_main(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            '--gzip',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, suffix='.fq.gz')

def test_main_multiProc(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '-n', '2',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix)

def test_main_unpaired(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-mflen', '0',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, paired=False)

def test_main_low_qual(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-paired',
                            '--art-maxQ', '25',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, paired=False)

def test_main_art_error_profiles(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    error_prof1 = os.path.join(data_dir, 'art_profiles/HiSeq2500L150R1.txt')
    error_prof2 = os.path.join(data_dir, 'art_profiles/HiSeq2500L150R2.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '--art-paired',
                            '--art-qprof1', error_prof1,                            
                            '--art-qprof2', error_prof2,
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1e4',
                            '--rndSeed', '8294',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, paired=False)
    
def test_main_pacbio(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '-n', '2',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1000',
                            '--pb-seq-depth', '100',
                            '--rndSeed', '8294',
                            '--gzip',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, suffix='.fq.gz')
    validate_fastq(output_prefix, paired=False, seq_type='pacbio', suffix='.fq.gz')

def test_main_nanopore(script_runner, tmp_path):
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    abund_table = os.path.join(data_dir, 'comm_wAbund.txt')
    temp_dir = os.path.join(str(tmp_path), 'temp_read_files')
    output_prefix = os.path.join(str(tmp_path), 'reads_out')
    ret = script_runner.run('MGSIM', 'reads',
                            '-n', '2',
                            '--art-paired',
                            '--tmp-dir', temp_dir,
                            '--sr-seq-depth', '1000',
                            '--np-seq-depth', '100',
                            '--rndSeed', '8294',
                            '--gzip',
                            genomeList, abund_table, output_prefix)
    chk_suc(ret)
    validate_fastq(output_prefix, suffix='.fq.gz')
    validate_fastq(output_prefix, paired=False, seq_type='nanopore', suffix='.fa.gz')
