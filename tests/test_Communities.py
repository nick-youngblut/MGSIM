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
def test_main():
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(data_dir, 'comm')
    args = [genomeList, output_prefix, '--n-comm', 2]
    Communities_CMD.opt_parse(args)

def test_main_ncomm1():
    genomeList = os.path.join(data_dir, 'genome_list.txt')
    output_prefix = os.path.join(data_dir, 'comm')
    args = [genomeList, output_prefix, '--n-comm', 1, '--richness', 1]
    Communities_CMD.opt_parse(args)
