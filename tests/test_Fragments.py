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
from MGSIM.Commands import Fragments as Fragments_CMD

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_cmd():
    genomeList = os.path.join(data_dir, 'genome_index.txt')
    fp = data_dir
    args = ['--fp', fp, genomeList]
    Fragments_CMD.opt_parse(args)
