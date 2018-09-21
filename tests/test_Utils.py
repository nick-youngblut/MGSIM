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

# data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')


# tests
def test_file_status():
    # real file
    f = os.path.join(data_dir, 'Methanosarcina_barkeri_MS.fna')
    Utils.is_file(f)
    Utils.checkExists(f)
    Utils.checkEmpty(f)
    # fake file
    f = os.path.join(data_dir, 'DOES_NOT_EXIST')
    with pytest.raises(IOError):
        Utils.is_file(f)
    with pytest.raises(IOError):
        Utils.checkExists(f)
    
