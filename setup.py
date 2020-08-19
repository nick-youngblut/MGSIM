#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import glob
import numpy

    
# dependencies
install_reqs = [
    'docopt>=0.6.2',
    'numpy>=1.10',
    'pandas>=0.18',
    'scipy>=0.17',
    'configobj>=5.0.6',
    'biopython>=1.68'
]

## install main application
desc = 'Metagenome simulation of multiple synthetic communities'
setup(
    name = 'MGSIM',
    version = '0.1.19',
    description = desc,
    long_description = desc + '\n See README for more information.',
    author = 'Nick Youngblut',
    author_email = 'nyoungb2@gmail.com',
    entry_points={
        'console_scripts': [
            'MGSIM = MGSIM.__main__:main'
        ]
    },
    install_requires = install_reqs,
    include_dirs = [numpy.get_include()],
    license = "MIT license",
    packages = find_packages(),
    package_dir={'MGSIM':
                 'MGSIM'},
    url = 'https://github.com/nick-youngblut/MGSIM'
)




