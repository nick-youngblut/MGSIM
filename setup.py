#!/usr/bin/env python
from setuptools import setup, find_packages

## install main application
desc = 'Metagenome simulation of multiple synthetic communities'
setup(
    name = 'MGSIM',
    version = '0.2.7',
    description = desc,
    long_description = desc + '\n See README for more information.',
    author = 'Nick Youngblut',
    author_email = 'nyoungb2@gmail.com',
    entry_points={
        'console_scripts': [
            'MGSIM = MGSIM.__main__:main'
        ]
    },
    license = "MIT license",
    packages = find_packages(),
    package_dir={'MGSIM':
                 'MGSIM'},
    url = 'https://github.com/nick-youngblut/MGSIM'
)




