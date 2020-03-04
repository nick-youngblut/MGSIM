[![Travis-CI Build Status](https://travis-ci.org/nick-youngblut/MGSIM.svg?branch=master)](https://travis-ci.org/nick-youngblut/MGSIM)

MGSIM
=====

Metagenome read simulation of multiple synthetic communities

#### Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIMULATION WORKFLOW](#simulation_workflow)
- [CHANGE LOG](#changelog)
- [LICENSE](#license)


# REFERENCE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696891.svg)](https://doi.org/10.5281/zenodo.3696891)


# INSTALLATION

## Dependencies

See the `.travis.yml` file for setup info.

## Install

It is best to used a conda environment, but you can just run the following:

`python setpy.py install`

> In theory, dependencies that are NOT found in your conda
environment will be downloaded during the `setup.py install`

## Testing

* conda-forge::pytest

In the MGSIM base directory, use the command `pytest` to
run all of the tests.

To run tests on a particular test file:

`pytest -s path/to/the/test/file`

Example:

`pytest -s ./tests/test_Genome_download.py`

# HOW-TO

See all subcommands:

`MGSIM --list`

## Download genomes

`MGSIM genome_download -h`

## Simulate communities

`MGSIM communities -h`

## Simulate reads for each genome in each community

### Simulating Illumina reads

`MGSIM reads -h`

### Simulating haplotagging reads

`MGSIM ht_reads -h`


# LICENSE

See [LICENSE](./LICENSE)


