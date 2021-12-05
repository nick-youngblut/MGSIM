[![MGSIM](https://github.com/nick-youngblut/MGSIM/actions/workflows/pythonpackage.yml/badge.svg)](https://github.com/nick-youngblut/MGSIM/actions/workflows/pythonpackage.yml)
[![Upload Python Package](https://github.com/nick-youngblut/MGSIM/actions/workflows/python-publish.yml/badge.svg)](https://github.com/nick-youngblut/MGSIM/actions/workflows/python-publish.yml)
[![PyPI version](https://badge.fury.io/py/MGSIM.svg)](https://badge.fury.io/py/MGSIM)

MGSIM
=====

Metagenome read simulation of multiple synthetic communities

#### Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIMULATION WORKFLOW](#simulation_workflow)
- [LICENSE](#license)


# REFERENCE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3696891.svg)](https://doi.org/10.5281/zenodo.3696891)

# DESCRIPTION

Straight-forward simulations of metagenome data from a
collection of reference bacterial/archaeal genomes. 

## Highlights

* Can simulate Illumina, PacBio, and/or Nanopore reads
  * For Illumina, synthetic long reads (read clouds) can also be simulated
* Generate communities differing in:
  * Sequencing depth
  * Richness
  * Beta diversity
  
The workflow:

* [optional] Download reference genomes
* Format reference genomes
  * e.g., rename contigs
* Simulate communities
* Simulate reads for each community 

# INSTALLATION

## Dependencies

See the `conda install` line in the [CI yaml](.github/workflows/pythonpackage.yml)

## Install

### via pip

`pip install MGSIM`

### via `setup.py`

`python setpy.py install`

## Testing

* conda-forge::pytest
* conda-forge::pytest-console-scripts

In the MGSIM base directory, use the command `pytest` to
run all of the tests.

To run tests on a particular test file:

`pytest -s --script-launch-mode=subprocess  path/to/the/test/file`

Example:

`pytest -s --script-launch-mode=subprocess ./tests/test_Reads.py`

# HOW-TO

See all subcommands:

`MGSIM --list`

## Download genomes

`MGSIM genome_download -h`

## Simulate communities

`MGSIM communities -h`

## Simulate reads for each genome in each community

### Simulating Illumina, PacBio, and/or Nanopore reads

`MGSIM reads -h`

### Simulating haplotagging reads (aka read-cloud data)

`MGSIM ht_reads -h`


# LICENSE

See [LICENSE](./LICENSE)


