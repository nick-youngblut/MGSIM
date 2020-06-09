![MGSIM](https://github.com/nick-youngblut/MGSIM/workflows/MGSIM/badge.svg)

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

See the `conda install` line in the [CI yaml](.github/workflows/pythonpackage.yml)

## Install

### via pip

`pip install MGSIM`

### via `setup.py`

`python setpy.py install`

## Testing

* conda-forge::pytest

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

### Simulating Illumina reads

`MGSIM reads -h`

### Simulating haplotagging reads

`MGSIM ht_reads -h`


# LICENSE

See [LICENSE](./LICENSE)


