MGSIM
=====

Metagenome read simulation for multiple synthetic communities

#### Sections

- [REFERENCE](#reference)
- [INSTALLATION](#installation)
- [TUTORIALS](#tutorials)
- [SIMULATION WORKFLOW](#simulation_workflow)
- [CHANGE LOG](#changelog)
- [LICENSE](#license)


# REFERENCE

[[top](#sections)]


# INSTALLATION

[[top](#sections)]

## Dependencies

* python 3
* conda-forge::numpy
* conda-forge::pandas
* conda-forge::docopt
* conda-forge::scipy
* conda-forge::biopython
* bioconda::art

## Install

`python setpy.py install`

> In theory, dependencies are NOT found in your conda
environment will be downloaded during the `setup.py install`

## Testing

* conda-forge::pytest

In the MGSIM base directory, use the command `pytest` to
run all of the tests.

To run tests on a particular test file:

`pytest -s path/to/the/test/file`

# HOW-TO

## Download genomes

`MGSIM genome_download -h`

## Simulate communities

`MGSIM communities -h`

## Simulate reads for each genome in each community

`MGSIM reads -h`


# TODO

* handle taxon names that include spaces
* include logging

# CHANGELOG

[[top](#sections)]


# LICENSE

[[top](#sections)]

