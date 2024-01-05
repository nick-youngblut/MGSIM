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

See `environment.yml` for a list of dependencies.

You can install via:

```bash
mamba env create -f environment.yml -n mgsim
```

> mamba is much faster than conda

## Install

### via pip

`pip install MGSIM`

### via `setup.py`

`python setpy.py install`

## Testing

* `conda-forge::pytest>=5.3`
* `conda-forge::pytest-console-scripts>=1.2`

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

# Tutorial

## Reference genome download

Create Taxon-accession table

```bash
mkdir -p tutorial

cat <<-EOF > tutorial/taxon_accession.tsv
Taxon	Accession
Escherichia coli O104-H4	NC_018658.1
Clostridium perfringens ATCC.13124	NC_008261
Methanosarcina barkeri [MS]	NZ_CP009528.1
EOF
```

Download genomes

```bash
MGSIM genome_download -d tutorial/ tutorial/taxon_accession.tsv > tutorial/genomes.tsv
```

## Simulate communities

Simulate 2 communities

```bash
MGSIM communities --n-comm 2 tutorial/genomes.tsv tutorial/communities
```

## Simulate reads

Illumina reads

```bash
MGSIM reads tutorial/genomes.tsv --sr-seq-depth 1e5 tutorial/communities_abund.txt tutorial/illumina_reads/
```

PacBio reads

```bash
MGSIM reads tutorial/genomes.tsv --pb-seq-depth 1e3 tutorial/communities_abund.txt tutorial/pacbio_reads/
```

Nanopore reads

```bash
MGSIM reads tutorial/genomes.tsv --np-seq-depth 1e3 tutorial/communities_abund.txt tutorial/nanopore_reads/
```


# LICENSE

See [LICENSE](./LICENSE)


