MGSIM
======

Metagenome read simulation for multiple synthetic communities

# algorthim

* download genomes
  * input:
    * table = taxon_name,accession
  * output:
    * table = taxon_name,genome_fasta
* simulate community
  * input:
    * table = taxon_name,genome_fasta
  * output:
    * table = taxon_name,sampleN_abs-abund,...
    * table = taxon_name,genome_fasta,genome_size
* community to DNA pool
  * description
    * genome relative fraction in the total DNA pool
      * bp_DNA = `genome_size * abs-abund`
      * scaling bp_DNA for all genomes in sample
  * input
    * table = taxon_name,sampleN_rel-abund
    * table = taxon_name,genome_fasta,genome_size
  * output
    * table = taxon_name,rel_DNA_pool_size
* sampling from DNA pool
  * for each genome (parallel?)
    * for each sample (parallel)
      * `genome_n_seqs = sequencing_depth * rel_DNA_pool_size`
      * simulate sequences (write to temp file)
        * Illumina: art
	* Nanopore: nanosim
      * reformat read headers (name by taxon)
  * for each sample
    * combine all genome reads
    * create table = sample,fasta_file
  * input
    * table = taxon_name,rel_DNA_pool_size
    * param = sequencing depth
  * output
    * table = sample,metagenome_fasta
    * metagenome_fasta files
 

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


# TUTORIALS


# CHANGELOG

[[top](#sections)]


# LICENSE

[[top](#sections)]

