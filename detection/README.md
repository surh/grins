# Detecting GRINS
This folder contains the scripts for detecting GRINS

## Key files

* **antismash2gff3.py**: Converts antiSMASH 5 BGC predictions to GFF3 format.

* **detect_in_genome.nf**: Pipeline to detect duplicated regions in
bacterial genomes in fasta format. Includes annotations via antiSMASH 5.

* **GRINS_detection_from_BOWTIE.py**: Takes duplicated regions identified by
`detect_in_genome.nf` and identifies GRINS regions.

* **GRINS_detection_from_BOWTIE.nf**: Pipeline that parallelizes the python
script of the same name.

* **GRINS_detection_pairwise_identity.py**: Script to detect GRINS within
PKS clusters.

* **pairwise_alignment_windows.nf**: Pipeline that parallelizes the python
script of the same name.

* **pairwise_alignment_windows.py**: Script that splits PKS cluster sequences
into windows and performs pairwise alignments which can be used to identify
duplicated regions and GRINS.

* **produce_windows_from_bam.py**: Reads bowtie2 BAM output and identifies
duplicated regions.

* **split_seq_into_windows.py**: Creates sliding windows for all analysis.


## Dependencies

The following software is required to run the current GRINS detection workflow.
The versions specified are those that were tested. It is possbile that the
workflow will run with similar versions but it has not been tested.

The recommended setup is to have conda installed as an
[Lmod](https://lmod.readthedocs.io/en/latest/index.html) module, and then
install dependencies as either other Lmod modules or conda environments.

* **antiSMASH (5.1.1 or 5.1.2)**: We recommend installing antiSMASH to its
own conda environment. See [their documentation](https://docs.antismash.secondarymetabolites.org/install/).
See also [our conda environment .yml file](../conda_envs/antismash.yml).

* **Python (3.6.13)**: Python 3 is required including the following modules:
argparse, biopython, json, os, pybedtools. We recommend to have a conda environment with
these dependencies. See our [conda environment .yml file](../conda_envs/fraserconda.yml)
for precise versions. Please note that this environment include many additional
Python 3 modules which are not required.

* **bowtie2 (2.4.1)**: Check [their website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for details.
We recommend compiling locally and installing bowtie2 as an Lmod module.

* **samtools (1.7 or 1.10)**: Check [their website](http://www.htslib.org/) for details.
We recommend compiling locally and installing samtools as an Lmod module.

* **bedtools (2.30)**: Check [their website](https://bedtools.readthedocs.io/en/latest/index.html) for details.
Since pybedtools requires this software, you might have installed it when you
created the Python 3 conda environment.
