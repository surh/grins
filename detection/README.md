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

The recommended setup is to have conda installed as an LMOD module, and then
install dependencies as either other LMOD modules or conda environments.

* **antiSMASH (5.1.1 or 5.1.2)**: We recommend installing antiSMASH to its
own conda environment. See [their documentation](https://docs.antismash.secondarymetabolites.org/install/).
See also [our conda environment .yml file](../conda_envs/antismash.yml).
