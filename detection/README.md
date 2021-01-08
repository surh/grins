This folder contains the scripts for detecting GRINS

# Key files

* **antismash2gff3.py**: Converts antiSMASH 5 BGC predictions to GFF3 format.

* **detect_in_genome.nf**: Pipeline to detect duplicated regions in
bacterial genomes in fasta format. Includes annotations via antiSMASH 5.

* **GRINS_detection_from_BOWTIE.py**: Takes duplicated regions identified by
`detect_in_genome.nf` and identifies GRINS regions.

* **GRINS_detection_from_BOWTIE.nf**: Pipeline that parallelizes the python
script of the same name.

* **GRINS_detection_pairwise_identity.py**: Script to detect GRINS within
PKS clusters.

* **produce_windows_from_bam.py**: Reads bowtie2 BAM output and identifies
duplicated regions.

* **split_seq_into_windows.py**: Creates sliding windows for all analysis.
