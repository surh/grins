# Phylogenetic analysis

Code for phylogeentic analysis of PKSs and other related sequences from
the GRINS project.

Most of the code was developed to run in fraser-server, though with the
appropriate modules installed it should be quite contained.

The best way to ensure you can run the software is by copying the full
contents of this directory to a given location.

# Files

* **seq2phylo.nf**: Nextflow pipeline that starts with DNA sequences and
produces a phylogeny. It can also start from aminoacid fasta sequences or
multiple sequence alignments. It is highly recommended to read the
documentation at the top of the script itself for more details.
The pipelines depends on lmod modules for RAxML (called `raxml`), and
either Clustal Omega (called `clustalo`) or MAFFT (called `mafft`).
Python 3 with Biopython must also be on the user PATH. For basic usage choose
one of the following options:

```bash
nextflow run seq2phylo.nf --nuc_dir <MY_NUC_DIR>
nextflow run seq2phylo.nf --faa_dir <MY_FAA_DIR>
nextflow run seq2phylo.nf --aln_dir <MY_ALN_DIR>
```

* **translate.py**: Python 3 script that takes a fasta file of nucleotide
sequences and translates them into aminoacid sequences. It requires BioPython.
For simple usage run:

```bash
translate.py --infile <dnaseqs.fasta> --outfile <aaseqs.fasta>
```

For full details and documentation use:

```bash
translate.py --help
```

* **rename_dist_matrix.py**: Script that reads tab-delimited file where the
first column ar row are coloumn and row IDs. It copies the file while renaming
both columns and rows. For usage details use:

```bash
rename_dist_matrix.py
```
