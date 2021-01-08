# Annotating genes

## Key files

* **run_roary/**: This is a git submodule linking to the
[surh/run_roary](giithub.com/surh/run_roary) repo.

* **run_roary/gbk2faa.nf***: Extracts CDS sequences from antiSMASH 5 GBK files,
and converts to aminoacid fasta sequences.

* **run_eggnog.nf**: Runs eggNOG mapper with eggNOG v4.5 to annotate aminoacid
sequences infasta format.
