# Phylogenetic analysis

Code for phylogeentic analysis of PKSs and other related sequences from
the GRINS project.

Most of the code was developed to run in fraser-server, though with the
appropriate modules installed it should be quite contained,

# Files

* **raxml_all_nuc.nf**: *DEPRECATED: Has become redundant with seq2phylo.nf
and won't be maintained*. Nextflow pipeline to run phylogeny in all alingnments
in a directory. For every file in the input directory (`--indir`) with
extension `.aln`. By default, it will call `raxmlHPC` (`--bin`). It assumes that
RAxML v8.2.X can be loaded with `module load raxml`. For each allignment,
the RAxML command  will first conduct a bootstrap search and a search for
the best–scoring maximum likelihood tree. One can specify the number of
bootstrap replicates or specify a convervegence criteria with option
`--botstrap`. An example call would be:

```bash
nextflow run raxml_all_nuc.nf --bootstrap autoMRE --indir myalns/ --seed 12345
```

* **seq2phylo.nf**: Nextflow pipeline that starts with DNA sequences and
produces a phylogeny. It can also start from aminoacid fasta sequences or
multiple sequence alignments. See the documentation on the top of the script
itself for more details. It requires RAxML (phylogeny), clustalo or mafft
(alignment) and python 3 Biopython (translation). For basic usage use one of
the following options:

```bash
nextflow run seq2phylo.nf --nuc_dir <MY_NUC_DIR>
nextflow run seq2phylo.nf --faa_dir <MY_FAA_DIR>
nextflow run seq2phylo.nf --aln_dir <MY_ALN_DIR>
```
