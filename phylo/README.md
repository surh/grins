# Phylogenetic analysis

Code for phylogeentic analysis of PKSs and other related sequences from
the GRINS project.

Most of the code was developed to run in fraser-server, though with the
appropriate modules installed it should be quite contained,

# Files

* **raxmk_all_nuc.nf**: Nextflow pipeline to run phylogeny in all alingnments
in a directory. For every file in the input directory (`--indir`) with
extension `.aln`. It will call `raxmlHPC` (`--bin`). It assumes that
RAxML v8.2.X can be loaded with `module load raxml`. For each allignment,
the RAxML command  will first conduct a bootstrap search and a search for
the best–scoring maximum likelihood tree. One can specify the number of
bootstrap replicates or specify a convervegence criteria with option
`--botstrap`. An example call would be:

```bash
nextflow run raxml_all_nuc.nf --bootstrap autoMRE --indir myalns/ --seed 12345
```

## TO DO LIST

1. Perform inference from aminoacide alignment
2. Include alignment in pipeline (clustalo or mafft)
3. Add alignment QC before inference
4. Automatic model selection for phylogeny?
