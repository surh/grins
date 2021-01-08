# Runnig checkM

Pipeline to run checkM on bacterial genomes.

## Key files

* **checkm_batch.nf**: Split genomes into batches and run checkM in parallel
in all batches.

* **clean_genome_dirs.py**: ad-hoc script to produce a directory structure
with genomes passing checkM params.

* **create_batches.py**: Utility script to define batches of genomes for
checkM.

* **process_checkm_results.r**: Reads checkM results and selects genomes
based on user defined thresholds.
