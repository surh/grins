#!/usr/bin/env nextflow
// Copyright (C) 2018 Sur Herrera Paredes

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Nextflow pipeline that runs a directory of genomes through checkm.
// Combines results in single table

// Parameters:
// --indir
// Input directory with GENOME files. Assumes genomes are organized in
// an an arbitrary number of subdirectories. And each genome is in a
// second level subdirectory, where there is an '.fna' file with the
// same name as the genome directory. For example, the following structure
// would be acceptable:
// genomes/group1/genome1/genome1.fna
// genomes/group1/genome2/genome2.fna
// genomes/group2/genome3/genome3.fna
// genomes/group2/genome4/genome4.fna

// --outdir
// Name of output directory

// --contamination, --completeness
// Thresholds for checkM statistics

// --batch_size
// Number of genomes to analyze by batch

// --threads, --memory, --time, --queue, --max_forks
// Parameters for checkm jobs.

// Configuration and dependencies:
// The processes here might require access to python 3, checkm, and R.
// Use the labeles 'py3', 'checkm', and 'r' in your nextflow.config
// file to provide access. See the end of this pipeline for an example file.

// Params
params.indir = 'genomes'
// params.outfile = 'checkm_results.txt'
params.outdir = 'output/'
params.contamination = 2
params.completeness = 98
params.heterogeneity = 0

params.batch_size = 200

// Process params
indir = file(params.indir)

process create_batch_map{
  label 'py3'

  input:
  file indir

  output:
  file('batch_map.txt') into create_batch
  file('checkm_batches/batch_*') into checkm_dirs

  """
  ${workflow.projectDir}/create_batches.py --indir ${params.indir} \
    --outdir checkm_batches \
    --outfile batch_map.txt \
    --batch_size ${params.batch_size}
  """
}

process run_checkm{
  label 'checkm'
  errorStrategy 'retry'
  maxRetries 2

  input:
  file checkm_dir from checkm_dirs.flatten()

  output:
  file "checkm_results.txt" into checkm_results
  // file "checkm_results_noheader.txt" into checkm_results_noheader

  """
  checkm lineage_wf \
    -t 8 \
    -f checkm_results.txt \
    --tab_table \
    ${checkm_dir} \
    results

  # tail -n +2 checkm_results.txt > checkm_results_noheader.txt
  """
}

process collect_results{
  label 'py3'
  publishDir params.outdir

  input:
  file "*.txt" from checkm_results.collect()

  output:
  file "checkm_results.txt" into CHECKM

  """
  ${workflow.projectDir}/../sutilspy/bin/cat_tables.py \
    *.txt \
    --outfile checkm_results.txt
  """
}

process filter_checkm{
  label 'r'
  publishDir params.outdir, pattern: "output/all_completeness_histogram.svg",
    saveAs: {"checkm/all_completeness_histogram.svg"}
  publishDir params.outdir, pattern: "output/all_contamination_histogram.svg",
    saveAs: {"checkm/all_contamination_histogram.svg"}
  publishDir params.outdir, pattern: "output/all_heterogeneity_histogram.svg",
    saveAs: {"checkm/all_heterogeneity_histogram.svg"}
  publishDir params.outdir,
    pattern: "output/all_completeness_vs_contamination.svg",
    saveAs: {"checkm/all_completeness_vs_contamination.svg"}
  publishDir params.outdir,
    pattern: "output/chosen_completeness_histogram.svg",
    saveAs: {"checkm/chosen_completeness_histogram.svg"}
  publishDir params.outdir,
    pattern: "output/chosen_contamination_histogram.svg",
    saveAs: {"checkm/chosen_contamination_histogram.svg"}
  publishDir params.outdir,
    pattern: "output/chosen_heterogeneity_histogram.svg",
    saveAs: {"checkm/chosen_heterogeneity_histogram.svg"}
  publishDir params.outdir,
    pattern: "output/chosen_completeness_vs_contamination.svg",
    saveAs: {"checkm/chosen_completeness_vs_contamination.svg"}
  publishDir params.outdir, pattern: "output/chosen_checkm_results.txt",
  saveAs: {"checkm/chosen_checkm_results.txt"}

  input:
  file "checkm_results.txt" from CHECKM

  output:
  file "output/all_completeness_histogram.svg"
  file "output/all_contamination_histogram.svg"
  file "output/all_heterogeneity_histogram.svg"
  file "output/all_completeness_vs_contamination.svg"
  file "output/chosen_completeness_histogram.svg"
  file "output/chosen_contamination_histogram.svg"
  file "output/chosen_heterogeneity_histogram.svg"
  file "output/chosen_completeness_vs_contamination.svg"
  file "output/chosen_checkm_results.txt" into CHOSEN

  """
  ${workflow.projectDir}/process_checkm_results.r \
    checkm_results.txt \
    --contamination ${params.contamination} \
    --completeness ${params.completeness} \
    --heterogeneity ${params.heterogeneity}
  """
}

process clean_genome_dirs{
  label 'py3'
  publishDir params.outdir,
    pattern: "clean",
    saveAs: {"clean"}

  input:
  file indir
  file 'checkm_results.txt' from CHOSEN

  output:
  file 'clean'

  """
  grep -v "# genomes" checkm_results.txt | cut -f 1 > keep.txt

  ${workflow.projectDir}/clean_genome_dirs.py ${params.indir} \
    --outdir clean/ \
    --keep keep.txt
  """
}


// Example nextflow.config
/*
process {
  queue = 'hbfraser,hns'
  maxForks = 100
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '2h'
  memory = '1G'
  withLabel: py3 {
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
  }
  withLabel: checkm {
    // module = 'prodigal:hmmer:pplacer:fraserconda'
    // conda = '/opt/modules/pkgs/anaconda/3.6/envs/python2'
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/checkm'
    cpus = 8
    memory = '40G'
    time = '2h'
  }
  withLabel: r {
    module = 'R/3.6.1'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}

*/
