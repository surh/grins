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

/* ############## DOCUMENTATION #################
This Nextflow pipeline takes a directory that has a collection of
sequence files and, for each file in the directory, it performs
translation, alignment and phylogenetic reconstruction.

The simplest way to utilize it is by simply specifying the directory
with the sequence files (either, nucleotide fasta, aminoacid fasta or
aligmnent) with one of the following commands:

$ nextflow run seq2phylo.nf --nuc_dir <MY_NUC_DIR>
$ nextflow run seq2phylo.nf --faa_dir <MY_FAA_DIR>
$ nextflow run seq2phylo.nf --aln_dir <MY_ALN_DIR>

The script creates an output directory with the results of every step. A
full description of the paramters follows below:
--nuc_dir, --faa_dir, --aln_dir
One and only one of these must be specified. Directory containing nucleotide
fasta, aminoacid fasta or alignment files.

--nuc_extension --aa_extension, --aln_extension
Extension for nucleotide fasta, aminoacid fasta or alignment files. Defaults
are '.fasta', '.faa' and '.aln' respectively.

--outdir
Directory to create for results. Default 'out/'

--aligner
Either 'mafft' or  'clustalo'. Program to make multiple sequence
alignments. Default 'mafft'.

--bootstrap
Bootstrap replicates to perform during phylogenetic reconstruction. See
RAxML manual for more info. Typically an integer number or 'autoMRE'.
Default: 100.

--table
Genetic code table for translation. Default 'Standard'.

--aln_mem, --aln_time, --aln_threads
--phylo_mem, --phylo_time, --phylo_threads
Parameters that control computational resource allocation for alignments
and phylogenetic reconstruction. Defaults are 24 hrs with 2GB of memory
single thread for alignment and three threads for phylogenetic
reconstruction.
################################################# */

// Parameters
params.nuc_dir = ''
params.faa_dir = ''
params.aln_dir = ''

params.nuc_extension = ".fasta"
params.aa_extension = ".faa"
params.aln_extension = ".aln"

params.outdir = "out/"
params.aligner = 'mafft'
// params.aligner = 'clustalo'

params.bootstrap = '100'
params.table = 'Standard'

params.aln_time = '24:00:00'
params.aln_mem = '2GB'
params.aln_threads = 1
params.phylo_time = '24:00:00'
params.phylo_mem = '2GB'
params.phylo_threads = 3

// ############## Process paramters #############
nuc_dir = params.nuc_dir
faa_dir = params.faa_dir
aln_dir = params.aln_dir

// ################# Run checks #################
// Check that only one dir is passed
n_dirs = 0
mydir = ''
if (nuc_dir != ''){
  n_dirs = n_dirs + 1
  mydir = nuc_dir
}
if (faa_dir != ''){
  n_dirs = n_dirs + 1
  mydir = faa_dir
}
if (aln_dir != ''){
  n_dirs = n_dirs + 1
  mydir = aln_dir
}
if(n_dirs != 1){
  error """ERROR: One and only one input directory (nuc_dir, faa_dir, \
    aln_dir) must be passed."""
}

// Check that passed directory is a directory and contains files
myfiles = file(mydir).listFiles()
if (myfiles == null){
  error """ERROR: The provided directory doesn't look like a directory."""
} else if(myfiles == []){
  error """ERROR: The provided directory is empty"""
}

// Check modules are installed
// Check RAxML
out = """
module () {
  eval \$($LMOD_CMD bash "$@") && eval $(\${LMOD_SETTARG_CMD:-:} -s sh)
}
module is-avail raxml""".execute()

if (out.exitValue() != 0){
  error """ERROR: there is no 'raxml' module"""
}

// Check aligner
cmd = "module is-avail ${params.aligner}"
out = cmd.execute()
if (out.exitValue() != 0){
  error """ERROR: there is no ${params.aligner} module"""
}

// ################ Run pipeline ################

// Determine which output is present
if (nuc_dir != ''){
  NUCS = Channel.
    fromPath("${nuc_dir}/*${params.nuc_extension}").
    map{file -> tuple(file.baseName, file)}

  // Processes
  process translate{
    publishDir "${params.outdir}/FAA/", mode: 'copy'

    input:
    set filename, file(seqs) from NUCS

    output:
    set filename, file("${filename}.faa") into FAAS
    file 'sequence_names_map.txt'

    """
    ${workflow.projectDir}/translate.py \
      --infile $seqs \
      --remove_stops \
      --outfile ${filename}.faa \
      --rename
    """
  }

  // Specify output directory
  faa_dir = "${params.outdir}/FAA/"
}

if (faa_dir != ''){
  if (params.faa_dir != ''){
    // If we are starting from here, read files
    FAAS = Channel.
      fromPath("${params.faa_dir}/*${params.aa_extension}").
      map{file -> tuple(file.baseName, file)}
  }

  process align{
    publishDir "${params.outdir}/ALN/", mode: 'copy'
    module params.aligner
    cpus params.aln_threads
    memory params.aln_mem
    time params.aln_time

    input:
    set filename, file(seqs) from FAAS

    output:
    set filename, file("${filename}.aln") into ALNS

    script:
    if(params.aligner == 'clustalo'){
      """
      clustalo \
        -i $seqs \
        -o ${filename}.aln \
        -t Protein \
        -v
      """
    }else if(params.aligner == 'mafft'){
      """
      mafft \
        --auto \
        --amino \
        $seqs \
        > ${filename}.aln
      """
    }else{
      error "Invalid aligner"
    }
  }

  // Specify output directory
  aln_dir = "${params.outdir}/ALN/"
}

if (aln_dir != ''){
  if (params.faa_dir != ''){
    // If we are starting from here, read files
    ALNS = Channel.
      fromPath("${params.aln_dir}/*${params.aln_extension}").
      map{file -> tuple(file.baseName, file)}
  }

  process raxml{
    module 'raxml'
    publishDir "${params.outdir}/TRE/", mode: 'copy'
    cpus params.phylo_threads
    memory params.phylo_mem
    time params.phylo_time

    input:
    set filename, file(aln) from ALNS

    output:
    file "RAxML_bestTree.${filename}"
    file "RAxML_bipartitionsBranchLabels.${filename}"
    file "RAxML_bipartitions.${filename}"
    file "RAxML_bootstrap.${filename}"

    """
    raxmlHPC-PTHREADS \
      -s $aln \
      -f a \
      -x 12345 \
      -p 12345 \
      -# autoMRE \
      -m PROTGAMMAAUTO \
      -T ${params.phylo_threads} \
      -n $filename
    """
  }
}
