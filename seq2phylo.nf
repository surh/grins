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


// Parameters
params.indir = 'seqs/'
params.nuc_dir = ''
params.faa_dir = ''
params.aln_dir = ''

params.nuc_extension = ".fasta"
params.aa_extension = ".faa"
params.aln_extension = ".aln"

params.outdir = "out/"
params.aligner = 'mafft'
params.aligner = 'clustalo'

params.bootstrap = '100'
params.table = 'Standard'
params.bindir = '~/micropopgen/src/grins/phylo/'

// Process paramters
// Determine which output is present
if (params.aln_dir != ''){
  ALNS = Channel.
    fromPath("${params.aln_dir}/*${params.aln_extesion}").
    map{file -> tuple(file.baseName, file)}
}else if (params.faa_dir != ''){
  FAAS = Channel.
    fromPath("${params.faa_dir}/*${params.aa_extension}").
    map{file -> tuple(file.baseName, file)}
}else if (params.nuc_dir != ''){
  NUCS = Channel.
    fromPath("${params.nuc_dir}/*${params.nuc_extension}").
    map{file -> tuple(file.baseName, file)}
}else{
  // Assume indir corresponds to nucleotide sequences
  NUCS = Channel.
    fromPath("${params.indir}/*${params.nuc_extension}").
    map{file -> tuple(file.baseName, file)}
}


// Processes
process translate{
  publishDir "${params.outdir}/FAA/", mode: 'copy'

  input:
  set filename, file(seqs) from NUCS

  output:
  set filename, file("${filename}.faa") into FAAS

  """
  ${params.bindir}/translate.py \
    --infile $seqs \
    --remove_stops \
    --outfile ${filename}.faa \
    --rename
  """
}


process translate{
  publishDir "${params.outdir}/ALN/", mode: 'copy'
  module params.aligner

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
