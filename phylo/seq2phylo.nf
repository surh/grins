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
params.bindir = '~/micropopgen/src/grins/phylo/'

// Process paramters
nuc_dir = params.nuc_dir
faa_dir = params.faa_dir
aln_dir = params.aln_dir

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
    ${params.bindir}/translate.py \
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
      fromPath("${params.aln_dir}/*${params.aln_extesion}").
      map{file -> tuple(file.baseName, file)}
  }

  process raxml{
    module 'raxml'
    publishDir "${params.outdir}/TRE/", mode: 'copy'
    cpus 3

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
      -T 3 \
      -n $filename
    """
  }
}
