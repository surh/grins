#!/usr/bin/env nextflow
// Copyright (C) 2019 Sur Herrera Paredes

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
// along with this program.  If not, see <http://www.gnu.org/licenses/>

params.names = ''
params.indir = ''
params.outdir = 'output'
params.min_size = 500
params.w_size = 150
params.s_size = 30

// Read file
names = file(params.names)
reader = names.newReader()
NAMES = []
INTERVALS = []
while(str = reader.readLine()){
  NAMES = NAMES + [tuple(str, file("${params.indir}/fasta/${str}.fasta"))]
  INTERVALS = INTERVALS + [tuple(str,
    file("${params.indir}/pGRINS.gff3/${str}.fasta.pgrins.gff3"))]
}
INTERVALS = Channel.from(INTERVALS)

process skews{
  label 'py3'
  publishDir "${params.outdir}/skews", mode: 'rellink'

  input:
  set name, file(fasta) from NAMES

  output:
  set name, file(fasta), file("${name}_skews.txt") into SKEWS

  """
  ${workflow.projectDir}/skews.py \
    --input $fasta \
    --prefix $name
  """
}

process cors{
  label 'py3'
  publishDir "${params.outdir}/cors", mode: 'rellink'

  input:
  set name, file(fasta), file(skews), file(gff3) from SKEWS.join(INTERVALS)

  output:
  file "${name}.fasta.pgrins_cors.bed"

  """
  ${workflow.projectDir}/skew_correlations.py \
    --input $gff3 \
    --skews $skews \
    --sequence $fasta \
    --include_complement
  """
}
