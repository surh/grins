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


// Params
params.indir = ''
params.outdir = 'output/'
params.format = 'genbank'
params.w_size = 150
params.s_size = 30

// Process params
if(params.format == 'genbank'){
  suffix = 'gbk'
}else if(params.format == 'fasta'){
  suffix = 'fasta'
}else{
  error "Wrong format\n"
}

SEQS = Channel.fromPath("${params.indir}/*.${suffix}")
if(params.format == 'genbank'){
  SEQS.into(GBKS)
}else if(params.format == 'fasta'){
  SEQS.into(FASTAS)
  GBKS = Channel.empty()
}else{
  error "Wrong format\n"
}

process gbk2fasta{
  label 'py3'

  input:
  file gbk_file from GBKS

  output:
  file '*.fasta' into FASTAS

  """
  ${workflow.projectDir}/gbk2fasta.py --input $gbk_file
  """
}

// Split fasta squences
FASTAS.into(FORWINDOWS; FORINDEX)

process split_in_windows{
  label 'py3'

  input:
  file seq from FORWINDWS
  val w_size from params.w_size
  val s_size from params.s_size

  output:
  file '*_windows.fasta' into WINDOWS

  """
  ${workflow.projectDir}/split_seq_into_windows.py \
    --input $seq \
    --format fasta \
    --w_size $w_size \
    --s_size $s_size
  """
}

process bowtie2{
  label 'bowtie2'
  label 'samtools'
  publishDir "${params.outdir}/bam/", mode: 'rellink'

  input:
  file ref from FORINDEX
  file windows from WINDOWS

  output:
  file '*.bam' into BAMS

  """
  bowtie2-build $ref
  bowtie2 \
    -f \
    --end-to-end \
    --sensitive \
    -a \
    --time \
    --threads 1 \
    -x $ref \
    -U windows | \
    samtools view -b - > ${ref}.bam

  """
}

// ~/micropopgen/src/SeqTools/gbk2fasta.pl -i AB469193.gbk -o AB469193.fasta
// ~/micropopgen/src/grins/detection/produce_windows_from_bam.py --input mapping.bam --output mapping_grins.gff3
