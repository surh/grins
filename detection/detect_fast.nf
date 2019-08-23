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
params.vsearch_id = 0.9

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
  SEQS.into{GBKS}
  FASTAS = Channel.empty()
}else if(params.format == 'fasta'){
  SEQS.into{FASTAS}
  GBKS = Channel.empty()
}else{
  error "Wrong format\n"
}

process gbk2fasta{
  label 'py3'
  publishDir "${params.outdir}/fasta/", mode: 'rellink'

  input:
  file gbk_file from GBKS

  output:
  file '*.fasta' into FASTASFROMGBK

  """
  ${workflow.projectDir}/gbk2fasta.py --input $gbk_file
  """
}

// Split fasta squences
FASTASFROMGBK.mix(FASTAS).into{FORWINDOWS; FORINDEX; FORGFFFASTA; FORPLOTS}

process split_in_windows{
  label 'py3'

  input:
  file seq from FORWINDOWS
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
  file '*.bam' into BAMS, BAMS_FOR_PLOT
  val "$ref" into REFNAMES_FROM_BOWTIE2

  """
  bowtie2-build $ref $ref
  bowtie2 \
    -f \
    --end-to-end \
    --sensitive \
    -a \
    --time \
    --threads 1 \
    -x $ref \
    -U $windows | \
    samtools view -b - > ${ref}.bam
  """
}

process merge_bam_windows{
  label 'py3'
  publishDir "${params.outdir}/pGRINS.gff3/", mode: 'rellink'

  input:
  file bam from BAMS
  val ref from REFNAMES_FROM_BOWTIE2

  output:
  file "${ref}.pgrins.gff3" into GFF3, GFF3_FOR_PLOT
  val ref into REFNAMES_FROM_MERGEBAM

  """
  ${workflow.projectDir}/produce_windows_from_bam.py \
    --input $bam \
    --output ${ref}.pgrins.gff3 \
    --w_size ${params.w_size}
  """
}

process gff_to_fasta{
  label 'py3'
  publishDir "${params.outdir}/pGRINS.fasta", mode: 'rellink'

  input:
  file gff from GFF3
  file fasta from FORGFFFASTA
  val ref from REFNAMES_FROM_MERGEBAM

  output:
  file "${ref}.pgrins.fasta" into PGRINS_FASTA
  val ref into REFNAMES_FROM_GFFFASTA

  """
  ${workflow.projectDir}/gff_to_fasta.py \
    --input $fasta \
    --gff3 $gff \
    --output ${ref}.pgrins.fasta \
    --format fasta
  """
}

process vsearch_pgrins{
  label 'vsearch'
  publishDir "${params.outdir}/vsearch.pgrins/", mode: 'rellink'

  input:
  file fasta from PGRINS_FASTA
  val ref from REFNAMES_FROM_GFFFASTA

  output:
  file "${ref}.pgrins.centroids.fasta" into CENTROIDS
  file "${ref}.clusters.uc" into PGRINS_UC
  val ref into REFNAMES_FROM_VSEARCH

  """
  vsearch \
    --threads 1 \
    --sizeout \
    --cluster_fast $fasta \
    --strand both \
    --centroids ${ref}.pgrins.centroids.fasta \
    --relabel centroid_ \
    --sizeorder \
    --uc ${ref}.clusters.uc \
    --id $params.vsearch_id
  """

}

process plot_fast_grins{
  label 'py3'
  publishDir "${params.outdir}/plots/", mode: 'rellink'

  input:
  file sequence from FORPLOTS
  file grins_gff3 from GFF3_FOR_PLOT
  file grins_clusters from PGRINS_UC
  file windows_bam from BAMS_FOR_PLOT
  val ref from REFNAMES_FROM_VSEARCH

  output:
  file "${ref}.png" into PLOTS

  """
  ${workflow.projectDir}/plot_fast_grins.py \
    --input $sequence \
    --grins_gff3 $grins_gff3 \
    --grins_clusters $grins_clusters \
    --windows_bam $windows_bam \
    --format fasta \
    --w_size $params.w_size \
    --s_size $params.s_size \
    --output ${ref}.png
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxFors = 300
  errorStrategy = 'finish'
  withLabel: 'py3'{
    module = 'conda/4.6.14'
    conda = '/home/groups/hbfraser/modules/packages/conda/4.6.14/envs/fraserconda'
  }
  withLabel: 'bowtie2'{
    module = 'bowtie2'
  }
  withLabel: 'samtools'{
    module = 'samtools'
  }
  withLabel: 'vsearch'{
    module = 'vsearch'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
