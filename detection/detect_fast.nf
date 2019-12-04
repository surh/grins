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
params.sensitivity = 'sensitive'
params.w_size = 150
params.s_size = 30
params.vsearch_id = 0.9

// Process params
SEQS = Channel.fromPath("${params.indir}/*", type: 'file')
if(params.format == 'genbank'){
  SEQS.into{GBKS}
  FASTAS = Channel.empty()
}else if(params.format == 'fasta'){
  SEQS.into{FORWINDOWS}
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
  file '*.fasta' into FORWINDOWS

  """
  ${workflow.projectDir}/gbk2fasta.py \
    --input $gbk_file
  """
}

process split_in_windows{
  label 'py3'

  input:
  file ref from FORWINDOWS
  val w_size from params.w_size
  val s_size from params.s_size

  output:
  set val("$ref"), file("$ref"), file("*_windows.fasta") into WINDOWS

  """
  ${workflow.projectDir}/split_seq_into_windows.py \
    --input $ref \
    --format fasta \
    --w_size $w_size \
    --s_size $s_size
  """
}

process bowtie2{
  label 'bowtie2'
  publishDir "${params.outdir}/bam/", mode: 'rellink'

  input:
  set ref, file(fasta), file(windows) from WINDOWS

  output:
  set val("$ref"), file("${ref}.bam"), file(fasta) into BOWTIE2_RES, BOWTIE_RES_FOR_PLOT

  script:
  if( params.sensitivity == 'sensitive' ){
    """
    bowtie2-build $fasta $fasta
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
  }else if(params.sensitivity == 'very-sensitive'){
    """
    bowtie2-build $fasta $fasta
    bowtie2 \
      -f \
      --end-to-end \
      --very-sensitive \
      -a \
      --time \
      --threads 1 \
      -x $ref \
      -U $windows | \
      samtools view -b - > ${ref}.bam
    """
    else
      error "Invalid sensitivity argument ($params.sensitivity)"
  }

}

process merge_bam_windows{
  label 'py3'
  publishDir "${params.outdir}/pGRINS.gff3/", mode: 'rellink'

  input:
  set ref, file(bam), file(fasta) from BOWTIE2_RES

  output:
  set val("$ref"), file("${ref}.pgrins.gff3"), file(fasta) into GFF3, GFF3_FOR_PLOT

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
  set ref, file("${ref}.pgrins.gff3"), file(fasta) from GFF3

  output:
  set ref, file("${ref}.pgrins.fasta") into PGRINS_FASTA

  """
  ${workflow.projectDir}/gff_to_fasta.py \
    --input $fasta \
    --gff3 ${ref}.pgrins.gff3 \
    --output ${ref}.pgrins.fasta \
    --format fasta
  """
}

process vsearch_pgrins{
  label 'vsearch'
  publishDir "${params.outdir}/vsearch.pgrins/", mode: 'rellink'

  input:
  set ref, file("${ref}.pgrins.fasta") from PGRINS_FASTA

  output:
  set ref, file("${ref}.pgrins.centroids.fasta"), file("${ref}.clusters.uc") into PGRINS_UC

  """
  vsearch \
    --threads 1 \
    --sizeout \
    --cluster_fast ${ref}.pgrins.fasta \
    --strand both \
    --centroids ${ref}.pgrins.centroids.fasta \
    --relabel centroid_ \
    --sizeorder \
    --uc ${ref}.clusters.uc \
    --id $params.vsearch_id \
    --maxseqlength 100000
  """

}

// println "==========================="
// GFF3_FOR_PLOT.join(PGRINS_UC).join(BOWTIE_RES_FOR_PLOT).subscribe{println it}

// FORPLOTS.cross(GFF3_FOR_PLOT).cross(PGRINS_UC).cross(BOWTIE_RES_FOR_PLOT).subscribe{println it}

process plot_fast_grins{
  label 'py3'
  publishDir "${params.outdir}/plots/", mode: 'rellink'

  input:
  set ref, file(grins_gff), file(fasta_ref1), file(grins_centroids), file(grins_clusters), file(windows_bam), file('temp') from GFF3_FOR_PLOT.join(PGRINS_UC).join(BOWTIE_RES_FOR_PLOT)
  // set ref2, file("${ref2}.pgrins.gff3"), file(fasta2) from GFF3_FOR_PLOT
  // set ref3, file("${ref3}.pgrins.centroids.fasta"), file("${ref3}.clusters.uc") from PGRINS_UC
  // set ref4, file("${ref4}.bam"), file(ref) from BOWTIE_RES_FOR_PLOT

  output:
  file "${ref}.png" into PLOTS

  """
  ${workflow.projectDir}/plot_fast_grins.py \
    --input $fasta_ref1 \
    --grins_gff3 $grins_gff \
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
    module = 'bowtie2:samtools'
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
*/
