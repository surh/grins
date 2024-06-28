#!/usr/bin/env nextflow

// Copyright (C) 2020-2021 Sur Herrera Paredes

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

// params
params.gbk = ''
params.gff_dups = ''
params.outdir = "output/"

//  Process parameters
gbk_dir = file(params.gbk)
gff_dups_dir = file(params.gff_dups)

GBKS = Channel.fromPath("${gbk_dir}/*", type: 'file')
  .map{gbk_file -> tuple(gbk_file.name.replaceAll(/\.gbk$/, ''),
    file(gbk_file))}

GFFS = Channel.fromPath("${gff_dups_dir}/*", type: 'file')
  .map{dup_file -> tuple(dup_file.name.replaceAll(/\.duplicated.gff3$/, ''),
    file(dup_file))}

process grinspred{
  tag "$genome"
  label "py3"
  publishDir "$params.outdir/outputs", mode: 'rellink',
    pattern:"output", saveAs: {"$genome"}
  publishDir "$params.outdir/res", mode: 'rellink',
    pattern:"GRINS_detected_in_genomes_and_BGCs.txt",
    saveAs: {"${genome}.txt"}

  input:
  tuple genome, file("gbk/${genome}.gbk"), file("dup/${genome}.duplicated.gff3") from GBKS.join(GFFS)

  output:
  tuple file("output"), file("GRINS_detected_in_genomes_and_BGCs.txt")

  """
  mkdir output
  python3 ${workflow.projectDir}/GRINS_detection_from_BOWTIE.py \
    --seq_input gbk/ \
    --dupl_input dup/ \
    --with_plots yes
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxFors = 300
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  withLabel: 'py3'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
  }
  withLabel: 'bedtools'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
  }
  withLabel: 'bowtie2'{
    module = 'bowtie2:samtools'
  }
  withLabel: 'antismash5'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/antismash'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
