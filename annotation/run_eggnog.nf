#!/usr/bin/env nextflow
// Copyright (C) 2018-2019 Sur Herrera Paredes

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

// Params
params.indir = "FAA/"
params.outdir = 'annots'

CDS = Channel.fromPath("${params.indir}/*.faa")
  .map{gbk_file -> tuple(gbk_file.name.replaceAll(/\.faa/, ""),
  file(gbk_file))}


process eggnog{
  label 'eggnog'
  tag "$genome"
  publishDir params.outdir, mode: 'rellink'

  input:
  tuple genome, file(faa_file) from CDS

  output:
  file "${genome}.emapper.annotations"

  script:
  """
  emapper.py --database bact \
    --data_dir /opt/pkgs/eggnog/1.0.3/data/ \
    --output_dir ./ \
    -i ${faa_file} \
    --cpu 1 \
    --output ${genome}
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxForks = 40
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  time = '48h'
  memory = '1G'
  withLabel: 'eggnog'{
    module = 'eggnog:anaconda'
    conda = "/opt/modules/pkgs/anaconda/3.6/envs/python2"
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
