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

// Parameters
params.indir = './'
params.format = 'genbank'
params.outdir = './output'
params.w_size = 150
params.s_size = 30


INPUTS = Channle.fromPath("${params.indir}/*")

process pairwise{
  label 'py3'
  publishDir params.outdir, mode: 'rellink'

  input:
  file sequence from INPUTS

  output:
  file "*.txt"

  """
  ${workflow.projectDir}.pairwise_alignment_windows.py \
    --input $sequence \
    --output ./ \
    --w_size ${params.w_size} \
    --s_size ${params.s_size} \
    --format ${params.format}
  """
}

// nextflow.config example
/*
process{
  maxForks = 300
  queue = 'hbfraser,hns'
  memory = '5G'
  time = 2d
  withLabel: 'py3'{
    module = 'conda'
    conda = '/home/groups/hbfraser/modules/packages/conda/4.6.14/envs/fraserconda'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  sumitRateLimit = '1 sec'
}
*/
