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
params.indir = 'genomes/'
params.outdir = 'output/'
params.antismash = "antismash/"
params.format = 'fasta'
params.sensitivity = 'sensitive'
params.w_size = 150
params.s_size = 30
params.vsearch_id = 0.9




// Process params
GENOMEFA = Channel.fromPath("${params.indir}/*", type: 'file')
  .map{genomefa -> tuple(genomefa.name.replaceAll(/\.fa(sta)?$/, ""),
    file(genomefa))}
ANTISMASH = Channel.fromPath("${params.antismash}/*/txt/*_BGC.txt",
  type: 'file')
  .map{bgcpreds -> tuple(bgcpreds.getParent().getParent().name,
    bgcpreds.name.replaceAll(/_BGC\.txt$/,''),
    file(bgcpreds))}


// GENOMEFA.join(ANTISMASH).subscribe{println it}
process intersect{
  label 'py3'


  input:
  tuple genome, file(genomefa), record, file(bgcpreds) from GENOMEFA.join(ANTISMASH)

  """
  # Convert to BED
  cut -f 4 $bgcpreds |  \
    grep -v BGC_range | \
    sed 's/;/\\t/' | \
    awk '{print "$record\\t" \$0}' > ${record}.bed

  # Intersect bed and get fasta2
  bedtools getfasta -fi $genomefa -bed ${record.bed} > ${record}_bgcs.fasta

  """



}


// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxFors = 300
  errorStrategy = 'finish'
  withLabel: 'py3'{
    module = 'fraserconda'
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/fraserconda'
  }
}
executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
