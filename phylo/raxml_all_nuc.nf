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

// Nextflow pipeline that takes a matrix of SNP frequencies, a file
// with phenotypes and a matrix of covariates and performs association
// for all SNPs and testing via permutation.

// Script to run phylogeny in a alingnments in folder_name

// Parameters
params.indir = './'
params.bin = 'raxmlHPC'
params.seed = 12345
params.bootstrap = 100
// params.bootstrap = 'autoMRE'

// Process inputs
ALNS = Channel.
  fromPath("${params.indir}/*.aln").
  map{file -> tuple(file.baseName, file)}

process raxml{
  cpus 1
  maxForks 10
  module 'raxml'
  time '100:00:00'

  input:
  set filename, file(aln) from ALNS

  output:
  file "RAxML_bestTree.${filename}"
  file "RAxML_bipartitionsBranchLabels.${filename}"
  file "RAxML_bipartitions.${filename}"
  file "RAxML_bootstrap.${filename}"

  """
  ${params.bin} \
    -s $aln \
    -f a \
    -x ${params.seed} \
    -p ${params.seed} \
    -# ${params.bootstrap} \
    -m GTRGAMMA \
    -n $filename
  """
}
