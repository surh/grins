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

// Nextflow pipeline to run phylogeny in all alingnments in folder_name
// It will look for all files in the input directory with extension `.aln`
// and call the specified raxml executable (default is --bin raxmlHPC).
// If an executable that allows for threads is passed, then the number of
// threads can be specified with --threads <n>. Nextflow assumes that
// there is a module called raxml that can be loaded with
// module load raxml, and that by default loads version 8.2.X.

// For each allignment, the RAxML command  will first conduct a bootstrap
// search and once that is done a search for the best–scoring
// maximum likelihood tree.

// It will return the bootstrapped trees, the best scoring maximum likelihood
// tree, and the bootstrap support values drawn on the best-scoring tree
// as node labels and branch labels.

// By increasing the number of bootstrap  replicates (--bootstrap parameter)
// the maximum lilkelihood will also be more thourough. The option
// --bootstrap autoMRE can also be passed and RAxML will execute a maximum
// of 1000 bootstrap though it might converge earlier.

// Example run:
// nextflow run raxml_all_nuc.nf --bootstrap autoMRE --indir myalns/ --seed 12345

// Parameters
params.indir = './'
params.bin = 'raxmlHPC'
params.seed = 12345
params.bootstrap = 100
params.threads = 1
// params.bootstrap = 'autoMRE'

// Process inputs
ALNS = Channel.
  fromPath("${params.indir}/*.aln").
  map{file -> tuple(file.baseName, file)}

process raxml{
  cpus params.threads
  maxForks 10
  module 'raxml'
  time 100.h
  memory 2.G

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
    -T ${params.threads}
    -n $filename
  """
}
