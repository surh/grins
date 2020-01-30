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


// ANTISMASH.join(GENOMEFA).subscribe{println it}
// GENOMEFA.cross(ANTISMASH).flatten().collate(5).subscribe{println it}


process intersect{
  label 'py3'
  tag "$genome-$record"
  publishDir params.outdir, mode: 'rellink'

  input:
  tuple genome, file(genomefa), genome2, record, file(bgcpreds) from GENOMEFA.cross(ANTISMASH).flatten().collate(5)

  output:
  tuple genome, file("${record}_bgcs.fasta") into BGCFASTAS

  """
  # Convert to BED
  cut -f 1,4 $bgcpreds | \
    grep -v BGC_range | \
    sed 's/;/\\t/' | \
    awk '{print "$record\\t" \$2 "\\t" \$3 "\\t" \$1}' > ${record}.bed

  # Remove version from accession ID
  cat $genomefa | \
    perl -e 'while(<>){chomp; \
      if(\$_ =~ /^>/){ \
        (\$id, @head) = split(/\\s/, \$_); \
        \$id =~ s/>//; \
        \$id =~ s/\\.[\\d]+\$//; \
        print ">\$id\\n" \
      }else{print "\$_\\n"}}' > versionless.fa

  # Create index
  samtools faidx versionless.fa

  # Intersect bed and get fasta2
  bedtools getfasta \
    -fi versionless.fa \
    -bed ${record}.bed \
    -name > ${record}_bgcs.fasta
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
