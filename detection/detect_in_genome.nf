#!/usr/bin/env nextflow
// Copyright (C) 2019-2020 Sur Herrera Paredes

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
/*
--indir
Directory containing genomes. One file per genome with extension .fa or
.fasta.

--outdir
Name for output directory. [output]

--sensitivity
Alignment sensitivity parameter for bowtie2. Either 'sensitive' or
'very-sensitive'. [sensitive]

--w_size
Size in bp for the sliding window to search for duplications. [150]

--s_size
Step size in bp for the sliding windows to search for duplications. [30]

--min_size
Minimum size of duplicated regions to report in bp. [500]

--antismash_cpus
Number of threads to use for antiSMASH 5 predictions. [8]
*/

// Params
params.indir = 'genomes/'
params.outdir = 'output/'
params.sensitivity = 'sensitive'
params.w_size = 150
params.s_size = 30
params.min_size = 500
params.antismash_cpus = 8



// Process params
GENOMEFA = Channel.fromPath("${params.indir}/*", type: 'file')
  .map{genomefa -> tuple(genomefa.name.replaceAll(/\.fa(sta)?$/, ""),
    file(genomefa))}
GENOMEFA.into{GENOME4DUPS; GENOME4ANTISMASH}

process antismash5{
  label 'antismash5'
  cpus params.antismash_cpus
  tag "$acc"
  publishDir "$params.outdir/antismash", mode: 'rellink',
    saveAs: {"$acc"}

  input:
  tuple acc, file(genome_file) from GENOME4ANTISMASH
  val threads from params.antismash_cpus

  output:
  tuple acc, file("output") into ANTISMASH

  """
  antismash \
    -c $threads \
    --taxon bacteria \
    --output-dir output/ \
    --verbose \
    --genefinding-tool prodigal \
    $genome_file
  """
}

process antismash2gff3{
  label 'py3'
  tag "$acc"
  publishDir "${params.outdir}/antismash.gff3/", mode: 'rellink'

  input:
  tuple acc, file('antismash') from ANTISMASH

  output:
  tuple acc, file("${acc}.gff3") into ASMASHREGIONS

  """
  ${workflow.projectDir}/antismash2gff3.py \
    --input antismash/${acc}.json \
    --output ${acc}.gff3
  """
}

process split_in_windows{
  label 'py3'
  tag "$acc"

  input:
  tuple acc, file(ref) from GENOME4DUPS
  val w_size from params.w_size
  val s_size from params.s_size

  output:
  tuple acc, file("$ref"), file("*_windows.fasta") into WINDOWS

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
  tag "$acc"

  input:
  tuple acc, file(fasta), file(windows) from WINDOWS

  output:
  tuple acc, file("${acc}.bam"), file(fasta) into BOWTIE2_RES, BOWTIE_RES_FOR_PLOT

  script:
  if( params.sensitivity == 'sensitive' )
    """
    bowtie2-build $fasta $fasta
    bowtie2 \
      -f \
      --end-to-end \
      --sensitive \
      -a \
      --time \
      --threads 1 \
      -x $fasta \
      -U $windows | \
      samtools view -b - > ${acc}.bam
    """
  else if(params.sensitivity == 'very-sensitive')
    """
    bowtie2-build $fasta $fasta
    bowtie2 \
      -f \
      --end-to-end \
      --very-sensitive \
      -a \
      --time \
      --threads 1 \
      -x $fasta \
      -U $windows | \
      samtools view -b - > ${acc}.bam
    """
    else
      error "Invalid sensitivity argument ($params.sensitivity)"
}

process merge_bam_windows{
  label 'py3'
  publishDir "${params.outdir}/duplicated.gff3/", mode: 'rellink'
  tag "$acc"

  input:
  tuple acc, file(bam), file(fasta) from BOWTIE2_RES
  val min_size from params.min_size
  val w_size from params.w_size

  output:
  tuple acc, file("${acc}.duplicated.gff3") into DUPGFF3

  """
  ${workflow.projectDir}/produce_windows_from_bam.py \
    --input $bam \
    --output ${acc}.duplicated.gff3 \
    --w_size $w_size \
    --min_size $min_size
  """
}

process intersect_asmash_dups{
  label 'bedtools'
  publishDir "${params.outdir}/bgcdups.gff3", mode: 'rellink'
  tag "$acc"

  input:
  tuple acc, file("duplicated.gff3"),
    file('regions.gff3') from DUPGFF3.join(ASMASHREGIONS)

  output:
  tuple acc, file("${acc}.bgcdups.gff3") into BGCDUPS

  """
  bedtools intersect -wo -a duplicated.gff3 -b regions.gff3 | \
    cut -f 1-9,18 | \
    sed 's/ID=region/Region=region/' | \
    perl -e 'while(<>){
      chomp;
      @line = split(/\\t/, \$_);
      \$line[8] = "\$line[8];\$line[9]";
      pop @line;
      print join("\\t", @line) . "\\n"
    }' > ${acc}.bgcdups.gff3
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
