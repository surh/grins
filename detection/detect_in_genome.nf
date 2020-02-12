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
*/

// Params
params.indir = 'genomes/'
params.outdir = 'output/'
// params.antismash = "antismash/"
// params.format = 'fasta'
params.sensitivity = 'sensitive'
params.w_size = 150
params.s_size = 30
// params.vsearch_id = 0.9
params.min_size = 500
params.antismash_cpus = 8



// Process params
GENOMEFA = Channel.fromPath("${params.indir}/*", type: 'file')
  .map{genomefa -> tuple(genomefa.name.replaceAll(/\.fa(sta)?$/, ""),
    file(genomefa))}
// ANTISMASH = Channel.fromPath("${params.antismash}/*/txt/*_BGC.txt",
//   type: 'file')
//   .map{bgcpreds -> tuple(bgcpreds.getParent().getParent().name,
//     bgcpreds.name.replaceAll(/_BGC\.txt$/,''),
//     file(bgcpreds))}


// ANTISMASH.join(GENOMEFA).subscribe{println it}
// GENOMEFA.cross(ANTISMASH).flatten().collate(5).subscribe{println it}
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
  publishDir "${params.outdir}/antismash.regions/", mode: 'rellink'

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
  tag "$genome"

  input:
  tuple genome, file(ref) from GENOME4DUPS
  val w_size from params.w_size
  val s_size from params.s_size

  output:
  tuple genome, file("$ref"), file("*_windows.fasta") into WINDOWS

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
  tag "$genome"

  input:
  tuple genome, file(fasta), file(windows) from WINDOWS

  output:
  tuple genome, file("${genome}.bam"), file(fasta) into BOWTIE2_RES, BOWTIE_RES_FOR_PLOT

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
      samtools view -b - > ${genome}.bam
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
      samtools view -b - > ${genome}.bam
    """
    else
      error "Invalid sensitivity argument ($params.sensitivity)"
}

process merge_bam_windows{
  label 'py3'
  publishDir "${params.outdir}/duplicated.gff3/", mode: 'rellink'
  tag "$genome"

  input:
  tuple genome, file(bam), file(fasta) from BOWTIE2_RES
  val min_size from params.min_size

  output:
  // set genome, file("${genome}.pgrins.gff3"), file(fasta) into GFF3, GFF3_FOR_PLOT
  tuple genome, file("${genome}.pgrins.gff3") into GFF3, GFF3_FOR_PLOT

  """
  ${workflow.projectDir}/produce_windows_from_bam.py \
    --input $bam \
    --output ${genome}.pgrins.gff3 \
    --w_size ${params.w_size} \
    --min_size $min_size
  """
}

// process intersect{
//   label 'py3'
//   tag "$genome-$record"
//   publishDir "${params.outdir}/bgcfastas", mode: 'rellink'
//
//   input:
//   tuple genome, file(genomefa), genome2, record, file(bgcpreds) from GENOMEFA.cross(ANTISMASH).flatten().collate(5)
//
//   output:
//   tuple genome, file("${record}_bgcs.fasta") into BGCFASTAS
//
//   """
//   # Convert to BED
//   cut -f 1,4 $bgcpreds | \
//     grep -v BGC_range | \
//     sed 's/;/\\t/' | \
//     awk '{print "$record\\t" \$2 "\\t" \$3 "\\t" \$1}' > ${record}.bed
//
//   # Remove version from accession ID
//   cat $genomefa | \
//     perl -e 'while(<>){chomp; \
//       if(\$_ =~ /^>/){ \
//         (\$id, @head) = split(/\\s/, \$_); \
//         \$id =~ s/>//; \
//         \$id =~ s/\\.[\\d]+\$//; \
//         print ">\$id\\n" \
//       }else{print "\$_\\n"}}' > versionless.fa
//
//   # Create index
//   samtools faidx versionless.fa
//
//   # Intersect bed and get fasta2
//   bedtools getfasta \
//     -fi versionless.fa \
//     -bed ${record}.bed \
//     -name > ${record}_bgcs.fasta
//   """
// }
//
// process cat_records{
//   tag "$genome"
//   publishDir "${params.outdir}/genomebgcs", mode: 'rellink'
//
//   input:
//   tuple genome, file(records) from BGCFASTAS.groupTuple()
//
//   output:
//   tuple genome, file("${genome}_bgcs.fasta") into GENOMEBGCS
//
//   """
//   cat $records > ${genome}_bgcs.fasta
//   """
//
// }



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
  withLabel: 'bowtie2'{
    module = 'bowtie2:samtools'
  }
  withLabel: 'antismash5'{
    module = 'anaconda'
    conda = '/opt/modules/pkgs/anaconda/3.6/envs/antismash5'
  }
}

executor{
  name = 'slurm'
  queueSize = 500
  submitRateLitmit = '1 sec'
}
*/
