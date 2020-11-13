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
  python ${workflow.projectDir}/GRINS_detection_from_BOWTIE.py \
    --seq_input gbk/ \
    --dupl_input dup/ \
    --with_plots yes
  """
}

// Example nextflow.config
/*
process{
  queue = 'hbfraser,hns'
  maxFors = 10
  errorStrategy = 'finish'
  stageInMode = 'rellink'
  withLabel: 'py3'{
    conda = '/opt/modules/pkgs/anaconda/4.8/envs/fraserconda'
  }
}

executor{
  name = 'slurm'
  queueSize = 50
  submitRateLitmit = '1 sec'
}
*/
