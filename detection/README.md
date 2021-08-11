# Detecting GRINS
This folder contains the scripts for detecting GRINS

## Key files

* **antismash2gff3.py**: Converts antiSMASH 5 BGC predictions to GFF3 format.

* **detect_in_genome.nf**: Pipeline to detect duplicated regions in
bacterial genomes in fasta format. Includes annotations via antiSMASH 5.

* **GRINS_detection_from_BOWTIE.py**: Takes duplicated regions identified by
`detect_in_genome.nf` and identifies GRINS regions.

* **GRINS_detection_from_BOWTIE.nf**: Pipeline that parallelizes the python
script of the same name.

* **GRINS_detection_pairwise_identity.py**: Script to detect GRINS within
PKS clusters.

* **pairwise_alignment_windows.nf**: Pipeline that parallelizes the python
script of the same name.

* **pairwise_alignment_windows.py**: Script that splits PKS cluster sequences
into windows and performs pairwise alignments which can be used to identify
duplicated regions and GRINS.

* **produce_windows_from_bam.py**: Reads bowtie2 BAM output and identifies
duplicated regions.

* **split_seq_into_windows.py**: Creates sliding windows for all analysis.


## Dependencies

The following software is required to run the current GRINS detection workflow.
The versions specified are those that were tested. It is possbile that the
workflow will run with similar versions but it has not been tested.

The recommended setup is to have conda installed as an
[Lmod](https://lmod.readthedocs.io/en/latest/index.html) module, and then
install dependencies as either other Lmod modules or conda environments.

* **antiSMASH (5.1.1 or 5.1.2)**: We recommend installing antiSMASH to its
own conda environment. See [their documentation](https://docs.antismash.secondarymetabolites.org/install/).
See also [our conda environment .yml file](../conda_envs/antismash.yml).

* **Python (3.6.13)**: Python 3 is required including the following modules:
argparse, biopython, json, matplotlib, numpy, os, and pybedtools. We recommend
to have a conda environment with these dependencies.
See our [conda environment .yml file](../conda_envs/fraserconda.yml)
for precise versions. Please note that this environment include many additional
Python 3 modules which are not required.

* **bowtie2 (2.4.1)**: Check [their website](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for details.
We recommend compiling locally and installing bowtie2 as an Lmod module.

* **samtools (1.7 or 1.10)**: Check [their website](http://www.htslib.org/) for details.
We recommend compiling locally and installing samtools as an Lmod module.

* **bedtools (2.30)**: Check [their website](https://bedtools.readthedocs.io/en/latest/index.html) for details.
Since pybedtools requires this software, you might have installed it when you
created the Python 3 conda environment. If not, we recommend compiling locally
and installing asn an Lmod module.

* **nextflow (20.10.0 build 5430)**: Check [their website](https://www.nextflow.io/) for details.
We recommend placing the nextflow executable directly in your PATH.

## Nextflow configuration

The current pipeline is divided in two nextflow sub-pipelines. Depending on
the type of installations you might need to use a `nextflow.config` file
to indicate where the different dependencies are located.

See the [nextflow documentation](https://www.nextflow.io/docs/latest/index.html)
for a full explanation of details and options. At the end of
[each](https://github.com/surh/grins/blob/ad4259bd86a5c587bb2e5ede7ba895192828b8f3/detection/detect_in_genome.nf#L202-L231)
[nextflow](https://github.com/surh/grins/blob/005d9dae4a3a619f2aecaee1bb24b81ea869883b/detection/GRINS_detection_from_BOWTIE.nf#L59-L88)
(**.nf**) script there is a commented section titled `Example nextflow.config`
which shows the `nextflow.config` file that we utilized for each step.

In brief, we use nextlow labels to indicate to each process which dependencies
it needs and where to find them. So if some dependencies are not readily
available, your nextflow configuration must include something like:

```nextflow
process{
  // ...
  withLabel: 'py3'{
    // Tell nextflow to load the Python 3 environment with the required modules.
  }
  withLabel: 'bedtools'{
    // Tell nextflow to load bedtools
  }
  withLabel: 'bowtie2'{
    // Tell nextflow to load bowtie2 & samtools
  }
  withLabel: 'antismash5'{
    // Tell nextflow to load antiSMASH 5
  }
}
```

## Workflow example

Below is an example of how to run the pipeline for a single genome. We will
use as example the RefSeq representative genome of *Mycobacterium leprae*
(GCF_003253775.1_ASM325377).

You will need to place the fasta genome assembly file of this  genome
(and all genomes you want to process) in a dedicated directory. The fasta
files need to have `.fasta` or `.fa` extensions to be recognized.

### Getting example genome
The following comands will: 1) create a new directory `genomes_fasta` for
our input files, 2) download the compressed fasta assembly into the directory
above, 3) uncompress the genome, 4) change the file extension so that it
works with the pipeline:

```bash
mkdir genomes_fasta
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/253/775/GCF_003253775.1_ASM325377v1/GCF_003253775.1_ASM325377v1_genomic.fna.gz -P genomes_fasta/
gzip -d genomes_fasta/GCF_003253775.1_ASM325377v1_genomic.fna.gz
rename 's/.fna/.fasta/' genomes_fasta/*.fna
```

That last line is overkill for a single file, but it can be useful to change
the extension of large numbers of files with one command.

You can list the contents of the newly created directory with
`ls -l genomes_fasta`, and you should see something like:

```
total 3152
drwxrwsr-x 2 sur science      62 ago 11 14:08 ./
drwxrwsr-x 3 sur science      34 ago 11 14:08 ../
-rw-rw-r-- 1 sur science 3227034 jun 22  2018 GCF_003253775.1_ASM325377v1_genomic.fasta
```

### Downloading the code for predicting GRINS

If you haven't, download the code with the following command.
Feel free to use your path of choice.

```bash
git clone git@github.com:surh/grins /path/to/grins/repo
```

### Annotating genomes and finding duplicated regions

The first step of the pipeline is running the `detect_in_genome.nf` nextflow
script. The full set of parameters is
[documented in the script itself](https://github.com/surh/grins/blob/10753e80b0bfca2a26e6eb7c19824e82ef74b2ac/detection/detect_in_genome.nf#L17-L41).
But you can start with the default parameters by running the following command.
Make sure you use the path where you downloaded the repository.

```bash
nextflow run /path/to/grins/repo/detection/detect_in_genome.nf --indir genomes_fasta --outdir first_step --antismash_cpus 2
```

This step can take a while depending on your system, specific genome, and
number of antiSMASH 5 cpus. When it finishes, you should see something like:

```
N E X T F L O W  ~  version 20.10.0
Launching `/home/sur/micropopgen/src/grins/detection/detect_in_genome.nf` [exotic_linnaeus] - revision: 64f432c5ae
executor >  slurm (6)
[c0/881fae] process > antismash5 (GCF_003253775.1_ASM325377v1_genomic)    [100%] 1 of 1 ✔
[0c/654fde] process > antismash2gff3 (GCF_003253775.1_ASM325377v1_geno... [100%] 1 of 1 ✔
[2f/1ddd03] process > split_in_windows (GCF_003253775.1_ASM325377v1_ge... [100%] 1 of 1 ✔
[53/a2c15a] process > bowtie2 (GCF_003253775.1_ASM325377v1_genomic)       [100%] 1 of 1 ✔
[60/970e77] process > merge_bam_windows (GCF_003253775.1_ASM325377v1_g... [100%] 1 of 1 ✔
[df/7d7ec4] process > intersect_asmash_dups (GCF_003253775.1_ASM325377... [100%] 1 of 1 ✔
Completed at: 11-Aug-2021 14:40:23
Duration    : 20m 2s
CPU hours   : 0.6
Succeeded   : 6
```
And you should be able to see the results in the newly created directory with
`ls -l first_step/`, which should show something like:

```
total 0                                                                                              │
drwxrwsr-x 2 sur science  56 ago 11 14:25 antismash                                                  │
drwxrwsr-x 2 sur science  61 ago 11 14:31 antismash.gff3                                             │
drwxrwsr-x 2 sur science 112 ago 11 14:40 bam                                                        │
drwxrwsr-x 2 sur science  69 ago 11 14:40 bgcdups.gff3                                               │
drwxrwsr-x 2 sur science  72 ago 11 14:40 duplicated.gff3   
```

### Identifying GRINS

The second step of the pipeline consists in running the nextflow script
`GRINS_detection_from_BOWTIE.nf`.

For running this second part we must first create a directory with the
antiSMASH-annotated genomes in genbank format. This is easy to do for
one genome, but can become unwieldy as their number increases. The following
bash commands will create a directory with all the genomes from the previous
step. It must be run from the **same directory** where you have performed all
previous operations

```bash
mkdir asmash_genbanks
cd asmash_genbanks/
find ../first_step/antismash/ -mindepth 1 -maxdepth 1 -exec basename {} \; | while read g; do ln -s ../first_step/antismash/$g/$g.gbk $g.gbk; done
cd ../
```
You can check the results of the previous commands with `ls -l asmash_genbanks`,
and you should see something like:

```
total 0
lrwxrwxrwx 1 sur science 99 ago 11 14:53 GCF_003253775.1_ASM325377v1_genomic.gbk -> ../first_step/antismash/GCF_003253775.1_ASM325377v1_genomic/GCF_003253775.1_ASM325377v1_genomic.gbk
```
Now you are ready to run the second part of the pipeline. There are only three
options and you can run it with the command below.
Make sure you use the path where you downloaded the repository.

```bash
nextflow run /path/to/grins/repo/detection/GRINS_detection_from_BOWTIE.nf --gbk asmash_genbanks/ --gff_dups first_step/duplicated.gff3/ --outdir second_step
```

When it is done you should see something like:

```
N E X T F L O W  ~  version 20.10.0
Launching `/home/sur/micropopgen/src/grins/detection/GRINS_detection_from_BOWTIE.nf` [focused_payne] - revision: d75e29b08b
executor >  slurm (1)
[db/3dc949] process > grinspred (GCF_003253775.1_ASM325377v1_genomic) [100%] 1 of 1 ✔
```

### Understanding the output
