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
[nextflow](https://github.com/surh/grins/blob/ad4259bd86a5c587bb2e5ede7ba895192828b8f3/detection/GRINS_detection_from_BOWTIE.nf#L42-L59)
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

Now you are ready to run the first step of the pipeline. If you haven't,
download the code with the following command. Feel free to use your path
of choice.

```bash
git clone git@github.com:surh/grins /path/to/grins/repo
```

The first step of the pipeline is running the `detect_in_genome.nf` nextflow
script. The full set of parameters is
[documented in the script itself](https://github.com/surh/grins/blob/10753e80b0bfca2a26e6eb7c19824e82ef74b2ac/detection/detect_in_genome.nf#L17-L41).
But you can start with the default parameters by running:

```bash
nextflow run /path/to/grins/repo/detection/detect_in_genome.nf --indir genomes_fasta --outdir first_step
```
