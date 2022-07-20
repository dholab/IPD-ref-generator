# Immuno Polymorphism Database - Reference Sequence Generator

_Latest Release:_ 3.9.0.0 (2022-07) build 209 (July 2022)

## Overview

This workflow pulls the latest allele and protein sequences from [IPD](https://www.ebi.ac.uk/ipd/) for the MHC and KIR regions. It also prepares a series of species-specific reference databases. These databases will contain sequences that are more than 100 base pairs long for Rhesus macaque (_Macaca mulatta_, a.k.a. Mamu), Cynomolgus macaque (_Macaca fascicularis_ a.k.a. Mafa), and Southern pig-tailed macaque (_Macaca nemestrina_, a.k.a. Mame), as well as for all non-human primates included in IPD. At minimum, this workflow should be run once after each IPD release of new MHC reference alleles.

The workflow also prepares context-specific reference databases for use with data generated at [AVRL](https://dholk.primate.wisc.edu/project/home/begin.view?):

1. Reference allele sequences for use with high depth-of-coverage whole exome sequence data enriched for immuno-genes (iWES)
2. Reference allele sequences for use with Illumina MiSeq amplicon libraries

## Getting Started

To run this workflow, simply `git clone` it into your working directory of choice, like so:

```
git clone https://github.com/nrminor/IPD-ref-generator.git .
```

Once the workflow bundle is in place, first ensure that the workflow scripts are executable, like so:

```
chmod +x bin/*.py
```

Next, build the Docker image that contains the workflow's dependencies:

```
docker build --tag ipd-ref-generator:v1_0_6 config/
```

Note that to build the above docker container, you may need to increase the amount of memory allotted to Docker in the Docker Engine preferences.

### Nextflow Installation

This workflow uses the [NextFlow](https://www.nextflow.io/) workflow manager. We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are set and ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in the BASH terminal:

```
nextflow run ipd-ref-generator.nf
```

If the workflow runs partway, but a computer outage or other issue interrupts its progress, no need to start over! Instead, run:

```
nextflow run ipd-ref-generator.nf -resume
```

The workflow's configurations (see below) tell NextFlow to plot the workflow and record run statistics. However, the plot the workflow, note that NextFlow requires the package GraphViz, which is easiest to install via the intructions on [GraphViz's website](https://graphviz.org/download/).

## Configuration

The following runtime parameters have been set for the whole workflow and are available for modification:

- `pull_mhc`, `pull_kir`, `pull_mhc_proteins`, `pull_kir_proteins` - These settings must be set to either _true_ or _false_, and tell the workflow which datasets to download. By default, in the `nextflow.config` file included in this repo, the workflow downloads all four datasets: all non-human primate MHC sequences, all non-human KIR sequences, all non-human primate MHC protein sequences, and all non-human KIR protein sequences. To download only one of these datasets, set that parameter to _true_ and the others to _false_.
- `mhc_allele_count`, `kir_allele_count`, `mhc_protein_count`, `kir_protein_count` - These are the most important settings for the workflow to generate up-to-date databases. To determine how many alleles the workflow should gather, scroll to the bottom right of the MHC allele list at [https://www.ebi.ac.uk/ipd/mhc/allele/list/](https://www.ebi.ac.uk/ipd/mhc/allele/list/) or the KIR allele list at [https://www.ebi.ac.uk/ipd/kir/alleles/](https://www.ebi.ac.uk/ipd/kir/alleles/)
  - Future versions of this workflow will allow users to download MHC alleles _and_ KIR alleles, though at this time only MHC alleles can be pulled.
- `iwes_exemplar` - This specifies the path to an exemplar for iWES data, which in most cases need not be changed.
- `miseq_exemplar` - This specifies the path to an exemplar for MiSeq data, which in most cases need not be changed.
- `results` - path to where the workflow's outputs should be placed. The default is `results/` within the workflow directory.

These parameters can be altered in the command line with a double-dash flag, like so:

```
nextflow run ipd-ref-generator.nf --results ~/ipd_results/
```

### Bundled Data

This workflow requires minimal pre-bundled data, as mentioned above: the iWES exemplar (Mamu-exon2-exemplar.fasta) and the MiSeq exemplar (Mamu_MiSeq_representative_alleles.fasta), both of which are stored in `resources/`.

## Workflow Summary

- First, the workflow downloads the full list of non-human primate MHC alleles from EMBL, as specified with the above mentioned `allele_count` parameter in `nextflow.config`. This takes place in the process PULL_IPD
- Next, each file generated in PULL_IPD is cleaned in the process CLEAN_IPD, which removes "X's" that IPD places at the end of each sequence. The process also filters out alleles that are less than 100 base pairs long.
- Each cleaned, genbank-formatted reference file is then trimmed in 2 ways independently and in parallel: Once for iWES data in the process IWES_TRIMMING, and once for MiSeq amplicon data in the process MISEQ_TRIMMING.

### A note on parallel computation

NextFlow automatically allocates cores to each process, as if it is a job on a compute node. Depending on how many cores you make available to this workflow, the four databases it outputs can each be downloaded, cleaned, and trimmed at the same time as one another.

### Output Files

For the following files, there should be 4 of each: one for mamu, mafa, mane, and nhp. In total, this means there will be 20 output files.

- `ipd-{locus}-{animal acronym}-{run date}_cleaned.gbk` is a genbank-format list of all downloaded alleles for each animal that are longer than 100 bases.
- `ipd-{locus}-{animal acronym}-{run date}_cleaned.miseq.trimmed.deduplicated.fasta` is FASTA-formatted list of trimmed alleles with no duplicates within the MiSeq amplicon regions.
- `ipd-{locus}-{animal acronym}-{run date}_cleaned.cdna.fasta` is a FASTA-formatted list of cDNA alleles.
- `ipd-{locus}-{animal acronym}-{run date}_cleaned.gdna.fasta` is a FASTA-formatted list of gDNA alleles.
- `ipd-{locus}-{animal acronym}-{run date}_cleaned.immunowes.fasta` is a FASTA-formatted list of alleles from gDNA or MHC exon 2.

## Acknowledgements

This workflow was developed with support from NIH/NIAID contract HHSN272201100026C. For more information, visit our [public MHC Contract portal](https://dholk.primate.wisc.edu/_webdav/dho/grants/mhc_contract/web_portal/@files/prototype/index.html).
