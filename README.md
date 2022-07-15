# Immuno Polymorphism Data Base - Reference Sequence Generator

### Overview

This workflow pulls the latest MHC allele sequences from [IPD](https://www.ebi.ac.uk/ipd/) and prepares a series of species-specific reference databases. These databases will contain MHC alleles that are more than 100 base pairs long for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus macaque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nemestrina, a.k.a. Mame), as well as for all non-human primates included in IPD. At minimum, this workflow should be run after each IPD release of new MHC reference alleles.

The workflow also prepares context-specific reference databases for use with data generated at [AVRL](https://dholk.primate.wisc.edu/project/home/begin.view?):

1. Reference allele sequences for use with high depth-of-coverage whole exome sequence data enriched for immuno-genes (iWES)
2. Reference allele sequences for use with Illumina MiSeq amplicon libraries

### Getting Started

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
docker build --tag ipd-ref-generator:v1_0_5 config/
```

Note that to build the above docker container, you may need to increase the amount of memory alloted to Docker in the Docker Engine preferences.

#### Nextflow Installation

This workflow uses the [NextFlow](https://www.nextflow.io/) workflow manager. We recommend you install NextFlow to your system in one of the two following ways:

##### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

##### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are set and ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in the BASH terminal:

```
nextflow run IPD-ref-generator.nf
```

If the workflow runs partway, but a computer outage or other issue interrupts its progress, no need to start over! Instead, run:

```
nextflow run IPD-ref-generator.nf -resume
```

The workflow's configurations (see below) tell NextFlow to plot the workflow and record run statistics. However, the plot the workflow, note that NextFlow requires the package GraphViz, which is easiest to install via the intructions on [GraphViz's website](https://graphviz.org/download/).

### Acknowledgements
