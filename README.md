# phastSim
Fast sequence evolution simulation for SARS-CoV-2 phylogenies and other genomic epidemiological datasets.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[TOC]: #

## Table of contents
- [Installation and dependencies](#installation-and-dependencies)
- [Usage](#usage)
- [Other included scripts](#other-included-scripts)

## Installation and dependencies

### Installation

phastSim is available either through [PyPi](https://pypi.org/) or by cloning this repository directly:

```sh
pip install phastSim
```

or

```sh
git clone https://github.com/NicolaDM/phastSim
```

### Dependencies

phastSim requires the Python packages [numpy](https://numpy.org/), [importlib_resources](https://importlib-resources.readthedocs.io/en/latest/), [ete3](http://etetoolkit.org/), and [biopython](https://biopython.org/). These packages are all available through both `PyPi`:

```sh
pip install numpy importlib-resources numpy ete3 biopython
```

and `conda`:

```sh
conda install -c defaults -c etetoolkit -c conda-forge numpy importlib_resources ete3 biopython
```

Versions used for testing were: python=3.9.1, numpy=1.19.2, importlib_resources=5.1.0, ete3=3.1.2, and biopython=1.78.


## Usage

If installation was performed using `pip`, you can run `phastSim` using:
```sh
mkdir simulation_output
phastSim --path simulation_output/ --seed 7 --createFasta --createInfo \
         --createNewick --createPhylip --treeFile [tree_name.newick] \
         --scale 3.0 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
         --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
         --reference [ref_genome.fasta]
```
where `--treeFile [tree_name.newick]` specifies a tree file in Newick format, and `--reference [ref_genome.fasta]` specifies a genome sequence provided as a single-record FASTA file.

Alternatively, if you have cloned this repository, an example command to run the main `phastSim` script from the terminal is:

```sh
mkdir simulation_output
python bin/phastSimulate --path simulation_output/ --seed 7 --createFasta --createInfo \
                         --createNewick --createPhylip --treeFile phastSim/example/example_sarscov2_tree.newick \
                         --scale 3.0 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
                         --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
                         --reference phastSim/example/MN908947.3.fasta
```
which uses an example SARS-CoV-2 tree included here, specified by `--treeFile phastSim/example/example_sarscov2_tree.newick`, and a reference SARS-CoV-2 genome in FASTA format, provided here and specified by `--reference phastSim/example/MN908947.3.fasta`.


As input, `phastSim` requires a tree in Newick format (specified using `--treeFile`, as shown in the examples above). This script will then simulate sequence evolution along the tree using a Gillespie algorithm, and will output efficiently represented genome sequences.
This approach is more efficient than alternatives when branch lengths are short, such as in SARS-CoV-2 and other genomic epidemiological datasets.
When a large proportion of the genome is expected to mutate on each branch, however, this approach might be slower than traditional methods.
Given a large phylogeny (say, >10,000 tips) and divergence levels typical for genomic epidemiology, this approach will be faster than other methods.
Even when simulating under a codon model, whole genome, and 100,000 sequences, simulations should only take a few seconds.
For small phylogenies (<1,000 tips) or long branches (many mutation events per branch on average), and stationary nucleotide models, the software [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/) by Rambaut and Grassly will typically be more efficient.
This script does not yet allow simulations of indels - for this you will probably need [INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/) by Fletcher and Yang.

Output can be created in FASTA format (`--createFasta`) and/or PHYLIP format (`--createPhylip`).
The standard output is however given as a list of modified bases w.r.t. the reference for each sample.
The history of mutation events can also be created and stored in an annotated Newick format (`--createNewick`).
A .info text file, describing the evolutionoary rates for each site of the genome, can also be created ('--createInfo').

Nucleotide or codon evolution models are allowed.
A nuclotide model is defined by a set of 12 rates, each for a non-diagonal entry of an UNREST substitution matrix https://doi.org/10.1089/1066527041410472.
By default, neutral mutation rates as inferred from SARS-CoV-2 https://doi.org/10.1101/2021.01.14.426705 are assumed.
Any other simpler nucleotide model can be specified by defining the corresponding entries in the UNREST matrix.

We allow users to specify mutation rate variation in multiple ways:
1) With a finite number of discrete mutation rate categories, or a continuous gamma distribution.
2) Jointly with 1., a proportion of non-mutable sites can also be specified.
3) Jointly with both 1. and 2., a new model of rate variation describing hypermutability, as observed for example in SARS-CoV-2, can also be specified. This model is triggered by options `--hyperMutProbs` and `--hyperMutRates`. In this model, small proportions of sites are given a (possibly much higher) mutation rate. For any of these sites, only one specific mutation rate is enhanced. For example, one site might have only the G→T mutation rate increased 100-fold, while all other rates at that site remain the same. This is to model the effects observed in SARS-CoV-2 which are possibly attributable to APOBEC, ROS activity, or other mechanisms which remain unclear (see https://doi.org/10.1101/2021.01.14.426705).

A codon model is built on top of a model of nucleotide mutation and mutation rate variation, as explained in the two paragraphs above.
See also https://doi.org/10.1093/molbev/mss266.
In addition, a codon model has a codon-specific omega parameter that increases or decreases the rate of nonsynonymous changes at the given codon.
Variation in omega can be defined with a finite number of omega classes, or with a continuous gamma distribution.

Thanks to the new algorithm in phastSim, the computational demand is almost not affected by the complexity of chosen evolution model.

Rates in phastSim are normalized so that the expected number of substitutions per nucleotide site per unit branch length at the root is 1.
This is true even when simulating under a codon model. If instead one wants branch lengths to be interpreted as number of expected substitutions per codon per unit time, then one needs to rescale the input tree by a factor of 3 using, for example, option `--scale 0.3333333`. Finally, since the substitution model will usually not  be assumed to be at equilibrium, the total evolution rate might typically decrease a bit moving downstream from the root to the tips of the phylogenetic tree.

## Other included scripts
Other scripts (in the "scripts" subdirectory) which are not required to run `phastSim`, but are useful for support, are also included in this repository. These are:

- **random_tree.py**, which efficiently generates a random birth-death tree with many tips, and

- **compareSimulators.py**, which is used to compare the running time of phastSim with other simulators (Seq-Gen, INDELible and pyvolve).
