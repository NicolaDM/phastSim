# phastSim
Fast sequence evolution simulation for SARS-CoV-2 phylogenies and other genomic epidemiological datasets.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


[TOC]: #

## Table of contents
- [Installation and dependencies](#installation-and-dependencies)
- [Usage](#usage)
- [Other included scripts](#other-included-scripts)
- [Tests](#tests)
- [Tutorial](#tutorial)

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

After cloning the repository, local installation with pip might still be require, for example executing
```sh
pip install -e .
```
within the phastSim directory.


### Dependencies

phastSim requires the Python packages [numpy](https://numpy.org/), [importlib_resources](https://importlib-resources.readthedocs.io/en/latest/), [ete3](http://etetoolkit.org/), [biopython](https://biopython.org/), [protobuf](https://pypi.org/project/protobuf/), and [google](https://pypi.org/project/google/). These packages are all available through both `PyPi`:

```sh
pip install numpy importlib-resources six ete3 biopython protobuf google
```

and `conda`:

```sh
conda install -c defaults -c etetoolkit -c conda-forge numpy importlib_resources six ete3 biopython protobuf google
```

Versions used for testing were: python=3.9.1, numpy=1.19.2, importlib_resources=5.1.0, ete3=3.1.2, biopython=1.78, protobuf=3.15.7 and google=3.0.0.


## Usage

If installation was performed using `pip`, you can run `phastSim` using:
```sh
mkdir simulation_output
phastSim --outpath simulation_output/ --seed 7 --createFasta --createInfo \
         --createNewick --createPhylip --treeFile [tree_name.newick] \
         --scale 0.333333333 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
         --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
         --reference [ref_genome.fasta]
```
where `--treeFile [tree_name.newick]` specifies a tree file in Newick format, and `--reference [ref_genome.fasta]` specifies a genome sequence provided as a single-record FASTA file.

Alternatively, if you have cloned this repository, an example command to run the main `phastSim` script from the terminal is:

```sh
mkdir simulation_output
python bin/phastSim --outpath simulation_output/ --seed 7 --createFasta --createInfo \
                    --createNewick --createPhylip --treeFile phastSim/example/example_sarscov2_tree.newick \
                    --scale 0.333333333 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 \
                    --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon \
                    --reference phastSim/example/MN908947.3.fasta
```
which uses an example SARS-CoV-2 tree included here, specified by `--treeFile phastSim/example/example_sarscov2_tree.newick`, and a reference SARS-CoV-2 genome in FASTA format, provided here and specified by `--reference phastSim/example/MN908947.3.fasta`.


As input, `phastSim` requires a tree in Newick format (specified using `--treeFile`, as shown in the examples above). This script will then simulate sequence evolution along the tree using a Gillespie algorithm, and will output efficiently represented genome sequences.
This approach is more efficient than alternatives when branch lengths are short, such as in SARS-CoV-2 and other genomic epidemiological datasets.
When a large proportion of the genome is expected to mutate on each branch, however, this approach might be slower than traditional methods.
Given a large phylogeny (say, >10,000 tips) and divergence levels typical for genomic epidemiology, this approach will be faster than other methods.
Even when simulating under a codon model, whole genome, and 100,000 sequences, simulations should only take a few seconds.
For small phylogenies (<1,000 tips) or long branches (many mutation events per branch on average), the software [seq-gen](http://tree.bio.ed.ac.uk/software/seqgen/) by Rambaut and Grassly will typically be more efficient.
PhastSim offers a broad choice of sequence evolution models, including indels (similarly to [INDELible](http://abacus.gene.ucl.ac.uk/software/indelible/) by Fletcher and Yang), non-reversible non-stationary substitution models, and codon models.

Output can be created in FASTA format (`--createFasta`) and/or PHYLIP format (`--createPhylip`); in the case of simulations with indels, the FASTA and PHYLIP formats will contain unaligned sequences.
The standard output is given as a list of modified bases w.r.t. the reference for each sample.
The history of mutation events can also be created and stored in an annotated Newick format (`--createNewick`).
A .info text file, describing the evolutionoary rates for each site of the genome, can also be created ('--createInfo').

Nucleotide or codon evolution models are allowed.
A nuclotide model can defined by a set of 12 rates, each for a non-diagonal entry of an UNREST substitution matrix https://doi.org/10.1089/1066527041410472; this allows the inderct specification of any nucleotide substitution model, but the user can alternatively also directly select a GTR, HKY, oor JC model.
By default, neutral mutation rates as inferred from SARS-CoV-2 https://doi.org/10.1101/2021.01.14.426705 are assumed.

We allow users to specify mutation rate variation in multiple ways:
1) With a finite number of discrete mutation rate categories, or a continuous gamma distribution.
2) Jointly with 1., a proportion of non-mutable sites can also be specified.
3) Jointly with both 1. and 2., a new model of rate variation describing hypermutability, as observed for example in SARS-CoV-2, can also be specified. This model is triggered by options `--hyperMutProbs` and `--hyperMutRates`. In this model, small proportions of sites are given a (possibly much higher) mutation rate. For any of these sites, only one specific mutation rate is enhanced. For example, one site might have only the G→T mutation rate increased 100-fold, while all other rates at that site remain the same. This is to model the effects observed in SARS-CoV-2 which are possibly attributable to APOBEC, ROS activity, or other mechanisms which remain unclear (see https://doi.org/10.1101/2021.01.14.426705).

A codon model is built on top of a model of nucleotide mutation and mutation rate variation, as explained in the two paragraphs above.
See also https://doi.org/10.1093/molbev/mss266.
In addition, a codon model has a codon-specific omega parameter that increases or decreases the rate of nonsynonymous changes at the given codon.
Variation in omega can be defined with a finite number of omega classes, or with a continuous gamma distribution.

Thanks to the new algorithm in phastSim, the computational demand is almost unaffected by the complexity of chosen evolution model.

Rates in phastSim are normalized so that the expected number of substitutions per nucleotide site per unit branch length at the root is 1.
This is true even when simulating under a codon model. If instead one wants branch lengths to be interpreted as number of expected substitutions per codon per unit time, then one needs to rescale the input tree by a factor of 3 using, for example, option `--scale 0.3333333`. Finally, since the substitution model will usually not  be assumed to be at equilibrium, the total evolution rate might typically decrease a bit moving downstream from the root to the tips of the phylogenetic tree.

Indels can also be simulated, under a similar model as INDELible, using option `--indels`. Insertion rate and deletion rate can separately bespecified with options `--insertionRate` and `--deletionRate`. Insertion length distribution and deletion length distribution can be specified with options `--insertionLength` and `--deletionLength`, which offer a number of possible types of distributions.


## Other included scripts
Other scripts (in the "scripts" subdirectory) which are not required to run `phastSim`, but are useful for support, are also included in this repository. These are:

- **random_tree.py**, which efficiently generates a random birth-death tree with many tips, and

- **compareSimulators.py**, which is used to compare the running time of phastSim with other simulators (Seq-Gen, INDELible and pyvolve, you will need to install these yourself and provide some command line arguments telling the script where these are installed).

- **some_experiment.sh**, a few bash scripts which automate the experiments used to produce the figures in the paper. 


## Tests
Tests are not required to run `phastSim` but can be helpful for development, and are also included in this repository. Tests make use of the [pytest](https://docs.pytest.org/en/) package. During development, it can be useful to install `phastSim` in editable mode; to do this run

```sh
pip install -e .
```

from the base directory of the project. 

To actually run the tests, do `pytest <name of test>`, and use the `-s` flag if you want the print statements to get printed to stdout. 

## Tutorial

This tutorial covers the most common use cases for phastSim, as well as some edge cases which might otherwise be quite subtle or ambiguous. 

### Basic usage

After installation (perhaps simplest through PyPi: `pip install phastSim`), all of the options available in phastSim are accessed as command-line 
arguments. 
```sh
phastSim --outpath mydirectory/ 
```
Run phastSim using all the default values, and output a succinct txt output file in the output directory. 
The default tree is found in phastSim/example/exampleTree.tree, while the default reference is the SARS-CoV-2 Wuhan reference genome, in
the same folder. The default model is an UNREST nucleotide substitution model, with substitution rates tuned for the SARS-CoV-2 genome. 

```sh
phastSim --outpath mydirectory/ --reference myreference.fasta --mutationRates HKY85 r_transition r_transversion pi_A pi_C pi_G pi_T --seed 10
```
Same as before but using a specified reference, substitution model, and random seed. 

```sh
phastSim --outpath mydirectory/ --reference myreference.fasta --mutationRates HKY85 r_transition r_transversion pi_A pi_C pi_G pi_T \
--seed 10 --createInfo --createFasta --createNewick --createMAT --createPhylip --verbose
```
Same as before but creating more output files:
- FASTA and PHYLIP files
- an info file specifying the parameters at each site in the genome (somewhat redundant in this case as they are all the same)
- a mutation annotated Newick tree
- a MAT protobuf file
 

```sh
phastSim --outpath mydirectory/ --reference myreference.fasta --treeFile mytree.tree --codon --scale 0.3333333333 \
        --createPhylip --createFasta --omegaCategoryProbs 0.1 0.2 0.3 0.4 --omegaCategoryRates 1.0 5.0 2.0 2.3 \
        --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --mutationRate JC69
```

Simulate under a codon model. The software assumes that branch lengths are given in substitutions per nucleotide, so
if they are actually given in substitutions per codon then all the rates should be scaled by a factor of 1/3 to account for this, 
hence the --scale 0.333333333. 
The underlying model for codon mutation is a GY94 model, meaning that each codon has a different value of omega (the ratio of
non-synonymous/synonymous mutations) speficied by the values given in --omegaCategoryRates, which appear with frequencies
given by --omegaCategoryProbs. Furthermore, a total of 5% of nucleotides are 'hypermutable', meaning they have a much inflated
value of mutation from their nucleotide to a randomly selected nucleotide. (In this case, the multipliers for these rates are 100 and 10,
values of hyperMutRates less than 1 are not allowed.)


```sh
phastSim --outpath mydirectory/ --reference myreference.fasta --treeFile mytree.tree --codon --scale 0.3333333333 \
        --createPhylip --createFasta --omegaCategoryProbs 0.1 0.2 0.3 0.4 --omegaCategoryRates 1.0 5.0 2.0 2.3 \
        --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --mutationRate JC69 \
        --indels --insertionRate GAMMA 0.1 1.0 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.9 --deletionLength NEGBINOMIAL 2 0.95 \
        --rootGenomeFrequencies 0
```

The same simulation as before but also including an indel model. The insertion rate for each codon is drawn from a gamma distribution with
parameters alpha=0.1 and beta=1.0 (hence a mean of 0.1), whereas the deletion rate for each codon is 0.1 uniformly. The insertion length
is drawn from a geometric distribution with p=0.9, and hence the mean insertion length is 10/9=1.11111 codons (since this is a codon simulation). 
Because --rootGenomeFrequencies has been set to 0 (rather than left at the default value which is tuned for SARS-CoV-2), phastSim will
count the frequencies of each codon during its initialisation, and draw from these frequencies when generating indels. 
This count assumes the whole genome is in frame 0, which may not be accurate. To override this, it would instead be possible to supply
frequencies for each of the 64 codons (in alphabetical order) to --rootGenomeFrequencies. 

```sh
phastSim --outpath mydirectory/ --reference myreference.fasta --treeFile mytree.tree \
        --categoryRates 1.0 1.1 1.5 2.0 --categoryProbs 0.9 0.01 0.04 0.05 --noHierarchy --createInfo
```

Use a nucleotide model with the default mutation rates coming from SARS-CoV-2 (since --mutationRates has been left unspecified). 
Each site in the genome has a rate coming from one of the 4 categories. The --noHierarchy option is compatible with these settings and should
improve performance. The only outputs are a concise text file containing a list of mutations at each tip of the tree, and an info file
specifying the category at each nucleotide. 

### How are genome positions indexed?

- In all output files, genomes are 1 indexed (the first nucleotide of the genome is at position 1). 

### Can I run different models for different parts of the genome?

- No, currently phastSim does not a configuration parameter analogous to the partitions in INDELible and Seq-Gen, though if 
we are sent multiple requests for this feature then we will add it. 

### How can I run phastSim directly in a python terminal or jupyter notebook?

- This isn't currently recommended, though it is possible, essentially by copying the code in the phastSim/bin/phastSim directory. 

### Should I change the mutation rates if I am using a codon model?

- phastSim assumes that branch lengths are given as substitutions per nucleotide, even if using a codon model. If
they are given as substitutions per codon, you should set --scale 0.33333333 to correctly adjust for this. 
The effect of --scale is to multiply the rate of each site (nucleotide or codon, depending on the model chosen) by --scale. 

### When can I use the --noHierarchy option?

- the --noHierarchy option is intended for simpler models (where it may result in a substantial increase in performance). It cannot be used with a codon model, with indels, or with non-homogeneous substitution rates (--alpha), though it can be used with --categoryRates. 

### What is the alternative output format?

- phastSim produces a text output file which is far more succinct than full FASTA files. This text file contains a list of mutations for each leaf of the tree, one on each time. These mutations can be written in one of two ways, either tab delimited as nucleotides_from \t position \t nucleotides_to, or using control characters i for insertions and d for deletions. 

### How do indels work in the codon model?

- Insertions and deletions operate on whole codons, and insertions are generated by randomly choosing codons according to the parameter 
--rootGenomeFrequencies. If this parameter is not specified, but a reference genome is supplied, then the root genome frequencies are calculated
by counting codons in the reference (using frame 0) - however to turn this behaviour on you need to specifically set --rootGenomeFrequencies to 0, as otherwise it will use a default value, which has been tuned for SARS-CoV-2. 

### I want to create a large number of small trees, should I use phastSim?

- phastSim has a relatively high startup time cost, due to creating a genome-tree data structure in memory, so it is likely to be slower than other simulators when creating a large number of small trees (e.g. 10000 trees each with 1000 tips). Because this startup cost only happens once, modifying the source code in the bin directory to repeatedly run several simulations using different random seeds should solve this issue, and if we are sent several requests for this feature then we will implement it. 

### What is insertionPos in the info file?

- If you use the --indels option then every nucleotide has 2 indices determining its position in the genome. The first of these is the insertionPos, which is an integer, where 0 refers to the root genome and n > 0 refers to the nth insertion created during the simulation. The 2nd index is simply the position of the nucleotide within its insertion, starting at 1. For example, the indices 0, 1049 would refer to the 1049th nucleotide on the reference genome. The indices 4, 3 would refer to the 3rd nucleotide on insertion number 4. 
