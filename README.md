# phastSim
Fast sequence evolution simulation for SARS-CoV-2 phylogenies and other genomic epidemiological datasets.

The required input is a newick tree. This script will then simulate sequence evolution along the tree using a Gillespie algorithm and efficiently representing genome sequences.
This approach is more efficient than alternative ones when branch lengths are short, such as in sars-cov-2 and other genomic epidelimiological datasets.
When a large proportion of the genome is expected to mutate on each branch, however, this approach might be slower than traditional methods.
Given a large phylogeny (say, >10,000 tips) and divergence levels tyical for genomic epidemiology, this approach will be fater than other methods.
Even when simulating under a codon model, whole genome, and 100,000 sequences, simulation should take only a few seconds.
For small phylogenies (<1,000 tips) and nucleotide models, the software seq-gen by Rambaut and Grassly is instead be tipically more efficient (due to fixed costs, for example associated with loading Python packages).
This script does not allow yet simulations of indels - for this you will probably need INDELible by Fletcher and Yang.

Output can be created in fasta and/or phylip format.
The standard output is however given as a list of modified bases wrt the reference for each sample.
The history of mutation events can also be created and stored in an annotated newick format.
An .info text file is also created describing the rate catogries of each site of the genome.

Nucleotide or codon evolution are allowed.
A nuclotide model is defined by a set of 12 rates, each for a non-diagonal entry of a UNREST substitution matrix https://doi.org/10.1089/1066527041410472 .
By default, neutral mutation rates as inferred from SARS-CoV-2 https://doi.org/10.1101/2021.01.14.426705 are assumed.
Any other simpler nucleotide model can be specified by defining the corresponding entries in the UNREST matrix.

We allow to specify mutation rate variation in multiple ways:
1) With a finite number of discrete mutation rate categories, or a continuous gamma distribution.
2) possibly jointly, with a proportion of un mutable sites.
3) Possibly jointly with the 2 above, a new model of rate variation describing the hypermutability observed in SARS-CoV-2. This model is triggered by options --hyperMutProbs and --hyperMutRates . In this model, small proportions of sites are given a (possibly much higher) mutation rate. For any of these sites, only one specific mutation rate is enhanced. For example, one site might have only the G->T mutation rate increased 100-fold, while all other rates at that site remain the same. This is to model the effects observed in SARS-CoV-2 and possibly attributable to APOBEC and ROS activity (or other still unclear mechanisms) https://doi.org/10.1089/1066527041410472 .

A codon model is built on top of a model of nucleotide mutation and mutation rate variation, as explained in the two paragraphs above.
See also https://doi.org/10.1093/molbev/mss266 .
In addition, a codon model has a codon-specific omega parameter that increases or decreases the rate of nonsynonymous changes at the given codon.
Variation in omega can be define with a finit number of omega classes, or with a continuous gamma distribution.

Thanks to algorithm in phastSim, the computational demand is not considerably affected by the chosen evolution model.

Rates in phastSim are normalized so that the expected number of substitutions per nucleotide site per unit branch length at the root is 1.
This is true even when simulating under a codon model; in instead one wants branch lengths to be interpreted as number of expected substitutions per codon er unit time, then one needs to rescale the input tree by a factor of 3 using for example option --scale 0.3333333 . Finally, since the substitution model will usually not  be assumed to be at equilibrium, the total evolution rate might typically decrease a bit moving downstream from the root to the tips of the phylogenetic tree.

An example command to run the script from terminal is:

python3 efficientSimuSARS2.py --path /simulations_folder/ --seed 7 --createFasta --createNewick --createPhylip --treeFile rob-11-7-20_new.newick --scale 3.0 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon --reference MN908947.3.fasta

In order to work it will need installation of several python modules including ete3 and numpy.

Other scripts, not required to run phastSim, but useful for support, are also included in this distribution.

random_tree.py is a script for efficiently generating a random birth-death tree with many tips.

compareSimulators.py is a custom script used to compare the running time of phastSim with other simulators (Seq-Gen, INDELible and pyvolve).
