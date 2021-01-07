# phastSim
Fast sequence evolution simulation for SARS-CoV-2 and similar.

Input a reference genome to be used as root (sars-cov-2 example provided), and a newick tree, and this script will simulate sequence evolution using a Gillespie algorithm and efficiently storing ancestral genomes.
This approach is more efficient than alternative ones when branch lengths are short, such as in sars-cov-2 and other genomic epidelimiological datasets.
When a large proportion of the genome is expected to mutate on each branch, however, this approach might be slower than traditional methods.

Output can be created in fasta and/or phylip format.
The standard output is however given as a list of modified bases wrt the reference for each sample.
The history of mutation events can also be created and stored in an annotated newick format.
For the case in which one wants to simulate rate variation along the genome, a .info text file is also created describing the catogries of each site.

For now, only nucleotide mutation models are allowed.
Also, only two types of rate variation are allowed:
1) Standard rate variation, specified by options --categoryProbs and --categoryRates , specifying respectively the frequencies of different categories and their rates; so far only a discrete number of categories is allowed with pre-specified rates.
2) A new model of rate variation describing the hypermutability observed in SARS-CoV-2. This model is triggered by options --hyperMutProbs and --hyperMutRates . In this model, small proportions of sites are given a (possibly much higher) mutation rate. However,for any of these sites, only one specific mutation rate is enhanced. For example, one site might have only the G->T mutation rate increased 100-fold, while all other rates at that site remain the same. This is to model the effects observed in SARS-CoV-2 and possibly attributable to APOBEC and ROS activity (or other still unclear mechanisms).

No matter what model of rate variation is used, the rates are normalized so that the expected number of substitutions per site per unit branch length at the root is 1. Since the substitution model might be not at equilibrium, the total rate might typically slightly decrease downstream of the root.

An example command line is:

python3 efficientSimuSARS2.py --path /simulations_folder/ --seed 7 --createFasta --createNewick --createPhylip --treeFile rob-11-7-20_new.newick --scale 3.0 --scale --categoryProbs 0.2 0.3 0.5 --categoryRates 0.01 0.2 1.0 --hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0

In order to work it will need installation of several python modules including ete3 and numpy.
