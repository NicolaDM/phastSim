# phastSim
Fast sequence evolution simulation for SARS-CoV-2 and similar.

Input a reference genome to be used as root (sars-cov-2 example provided), and a newick tree, and this script will simulate sequence evolution using a Gillespie algorithm and efficiently storing ancestral genomes.
This approach is more efficient than alternative ones when branch lengths are short, such as in sars-cov-2 and other genomic epidelimiological datasets.
When a large proportion of the genome is expected to mutate on each branch, however, this approach might be slower than traditional methods.
Output can be created in fasta and/or phylip format. 
The standard output is however given as a list of modified bases wrt the reference for each sample.
The history of mutation events can also be created and stored in an annotated newick format.

An example command line is:

python3 efficientSimuSARS2.py --path /simulations_folder/ --seed 7 --createFasta --createNewick --createPhylip --treeFile rob-11-7-20_new.newick --scale 3.0

In order to work it will need installation of several python modules including ete3 and numpy.
