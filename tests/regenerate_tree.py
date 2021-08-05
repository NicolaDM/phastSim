import os
from ete3 import Tree

# generate a balanced tree with 2**10 leaves using ete3
def grow_tree(my_tree):
    for c in my_tree.get_leaves():
        c.add_child(name=f"{c.name + 'a'}", dist=12.0/30000.0)
        c.add_child(name=f"{c.name + 'b'}", dist=12.0/30000.0)

t = Tree()
for i in range(10):
    grow_tree(t)
    
t.write(format=5, outfile="../phastSim/simulation_output_3/sars-cov-2_test_tree.tree")

#t = Tree()
#t.populate(2**10, random_branches=True, branch_range=(5,10))
#t.write(format=1, outfile="test_tree.tree")


# run phastSim using the example tree
command_string = """phastSim --outpath simulation_output_3/ 
--seed 7 --createInfo --createNewick --createPhylip --createMAT 
--treeFile ./simulation_output_3/sars-cov-2_test_tree.tree 
--scale 0.333333333 --invariable 0.1 --alpha 1.0 --omegaAlpha 1.0 
--hyperMutProbs 0.01 0.01 --hyperMutRates 20.0 200.0 --codon 
--reference ./phastSim/example/MN908947.3.fasta --indels 
--insertionRate CONSTANT 0.01 --deletionRate CONSTANT 0.01 
--insertionLength GEOMETRIC 0.1 --deletionLength GEOMETRIC 0.1 --createFasta"""

os.system(command_string)

# run mafft
fasta_output_file = ""
alignment_file = ""
os.system(f"mafft --auto {fasta_output_file} > {alignment_file}")

# run iqtree
os.system(f"iqtree -s {alignment_file} -m GTR+F+R4")

# compare the original tree with the result from iqtree
os.system(f"ete3 compare -t sars-cov-2_test_tree.tree -r sars-cov-2_simulation_alignment.fasta.treefile --unrooted")