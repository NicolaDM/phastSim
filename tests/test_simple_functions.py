import pytest
from types import SimpleNamespace
import numpy as np
from ete3 import Tree
from phastSim import phastSim




def test_phastSim_imported_correctly():
    
    # if this test failed something has gone very wrong and most likely it is due to an import failing. 
    c = phastSim.Constants()
    print(c.alleles)
    assert c.alleles == {'A':0, 'C':1, 'G':2, 'T':3}

def create_standard_mock_arguments():
    return SimpleNamespace(
        outpath="_",
        reference="_",
        rootGenomeLength=10,
        rootGenomeFrequencies=[],
        treeFile="_",
        scale=1,
        seed=1,
        alpha=0,
        invariable=0.0,
        mutationRates=[0.039, 0.310, 0.123, 0.140, 0.022, 3.028, 0.747, 0.113, 2.953, 0.056, 0.261, 0.036],
        codon=False,
        omegaAlpha=0.0,
        omegaCategoryProbs=[1.0],
        omegaCategoryRates=[1.0],
        categoryProbs=[1.0],
        categoryRates=[1.0],
        hyperMutProbs=[],
        hyperMutRates=[],
        noHierarchy=False,
        verbose=True,
        outputFile="_",
        createNewick="_",
        createFasta="_",
        createPhylip="_",
        createInfo="_"
    )

def test_phastSimRun_initialises():
    
    mock_args = create_standard_mock_arguments()

    sim_run = phastSim.phastSimRun(mock_args)

    assert len(sim_run.create_rootGenome_nuc()[1]) == 10
    assert len(sim_run.create_rootGenome_codon()[1]) == 9

    assert sim_run.init_substitution_rates()[0][0] == -0.472
    


def test_GenomeTree_hierarchical_initialises():

    mock_args = create_standard_mock_arguments()
    sim_run = phastSim.phastSimRun(mock_args)

    t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);")
    print(t)

    #genome_tree = phastSim.GenomeTree_hierarchical(
    #    nCodons=
    #)

    # TODO - the rest of this test. 