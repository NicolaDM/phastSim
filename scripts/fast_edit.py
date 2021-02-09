import subprocess
import time
import mmap
import numpy as np
import random


def execute(command):
    # create the unix process
    running = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               encoding='utf-8', shell=True)
    # run on a shell and wait until it finishes
    stdout, stderr = running.communicate()


def edit_sites(filename, sites, new_characters):
    # open the file
    with open(filename, mode="r+") as f:
        # instantiate a memory map with writing acess
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_WRITE)
        # loop through the sites we want to change and assign new characters
        for s in range(sites.shape[0]):
            mm[sites[s]: sites[s] + 1] = new_characters[s].encode('ascii')
        # flush the changes from memory to disk
        mm.flush()


if __name__ == "__main__":
    # number of tips in the tree
    n_copies = 10000
    reference_file = "phastSim/example/MN908947.3.fasta"
    # unix command that appends a file to itself "n_copies" times
    command= f"yes {reference_file} | head -n {n_copies} | xargs cat >ref.fa"

    # generate random sites and random characters
    n_sites = 10000
    sites = np.random.choice(int(1e5), size=n_sites, replace=False)
    new_characters = random.choices(["W", "X", "Y", "Z"], k=n_sites)
    new_characters = np.array(new_characters)

    tic = time.time()
    execute(command)
    edit_sites("ref.fa", sites, new_characters)
    toc = time.time()
    print(toc - tic)




