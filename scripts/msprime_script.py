# Simulation of Michael Jardine cage sequencing experiment. Reuter lab 2021
import msprime, tskit, time
import matplotlib.pyplot as plt
import numpy as np
from random import randint

# Population parameters (note these might have to be lower if running from local machine)
pop_size = 5000; seq_length = 32079331; sample_size = 90
demography = msprime.Demography()

# Set parameters from cage populations
demography.add_population(name='Snine_Lone_one', initial_size=500)
demography.add_population(name='Snine_Lone_two', initial_size=500)
demography.add_population(name='Snine_Lone_three', initial_size=500)
demography.add_population(name='Snine_Lone_four', initial_size=500)
demography.add_population(name='Snine_Lone_five', initial_size=500)
demography.add_population(name='Sone_Lnine_one', initial_size=500)
demography.add_population(name='Sone_Lnine_two', initial_size=500)
demography.add_population(name='Sone_Lnine_three', initial_size=500)
demography.add_population(name='Sone_Lnine_four', initial_size=500)
demography.add_population(name='Sone_Lnine_five', initial_size=500)

# Specify parameters for ancestral model
mut_model = msprime.JC69()

# Begin with the simulation of the ancestry for tree sequence
ts = msprime.sim_ancestry(
    samples=sample_size,
    population_size=pop_size, sequence_length=seq_length,
    recombination_rate=1.71642e-08, random_seed=12345
)

# Add mutations
mts = msprime.sim_mutations(ts, rate=5.49e-09,
    random_seed=54321, model=mut_model
)

# Parameters count check
print('Number of Trees:', mts.num_trees); print('Number of mutations:', mts.num_mutations)
print('Number of external nodes:', mts.num_samples); print('Number of samples:', mts.num_individuals)

# Tree count
for tree in ts.trees():
    print(f'Tree {tree.index} covers {tree.interval}')
    if tree.index >= 5:
        print('...')
        break
print(f"Tree {ts.last().index} covers {ts.last().interval}") # Jump to last tree

# Processing time calculation
elapsed = time.time()
for tree in ts.trees():
    if tree.has_multiple_roots:
        print("Tree {tree.index} has not coalesced")
        break
else:
    elapsed = time.time() - elapsed
    print(f"All {ts.num_trees} trees coalesced")
    print(f"Completed in {elapsed:.6g} secs")

# Display genotypes for external nodes
np.set_printoptions(linewidth=200)
for v in mts.variants():
    print(f'Site {v.site.id}: {v.genotypes}')
    if v.site.id >= 5:
        print('...')
        break
print('First 5 tree genotypes available')
print('Creating MRCA image...')

#  Display time to MRCA as graph image
kb = [0]  # Starting genomic position
mrca_time = []
for tree in ts.trees():
    kb.append(tree.interval.right/1000)
    mrca = ts.node(tree.root)
    mrca_time.append(mrca.time)
plt.stairs(mrca_time, kb, baseline=None)
plt.xlabel("Genome position (kb)")
plt.ylabel("Time of root (or MRCA) in generations")
plt.yscale("log")
plt.show()
print('Creating image.svg file')

# Output 0-5208 genome region to a file (image.svg) to check trees, limit external nodes to 10 for computational time
f = open('image.svg', 'w')
simple_ts = ts.simplify([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
image = simple_ts.draw_svg(x_lim=(0, 5208))
print(image, file=f)
print('Process complete')