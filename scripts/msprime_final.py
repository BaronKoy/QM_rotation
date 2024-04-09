# Simulation of Michael Jardine cage sequencing experiment. Reuter lab 2021
# Written with the purpose of changing random seed and rerunning the simulation from scratch every run
import msprime, tskit, time, json, sys, random
import matplotlib.pyplot as plt
import numpy as np

# Population, sampling, & sequence parameters
pop_size = 5000; seq_length = 32079331; sample_size = 5000 # Based on Michael Jardine cage experiments & Drosophila chromosome arm 3R

# Use block code to simulate cage populations
# NOTE current simulations runs on a single population of Ne = 5000
'''demography = msprime.Demography()
demography.add_population(name='Snine_Lone_one', initial_size=500)
demography.add_population(name='Snine_Lone_two', initial_size=500)
demography.add_population(name='Snine_Lone_three', initial_size=500)
demography.add_population(name='Snine_Lone_four', initial_size=500)
demography.add_population(name='Snine_Lone_five', initial_size=500)
demography.add_population(name='Sone_Lnine_one', initial_size=500)
demography.add_population(name='Sone_Lnine_two', initial_size=500)
demography.add_population(name='Sone_Lnine_three', initial_size=500)
demography.add_population(name='Sone_Lnine_four', initial_size=500)
demography.add_population(name='Sone_Lnine_five', initial_size=500)'''

# Assign ancestral & mutation models
# "duration" can be specified to msprime.StandardCoalescent object if generation time should be included in simulation
# Current model simulates until all trees have coalesced
anc_model = msprime.StandardCoalescent()
mut_model = msprime.BinaryMutationModel(state_independent=True)

# Assign random seed value - this will generate a different random seed string value everytime the simulation is run
ancestry_seed = random.randint(1, 10000)
mutation_seed = random.randint(10001, 20000)

# Ancestry simulation. Recombination rate for 3R retrieved from: https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_DroMel
ts = msprime.sim_ancestry(
    model=anc_model, samples=sample_size,
    population_size=pop_size, sequence_length=seq_length,
    recombination_rate=1.71642e-08, random_seed=ancestry_seed #random seed point
)

# Add mutations. Mutation rate for 3R retrieved from: https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_DroMel
mts = msprime.sim_mutations(ts, rate=5.49e-09,
    model=mut_model, random_seed=mutation_seed # random seed point
)

# Tree & external nodes count check 
print('Number of Trees:', mts.num_trees); print('Number of mutations:', mts.num_mutations)
print('Number of external nodes:', mts.num_samples); print('Number of samples:', mts.num_individuals)

# Tree count
for tree in mts.trees():
    print(f'Tree {tree.index} covers {tree.interval}')
    if tree.index >= 5:
        print('...')
        break # Output first 5 trees to save time
print(f"Tree {mts.last().index} covers {mts.last().interval}") # Jump to last tree

# Processing time calculation
# Check to determine if all trees have coalesced
elapsed = time.time()
for tree in mts.trees():
    if tree.has_multiple_roots:
        print("Tree {tree.index} has not coalesced")
        break
else:
    elapsed = time.time() - elapsed
    print(f"All {mts.num_trees} trees coalesced")
    print(f"Completed in {elapsed:.6g} secs")

# Check & display genotypes for the first 5 external nodes
print('Generating genotypes for the first 5 external nodes...')
print('0-5 site genotypes:')
np.set_printoptions(linewidth=200)
for v in mts.variants():
    alleles = np.array(v.position)
    print(f'Site {v.site.id}: {v.genotypes}')
    if v.site.id >= 5:
        print('...')
        break

# AF calculation per mutation site
# Ensure this step is done to output AF values for SLiM input
complete_file = open('All_mut.csv', 'a') # change this to 'a' to append new outputs if file already exists. 'w' to create new file
allele_file = open('SLiM_AF_count.csv', 'a') # change this to 'a' to append new outputs if file already exists. 'w' to create new file
print('Position,Ref_allele', file=complete_file)
print(file=allele_file)
for var in mts.variants():
    locat = int(var.position)
    af = np.array(var.genotypes)
    af_count = np.count_nonzero(af == 1) # Count alternate allele
    af_per = (af_count/100)
    if 5 <= af_per <= 15 or 85 <= af_per <= 95:
        print(f'{locat},{af_per}', file=allele_file)
    print(f'{locat},{af_per}', file=complete_file)
    '''if var.index >= 5:
        break''' # Include to save script having to go through all mutations!

# The next block of code creates a MRCA graph image
# TODO - Current error in which graph can not be generated
'''
print('Creating MRCA image...')
kb = [0]  # Starting genomic position
mrca_time = []
for tree in ts.trees():
    kb.append(tree.interval.right/1000)
    mrca = ts.node(tree.roots)
    mrca_time.append(mrca.time)
plt.stairs(mrca_time, kb, baseline=None)
plt.xlabel("Genome position (kb)")
plt.ylabel("Time of root (or MRCA) in generations")
plt.yscale("log")
plt.show()
print('Creating image.svg file')
'''

print('Creating simple svg image...')
# Output 0-10000 genome region to a file (image.svg) to check trees, limit external nodes to 10 for computational load
f = open('image.svg', 'w')
simple_mts = mts.simplify([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
image = simple_mts.draw_svg(y_axis=True, x_lim=(0, 10000))
print(image, file=f)
print('image.svg successfully created')

# Create tree output file for SLiM
mts.dump('msprime_tree_sequence.trees'); print('.trees file created')

# Create vcf output file for SLiM
vcf_file = open('msprime_out.vcf', 'w')
create_vcf = mts.write_vcf(output=vcf_file); print('.vcf created')
print(f'Random seeds: Anc: {ancestry_seed}, Mut: {mutation_seed}') # Check random seed strings used
print('Process complete')