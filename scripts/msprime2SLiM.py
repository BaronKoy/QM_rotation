# Simulation of Michael Jardine cage sequencing experiment. Reuter lab 2021
# Written with the purpose of changing random seed and rerunning the simulation from scratch every run
import msprime, tskit, time, json, sys, random, pyslim
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display

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
mut_model = msprime.SLiMMutationModel(type = 1)
# this model is not required as we will be using the pyslim mutation metadata model
# mut_model = msprime.BinaryMutationModel(state_independent=True)

# Assign random seed value - this will generate a different random seed string value everytime the simulation is run
ancestry_seed = random.randint(1, 10000)
mutation_seed = random.randint(10001, 20000)
rng_seed = random.randint(20001, 30000)

# Ancestry simulation. Recombination rate for 3R retrieved from: https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_DroMel
ts = msprime.sim_ancestry(
    model = anc_model, samples = sample_size,
    population_size = pop_size, sequence_length = seq_length,
    recombination_rate = 1.71642e-08, random_seed = ancestry_seed #random seed point
)

# Annotate metadata with pyslim
ts = pyslim.annotate(ts, model_type = 'WF', tick = 1, stage = 'late')

# Add mutations. Mutation rate for 3R retrieved from: https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html#sec_catalog_DroMel
mut = msprime.RateMap.uniform(sequence_length = seq_length, rate = 5.49e-09)
mts = msprime.sim_mutations(ts, rate = mut,
    model = mut_model, random_seed = mutation_seed # random seed point
)

#Build map of SLiM mutations & assign selection coefficients
# TODO fix issue of assigning selection coefficients
'''
rng = np.random.default_rng(seed=rng_seed)
tables = mts.tables
tables.mutations.clear()
mut = {}
for m in mts.mutations():
    md_list = m.metadata['mutation_list']
    slim_ids = m.derived_state.split(',')
    assert len(slim_ids) == len(md_list)
    for sid, md in zip(slim_ids, md_list):
        md['selection_coeff'] = mut[sid]
    _ = tables.mutations.append(m.replace(metadata = {'mutation_list': md_list}))
assert tables.mutations.num_rows == mts.num_mutations
'''

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
complete_file = open('All_mut_2.csv', 'w') # change this to 'a' to append new outputs if file already exists. 'w' to create new file
allele_file = open('SLiM_AF_count_2.csv', 'w') # change this to 'a' to append new outputs if file already exists. 'w' to create new file
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
for tree in mts.trees():
    kb.append(tree.interval.right/1000)
    mrca = mts.node(tree.roots)
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
f = open('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/Image.svg', 'w')
svg_size = (2000, 1000)
css_string = ('{background-color: white}')
simple_mts = mts.simplify([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
image = simple_mts.draw_svg(y_axis=True,
    x_scale = 'treewise',
    x_lim=(0, 10000), size = svg_size, style = css_string)
print(image, file=f)
print('image.svg successfully created')

# Create tree output file for SLiM
mts.dump('msprime_tree_sequence.trees'); print('.trees file created')

# Create vcf output file for SLiM
vcf_file = open('msprime_out_3.vcf', 'w')
create_vcf = mts.write_vcf(output=vcf_file); print('.vcf created')

# Print end process messages
print(f'Random seeds: Anc: {ancestry_seed}, Mut: {mutation_seed}') # Check random seed strings used
print('Process complete')
print(f'Max root time: {mts.max_root_time}')
print(f'Overall Diversity (pi) = {mts.diversity(mode = "site")}')
afs = mts.allele_frequency_spectrum(sample_sets = None, mode = 'site', polarised = True, span_normalise = False)
print(f'Site frequency spectrum sample: {afs}')

# Figure for SFS
fig, ax = plt.subplots()
ax.bar(np.arange(mts.num_samples + 1), afs, color = 'steelblue') # Plot sfs
# Set labels for title & axis
ax.set_title('Site frequency spectrum (SFS)\nSimulation_4', fontsize = '30',
 fontweight = 'bold', ); ax.set_xlabel('Alt Allele frequency class',
    fontsize = '20'); ax.set_ylabel('Proportion of Alt allele', fontsize = '20')
'''NOTE the lines below limits the x&y axis' to 1000, which in result will remove any bins with frequency
>1000, which is more than likely with this simulation. Only include to produce a cleaner plot
which can be used to demonstrate exponential distribution
which is expected with neutral model simulation'''
ax.set_xlim([0, 1000]) # limit x axis to 1000
ax.set_ylim([0, 1000]) # limit y axis to 1000
plt.savefig('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/SFS_4.png', 
    dpi = 1000, format = 'png',
    bbox_inches = 'tight')
plt.show()