import msprime, tskit
from IPython.display import SVG, display

# Simulation for cage experiments from Michael Jardine, Reuter lab 2021
# Neutral simulation
# Parameters are low at the moment to avoid computational heavy process on local machine
ts = msprime.sim_ancestry(samples=45, sequence_length=13500,
    random_seed=12345, recombination_rate=1.71642e-08
)
mts = msprime.sim_mutations(ts, rate=5.49e-09, random_seed=5678)

# Displays the chromosomal position, ALT & REF, and genotypes
for var in mts.variants():
    print(var.site.position, var.alleles, var.genotypes, sep='\t')

# Display the number of generated mutations
print('Number of mutations:', mts.num_mutations)

# Visualisation for the simulation
t_image = mts.draw_svg(y_axis=True, y_gridlines=True, tree_height_scale=10,
    time_scale='log_time', y_ticks=[2, 4, 8, 12, 20, 28, 36, 44, 56],
    x_label='Genome sequence'
)

# Output to svg image file
f = open('output.svg', 'w')
print(t_image, file=f)