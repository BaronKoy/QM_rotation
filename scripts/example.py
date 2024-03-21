# Running the example tutorial found here: https://tskit.dev/tutorials/terminology_and_concepts.html#sec-terminology-nodes
import msprime, tskit
from IPython.display import SVG, display
ts = tskit.load('/home/baron/Documents/rotation_2/tskit_tutorial/data/basics.trees')

# Length of spanning genome
print('Genome length:', int(ts.sequence_length))

# Draw tree sequence from loaded data
t_image = ts.draw_svg(
    y_axis=True, y_gridlines=True,
    time_scale='log_time', y_ticks=[0, 3, 10, 30, 100, 300, 1000]
)

# Output tree image to file 'image.svg'
f = open('image.svg', 'w')
print(t_image, file=f)