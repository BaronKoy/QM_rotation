# import modules
library(readr)
library(ggplot2)
options(scipen = 999) # prevent scientific notation

# Read in output from msprime simulation
all_mut <- read_csv("Documents/rotation_2/QM_rotation/outputs/All_mut.csv", 
                     col_types = cols(Position = col_number(), 
                                      Ref_allele = col_number()))
# Include to view mutations table
# View(all_mut)

# ggplot genomic position vs reference allele count
initial_plot <- ggplot(all_mut, aes(x = Position, y = Ref_allele)) + 
  geom_point(shape = 19, size = 2, color = 'deepskyblue3') +
  labs(title = '% of Ancestral alleles across chromosome arm 3R', y = 'Ancestral allele(0) %',
       x = 'Genomic position (Chromosome 3R)') +
  # Intercept lines for 5, 15, 85, & 95%
  theme_grey() + geom_hline(yintercept = 5, color = 'red') +
  geom_hline(yintercept = 15, color = 'red') +
  geom_hline(yintercept = 85, color = 'red') +
  geom_hline(yintercept = 95, color = 'red') +
  geom_hline(yintercept = 50, color = 'black') # Black line added for 50%
# Title & axis text size adjustments
final_plot <- initial_plot + theme(plot.title = element_text(size = 50, face = 'bold'),
                     axis.text.x = element_text(size = 20),
                     axis.title.x = element_text(size = 40),
                     axis.text.y = element_text(size = 20),
                     axis.title.y = element_text(size = 40)) +
  coord_cartesian(xlim = c(0, 32079331), ylim = c(0, 100), expand = FALSE)
final_plot