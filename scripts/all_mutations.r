# import modules
library(readr)
library(ggplot2)
options(scipen = 999) # prevent scientific notation

# Read in output from msprime simulation
all_mut <- read_csv("Documents/rotation_2/QM_rotation/scripts/All_mut.csv", 
                     col_types = cols(Position = col_number(), 
                                      Ref_allele = col_number()))
# Include to view mutations table
# View(all_mut)

# ggplot genomic position vs reference allele count
initial_plot <- ggplot(all_mut, aes(x = Position, y = Ref_allele)) + 
  geom_point(shape = 19, size = 1, color = 'blue') +
  labs(title = '% of Ancestral alleles across chromosome arm 3R', y = 'Ancestral allele(0) %',
       x = 'Genomic position (Chromosome 3R)') +
  theme_grey() + geom_hline(yintercept = 5, color = 'red') +
  geom_hline(yintercept = 15, color = 'red') +
  geom_hline(yintercept = 85, color = 'red') +
  geom_hline(yintercept = 95, color = 'red') +
  geom_hline(yintercept = 50, color = 'black')
final_plot <- initial_plot + theme(plot.title = element_text(size = 50, face = 'bold'),
                     axis.text.x = element_text(size = 20),
                     axis.title.x = element_text(size = 40),
                     axis.text.y = element_text(size = 20),
                     axis.title.y = element_text(size = 40))
final_plot