# Modules
library(readr)
library(ggplot2)
library(reshape2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

population <- read_csv('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/vcfs/all.csv',
                   col_types = cols(Generation = col_number(), Minor_allele_frequency = col_number(),
                                    Cage = col_character()))
#View file if required
#View(allele)

cage_plot <- ggplot(population, aes(x = Generation, y = Minor_allele_frequency, color = Cage)) + geom_line() +
  theme_grey() +
  labs(title = 'Alternate allele frequency per generation \nOverdominance (h=20.0, s=0.5)', y = 'Alternate allele frequency', x = 'Generations') +
  scale_x_continuous(breaks = c(0,2,4,8,12,20,28,36,44,56)) +
  scale_color_discrete(breaks = c(1,2,3,4,5,6,7,8,9,10))
  
cage_plot