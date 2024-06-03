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

cage_plot <- ggplot(population, aes(x = Generation, y = Minor_allele_frequency, color = Cage)) + geom_line()
cage_plot