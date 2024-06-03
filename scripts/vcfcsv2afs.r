# Modules
library(readr)
library(ggplot2)
options(scipen = 999) # prevent scientific notation for genomic position on x axis

alleles <- read_csv('/home/baron/Documents/rotation_2/QM_rotation/scripts/outputs/vcfs/cage_1/afs/gen2.csv', # load input csv
                col_types = cols(Location = col_number(),
                                 Minor_allele_frequency = col_number()))
View(alleles)

plot <- ggplot(alleles, aes(x = Location, y = Minor_allele_frequency)) + geom_line()
plot