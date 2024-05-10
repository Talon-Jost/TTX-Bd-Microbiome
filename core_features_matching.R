rm(list=ls())


#BiocManager
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
require('BiocManager')

#different load in technique <3
packages <- c('qiime2R', 'magrittr', 'purr', 'stringr', 'devtools')


install_if_missing <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE)
    }
  }
}

#loading in libraries
library(qiime2R) #(jbisanz/qiime2R)
library(magrittr)
library(purrr)
library(stringr)
library(devtools)
library(phyloseq) #(joey711/phyloseq)
library(dplyr)
library(readr)


install_if_missing(packages)

setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/core_features/')

cf_table_group = qiime2R::read_qza("core_features_matches/clustered_table.qza")
cf_table_group = as.data.frame(cf_table_group$data)

sum_sequences = cf_table_group %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

sum_sequences = as.data.frame(sum_sequences)

sum_sequences = cbind(Indiv_ID = rownames(sum_sequences), sum_sequences)

rownames(sum_sequences) <- NULL

sum_sequences = sum_sequences[order(sum_sequences$Indiv_ID),]

sum_sequences = as.data.frame(sum_sequences)

sum_data = sum_sequences %>%
  mutate(relative_abundance = (sum_sequences/2200)*100)


pres_ab = qiime2R::read_qza("core_features_matches/core_features_PresAbs_clustered_table.qza")

pres_ab = as.data.frame(pres_ab$data)

rep_samp_richness = pres_ab %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

rep_samp_richness = as.data.frame(rep_samp_richness)

sum_data_2 = cbind(sum_data, rep_samp_richness)
colnames(sum_data_2)[1] = 'SampleID'

write.csv(sum_data_2, 
          'core_microbiome_results.csv', 
          row.names = FALSE)


sum(sum_data_2$sum_sequences)

#        hash                       freq    #samples
# 4cf4fd02fdd37f4390cfc5c57a91305c 	24,831 	59
# dd27df27e55b4607445f2acf1d26b23e 	10,811 	60
# d0099ee91dbd9830ee9ca81d0ee85747 	7,542 	71
# 7c680f7e9da7731a12a0751a0221fc67 	2,111 	59
# 34aa31861212a9f72ffd8bf087ca47d2 	1,664 	64



metadata_2 <- read.csv('2200_full_metadata.csv', header = T)

new_dat <- sum_data_2 %>%
  left_join(metadata_2, by = 'SampleID')
percentage <- new_dat %>%
  mutate((sum_sequences/2200)*100)
