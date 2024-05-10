#clearing the environment
rm(list = ls())

#BiocManager
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
require('BiocManager')

#different load in technique <3
packages <- c(
  'tidyverse', 'readr', 'ggplot2', 'car',
  'MuMIn', 'vegan', 'effects', 'ggtext',
  'lme4', 'mgcv'
)

install_if_missing <- function(packages) {
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE)
    }
  }
}

#loading in libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(car)
library(MuMIn)
library(vegan)
library(effects)
library(ggtext)
library(lme4)
library(mgcv)
library(knitr)


install_if_missing(packages)

check_installed <- function(packages) {
  missing_packages <- packages[!sapply(packages, requireNamespace, dependencies = TRUE, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    cat("The following packages are not installed:\n")
    cat(paste(missing_packages, collapse = ", "), "\n")
  } else {
    cat("All required packages are installed.\n")
  }
}

check_installed(packages)

setwd('/home/toast/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/biom/')
list.files()

feature_table <- read.csv('core_feature_table.csv', header = TRUE)
colnames(feature_table)[1] <- 'id'

sequence_table <- read_tsv('sequence_table_merged_rmvd.tsv')
sequence_table <- sequence_table[-1,]
sequence_table <- sequence_table[,-3,4]
sequence_table <- separate(sequence_table, Taxon, into = c('Kingdom', 'Phylum', 'Class',
                                                           'Order', 'Family', 'Genus',
                                                           'Species'), sep = '; ')
sequence_table <- sequence_table[,-9]
write.csv(sequence_table, 'family_level_coracle.csv', row.names = FALSE)

joined_dataset <- right_join(sequence_table, feature_table, by = 'id')
joined_dataset <- joined_dataset[,-2]
colnames(joined_dataset)[1] <- 'ASV'
joined_dataset$ASV <- lapply(strsplit(joined_dataset$ASV, '; '), unlist)
print(joined_dataset$ASV)

core_features <-read.csv('core_feature_table.csv', header = TRUE)
core_features <- as.data.frame(core_features)
colnames(core_features) <- core_features[1,]
core_features <- core_features[-1,]
colnames(core_features)[1] <- 'id'
metadata_coracle <- read_tsv('newt_metadata.tsv')

core_features <- core_features %>%
  mutate(id = gsub("\\.", "-", id))
core_features <- core_features[-84,]

intersect_samples <- intersect(core_features$SampleID, metadata_coracle$SampleID)
filtered_dataset <- metadata_coracle %>% filter(SampleID %in% intersect_samples)
dim(filtered_dataset)
dim(core_features)
merged <- bind_cols(filtered_dataset, core_features)
merged_dataframe <- merged[,-3]
write.csv(merged_dataframe, 'x_file.csv', row.names = FALSE)

y <- merged_dataframe[,1:2]
write.csv(y, 'y_file.csv', row.names = FALSE)
