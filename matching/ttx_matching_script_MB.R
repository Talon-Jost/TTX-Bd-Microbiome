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

setwd('/home/toast/Documents/Thesis_Part_2/Newt_TTX_Project/fasta_files/vsearch/')

TTX_table = qiime2R::read_qza("TTX_Matches_99/clustered_table.qza")
TTX_table = as.data.frame(TTX_table$data)

group_table = qiime2R::read_qza('group_TTX_matches/clustered_table.qza')
group_table = as.data.frame(group_table$data)

tt_group_table = qiime2R::read_qza('taricha_matches/clustered_table.qza')
tt_group_table = as.data.frame(tt_group_table$data)


TotalTTX = TTX_table %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

GroupTotalTTX = group_table %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

tt_group_ttx = tt_group_table %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

TotalTTX = as.data.frame(TotalTTX)
GroupTotalTTX = as.data.frame(GroupTotalTTX)
group_Totalttx = as.data.frame(tt_group_ttx)

TotalTTX = cbind(SampleID = rownames(TotalTTX), TotalTTX)
GroupTotalTTX = cbind(SampleID = rownames(GroupTotalTTX), GroupTotalTTX)
taricha_Group_TTX = cbind(SampleID = rownames(group_Totalttx), group_Totalttx)

rownames(TotalTTX) <- NULL
rownames(GroupTotalTTX) <- NULL
rownames(taricha_Group_TTX) <- NULL

# cat("Enter rarefying depth:");
# rardepth <- readLines("stdin",n=1);
# cat("You entered")
# str(rardepth);
# cat( "\n" )

rardepth <- 5200
rardepth <- as.numeric(rardepth)
tt_rardepth <- 3800

TotalTTX_data = TotalTTX %>%
  mutate(Propor_TotalTTX = TotalTTX/rardepth)

Group_TTX_data = GroupTotalTTX %>%
  mutate(Group_Propor_TotalTTX = GroupTotalTTX/rardepth)

tt_total_TTX_data = taricha_Group_TTX %>%
  mutate(taricha_proport_totalTTX = tt_group_ttx/tt_rardepth)

# cat("Enter Metadata filepath:");
# Metapath <- readLines("stdin",n=1);
# cat("You entered")
# str(Metapath);
# cat( "\n" )

Metapath <- read_tsv('new_newt_metadata.tsv')
names(Metapath)[1] <- 'SampleID'
Metapath <- str_trim(Metapath)
group_metadata <- Metapath

tt_metadata <- read.csv('Ttorosa_mapping_Alpha_AntiBd_TTX.csv')
tt_metadata <- str_trim(tt_metadata)

Metadata = readr::read_tsv(Metapath)

Metadata_TTX = TotalTTX_data %>%
  left_join(Metadata, by = "SampleID")

group_metadata = Group_TTX_data %>%
  left_join(Metadata, by = 'SampleID')

tt_metadata = tt_total_TTX_data %>%
  left_join(tt_metadata, by = 'SampleID')

TTX_tablePA = qiime2R::read_qza("TTX_Matches_99/PresAbs_table_clustered.qza")
group_tablePA = qiime2R::read_qza('group_TTX_matches/group_presabs_clustered.qza')
tt_group_tablePA = qiime2R::read_qza('taricha_matches/taricha_Grouped_PresAbs_clustered_table.qza')

TTX_tablePA = as.data.frame(TTX_tablePA$data)
group_tablePA = as.data.frame(group_tablePA$data)
tt_group_tablePA = as.data.frame(tt_group_tablePA$data)

TTX_Richness = TTX_tablePA %>%
  select_if(is.numeric) %>%
  map_dbl(sum)
group_TTX_richness = group_tablePA %>%
  select_if(is.numeric) %>%
  map_dbl(sum)
TTXRichness = tt_group_tablePA %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

TTX_Richness = as.data.frame(TTX_Richness)
group_TTX_richness = as.data.frame(group_TTX_richness)
TTXRichness = as.data.frame(TTXRichness)

TTX_Richness = cbind(SampleID = rownames(TTX_Richness), TTX_Richness)
group_TTX_richness = cbind(SampleID = rownames(group_TTX_richness), group_TTX_richness)
TTXRichness = cbind(SampleID = rownames(TTXRichness), TTXRichness)

rownames(TTX_Richness) <- NULL
rownames(group_TTX_richness) <- NULL
rownames(TTXRichness) <- NULL

Metadata_TTX_Richness = MetadatNULLMetadata_TTX_Richness = Metadata_TTX %>%
  left_join(TTX_Richness, by = "SampleID")
group_TTX_richness = group_metadata %>%
  left_join(group_TTX_richness, by = 'SampleID')
TTXRichness = tt_metadata %>%
  left_join(TTXRichness, by = 'SampleID')

Metadata_TTX_Richness = Metadata_TTX_Richness %>%
  select(-TotalTTX,everything())  %>%
  select(-Propor_TotalTTX, everything())

readr::write_tsv(Metadata_TTX_Richness,"Metadata_TTX_Predictions.txt")
write.csv(Metadata_TTX_Richness, 'ttx_richness_metadata.csv', row.names = FALSE)

group_metadata = group_metadata %>%
  select(-GroupTotalTTX,everything())  %>%
  select(-Group_Propor_TotalTTX, everything())
write.csv(group_metadata, 'group_ttx_richness_metadata.csv', row.names = FALSE)
  
Metadat_richness_ttx = MetadataNULLmetadata_TTX_Richness = tt_metadata %>%
  select(-tt_total_TTX_data, everything()) %>%
  select(-taricha_proport_totalTTX, everything())
final_tt <- TTXRichness[,-c(4:34)]
write.csv(final_tt, 'tt_group_ttx_richness_metadata.csv', row.names = FALSE)



tg_amphibac_table = qiime2R::read_qza('tg_amphibac/clustered_table.qza')
tg_amphibac_table = as.data.frame(tg_amphibac_table$data)
totalAntifungal = tg_amphibac_table %>%
  select_if(is.numeric) %>%
  map_dbl(sum)
TotalAntifungal = as.data.frame(totalAntifungal)
TotalAntifungal = cbind(SampleID = rownames(TotalAntifungal), TotalAntifungal)
rownames(TotalAntifungal) <- NULL
antifungal_metadata <- read.table('t_granulosa/Metadata_TTX_Predictions.txt', header = T)
TotalAntifungal <- antifungal_metadata %>%
  left_join(TotalAntifungal, by = 'SampleID')
ProportionAntifungal = TotalAntifungal %>%
  mutate(ProportionAntifungal = totalAntifungal/2200)
# ANtifungalRichness = TotalAntifungal %>%
#   select(-TotalAntifungal,everything())  %>%
#   select(-ProportionAntifungal, everything())
antifungal_presab <- read_qza('tg_amphibac/tg_amphibac_Presence_absence.qza')
antifungal_presab <- as.data.frame(antifungal_presab$data)
Antifungal_Richness <- antifungal_presab %>%
  select_if(is.numeric) %>%
  map_dbl(sum)
Antifungal_Richness <- as.data.frame(Antifungal_Richness)
AntifungalRichness = cbind(SampleID = rownames(Antifungal_Richness), Antifungal_Richness)
AntifungalRichness <- TotalAntifungal %>%
  left_join(AntifungalRichness, by = 'SampleID')
write.csv(AntifungalRichness, 'amphibac_ttx_variables.csv', row.names = F)

