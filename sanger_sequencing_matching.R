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

setwd('/home/toast/Documents/Thesis_Part_2/Newt_TTX_Project/sanger/vsearch/')

sanger_table_group = qiime2R::read_qza("sanger_matches_99id/clustered_table.qza")
sanger_table_group = as.data.frame(sanger_table_group$data)

sum_sequences = sanger_table_group %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

sum_sequences = as.data.frame(sum_sequences)

sum_sequences = cbind(Indiv_ID = rownames(sum_sequences), sum_sequences)

rownames(sum_sequences) <- NULL

sum_sequences = sum_sequences[order(sum_sequences$Indiv_ID),]

sum_sequences = as.data.frame(sum_sequences)

sum_data = sum_sequences %>%
  mutate(relative_abundance = (sum_sequences/2200)*100)


pres_ab = qiime2R::read_qza("sanger_matches_99id/sanger_PresAbs_clustered_table.qza")

pres_ab = as.data.frame(pres_ab$data)

rep_samp_richness = pres_ab %>%
  select_if(is.numeric) %>%
  map_dbl(sum)

rep_samp_richness = as.data.frame(rep_samp_richness)

sum_data_2 = cbind(sum_data, rep_samp_richness)
colnames(sum_data_2)[1] = 'SampleID'

#TMC02B
sanger_tab2 <- sanger_table_group
tab <- t(sanger_tab2)
tab <- tab[-84,]
tab <- as.data.frame(tab)
write.csv(tab, 'sanger_results_table.csv', row.names = T)
sanger_tab2$ASV <- rownames(sanger_tab2)
tmc_dat <- sanger_tab2 %>%
  filter(ASV %in% 'TMC02B')

tmc_dat2 <- t(sanger_table_group)
tmc_dat2 <- data.frame(tmc_dat2)
tmc_dat2$SampleID <- rownames(tmc_dat2)
tmc_dat2 <- tmc_dat2[-84,]
tmc_dat2 <- as.numeric(tmc_dat2)
tmc_dat2$TMC02B <- as.numeric(tmc_dat2$TMC02B)

sum_data_3 <- cbind(sum_data_2, tmc_dat2)
sum_data_3 <- sum_data_3[,-6]
sum_data_3$TMC02B <- as.numeric(sum_data_3$TMC02B)

write.csv(sum_data_2, 
          'sanger_sequence_results.csv', 
          row.names = FALSE)
metadata_2 <- read.csv('2200_full_metadata.csv', header = T)

new_dat <- sum_data_3 %>%
  left_join(metadata_2, by = 'SampleID')
percentage <- new_dat %>%
  mutate((sum_sequences/2200)*100)

new_dat <- new_dat %>%
  mutate(tmc_relab, (TMC02B/2200)*100)
colnames(new_dat)[38] <- 'tmc_relab'

#_______exploring_________
max(new_dat$relative_abundance)

rep_kw <- kruskal.test(new_dat$rep_samp_richness, new_dat$Location)
rep_kw

tmc_kw <- kruskal.test(new_dat$tmc_relab, new_dat$Location)
tmc_kw

tmc_lm <- lm(new_dat$TTX_mg ~ new_dat$tmc_relab)
summary(tmc_lm)

reads_kw <- kruskal.test(new_dat$relative_abundance, new_dat$Location)
reads_kw

ttx_rep_kw <- kruskal.test(new_dat$rep_samp_richness, new_dat$TTX_Real)
ttx_rep_kw
ttx_reL_kw <- kruskal.test(new_dat$relative_abundance, new_dat$TTX_Real)
ttx_reL_kw

ttx_lm_rep <- lm(new_dat$rep_samp_richness ~ new_dat$TTX_mg)
summary(ttx_lm_rep)
ttx_lm_rel <- lm(new_dat$relative_abundance ~ new_dat$TTX_mg)
summary(ttx_lm_rel)

inten_lm <- lm(new_dat$Log_bd ~ new_dat$rep_samp_richness)
summary(inten_lm)
ttx_glm_rep <- glm(new_dat$infected ~ new_dat$rep_samp_richness)
summary(ttx_lm_rep)

ttx_glm_rel <- glm(new_dat$ infected ~ new_dat$relative_abundance)
summary(ttx_lm_rel)
rel_inten_lm <- lm(new_dat$Log_bd ~ new_dat$relative_abundance)
summary(rel_inten_lm)


#bd_pos_sites
bd_sites_data <- new_dat %>%
  filter(!Location %in% c('FCSP', 'LSSP'))

site_kw <- kruskal.test(bd_sites_data$relative_abundance, bd_sites_data$Location)
site_kw

site_rep_kw <- kruskal.test(bd_sites_data$rep_samp_richness, bd_sites_data$Location)
site_rep_kw


pos_inten_lm <- lm(bd_sites_data$Log_bd ~ bd_sites_data$rep_samp_richness)
summary(pos_inten_lm)
pos_ttx_glm_rep <- glm(bd_sites_data$infected ~ bd_sites_data$rep_samp_richness)
summary(pos_ttx_glm_rep)

pos_ttx_glm_rel <- glm(bd_sites_data$infected ~ bd_sites_data$relative_abundance)
summary(pos_ttx_glm_rel)
pos_rel_inten_lm <- lm(bd_sites_data$Log_bd ~ bd_sites_data$relative_abundance)
summary(pos_rel_inten_lm)

#bd_pos_individuals
bd_indiv_data <- new_dat %>%
  subset(infected == 1)

indiv_kw <- kruskal.test(bd_indiv_data$relative_abundance, bd_indiv_data$Location)
indiv_kw

indiv_site_rep_kw <- kruskal.test(bd_indiv_data$rep_samp_richness, bd_indiv_data$Location)
indiv_site_rep_kw

indiv_lm <- lm(bd_indiv_data$Log_bd ~ bd_indiv_data$relative_abundance)
summary(indiv_lm)
indiv_lm2 <- lm(bd_indiv_data$Log_bd ~ bd_indiv_data$rep_samp_richness)
summary(indiv_lm2)


