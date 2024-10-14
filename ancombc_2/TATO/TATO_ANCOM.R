rm(list=ls())


library(tidyverse)
library(readr)
library(ggplot2)
library(car)
library(MuMIn)
library(vegan)
library(effects)
library(qiime2R)
set.seed(001)




meta.data <- read.csv('ancombc_2/Ttorosa_mapping_Alpha_AntiBd_TTX.csv', sep =',') 
#Taxonomy ####
taxonomy <- as.data.frame(read_qza('ancombc_2/classification.qza')$data)
taxonomy <- taxonomy %>% 
  rename(id = 'Feature.ID')
taxa <- taxonomy %>% 
  mutate(taxonomy_parts = str_split(Taxon, "; ")) %>% 
  unnest(taxonomy_parts) %>% 
  separate(taxonomy_parts, into =c('Rank', 'Value'), sep = '__') %>% 
  pivot_wider(names_from = Rank, values_from = Value) %>% 
  rename(
    Domain = d, Phylum = p, Class = c,
    Order = o, Family = f, Genus = g,
    Species = s
  ) %>% 
  select(-Taxon) %>% 
  select("id", 'Family')

#differentials ####
differentails_p_val <- read.table('ancombc_2/ancom/exp_difs/p_val_slice.csv', header = T, sep = ',')%>% 
  rename(intercept = 'X.Intercept.') %>% 
  rename(p_value = 'TTX_WholeLow')

qval <- read.table('ancombc_2/ancom/exp_difs/q_val_slice.csv', header = T, sep = ',')  %>% 
  rename(q_value = 'TTX_WholeLow')

matched <- differentails_p_val %>% 
  left_join(qval, by = 'id') %>% 
  filter(q_value < 0.05)

#there are no individual taxa that have a q_value less than 0.05, meaning there are no significantly abundant taxa and there is no need to continue to go any further with this analysis.
