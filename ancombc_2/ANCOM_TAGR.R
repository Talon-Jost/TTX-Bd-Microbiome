rm(list=ls())


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
library(Hmisc)
library(qiime2R)
library(ggords)
library(ggpcoa)
library(ggforce)
library(ggpubr)
library(gridExtra)
library(ecole) #remotes::install_github('phytomosaic/ecole')
library(spearmanCI)
library(rstatix)
library(sjstats)
library(dunn.test)
library(multcompView)
library(FSA)
library(rcompanion)
library(broom.mixed)
library(lmerTest)
library(sjPlot)
library(qiime2R)
set.seed(12345)

location_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                     'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')
tg_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan')
tt_colors <- c('salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')

setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/')


meta.data <- read.csv('newt_metadata.tsv', sep ='\t') %>% 
  rename(SampleID = 'X.SampleID') %>% 
  filter(!SampleID %in% c('Pos',
                          'Neg',
                          'FCSP-05'))

#Taxonomy ####
taxonomy <- as.data.frame(read_qza('taxonomy/classification.qza')$data)
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
differentails <- read.table('ancombc_2/exported_difs/p_val_slice.csv', header = T, sep = ',') %>% 
  rename(intercept = 'X.Intercept.') %>% 
  rename(p_value = 'TTX_RealLow')

qval <- read.table('ancombc_2/exported_difs/q_val_slice.csv', header = T, sep = ',')  %>% 
  rename(q_value = 'TTX_RealLow')

matched <- differentails %>% 
  left_join(qval, by = 'id') %>% 
  filter(q_value < 0.05)

lfc <- read.csv('ancombc_2/exported_difs/lfc_slice.csv', header = T, sep = ',')

full_match <- matched %>% 
  left_join(taxa, by = 'id')

#you ended wanting to join all the datasets. all the folders in exported difs and differentails are the same. just go join them because you need the beta and standard error to make the plot
