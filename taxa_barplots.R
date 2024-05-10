#clearing the environment
rm(list = ls())


#BiocManager
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
require('BiocManager')

#different load in technique <3
packages <- c(
  "usethis", "dplyr", "MASS", "ggplot2", "vegan", "readr", "plotrix",
  "reticulate", "effects", "remotes", "MuMIn", "ggtext", "car", "minpack.lm",
  "nlcor", "nlraa", "mgcv", "tidymv", "gratia", "usedist", "lme4", "vegan",
  "ape", "qiime2R", "tidyverse", "performance", "ggfortify", "plotly", "ggord",
  "ggpubr", "ggords", "sjPlot", "sjmisc", "emmeans", "rstatix", "patchwork",
  "gridExtra", "scales", "ggthemes", "egg"
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
library(Hmisc)
library(qiime2R)
library(ggords)
library(ggpcoa)
library(ggforce)
library(ggpubr)
library(gridExtra)
library(ecole)
library(RColorBrewer)
library(viridis)



setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/')
relfreq_1 <- read.table('relfreq.tsv', sep = '\t', header = F)
metadata <- read.csv('2200_full_metadata.csv', header = T)
colnames(relfreq_1) <- relfreq_1[1,]
relfreq_1 <- relfreq_1[-1,-1]

grouped <- relfreq_1 %>%
  group_by(Species) %>%
  slice(-1) %>%
  summarise_all(rel_abun = colSums(,2:48))

pasetel <- brewer.pal(9, 'Pastel1')

print(pasetel)
palette1 <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#E5D8BD", "#FDDAEC",
              "#F2F2F2", 'red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                 'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033',
                 "#8C96C6", "#78C679", "#41B6C4", "#FE9929", "#FD8D3C", "#018571", 
                 "#4DAC26","#008837", "#5E3C99", "#0571B0", "#404040", "#2C7BB6", 
                 "#1A9641","#FFFF99", "#E7298A", "#33A02C", "#DECBE4", "#F4CAE4", 
                 "#984EA3", "#E78AC3", "#FB8072", "#2171B5", "#238B45", "#88419D", 
                 "#2B8CBE", "#238B45", "#525252", "#D94701", "#D7301F", "#0570B0", 
                 "#02818A", "#CE1256", "#6A51A3", "#AE017E","#225EA8", "#CC4C02", 
                 "#E31A1C", "#FDC086", "#7570B3")
print(palette1)

tato_colors <- c("#FBB4AE", "#0571B0", "#008837", "#DECBE4", "#F39929", "#7570B3", "#D7301F", "#FFFF99",
                 "#F2F2F2", 'red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                 'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033',
                 "#8C96C6", "#78C679", "#41B6C4", "#FE9929", "#FD8D3C", "#018571", 
                 "#4DAC26","#008837", "#5E3C99", "#0571B0", "#404040", "#2C7BB6", 
                 "#1A9641","#FFFF99", "#E7298A", "#33A02C", "#DECBE4", "#F4CAE4", 
                 "#984EA3", "#E78AC3", "#FB8072", "#2171B5", "#238B45", "#88419D", 
                 "#2B8CBE", "#238B45", "#525252", "#D94701", "#D7301F", "#0570B0", 
                 "#02818A", "#CE1256", "#6A51A3", "#AE017E","#225EA8", "#CC4C02", 
                 "#E31A1C", "#FDC086", "#7570B3")
tagr_colors <- c("#FBB4AE", "#008837", "#FE9929", "#7570B3", "#DECBE4", "#FFFF99", "#0571B0", "#D7301F",
                 "#F2F2F2", 'red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                 'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033',
                 "#8C96C6", "#78C679", "#41B6C4", "#FE9929", "#FD8D3C", "#018571", 
                 "#4DAC26","#008837", "#5E3C99", "#0571B0", "#404040", "#2C7BB6", 
                 "#1A9641","#FFFF99", "#E7298A", "#33A02C", "#DECBE4", "#F4CAE4", 
                 "#984EA3", "#E78AC3", "#FB8072", "#2171B5", "#238B45", "#88419D", 
                 "#2B8CBE", "#238B45", "#525252", "#D94701", "#D7301F", "#0570B0", 
                 "#02818A", "#CE1256", "#6A51A3", "#AE017E","#225EA8", "#CC4C02", 
                 "#E31A1C", "#FDC086", "#7570B3")

gamma <- #fbb4ae
bacterioda <- #008837
alpha <- #fe9929
verrunco <- #7570b3
actino
clostridia
bacilli
gracil

new_data <- read.csv('taxa_barplot.csv', header = T, sep = '')
col_sum <- sum(new_data$Frequency)
new_data$rel_abundance <- (new_data$Frequency/col_sum)*100
new_data$Class <- factor(new_data$Class, 
                         levels = new_data$Class[order(new_data$rel_abundance, decreasing = TRUE)])
final_tato <- new_data %>%
  filter(Class %in% unique(Class[1:8]))
taxa_plot <- final_tato %>%
  ggplot(aes(x = 1, y = rel_abundance, fill = Class))+ 
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = tato_colors)+
  ylab('Relative Frequency')+
  xlab('Taricha torosa')+
  theme_minimal()+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(face = 'italic',
                                size = 15),
    axis.title.y = element_text(face = 'bold',
                                size = 15),
    legend.title = element_text(face = 'bold',
                                size = 20)
  )+
  scale_x_continuous(breaks = NULL)+
  labs(tag = '')
taxa_plot

tagr_data <- read.csv('tagr_taxa_barplot.csv', header = T, sep = '')
tagr_col_sum <-sum(tagr_data$Frequency)
tagr_data$rel_abundance <- (tagr_data$Frequency/tagr_col_sum)*100
tagr_data <- tagr_data[-21,]

tagr_data <- tagr_data |>
  separate(Taxonomy, into = c('Domain', 'Phylum', 'Class'), 
                      sep = ';') |>
  mutate(Class = gsub('c__', '', Class)) |>
  select(-c(1,2)) |>
  filter(Class %in% unique(Class[1:8]))

tagr_data$Class <- factor(tagr_data$Class, 
                          levels = tagr_data$Class[order(tagr_data$rel_abundance, decreasing = TRUE)])


tagr_barplot <- tagr_data %>%
  ggplot(aes(x = 1, y = rel_abundance, fill = Class))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = tagr_colors)+
  ylab('Relative Frequency')+
  xlab('Taricha granulosa')+
  theme_minimal()+
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(face = 'italic',
                                size = 15),
    axis.title.y = element_text(face = 'bold',
                                size = 15),
    legend.title = element_text(face = 'bold',
                                size = 20)
  )+
  scale_x_continuous(breaks = NULL) +
  labs(tag = 'b)')
tagr_barplot


combined_taxonomy_barplot <- ggarrange(tagr_barplot, taxa_plot,
                                       common.legend = T, legend = 'right')

#______graveyard
# relfreq_1 <- relfreq_1[-30,]
# relfreq_2 <- t(relfreq_1)
# colnames(relfreq_2) <- relfreq_2[1,]
# colnames(relfreq_2)[1] <- 'SampleID'
# relfreq_2 <- relfreq_2[-1,]
# rownames(relfreq_2) <- relfreq_2[,1]
# relfreq_2 <- relfreq_2[-119,]
# relfreq_2 <- as.data.frame(relfreq_2)
# 
# trial_barplot <- taxa_barplot(relfreq_2, metadata, 'SampleID')
