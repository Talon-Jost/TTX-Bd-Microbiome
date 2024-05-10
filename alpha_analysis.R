
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
library(ecole) > #remotes::install_github('phytomosaic/ecole')
  
  
setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/')

metadata <- read.csv('2200_full_metadata.csv', header = TRUE)
tg_data2 <- subset(metadata, Species == 'Taricha_granulosa')
tt_data2 <- subset(metadata, Species == 'Taricha_torosa')

species_list <- list(metadata = metadata, tg_data2 = tg_data2, tt_data2 = tt_data2)
species_names <- c("metadata", "tg_data2", "tt_data2")
dataframe_names <- c("metadata", "tg_data2", "tt_data2")  # Vector to store dataframe names

alpha_div <- c('pielou_evenness', 'shannon_entropy', 'observed_features','faith_pd', 
               'TTX_mg', 'mass_g', 'Log_bd')
variables <- c('TTX_mg', 'Log_bd', 'mass_g',
               'pielou_evenness', 'shannon_entropy', 'observed_features', 'faith_pd')

spear_results_df <- data.frame(
  matrix = character(),
  predictor = character(),
  response = character(),
  pseudo_F = numeric(),
  p_value = numeric(),
  dataframe_name = character(),  # Change column name to dataframe_name
  species_name = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(species_names)) {
  species_name <- species_names[i]
  dataframe_name <- dataframe_names[i]
  species_data <- species_list[[species_name]]
  for (predictor in alpha_div) {
    for (var in variables) {
      spear_result <- cor.test(species_data[[predictor]], species_data[[var]], method = "spearman")
      spear_results_df <- rbind(spear_results_df, data.frame(
        matrix = species_name,
        predictor = predictor,
        response = var,
        pseudo_F = spear_result$statistic,
        p_value = spear_result$p.value,
        dataframe_name = dataframe_name,  # Assign the dataframe name
        species_name = unique(species_data$Species)
      ))
    }
  }
}

write.csv(spear_results_df, file = "spear_results.csv", row.names = F)

cor.test(tg_data2$TTX_mg, tg_data2$faith_pd, method = 'spearman')
cor.test(tg_data$TTX_mg, tg_data$shannon_entropy, method = 'spearman')
cor.test(tg_data2$TTX_mg, tg_data2$observed_features, method = 'spearman')
cor.test(tg_data$TTX_mg, tg_data$pielou_evenness, method = 'spearman')

cor.test(tg_data$Log_bd, tg_data$faith_pd, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$shannon_entropy, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$observed_features, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$pielou_evenness, method = 'spearman')


cor.test(tg_data$mass_g, tg_data$faith_pd, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$shannon_entropy, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$observed_features, method = 'spearman')
cor.test(tg_data$Log_bd, tg_data$pielou_evenness, method = 'spearman')

cor.test(tt_data2$TTX_mg, tt_data2$faith_pd, method = 'spearman')
cor.test(tt_data$TTX_mg, tt_data$shannon_entropy, method = 'spearman')
cor.test(tt_data$TTX_mg, tt_data$observed_features, method = 'spearman')
cor.test(tt_data$TTX_mg, tt_data$pielou_evenness, method = 'spearman')

cor.test(tt_data$Log_bd, tt_data$faith_pd, method = 'spearman')
cor.test(tt_data$Log_bd, tt_data$shannon_entropy, method = 'spearman')
cor.test(tt_data$Log_bd, tt_data$observed_features, method = 'spearman')
cor.test(tt_data$Log_bd, tt_data$pielou_evenness, method = 'spearman')

cor.test(tt_data$mass_g, tt_data$faith_pd, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$shannon_entropy, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$observed_features, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$pielou_evenness, method = 'spearman')

new_lm <- lm(tt_data2$faith_pd ~ tt_data2$Log_bd)
summary(new_lm)
tt_data23 <- tt_data2 %>%
  filter(!Location %in% c('Bolinger', 'Crocker'))
plot(tt_data23$Log_bd, tt_data23$faith_pd)
plot(tt_data23$TTX_mg, tt_data23$faith_pd)


#_____-bd___pos____________
alpha_div <- c('Log_bd')
variables <- c('pielou_evenness', 'shannon_entropy', 'observed_features', 'faith_pd')

tg_pos <- tg_data2 %>%
  subset(infected == 1)

spear_results_df <- data.frame(
  matrix = character(),
  predictor = character(),
  response = character(),
  pseudo_F = numeric(),
  p_value = numeric(),
  dataframe_name = character(),
  stringsAsFactors = FALSE,
  rho = numeric()
)

for (predictor in alpha_div) {
  for (var in variables) {
    spear_result <- cor.test(alpha_div[[predictor]], variables[[var]], method = "spearman")
    spear_results_df <- rbind(spear_results_df, data.frame(
      predictor = predictor,
      response = var,
      pseudo_F = spear_result$statistic,
      p_value = spear_result$p.value,
      rho = spear_result$estimate,
    ))
  }
}

write.csv(spear_results_df, file = "spear_results.csv", row.names = FALSE)

cor.test(tg_pos$Log_bd, tg_pos$pielou_evenness, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$shannon_entropy, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$observed_features, method = 'spearman')
cor.test(tt_data$mass_g, tt_data$pielou_evenness, method = 'spearman')

kruskal.test(tg_data2$faith_pd, tg_data2$TTX_Real)
cor.test(tg_data2$faith_pd, tg_data2$TTX_mg, method = 'spearman')
cor.test(tt_data2$faith_pd, tt_data2$TTX_mg, method = 'spearman')

new_lin <- lm(tg_data2$TTX_RelativeAbundance ~ tg_data2$TTX_mg)
summary(new_lin)
