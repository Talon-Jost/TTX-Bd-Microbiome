#rm(list=ls())

library(usethis)
library(dplyr)
library(MASS)
library(ggplot2)
library(lme4) 
library(nlme)
library(vegan)
library(readr)
library(plotrix)
library(reticulate)
library(effects)
library(remotes)
library(MuMIn)
library(ggtext)
library(car)
library(MuMIn)
library(minpack.lm)
#remotes::install_github('ProcessMiner/nlcor')
library(nlcor)
#remotes::install_github('femiguez/nlraa')
library(nlraa)
library(mgcv)
library(tidymv)
#remotes::install_github("gavinsimpson/gratia")
library(gratia)
library(usedist) #tools for distance matrix
library(lme4) #statistics package for more complex models
library(vegan) #used for microbiome statistics
#library(pairwiseAdonis) #also aids in microbiome stats
library(ape) #transforms distance matrix for plotting
library(qiime2R) #imports QIIME files for use in R
library(tidyverse) #combination of useful packages that aid and perform different functions
library(performance)
library(ggfortify)
library(plotly)
#library(mdthemes)

# Enable the r-universe repo
#options(repos = c(
#  fawda123 = 'https://fawda123.r-universe.dev',
#  CRAN = 'https://cloud.r-project.org'))

# Install ggord
#install.packages('ggord')
library(ggord)
library(ggpubr)

#remotes::install_github("wdy91617/ggords")
library(ggords)
library(sjPlot)
library(sjmisc)
library(corrr)
library(ggcorrplot)
library(FactoMineR)

#loading in!
setwd("~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis")
df1 <- read.csv('TBJALH_TTX_DATA_CLEANED.csv', header = TRUE)
df2 <- read.csv('phylogenetic_torosa_granulosa_metadata.csv', header = T)
df <- subset(df2, Species == 'Taricha_granulosa')
rownames(df) <- df[,1]
weighted_unifrac_tg <- read_qza('microbiome/tg/weighted_unifrac/weighted_unifrac_distance_matrix.qza')$data
wu_dm <- as.matrix(weighted_unifrac_tg)

#df exploration
colSums(is.na(df))
#looks like we have no NAs, neat.
#making it **numerical**
df_numerical <- select_if(df, is.numeric)
rownames(df_numerical)
merged <- merge(wu_dm, df_numerical)
merged <- t(merged)

ad <- adonis2(wu_dm ~ df$TTX.mg.)
df2 <- adonis2(df$TTX.mg. ~ wu_dm)
df2
ad
an <- anosim(wu_dm ~ df$TTX.mg.)
