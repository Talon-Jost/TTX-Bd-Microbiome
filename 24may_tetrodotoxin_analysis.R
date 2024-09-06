#clearing the environment
rm(list = ls())

#BiocManager
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
require('BiocManager')

#different load in technique <3
packages <- c(
  "tidyverse", "readr", "ggplot2", "car", "MuMIn", "vegan",
  "effects", "ggtext", "lme4", "mgcv", "Hmisc",
  "qiime2R", "ggords", "ggpcoa", "ggforce", "ggpubr",
  "gridExtra", "ecole", "spearmanCI", "rstatix", "sjstats",
  "dunn.test", "multcompView", "FSA", "rcompanion", "LDM"
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
set.seed(000)

location_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                     'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')
tg_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan')
tt_colors <- c('salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')

setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/')

metadata <- read.csv('2200_full_metadata.csv', header = TRUE)
metadata$TTX_Real <- ifelse(metadata$TTX_mg > 1, 'High', 
                     ifelse(metadata$TTX_mg < 1, 'Low', 'Low'))
metadata_2 <- metadata %>%
  filter(!Location %in% c('Bolinger', 'Crocker')) %>%
  filter(!SampleID %in% c('FCSP-05', 'Tt.10B2'))
metadata_3 <- metadata %>%
  filter(!Location %in% c('Bolinger', 'Crocker')) %>%
  filter(!SampleID %in% 'Tt.10B2')
tg_data2 <- subset(metadata, Species == 'Taricha_granulosa')
tt_data2 <- subset(metadata, Species == 'Taricha_torosa')
tt_data2 <- tt_data2 %>%
  filter(!Location %in% c('Bolinger', 'Crocker') & !TTX_mg > 2)
fcsp_neg <- tg_data2 %>%
  filter(!SampleID %in% 'FCSP-05')

#TTX between species 
species_ttx <- kruskal.test(metadata_2$TTX_mg, metadata_2$Species)
species_ttx

#see if TTX varies across locations for TAGR
tg_ttx_kruskal <- kruskal.test(tg_data2$TTX_mg, tg_data2$Location) 
tg_ttx_kruskal
# big difference

#which sites are different?
dunn1 <- dunn.test(tg_data2$TTX_mg, tg_data2$Location, method = 'bonferroni')
cld <- cldList(comparison = dunn1$comparisons,
               p.value = dunn1$P.adjusted,
               threshold = 0.05)
names(cld)[1] <- 'Location'
cld 

dt1 <- fcsp_neg %>%
  group_by(Location) %>%
  summarise(ttx = max(TTX_mg))
dt1$Letter <- cld$Letter
dt1

#plot and add significance
tg_ttx_location <- fcsp_neg %>%
  ggplot(aes(x = Location, y = TTX_mg, fill = Location)) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = tg_colors) +
  scale_color_manual(values = tg_colors) +
  ylab('Tetrodotoxin Concentration (mg)') +
  xlab('') +
  ggtitle('Taricha granulosa') +
  scale_x_discrete(labels = c('FCSP', "GRP","LSSP", "SLNF", "TMC", 'WBC')) +
  theme(
    plot.title = element_text(face = 'italic', hjust = 0.5, size = 20),
    axis.text.x = element_text(face = 'bold', size = 14),
    axis.text.y = element_text(face = 'bold', size = 14),
    axis.title = element_text(face = 'bold', size = 14),
    axis.line = element_line(),
    panel.grid = element_blank()
  ) +
  guides(fill = FALSE) +
  geom_text(data = dt1, aes(label = Letter, x = Location, y = ttx + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  ylim(0, 8)+
  labs(tag = 'a)')
tg_ttx_location

#differences between sexes?
tg_sex_kw <- kruskal.test(tg_data2$TTX_mg, tg_data2$sex)
tg_sex_kw
#not different between sexes


#TATO time!
tt_ttx_loc_kw <- kruskal.test(tt_data2$TTX_mg, tt_data2$Location)
tt_ttx_loc_kw

#dunn test time
dunn2 <- dunn.test(tt_data2$TTX_mg, tt_data2$Location, method = 'bonferroni')
cld2 <- cldList(comparison = dunn2$comparisons,
               p.value = dunn2$P.adjusted,
               threshold = 0.05)
names(cld2)[1] <- 'Location'
cld2 

dt2 <- tt_data2 %>%
  group_by(Location) %>%
  summarise(ttx = max(TTX_mg))
dt2$Letter <- cld2$Letter
dt2

tt_ttx_location <- tt_data2 %>%
  ggplot(aes(x = Location, y = TTX_mg, fill = Location))+
  geom_boxplot()+
  theme_minimal()+
  scale_fill_manual(values = tt_colors)+
  scale_color_manual(values = tt_colors)+
  ylab('Tetrodotoxin Concentration (mg)')+
  xlab('')+
  ggtitle('Taricha torosa')+
  scale_x_discrete(labels = c("Cold Creek", "Grunstock",
                              "Madonna","Muir"))+
  theme(
    plot.title = element_text(face = 'italic',
                              hjust = 0.5,
                              size = 20),
    axis.text.x = element_text(face = 'bold',
                               size = 14),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    axis.title = element_text(face = 'bold',
                              size = 14),
    axis.line = element_line(),
    panel.grid = element_blank()
  )+
  guides(fill = FALSE)+
  geom_text(data = dt2, aes(label = Letter, x = Location, y = ttx + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag = 'b)')
tt_ttx_location

#combine the plots
ttx_location_arranged <- ggarrange(tg_ttx_location, tt_ttx_location)
ttx_location_arranged

#TATO sex kw
tt_sex_kw <- kruskal.test(tt_data2$TTX_mg, tt_data2$sex)
tt_sex_kw


#Bd stuff
Bd_pos_total <- metadata_3 %>%
  subset(infected == 1)
bd_tg_df <- Bd_pos_total %>%
  subset(Species == 'Taricha_granulosa')
tato_bd_df <- metadata_3 %>%
  subset(Species == 'Taricha_torosa')

#species bd
species_bd_kw <- kruskal.test(metadata_3$Log_bd, metadata_3$Species)
species_bd_kw

infected_bd_kw <- kruskal.test(Bd_pos_total$Log_bd, Bd_pos_total$Species)
infected_bd_kw

infected <- kruskal.test(metadata_3$infected, metadata_3$Species)
infected

#TAGR Bd intentisy
tg_bd_kw <- kruskal.test(tg_data2$Log_bd, tg_data2$Location)
tg_bd_kw

dunn_bd_kw <- dunn.test(tg_data2$Log_bd, tg_data2$Location, method = 'bonferroni')
cld_bd <- cldList(comparison = dunn_bd_kw$comparisons,
                p.value = dunn_bd_kw$P.adjusted,
                threshold = 0.05)
names(cld_bd)[1] <- 'Location'
cld_bd 

dt3 <- tg_data2 %>%
  group_by(Location) %>%
  summarise(log_bd = max(Log_bd))
dt3$Letter <- cld_bd$Letter
dt3


#TAGR infection status
tg_bd_kw_pres <- kruskal.test(tg_data2$infected, tg_data2$Location)
tg_bd_kw_pres

dunn_bd_kw_inten <- dunn.test(tg_data2$infected, tg_data2$Location, method = 'bonferroni')
cld_bd_infected <- cldList(comparison = dunn_bd_kw_inten$comparisons,
                  p.value = dunn_bd_kw_inten$P.adjusted,
                  threshold = 0.05)
names(cld_bd_infected)[1] <- 'Location'
cld_bd_infected

dt31 <- tg_data2 %>%
  group_by(Location) %>%
  summarise(infected = mean(infected)*100)
dt31$Letter <- cld_bd_infected$Letter
dt31


#TATO Bd intensity
tt_bd_kw <- kruskal.test(tt_data2$Log_bd, tt_data2$Location)
tt_bd_kw

tato_dunn_bd_kw <- dunn.test(tt_data2$Log_bd, tt_data2$Location, method = 'bonferroni')
cld_bd2 <- cldList(comparison = tato_dunn_bd_kw$comparisons,
                  p.value = tato_dunn_bd_kw$P.adjusted,
                  threshold = 0.05)
names(cld_bd2)[1] <- 'Location'
cld_bd2

dt4 <- tt_data2 %>%
  group_by(Location) %>%
  summarise(log_bd = max(Log_bd))
dt4$Letter <- cld_bd2$Letter
dt4

#TATO infection status
tt_bd_kw_infected <- kruskal.test(tt_data2$Log_bd, tt_data2$Location)
tt_bd_kw_infected

tt_dun_status <- dunn.test(tt_data2$infected, tt_data2$Location, method = 'bonferroni')
tt_status_cld <- cldList(comparison = tt_dun_status$comparisons,
                           p.value = tt_dun_status$P.adjusted,
                           threshold = 0.05)
names(tt_status_cld)[1] <- 'Location'
tt_status_cld

dt_status <- tt_data2 %>%
  group_by(Location) %>%
  summarise(infected = mean(infected)*100)
dt_status$Letter <- tt_status_cld$Letter
dt_status

#bd infected locations kw TAGR
loc_bd_pos_tg <- tg_data2 %>%
  filter(!Location %in% 'FCSP')

tg_loc <- kruskal.test(loc_bd_pos_tg$Log_bd, loc_bd_pos_tg$Location)
tg_loc

tg_loc_status <- kruskal.test(loc_bd_pos_tg$infected, loc_bd_pos_tg$Location)
tg_loc_status

#bd infected locations TATO
loc_tt <- tt_data2 %>%
  filter(!Location %in% c('Gunstock', 'Madonna'))

tt_loc_inten <- kruskal.test(loc_tt$Log_bd, loc_tt$Location)
tt_loc_inten

tt_status_bd <- kruskal.test(loc_tt$infected, loc_tt$Location)
tt_status_bd

#bd pos indivs
indiv_bd <- kruskal.test(bd_tg_df$Log_bd, bd_tg_df$Location)
indiv_bd

indiv_tt <- kruskal.test(tato_bd_df$Log_bd, tato_bd_df$Location)
indiv_tt

;#Bd infection plots______________________________________________________________________________________
bd_loc_plot2 <- tg_data2 %>%
  ggplot(aes(x = Location, y = Log_bd, fill = Location))+
  geom_boxplot()+
  theme_classic()+
  xlab('')+
  ylab('ITS Copies')+
  ggtitle('Taricha granulosa')+
  theme(axis.text.x = element_text(face = 'bold',
                                   size = 14),
        axis.text.y = element_text(size = 16,
                                   face = 'bold'),
        axis.title.y = element_text(face = 'bold',
                                    size = 14),
        plot.title = element_text(face = 'italic',
                                  size = 20,
                                  hjust = 0.5),
        panel.grid = element_blank())+
  scale_x_discrete(labels = c(
    'FCSP', "GRP",
    "LSSP","SLNF", 
    "TMC", 'WBC'))+
  guides(fill = FALSE)+
  scale_fill_manual(values = tg_colors)+
  coord_cartesian(ylim=c(0,NA))+
  expand_limits(y = 0)+
  geom_text(data = dt3, aes(label = Letter, x = Location, y = log_bd + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag='a)')
bd_loc_plot2

bd_loc_plot3 <- tt_data2 %>%
  ggplot(aes(x = Location, y = Log_bd, fill = Location))+
  geom_boxplot()+
  theme_classic()+
  xlab('')+
  ylab('ITS Copies')+
  ggtitle('Taricha torosa')+
  theme(axis.text.x = element_text(face = 'bold',
                                   size = 14),
        axis.text.y = element_text(size = 16,
                                   face = 'bold'),
        axis.title.y = element_text(face = 'bold',
                                    size = 14),
        plot.title = element_text(face = 'italic',
                                  size = 20,
                                  hjust = 0.5),
        panel.grid = element_blank())+
  scale_x_discrete(labels = c("Cold Creek","Grunstock","Madonna","Muir"))+
  guides(fill = FALSE)+
  scale_fill_manual(values = tt_colors)+
  coord_cartesian(ylim=c(0,NA))+
  expand_limits(y = 0)+
  geom_text(data = dt4, aes(label = Letter, x = Location, y = log_bd + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag = 'b)')
bd_loc_plot3

bd_loc_plot_joined <- ggarrange(bd_loc_plot2, bd_loc_plot3)
bd_loc_plot_joined

tg_prev_prop_plot <- tg_data2 %>%
  group_by(Location) %>%
  summarise(Proportion = mean(infected == 1) * 100)

tg_prev_plot <- ggplot(
  tg_prev_prop_plot, aes(x = Location, y = Proportion))+
  geom_bar(stat = 'identity', fill = '#F8766D')+
  theme_minimal()+
  xlab('')+
  ylab('Infection Prevalence (%)')+
  ggtitle('Taricha granulosa')+
  theme(
    plot.title = element_text(face = 'italic',
                              size = 20,
                              hjust = 0.5),
    axis.text.x = element_text(face = 'bold',
                               size = 14),
    axis.title.y = element_text(face = 'bold',
                                size = 14),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    axis.line = element_line(),
    panel.grid = element_blank())+
  geom_text(data = dt31, aes(label = Letter, x = Location, y = infected + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag = 'c)')

tg_prev_plot



tt_prev <- tt_data2 %>%
  group_by(Location) %>%
  summarise(Proportion = mean(infected == 1) * 100)

tt_prev_plot <- ggplot(
  tt_prev, aes(x = Location, y = Proportion))+
  geom_bar(stat = 'identity', fill = '#00BFC4')+
  theme_minimal()+
  xlab('')+
  ylab('Infection Prevalence (%)')+
  ggtitle('Taricha torosa')+
  theme(
    axis.text.x = element_text(face = 'bold',
                               size = 14),
    axis.title.y = element_text(face = 'bold',
                                size = 14),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    plot.title = element_text(face = 'italic',
                              size = 20,
                              hjust = 0.5),
    axis.line = element_line(),
    panel.grid = element_blank())+
  scale_x_discrete(labels = c("Cold Creek","Grunstock","Madonna","Muir"))+
  geom_text(data = dt_status, aes(label = Letter, x = Location, y = infected + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag = 'd)')
tt_prev_plot

prev_plots <- ggarrange(tg_prev_plot, tt_prev_plot)
prev_plots

tt_bd_prev_plot <- tt_data2 %>%
  ggplot(aes(x = Location, fill = infected)) +
  geom_bar(position = position_dodge()) +
  theme_minimal() +
  xlab('')+
  ylab('Infected Individuals')+
  labs(fill = 'Infection Status')+
  ggtitle('Bd Infection Prevalence')+
  theme(axis.text.x = element_text(face = 'bold',
                                   size = 16),
        axis.text.y = element_text(face = 'bold',
                                   size = 14),
        axis.title.x = element_text(face = 'bold',
                                    size = 14,
                                    angle = 45,
                                    hjust = 0.05),
        plot.title = element_text(face = 'bold',
                                  hjust = 0.5,
                                  size = 20),
        axis.title.y = element_text(face = 'bold',
                                    size = 14))+
  scale_x_discrete(labels = c('Bolinger', "Cold Creek",'Crocker',
                              "Grunstock","Madonna","Muir"))+
  geom_vline(xintercept = 0, 
             color = "black", 
             size = 0.5)+
  geom_hline(yintercept = 0, 
             color = "black", 
             size = 0.5)+
  labs(tag = 'b)')


multi_plots <- ggarrange(bd_loc_plot_joined, prev_plots,
                         nrow = 2)
multi_plots


#__________ttx-bd__________________________________________________________________________________

tg_ttx_bd <- glmer(infected ~ TTX_mg * mass_g * svl_mm + (1 | Location), data = tg_data2, family = 'binomial')
tg_ttx_bd2 <- glmer(infected ~ TTX_mg + (1 | Location), data = tg_data2, family = 'binomial')
summary(tg_ttx_bd)
anova(tg_ttx_bd)
plot(allEffects(tg_ttx_bd))
anova(tg_ttx_bd, tg_ttx_bd2, test = 'Chisq')

tg_sex_glmm <- glmer(infected ~ sex + (1 | Location), data = tg_data2, family = 'binomial')
summary(tg_sex_glmm)
anova(tg_sex_glmm)

tt_ttx_bd <- glmer(infected ~ TTX_mg + (1 | Location), data = tt_data2, family = 'binomial')
tt_ttx_bd2 <- glmer(infected ~ TTX_mg * mass_g + (1 | Location), data = tt_data2, family = 'binomial')
summary(tt_ttx_bd)
anova(tt_ttx_bd)
anova(tt_ttx_bd, tt_ttx_bd2, test = 'Chisq')

tt_sex_glmm <- glmer(infected ~ sex + (1 | Location), data = tt_data2, family = 'binomial')
summary(tt_sex_glmm)
anova(tt_sex_glmm)

#Plot
unique_locations <- unique(tg_data2$Location)
ttx_weight <- seq(0.05, 7, 0.01)
location_colors <- data.frame(Location = unique_locations, Color = tg_colors[1:length(unique_locations)])
newdata <- merge(tg_data2, location_colors, by = "Location")
newdata <- expand.grid(TTX_mg = ttx_weight, Location = unique_locations)
predicted_responses <- predict(tg_ttx_bd, newdata = newdata, type = 'response')
newdata$Predicted_Response <- predicted_responses


tg_glmm <- ggplot(newdata, aes(x = TTX_mg, y = Predicted_Response, 
                    group = Location, color = Location)) +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values = tg_colors) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = 'bold',
                               size = 16),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    axis.title.x = element_text(face = 'bold',
                                size = 14,
                                hjust = 0.5),
    plot.title = element_text(face = 'italic',
                              hjust = 0.5,
                              size = 20),
    axis.title.y = element_text(face = 'bold',
                                size = 14)
  ) +
  labs(title = "Taricha granulosa",
       x = "TTX Weight (mg)",
       y = "Infection Status")+  
  geom_vline(xintercept = -0.1, 
             color = "black", 
             linewidth = 0.5)+
  geom_hline(yintercept = 0, 
             color = "black", 
             linewidth = 0.5)+
  expand_limits(x = c(0, max(newdata$TTX_mg)), y = c(0, max(newdata$Predicted_Response)))+
  labs(tag = 'a)')

unique_locations_tt <- unique(tt_data2$Location)
tt_location_colors <- data.frame(Location = unique_locations_tt, Color = tt_colors[1:length(unique_locations_tt)])
tt_data2 <- merge(tt_data2, tt_location_colors, by = "Location")
tt_ttx_weight <- seq(0.002, 2, 0.01)
tt_newdata <- expand.grid(TTX_mg = tt_ttx_weight, Location = unique_locations_tt)
tt_predicted_responses <- predict(tt_ttx_bd, newdata = tt_newdata, type = 'response')
tt_newdata$Predicted_Response <- tt_predicted_responses


tt_plot_glmm <- ggplot(tt_newdata, aes(x = TTX_mg, y = Predicted_Response, 
                    group = Location, color = Location)) +
  geom_line() +
  theme_minimal() +
  scale_color_manual(values = tt_colors) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = 'bold',
                               size = 16),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    axis.title.x = element_text(face = 'bold',
                                size = 14,
                                hjust = 0.5),
    plot.title = element_text(face = 'italic',
                              hjust = 0.5,
                              size = 20),
    axis.title.y = element_text(face = 'bold',
                                size = 14)
  ) +
  labs(title = "Taricha torosa",
       x = "TTX Weight (mg)",
       y = "Infection Status")+  
  geom_vline(xintercept = -0.1, 
             color = "black", 
             linewidth = 0.5)+
  geom_hline(yintercept = 0, 
             color = "black", 
             linewidth = 0.5)+
  expand_limits(x = c(0, max(tt_newdata$TTX_mg)), y = c(0, max(tt_newdata$Predicted_Response)))+
  labs(tag = 'b)')

new_glmm <- ggarrange(tg_glmm, tt_plot_glmm)
new_glmm
#pos loc

tg_loc_glmm <- glmer(infected ~ TTX_mg + (1 | Location), data = loc_bd_pos_tg, family = 'binomial')
summary(tg_loc_glmm)
anova(tg_loc_glmm)

tt_loc_glmm <- glmer(infected ~ TTX_mg + (1 | Location), data = loc_tt, family = 'binomial')
summary(tt_loc_glmm)
anova(tt_loc_glmm)


#_________________________alpha div______________________
shapiro.test(tt_data2$faith_pd)
shapiro.test(tt_data2$observed_features)
tt_data2$log_asv <- log(tt_data2$observed_features + 1)
shapiro.test(tt_data2$log_asv)
shapiro.test(tt_data2$pielou_evenness)
shapiro.test(tt_data2$shannon_entropy)

tg_data2$log_faith <- log(tg_data2$faith_pd + 1)
shapiro.test(tg_data2$log_faith)
tg_data2$log_asv <- log(tg_data2$observed_features + 1)
shapiro.test(tg_data2$log_asv)
shapiro.test(tg_data2$log_pielou)
tg_data2$log_pielou <- log(tg_data2$pielou_evenness + 1)
shapiro.test(tg_data2$shannon_entropy)

tg_faith_kw <- kruskal.test(tg_data2$faith_pd, tg_data2$Location)
tg_faith_kw

tg_asv_kw <- kruskal.test(tg_data2$observed_features, tg_data2$Location)
tg_asv_kw

tg_shan_kw <- kruskal.test(tg_data2$shannon_entropy, tg_data2$Location)
tg_shan_kw

tg_pielou_kw <- kruskal.test(tg_data2$pielou_evenness, tg_data2$Location)
tg_pielou_kw

tg_observed_features_lm <- lmer(observed_features ~ TTX_mg + (1 | Location), data = tg_data2)
summary(tg_observed_features_lm)
anova(tg_observed_features_lm)

tg_shannon_entropy_lm <- lmer(shannon_entropy ~ TTX_mg + (1 | Location), data = tg_data2)
summary(tg_shannon_entropy_lm)
anova(tg_shannon_entropy_lm)

faith_ttx <- lmer(faith_pd ~ TTX_mg + (1 | Location), tg_data2)
summary(faith_ttx)
anova(faith_ttx)

tg_pielou_evenness_lm <- lmer(pielou_evenness ~ TTX_mg + (1 | Location), data = tg_data2)
summary(tg_pielou_evenness_lm)
anova(tg_pielou_evenness_lm)
pe <- tg_data2$pielou_evenness
lin_1 <- lm(pe ~ TTX_mg, data = tg_data2)
box_res_pielou <- boxCox(lin_1)
best_lambda_pe <- box_res_pielou$x[which.max(box_res_pielou$y)]
transformed_variable <- bcPower(pe, best_lambda_pe)
tg_data2$transformed_pielou <- transformed_variable
shapiro_test(tg_data2$transformed_pielou)
tg_pielou_evenness_lm2 <- lmer(pielou_evenness ~ TTX_mg + (1 | Location), data = tg_data2)
summary(tg_pielou_evenness_lm2)
anova(tg_pielou_evenness_lm2)

tt_observed_features_lm <- lmer(observed_features ~ TTX_mg + (1 | Location), data = tt_data2)
summary(tt_observed_features_lm)
anova(tt_observed_features_lm)
plot(tt_observed_features_lm)

tt_shannon_entropy_lm <- lmer(shannon_entropy ~ TTX_mg + (1 | Location), data = tt_data2)
summary(tt_shannon_entropy_lm)
anova(tt_shannon_entropy_lm)


tt_faith_lm <- lmer(faith_pd ~ TTX_mg + (1 | Location), data = tt_data2)
summary(tt_faith_lm)
anova(tt_faith_lm)
tt_faith_lm

tt_pielou_evenness_lm <- lmer(pielou_evenness ~ TTX_mg + (1 | Location), data = tt_data2)
summary(tt_pielou_evenness_lm)
anova(tt_pielou_evenness_lm)


#Bd______________
tg_observed_features_lm_bd <- lmer(observed_features ~ Log_bd + (1 | Location), data = tg_data2)
summary(tg_observed_features_lm_bd)
anova(tg_observed_features_lm_bd)

tg_shannon_entropy_lm_bd <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = tg_data2)
summary(tg_shannon_entropy_lm_bd)
anova(tg_shannon_entropy_lm_bd)

tg_faith_lm_bd <- lmer(faith_pd ~ Log_bd + (1 | Location), data = tg_data2)
summary(tg_faith_lm_bd)
anova(tg_faith_lm_bd)

tg_pielou_evenness_lm_bd <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = tg_data2)
summary(tg_pielou_evenness_lm_bd)
anova(tg_pielou_evenness_lm_bd)





tt_observed_features_lm_bd <- lmer(observed_features ~ Log_bd + (1 | Location), data = tt_data2)
summary(tt_observed_features_lm_bd)
anova(tt_observed_features_lm_bd)

tt_shannon_entropy_lbd <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = tt_data2)
summary(tt_shannon_entropy_lbd)
anova(tt_shannon_entropy_lbd)

tt_faith_lm_bd <- lmer(faith_pd ~ Log_bd + (1 | Location), data = tt_data2)
summary(tt_faith_lm_bd)
anova(tt_faith_lm_bd)

tt_pielou_evenness_lm_bd <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = tt_data2)
summary(tt_pielou_evenness_lm_bd)
anova(tt_pielou_evenness_lm_bd)

#Bd locations______________
tg_obs_bdloc <- lmer(observed_features ~ Log_bd + (1 | Location), data = loc_bd_pos_tg)
summary(tg_obs_bdloc)
anova(tg_obs_bdloc)

tg_shan_bdloc <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = loc_bd_pos_tg)
summary(tg_shan_bdloc)
anova(tg_shan_bdloc)

tg_faith_lm_bdloc <- lmer(faith_pd ~ Log_bd + (1 | Location), data = loc_bd_pos_tg)
summary(tg_faith_lm_bdloc)
anova(tg_faith_lm_bdloc)


tg_pielou_bdloc <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = loc_bd_pos_tg)
summary(tg_pielou_bdloc)
anova(tg_pielou_bdloc)




tt_df3 <- tt_data2 %>%
  filter(Location %in% c('Cold_Creek', 'Muir'))

tt_observed_features_lm_bd_bdloc <- lmer(observed_features ~ Log_bd + (1 | Location), data = tt_df3)
summary(tt_observed_features_lm_bd_bdloc)
anova(tt_observed_features_lm_bd_bdloc)

tt_shannon_entropy_lbd_bdloc <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = tt_df3)
summary(tt_shannon_entropy_lbd_bdloc)
anova(tt_shannon_entropy_lbd_bdloc)

tt_faith_lm_bd_bdloc <- lmer(faith_pd ~ Log_bd + (1 | Location), data = tt_df3)
summary(tt_faith_lm_bd_bdloc)
anova(tt_faith_lm_bd_bdloc)

tt_pielou_evenness_lm_bd_bdloc <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = tt_df3)
summary(tt_pielou_evenness_lm_bd)
anova(tt_pielou_evenness_lm_bd_bdloc)

#bd_pos indiv________________________________________
tg_data4 <- tg_data2 %>%
  subset(infected == 1)

tg_obs_indiv <- lmer(observed_features ~ Log_bd + (1 | Location), data = tg_data4)
summary(tg_obs_indiv)
anova(tg_obs_indiv)

tg_shan_indiv <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = tg_data4)
summary(tg_shan_indiv)
anova(tg_shan_indiv)

tg_faith_lm_indiv <- lmer(faith_pd ~ Log_bd + (1 | Location), data = tg_data4)
summary(tg_faith_lm_indiv)
anova(tg_faith_lm_indiv)

tg_pielou_indiv <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = tg_data4)
summary(tg_pielou_indiv)
anova(tg_pielou_indiv)


tt_df4 <- tt_data2 %>%
  subset(infected == 1)

tt_obs_indiv <- lmer(observed_features ~ Log_bd + (1 | Location), data = tt_df4)
summary(tt_obs_indiv)
anova(tt_obs_indiv)

tt_shannon_indiv <- lmer(shannon_entropy ~ Log_bd + (1 | Location), data = tt_df4)
summary(tt_shannon_indiv)
anova(tt_shannon_indiv)

tt_faith_indiv <- lmer(faith_pd ~ Log_bd + (1 | Location), data = tt_df4)
summary(tt_faith_indiv)
anova(tt_faith_indiv)

tt_pielou_indiv <- lmer(pielou_evenness ~ Log_bd + (1 | Location), data = tt_df4)
summary(tt_pielou_indiv)
anova(tt_pielou_indiv)

#Figures for this
labs_for_stuff <- c('Taricha_granulosa' = '',
                    'Taricha_torosa' = 'Taricha torosa')
tg_obs <- kruskal.test(fcsp_neg$observed_features, fcsp_neg$TTX_Real)
tg_faith <- kruskal.test(fcsp_neg$faith_pd, fcsp_neg$TTX_Real)
tg_faith 

faith_data <- na.omit(fcsp_neg$faith_pd)
mean(faith_data$faith_pd[faith_data$TTX_Real == 'High'])
alpha_boxplot <- fcsp_neg %>%
  group_by(Species) %>%
  filter(Species == 'Taricha_granulosa' & !is.na(TTX_Real) & !is.na(faith_pd)) %>%
  ggplot(aes(x = TTX_Real, y = faith_pd, fill = TTX_Real), na.rm = T) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  labs(
    x = '',
    y = "Faith's Phylogenetic Diversity",
    title = 'Taricha granulosa'
  )+
  theme(
    axis.text.y = element_text(face = 'bold', size = 20),
    axis.title.y = element_text(face = 'bold', size = 26),
    axis.text.x = element_text(face = 'bold', size = 26),
    strip.text = element_text(face = 'italic', size = 25),
    axis.line = element_line(),
    panel.grid = element_blank(),
    plot.title = element_text(face = 'italic',
                         size = 25,
                         hjust = 0.5)
  ) +
  scale_x_discrete(labels = c('High (TTX >1mg)', 'Low (TTX < 1mg)'))+
  facet_wrap(~Species, labeller = labeller(Species = labs_for_stuff)) +
  stat_compare_means(method = 'kruskal.test', label = "p.signif", label.sep = "") +
  annotate("text", x = 1.5, y = max(fcsp_neg$faith_pd, na.rm = TRUE) + 5, 
           label = paste("Kruskal-Wallis chi-squared =", round(tg_faith$statistic, 2), 
                         ", p-value =", format.pval(tg_faith$p.value, digits = 4)), size = 5, hjust = 0.5, fontface = "bold")+
  labs(tag = 'a)')

alpha_boxplot


tt_faith <- kruskal.test(tt_data2$faith_pd, tt_data2$TTX_Real)
tt_faith
tt_alpha_boxplot <- tt_data2 %>%
  group_by(Species) %>%
  filter(Species == 'Taricha_torosa' & !is.na(TTX_Real) & !is.na(faith_pd)) %>%
  ggplot(aes(x = TTX_Real, y = faith_pd, fill = TTX_Real), na.rm = T) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  xlab('') +
  ylab("Faith's Phylogenetic Diversity") +
  theme(
    axis.text.y = element_text(face = 'bold', size = 20),
    axis.title.y = element_text(face = 'bold', size = 26),
    axis.text.x = element_text(face = 'bold', size = 26),
    strip.text = element_text(face = 'italic', size = 25),
    axis.line = element_line(),
    plot.title = element_text(size = 25,
                              face = 'italic')
  ) +
  scale_x_discrete(labels = c('High (TTX >1mg)', 'Low (TTX < 1mg)'))+
  facet_wrap(~Species, labeller = labeller(Species = labs_for_stuff)) +
  stat_compare_means(method = 'kruskal.test', label = "p.signif", label.sep = "") +
  annotate("text", x = 1.5, y = max(tt_data2$faith_pd, na.rm = TRUE) + 5, 
           label = paste("Kruskal-Wallis chi-squared =", round(tt_faith$statistic, 2), 
                         ", p-value =", format.pval(tt_faith$p.value, digits = 4)), size = 5, hjust = 0.5, fontface = "bold")

tt_alpha_boxplot

faith_arranged <- ggarrange(alpha_boxplot, tt_alpha_boxplot)
faith_arranged

#################
tg_obs <- kruskal.test(fcsp_neg$observed_features, fcsp_neg$TTX_Real)
tg_obs
obs_data <- na.omit(fcsp_neg$observed_features)
mean(obs_data$tg_obs[tg_obs$TTX_Real == 'High'])

alpha_boxplot_obs <- fcsp_neg %>%
  group_by(Species) %>%
  filter(Species == 'Taricha_granulosa' & !is.na(TTX_Real) & !is.na(observed_features)) %>%
  ggplot(aes(x = TTX_Real, y = observed_features, fill = TTX_Real), na.rm = T) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  labs(
    x = '',
    y = "ASV Richness",
    title = 'Taricha granulosa'
  )+
  theme(
    axis.text.y = element_text(face = 'bold', size = 20),
    axis.title.y = element_text(face = 'bold', size = 26),
    axis.text.x = element_text(face = 'bold', size = 26),
    strip.text = element_text(face = 'italic', size = 25),
    axis.line = element_line(),
    panel.grid = element_blank(),
    plot.title = element_text(face = 'italic',
                              size = 25,
                              hjust = 0.5)
  ) +
  scale_x_discrete(labels = c('High (TTX >1mg)', 'Low (TTX < 1mg)'))+
  facet_wrap(~Species, labeller = labeller(Species = labs_for_stuff)) +
  stat_compare_means(method = 'kruskal.test', label = "p.signif", label.sep = "") +
  annotate("text", x = 1.5, y = max(fcsp_neg$observed_features, na.rm = TRUE) + 5, 
           label = paste("Kruskal-Wallis chi-squared =", round(tg_obs$statistic, 2), 
                         ", p-value =", format.pval(tg_obs$p.value, digits = 4)), size = 5, hjust = 0.5, fontface = "bold")+
  labs(tag = 'a)')

alpha_boxplot_obs


tt_obs <- kruskal.test(tt_data2$observed_features, tt_data2$TTX_Real)
tt_obs
tt_alpha_boxplot_obs <- tt_data2 %>%
  group_by(Species) %>%
  filter(Species == 'Taricha_torosa' & !is.na(TTX_Real) & !is.na(observed_features)) %>%
  ggplot(aes(x = TTX_Real, y = observed_features, fill = TTX_Real), na.rm = T) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal() +
  xlab('') +
  ylab("ASV Richness") +
  theme(
    axis.text.y = element_text(face = 'bold', size = 20),
    axis.title.y = element_text(face = 'bold', size = 26),
    axis.text.x = element_text(face = 'bold', size = 26),
    strip.text = element_text(face = 'italic', size = 25),
    axis.line = element_line(),
    plot.title = element_text(size = 25,
                              face = 'italic')
  ) +
  scale_x_discrete(labels = c('High (TTX >1mg)', 'Low (TTX < 1mg)'))+
  facet_wrap(~Species, labeller = labeller(Species = labs_for_stuff)) +
  stat_compare_means(method = 'kruskal.test', label = "p.signif", label.sep = "") +
  annotate("text", x = 1.5, y = max(tt_data2$observed_features, na.rm = TRUE) + 5, 
           label = paste("Kruskal-Wallis chi-squared =", round(tt_obs$statistic, 2), 
                         ", p-value =", format.pval(tt_obs$p.value, digits = 4)), size = 5, hjust = 0.5, fontface = "bold")

tt_alpha_boxplot_obs

obs_arranged <- ggarrange(alpha_boxplot_obs, tt_alpha_boxplot_obs)
obs_arranged

#_______mixed alpha

faith_tg <- lmer(faith_pd ~ TTX_Real + (1 | Location), data = tg_data2)
summary(faith_tg)
anova(faith_tg)
plot(faith_tg)

obs_tg <- lmer(observed_features ~ TTX_Real + (1 | Location), data = tg_data2)
summary(obs_tg)
anova(obs_tg)

shan_tg <- lmer(shannon_entropy ~ TTX_Real + (1 | Location), data = tg_data2)
summary(shan_tg)
anova(shan_tg)

piel_tg <- lmer(pielou_evenness ~ TTX_Real + (1 | Location), data = tg_data2)
summary(piel_tg)
anova(piel_tg)



faith_1 <- lmer(faith_pd ~ TTX_Real + (1 | Location), data = tt_data2)
summary(faith_1)
anova(faith_1)

obs_1 <- lmer(observed_features ~ TTX_Real + (1 | Location), data = tt_data2)
summary(obs_1)
anova(obs_1)

shan_tt <- lmer(shannon_entropy ~ TTX_Real + (1 | Location), data = tt_data2)
summary(shan_tt)
anova(shan_tt)


piel_tt <- lmer(pielou_evenness ~ TTX_Real + (1 | Location), data = tt_data2)
summary(piel_tt)
anova(piel_tt)

#sex mixed 
faith_tg_sex <- lmer(faith_pd ~ sex + (1 | Location), data = tg_data2)
summary(faith_tg_sex)
anova(faith_tg_sex)
boxplot(tg_data2$faith_pd ~ tg_data2$sex)

obs_tg_sex <- lmer(observed_features ~ sex + (1 | Location), data = tg_data2)
summary(obs_tg_sex)
anova(obs_tg_sex)

shan_tg_sex <- lmer(shannon_entropy ~ sex + (1 | Location), data = tg_data2)
summary(shan_tg_sex)
anova(shan_tg_sex)

piel_tg_sex <- lmer(pielou_evenness ~ sex + (1 | Location), data = tg_data2)
summary(piel_tg_sex)
anova(piel_tg_sex)



faith_1_sex <- lmer(faith_pd ~ sex + (1 | Location), data = tt_data2)
summary(faith_1_sex)
anova(faith_1_sex)

obs_1_sex <- lmer(observed_features ~ sex + (1 | Location), data = tt_data2)
summary(obs_1_sex)
anova(obs_1_sex)

shan_tt_sex <- lmer(shannon_entropy ~ sex + (1 | Location), data = tt_data2)
summary(shan_tt_sex)
anova(shan_tt_sex)
boxplot(tt_data2$shannon_entropy ~ tt_data2$sex)

piel_tt_sex <- lmer(pielou_evenness ~ sex + (1 | Location), data = tt_data2)
summary(piel_tt_sex)
anova(piel_tt_sex)
  
#mass lmers
tg_mass_lm <- lmer(observed_features ~ mass_g + (1 | Location), data = tg_data2)
summary(tg_mass_lm)
anova(tg_mass_lm)

tg_mass_shannon <- lmer(shannon_entropy ~ mass_g + (1 | Location), data = tg_data2)
summary(tg_mass_shannon)
anova(tg_mass_shannon)

tg_fatih_mass <- lmer(faith_pd ~ mass_g + (1 | Location), data = tg_data2)
summary(tg_fatih_mass)
anova(tg_fatih_mass)

tg_pielou_mass <- lmer(pielou_evenness ~ mass_g + (1 | Location), data = tg_data2)
anova(tg_pielou_mass)
summary(tg_pielou_mass)



tt_obs_mass <- lmer(observed_features ~ mass_g + (1 | Location), data = tt_data2)
summary(tt_obs_mass)
car::Anova(tt_obs_mass, type = 'III')

tt_shan_mass <- lmer(shannon_entropy ~ mass_g + (1 | Location), data = tt_data2)
summary(tt_shan_mass)
car::Anova(tt_shan_mass, type = 'III')
tt_faith_mass <- lmer(faith_pd ~ mass_g + (1 | Location), data = tt_data2)
summary(tt_faith_mass)
car::Anova(tt_faith_mass, type = 'III')
tt_pielou_mass <- lmer(pielou_evenness ~ mass_g + (1 | Location), data = tt_data2)
car::Anova(tt_pielou_mass, type = 'III')
summary(tt_pielou_mass)

#Plots for these
require(plyr)
fit <- faith_tg
anova(faith_tg)
sjPlot::tab_model(faith_tg)
parameters::model_parameters(faith_tg)
plot_model(fit, type = 'pred', terms = c('High', 'Low'))
faith_beta <- data.frame(summary(faith_tg)$coef)
betas <- mutate(faith_beta, name = row.names(faith_beta), 
                beta = Estimate, lo  = beta - 2*Std..Error, hi = beta + 2*Std..Error)
betas <- subset(betas, grepl("Chick", name))
ggplot(betas, aes(y=beta, x=name, ymin=lo, ymax=hi)) + geom_pointrange() + coord_flip()
tg_alpha <- fcsp_neg %>%
  ggplot(aes(x=TTX_mg,y=faith_pd))+
  geom_point(aes(color=Location),
             size = 3)+
  theme_minimal()+
  geom_smooth(method = lm,colour="black")+
  scale_colour_manual(values = tg_colors,
                      labels = c("FCSP", "GRP", "LSSP", "SLNF", "TMC", "WBC"))+
  xlab('TTX (mg)')+
  ylab("Faith's Phylogenetic Diversity")+
  theme(axis.line = element_line(colour = 'black'),
        panel.border = element_blank(),
        axis.text = element_text(face = 'bold',
                                 size = 14),
        axis.title = element_text(face = 'bold',
                                  size = 14),
        legend.text = element_text(face = 'bold',
                                   size = 14),
        legend.title = element_text(face = 'bold',
                                    size = 20),
        panel.grid = element_blank())+
  coord_cartesian(ylim = c(0,NA))+
  labs(tag='a)')
tg_alpha

#Beta_analysis
setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/')

#Weighted_unifrac__________________________________________________
weighted_unifram <- qiime2R::read_qza('core_metrics/all_samples/beta/weighted_unifrac/weighted_unifrac_distance_matrix.qza')
unifrac_matrix <- as(weighted_unifram$data, "matrix")
unifrac_matrix <- as.data.frame(unifrac_matrix)
met_2200 <- read.csv('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/2200_full_metadata.csv', header = T)
row_to_remove <- which(rownames(unifrac_matrix) == 'FCSP-05')
col_to_rmv <- which(colnames(unifrac_matrix) == 'FCSP-05')
unifrac_matrix <- unifrac_matrix[-row_to_remove,-col_to_rmv]
ids <- rownames(unifrac_matrix)
filtered_metadata <- met_2200[met_2200$SampleID %in% ids,]
rownames(filtered_metadata) <- filtered_metadata$SampleID
new_order <- order(rownames(unifrac_matrix))

new_matrix <- unifrac_matrix[new_order,]
new_matrix1 <- new_matrix[,new_order]
new_matrix1 <- as.matrix(new_matrix1)
uni_matrix <- as.matrix(unifrac_matrix)
tg_weighted_unifrac_matrix <- new_matrix1

dist1 <- dist(tg_weighted_unifrac_matrix)
dist2 <- dist(filtered_metadata$TTX_mg)
mantel_result <- mantel(dist1, dist2, method = 'spearman')
print(mantel_result)
#jaccard matrix__________________________________________________________________________________________--
jac_import <- qiime2R::read_qza('core_metrics/all_samples/beta/jaccard/jaccard_distance_matrix.qza')
jac_mat1 <- as(jac_import$data, 'matrix')
jac_mat1 <- as.data.frame(jac_mat1)
jac_row_to_remove <- which(rownames(jac_mat1) == 'FCSP-05')
jac_col_to_rmv <- which(colnames(jac_mat1) == 'FCSP-05')
jac_mat2 <- jac_mat1[-jac_row_to_remove, -jac_col_to_rmv]
jac_ids <- rownames(jac_mat2)
jac_new_order <- order(rownames(jac_mat2))
jac_ordered <- jac_mat2[jac_new_order,jac_new_order]
tg_jaccard_matrix <- as.matrix(jac_ordered)

dist3 <- dist(tg_jaccard_matrix)
mantel_result2 <- mantel(dist3, dist2, method = 'spearman')
print(mantel_result2)
#Unweighted-Unifrac___________________________________________________________________________________________--
uwu_import <- qiime2R::read_qza('core_metrics/all_samples/beta/unweighted/unweighted_unifrac_distance_matrix.qza')
uwu_mat1 <- as(uwu_import$data, 'matrix')
uwu_mat1 <- as.data.frame(uwu_mat1)
uwu_row_to_remove <- which(rownames(uwu_mat1) == 'FCSP-05')
uwu_col_to_rmv <- which(colnames(uwu_mat1) == 'FCSP-05')
uwu_mat2 <- uwu_mat1[-uwu_row_to_remove, -uwu_col_to_rmv]
uwu_ids <- rownames(uwu_mat2)
uwu_new_order <- order(rownames(uwu_mat2))
uwu_ordered <- uwu_mat2[uwu_new_order,uwu_new_order]
tg_unweighted_unifrac_matrix <- as.matrix(uwu_ordered)

dist4 <- dist(tg_unweighted_unifrac_matrix)
mantel_result3 <- mantel(dist4, dist2, method = 'spearman')
print(mantel_result3)
#_Bray-Curtis___________________________________________________________________________________________--
bc_import <- qiime2R::read_qza('core_metrics/all_samples/beta/bray_curtis/bray_curtis_distance_matrix.qza')
bc_mat1 <- as(bc_import$data, 'matrix')
bc_mat1 <- as.data.frame(bc_mat1)
bc_row_to_remove <- which(rownames(bc_mat1) == 'FCSP-05')
bc_col_to_rmv <- which(colnames(bc_mat1) == 'FCSP-05')
bc_mat2 <- bc_mat1[-bc_row_to_remove, -bc_col_to_rmv]
bc_ids <- rownames(bc_mat2)
bc_new_order <- order(rownames(bc_mat2))
bc_ordered <- bc_mat2[bc_new_order,bc_new_order]
tg_bray_curtis_matrix <- as.matrix(bc_ordered)

dist5 <- dist(tg_bray_curtis_matrix)
mantel_result4 <- mantel(dist5, dist2, method = 'spearman')
print(mantel_result4)
#_______________________________________________________________________________________________________

permanova_1 <- adonis2(formula = as.dist(new_matrix1) ~ filtered_metadata$TTX_mg * filtered_metadata$Location, na.rm = T)
print(permanova_1)
nm_plot_mat1 <- metaMDS(new_matrix1, k = 3)


perm_1.5 <- adonis2(formula = new_matrix1 ~ filtered_metadata$TTX_mg * filtered_metadata$Location)
print(perm_1.5)
permanova_2 <- adonis2(uni_matrix ~ filtered_metadata$shannon_entropy, na.rm = T)
print(permanova_2)

perm_3 <- adonis2(formula = as.dist(new_matrix1) ~ filtered_metadata$Location)
summary(perm_3)
perm_3
pairwise_perm_test <- permanova_pairwise(as.dist(new_matrix1), filtered_metadata$Location)
pairwise_perm_test1 <- permanova_pairwise(as.dist(new_matrix1), filtered_metadata$Location)
write.csv(pairwise_perm_test1, 'pairwise_perm_test_unifraccsv', row.names = F)
pcoa_res <- cmdscale(uni_matrix)
df_pc <- as.data.frame(pcoa_res)
plot(df_pc)

adonis2(tg_bray_curtis_matrix ~ filtered_metadata$TTX_mg, data = filtered_metadata)
adonis2(tg_weighted_unifrac_matrix ~ filtered_metadata$TTX_mg, data = filtered_metadata)
adonis2(tg_unweighted_unifrac_matrix ~ filtered_metadata$TTX_mg, data = filtered_metadata)
adonis2(tg_unweighted_unifrac_matrix ~ filtered_metadata$mass_g, data = filtered_metadata)
jacadon <- adonis2(tg_jaccard_matrix ~ TTX_mg, data = filtered_metadata)
jacadon


#________________________________________________________________________________________
perm_results_df <- data.frame(matrix = character(),
                              predictor = character(),
                              response = character(),
                              pseudo_F = numeric(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)

matrix_list <- list(tg_weighted_unifrac_matrix = tg_weighted_unifrac_matrix, 
                    tg_unweighted_unifrac_matrix = tg_unweighted_unifrac_matrix,
                    tg_jaccard_matrix = tg_jaccard_matrix,
                    tg_bray_curtis_matrix = tg_bray_curtis_matrix)

predictor_list <- c('TTX_mg', 'Location', 
                    'Log_bd', 'mass_g')

# Iterate over distance matrices and predictor variables
for(matrix_name in names(matrix_list)) {
  for(predictor_var in predictor_list) {
    perm_result <- adonis2(as.matrix(matrix_list[[matrix_name]]) ~ filtered_metadata[[predictor_var]], data = filtered_metadata)
    perm_results_df <- rbind(perm_results_df, 
                             data.frame(matrix = matrix_name,
                                        predictor = predictor_var,
                                        pseudo_F = perm_result$F[1],
                                        p_value = perm_result$`Pr(>F)`[1]))
  }
}


# Write the results to a CSV file
write.csv(perm_results_df, 'tg_perm_results.csv', row.names = FALSE)

# Optional: Print the results in a table format
library(knitr)
table_txt <- kable(perm_results_df, format = "markdown", 
                   col.names = c("Distance Matrix", 
                                 "Predictor Variable", 
                                 "Response Variable", 
                                 "Pseudo-F", 
                                 "p-value"))
print(table_txt)

#TG mantel____________________________________________________________________________________

library(vegan)

matrix_list <- list(tg_weighted_unifrac_matrix = tg_weighted_unifrac_matrix, 
                    tg_unweighted_unifrac_matrix = tg_unweighted_unifrac_matrix,
                    tg_jaccard_matrix = tg_jaccard_matrix,
                    tg_bray_curtis_matrix = tg_bray_curtis_matrix)

variables <- list(as.matrix(filtered_metadata$TTX_mg), 
                  as.matrix(filtered_metadata$Log_bd), 
                  as.matrix(filtered_metadata$mass_g))

library(vegan)


mantel_results_df <- data.frame(
  matrix = character(),
  variable = character(),
  statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (predictor_name in names(matrix_list)) {
  predictor_matrix <- matrix_list[[predictor_name]]
  
  for (var_name in names(variables)) {
    var_matrix <- variables[[var_name]]
    
    num_columns <- ncol(predictor_matrix)
    
    for (i in 1:num_columns) {
      mantel_result <- mantel(predictor_matrix[, i], 
                              var_matrix,
                              permutations = 999,
                              method = "spearman")
      
      mantel_results_df <- rbind(mantel_results_df, data.frame(
        matrix = predictor_name,
        matrix_column = paste0("column_", i),
        variable = var_name,
        statistic = mantel_result$statistic,
        p_value = mantel_result$signif
      ))
    }
  }
}

print(mantel_results_df)
write.csv(mantel_results_df, file = "mantel_results.csv", row.names = FALSE)


#___________________________________________________________________
tt_perm_results_df <- data.frame(matrix = character(),
                                 predictor = character(),
                                 response = character(),
                                 pseudo_F = numeric(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

tt_matrix_list <- list(tt_unweighted_unifrac_matrix = tt_unweighted_unifrac_matrix,
                       tt_weighted_unifrac_matrix = tt_weighted_unifrac_matrix,
                       tt_jaccard_matrix = tt_jaccard_matrix,
                       tt_bray_curtis_matrix = tt_bray_curtis_matrix)

tt_predictor_list <- c('TTX_mg', 'Location', 'Log_bd')

# Iterate over distance matrices and predictor variables
for(matrix_name in names(tt_matrix_list)) {
  for(predictor_var in tt_predictor_list) {
    tt_perm_result <- adonis2(as.matrix(tt_matrix_list[[matrix_name]]) ~ tt_filtered_metadata[[predictor_var]], data = tt_filtered_metadata)
    tt_perm_results_df <- rbind(tt_perm_results_df, 
                                data.frame(matrix = matrix_name,
                                           predictor = predictor_var,
                                           pseudo_F = tt_perm_result$F[1],
                                           p_value = tt_perm_result$`Pr(>F)`[1]))
  }
}


# Write the results to a CSV file
write.csv(tt_perm_results_df, 'tt_perm_results.csv', row.names = FALSE)

# Optional: Print the results in a table format
library(knitr)
table_txt <- kable(perm_results_df, format = "markdown", 
                   col.names = c("Distance Matrix", 
                                 "Predictor Variable", 
                                 "Response Variable", 
                                 "Pseudo-F", 
                                 "p-value"))
print(table_txt)
#_________________________________________________________________________________________-

nmds_t <- metaMDS(uni_matrix, k = 2)
cords <- scores(nmds_t)
ordiplot(cords, type = 'n')
locs <- filtered_metadata$Location
colorss <- rainbow(length(unique(locs)))
points(cords, col = colorss[as.numeric(factor(locs))], pch = 16)

#_____________________________________________________________________________
location_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                     'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')
tg_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan')
tt_colors <- c('salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')


tg_prcomp3 <- wcmdscale(new_matrix1, eig = TRUE)
print(tg_prcomp3)
require(ggplot2)
tg_part_1 <- ggpcoa(ord = tg_prcomp3, 
                    ordata = new_matrix1,
                    spearrow = NULL, 
                    spe = FALSE, 
                    cirline = 3,
                    ellprob = 0.95,
                    ellipse = TRUE,
                    groups = filtered_metadata$Location
)+ 
  ggtitle('Taricha granulosa')+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = factor(filtered_metadata$Location)),
             size = 2)+
  guides(color = guide_legend('Location'),
         size = FALSE,
         shape = FALSE)+
  scale_shape_manual(values = c(16, 16, 16,
                                16, 16, 16))+
  scale_fill_manual(values = tg_colors)+
  scale_color_manual(values = tg_colors)+
  labs(tag = 'a)')

tg_part_1


#_________________________________________________________________________

tg_prcomp2 <- wcmdscale(new_matrix1, eig = TRUE)
require(ggplot2)
tg_part_2 <- ggpcoa(ord = tg_prcomp2, 
                    ordata = new_matrix1,
                    spearrow = NULL, 
                    spe = FALSE, 
                    cirline = 3,
                    ellprob = 0.95,
                    ellipse = TRUE
)+ 
  ggtitle('Taricha granulosa') +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = filtered_metadata$TTX_mg,
                 fill = filtered_metadata$TTX_mg,
                 shape = filtered_metadata$Location), size = 4) +
  scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                        values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                   4.5, 5, 6.1, 6.2))) +
  scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                       values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                  4.5, 5, 6.1, 6.2)))+
  guides(size = FALSE,
         color = guide_colorbar('TTX (mg)', order = 1),
         shape = guide_legend('Location', order = 2),
         fill = FALSE)+
  scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
  labs(tag = 'b)')


tg_part_2
#_________________________________________

tg_part_3 <- ggpcoa(ord = tg_prcomp2, 
                    ordata = new_matrix1,
                    spearrow = NULL, 
                    spe = FALSE, 
                    cirline = 3,
                    ellprob = 0.95,
                    ellipse = TRUE
)+ 
  ggtitle('Taricha granulosa') +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = as.factor(filtered_metadata$infected),
                 fill = as.factor(filtered_metadata$infected),
                 shape = filtered_metadata$Location), size = 4) +
  scale_color_manual(labels = c('Negative', 'Positive'),
                     values = c('#00BFC4',
                                "#F8766D"))+
  scale_fill_manual(values = c('#00BFC4',
                               "#F8766D"))+
  guides(size = FALSE,
         shape = guide_legend('Location', 
                              order = 2),
         color = guide_legend('Infection Status',
                              order = 1),
         fill = FALSE)+
  scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
  labs(tag = 'c)')


tg_part_3



nmds_plots <- ggarrange(tg_part_1, tg_part_2, tg_part_3, ncol = 3,
                        tt_part_1, part_2, part_3, nrow = 2)
nmds_plots
scaled <- arrangeGrob(nmds_plots, scale = 0.8)
scaled
ggsave('ndms_plots.pdf', 
       plot = nmds_plots, 
       device = 'pdf',
       width = 16)



#_________________________________________________________________________________

tt_wu_import <- qiime2R::read_qza('tato_new/core_metrics/weighted_unifrac/weighted_unifrac_distance_matrix.qza')
tt_wu_mat1 <- as(tt_wu_import$data, 'matrix')
tt_wu_mat1 <- as.data.frame(tt_wu_mat1)
tt_wu_ids <- rownames(tt_wu_mat1)
tt_wu_new_order <- order(rownames(tt_wu_mat1))
tt_wu_ordered <- tt_wu_mat1[tt_wu_new_order,tt_wu_new_order]
tt_weighted_unifrac_matrix <- as.matrix(tt_wu_ordered)

tt_met_2200 <- read.table('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/new_tt/new_metadata.tsv', header = T, sep = '\t')
tt_filtered_metadata <- tt_met_2200 #[tt_met_2200$SampleID %in% tt_wu_ids,]
rownames(tt_filtered_metadata) <- tt_filtered_metadata$SampleID
tt_order <- order(rownames(tt_filtered_metadata))
tt_filtered_metadata <- tt_filtered_metadata[tt_order,]
tt_new_matrix1 <- as.matrix(tt_weighted_unifrac_matrix)

tt_dist1 <- dist(tt_new_matrix1)
tt_met <- dist(tt_filtered_metadata$TTX_mg)
tt_mantel_result1 <- mantel(tt_dist1, tt_met, method = 'spearman')
print(tt_mantel_result1)
#___________________________________________________________
tt_jac_import <- qiime2R::read_qza('tato_new/core_metrics/jaccard/jaccard_distance_matrix.qza')
tt_jac_mat1 <- as(tt_jac_import$data, 'matrix')
tt_jac_mat1 <- as.data.frame(tt_jac_mat1)
tt_jac_ids <- rownames(tt_jac_mat1)
tt_jac_new_order <- order(rownames(tt_jac_mat1))
tt_jac_ordered <- tt_jac_mat1[tt_jac_new_order,tt_jac_new_order]
tt_jaccard_matrix <- as.matrix(tt_jac_ordered)

tt_dist2 <- dist(tt_jaccard_matrix)
tt_mantel_result2 <- mantel(tt_dist2, tt_met, method = 'spearman')
print(tt_mantel_result2)
#___________________________________________________________
tt_uwu_import <- qiime2R::read_qza('tato_new/core_metrics/unweighted/unweighted_unifrac_distance_matrix.qza')
tt_uwu_mat1 <- as(tt_uwu_import$data, 'matrix')
tt_uwu_mat1 <- as.data.frame(tt_uwu_mat1)
tt_uwu_ids <- rownames(tt_uwu_mat1)
tt_uwu_new_order <- order(rownames(tt_uwu_mat1))
tt_uwu_ordered <- tt_uwu_mat1[tt_uwu_new_order,tt_uwu_new_order]
tt_unweighted_unifrac_matrix <- as.matrix(tt_uwu_ordered)

tt_dist3 <- dist(tt_unweighted_unifrac_matrix)
tt_mantel_result3 <- mantel(tt_dist3, tt_met, method = 'spearman')
print(tt_mantel_result3)
#_Bray-Curtis___________________________________________________________________________________________--
tt_bc_import <- qiime2R::read_qza('tato_new/core_metrics/bray_curtis/bray_curtis_distance_matrix.qza')
tt_bc_mat1 <- as(tt_bc_import$data, 'matrix')
tt_bc_mat1 <- as.data.frame(tt_bc_mat1)
tt_bc_ids <- rownames(tt_bc_mat1)
tt_bc_new_order <- order(rownames(tt_bc_mat1))
tt_bc_ordered <- tt_bc_mat1[tt_bc_new_order,tt_bc_new_order]
tt_bray_curtis_matrix <- as.matrix(tt_bc_ordered)

tt_dist4 <- dist(tt_bray_curtis_matrix)
tt_mantel_result4 <- mantel(tt_dist4, tt_met, method = 'spearman')
print(tt_mantel_result4)

adonis2(tt_bray_curtis_matrix ~ Log_bd, data = tt_filtered_metadata)
adonis2(tt_weighted_unifrac_matrix ~ Log_bd, data = tt_filtered_metadata)

#___________________________________________________________________
tt_perm_results_df <- data.frame(matrix = character(),
                                 predictor = character(),
                                 response = character(),
                                 pseudo_F = numeric(),
                                 p_value = numeric(),
                                 stringsAsFactors = FALSE)

tt_matrix_list <- list(tt_unweighted_unifrac_matrix = tt_unweighted_unifrac_matrix,
                       tt_weighted_unifrac_matrix = tt_weighted_unifrac_matrix,
                       tt_jaccard_matrix = tt_jaccard_matrix,
                       tt_bray_curtis_matrix = tt_bray_curtis_matrix)

tt_predictor_list <- c('TTX_mg', 'Log_bd')

# Iterate over distance matrices and predictor variables
for(matrix_name in names(tt_matrix_list)) {
  for(predictor_var in tt_predictor_list) {
    tt_perm_result <- adonis2(as.matrix(tt_matrix_list[[matrix_name]]) ~ tt_filtered_metadata[[predictor_var]], data = tt_filtered_metadata)
    tt_perm_results_df <- rbind(tt_perm_results_df, 
                                data.frame(matrix = matrix_name,
                                           predictor = predictor_var,
                                           pseudo_F = tt_perm_result$F[1],
                                           p_value = tt_perm_result$`Pr(>F)`[1]))
  }
}


# Write the results to a CSV file
write.csv(tt_perm_results_df, 'tt_adonis_results.csv', row.names = FALSE)

# Optional: Print the results in a table format


#___________________________________________________________________


tt_prcomp2 <- wcmdscale(tt_new_matrix1, eig = TRUE)
part_1 <- ggpcoa(ord = tt_prcomp2, 
                 ordata = tt_new_matrix1,
                 spearrow = NULL, 
                 spe = FALSE, 
                 cirline = 3,
                 ellprob = 0.95,
                 ellipse = TRUE,
                 groups = tt_filtered_metadata$Location
)+ 
  ggtitle('Taricha torosa')+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic'))+
  geom_point(aes(color = factor(tt_filtered_metadata$Location)),
             size = 2)+
  guides(color = guide_legend('Location'),
         size = FALSE,
         shape = FALSE)+
  scale_shape_manual(values = c(16, 16, 16,
                                16, 16, 16))+
  scale_fill_manual(values = tt_colors)+
  scale_color_manual(values = tt_colors,
                     labels = c(
                       'Bolinger', 'Cold Creek', 'Crocker',
                       'Gunstock', 'Madonna', 'Muir'
                     ))+
  labs(tag = 'd)')

part_1
#____________________________________________________________

part_2 <- ggpcoa(ord = tt_prcomp2, 
                 ordata = tt_new_matrix1,
                 spearrow = NULL, 
                 spe = FALSE, 
                 cirline = 3,
                 ellprob = 0.95,
                 ellipse = TRUE
)+ 
  ggtitle('Taricha torosa') +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = tt_filtered_metadata$TTX_mg,
                 fill = tt_filtered_metadata$TTX_mg,
                 shape = tt_filtered_metadata$Location), size = 4) +
  scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                        values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                   4.5, 5, 6.1, 6.2))) +
  scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                       values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                  4.5, 5, 6.1, 6.2)))+
  guides(size = FALSE,
         shape = guide_legend('Location'),
         color = guide_colorbar('TTX (mg)'),
         fill = FALSE)+
  scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
  labs(tag = 'e)')


part_2




#____________________________________________________________

part_3 <- ggpcoa(ord = tt_prcomp2, 
                 ordata = tt_new_matrix1,
                 spearrow = NULL, 
                 spe = FALSE, 
                 cirline = 3,
                 ellprob = 0.95,
                 ellipse = TRUE
)+ 
  ggtitle('Taricha torosa') +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = as.factor(tt_filtered_metadata$infected),
                 fill = as.factor(tt_filtered_metadata$infected),
                 shape = tt_filtered_metadata$Location), size = 4) +
  scale_color_manual(labels = c('Negative', 'Positive'),
                     values = c('#00BFC4',
                                "#F8766D"))+
  scale_fill_manual(values = c('#00BFC4',
                               "#F8766D"))+
  guides(size = FALSE,
         shape = guide_legend('Location', order = 2),
         color = guide_legend('Infection Status', order = 1),
         fill = FALSE)+
  scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
  labs(tag = 'f)')


part_3

ggarrange(part_1, part_2 ,part_3)




nmds_plots <- ggarrange(tg_part_1, tg_part_2, tg_part_3, ncol = 3,
                        part_1, part_2, part_3, nrow = 2)

nmds_plots    









#_____________Merged_PCoA_Files_____________________________________________________________________________________-
#___________________________________________
merged_weighted_unifram <- qiime2R::read_qza('tt_tg_merged/merged_core_metrics/weighted_unifrac_distance_matrix.qza')
merged_unifrac_matrix <- as(merged_weighted_unifram$data, "matrix")
merged_unifrac_matrix <- as.data.frame(merged_unifrac_matrix)
met_2200 <- read.csv('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/2200_full_metadata.csv', header = T)
row_to_remove <- which(rownames(merged_unifrac_matrix) == 'FCSP-05')
col_to_rmv <- which(colnames(merged_unifrac_matrix) == 'FCSP-05')
merged_unifrac_matrix <- merged_unifrac_matrix[-row_to_remove,-col_to_rmv]
merged_ids <- rownames(merged_unifrac_matrix)
merged_filtered_metadata <- met_2200[met_2200$SampleID %in% merged_ids,]
rownames(merged_filtered_metadata) <- merged_filtered_metadata$SampleID
new_order <- order(rownames(merged_unifrac_matrix))

merged_new_matrix <- merged_unifrac_matrix[new_order,]
merged_new_matrix1 <- merged_new_matrix[,new_order]
merged_new_matrix1 <- as.matrix(merged_new_matrix1)
merged_weighted_unifrac_matrix <- merged_new_matrix1

ttx_perm_merged <- adonis2(formula = as.dist(merged_weighted_unifrac_matrix) ~ merged_filtered_metadata$TTX_mg, na.action = na.omit)
print(ttx_perm_merged)

merged_prcomp <- wcmdscale(merged_weighted_unifrac_matrix, eig = TRUE)
merged_pcoa <- ggpcoa(ord = merged_prcomp, 
                      ordata = merged_weighted_unifrac_matrix,
                      spearrow = NULL, 
                      spe = FALSE, 
                      cirline = 1,
                      ellprob = 0.95,
                      ellipse = TRUE,
                      groups = merged_filtered_metadata$Species
)+
  geom_point(aes(color = factor(merged_filtered_metadata$Species)),
             size = 4)+
  guides(color = guide_legend('Species',
                              title.theme = element_text(size = 20),
                              label.theme = element_text(face = 'italic',
                                                         size = 15)),
         size = FALSE,
         shape = FALSE)+
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16,
                                16, 16, 16, 16, 16, 16))+
  scale_fill_manual(values = c("#F8766D", '#00BFC4'))+
  scale_color_manual(values = c("#F8766D", '#00BFC4'),
                     labels = c('Taricha granulosa', 'Taricha torosa')) +
  labs(tag = 'a)')

merged_pcoa
#merged file for manuscript. other graph has to be loaded in and is named 'taxa_barplots'

merged_ggrange <- ggarrange(merged_pcoa, combined_taxonomy_barplot)
merged_ggrange





#___________________________________________________________
merged_jac_import <- qiime2R::read_qza('tt_tg_merged/merged_core_metrics/jaccard_distance_matrix.qza')
merged_jac_mat1 <- as(merged_jac_import$data, 'matrix')
merged_jac_mat1 <- as.data.frame(merged_jac_mat1)
row_to_remove <- which(rownames(merged_jac_mat1) == 'FCSP-05')
col_to_rmv <- which(colnames(merged_jac_mat1) == 'FCSP-05')
merged_jac_mat1 <- merged_jac_mat1[-row_to_remove,-col_to_rmv]
merged_ids <- rownames(merged_jac_mat1)
merged_jac_ids <- rownames(merged_jac_mat1)
merged_jac_new_order <- order(rownames(merged_jac_mat1))
merged_jac_ordered <- merged_jac_mat1[merged_jac_new_order,merged_jac_new_order]
merged_jaccard_matrix <- as.matrix(merged_jac_ordered)
#___________________________________________________________
merged_uwu_import <- qiime2R::read_qza('tt_tg_merged/merged_core_metrics/unweighted_unifrac_distance_matrix.qza')
merged_uwu_mat1 <- as(merged_uwu_import$data, 'matrix')
merged_uwu_mat1 <- as.data.frame(merged_uwu_mat1)
row_to_remove <- which(rownames(merged_uwu_mat1) == 'FCSP-05')
col_to_rmv <- which(colnames(merged_uwu_mat1) == 'FCSP-05')
merged_uwu_mat1 <- merged_uwu_mat1[-row_to_remove,-col_to_rmv]
merged_ids <- rownames(merged_uwu_mat1)
merged_uwu_ids <- rownames(merged_uwu_mat1)
merged_uwu_new_order <- order(rownames(merged_uwu_mat1))
merged_uwu_ordered <- merged_uwu_mat1[merged_uwu_new_order,merged_uwu_new_order]
merged_unweighted_unifrac_matrix <- as.matrix(merged_uwu_ordered)
#_Bray-Curtis___________________________________________________________________________________________--
merged_bc_import <- qiime2R::read_qza('tt_tg_merged/merged_core_metrics/bray_curtis_distance_matrix.qza')
merged_bc_mat1 <- as(merged_bc_import$data, 'matrix')
merged_bc_mat1 <- as.data.frame(merged_bc_mat1)
row_to_remove <- which(rownames(merged_bc_mat1) == 'FCSP-05')
col_to_rmv <- which(colnames(merged_bc_mat1) == 'FCSP-05')
merged_bc_mat1 <- merged_bc_mat1[-row_to_remove,-col_to_rmv]
merged_ids <- rownames(merged_bc_mat1)

merged_bc_ids <- rownames(merged_bc_mat1)
merged_bc_new_order <- order(rownames(merged_bc_mat1))
merged_bc_ordered <- merged_bc_mat1[merged_bc_new_order,merged_bc_new_order]
merged_bray_curtis_matrix <- as.matrix(merged_bc_ordered)
#_______________________________________________________________________________________________________

merged_perm_results_df <- data.frame(matrix = character(),
                                     predictor = character(),
                                     response = character(),
                                     pseudo_F = numeric(),
                                     p_value = numeric(),
                                     stringsAsFactors = FALSE)

merged_matrix_list <- list(merged_unweighted_unifrac_matrix = merged_unweighted_unifrac_matrix,
                           merged_weighted_unifrac_matrix = merged_weighted_unifrac_matrix,
                           merged_jaccard_matrix = merged_jaccard_matrix,
                           merged_bray_curtis_matrix = merged_bray_curtis_matrix)

merged_predictor_list <- c('TTX_mg', 'Location', 'Log_bd', 'Species')
merged_filtered_metadata <- met_2200[met_2200$SampleID %in% merged_ids,]


# Iterate over distance matrices and predictor variables
for(matrix_name in names(merged_matrix_list)) {
  for(predictor_var in merged_predictor_list) {
    merged_perm_result <- adonis2(as.matrix(merged_matrix_list[[matrix_name]]) ~ merged_filtered_metadata[[predictor_var]], data = merged_filtered_metadata)
    merged_perm_results_df <- rbind(merged_perm_results_df, 
                                    data.frame(matrix = matrix_name,
                                               predictor = predictor_var,
                                               pseudo_F = merged_perm_result$F[1],
                                               p_value = merged_perm_result$`Pr(>F)`[1]))
  }
}


# Write the results to a CSV file
write.csv(merged_perm_results_df, 'merged_perm_results.csv', row.names = FALSE)

# Optional: Print the results in a table format
library(knitr)
table_txt <- kable(merged_perm_results_df, format = "markdown", 
                   col.names = c("Distance Matrix", 
                                 "Predictor Variable", 
                                 "Response Variable", 
                                 "Pseudo-F", 
                                 "p-value"))
print(table_txt)


#__species pcoa________________________________________________________________________________
merged_prcomp <- wcmdscale(merged_weighted_unifrac_matrix, eig = TRUE)
merged_pcoa <- ggpcoa(ord = merged_prcomp, 
                      ordata = merged_weighted_unifrac_matrix,
                      spearrow = NULL, 
                      spe = FALSE, 
                      cirline = 1,
                      ellprob = 0.95,
                      ellipse = TRUE,
                      groups = merged_filtered_metadata$Species
)+ 
  ggtitle('Genus Taricha')+
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 20, 
                                  face = 'italic')) +
  geom_point(aes(color = factor(merged_filtered_metadata$Species)),
             size = 4)+
  guides(color = guide_legend('Location'),
         size = FALSE,
         shape = FALSE)+
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16,
                                16, 16, 16, 16, 16, 16))+
  scale_fill_manual(values = c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                               'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033'))+
  scale_color_manual(values = c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                                'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033'))+
  labs(tag = 'd)')

merged_pcoa

#_____________________________________________________________________________________________________-





#color_pallete
colorblind_2 <- c('red', 'blue', 'green', 'orange', 'black', 'cyan',
                  'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')
pie(rep(1,12), col = colorblind_2)
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", '#be0032', '#FFFF99', "#008837",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7", '#fa8072', '#9FE2BF')
pie(rep(1, 13), col = colorBlindBlack8)

values <- structure(c("#A6611A", "#D01C8B", "#7B3294", "#E66101", "#CA0020", 
                      "#CA0020", "#D7191C", "#D7191C", "#D7191C", "#7FC97F", "#1B9E77", 
                      "#A6CEE3", "#FBB4AE", "#B3E2CD", "#E41A1C", "#66C2A5", "#8DD3C7", 
                      "#EFF3FF", "#EDF8FB", "#EDF8FB", "#F0F9E8", "#EDF8E9", "#F7F7F7", 
                      "#FEEDDE", "#FEF0D9", "#F1EEF6", "#F6EFF7", "#F1EEF6", "#F2F0F7", 
                      "#FEEBE2", "#FEE5D9", "#FFFFCC", "#FFFFCC", "#FFFFD4", "#FFFFB2", 
                      "#DFC27D", "#F1B6DA", "#C2A5CF", "#FDB863", "#F4A582", "#F4A582", 
                      "#FDAE61", "#FDAE61", "#FDAE61", "#BEAED4", "#D95F02", "#1F78B4", 
                      "#B3CDE3", "#FDCDAC", "#377EB8", "#FC8D62", "#FFFFB3", "#BDD7E7", 
                      "#B2E2E2", "#B3CDE3", "#BAE4BC", "#BAE4B3", "#CCCCCC", "#FDBE85", 
                      "#FDCC8A", "#BDC9E1", "#BDC9E1", "#D7B5D8", "#CBC9E2", "#FBB4B9", 
                      "#FCAE91", "#C2E699", "#A1DAB4", "#FED98E", "#FECC5C", "#80CDC1", 
                      "#B8E186", "#A6DBA0", "#B2ABD2", "#92C5DE", "#BABABA", "#ABD9E9", 
                      "#A6D96A", "#ABDDA4", "#FDC086", "#7570B3", "#B2DF8A", "#CCEBC5", 
                      "#CBD5E8", "#4DAF4A", "#8DA0CB", "#BEBADA", "#6BAED6", "#66C2A4", 
                      "#8C96C6", "#7BCCC4", "#74C476", "#969696", "#FD8D3C", "#FC8D59", 
                      "#74A9CF", "#67A9CF", "#DF65B0", "#9E9AC8", "#F768A1", "#FB6A4A", 
                      "#78C679", "#41B6C4", "#FE9929", "#FD8D3C", "#018571", "#4DAC26", 
                      "#008837", "#5E3C99", "#0571B0", "#404040", "#2C7BB6", "#1A9641", 
                      "#2B83BA", "#FFFF99", "#E7298A", "#33A02C", "#DECBE4", "#F4CAE4", 
                      "#984EA3", "#E78AC3", "#FB8072", "#2171B5", "#238B45", "#88419D", 
                      "#2B8CBE", "#238B45", "#525252", "#D94701", "#D7301F", "#0570B0", 
                      "#02818A", "#CE1256", "#6A51A3", "#AE017E", "#CB181D", "#238443", 
                      "#225EA8", "#CC4C02", "#E31A1C")
                    
                    
                    #___tg________-bd_pos_sites
                    
                    #_____________-bray_curtis___________________
                    tg_prcomp4 <- wcmdscale(tg_bray_curtis_matrix, eig = TRUE)
                    print(tg_prcomp4)
                    require(ggplot2)
                    tg_part_12 <- ggpcoa(ord = tg_prcomp4, 
                                         ordata = tg_bray_curtis_matrix,
                                         spearrow = NULL, 
                                         spe = FALSE, 
                                         cirline = 3,
                                         ellprob = 0.95,
                                         ellipse = TRUE,
                                         groups = filtered_metadata$Location
                    )+ 
                      ggtitle('Taricha granulosa')+
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic')) +
                      geom_point(aes(color = factor(filtered_metadata$Location)),
                                 size = 2)+
                      guides(color = guide_legend('Location'),
                             size = FALSE,
                             shape = FALSE)+
                      scale_shape_manual(values = c(16, 16, 16,
                                                    16, 16, 16))+
                      scale_fill_manual(values = tg_colors)+
                      scale_color_manual(values = tg_colors)+
                      labs(tag = 'a)')
                    
                    tg_part_12
                    
                    
                    #_________________________________________________________________________
                    
                    tg_prcomp22 <- wcmdscale(tg_bray_curtis_matrix, eig = TRUE)
                    require(ggplot2)
                    tg_part_22 <- ggpcoa(ord = tg_prcomp22, 
                                         ordata = tg_bray_curtis_matrix,
                                         spearrow = NULL, 
                                         spe = FALSE, 
                                         cirline = 3,
                                         ellprob = 0.95,
                                         ellipse = TRUE
                    )+ 
                      ggtitle('Taricha granulosa') +
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic')) +
                      geom_point(aes(color = filtered_metadata$TTX_mg,
                                     fill = filtered_metadata$TTX_mg,
                                     shape = filtered_metadata$Location), size = 4) +
                      scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                            values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                       4.5, 5, 6.1, 6.2))) +
                      scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                           values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                      4.5, 5, 6.1, 6.2)))+
                      guides(size = FALSE,
                             color = guide_colorbar('TTX (mg)', order = 1),
                             shape = guide_legend('Location', order = 2),
                             fill = FALSE)+
                      scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
                      labs(tag = 'b)')
                    
                    
                    tg_part_22
                    #_________________________________________
                    
                    tg_part_32 <- ggpcoa(ord = tg_prcomp22, 
                                         ordata = tg_bray_curtis_matrix,
                                         spearrow = NULL, 
                                         spe = FALSE, 
                                         cirline = 3,
                                         ellprob = 0.95,
                                         ellipse = TRUE
                    )+ 
                      ggtitle('Taricha granulosa') +
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic')) +
                      geom_point(aes(color = filtered_metadata$Log_bd,
                                     fill = filtered_metadata$Log_bd,
                                     shape = filtered_metadata$Location), size = 4) +
                      scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                            values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                       4.5, 5, 6.1, 6.2))) +
                      scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                           values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                      4.5, 5, 6.1, 6.2)))+
                      guides(size = FALSE,
                             shape = guide_legend('Location', 
                                                  order = 2),
                             color = guide_colorbar('Infection Intensity',
                                                    order = 1),
                             fill = FALSE,
                             size = FALSE)+
                      scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
                      labs(tag = 'c)')
                    
                    
                    tg_part_32
                    
                    tg_bc <- ggarrange(tg_part_12, tg_part_22, tg_part_32, nrow = 1)
                    tg_bc
                    
                    #____________________________-torosa
                    tt_prcomp22 <- wcmdscale(tt_bray_curtis_matrix, eig = TRUE)
                    part_12 <- ggpcoa(ord = tt_prcomp22, 
                                      ordata = tt_bray_curtis_matrix,
                                      spearrow = NULL, 
                                      spe = FALSE, 
                                      cirline = 3,
                                      ellprob = 0.95,
                                      ellipse = TRUE,
                                      groups = tt_filtered_metadata$Location
                    )+ 
                      ggtitle('Taricha torosa')+
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic'))+
                      geom_point(aes(color = factor(tt_filtered_metadata$Location)),
                                 size = 2)+
                      guides(color = guide_legend('Location'),
                             size = FALSE,
                             shape = FALSE)+
                      scale_shape_manual(values = c(16, 16, 16,
                                                    16, 16, 16))+
                      scale_fill_manual(values = tt_colors)+
                      scale_color_manual(values = tt_colors,
                                         labels = c(
                                           'Bolinger', 'Cold Creek', 'Crocker',
                                           'Gunstock', 'Madonna', 'Muir'
                                         ))+
                      labs(tag = 'd)')
                    
                    part_12
                    #____________________________________________________________
                    
                    part_22 <- ggpcoa(ord = tt_prcomp22, 
                                      ordata = tt_bray_curtis_matrix,
                                      spearrow = NULL, 
                                      spe = FALSE, 
                                      cirline = 3,
                                      ellprob = 0.95,
                                      ellipse = TRUE
                    )+ 
                      ggtitle('Taricha torosa') +
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic')) +
                      geom_point(aes(color = tt_filtered_metadata$TTX_mg,
                                     fill = tt_filtered_metadata$TTX_mg,
                                     shape = tt_filtered_metadata$Location), size = 4) +
                      scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                            values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                       4.5, 5, 6.1, 6.2))) +
                      scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                           values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                      4.5, 5, 6.1, 6.2)))+
                      guides(size = FALSE,
                             shape = guide_legend('Location'),
                             color = guide_colorbar('TTX (mg)'),
                             fill = FALSE)+
                      scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
                      labs(tag = 'e)')
                    
                    
                    part_22
                    
                    
                    
                    
                    #____________________________________________________________
                    
                    part_32 <- ggpcoa(ord = tt_prcomp22, 
                                      ordata = tt_bray_curtis_matrix,
                                      spearrow = NULL, 
                                      spe = FALSE, 
                                      cirline = 3,
                                      ellprob = 0.95,
                                      ellipse = TRUE
                    )+ 
                      ggtitle('Taricha torosa') +
                      theme(plot.title = element_text(hjust = 0.5, 
                                                      size = 20, 
                                                      face = 'italic')) +
                      geom_point(aes(color = tt_filtered_metadata$Log_bd,
                                     fill = tt_filtered_metadata$Log_bd,
                                     shape = tt_filtered_metadata$Location), size = 4) +
                      scale_color_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                            values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                       4.5, 5, 6.1, 6.2))) +
                      scale_fill_gradientn(colours = c("#40E0D0", '#ff2626', 'black'),
                                           values = scales::rescale(c(0.25, 1, 1.5, 3, 4, 
                                                                      4.5, 5, 6.1, 6.2)))+
                      guides(size = FALSE,
                             shape = guide_legend('Location', 
                                                  order = 2),
                             color = guide_colorbar('Infection Intensity',
                                                    order = 1),
                             fill = FALSE)+
                      scale_shape_manual(values = c(16, 17, 18, 15, 20, 25))+
                      labs(tag = 'f)')
                    
                    
                    part_32
                    
                    tt_bc <- ggarrange(part_12, part_22 ,part_32, nrow = 1)
                    tt_bc
                    
                    final_bc_plots <- ggarrange(tg_bc, tt_bc,
                                                nrow = 2)
                    final_bc_plots
                    
                    
                    
                    #LDM
                    
                    #check to see if they're different
                    
                    BC_LDM1 <- permanovaFL(tg_bray_curtis_matrix ~ filtered_metadata$Location, 
                                           perm.between.type = 'free',
                                           perm.within.type = 'free',
                                           seed = 999,
                                           square.dist = FALSE,
                                           n.cores = 8,
                                           n.perm.max = 5000)
                    BC_LDM1$F.statistics
                    BC_LDM1$R.squared
                    BC_LDM1$p.permanova
                    
                    #how are they different
                    
                    BC_disp<-betadisper(as.dist(tg_bray_curtis_matrix), 
                                        filtered_metadata$Location,
                                        type = c("centroid"),
                                        bias.adjust = TRUE,
                                        sqrt.dist = FALSE,
                                        add = FALSE)
                    anova(BC_disp)
                    
                    
                    Taxa <-parse_taxonomy(read_qza("./taxonomy/classification.qza")$data)%>%
                      rownames_to_column("featureid")%>%
                      dplyr::select(-Kingdom)%>%
                      mutate(across(everything(),~gsub("[][]","",.)),
                             across(everything(),~replace_na(.,"")))
                    
                    for (i in 1:nrow(Taxa)){
                      if (Taxa[i,7] != ""){
                        Taxa$Species[i] <- paste(Taxa$Genus[i], Taxa$Species[i], sep=" ")
                      } else if (Taxa[i,2] == ""){
                        Phylum <- paste("Unclassified", "Bacteria", sep="-")
                        Taxa[i, 2:7] <- Phylum
                      } else if (Taxa[i,3] == ""){
                        Class <- paste(Taxa[i,2],"(P)", sep=" ")
                        Taxa[i, 3:7] <- Class
                      } else if (Taxa[i,4] == ""){
                        Order <- paste(Taxa[i,3],"(C)", sep=" ")
                        Taxa[i, 4:7] <- Order
                      } else if (Taxa[i,5] == ""){
                        Family <- paste(Taxa[i,4],"(O)", sep=" ")
                        Taxa[i, 5:7] <- Family
                      } else if (Taxa[i,6] == ""){
                        Genus <- paste(Taxa[i,5],"(F)", sep=" ")
                        Taxa[i, 6:7] <- Genus
                      } else if (Taxa[i,7] == ""){
                        Taxa$Species[i] <- paste(Taxa$Genus[i],"sp.", sep=" ")
                      }
                    }
                    
                    rar_tab <- as.data.frame(read_qza('./core_metrics/all_samples/2200/all_samples_table.qza')$data)
                    tg_2200_rare <- as.data.frame(t(rar_tab)) %>%
                      rownames_to_column("SampleID")
                    
                    tg_OTU_table <- rar_tab %>%
                      rownames_to_column('featureID')
                    
                    tg_2200_rare <- tg_2200_rare %>%
                      filter(!SampleID == 'FCSP-05') %>%
                      arrange("SampleID")
                    
                    matched_1 <- tg_2200_rare %>%
                      filter(SampleID %in% filtered_metadata$SampleID)%>%
                      column_to_rownames("SampleID")%>%
                      .[colSums(.[]) !=0] %>%
                      rownames_to_column("SampleID")
                    
                    matched_tab <- matched_1 %>%
                      .[ order(match(matched_1$SampleID, metadata_frame$SampleID)),] %>%
                      select(.,-SampleID)
                    
                    metadata_frame <- read.csv('~/Documents/Thesis_Part_2/Newt_TTX_Project/newt_ttx/newt_metadata.tsv', 
                                               sep = '\t',
                                               row.names = NULL) %>%
                      rename(SampleID = 1) %>%
                      filter(!SampleID == 'FCSP-05') %>%
                      as.data.frame() %>%
                      arrange(SampleID) %>%
                      filter(SampleID %in% filtered_metadata$SampleID)
                    
                    metadata_new <- metadata_frame %>%
                      select(SampleID, TTX_mg) %>%
                      column_to_rownames("SampleID")
                    
                    ttx_mg <- as.data.frame(metadata_new$TTX_mg)
                    colnames(ttx_mg)[1] <- "TTX_mg"
                    rownames_to_column(ttx_mg)
                    rownames(ttx_mg) <- NULL
                    any(duplicated(rownames(ttx_mg)))
                    any(duplicated(rownames(matched_tab)))
                    rownames(metadata_frame) <- seq_len(nrow(metadata_frame))
                    
                    tg_ttx_location_ldm <- ldm(matched_tab ~ metadata_frame$Site,
                                               perm.within.type="free", 
                                               perm.between.type="none",
                                               seed=999,dist.method = "bray", 
                                               fdr.nominal = 0.01, 
                                               n.cores = 10,
                                               n.perm.max = 10000)
