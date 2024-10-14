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
loc_bd_pos_tg <- tg_data2 %>%
  filter(!Location %in% 'FCSP')
loc_tt <- tt_data2 %>%
  filter(!Location %in% c('Gunstock', 'Madonna'))
tt_df3 <- tt_data2 %>%
  filter(Location %in% c('Cold_Creek', 'Muir'))
tt.data3 <- tt_data2 %>% 
  filter(!is.na(observed_features))
tg.data3 <- tg_data2 %>% 
  filter(!is.na(observed_features))
tg_data4 <- tg_data2 %>%
  subset(infected == 1)
tt_df4 <- tt_data2 %>%
  subset(infected == 1)

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
ggsave("figure_2.jpeg", ttx_location_arranged, dpi = 300,
        units = 'in', width = 15, height = 12)
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


#TAGR infection status #######
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

#Bd infection plots______________________________________________________________________________________
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
  expand_limits(y = c(0,5))+
  geom_text(data = dt3, aes(label = Letter, x = Location, y = log_bd + 0.1),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag='(A)')
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
  expand_limits(y = c(0,5))+
  geom_text(data = dt4, aes(label = Letter, x = Location, y = log_bd + 0.25),
            vjust = -0.5, hjust = 0.5, fontface = 'bold', size = 4, na.rm = TRUE) +
  labs(tag = '(B)')
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
  expand_limits(y = c(0,100))+
  labs(tag = '(C)')

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
  labs(tag = '(D)')
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
  expand_limits(y = c(0,5))+
  labs(tag = '(B)')
tt_bd_prev_plot

multi_plots <- ggarrange(bd_loc_plot_joined, prev_plots,
                         nrow = 2)
multi_plots

ggsave("figure_3.pdf", multi_plots, units = 'in', width = 11, height = 8, dpi = 500)
#infection as a function of TTX #########
tg_ttx_bd2 <- glmer(infected ~ TTX_mg + (1 | Location), data = fcsp_neg, family = 'binomial')
tg_ttx_bd_null <- glmer(infected ~ (1 | Location), data = fcsp_neg, family = 'binomial')
summary(tg_ttx_bd2)
drop1(tg_ttx_bd2, test = 'Chisq', na.action = 'na.omit')
lmtest::lrtest(tg_ttx_bd2, tg_ttx_bd_null)

tt_ttx_bd2 <- glmer(infected ~ TTX_mg + (1 | Location), data = tt_data2, family = 'binomial')
tt_ttx_bd_null <- glmer(infected ~ (1 | Location), data = tt_data2, family = 'binomial')
summary(tt_ttx_bd2)
drop1(tt_ttx_bd2, test = 'Chisq')
lmtest::lrtest(tt_ttx_bd2, tt_ttx_bd_null)
# neither are significant
################################################____________________________




#testing infection status as a factor of sex/location ####
tg_sex_glmm <- glmer(infected ~ sex + (1 | Location), data = tg_data2, family = 'binomial')
tg_sex_null <- glmer(infected ~ (1 | Location), data = tg_data2, family = 'binomial')
summary(tg_sex_glmm)
drop1(tg_sex_glmm, test = 'Chisq')
lmtest::lrtest(tg_sex_glmm, tg_sex_null)
#not significant


tt_sex_glmm <- glmer(infected ~ sex + (1 | Location), data = tt_data2, family = 'binomial')
tt_sex_glmm_null <- glmer(infected ~ (1 | Location), data = tt_data2, family = 'binomial')
summary(tt_sex_glmm)
drop1(tt_sex_glmm, test = 'Chisq')
lmtest::lrtest(tt_sex_glmm, tt_sex_glmm_null)
#not significant
################################################____________________________


  #Plot
mass_values <- seq(min(tg_data2$mass_g), max(tg_data2$mass_g), length.out = 100)  # Adjust as needed.
svl_values <- seq(min(tg_data2$svl_mm), max(tg_data2$svl_mm), length.out = 100)   # Adjust as needed.
unique_locations <- unique(tg_data2$Location)
ttx_weight <- seq(0.05, 7, 0.01)
location_colors <- data.frame(Location = unique_locations, Color = tg_colors[1:length(unique_locations)])
newdata <- merge(tg_data2, location_colors, by = "Location")
newdata <- expand.grid(TTX_mg = ttx_weight, Location = unique_locations)
# Expand the grid to include these values.
newdata <- expand.grid(
  TTX_mg = ttx_weight,
  Location = unique_locations,
  mass_g = mean(mass_values),   # Use mean or other appropriate value.
  svl_mm = mean(svl_values)     # Use mean or other appropriate value.
)

predicted_responses <- predict(tg_ttx_bd2, newdata = newdata, type = 'response')
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
    plot.title = element_blank(),
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
  labs(tag = '(A)')

tg_glmm

unique_locations_tt <- unique(tt_data2$Location)
tt_location_colors <- data.frame(Location = unique_locations_tt, Color = tt_colors[1:length(unique_locations_tt)])
tt_data2 <- merge(tt_data2, tt_location_colors, by = "Location")
tt_ttx_weight <- seq(0.002, 2, 0.01)
tt_newdata <- expand.grid(TTX_mg = tt_ttx_weight, Location = unique_locations_tt)
tt_predicted_responses <- predict(tt_ttx_bd2, newdata = tt_newdata, type = 'response')
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
    plot.title = element_blank(),
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
  labs(tag = '(B)')

new_glmm <- ggarrange(tg_glmm, tt_plot_glmm)
new_glmm

#tg_tt ttx-bd_prev plot ##########################################
tg.loc.prev <- fcsp_neg %>%
  group_by(Location) %>%
  dplyr::summarise(Proportion = mean(infected == 1) * 100,
                   mean_ttx = mean(TTX_mg),
                   ttx_se = sd(TTX_mg, na.rm = TRUE) / sqrt(n()),
                   sample_size = n(),
                   prev_sd = sqrt(sample_size * (Proportion / 100) * (1 - Proportion / 100))
                   )

tg_bd_prev_ttx_plot <- tg.loc.prev %>%  
  ggplot(aes(mean_ttx, Proportion, color = Location)) + 
  geom_errorbarh(aes(x = mean_ttx, xmin = mean_ttx - ttx_se, 
                     xmax = mean_ttx + ttx_se),
                size = 1,
                height = 0,
                color = 'black') +
  geom_errorbar(aes(y = Proportion,
                    ymin = Proportion - prev_sd,
                    ymax = Proportion + prev_sd),
                size = 1,
                width = 0,
                color = 'black') +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = 'Mean TTX (mg)',
       y = 'Infection Prevalence') +
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
    plot.title = element_blank(),
    axis.title.y = element_text(face = 'bold',
                                size = 14),
    axis.line = element_line(),
    legend.position = c(1,0.5),
    legend.justification = c('right', 'top'),
    legend.box.background = element_rect(),
    legend.title = element_text(face = "bold",
                                size = 20,
                                hjust = 0.5),
    legend.text = element_text(size = 15)
  ) + 
  labs(tag = '(C)') +
  scale_y_continuous(limits = c(0,100))
  

tg_bd_prev_ttx_plot


#tt plot
tt.loc.prev <- tt_data2 %>%
  group_by(Location) %>%
  dplyr::summarise(Proportion = mean(infected == 1) * 100,
                   mean_ttx = mean(TTX_mg),
                   ttx_se = sd(TTX_mg, na.rm = TRUE) / sqrt(n()),
                   sample_size = n(),
                   prev_sd = sqrt(sample_size * (Proportion / 100) * (1 - Proportion / 100))
  )

tt_bd_prev_ttx_plot <- tt.loc.prev %>%  
  ggplot(aes(mean_ttx, Proportion, color = Location)) + 
  geom_errorbarh(aes(x = mean_ttx, xmin = mean_ttx - ttx_se, 
                     xmax = mean_ttx + ttx_se),
                 size = 1,
                 height = 0,
                 color = 'black') +
  geom_errorbar(aes(y = Proportion,
                    ymin = Proportion - prev_sd,
                    ymax = Proportion + prev_sd),
                size = 1,
                width = 0,
                color = 'black') +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = 'Taricha torosa',
       x = 'Mean TTX (mg)',
       y = 'Infection Prevalence') +
  scale_color_manual(values = tt_colors,
                     labels = c('Cold Creek', 'Gunstock', 'Madonna', 'Muir')) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = 'bold',
                               size = 16),
    axis.text.y = element_text(face = 'bold',
                               size = 14),
    axis.title.x = element_text(face = 'bold',
                                size = 14,
                                hjust = 0.5),
    plot.title = element_blank(),
    axis.title.y = element_text(face = 'bold',
                                size = 14),
    axis.line = element_line(),
    legend.position = c(1,0.5),
    legend.justification = c('right', 'top'),
    legend.box.background = element_rect(),
    legend.title = element_text(face = "bold",
                                size = 20,
                                hjust = 0.5),
    legend.text = element_text(size = 15)
  ) +
  guides(color = guide_legend(title = "Location")) +
  labs(tag = "(D)")

tt_bd_prev_ttx_plot

glmms <- ggarrange(tg_glmm, tt_plot_glmm,
                   ncol = 2,
                   nrow = 1,
                   widths = c(2,2)) 
glmms
ggarrange(glmms, glmms2)

glmm_means <- ggaglmm_means <- ggaglmm_means <- ggarrange(tg_bd_prev_ttx_plot, tt_bd_prev_ttx_plot)
glmm_means
ggarrange(glmms, glmm_means,
          nrow = 2,
          ncol = 2,
          heights = c(0.5, 0.5, 0.5, 0.5),
          widths = c(0.5, 0.5, 0.5, 0.5)
          )

glmms2 <- ggarrange(tt_plot_glmm, tt_bd_prev_ttx_plot,
                    common.legend = T,
                    legend = 'bottom',
                    ncol = 2,
                    nrow = 1,
                    heights = c(2,1))

glmms2
final.plot <- ggarrange(glmms.update, glmms2,
          nrow = 1,
          ncol = 2,
          heights = c(1, 1),
          widths = c(0.5, 0.5))
final.plot <- ggarrange(glmms, glmm_means,
          nrow = 2)
ggsave("figure_4.pdf", plot = final.plot, width = 13, height = 10,
       units = 'in', dpi = 1000)


f# Bd infection status as TTX with location ######______________
tg_loc_glmm <- glmer(infected ~ TTX_mg + (1 | Location), data = loc_bd_pos_tg, family = 'binomial')
tg_loc_glmm_null1 <- glmer(infected ~ (1 | Location), data = loc_bd_pos_tg, family = 'binomial')
summary(tg_loc_glmm)
drop1(tg_loc_glmm, test = 'Chisq')
lmtest::lrtest(tg_loc_glmm, tg_loc_glmm_null1)

tt_loc_glmm <- glmer(infected ~ TTX_mg + (1 | Location), data = loc_tt, family = 'binomial')
tt_loc_glmm_null <- glmer(infected ~ (1 | Location), data = loc_tt, family = 'binomial')
summary(tt_loc_glmm)
drop1(tt_loc_glmm, test = 'Chisq')
lrtest(tt_loc_glmm, tt_loc_glmm_null)
#neither are significant ####################3


#__________________tt_loc_glmm_full#_________________________alpha div______________________
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

tt_pielou_evenness_lm <- lmer(pielou_evenness ~ TTX_mg + (1 | Location), data = tt_data2)
summary(tt_pielou_evenness_lm)
anova(tt_pielou_evenness_lm)
################################################################


tt_asv_kw <- kruskal.test(tt_data2$observed_features, tt_data2$Location)
tt_asv_kw

tt_shan_kw <- kruskal.test(tt_data2$shannon_entropy, tt_data2$Location)
tt_shan_kw

tt_faith_kw <- kruskal.test(tt_data2$faith_pd, tt_data2$Location)
tt_faith_kw

tt_pielou_kw <- kruskal.test(tt_data2$pielou_evenness, tt_data2$Location)
tt_pielou_kw


#Bd______________ #####################################33
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

##########################################
#removing na rows because they mess up the drop1 model
tg.data3 <- tg_data2 %>% 
  filter(!is.na(observed_features))
tg_obs_feat_inf <- glmer(infected ~ observed_features + (1 | Location), data = tg.data3, family = 'binomial')
tg_obs_feat_inf_null <- glmer(infected ~ (1 | Location), data = tg.data3, family = 'binomial')
summary(tg_obs_feat_inf)
drop1(tg_obs_feat_inf, test = 'Chisq')
lmtest::lrtest(tg_obs_feat_inf, tg_obs_feat_inf_null)

tg_shan_feat_inf <- glmer(infected ~ shannon_entropy + (1 | Location), data = tg.data3, family = 'binomial')
tg_shan_feat_inf_null <- glmer(infected ~ (1 | Location), data = tg.data3, family = 'binomial')
summary(tg_shan_feat_inf)
drop1(tg_shan_feat_inf, test = "Chisq")
lmtest::lrtest(tg_shan_feat_inf, tg_shan_feat_inf_null)


tg_faith_feat_inf <- glmer(infected ~ faith_pd + (1 | Location), data = tg.data3, family = 'binomial')
tg_faith_feat_inf_null <- glmer(infected ~ (1 | Location), data = tg.data3, family = 'binomial')
summary(tg_faith_feat_inf)
drop1(tg_faith_feat_inf, test = "Chisq")
lmtest::lrtest(tg_faith_feat_inf, tg_faith_feat_inf_null)

tg_pielou_feat_inf <- glmer(infected ~ pielou_evenness + (1 | Location), data = tg.data3, family = 'binomial')
tg_pielou_feat_inf_null <- glmer(infected ~ (1 | Location), data = tg.data3, family = 'binomial')
summary(tg_pielou_feat_inf)
drop1(tg_pielou_feat_inf, test = "Chisq")
lmtest::lrtest(tg_pielou_feat_inf, tg_pielou_feat_inf_null)

#remove nas from data for TATO
tt.data3 <- tt_data2 %>% 
  filter(!is.na(observed_features))
tt_obs_feat_inf <- glmer(infected ~ observed_features + (1 | Location), data = tt.data3, family = 'binomial')
tt_obs_feat_inf_null <- glmer(infected ~ (1 | Location), data = tt.data3, family = 'binomial')
summary(tt_obs_feat_inf)
drop1(tt_obs_feat_inf, test = "Chisq")
lmtest::lrtest(tt_obs_feat_inf, tt_obs_feat_inf_null)

tt_shan_feat_inf <- glmer(infected ~ shannon_entropy + (1 | Location), data = tt.data3, family = 'binomial')
tt_shan_feat_inf_null <- glmer(infected ~ (1 | Location), data = tt.data3, family = 'binomial')
summary(tt_shan_feat_inf)
drop1(tt_shan_feat_inf, test = "Chisq")
lmtest::lrtest(tt_shan_feat_inf, tt_shan_feat_inf_null)

tt_faith_feat_inf <- glmer(infected ~ faith_pd + (1 | Location), data = tt.data3, family = 'binomial')
tt_faith_feat_inf_null <- glmer(infected ~ (1 | Location), data = tt.data3, family = 'binomial')
summary(tt_faith_feat_inf)
drop1(tt_faith_feat_inf, test = "Chisq")
lmtest::lrtest(tt_faith_feat_inf, tt_faith_feat_inf_null)

tt_pielou_feat_inf <- glmer(infected ~ pielou_evenness + (1 | Location), data = tt.data3, family = 'binomial')
tt_pielou_feat_inf_null <- glmer(infected ~ (1 | Location), data = tt.data3, family = 'binomial')
summary(tt_pielou_feat_inf)
drop1(tt_pielou_feat_inf, test = "Chisq")
lmtest::lrtest(tt_pielou_feat_inf, tt_pielou_feat_inf_null)



#Bd locations______________ ############################
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



#bd_pos indiv________________________________________ ############################
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




#Figures for this #####
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

#alpha and the rest ############

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




#sex mixed ##### 
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


  
#mass lmers #################33
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
anova(tt_obs_mass, type = 'III')


tt_shan_mass <- lmer(shannon_entropy ~ mass_g + (1 | Location), data = tt_data2)
summary(tt_shan_mass)
anova(tt_shan_mass, type = 'III')

tt_faith_mass <- lmer(faith_pd ~ mass_g + (1 | Location), data = tt_data2)
summary(tt_faith_mass)
anova(tt_faith_mass, type = 'III')

tt_pielou_mass <- lmer(pielou_evenness ~ mass_g + (1 | Location), data = tt_data2)
anova(tt_pielou_mass, type = 'III')
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
  labs(tag = '(A)')

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
  labs(tag = '(B)')


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
  labs(tag = '(D)')


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
  labs(tag = '(D)')

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
  labs(tag = '(E)')


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
  labs(tag = '(F)')


part_3

ggarrange(part_1, part_2 ,part_3)




nmds_plots <- ggarrange(tg_part_1, tg_part_2, tg_part_3, ncol = 3,
                        part_1, part_2, part_3, nrow = 2)

nmds_plots    

ggsave('Figure_7.pdf', nmds_plots, units = 'in', dpi = 1000, height = 12, width = 18)







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
  labs(tag = '(A)')

merged_pcoa

#################################################
#Barplots
setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/')
relfreq_1 <- read.table('relfreq.tsv', sep = '\t', header = F)
relfreq_1 <- relfreq_1 %>% 
  column_to_rownames(var = 'SampleID')
metadata <- read.csv('2200_full_metadata.csv', header = T)
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
  labs(tag = '(B)')
tagr_barplot

barplots <- ggarrange(tagr_barplot, taxa_plot)
merged_fig_5 <- ggarrange(merged_pcoa, barplots)
merged_fig_5
ggsave('figure_5.pdf', merged_fig_5, units = 'in', dpi = 1000, height = 10, width = 18)
