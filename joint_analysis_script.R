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




#load datasets
setwd('~/Documents/Thesis_Part_2/Newt_TTX_Project/Joint_analysis/')

#This dataset includes all alpha diversity metrics I will work with as time goes on, it includes no rows that had NAs (in either the diversity or anywhere else)
df <- read.csv('2200_full_metadata.csv',
               header = TRUE)


df$Bd_Real <- ifelse(df$infected > 0, 'Positive', 
                     ifelse(df$infected == 0, 'Negative', 'Negative'))
#species subsets
df_torosa <- subset(df2, Species == 'Taricha_torosa')
df_granulosa <- subset(df, Species == 'Taricha_granulosa')

df2<- df %>%
  filter(!Location%in%c("Bolinger","Crocker"))

tt_df2<- df2 %>%
  subset(Species == 'Taricha_torosa')
tg_df2 <- df2 %>%
  subset(Species == 'Taricha_granulosa')

#____________________________________________________________________________________________________________________________________________________________
#normal distributions of our diversity metrics?
#exploratory statistics

#what this code does is defines variables from the datasets I want to work with and tests them all in one go from all the datasets rather than having to do them individually. this can evven be expanded
variables_to_test <- c('infected','Log_bd','TTX_mg','mass_g',
                       'TTX_Richness', 'TotalTTX', 'faith_pd',
                       'shannon_entropy','pielou_evenness','observed_features',
                       'TTX_RelativeAbundance')
data_frames <- list(df2 = df2, 
                    tt_df2 = tt_df2, 
                    tg_df2 = tg_df2)
results_df <- data.frame(variable = character(),
                         dataframe = character(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

for (variable_name in variables_to_test) {
  cat("Shapiro-Wilk test for", variable_name, ":\n")
  
  for (df_name in names(data_frames)) {
    result <- shapiro.test(data_frames[[df_name]][[variable_name]])
    cat("From", df_name, ":\n")
    print(result)
    cat("\n")
    results_df <- rbind(results_df, data.frame(variable = variable_name,
                                               dataframe = df_name,
                                               p_value = result$p.value))
  }
}
write.csv(results_df, 
          "shapiro_wilk_results.csv", 
          row.names = FALSE)
#significant from the shapiro-wilks test
#infected: df2, tt_df2, tg_df2
#log_bd: df2, tt_df2, tg_df2
#ttx_mg: df2, tt_df2, tg_df2
#mass_g: df2
#observed_features: df2, tt_df2, tg_df2
#TTX_asvs
#antifungal_asvs:

#not significant
#mass_g: tt_df2, tg_df2
tt_bd_pos <- tt_df2 %>%
  subset(infected == 1)
tg_bd_pos <- tg_df2 %>%
  subset(infected ==1)
tg_bd_neg <- tg_df2 %>%
  subset(infected == 0 & !SampleID == 'FCSP-05')
tt_sites <- tt_df2 %>%
  filter(Location%in%c('Cold_Creek', 'Muir'))
tg_sites <- tg_df2 %>%
  filter(!Location%in%c('LSSP', 'FCSP'))

print(mean(tg_bd_pos$TTX_mg))
print(mean(tg_bd_neg$TTX_mg))
ttxttx <- lm(tg_df2$Log_bd ~ tg_df2$TTX_mg)
summary(ttxttx)

kruskal.test(tg_bd_pos$Log_bd, tg_bd_pos$Location)

bdposlin <- lm(tg_bd_pos$Log_bd ~ tg_bd_pos$TTX_mg)
bdposlin2 <- lm(tg_sites$Log_bd ~ tg_sites$TTX_mg)
summary(bdposlin2)
summary(bdposlin)
bdposlin2 <- lm(tt_bd_pos$Log_bd ~ tt_bd_pos$TTX_mg)
new_enw_new <- lm(tt_sites$Log_bd ~ tt_sites$TTX_mg)
summary(new_enw_new)
summary(bdposlin2)

kruskal.test(tt_sites$Log_bd, tt_sites$Location)
kruskal.test(tg_sites$Log_bd, tg_sites$Location)
kruskal.test(tg_df2$Log_bd, tg_df2$Location)
kruskal.test(tt_sites$infected, tt_sites$Location)
kruskal.test(tt_df2$infected, tt_df2$Location)

kruskal.test(tt_bd_pos$Log_bd, tt_bd_pos$Location)
kruskal.test(tt_bd_pos$infected)

kruskal.test(tg_df2$Log_bd, tg_df2$Location)
kruskal.test(tt_df2$Log_bd, tt_df2$Location)
kruskal.test(tg_df2$infected, tg_df2$Location)
kruskal.test(tt_df2$infected, tt_df2$Location)
kruskal.test(tg_df2$TTX_mg, tg_df2$sex)
kruskal.test(tt_df2$TTX_mg, tt_df2$sex)
kruskal.test(df2$Log_bd, df2$Location)
kruskal.test(df2$infected, df2$Location)
kruskal.test(df2$TTX_mg, df2$Species)
print(max(df$TTX_mg, na.rm = T))
print(mean(tg_df2$TTX_mg, na.rm = T))
print(min(df$TTX_mg))
new_tt <- df %>%
  subset(Species == 'Taricha_torosa')
print(mean(new_tt$TTX_mg))

glm1 <- glm(data = tg_df2, infected ~ TTX_mg, family = 'binomial')
summary(glm1)
glm12 <- glm(data = tg_sites, infected ~ TTX_mg, family = 'binomial')
summary(glm12)
glm122 <- glm(data = tg_bd_pos, infected ~ TTX_mg, family = 'binomial')
summary(glm122)
glm2 <- glm(data = tt_df2, infected ~ TTX_mg, family = 'binomial')
summary(glm2)
glm22 <- glm(data = tt_sites, infected ~ TTX_mg, family = 'binomial')
summary(glm22)
#Shadf2#Shannnon______________________________________________________________________________
response <- c('shannon_entropy', 'faith_pd', 
              'pielou_evenness', 'observed_features', 
              'mass_g', 'TTX_mg', 'Log_bd',
              'TTX_RelativeAbundance')
predictor <- c('TTX_Real', 'Location', 
               'Bd_Real', 'sex', 'Species', 'mass_g')
data_frames_1 <- list(df2 = df2, tt_df2 = tt_df2, tg_df2 = tg_df2)

multi_kw_results_df <- data.frame(predictor = character(),
                                  response = character(),
                                  dataframe = character(),
                                  chi_squared = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (predictor_var in predictor) {
  for (response_var in response) {
    for (df_name in names(data_frames_1)) {
      predictor_data <- data_frames_1[[df_name]][[predictor_var]]
      response_data <- data_frames_1[[df_name]][[response_var]]
      
      multi_kw_result <- kruskal.test(response_data ~ predictor_data)
      
      multi_kw_results_df <- rbind(multi_kw_results_df, data.frame(predictor = predictor_var,
                                                                   response = response_var,
                                                                   dataframe = df_name,
                                                                   chi_squared = multi_kw_result$statistic,
                                                                   p_value = multi_kw_result$p.value))
    }
  }
}

write.csv(multi_kw_results_df, "multi_kruskal_wallis_results.csv", row.names = FALSE)

multi_kw_results_df <- multi_kw_results_df[order(multi_kw_results_df$response, 
                                                 multi_kw_results_df$predictor), ]
table1_txt <- knitr::kable(multi_kw_results_df, format = "markdown", 
                           col.names = c("Response Variable", 
                                         "Predictor Variable", 
                                         "Data Frame",
                                         'Chi-Squared',
                                         "p-value"))
write.table(table1_txt, 
            "multi_overall_kruskal_wallis_results.txt", 
            quote = FALSE, 
            row.names = FALSE)

#_______________________________________________________________________________--

response2 <- c('Log_bd','TTX_RelativeAbundance', 'infected')
predictor2 <- c('Species')
data_frames_2 <- list(df2 = df2, bd_pos = bd_pos)

multi_kw_results_df2 <- data.frame(predictor = character(),
                                  response = character(),
                                  dataframe = character(),
                                  chi_squared = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (predictor_var in predictor2) {
  for (response_var in response2) {
    for (df_name in names(data_frames_2)) {
      predictor_data2 <- data_frames_2[[df_name]][[predictor_var]]
      response_data2 <- data_frames_2[[df_name]][[response_var]]
      
      multi_kw_result2 <- kruskal.test(response_data2 ~ predictor_data2)
      
      multi_kw_results_df2 <- rbind(multi_kw_results_df2, data.frame(predictor = predictor_var,
                                                                   response = response_var,
                                                                   dataframe = df_name,
                                                                   chi_squared = multi_kw_result2$statistic,
                                                                   p_value = multi_kw_result2$p.value))
    }
  }
}

write.csv(multi_kw_results_df2, "spp_kruskal_wallis_results.csv", row.names = FALSE)

multi_kw_results_df <- multi_kw_results_df[order(multi_kw_results_df$response, 
                                                 multi_kw_results_df$predictor), ]
table1_txt <- knitr::kable(multi_kw_results_df, format = "markdown", 
                           col.names = c("Response Variable", 
                                         "Predictor Variable", 
                                         "Data Frame",
                                         'Chi-Squared',
                                         "p-value"))
write.table(table1_txt, 
            "multi_overall_kruskal_wallis_results.txt", 
            quote = FALSE, 
            row.names = FALSE)

#_______________________________________________________________________________
kw_2_results_df <- data.frame(predictor = character(),
                            response = numeric(),
                            chi_sqquared = numeric(),
                            p_value = numeric(),
                            stringsAsFactors = FALSE)
kw_predictor_2 <- 'Species'
kw_response_2 <- c('shannon_entropy', 'faith_pd', 
                'pielou_evenness', 'observed_features', 
                'mass_g', 'TTX_mg')
kw_species <- for(response_vari in kw_response_2){
  species_kw <- kruskal.test(df2[[response_vari]] ~ df2$Species, na.action = na.omit)
  print(species_kw)
  kw_2_results_df <- rbind(kw_2_results_df, data.frame(predictor = kw_predictor_2,
                                                   response = response_vari,
                                                   dataframe = 'df2',
                                                   p_value = species_kw$p.value))
}
write.csv(kw_2_results_df, 'species_kw_results.csv', row.names = FALSE)
kw_2_results_df <- kw_2_results_df[order(kw_2_results_df$response, 
                                         kw_2_results_df$predictor), ]
table_txt <- kable(kw_2_results_df, format = "markdown", 
                   col.names = c("Response Variable", 
                                 "Predictor Variable", 
                                 "Data Frame", 
                                 "p-value"))
write.table(table_txt, 
            "species_kruskal_wallis_results.txt", 
            quote = FALSE, 
            row.names = FALSE)

print(mean(tt_df2$TTX_mg))
print(max(tt_df2$TTX_mg, na.rm = T))
print(min(tt_df2$TTX_mg, na.rm = T))
mean_tg <- mean(tg_df2$TTX_mg, na.rm = TRUE)
print(mean_tg)
min_tg <- min(tg_df2$TTX_mg, na.rm = T)
print(min_tg)
max_tg <- max(tg_df2$TTX_mg, na.rm = T)
print(max_tg)


df2$Species <- as.factor(df2$Species)
ktest <- kruskal.test(data = tt_df2, TTX_mg ~ Location, na.action = na.omit)
ktest
tg_ktest <- kruskal.test(data = tg_df2, TTX_mg ~ Location, na.action = na.omit)
tg_ktest

print(sum(df$infected))
print(sum(tt_df2$infected))
print(sum(tg_df2$infected))


basic_plot <- df2 %>%
  ggplot(aes(Species, TTX_mg, fill = Species))+
  geom_boxplot()+
  theme_minimal()+
  guides(fill=FALSE)

basic_plot



#TTX linear models_____________________________________________________________
bd_pos <- df2 %>%
  subset(infected == 1)

tg_bdpos <- bd_pos %>%
  subset(Species == 'Taricha_granulosa')

tt_bdpos <- bd_pos %>%
  subset(Species == 'Taricha_torosa')

response_variables <- c(
  'shannon_entropy', 'faith_pd', 
  'pielou_evenness', 'observed_features', 
  'mass_g', 'TTX_mg', 'Log_bd',
  'antifungal_ASVs', 'antifungal_reads',
  'Propor_TotalAntiFungal', 'TTX_RelativeAbundance',
  'TTX_Richness', 'TotalTTX'
)

predictor_variables <- c(
  'shannon_entropy', 'faith_pd', 
  'pielou_evenness', 'observed_features', 
  'mass_g', 'TTX_mg', 'Log_bd',
  'antifungal_ASVs', 'antifungal_reads',
  'Propor_TotalAntiFungal', 'TTX_RelativeAbundance',
  'TTX_Richness', 'TotalTTX'
)
variables_indf <- response_variables %in% colnames(df)
print(variables_indf)

lm_data_frames <- list(df2 = df2, tt_df2 = tt_df2, tg_df2 = tg_df2,
                       bd_pos = bd_pos, tg_bdpos = tg_bdpos, tt_bdpos = tt_bdpos)

model_results_df <- data.frame(response = character(),
                               predictor = character(),
                               dataframe = character(),
                               coefficients = numeric(),
                               t_value = numeric(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)

for(lm_df_names in names(lm_data_frames)){
  lm_data_frame <- lm_data_frames[[lm_df_names]]
  for(lm_response_var in response_variables){
    for(lm_predictor_var in predictor_variables){
      linear_model <- lm(lm_data_frame[[lm_response_var]] ~ lm_data_frame[[lm_predictor_var]])
      coeficient_value <- coef(summary(linear_model))[2, 1]
      t_value <- coef(summary(linear_model))[2, 3]
      p_value <- coef(summary(linear_model))[2, 4]
      model_results_df <- rbind(model_results_df, data.frame(response = lm_response_var,
                                                             predictor = lm_predictor_var,
                                                             dataframe = lm_df_names,
                                                             coefficients = coeficient_value,
                                                             t_value = t_value,
                                                             p_value = p_value))
    }
  }
}
write.csv(model_results_df, "linear_regression_results.csv", row.names = FALSE)
lm_table <- kable(model_results_df, format = "markdown", 
                         col.names = c("Response Variable", "Predictor Variable", "Data Frame", "Coefficients", "p-value"))
writeLines(c(lm_table), "linear_regression_results.txt")

#zoopsore load graph_______________________________________________________________
yexpres <- expression(bold("Zoospore" ~ "Equivalents" ~ "(log"['10']*"+1)"))
zooplot <- df2 %>%
  ggplot(aes(TTX_mg, Log_bd))+
  geom_smooth(method = 'lm',
                  fill = 'grey',
                  aes(color = Species))+
  geom_point(aes(fill = Species,
                 color = Species),
             size = 2,
             show.legend = TRUE)+
  scale_color_manual(values = c(
    "#CC79A7", "#56B4E9")) + 
  theme_classic() +
  xlab('TTX (mg)') +
  ylab(yexpres) +
  #labs(tag = 'c)') +
  ggtitle('Linear Model of Infection Intensity and TTX Concentration')+
  theme(
    axis.title = element_text(face = 'bold', size = 14),
    axis.title.x = element_text(face = 'bold'),
    axis.text = element_text(face = 'bold', size = 14),
    plot.title = element_text(face = 'bold', size = 20, hjust = 0.5),
    title = element_text(face = "bold", size = 12)
  )

zooplot





zoospore_ttx_lm <- df2 %>%
  ggplot(aes(TTX_mg, Log_bd,
             color = Species)) +
  geom_smooth(method = 'lm', 
              fill = 'grey', 
              color = 'red') +
  geom_point(aes(fill = Species, color = Species), 
             size = 4, 
             show.legend = TRUE)+
  scale_color_manual(values = c(
    "green", "#56B4E9")) +
    theme_classic() +
    xlab('TTX (mg)') +
    ylab('Zoospore Equivalents') +
    labs(tag = 'c)') +
    theme(
      axis.title = element_text(face = 'bold', size = 14),
      axis.title.x = element_text(face = 'bold'),
      axis.text = element_text(face = 'bold', size = 14),
      plot.title = element_text(face = 'bold', size = 20, hjust = 0.5)
    )

zoospore_ttx_lm
  
theme_classic() +
  xlab('TTX (mg)') +
  ylab('Zoospore Equivalents') +
  labs(tag = 'c)') +
  theme(
    axis.title = element_text(face = 'bold', size = 14),
    axis.title.x = element_text(face = 'bold'),
    axis.text = element_text(face = 'bold', size = 14),
    plot.title = element_text(face = 'bold', size = 20, hjust = 0.5)
  )+
  guides(size = FALSE)

zoospore_ttx_lm




#_GLM and Graph______________________________________________________________
overall_infection_glm <- glm(infected ~ TTX_mg, data = df2, family = binomial)
summary(overall_infection_glm)
tt_infection_glm <- glm(infected ~ TTX_mg, data = tt_df2, family = binomial)
summary(tt_infection_glm)
tg_infection_glm <- glm(infected ~ TTX_mg, data = tg_df2, family = binomial)
summary(tg_infection_glm)

overall_infection_glm <- glm(infected ~ TTX_mg, data = bd_pos_sites, family = binomial)
summary(overall_infection_glm)
tt_infection_glm <- glm(infected ~ TTX_mg, data = bd_tt_df2, family = binomial)
summary(tt_infection_glm)
tg_infection_glm <- glm(infected ~ TTX_mg, data = bd_tg_df2, family = binomial)
summary(tg_infection_glm)

overall_infection_lm <- lm(Log_bd ~ TTX_mg, data = bd_pos_sites)
summary(overall_infection_glm)
tt_infection_lm <- lm(Log_bd ~ TTX_mg, data = bd_tt_df2)
summary(tt_infection_lm)
tg_infection_lm <- lm(Log_bd ~ TTX_mg, data = bd_tg_df2)
summary(tg_infection_lm)



















#torosa
shapiro.test(torosa_full$mass.g.)
mean(torosa_full$mass.g.)
print(std.error(torosa_full$mass.g.))
torosa_mass_aov <- aov(data = torosa_full, mass.g. ~ Location)
summary(torosa_mass_aov)

bd_zoo <- kruskal.test(data = df2, Log_bd ~ Location)
bd_zoo
tt_zoo <- kruskal.test(data = tt_df2, Log_bd ~ Location)
tt_zoo
tg_zoo <- kruskal.test(data = tg_df2, Log_bd ~ Location)
tg_zoo

shapiro.test(df_base$mass.g.)
mean(df_base$mass.g.)
print(std.error(df_base$mass.g.))
mass_aov <- aov(df_base$mass.g. ~ df_base$Location)
summary(mass_aov)

torosa_males <- subset(torosa_full, sex == 'M')
torosa_females <- subset(torosa_full, sex == 'F')

print(mean(torosa_males$mass.g.))
print(std.error(torosa_males$mass.g.))
print(mean(torosa_females$mass.g.))
print(std.error(torosa_females$mass.g.))

shapiro.test(torosa_full$mass.g.)
sex_mass_aov <- aov(data = torosa_full, mass.g. ~ sex)
summary(sex_mass_aov)

# more LMs
torosa_lm_TTX_mass <- lm(data = torosa_alpha, TTX.mg. ~ mass.g. * sex)
summary(torosa_lm_TTX_mass)

granulosa_lm_ttx_mass <- lm(data = granulosa_alpha, TTX.mg. ~ mass.g. * sex)
summary(granulosa_lm_ttx_mass)

full_lm_ttx_mass <- lm(data = overall_dataset, TTX.mg. ~ mass.g. * sex)
summary(full_lm_ttx_mass)

#lms 
torosa_lm_TTX_mass_location <- lm(data = torosa_alpha, TTX.mg. ~ mass.g. * sex * Location)
summary(torosa_lm_TTX_mass_location)

granulosa_lm_ttx_mass_locaion <- lm(data = granulosa_alpha, TTX.mg. ~ mass.g. * sex * Location)
summary(granulosa_lm_ttx_mass_locaion)

full_lm_ttx_mass_location <- lm(data = overall_dataset, TTX.mg. ~ mass.g. * sex * Location)
summary(full_lm_ttx_mass_location)
plot(full_lm_ttx_mass_location)

#dredge
#mine_five <- glm(Bd_Presence ~ TTX.mg. * Location * Mass.g. * SVL.mm., data = overall_dataset, family = binomial)
#summary(mine_five)


#tetrodotoxin analysis
shapiro.test(torosa_full$TTX.mg.)
shapiro.test(torosa_edited_df$TTX.mg.)
#looks like TTX concentration is not normally distributed, therefore kw?

ttx_torosa_site <- kruskal.test(data = torosa_edited_df, TTX.mg. ~ Location)
ttx_torosa_site
ttx_torosa_sex <- kruskal.test(data = torosa_edited_df, TTX.mg. ~ sex)
ttx_torosa_sex

granulosa_bd_dataset <- read.csv('TBJALH_TTX_DATA_CLEANED.csv', header = TRUE)

shapiro.test(overall_dataset$TTX.mg.)
ttx_overall_site <- kruskal.test(data = overall_dataset, TTX.mg. ~ Location)
ttx_overall_site
pairwise.wilcox.test(overall_dataset$TTX.mg., overall_dataset$Location, p.adjust.method = 'BH')
ttx_torosa_site <- kruskal.test(data = torosa_edited_df, TTX.mg. ~ Location)
ttx_granulosa_site <- kruskal.test(data = granulosa_full, TTX.mg. ~ Location)

overall_ttx_plot <- ggplot(overall_dataset, aes(Location, TTX.mg., fill = Location))+
  geom_boxplot()+
  theme_classic()+
  ylab('TTX (mg)')+
  xlab('')+
  ggtitle('Tetrodotoxin concentration by location')+
  scale_fill_brewer(palette = 'Spectral')+
  scale_x_discrete(labels=c("Cold Creek","Gunstock","Madanna","Muir","FCSP","TMC","LSSP","GRP", "WBC", "SLNF"))+
  theme(axis.title = element_text(face = 'bold', size = 14), 
        plot.title = element_text(face = 'bold', size = 25, hjust = 0.5),
        axis.text.x = element_text(face = 'bold', size = 14))

overall_ttx_plot


speciated_plot <- overall_dataset%>%
  group_by(Species)%>%
  ggplot(aes(Location, TTX.mg., fill = Location))+
  geom_boxplot()+
  theme_classic()+
  ylab('TTX (mg)')+
  xlab('')+
  ggtitle('Tetrodotoxin concentration by location')+
  scale_fill_brewer(palette = 'Spectral')+
  theme(axis.title = element_text(face = 'bold', size = 14), 
        plot.title = element_text(face = 'bold', size = 25, hjust = 0.5),
        axis.text.x = element_text(face = 'bold', size = 14))+
  facet_wrap('Species', scales = 'free_x')
speciated_plot

tt_cmprs <- list(c('Cold_Creek', 'Gunstock', 'Madanna', 'Muir'),
                 c('Gunstock', 'Madanna', 'Muir'), 
                 c('Madanna', 'Muir'))
Tt_speciated_ttx <- torosa_edited_df%>%
  ggplot(aes(x=Location,y=TTX.mg.,fill=Location))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_brewer(palette="Spectral")+
  theme_bw()+
  ylab("TTX (mg)")+
  xlab('')+
  scale_x_discrete(labels=c("Cold Creek","Gunstock","Madanna","Muir"))+
  ggtitle('Taricha torosa')+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  geom_signif(comparisons = list(c('Cold_Creek', 'Gunstock'),
                                 c('Cold_Creek', 'Madonna'),
                                 c('Cold_Creek', 'Muir'),
                                 c('Gunstock', 'Madonna'),
                                 c('Gunstock', 'Muir'),
                                 c('Madonna', 'Muir')),
              map_signif_level = TRUE,
              y_position = c(2, 
                             2.2, 
                             2.4, 
                             2.6, 
                             2.8, 
                             3.0))+
  stat_compare_means(label.y = 3.5, size = 6)
Tt_speciated_ttx

#There's a type here and I'm not sure what to do about it but it works all the same
tg_speciated_ttx <- granulosa_full %>%
  ggplot(aes(Location, TTX.mg., fill = Location))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_brewer(palette="Accent")+
  theme_bw()+
  ylab("TTX (mg)")+
  xlab('')+
  scale_x_discrete(labels=c("FCSP","LSSP","SLNF","GRP", "WBC", "TMC"))+
  ggtitle('Taricha granulosa')+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  geom_signif(comparisons = list(c('FSCP', 'LSSP'),
                                 c('FSCP', 'SLNF'),
                                 c('FSCP', 'GRP'),
                                 c('FSCP', 'WBC'),
                                 c('FSCP', 'TMC'),
                                 c('LSSP', 'SLNF'),
                                 c('LSSP', 'GRP'),
                                 c('LSSP', 'WBC'),
                                 c('LSSP', 'TMC'),
                                 c('SLNF', 'GRP'),
                                 c('SLNF', 'WBC'),
                                 c('SLNF', 'TMC'),
                                 c('GRP', 'WBC'),
                                 c('GRP', 'TMC'),
                                 c('WBC', 'TMC')),
              map_signif_level = TRUE,
              y_position = c(7,
                             7.5,
                             8,
                             8.5,
                             9,
                             9.5,
                             10,
                             10.5,
                             11,
                             11.5,
                             12,
                             12.5,
                             13,
                             13.5,
                             14))+
  stat_compare_means(label.y = 16, label.x = 1.5, size = 6)
tg_speciated_ttx



ggarrange(Tt_speciated_ttx, tg_speciated_ttx)
#annotate_figure(arranged, bottom = text_grob('Location', face = 'bold', size = 20, vjust = 0.5), top = text_grob('Phylogenetic Diversity Between Species', size = 20, face = 'bold'))




pairwise.wilcox.test(torosa_edited_df$TTX.mg., torosa_edited_df$Location)
wilcox.test(torosa_edited_df$Location, torosa_edited_df$TTX.mg.)
wilcox.test(torosa_edited_df$TTX.mg.[which(torosa_edited_df$Location == 'Cold_Creek')], 
            torosa_edited_df$TTX.mg.[which(torosa_edited_df$Location == 'Gunstock')])
# Tt_speciated_ttx


print(range(torosa_full$TTX.mg.))

library(sjPlot)
library(sjmisc)
#Is bd infection prevalence correlated to TTX concentration?

new_trial2 <- glm(Bd_Presence ~ Whole_TTX, data = df_diversity, family = binomial)
summary(new_trial2)
ggplot(df_diversity, aes(Whole_TTX, Bd_Presence))+
  geom_smooth(method = 'glm', method.args = list(family = binomial), fill = 'grey', color = 'black')+
  theme_classic()+
  xlab('TTX (mg)')+
  ylab('Bd Presence')+
  theme(axis.title = element_text(face = 'bold', size = 14))+
  ggtitle('Tg GLM, P = 0.0641')+
  theme(plot.title = element_text(face = 'bold', size = 20, hjust = 0.5))
#it does not appear so

ggarrange()

#zoospore load and TTX concentration
zoospore_estimates <- log(granulosa_bd_dataset$ZE.avg.Bd + 1)
zoospore_lm <- lm(zoospore_estimates ~ Whole_TTX, data = granulosa_bd_dataset)
summary(zoospore_lm)
zoospore_lm_fr <- abline(lm(zoospore_estimates ~ Whole_TTX, data = granulosa_bd_dataset))
plot(zoospore_estimates, granulosa_bd_dataset$Whole_TTX)

zoospore_ttx_lm <- ggplot(granulosa_bd_dataset, aes(Whole_TTX, zoospore_estimates))+
  geom_smooth(method = 'lm', fill = 'grey', color = 'red')+
  geom_point(aes(shape = Location, fill = Location, color = Location), size = 2)+
  theme_classic()+
  scale_shape_manual(values = c(15,16,17,18,19,25))+
  xlab('TTX (mg)')+
  ylab('Zoospore Equivalents')+
  theme(axis.title = element_text(face = 'bold', size = 14))+
  ggtitle('P = 0.0301')+
  theme(plot.title = element_text(face = 'bold', size = 20, hjust = 0.5))+
  guides(size = FALSE, color = guide_legend('Location'))

granulosa_bd_dataset$Bd_load <- log(granulosa_bd_dataset$ZE.avg.Bd + 1)
bd_loc_plot <- ggplot(data = granulosa_bd_dataset, aes(x = Location, y = Bd_load, fill = Location))+
  geom_boxplot()+
  theme_classic()+
  scale_fill_brewer(palette = 'Accent')
bd_loc_plot
library(RColorBrewer)
library(colorspace)

ggarrange(bd_loc_plot, tg_speciated_ttx, zoospore_ttx_lm)


#Bd positive sites
bd_data$zoospore_estimate <- log(bd_data$ZE.avg.Bd + 1)
Bd_pos_df <- subset(bd_data, Site %in% (c('GRP', 'SLNF', 'TMC', 'WBC')))

#is Bd load different among sites with only sites with positive bd infection confirmed?
shapiro.test(bd_data$zoospore_estimate)
shapiro.test(Bd_pos_df$zoospore_estimate)

kruskal.test(bd_data$zoospore_estimate, bd_data$Site)
kruskal.test(Bd_pos_df$zoospore_estimate, Bd_pos_df$Site)
overall_bd_anova <- aov(data = bd_data, zoospore_estimate ~ Site)
summary(overall_bd_anova)
TukeyHSD(overall_bd_anova)

#it appears that they are significantly different but not quite as much with tukeysHSD

#is TTX different among the Bd positive only sites?
site_ttx_anova <- aov(Bd_pos_df$Whole_TTX ~ Bd_pos_df$Site)
summary(site_ttx_anova)
TukeyHSD(site_ttx_anova)
#TTX concentration is significantly different among Bd positive sites

tt_prcomp1 <- wcmdscale(w_unifrac, eig = TRUE)
ggpcoa(ord = tt_prcomp1, ordata = w_unifrac, groups = torosa_alpha$Location, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  ggtitle('Taricha torosa')+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'italic'))+
  scale_color_manual(labels = torosa_alpha$Location, values = c('red', 'blue', 'green', 'orange', 'black', 'yellow'))+
  guides(shape = FALSE, color = guide_legend('Location'))
plot(prcomp2)

ttx_pcoa <- wcmdscale(tg_w_unifrac_dm, ttx_vec, eig = TRUE)
ttx_vec <- as.numeric(unlist(granulosa_alpha$TTX.mg.))
ttx_vec2 <- as.matrix(ttx_vec)
ggpcoa(ord = ttx_pcoa, ordata = w_unifrac, groups = ttx_vec, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+
  ggtitle('granulosa PCoA TTX')+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'))

prcomp(tg_w_unifrac_dm, scale = TRUE)

full_granulosa_dataset$zoospore_estimate <- log(full_granulosa_dataset$ZE.avg.Bd + 1)
bd_lm <- lm(formula = zoospore_estimate ~ Whole_TTX * Mass_g * Sex, data = full_granulosa_dataset)
summary(bd_lm)


#Microbiome analysis__________________________________________________________________________________________
setwd('microbiome/')
#________________________
#torosa files
bc_qza <- read_qza('tt/bray_curtis/bray_curtis_distance_matrix.qza')$data
jaccard_qza <- read_qza('tt/Jaccard/jaccard_distance_matrix.qza')$data
uwu_qza <- read_qza('tt/unweighted_unifrac/unweighted_unifrac_distance_matrix.qza')$data
w_unifrac_qza <- read_qza('tt/weighted_unifrac/weighted_unifrac_distance_matrix.qza')$data

bc_dm <- as.matrix(bc_qza)
jaccard_dm <- as.matrix(jaccard_qza)
uwu_dm <- as.matrix(uwu_qza)
w_unifrac <- as.matrix(w_unifrac_qza)

tt_taxonomy <- read_qza('tt/rarefied3800_table.qza')$data

Taxa <- tt_taxonomy%>%
  mutate(Taxon=gsub("[][]", "", Taxon))%>%
  parse_taxonomy(.)%>%
  rownames_to_column("featureid")%>%
  select(-c(Kingdom,Phylum))

tt_prcomp1 <- wcmdscale(w_unifrac, eig = TRUE)
ggpcoa(ord = tt_prcomp1, ordata = w_unifrac, groups = torosa_alpha$Location, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  ggtitle('Taricha torosa')+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'italic'))+
  scale_color_manual(labels = torosa_alpha$Location, values = c('red', 'blue', 'green', 'orange', 'black', 'yellow'))+
  guides(shape = FALSE, color = guide_legend('Location'))
plot(prcomp2)

#granulosa files
tg_bc_qza <- read_qza('tg/bray_curtis/bray_curtis_distance_matrix.qza')$data
tg_jaccard_qza <- read_qza('tg/jaccard/jaccard_distance_matrix.qza')$data
tg_uwu_qza <- read_qza('tg/unweighted_unifrac/unweighted_unifrac_distance_matrix.qza')$data
tg_w_unifrac_qza <- read_qza('tg/weighted_unifrac/weighted_unifrac_distance_matrix.qza')$data

tg_bc_dm <- as.matrix(tg_bc_qza)
tg_jaccard_dm <- as.matrix(tg_jaccard_qza)
tg_uwu_dm <- as.matrix(tg_uwu_qza)
tg_w_unifrac_dm <- as.matrix(tg_w_unifrac_qza)

tg_taxonomy <- read_qza('tg/rarefied_table.qza')$data


trial_prcomp <- prcomp(tg_w_unifrac_dm, scale = TRUE)
summary(trial_prcomp)

prcomp2 <- wcmdscale(tg_w_unifrac_dm, eig = TRUE)
tg_wunifrac_pcoa <- ggpcoa(ord = prcomp2, ordata = tg_w_unifrac_dm, groups = granulosa_alpha$Location, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  ggtitle('Weighted Unifrac PCoA')+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  scale_color_manual(labels = granulosa_alpha$Location, values = c('red', 'blue', 'green', 'orange', 'black', 'yellow'))+
  guides(shape = FALSE, color = guide_legend('Location'))
plot(prcomp2)


uwunifrac_prcomp <- wcmdscale(tg_uwu_dm, eig = TRUE)
tg_uwu_pcoa <- ggpcoa(ord = uwunifrac_prcomp, ordata = tg_uwu_dm, groups = granulosa_alpha$Location, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(16,16,16,16,16,16))+
  ggtitle('Unweighted Unifrac PCoA')+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  scale_color_manual(labels = granulosa_alpha$Location, values = c('red', 'blue', 'green', 'orange', 'black', 'yellow'))+
  guides(shape = FALSE, color = guide_legend('Location'))
plot(uwunifrac_prcomp)

ggarrange(tg_wunifrac_pcoa, tg_uwu_pcoa, common.legend = T, legend = 'right')

#Bd PCoA Plot weighted_unifrac
bd_data <- read.csv('final_metadata.csv', header = TRUE)
bd_data$Bd_Presence <- as.character(bd_data$Bd_Presence)
tg_weighted_prcomp <- ggpcoa(ord = prcomp2, ordata = tg_w_unifrac_dm, groups = bd_data$Bd_Presence, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(0,1,24,6,18,16))+
  ggtitle('Weighted Unifrac PCoA Bd Infection')+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))+
  scale_color_manual(labels = c('Negative Infection', 'Positive Infection'), values = c('blue', 'red'))+
  theme_minimal()+
  guides(shape = FALSE, color = guide_legend("Infection Status"))
plot(prcomp2)

#Bd PcoA plot unweighted_unifrac
tg_prcomp_graph <- ggpcoa(ord = uwunifrac_prcomp, ordata = tg_uwu_dm, groups = bd_data$Bd_Presence, spearrow = NULL, spe = FALSE, ellipse = TRUE, cirline = 1)+ 
  scale_shape_manual(values = c(0,1,24,6,18,16))+
  ggtitle('Unweighted Unifrac PCoA')+
  theme(plot.title = element_text(hjust = 1, face = 'bold', size = 20))+
  scale_color_manual(labels = c('Negative Infection', 'Positive Infection'), values = c('blue', 'red'))+
  theme_bw()+
  guides(shape = FALSE, color = guide_legend("Infection Status"))
tg_prcomp_graph

#Tt and Tg pcoa plots
title <- '<i>Taricha granulosa<i> Bd Infection PCoA'
tgrob <- text_grob(title, size = 20) + theme(plot.margin = margin(0,3,0,0, "cm"))
bd_pcoa_tg_plot <- ggarrange(tg_prcomp_graph, tg_weighted_prcomp, common.legend = TRUE, legend = 'right')

annotate_figure(bd_pcoa_tg_plot) 


#_____________________________________________
#TTX and anti-microbial import

packages <- c('qiime2R', 'magrittr', 'purr', 'stringr', 'devtools')
library(qiime2R) #(jbisanz/qiime2R)
library(magrittr)
library(purr)
library(stringr)
library(devtools)


mean_ttx <- aggregate(TTX_mg ~ Location, data = df, FUN = mean)
print(mean_ttx)





















#bd_pos_sites___________________________________________________________________________________

bd_infect_table <- subset(df, infected > 0)
print(table(bd_infect_table$Location))
print(min(bd_infect_table$zoospore_estimates))
print(max(bd_infect_table$zoospore_estimates))
print(mean(bd_infect_table$zoospore_estimates))
sites <- c('Cold_Creek', 'GRP', 'LSSP',
           'Muir', 'SLNF', 'TMC', 'WBC')
bd_pos_sites <- df %>%
  filter(Location %in% sites)
bd_tt_df2 <- no_lssp_bd %>%
  subset(Species == 'Taricha_torosa')
bd_tg_df2 <- no_lssp_bd %>%
  subset(Species == 'Taricha_granulosa')
print(table(no_lssp_bd$Location))

gen_plot <- bd_pos_sites %>%
  ggplot(aes(Location, Log_bd)) +
  geom_boxplot() +
  geom_point()

gen_plot

sites_2 <- c('Cold_Creek', 'GRP','Muir', 
             'SLNF', 'TMC', 'WBC')
no_lssp_bd <- df %>%
  filter(Location %in% sites_2)
gen_plot2 <- no_lssp_bd %>%
  ggplot(aes(Location, Log_bd)) +
  geom_boxplot() +
  geom_point()

gen_plot2


kruskal.test(data = no_lssp_bd, Log_bd ~ Location)
contin_table <- table(no_lssp_bd$Location, no_lssp_bd$infected)
print(contin_table)
chisq.test(contin_table)

kruskal.test(data = bd_tg_df2, Log_bd ~ Location)
tg_contin <- table(bd_tg_df2$Location, bd_tg_df2$infected)
chisq.test(tg_contin)

kruskal.test(data = bd_tt_df2, Log_bd ~ Location)
tt_contin <- table(bd_tt_df2$Location, bd_tt_df2$infected)
chisq.test(tt_contin)

overall_infection_glm






#normal distributions of our diversity metrics?
#exploratory statistics

#what this code does is defines variables from the datasets I want to work with and tests them all in one go from all the datasets rather than having to do them individually. this can evven be expanded
variables_to_test <- c('Log_bd','TTX_mg','mass_g',
                       'TTX_Richness', 'TotalTTX', 'faith_pd',
                       'shannon_entropy','pielou_evenness','observed_features')
data_frames <- list(bd_infect_table = bd_infect_table, 
                    bd_tt_df2 = bd_tt_df2, 
                    bd_tg_df2 = bd_tg_df2)
results_df <- data.frame(variable = character(),
                         dataframe = character(),
                         p_value = numeric(),
                         stringsAsFactors = FALSE)

for (variable_name in variables_to_test) {
  cat("Shapiro-Wilk test for", variable_name, ":\n")
  
  for (df_name in names(data_frames)) {
    result <- shapiro.test(data_frames[[df_name]][[variable_name]])
    cat("From", df_name, ":\n")
    print(result)
    cat("\n")
    results_df <- rbind(results_df, data.frame(variable = variable_name,
                                               dataframe = df_name,
                                               p_value = result$p.value))
  }
}
write.csv(results_df, 
          "bd_indiv_shapiro_wilk_results.csv", 
          row.names = FALSE)
#significant from the shapiro-wilks test
#infected: df2, tt_df2, tg_df2
#log_bd: df2, tt_df2, tg_df2
#ttx_mg: df2, tt_df2, tg_df2
#mass_g: df2
#observed_features: df2, tt_df2, tg_df2
#TTX_asvs
#antifungal_asvs:

#not significant
#mass_g: tt_df2, tg_df2





#Shannnon______________________________________________________________________________
response <- c('shannon_entropy', 'faith_pd', 
              'pielou_evenness', 'observed_features', 
              'mass_g', 'TTX_mg', 'Log_bd')
predictor <- c('TTX_Real')
data_frames_1 <- list(bd_infect_table = bd_infect_table, 
                      bd_tt_df2 = bd_tt_df2, 
                      bd_tg_df2 = bd_tg_df2)

multi_kw_results_df <- data.frame(predictor = character(),
                                  response = character(),
                                  dataframe = character(),
                                  chi_squared = numeric(),
                                  p_value = numeric(),
                                  stringsAsFactors = FALSE)

for (predictor_var in predictor) {
  for (response_var in response) {
    for (df_name in names(data_frames_1)) {
      predictor_data <- data_frames_1[[df_name]][[predictor_var]]
      response_data <- data_frames_1[[df_name]][[response_var]]
      
      multi_kw_result <- kruskal.test(response_data ~ predictor_data)
      
      bd_multi_kw_results_df <- rbind(multi_kw_results_df, data.frame(predictor = predictor_var,
                                                                   response = response_var,
                                                                   dataframe = df_name,
                                                                   chi_squared = multi_kw_result$statistic,
                                                                   p_value = multi_kw_result$p.value))
    }
  }
}

write.csv(multi_kw_results_df, "bd_multi_kruskal_wallis_results.csv", row.names = FALSE)

bd_multi_kw_results_df <- bd_multi_kw_results_df[order(multi_kw_results_df$response, 
                                                 multi_kw_results_df$predictor), ]
table2_txt <- knitr::kable(bd_multi_kw_results_df, format = "markdown", 
                           col.names = c("Response Variable", 
                                         "Predictor Variable", 
                                         "Data Frame",
                                         'Chi-Squared',
                                         "p-value"))
write.table(table2_txt, 
            "bd_multi_overall_kruskal_wallis_results.txt", 
            quote = FALSE, 
            row.names = FALSE)





#_______________________________________________________________________________
kw_2_results_df <- data.frame(predictor = character(),
                              response = numeric(),
                              chi_sqquared = numeric(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)
kw_predictor_2 <- 'Species'
kw_response_2 <- c('shannon_entropy', 'faith_pd', 
                   'pielou_evenness', 'observed_features', 
                   'mass_g', 'TTX_mg')
kw_species <- for(response_vari in kw_response_2){
  species_kw <- kruskal.test(df2[[response_vari]] ~ df2$Species, na.action = na.omit)
  print(species_kw)
  kw_2_results_df <- rbind(kw_2_results_df, data.frame(predictor = kw_predictor_2,
                                                       response = response_vari,
                                                       dataframe = 'df2',
                                                       p_value = species_kw$p.value))
}
write.csv(kw_2_results_df, 'species_kw_results.csv', row.names = FALSE)
kw_2_results_df <- kw_2_results_df[order(kw_2_results_df$response, 
                                         kw_2_results_df$predictor), ]
table_txt <- kable(kw_2_results_df, format = "markdown", 
                   col.names = c("Response Variable", 
                                 "Predictor Variable", 
                                 "Data Frame", 
                                 "p-value"))
write.table(table_txt, 
            "kruskal_wallis_results.txt", 
            quote = FALSE, 
            row.names = FALSE)

print(mean(tt_df2$TTX_mg))
mean_tg <- mean(tg_df2$TTX_mg, na.rm = TRUE)
print(mean_tg)
df2$Species <- as.factor(df2$Species)
ktest <- kruskal.test(data = df2, TTX_mg ~ Location, na.action = na.omit)
ktest

basic_plot <- df2 %>%
  ggplot(aes(Species, TTX_mg, fill = Species))+
  geom_boxplot()+
  theme_minimal()+
  guides(fill=FALSE)

basic_plot



#TTX linear models_____________________________________________________________
response_variables <- c(
  'shannon_entropy', 'faith_pd', 
  'pielou_evenness', 'observed_features', 
  'mass_g', 'TTX_mg', 'Log_bd',
  'antifungal_ASVs', 'antifungal_reads',
  'Propor_TotalAntiFungal', 'Propor_TotalTTX',
  'TTX_Richness', 'TotalTTX'
)

predictor_variables <- c(
  'shannon_entropy', 'faith_pd', 
  'pielou_evenness', 'observed_features', 
  'mass_g', 'TTX_mg', 'Log_bd',
  'antifungal_ASVs', 'antifungal_reads',
  'Propor_TotalAntiFungal', 'Propor_TotalTTX',
  'TTX_Richness', 'TotalTTX'
)
variables_indf <- response_variables %in% colnames(df)
print(variables_indf)

lm_data_frames <- list(df2 = df2, tt_df2 = tt_df2, tg_df2 = tg_df2)

model_results_df <- data.frame(response = character(),
                               predictor = character(),
                               dataframe = character(),
                               coefficients = numeric(),
                               p_value = numeric(),
                               stringsAsFactors = FALSE)

for(lm_df_names in names(lm_data_frames)){
  lm_data_frame <- lm_data_frames[[lm_df_names]]
  for(lm_response_var in response_variables){
    for(lm_predictor_var in predictor_variables){
      linear_model <- lm(lm_data_frame[[lm_response_var]] ~ lm_data_frame[[lm_predictor_var]])
      coeficient_value <- coef(summary(linear_model))[2, 1]
      p_value <- coef(summary(linear_model))[2, 4]
      model_results_df <- rbind(model_results_df, data.frame(response = lm_response_var,
                                                             predictor = lm_predictor_var,
                                                             dataframe = lm_df_names,
                                                             coefficients = coeficient_value,
                                                             p_value = p_value))
    }
  }
}
write.csv(model_results_df, "linear_regression_results.csv", row.names = FALSE)
lm_table <- kable(model_results_df, format = "markdown", 
                  col.names = c("Response Variable", "Predictor Variable", "Data Frame", "Coefficients", "p-value"))
writeLines(c(lm_table), "linear_regression_results.txt")