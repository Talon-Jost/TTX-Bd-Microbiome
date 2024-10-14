rm(list=ls())


library(tidyverse)
library(readr)
library(ggplot2)
library(car)
library(MuMIn)
library(vegan)
library(effects)
library(qiime2R)
set.seed(12345)

location_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan',
                     'salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')
tg_colors <- c('red', 'blue', 'green', 'orange', '#8DD3C7', 'cyan')
tt_colors <- c('salmon', 'violet', '#009933', '#ffcc00', '#330099', '#660033')


meta.data <- read.csv('newt_metadata.tsv', sep ='\t') %>% 
  rename(SampleID = 'X.SampleID') %>% 
  filter(!SampleID %in% c('Pos',
                          'Neg',
                          'FCSP-05'))

#Taxonomy ####
taxonomy <- as.data.frame(read_qza('classification.qza')$data)
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
differentails <- read.table('exported_difs/p_val_slice.csv', header = T, sep = ',') %>% 
  rename(intercept = 'X.Intercept.') %>% 
  rename(p_value = 'TTX_RealLow')

qval <- read.table('exported_difs/q_val_slice.csv', header = T, sep = ',')  %>% 
  rename(q_value = 'TTX_RealLow')

matched <- differentails %>% 
  left_join(qval, by = 'id') %>% 
  filter(q_value < 0.05)

lfc <- read.csv('exported_difs/lfc_slice.csv', header = T, sep = ',') %>% 
  rename(lfc = 'TTX_RealLow')

se <- read.csv('exported_difs/se_slice.csv', header = T, sep = ',') %>% 
  rename(std.err = 'TTX_RealLow')

full_match <- matched %>% 
  left_join(taxa, by = 'id') %>%
  left_join(lfc, by = 'id') %>% 
  left_join(se, by = 'id') %>% 
  arrange(desc(lfc)) %>% 
  mutate(asv_id = paste0("ASV", sprintf("%03d", seq(1, nrow(matched))))) %>% 
  mutate(label = paste0(asv_id, " (", Family, ")"))

asvs_to_combine <- c("ASV001", "ASV003")
asvs_to_combine2 <- c("ASV002", "ASV005")

combined_data <- full_match %>%
  mutate(Family = ifelse(asv_id %in% asvs_to_combine, "Pseudomonadaceae 1", Family)) %>% 
  mutate(Family = ifelse(asv_id %in% asvs_to_combine2, "Pseudomonadaceae 2", Family)) %>% 
  group_by(Family) %>%
  summarise(lfc = sum(lfc), .groups = 'drop')

combined_data <- combined_data %>%
  mutate(fill_color = case_when(
    Family == "Pseudomonadaceae 1" ~ "red",  
    Family == "Pseudomonadaceae 2" ~ "red",
    Family == 'Sphingomonadaceae' ~ 'red',
    TRUE ~ ifelse(lfc > 0, "#61D04F", "#28E2E5")
  )) %>% 
  mutate(Family = ifelse(Family == 'uncultured', 'Uncultured Verrucomicrobiae', Family))

#you ended wanting to join all the datasets. all the folders in exported difs and differentails are the same. just go join them because you need the beta and standard error to make the plot

library(forcats)
da.fig <-ggplot(combined_data, aes(fct_reorder(Family, lfc), 
                          y = lfc, 
                          fill = fill_color)) +
  geom_bar(stat = 'identity', width = 0.7) +
  coord_flip() +
  scale_fill_identity() +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold',
                                  hjust = 0.5),
        axis.title.x = element_text(face = 'bold'), 
        axis.title.y = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        axis.line.x = element_line(), 
        axis.line.y = element_line(),
        plot.margin = margin(t = 8, r = 8, b = 8, l = 8),  
        legend.position = 'none') +
  labs(title = 'Differential Abundance of Bacterial Families',
       y = 'Log Fold Change (LFC)',
       x = '')
da.fig

ggsave('differential_abundance_figure.pdf', da.fig, dpi = 1000, unit = 'in', height = 10, width = 8)
ggsave('differential_abundance_figure.jpg', da.fig, dpi = 1000, unit = 'in', height = 10, width = 8)

