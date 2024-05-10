
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
#___________________________________________________________
tt_jac_import <- qiime2R::read_qza('tato_new/core_metrics/jaccard/jaccard_distance_matrix.qza')
tt_jac_mat1 <- as(tt_jac_import$data, 'matrix')
tt_jac_mat1 <- as.data.frame(tt_jac_mat1)
tt_jac_ids <- rownames(tt_jac_mat1)
tt_jac_new_order <- order(rownames(tt_jac_mat1))
tt_jac_ordered <- tt_jac_mat1[tt_jac_new_order,tt_jac_new_order]
tt_jaccard_matrix <- as.matrix(tt_jac_ordered)
#___________________________________________________________
tt_uwu_import <- qiime2R::read_qza('tato_new/core_metrics/unweighted/unweighted_unifrac_distance_matrix.qza')
tt_uwu_mat1 <- as(tt_uwu_import$data, 'matrix')
tt_uwu_mat1 <- as.data.frame(tt_uwu_mat1)
tt_uwu_ids <- rownames(tt_uwu_mat1)
tt_uwu_new_order <- order(rownames(tt_uwu_mat1))
tt_uwu_ordered <- tt_uwu_mat1[tt_uwu_new_order,tt_uwu_new_order]
tt_unweighted_unifrac_matrix <- as.matrix(tt_uwu_ordered)
#_Bray-Curtis___________________________________________________________________________________________--
tt_bc_import <- qiime2R::read_qza('tato_new/core_metrics/bray_curtis/bray_curtis_distance_matrix.qza')
tt_bc_mat1 <- as(tt_bc_import$data, 'matrix')
tt_bc_mat1 <- as.data.frame(tt_bc_mat1)
tt_bc_ids <- rownames(tt_bc_mat1)
tt_bc_new_order <- order(rownames(tt_bc_mat1))
tt_bc_ordered <- tt_bc_mat1[tt_bc_new_order,tt_bc_new_order]
tt_bray_curtis_matrix <- as.matrix(tt_bc_ordered)

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
