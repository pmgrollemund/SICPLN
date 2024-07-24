################################################################################
######### Application sur les données GENUS2
#########    (données GENUS2)
######### Guy Darcy Remesha et Paul-Marie Grollemund
######### 2024-04-09
################################################################################

# Nettoyage de l'environnement de travail
rm(list = ls())

# Récupération des arguments passés au script
args <- commandArgs(trailingOnly = TRUE)

# Vérification du nombre d'arguments passés
if (length(args) < 3) {
  stop("Veuillez spécifier le chemin de la base de données, le répertoire de sortie en arguments et le répertoire des fonctions personnalisées.")
}

# général options
data_path <- args[1]
output <- args[2]
code_path <- args[3]

# Créer le répertoire de sortie
temps <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path(output, paste("SICPLN_execution", temps, sep = "_"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Importer les fonctions
source(code_path)

# Chargement de la base de données à partir du premier argument
load(data_path)
data <- genus2

# Séparation des données d'abondance et des covariables
Abundance <- genus2[, 1:15]
Covariate <- genus2[, -c(1:15)]

# Génération des variables environnementales catégorielles
#var_cat <- as.factor(matrix(sample(letters[1:4], nrow(Covariate), replace = TRUE), ncol = 1))
# cat1 = as.factor(matrix(sample(letters[1:3], nrow(Covariate), replace = TRUE), ncol = 1))
# cat2 = as.factor(matrix(sample(LETTERS[1:2], nrow(Covariate), replace = TRUE), ncol = 1))
# cat3 = as.factor(ifelse(cat2 == "A", "Oui", "Non"))
# 
# Covariate <- cbind(Covariate, var_cat, cat1, cat2, cat3)

Covariate$rainfall <- rowSums(Covariate[, 27:38])
Covariate$vegetation_index <- rowSums(Covariate[, 1:23])

# Préparation des données pour l'analyse
data_to_PLN <- PLNmodels::prepare_data(Abundance, Covariate)

prep_Covar_environement <- data_to_PLN[ , c("vegetation_index","rainfall",
                                            "altitude","awd","mcwd", "wetness",
                                            "CWD","center_x", "center_y",
                                            "pluvio_an", "surface")]  # , "var_cat", "cat1", "cat2", "cat3"

# Calcul des VIF (variance inflation factors)
Res_vif <- calculate_VIF(prep_Covar_environement, 5)


# Variables retenues après étude des VIF
retain_Covar <- prep_Covar_environement[, Res_vif$selected_variables]
rownames(retain_Covar) <- rownames(data_to_PLN$Abundance)

data_tmp <- list(Abundance = data_to_PLN$Abundance, Covariate = retain_Covar)
data_to_PLN <- prepare_data(data_tmp$Abundance, data_tmp$Covariate)
data_to_PLN <- data_to_PLN[, -ncol(data_to_PLN)]



# Appeler la fonction pour générer les graphiques
data_reg <- data.frame(
  Abundance = data_to_PLN$Abundance,
  vegetation_index = data_to_PLN$vegetation_index,
  rainfall = data_to_PLN$rainfall,
  wetness = data_to_PLN$wetness,
  center_x = data_to_PLN$center_x,
  center_y = data_to_PLN$center_y,
  surface = data_to_PLN$surface
)
row.names(data_reg) <- as.character(1:990)

regression_species <- plot_abundance_regressions(data = data_reg, file_path = output_dir)


# Spécifiez les paramètres de contrôle
control <- PLN_param(backend = "nlopt", covariance = "full")

###### Cas avec PLN : ----
cat("\n\n\tlancement du modèle PLN.\n")
cat("        ==========================\n")

offset_variable <- "surface"
offset_values <- data_to_PLN[[offset_variable]]
genus_to_PLN <- data_to_PLN[, colnames(data_to_PLN) != offset_variable]

genus2_pln_offset <- PLNmodels::PLN(Abundance ~ . + offset(offset_values), control = control, data = genus_to_PLN)
genus2_pln_offset$model_par$B

cat("\n\n\tlancement du modèle SICPLN.\n")
cat("        ============================\n")

formule <- Abundance ~ .
data <- data_to_PLN

genus2_sicpln <- SICPLN(formule, offset_column_name = "surface", data = data_to_PLN, numb_eps = 100, maxIter = 300, control = control)
genus2_sicpln$res_fisher$model_par$B


# Extract abundance and covariate matrices
Abundance_matrix <- data$Abundance
Covariate_matrix <- data[, colnames(retain_Covar)]
Covariate_matrix <- Covariate_matrix[, colnames(Covariate_matrix) != offset_variable]
# Check the structure of the matrices and vectors
str(Abundance_matrix)
str(Covariate_matrix)

genus_glmnet <- glmnet_lasso_poisson(Abundance = data$Abundance,
                                     Covariate = Covariate_matrix,
                                     offset_column_name = "surface",
                                     data = data,
                                     verbose = TRUE)

genus_glmnet$hat_B
genus_glmnet$prediction


# Convertir la matrice en un format adapté à ggplot2
data_df <- melt(data$Abundance)
data_sic <- melt(genus2_sicpln$fitted)
data_glmnet <- melt(genus_glmnet)

# Renommer les colonnes
colnames(data_df) <- c("Site", "Species", "Comptage")
colnames(data_sic) <- c("Site", "Species", "Comptage")
colnames(data_glmnet) <- c("Site", "Species", "Comptage")
# Graphiques
plot_coef_pln <- coef_genus_barre(genus2_pln_offset$model_par$B[-1,], 
                                  axis_names = c("Genus", "Variables", "Coefficients value", "PLN"), data = data)
plot_coef_sicpln <- coef_genus_barre(genus2_sicpln$res_fisher$model_par$B[-1,],
                                     axis_names = c("Genus species", "Variables", "Coefficients value", "SICPLN"), data = data)
plot_coef_glmnet <- coef_genus_barre(genus_glmnet$hat_B[-1,],
                                     axis_names = c("Genus species", "Variables", "Coefficients value", "GLMNET"), data = data)


# plot_coef_pln + geom_hline(yintercept = 3.5) + geom_hline(yintercept = 5.5)

plt_abund <- boxplot_graph(data_df = data_df, x = "Species", y = "Comptage", title = "Boxplot of Abundances by Species")
plt_abund_sic <- boxplot_graph(data_df = data_sic, x = "Species", y = "Comptage", title = "Boxplot of Abundances fitted by Species")
plt_abund_glmnet <- boxplot_graph(data_df = data_glmnet, x = "Species", y = "Comptage", title = "Boxplot of Abundances fitted by Species")

# Utilisation de la fonction avec les données d'abondance et les valeurs prédites
resid_pln <- plot_residus_density(data$Abundance, genus2_pln_offset$fitted, output_dir) #, output_dir
resid_sicpln <- plot_residus_density(data$Abundance, genus2_sicpln$fitted, output_dir) #, output_dir
resid_glmnet <- plot_residus_density(data$Abundance, genus_glmnet$fitted, output_dir) #, output_dir



# # density des y et densité des y prédit
density_abundance_pln <- plot_density(data$Abundance, genus2_pln_offset$fitted, main = "Density Plot of Abundance from PLN", xlab = "Abundance")
density_abundance_sicpln <- plot_density(data$Abundance, genus2_sicpln$fitted, main = "Density Plot of Abundance from SICPLN", xlab = "Abundance")
density_abundance_glmnet <- plot_density(data$Abundance, genus_glmnet$fitted, main = "Density Plot of Abundance from GLMNET", xlab = "Abundance")


# image de la matrice sigma (precision matrix)
tmp <- solve(genus2_sicpln$res_fisher$model_par$Sigma)
norm <- sqrt(diag(tmp)) %*% t(sqrt(diag(tmp)))
cor_mat <- tmp / norm
precision_mat_sicpln <- plot_precision_matrix(cor_mat, axis_names = c("", "", "", "SICPLN_Genus precision matrix of Sigma"), limit = NULL)

# Appel de la fonction plot_mse
label_order <- paste0("gen",1:15,sep = "")
MSE_pln_graph <- plot_mse_barplot(data$Abundance, genus2_pln_offset$fitted, label_order = label_order)
MSE_sicpln_graph <- plot_mse_barplot(data$Abundance, genus2_sicpln$fitted, label_order = label_order)
MSE_glmnet_graph <- plot_mse_barplot(data$Abundance, genus_glmnet$fitted, label_order = label_order)

# # Relation entre variable environnementale et Espèce
plotcovar_Abund_sic <- plot_abundance_vs_environment(genus2_sicpln$res_fisher$model_par$B, data, output_dir) #, output_dir
#plotcovar_Abund_glmnet <- plot_abundance_vs_environment(genus_glmnet$hat_B, data, output_dir) #, output_dir

combined_density_plot <- plot_combined_density(data$Abundance, genus2_pln_offset$fitted,genus2_sicpln$fitted, genus_glmnet$fitted,
                                               main = "Combined Density Plot of Abundance", xlab = "Abundance")

# Sauvegarde des graphiques dans le dossier de sortie spécifié en argument
ggsave(file.path(output_dir, "plot_pln.pdf"), plot_coef_pln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plot_sicpln.pdf"), plot_coef_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plot_glmnet.pdf"), plot_coef_glmnet, width = 10, height = 6, units = "in")


ggsave(file.path(output_dir, "plt_abund.pdf"), plt_abund, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plt_abund_sic.pdf"), plt_abund_sic, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plt_abund_glmnet.pdf"), plt_abund_glmnet, width = 10, height = 6, units = "in")

ggsave(file.path(output_dir, "density_abundance_pln.pdf"), density_abundance_pln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "density_abundance_sicpln.pdf"), density_abundance_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "density_abundance_glmnet.pdf"), density_abundance_glmnet, width = 10, height = 6, units = "in")

ggsave(file.path(output_dir, "precision_mat_sicpln.pdf"), precision_mat_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "MSE_pln_graph.pdf"), MSE_pln_graph, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "MSE_sicpln_graph.pdf"), MSE_sicpln_graph, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "MSE_glmnet_graph.pdf"), MSE_glmnet_graph, width = 10, height = 6, units = "in")

ggsave(file.path(output_dir, "combined_density_plot.pdf"), combined_density_plot, width = 10, height = 6, units = "in")

# ggsave(file.path(output_dir, "resid_pln.pdf"), resid_pln, width = 10, height = 6, units = "in")
# ggsave(file.path(output_dir, "resid_sicpln.pdf"), resid_sicpln, width = 10, height = 6, units = "in")

save.image(file.path(output, "Renv.RData"))
