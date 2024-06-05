################################################################################
######### Application sur les données GENUS2
#########    (données HSPIDER)
######### Guy Darcy Remesha et Paul-Marie Grollemund
######### 2024-04-09
################################################################################

# Nettoyage de l'environnement de travail
rm(list = ls())

# Récupération des arguments passés au script
args <- commandArgs(trailingOnly = TRUE)

# Vérification du nombre d'arguments passés
if (length(args) < 2) {
  stop("Veuillez spécifier le chemin de la base de données, le répertoire de sortie en arguments et le répertoire des fonctions personnalisées.")
}

# général options
output <- args[1]
code_path <- args

output <- "result"
code_path <- "../code/SICPLN_Pack.R"


# Créer le répertoire de sortie
temps <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path(output, paste("SICPLN_execution", temps, sep = "_"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Importer les fonctions
source(code_path)

#source("/home/gdremesh/Bureau/new_code/demo_spider/SICPLN_Pack.R")
suppressPackageStartupMessages(library(VGAM))
spider <- hspider

cat("\n\n\t Pretreament .... :\n")
# Séparation des données d'abondance et des covariables
Abundance <- spider[,-c(1:6)]
Covariate <- spider[,c(1:6)]

# enleve de ReflLux car fortement correllée
prep_spider <-  PLNmodels::prepare_data(Abundance,Covariate)

# calcul des VIF avec variable Reflux
Res_vif_A_Ref <- calculate_VIF(Covariate[,-3], 5)
cat("\n\n\t Variable selection using VIF :\n")
print(Res_vif_A_Ref$selected_variables)

# Variables retenues après étude des VIF
#retain_Covar <- Covariate[, Res_vif_A_Ref$selected_variables]
retain_Covarr <- prep_spider[, Res_vif_A_Ref$selected_variables]

rownames(retain_Covarr) <- rownames(prep_spider$Abundance)


data_tmp <- list(Abundance = prep_spider$Abundance, Covariate = retain_Covarr)
data_to_PLN <- prepare_data(data_tmp$Abundance, data_tmp$Covariate)
data_to_PLN <- data_to_PLN[, -ncol(data_to_PLN)]



# Spécifiez les paramètres de contrôle
control <- PLN_param(backend = "nlopt", covariance = "full")

###### Cas avec PLN (spider) : ----
# Ajustement du modèle PLN (Poisson-Lognormal) avec la surface en tant qu'offset
cat("\n\n\t launch of the PLN model .\n")
cat("        ==========================\n")

spider_pln <- PLNmodels::PLN(Abundance ~ . , control = control, data = data_to_PLN) 
print(spider_pln$model_par$B)


# Ajustement pour la selection de variable
cat("\n\n\t launch of the SIC-PLN Model. \n")
cat("        ============================ \n")

formule <- Abundance ~ .
data <- data_to_PLN

spider_sicpln <- SICPLN(formule , data = data_to_PLN, control = control)
spider_sicpln$res_fisher$model_par$B


# Convertir la matrice en un format adapté à ggplot2
data_df <- melt(data$Abundance)
data_sic <- melt(spider_sicpln$fitted)


# Renommer les colonnes
colnames(data_df) <- c("Site", "Species", "Comptage")
colnames(data_sic) <- c("Site", "Species", "Comptage")

# Graphiques
plot_coef_pln <- coef_genus_barre(spider_pln$model_par$B[-1,], 
                                  axis_names = c("Spider", "Variables", "Coefficients value", "PLN"), data)
plot_coef_sicpln <- coef_genus_barre(spider_sicpln$res_fisher$model_par$B[-1,],
                                     axis_names = c("Spider species", "Variables", "Coefficients value", "SICPLN"), data)



plt_abund <- boxplot_graph(data_df = data_df, x = "Species", y = "Comptage", title = "Boxplot of Abundances by Species")
plt_abund_sic <- boxplot_graph(data_df = data_sic, x = "Species", y = "Comptage", title = "Boxplot of Abundances fitted by Species")

# Utilisation de la fonction avec les données d'abondance et les valeurs prédites
resid_pln <- plot_residus_density(data$Abundance, spider_pln$fitted, output_dir) #, output_dir
resid_sicpln <- plot_residus_density(data$Abundance, spider_sicpln$fitted, output_dir) #, output_dir

# # density des y et densité des y prédit
density_abundance_pln <- plot_density(data$Abundance, spider_pln$fitted, main = "Density Plot of Abundance from PLN", xlab = "Abundance")
density_abundance_sicpln <- plot_density(data$Abundance, spider_sicpln$fitted, main = "Density Plot of Abundance from SICPLN", xlab = "Abundance")

# image de la matrice sigma (precision matrix)
tmp <- solve(spider_sicpln$res_fisher$model_par$Sigma)
norm <- sqrt(diag(tmp)) %*% t(sqrt(diag(tmp)))
cor_mat <- tmp / norm
precision_mat_sicpln <- plot_precision_matrix(cor_mat, axis_names = c("", "", "", "SICPLN_Spider precision matrix of Sigma"), limit = NULL)

# Appel de la fonction plot_mse
#label_order <- paste0("Species",1:ncol(data$Abundance),sep = "")
MSE_pln_graph <- plot_mse_barplot(data$Abundance, spider_pln$fitted)
MSE_sicpln_graph <- plot_mse_barplot(data$Abundance, spider_sicpln$fitted)

# # Relation entre variable environnementale et Espèce
plotcovar_Abund <- plot_abundance_vs_environment(spider_sicpln$res_fisher$model_par$B, data, output_dir) #, output_dir

# Sauvegarde des graphiques dans le dossier de sortie spécifié en argument
ggsave(file.path(output_dir, "plot_pln.pdf"), plot_coef_pln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plot_sicpln.pdf"), plot_coef_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plt_abund.pdf"), plt_abund, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "plt_abund_sic.pdf"), plt_abund_sic, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "density_abundance_pln.pdf"), density_abundance_pln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "density_abundance_sicpln.pdf"), density_abundance_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "precision_mat_sicpln.pdf"), precision_mat_sicpln, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "MSE_pln_graph.pdf"), MSE_pln_graph, width = 10, height = 6, units = "in")
ggsave(file.path(output_dir, "MSE_sicpln_graph.pdf"), MSE_sicpln_graph, width = 10, height = 6, units = "in")

save.image(file.path(output, "Renv.RData"))
