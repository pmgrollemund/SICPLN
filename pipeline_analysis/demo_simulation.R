################################################################################
######### Démo des fonction de moi
#########    (données simulées)
######### Guy Darcy Remesha et Paul-Marie Grollemund
######### 2024-03-25
################################################################################

# Nettoyage de l'environnement de travail
rm(list = ls())

# Récupération des arguments passés au script
args <- commandArgs(trailingOnly = TRUE)

# Vérification du nombre d'arguments passés
if (length(args) < 2) {
  stop("Veuillez spécifier le répertoire des fonctions personnalisées ainsi que le repertoire de sorie.")
}

# général options
code_path <- args[1]
output <- args[2]

# Créer le répertoire de sortie
temps <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
output_dir <- file.path(output, paste("SICPLN_execution", temps, sep = "_"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Chargement du fichier contenant les fonctions personnalisées
source(code_path)

# Definition des  parametres pour la génération des données d'exemple
set.seed(123)  # Pour la reproductibilité
n <- 10  # Nombre d'observations
p <- 2  # Nombre de covariables
d <- 3  # Nombre d'espèces

# Définition des paramètres pour le calcul des erreurs 
nbre_repetition <- 5
beta <- matrix(1:(p*d),nrow = d,ncol = p)
Sigma <- diag(p)

# Génération des covariables (matrice aléatoire)
Covariate <- matrix(runif(n * d), nrow = n, ncol = d)

# Génération des données d'exemple
gen_data <- simul_Abundance(B = beta, Sigma, Covariate)

rownames(gen_data$Covariate) <- rownames(gen_data$Abundance)

gen_data$Covariate <- cbind(gen_data$Covariate, o = gen_data$o)


data_tmp <- list(Abundance = gen_data$Abundance, Covariate = gen_data$Covariate)

data_to_PLN <- prepare_data(data_tmp$Abundance, data_tmp$Covariate)
data_to_PLN <- data_to_PLN[, -ncol(data_to_PLN)]


formule <- Abundance ~ .
data <- data_to_PLN
# Spécifiez les paramètres de contrôle
control <- PLN_param(backend = "nlopt", covariance = "full")


# Application de la fonction SICPLN avec les données générées
res_SICPLN <- SICPLN(formule , offset_column_name = "o", data = data_to_PLN, control = control)
print(res_SICPLN$res_fisher$model_par$B)



# Application de la fonction glmnet_lasso_poisson avec les données générées
res_glmpois <- glmnet_lasso_poisson(data_tmp$Abundance,data_tmp$Covariate[, -ncol(data_tmp$Covariate)], verbose = FALSE)
res_glmpois$hat_B
res_glmpois$prediction

# Application de la fonction error_cal avec les paramètres spécifiés
resultats_error <- error_cal(nombre_espece = p, taille_echant = n, nbre_repetition =10)
resultats_error$df_bp
# Calcul de l'erreur MSE entre les vraies valeurs et les valeurs prédites pour SICPLN
mse_sicpln <- mse_predict(gen_data$Abundance, res_SICPLN$res_pln$fitted)


# Calcul de l'erreur MSE entre les vraies valeurs et les valeurs prédites pour GLMNET
mse_glmnet <- mse_predict(gen_data$Abundance, res_glmpois$prediction)

# Calcul de l'erreur MSE entre les vraies valeurs et les valeurs prédites pour PLN
mse_pln <- mse_predict(gen_data$Abundance, res_SICPLN$fitted)



# Affichage des résultats
print("Résultats SICPLN :")
coef(res_SICPLN)
print("Résultats GLMNET :")
print(res_glmpois)

print("Erreurs de prédiction (MSE) :")
print("=============================")
print(paste("MSE SICPLN:", mse_sicpln))
print(paste("MSE GLMNET:", mse_glmnet))
print(paste("MSE PLN:", mse_pln))

str(resultats_error)
#=====================================================================================================================
# Definition des  parametres pour la génération des données d'exemple
set.seed(123)  # Pour la reproductibilité
n <- 30  # Nombre d'observations
p <- 4  # Nombre d'espèces
d <- 3  # Nombre de covariables

# Définition des paramètres pour le calcul des erreurs 
nbre_repetition <- 10

Error_n_30 <- error_cal(nombre_espece = p,taille_echant = n, nbre_repetition)

#====================================================================================
# Definition des  parametres pour la génération des données d'exemple
set.seed(123)  # Pour la reproductibilité
n <- 50  # Nombre d'observations
p <- 4  # Nombre d'espèces
d <- 3  # Nombre de covariables

# Définition des paramètres pour le calcul des erreurs 
nbre_repetition <- 100

Error_n_50 <- error_cal(nombre_espece = p,taille_echant = n,nbre_repetition)

#=====================================================================================
# Definition des  parametres pour la génération des données d'exemple

n <- 100  # Nombre d'observations
p <- 4  # Nombre d'espèces
d <- 3  # Nombre de covariables

# Définition des paramètres pour le calcul des erreurs 
nbre_repetition <- 100

Error_n_100 <- error_cal(nombre_espece = p,taille_echant = n ,nbre_repetition)

#=========================================================================================
# Definition des  parametres pour la génération des données d'exemple

n <- 1000  # Nombre d'observations
p <- 4  # Nombre d'espèces
d <- 3  # Nombre de covariables

# Définition des paramètres pour le calcul des erreurs 
nbre_repetition <- 100

Error_n_1000 <- error_cal(nombre_espece,taille_echant,nbre_repetition)

#==================================================================================================
##  GRAPHIQUE POUR LES ERREURs GENERES EN TENANT COMPTE DU CHANGEMENT DE TAILLE D"ECHANTILLON

# Graphique des coefficients des matrices B
boxplt_echantillon <- boxplot_error_norm(Error_n_30, Error_n_50, Error_n_100, Error_n_1000)
ggsave("E:/Master SAAD/Stage Pro/new_code/plot_simulation/boxplot_error_coef_ech.pdf", boxplt_echantillon, width = 10, height = 6, units = "in")


# erreur de prediction en fonction de la taille d'échantillon
plt_echantillon <- plot_error_prediction(Error_n_30, Error_n_50, Error_n_100, Error_n_1000)
ggsave("E:/Master SAAD/Stage Pro/new_code/plot_simulation/plot_error_pred.pdf", plt_echantillon, width = 10, height = 6, units = "in")
