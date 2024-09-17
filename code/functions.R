# fonction génération des données : simul_Abundance
######## Description  : ----
#   fonction génération des données
# Paramètre :
#   B : une matrice de coefficient, de taille p x d
#   Sigma_species : une matrice de covariance entre les espèces
#   X : une matrice de variable sans l'intercept
# Valeur renvoyée : les Y
## Exemple :   ----
#  p <- 2
#  d<- 3
#  n <- 10
#  B <- matrix(1:(p*d),nrow = p,ncol = d)
#  Sigma <- diag(d)
#  X <- matrix(runif(n*p,0,1),nrow = n,ncol = p)
#  # o=rep(1,nrow(X))
#  gen_data <- simul_Abundance(B,Sigma,X)
#  #Afficher les resultats
#  print(gen_data)
# 
# # Afficher les donnees
# print(gen_data)
##### Fonction simul_Abundance : ----
simul_Abundance <-function(B,Sigma = diag(p),Covariate,o=rep(1,nrow(Covariate))){
  #### Initialisation
  #nombre d'observations de la matrice des covariates
  n <- nrow(Covariate)
  #nombre de variables explicatives (Xi)
  d <- ncol(Covariate)
  #Nombre de colonnes
  p <- ncol(B)
  
  #matrices vide qui contiendra les données de comptages
  Abundance <- matrix(NA,nrow=n,ncol=p)
  
  ##### Simuler les données
  for(i in 1:n){
    # Generation de la couche latente
    mu_Z <- o[i] + crossprod(Covariate[i,],B) # GDR - 2024-03-20 : %*% -> crossprod
    Z_i <- MASS::mvrnorm(1,mu_Z,Sigma)
    dim(Z_i)
    # Calculer les Y à partir du modèle
    for(j in 1:p)
    {Abundance[i,j] <- rpois(1,exp(Z_i[j]))}
    
  }
  
  ##### Post-traitement
  dimnames(Abundance) <- list(paste0("obs", 1:n), paste0("Y", 1:p))
  res_simu <- list(Abundance = Abundance, Covariate = Covariate, o =o)
  ##### Renvoyer le résultat
  return(res_simu)
}
# ----

# Fonction glmnet_lasso_poisson  ----
####### Destription : ----
#   Applique le modele GLMNET sur les données(Abundance, Covariate)
#   les utilisent comme entrée pour la fonction glmnet_lasso_poisson
#   affiche les coefficients estimés et les prédictions résultantes
####Parametre:
# Abundance : matrice contenant les mesures d'Abundance des espèces et qui suit la loi de poisson
# Covariate : matrice contenant les variables explicatives et qui suit la loi normale
# Valeur renvoye : Retourne une liste contenant les coefficients estimés et les prédictions

### Exemmple à faire :----
# # Génération de données d'exemple : 100 observations, 3 covariables et 2 variables d'Abundance
# set.seed(123)                            # Pour la reproductibilité
#Abundance <- matrix(rpois(100*2, 5), ncol = 2)            # Deux variables d'Abundance
#Covariate <- matrix(rnorm(100*3), ncol = 3)               # Trois covariables
# # Appel de la fonction glmnet_lasso_poisson avec les données d'exemple
#resultats <- glmnet_lasso_poisson(Abundance, Covariate)
# # Affichage des coefficients estimés
#print("Coefficients estimés :")
#print(resultats$coefficient)




##### Fonction glmnet_lasso_poisson_withoffset : ----
##### Fonction glmnet_lasso_poisson : ----
glmnet_lasso_poisson <- function(Abundance, Covariate, offset_column_name = NULL, data, verbose = TRUE) {
  
  # Initialisations
  p <- ncol(Abundance)
  n <- nrow(Abundance)
  d <- ncol(Covariate) + 1  # Nombre de covariables + intercept
  
  hat_B <- matrix(NA, nrow = d, ncol = p)  # Matrice pour stocker les coefficients
  colnames(hat_B) <- colnames(Abundance)
  colnames <- c("(Intercept)", colnames(Covariate))
  
  
  prediction <- matrix(NA, ncol = p, nrow = n)
  colnames(prediction) <- colnames(Abundance)
  
  # Gestion de l'offset
  if (!is.null(offset_column_name)) {
    offset_values <- data[[offset_column_name]]
    data <- data[, !(names(data) %in% offset_column_name)]
  }
  
  # Boucle de modélisation
  for (j in 1:p) {
    if (is.null(offset_column_name)) {
      res_glmnet <- glmnet(x = Covariate, y = Abundance[, j], family = "poisson")
      if (verbose) {
        cat("Coefficients pour la variable de réponse", j, ":\n")
        print(coef(res_glmnet))
      }
      
      hat_B[, j] <- as.vector(coef(res_glmnet, s = 1)) 
      suppressWarnings(prediction[, j] <- predict(res_glmnet, newx = as.matrix(Covariate), s = 1, newoffset = ))
      
    } else {
      res_glmnet <- glmnet(x = Covariate, y = Abundance[, j], offset = offset_values, family = "poisson")
      if (verbose) {
        cat("Coefficients pour la variable de réponse", j, ":\n")
        print(coef(res_glmnet))
      }
      
      hat_B[, j] <- as.vector(coef(res_glmnet, s = 1)) 
      suppressWarnings(prediction[, j] <- predict(res_glmnet, newx = as.matrix(Covariate), s = 1, newoffset = offset_values))
    }
    
    
  }
  
  # Retourner une liste de hat_B et prediction
  return(list(hat_B = hat_B, fitted = prediction))
}

# ----

## fonction  calcul des parametres : param_init ----
#   Description : ----
#   param_init initialise les parametres de regressions estimes par lm.fit (Beta, s, M)
# Paramètre :
#   covariates : une matrice des variables  explicatives, de taille n x d
#   Abundance : une matrice Abondace(ncol = Nombre d'especes, nrow = Nombre de site)
# Valeur renvoyée : les coefficients de regressions

#### Exemple a faire : ----
# set.seed(123)                       # Pour la reproductibilité
# # Génération de données d'exemple : 100 observations, 3 covariables et 2 variables d'Abundance
# Abundance <- matrix(rpois(100*2, 5), ncol = 2)             # Deux variables d'Abundance
# Covariate <- matrix(rnorm(100*3), ncol = 3)                # Trois covariables
# 
# # Appel de la fonction param_init pour initialiser les paramètres du modèle
# params_init <- param_init(Covariate = Covariate, Abundance = Abundance)
# 
# # Affichage des coefficients estimés
# print("Coefficients estimés :")
# print(params_init$B)
# # Affichage des résidus
# print("Résidus :")
# print(params_init$M)L'objectif de cette étude était de vérifier si les espèces d'araignées chasseuses se trouvaient toujours dans les mêmes endroits, en comparant avec des études précédentes. En outre, elle visait à décrire où chaque espèce d'araignée était la plus abondante et dans quelles conditions environnementales elle pouvait survivre.


### Fonction : ---- 
param_init <- function(Covariate,Abundance,Offset = matrix(0, ncol = ncol(Abundance),
                                                           nrow = nrow(Abundance))){
  
  #### Initialisation
  #Nombre de site
  n <- nrow(Abundance)
  #Nombre d'especes
  p <- ncol(Abundance)
  
  # Matrice des Covariate + l'intercept ou ajout de l'intercept dans les données
  Covariate <- cbind(rep(1,n),as.matrix(Covariate))
  
  
  # Estimation du modèle de régression sur les Covariates
  res_lm.fit <- lm.fit(Covariate, log((1 + as.matrix(Abundance))/(1 + exp(as.matrix(Offset)))))
  
  # Matrice des coefficients de regression estimés
  B <- res_lm.fit$coefficients
  # Matrice des residus
  M <- matrix(res_lm.fit$residuals, n, p)
  # Matrice d'observation
  S <- matrix(1, n, p)
  
  
  # renvoie une liste de coefficients de regression estimés, des residus et d'observation 
  return(list(B=B,M=M,S=S))
}

# ----

#  Fonction pour formater les données pour le modèle : format_PLN ----
# Description :Fonction donner un format les données pour le modèle
## parametres :
#   covariates : une matrice des variables  explicatives, de taille n x d
#   Abundance : une matrice Abondace(ncol = Nombre d'especes, nrow = Nombre de site)
## Valeur renvoyée : une liste contenannt les donnees dans data, les parametres dans params

# Exemple a faire : ----
# set.seed(123)                       # Pour la reproductibilité
# # Génération de données d'exemple : 100 observations, 3 covariables et 2 variables d'Abundance
# Abundance <- matrix(rpois(100*2, 5), ncol = 2)             # Deux variables d'Abundance
# Covariate <- matrix(rnorm(100*3), ncol = 3)                # Trois covariables
# 
# data_formatted <- format_PLN(Abundance = Abundance, Covariate = Covariate)
# print("Données formatées pour le modèle :")
# print(data_formatted)
# Fonction format_PLN : ----
format_PLN <- function(Abundance, Covariate) {
  
  # Calcul du poids w
  n <- nrow(Abundance)
  w <- rep(1, n)  # poids égal à 1 pour chaque observation
  
  #### Appel des fonctions simul_Y et param_init pour obtenir les paramètres
  param_data <- param_init(Covariate, Abundance)
  
  # Les paramètres du modèle
  params <- list(B = param_data$B , Sigma = Sigma, M = param_data$M , S = param_data$S)
  
  
  # Données pour le modèle : Y, X, offset O, poids w
  data <- list(Abundance = Abundance, Covariate = Covariate, O = param_data$o, w = w)
  
  # Renvoie une liste contenant les données de comptage (Y, X, O, w) et les paramètres
  return(list(data = data, params = params))
}
# ----


##### Fonction glmnet_lasso_poisson_withoffset : ----
glmnet_lasso_poisson_withoffset <- function(Abundance, Covariate, offset_var, verbose=T){
  #### Initialisation
  #Nombre de colonnes dans la matrice Abundance, qui représente le nombre de variables à prédire
  p <- ncol(Abundance)
  
  #Nombre de lignes dans la matrice Abundance, qui représente le nombre d'observations
  n <- nrow(Abundance)
  
  # Nombre de colonnes dans la matrice Covariate, qui représente le nombre de Covariates plus un pour l'intercept
  d <- ncol(Covariate)+1
  
  # Initialisation d'une matrice vide pour stocker les coefficients estimés pour chaque variable réponse
  hat_B <- matrix(NA, ncol=p, nrow=d)
  
  # Initialisation d'une matrice vide pour stocker les prédictions pour chaque variable réponse
  prediction <- matrix(NA, ncol = p, nrow = n)
  
  #### Ajuste le modèle
  # Boucle sur chaque variable réponse
  for (j in 1:p) {
    
    # Estimation du modèle de régression LASSO avec glmnet pour la variable réponse j
    res_Pois_off <- glmnet(x = Covariate, y = Abundance[, j], offset = offset_var, family = ("poisson"))
    if(verbose)
      print(coef(res_Pois_off))
    
    # Stockage des coefficients estimés pour la variable réponse j
    hat_B[, j] <- as.matrix(coef(res_Pois_off, s = 1))
    
    # Prédiction des valeurs pour la variable réponse j à partir du modèle ajusté
    prediction[, j] <- predict(res_Pois_off, newx = Covariate, s = 1, newo = offset_var)
    
  }
  
  #### Résultat
  # Retourne une liste contenant les coefficients estimés et les prédictions
  return(list( hat_B = hat_B,prediction = prediction))
}
#----

#Description : ----
# Fonction pour calculer les facteurs d'inflation de la variance (VIF)
# parametres : Covariate
#renvoie : les VIF sur chaque variable

# Exemple a executer : ----
# Création d'un ensemble de données fictif
#set.seed(123)  # Pour la reproductibilité
#vegetation_index <- rnorm(100, mean = 50, sd = 10)
#wetness <- rnorm(100, mean = 12, sd = 3)
#center_xv <- rnorm(100, mean = 40, sd = 5)

# Créer un dataframe avec ces variables
#data <- data.frame(vegetation_index, wetness, center_xv)
# Utiliser la fonction calculate_VIF pour calculer les VIF
#VIF_results <- calculate_VIF(data)
# Afficher les résultats
#print(VIF_results)
