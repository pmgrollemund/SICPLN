######### Guy Darcy Remeasha et Paul-Marie Grollemund
######### 2024-03-25
################################################################################

################################################################################ ----
#######
# Required packages 
#######
################################################################################ ----

suppressPackageStartupMessages(library("torch"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("PLNmodels"))
suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("MASS"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggpattern"))
suppressPackageStartupMessages(library("svMisc"))
suppressPackageStartupMessages(library("rlang"))
################################################################################ ----
#######
# Required objects 
#######
################################################################################ ----




################################################################################ ----
#######
# Functions
#######
################################################################################

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


##### Fonction glmnet_lasso_poisson : ----
glmnet_lasso_poisson <- function(Abundance, Covariate,verbose=T){
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
    res_glmnet <- glmnet(x = Covariate, y = Abundance[, j], family = "poisson")
    if(verbose)
      print(coef(res_glmnet))
    
    # Stockage des coefficients estimés pour la variable réponse j
    hat_B[, j] <- as.matrix(coef(res_glmnet, s = 1))
    
    # Prédiction des valeurs pour la variable réponse j à partir du modèle ajusté
    prediction[, j] <- predict(res_glmnet, newx = Covariate, s = 1)
    
  }
  
  #### Résultat
  # Retourne une liste contenant les coefficients estimés et les prédictions
  return(list( hat_B = hat_B,prediction = prediction))
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
    # print(params_init$M)

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

  
## Fonction qui calcule l'offset ----
  #Description :
  #Parametre : variable dans covariate
  # Fonction pour calculer la matrice O en utilisant le nom de la colonne
calcul_Offset <- function(data, column_name = NULL) {
    if (!is.null(data)) {
      if (!"Abundance" %in% colnames(data)) {
        stop("La variable 'Abundance' n'existe pas dans les données.")
      }
      if (!is.null(column_name) && !column_name %in% colnames(data)) {
        stop("La colonne spécifiée n'existe pas dans les données.")
      }
    } else {
      stop("Les données ne sont pas spécifiées.")
    }
    
    if (!is.null(column_name)) {
      column_selected <- data[[column_name]]
      O <- matrix(rep(column_selected, ncol(data$Abundance)), ncol = ncol(data$Abundance), nrow = nrow(data$Abundance))
    } else {
      O <- matrix(rep(0, nrow(data$Abundance) * ncol(data$Abundance)), nrow = nrow(data$Abundance), ncol = ncol(data$Abundance))
    }
    
    return(O)
  }
  
  #----
predict_sicpln <- function(Covariate, Offset , params){       #GDR : 27-3-2024 : j'ai remplace le nom de la fonction fitted_pln par Error_var
  # data = data_to_predict
  # Offset = offsett
  # params = list(B = res_fs$model_par$B,S = res_fs$var_par$S , M = res_fs$var_par$M)
  ### Initialisation
  n <- nrow(data$Abundance)  #Nombre d'observations
  index <- torch_tensor(1:n)
  
  ## Conversion of parameters to torch tensors (pointers)
  params <- lapply(params, torch_tensor, requires_grad = F)
  
  
  #Calcul des Coeffs pour les donnees ajustees
  S2 <- torch_square(params$S[index])
  Z <- torch_tensor(Offset)[index] + params$M[index] + torch_matmul(torch_tensor(Covariate)[index], params$B)
  fitted_values <- (torch_exp(Z + .5 * S2))
  fitted_values <- as.matrix(fitted_values)
  colnames(fitted_values) <- colnames(data$Abundance)
  
  return(fitted_values)
}

compute_r_squared <- function(observed_values, predicted_values) {
  # Calcul de la moyenne des valeurs observées
  mean_observed <- mean(observed_values)
  
  # Calcul de la somme des carrés des différences entre les valeurs observées et la moyenne
  total_sum_squares <- sum((observed_values - mean_observed)^2)
  
  # Calcul de la somme des carrés des résidus
  residual_sum_squares <- sum((observed_values - predicted_values)^2)
  
  # Calcul du coefficient de détermination R²
  r_squared <- 1 - (residual_sum_squares / total_sum_squares)
  
  return(r_squared)
}

### Fonction compute_fisher_opt_pln : ----
## Description : ----
   # Introduit le calcul sur la selection de variables avec SIC

## Exemple a faire : ----
   #cfr en bas de la fonction SICPLN

## Fonction compute_fisher_opt_pln : ----
  # Introduit le calcul sur la selection de variables avec SIC dans le modèle PLN
compute_fisher_opt_pln <- function(Covariate, Abundance, O, numb_eps, maxIter, plngendata) {
  
  # Covariate = Covariate; Abundance = Abundance;
  # O = offsett; numb_eps = numb_eps; maxIter = maxIter; plngendata = plngendata
  
  p <- ncol(Abundance)                      # nombre d'especes
  n <- nrow(Covariate)                      # nombre d'observation

  O <- O
  d <- ncol(Covariate) 
  w <- rep(1,n)                     # le poids de chaque observation
  B <- plngendata$model_par$B        # Le parametre  B generé par la fontion PLN          GDR 05-04-2024  : A revoir le commentaire
  Sigma <- plngendata$model_par$Sigma # Le parametre sigma generé par la fontion PLN
  M   <- plngendata$var_par$M           # Le parametre M genere par la fontion PLN
  S <- plngendata$var_par$S           # Le parametre S genere par la fontion PLN
  
  
  ##################################
  lambda_B <- log(n)*((p*(d))+p)     #GDR 29-03-2024 : coeffficient de penalite
  E <- c()                             # Initialisation du vecteur E(Pour le stockage des epsilones)
  e1 <- 10
  E[1] <- e1                           #affectation du premier valeur au vecteur E
  
  for(t in 2:numb_eps){E[t] <- e1*(0.87)^(t-1)}
  #calcul des  100 autres elements du vecteur
  
  vecepsilon <- sort(E, decreasing = TRUE)
  B_tol <- 0.00001
  dd <- matrix(NA, ncol = p)
  
  
  soltest <- tryCatch({ 
    for(i in 1:numb_eps) {
      epsilon <- vecepsilon[i]
      
      # Creation de bar de progression 
      progress(i, max.value =  numb_eps, progress.bar = TRUE)
      
      # convergence newton
      B_oldf <- B                         # stockage de l'ancien B  
      for(k in 1:maxIter)
      {
        Z <- O+M+as.matrix(torch_matmul(Covariate, B))
        A <- torch_exp(Z + torch_pow(S, 2)/2)
        A <- as.matrix(A)
        
        # Calcul Matrice d'information de fisher
        fisher <- Matrix::bdiag(lapply(1:p, function(j)
        {crossprod(Covariate, A[, j] * Covariate) })) 
        fisher <- as.matrix(fisher)
        
        # Calcul de la matrice de score
        score <- Matrix::bdiag(lapply(1:p, function(j)
        {crossprod(Covariate, Abundance[, j]-A[, j])}))       
        score=as.matrix(score)
        
        # mise à jour
        Sigma_B <- diag(c(rep(0,ncol(Abundance)),
                          (lambda_B / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (B[-1,]^2))) / ((B[-1,]^2 + epsilon^2)^3))))
        
        v_B <- (lambda_B/2)*(((2 * (epsilon^2))*(B[-1,]))) / ((B[-1,]^2 + epsilon^2)^2)
        v_B <- rbind(rep(0,ncol(Abundance)),v_B)
        
        #%%%%%%%%%%%%%%%%%%% Transformation en bloc
        #%%%% v_B en bloc
        v_B <- (lambda_B/2)*(((2 * (epsilon^2))*(B[-1,]))) / ((B[-1,]^2 + epsilon^2)^2)
        v_B <- rbind(rep(0,ncol(Abundance)),v_B)
        
        v_bloc <- Matrix::bdiag(lapply(1:p, function(j)
        {v_B[,j] }))  # t(X) %*% diag(A[, i]) %*% X
        v_bloc=as.matrix(v_bloc)
        
        #%%%% Sigma en bloc
        Sigma_B <- (lambda_B / 2) * (((2 * (epsilon^2)) * (epsilon^2 - 3 * (B[-1,]^2))) / ((B[-1,]^2 + epsilon^2)^3))
        
        Sigma_B <- rbind(rep(0,p),Sigma_B)
        Sigma_Bt <- diag(c(Sigma_B))
        
        Sig_bloc <- Matrix::bdiag(lapply(1:p, function(j){Sigma_B[,j] }))          # t(X) %*% diag(A[, i]) %*% X
        Sig_bloc=as.matrix(Sig_bloc)
        
        #%%%%%%%%%%%%%%%%%% B en bloc
        B_bloc <- Matrix::bdiag(lapply(1:p, function(j)
        {B[,j] }))
        
        B_bloc <- as.matrix(B_bloc)
        
        fish_Sigma <- fisher+Sigma_Bt
        score_v_bloc <- score-v_bloc
        
        # QR decomposition de pour calculer l'inverse de l'information de fisher
        qr_decomposition <- qr(fish_Sigma)
        
        # Matrice Q
        Q <- qr.Q(qr_decomposition)
        
        # Matrice R
        R <- qr.R(qr_decomposition)
        
        # Calcul de l'inverse de R
        inverse_R <- solve(R)
        
        # Transposée de Q
        transpose_Q <- t(Q)
        
        # Calcul de l'inverse de l'information de fisher
        inverse_fish_Sigma <-inverse_R %*% transpose_Q
      
        # mise à jour de B
        B_new <- as.matrix(B_bloc+inverse_fish_Sigma%*%(score_v_bloc))
        
        #======================
        ki = 1
        d <- ncol(Covariate)       # nombre de covariables 
        B_newmat <- matrix(NA, ncol=ncol(Abundance), nrow=nrow(B))
        for(t in 1:ncol(Abundance)){
          
          B_newmat[,t]=B_new[ki:d,t]
          
          ki=d+1
          d=d+ncol(Covariate)
          
        }
        
        diff_B <- sum(abs(B_newmat-B_oldf))
        if(diff_B <= B_tol)
        {
          for(bi in 2:nrow(B_newmat))
          {
            for(bj in 1:ncol(B_newmat))
            {
              if(abs(B_newmat[bi,bj]) <= B_tol) {B_newmat[bi,bj]=0}
            }
          }
          
          B <- B_newmat
          break
        }
        if(k==maxIter)
        {
          B <- B_newmat
          B_oldf <- B
        }
        
        B <- B_newmat
        B_oldf <- B
      }
      dd <- rbind(dd,B[-1,])
      
      varpar <-  plngendata$optimize_vestep(Covariate, O, Abundance,weights=rep(1,n) ,B = B, Omega = plngendata$model_par$Omega,
                                            control = PLN_param(backend = "nlopt"))
      M <- varpar$M
      S <- varpar$S
      B <- B
      Sigma <- plngendata$model_par$Sigma
      Omega <- plngendata$model_par$Omega
      
    }
    model_par <- list(B = B, Sigma = Sigma,Omega = Omega)
    var_par <- list(M = M, S = plngendata$var_par$S, S2 = plngendata$var_par$S2^2)
    res1 <- list(model_par = model_par, var_par = var_par,vcov_model = plngendata$vcov_model, Offset = O, Covariate = Covariate, dd = dd)
    return(res1)
    
  } , error = function(e){ return(res1) }) #Spécification de la classe d'erreur et gestion de l'erreur
  cat("\t\n Execution of SICPLN is DONE! \n")
  return(soltest) 
}
# -----



# fonction SICPLN : ----
### Description : ----
    #    initialisation PLN + SIC
###   Paramètres de la fonction :
   #   - Covariate : Matrice des covariables
   #   - Abundance : Matrice des Abundances
   #   - O : Matrice des offsets (par défaut une matrice de zéros)

###  Exemple a faire : ----
    #cfr en bas
###  Fonction SICPLN : -----
#interactive=TRUE,covariance_structure="full"
SICPLN <- function(formule, offset_column_name , data , numb_eps, maxIter, control = control) {
  if (!is.null(offset_column_name)) {
    offsett <- calcul_Offset(data, column_name = offset_column_name)
    offset_values <- data[[offset_column_name]]
    
    data <- data[, colnames(data) != offset_column_name]
    
    Covariate <- model.matrix(formule, data)
    Abundance <- model.response(model.frame(formule, data))
    
    # Création du modèle PLN avec les spécifications fournies
    capture.output({plngendata <- PLN(Abundance ~ . + offset(offset_values), data = data, control = control)}, file = NULL)
    colnames(offsett) <- colnames(plngendata$model_par$B)
    
    res_fs <- compute_fisher_opt_pln(Covariate = Covariate, Abundance = Abundance,
                                     O = offsett, numb_eps = numb_eps, maxIter = maxIter, plngendata = plngendata)
    
    colnames(res_fs$model_par$B) <- colnames(plngendata$model_par$B)
    rownames(res_fs$model_par$B) <- rownames(plngendata$model_par$B)
    
    predictions <- predict_sicpln(Covariate = Covariate , Offset = offsett ,params = list(B = res_fs$model_par$B,S = res_fs$var_par$S , M = res_fs$var_par$M))
    predict_pln <- predict_sicpln(Covariate = Covariate , Offset = offsett ,params = list(B = plngendata$model_par$B,S = plngendata$var_par$S , M = plngendata$var_par$M))
    residuals <- sum(abs(data$Abundance - predictions))
    residuals_pln <- sum(abs(data$Abundance - predict_pln))
    Rsquared <- compute_r_squared(Abundance , predictions)
    # Résultats du modèle
    resf <- list(res_pln = plngendata, res_fisher = res_fs, fitted = predictions , residuals = residuals , residuals_pln = residuals_pln, R_sqared = Rsquared)
  } else {
    offsett <- calcul_Offset(data, column_name = offset_column_name)
    
    Covariate <- model.matrix(formule, data)
    Abundance <- model.response(model.frame(formule, data))
    
    # Création du modèle PLN avec les spécifications fournies
    
    capture.output({plngendata <- PLN(Abundance ~ ., data = data, control = control)}, file = NULL)
    colnames(offsett) <- colnames(plngendata$model_par$B)
    
    res_fs <- compute_fisher_opt_pln(Covariate = Covariate, Abundance = Abundance,
                                     O = offsett, numb_eps = numb_eps, maxIter = maxIter, plngendata = plngendata)
    
    colnames(res_fs$model_par$B) <- colnames(plngendata$model_par$B)
    rownames(res_fs$model_par$B) <- rownames(plngendata$model_par$B)
    
    predictions <- predict_sicpln(Covariate = Covariate , Offset = offsett ,params = list(B = res_fs$model_par$B,S = res_fs$var_par$S , M = res_fs$var_par$M))
    predict_pln <- predict_sicpln(Covariate = Covariate , Offset = offsett ,params = list(B = plngendata$model_par$B,S = plngendata$var_par$S , M = plngendata$var_par$M))
    residuals <- sum(abs(data$Abundance - predictions))
    residuals_pln <- sum(abs(data$Abundance - predict_pln))
    Rsquared <- compute_r_squared(Abundance , predictions)
    # Résultats du modèle
    resf <- list(res_pln = plngendata, res_fisher = res_fs, fitted = predictions , residuals = residuals , residuals_pln = residuals_pln, R_sqared = Rsquared)
  }
  cat("\n\t")
  return(resf)  # Retourne les résultats
}
  
  # ----
  
# Exemple a faire pour la fonction SICPLN &  Compute-fisher : ----
   # # Création de données factices pour l'exemple
    #set.seed(123)  # Pour la reproductibilité
    #n <- 100  # Nombre d'observations
    #p <- 3  # Nombre de covariables
    #d <- 2  # Nombre d'espèces
    
    # Génération des covariables (matrice aléatoire)
     #Covariate <- matrix(runif(n * p), nrow = n, ncol = p)
    
    # Génération des Abundances (matrice aléatoire)
    #Abundance <- matrix(rpois(n * d, lambda = 10), nrow = n, ncol = d)
    
    # Appel de la fonction SICPLN avec les données factices
    #result <- SICPLN(Covariate, Abundance)
    
    # Affichage des résultats
    #print(result)

   
   
# ----
# Fonction MSE (Mean Squared Error) ----
## Description : ----
   # Fonction qui calcule la moyenne des erreurs au carré entre les vraies valeurs et les valeurs prédites.
   
   # Paramètres initiaux :
   # vrai_comptage : les données de comptage réelles (data$Y)
   # comptage_predicte : les données de comptage prédites
   # renvoie le MSE
##  Exemple a executer : ----
   # # Données de comptage réelles et prédites
   # vrai_comptage <- c(10, 20, 30, 40)
   # comptage_predicte <- c(12, 18, 32, 38)
   # 
   # # Calcul de l'erreur quadratique moyenne
   # mse <- mse_predict(vrai_comptage, comptage_predicte)
   # 
   # # Affichage du résultat
   # print(paste("MSE:" , mse))
   
   
## fonction mse_predict : ----  
mse_predict <- function(vrai_comptage, comptage_predicte) {
  # Vérifie que les deux vecteurs ou matrices ont la même dimension
  if (!all(dim(vrai_comptage) == dim(comptage_predicte))) {
    stop("Les dimensions de vrai_comptage et comptage_predicte doivent être les mêmes")
  }
  # Calcul de l'erreur quadratique moyenne
  mse <- mean((vrai_comptage - comptage_predicte) ^ 2)
  
  return(mse)  # Retourne la valeur de l'erreur quadratique moyenne
}


# ----
   
# Erreur de prédiction en tenant compte des paramètres variationnels : ----
   
#GDR : 27-3-2024 :commentaire de Jean-Yves --->   # Erreur de prédiction en tenant compte des paramètres variationnels
### Description : ---- 
# Erreur de prédiction en tenant compte des paramètres variationnels
## Parametres initiales :
#   data : données du modèle
#   params : paramètres du modèle
# Renvoie l'erreur apres varitionnal EM
##  Exemple a executer : ----
   
## fonction Error_var : ----
Error_var<-function(data, params){       #GDR : 27-3-2024 : j'ai remplace le nom de la fonction fitted_pln par Error_var
     ### Initialisation
     p <- ncol(data$Abundance)  #nombre d'especes
     d <- ncol(data$Covariate)  #Nombre de variables explicatives(Covariates)
     n <- nrow(data$Covariate)  #Nombre d'observations
     index <- torch_tensor(1:n)
     
     ## Conversion of data and parameters to torch tensors (pointers)
     data   <- lapply(data, torch_tensor)
     params <- lapply(params, torch_tensor, requires_grad = F)
     
     
     #Calcul des Coeffs pour les donnees ajustees
     S2 <- torch_square(params$S[index])
     Z <- data$O[index] + params$M[index] + torch_matmul(data$Covariate[index], params$B)
     res_tmp <- (torch_exp(Z + .5 * S2))
     
     
     return(as.matrix(res_tmp))
   }
  # ----  
   
# calcul erreur de prediction ( Ou Taux de certitude au niveau de Beta) : ----
## Description : ----
   # Fonction qui calcule l'erreur de prediction.
## Paramètres initiaux :
   # vrai_beta : les données de comptage réelles (data$Y)
   # beta_estime : les données de comptage prédites
   # Valeur renvoyee : Taux de certitude pour les Beta
##  Exemple a executer : ----
   # # Données de vrais beta et beta estimés
   # vrai_beta <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2, byrow = TRUE)
   # beta_estime <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 2, byrow = TRUE)
   # 
   # # Calcul du taux de certitude pour les Beta
   # taux_certitude <- sparsity_recognition(vrai_beta, beta_estime)
   # 
   # # Affichage des résultats
   # print(paste("Taux de certitude pour les Beta - TNR:", taux_certitude$TNR))
   # print(paste("Taux de certitude pour les Beta - TPR:", taux_certitude$TPR))
   
## fonction sparsity_recognition : ----
sparsity_recognition=function(vrai_beta, beta_estime){
     TN <- 0            # Effectifs des vraies zeros
     TP <- 0            # Effectifs des elements non-zeros
     for(i in 1:nrow(vrai_beta))
     {
       for(j in 1:ncol(vrai_beta))
       {
         if(vrai_beta[i,j]!=0 & beta_estime[i,j]!=0){TP <- TP+1}
         
         if(vrai_beta[i,j]==0 & beta_estime[i,j]==0){TN <- TN+1}
       }
     }
     # Taux de non zero estimée estimée zero
     TPR <- TP/(length(which(vrai_beta!=0)))
     # Taux de zero estimée zéro
     TNR <- TN/(length(which(vrai_beta==0)))
     
     # stockage des taux de certitude dans une liste
     sol <- list(TNR=TNR,TPR=TPR)
     return(sol)
   }
# ---- 
   
# Avec pln et optimisation des paramètre variationnellle : ----
## Description : ----
   #
   
## Exemple a executer :----
   # # Définir les paramètres
   # nombre_espece <- 10
   # taille_echant <- 100
   # nbre_repetition <- 5
   # 
   # # Appeler la fonction error_cal avec les paramètres spécifiés
   # resultats <- error_cal(nombre_espece, taille_echant, nbre_repetition)
   # 
   # # Afficher les résultats
   # print(resultats)
   
   ## Fonction error_cal : ----
error_cal <- function(nombre_espece , taille_echant, nbre_repetition) {
     
     # Initialisation des paramètres
     nombre_espece <- p
     taille_echant <- n
     sp <- nombre_espece      # Nombre d'espèces
     n <- taille_echant
     
     nb_times <- nbre_repetition       # Nombre de répétitions pour la simulation
     
     # Vecteurs pour stocker les erreurs pour chaque itération
     err_sic_sp <- numeric(nb_times)
     err_pln_sp <- numeric(nb_times)
     err_glm_sp <- numeric(nb_times)
     
     # Stockage des MSE pour chaque méthode
     MSE_SIC <- numeric(nb_times)
     MSE_PLN <- numeric(nb_times)
     MSE_GLMNET <- numeric(nb_times)
     
     # Stockage des taux de faux positifs et faux négatifs
     TPR_SIC <- numeric(nb_times)
     TPR_PLN <- numeric(nb_times)
     TPR_GLMNET <- numeric(nb_times)
     TNR_SIC <- numeric(nb_times)
     TNR_PLN <- numeric(nb_times)
     TNR_GLMNET <- numeric(nb_times)
     
     # Initialisation des matrices pour stocker les coefficients estimés
     list_mat_B_sp <- list(SICPLN = matrix(NA, 1, sp + 2),
                           PLN = matrix(NA, 1, sp + 2),
                           GLMNET = matrix(NA, 1, sp + 2))
     
     # Boucle sur le nombre de répétitions
     for (t in 1:nb_times) {  
       # Code pour chaque répétition de simulation
       gen_data <- simul_Abundance(B = beta, Sigma = Sigma, Covariate = Covariate)
       Covariate <- gen_data$Covariate; Abundance <- gen_data$Abundance; p <- ncol(gen_data$Abundance); d <- ncol(gen_data$Covariate); n <- nrow(gen_data$Abundance)
       res_SICPLN <- SICPLN(Covariate , Abundance , O = matrix(0,ncol = p,nrow = n), interactive = FALSE, covariance_structure = "full")
       res_glmpois <- glmnet_lasso_poisson(Abundance, Covariate)
       
       # Calcul des erreurs pour chaque méthode
       err_sic_sp[t] <- norm(((beta)-(res_SICPLN$res_fisher$B[-1,])),"F")/norm(beta,"F")
       err_pln_sp[t] <- norm(((beta)-(res_SICPLN$res_pln$B[-1,])),"F")/norm(beta,"F")
       err_glm_sp[t] <- norm(((beta)-(res_glmpois$hat_B[-1,])),"F")/norm(beta,"F")
       
       # Calcul des MSE
       # Prediction SICPLN
       # params=params1
       data_forma_pred <- list(Covariate=cbind(rep(1,n),Covariate),Abundance=Abundance
                               ,O=res_SICPLN$res_fisher$O)
       data <- data_forma_pred
       params1 <- list(B=res_SICPLN$res_fisher$B,
                       S=res_SICPLN$res_fisher$S, M=res_SICPLN$res_fisher$M)
       n <- nrow(res_SICPLN$res_fisher$M)
       pred_sic <- Error_var(data=data_forma_pred, params=params1)  #Couche latente Z avec SIC
       
       # Prediction PLN
       # data_forma_pred$Y
       params2 <- list(B=res_SICPLN$res_pln$B,S=res_SICPLN$res_pln$S
                       ,M=res_SICPLN$res_pln$M)
       pred_pln <- Error_var(data=data_forma_pred, params=params2)  #couche latente Z avec PLNmodel 
       MSE_SIC[t] <- mse_predict(gen_data$Abundance,pred_sic)  #Erreur genere en utilisant SIC
       MSE_PLN[t] <- mse_predict(gen_data$Abundance,pred_pln)  #Erreur genere en utilisant PLNmodel
       MSE_GLMNET[t] <- mse_predict(gen_data$Abundance,res_glmpois$prediction)  #Erreur genere en utilisant GLMNET (la fonction calcule ses valeurs predicts)
       
       
       # Calcul des taux de faux positifs et faux négatifs
       SIC_sparsity <- sparsity_recognition(beta,res_SICPLN$res_fisher$B[-1,])
       GLMNET_sparsity <- sparsity_recognition(beta,res_glmpois$hat_B[-1,])
       PLN_sparsity <- sparsity_recognition(beta,res_SICPLN$res_pln$B[-1,])
       TPR_SIC[t] <- SIC_sparsity$TPR
       TPR_PLN[t] <- PLN_sparsity$TPR
       TPR_GLMNET[t] <- GLMNET_sparsity$TPR
       TNR_SIC[t] <- SIC_sparsity$TNR
       TNR_PLN[t] <- PLN_sparsity$TNR
       TNR_GLMNET[t] <- GLMNET_sparsity$TNR
       
       # Stockage des résultats dans les matrices
       
       list_mat_B_sp[["SICPLN"]] <- rbind(list_mat_B_sp[[1]],cbind(abs(beta-res_SICPLN$res_fisher$B[-1,]),rep(n,d),rep(sp,d)))
       list_mat_B_sp[["PLN"]] <- rbind(list_mat_B_sp[[2]],cbind(abs(beta-res_SICPLN$res_pln$B[-1,]),rep(n,d),rep(sp,d)))
       list_mat_B_sp[["GLMNET"]] <- rbind(list_mat_B_sp[[3]],cbind(abs(beta-res_glmpois$hat_B[-1,]),rep(n,d),rep(sp,d)))
       
     }
     
     
     df_sp_size=data.frame(nom_fb=c(err_sic_sp,err_pln_sp, err_glm_sp),
                           lab_methode = as.factor(rep(c("SICPLN","PLN","GLMNET"),each=length(err_sic_sp))),
                           ech_size=rep(c(n,n,n),each=length(err_sic_sp)),
                           sp_size=as.factor(rep(c(sp,sp,sp),each=length(err_sic_sp))))
     
     names(df_sp_size)=c("nom_fb","lab_methode","ech_size","sp_size")
     
     # Calcul de la moyenne des erreurs sur toutes les répétitions
     err_sic_sp_mean <- mean(err_sic_sp)
     err_pln_sp_mean <- mean(err_pln_sp)
     err_glm_sp_mean <- mean(err_glm_sp)
     
     # Calcul de la moyenne des taux de vrais positifs et vrais négatifs
     TPR_SIC_mean <- mean(TPR_SIC)
     TPR_PLN_mean <- mean(TPR_PLN)
     TPR_GLMNET_mean <- mean(TPR_GLMNET)
     TNR_SIC_mean <- mean(TNR_SIC)
     TNR_PLN_mean <- mean(TNR_PLN)
     TNR_GLMNET_mean <- mean(TNR_GLMNET)
     
     # Calcul de la moyenne des MSE sur toutes les répétitions
     MSE_SIC_mean <- mean(MSE_SIC)
     MSE_PLN_mean <- mean(MSE_PLN)
     MSE_GLMNET_mean <- mean(MSE_GLMNET)
     # 
     
     # Création des structures de données pour les résultats
     SICPLN_matb <- as.matrix(list_mat_B_sp$SICPLN)
     PLN_matb <- as.matrix(list_mat_B_sp$PLN)
     GLM_matb <- as.matrix(list_mat_B_sp$GLMNET)
     head(list_mat_B_sp$SICPLN)
     GLM_bpt <- as.numeric(GLM_matb[,1:sp])
     SIC_bpt <- as.numeric(SICPLN_matb[,1:sp])
     PLN_bpt <- as.numeric(PLN_matb[,1:sp])
     lab_methode <- rep(c("GLMNET","SICPLN","PLN"),each=length(PLN_bpt))
     lab_ech_size <- rep(n,length(lab_methode))
     lab_sp_size <- rep(sp,length(lab_methode))
     df_bp <- data.frame(coef=c(GLM_bpt,SIC_bpt,PLN_bpt),lab_methode=lab_methode,ech_size=lab_ech_size,sp_size=lab_sp_size)
     df_bp$lab_methode <- as.factor(df_bp$lab_methode)
     df_bp$ech_size <- as.factor(df_bp$ech_size)
     df_bp$sp_size <- as.factor(df_bp$sp_size)
     MSE <- list(MSE_SIC_mean = MSE_SIC_mean, MSE_PLN_mean = MSE_PLN_mean, MSE_GLMNET_mean = MSE_GLMNET_mean)
     MSE_replication <- list(MSE_SIC = MSE_SIC, MSE_PLN = MSE_PLN, MSE_GLMNET = MSE_GLMNET)
     Performance <- list(TPR_SIC_mean=TPR_SIC_mean, TPR_PLN_mean=TPR_PLN_mean, TPR_GLMNET_mean=TPR_GLMNET_mean,TNR_SIC_mean=TNR_SIC_mean,TNR_PLN_mean=TNR_PLN_mean,TNR_GLMNET_mean=TNR_GLMNET_mean)
     performance_replication <- list(TPR_SIC=TPR_SIC,TPR_PLN=TPR_PLN,TPR_GLMNET=TPR_GLMNET,TNR_SIC=TNR_SIC,TNR_PLN=TNR_PLN,TNR_GLMNET=TNR_GLMNET)
     
     
     # Retourner les résultats
     return(list(
       df_bp = df_bp,
       erreur_moy_sic = err_sic_sp,
       erreur_moy_pln = err_pln_sp,
       erreur_moy_glm = err_glm_sp,
       list_matrice_coef = list_mat_B_sp,
       MSE = MSE,
       Performance = Performance
     ))
   }
   # ----    
   
 #######################   GRAPHIQUE   ############################
# boxplot
# parametres : dataframe contenant les resultats d'erreur en variant l'echantillon ou le nbre d'especes
# renvoie les boxplots

# Fonction boxplot_Error_norm : ----   
boxplot_error_norm <- function(Error_n_30, Error_n_50, Error_n_100, Error_n_1000, graph_parm = c("ech_size", "nom_fb", "lab_methode")) {
     df_combine_echantillon <- rbind(Error_n_30$df_bp, Error_n_50$df_bp, Error_n_100$df_bp, Error_n_1000$df_bp)
     
     plt_echantillon <- ggplot(df_combine_echantillon, aes(x = as.factor(ech_size), y = nom_fb, fill = lab_methode, dodge.position)) + 
       geom_boxplot(aes(fill = lab_methode), notch = FALSE, position = position_dodge(0.9)) +
       scale_fill_manual(values = c("lightgray", "#d1ab75", "#617a89"), name = "Estimation method") +
       theme(
         plot.title = element_text(size = 20, hjust = 0.5),
         axis.title = element_text(size = 20, colour = "black"),
         legend.text = element_text(size = 10),
         axis.text.x = element_text(size = 8),
         strip.text.x = element_text(size = 20, colour = "black"),
         axis.text.y = element_text(size = 10, colour = "black")
       ) +
       labs(x = "Sample size", y = "Coefficients errors", title = "") +
       scale_x_discrete(labels = c("n = 30", "n = 50", "n = 100", "n = 1000"))
     
     return(plt_echantillon)
   }
# ----   

# Fonction  plot_error_prediction :   
# description :----
   # parametres : les differentes resultats sur le erreurs
   
# Fonction  plot_error_prediction : ----  
plot_error_prediction <- function(Error_n_30, Error_n_50, Error_n_100, Error_n_1000) {
     df_error_prediction_echantillon <- data.frame(
       SICPLN = c(Error_n_30$MSE$MSE_SIC, Error_n_50$MSE$MSE_SIC, Error_n_100$MSE$MSE_SIC, Error_n_1000$MSE$MSE_SIC),
       PLN = c(Error_n_30$MSE$MSE_PLN, Error_n_50$MSE$MSE_PLN, Error_n_100$MSE$MSE_PLN, Error_n_1000$MSE$MSE_PLN),
       GLMNET = c(Error_n_30$MSE$MSE_GLMNET, Error_n_50$MSE$MSE_GLMNET, Error_n_100$MSE$MSE_GLMNET, Error_n_1000$MSE$MSE_GLMNET)
     )
     
     df_error_prediction_echantillon$x <- c(30, 50, 100, 1000)
     
     plot_error_prediction_ech <- ggplot(df_error_prediction_echantillon, aes(x)) +
       geom_point(aes(y = PLN, color = "PLN"), size = 3) +   # Points pour y1 (en bleu)
       geom_line(aes(y = PLN, color = "PLN"), show.legend = TRUE) +              # Ligne pour y1 (en bleu)
       geom_point(aes(y = GLMNET, color = "GLMNET"), size = 3) +    # Points pour y2 (en rouge)
       geom_line(aes(y = GLMNET, color = "GLMNET")) +
       geom_point(aes(y = SICPLN, color = "SICPLN"), size = 3) +    # Points pour y3 (en vert)
       geom_line(aes(y = SICPLN, color = "SICPLN"), show.legend = TRUE) + # Ligne pour y3 (en vert)
       scale_y_log10() +
       labs(x = "Sample size", y = "Prediction error", title = "",
            color = "") +
       scale_color_manual(values = c("PLN" = "blue", "GLMNET" = "red", "SICPLN" = "green"), 
                          labels = c("GLMNET", "PLN", "SICPLN")) +
       theme_minimal() +
       theme(
         legend.position = "top",
         plot.title = element_text(size = 20, hjust = 0.5),
         axis.title = element_text(size = 20, colour = "black"),
         legend.text = element_text(size = 15),
         axis.text.x = element_text(size = 8),
         strip.text.x = element_text(size = 20, colour = "black"),
         axis.text.y = element_text(size = 10, colour = "black")
       ) +
       labs(x = " Sample size", y = "Prediction error", title = "")
     
     print(plot_error_prediction_ech)
   }
# ---- 

#################################  GENUS2 #####################################
#######################################
# Fonction glmnet_lasso_poisson_withoffset  ----
####### Destription : ----
#   Applique le modele GLMNET sur les données(Abundance, Covariate, offset)
#   les utilisent comme entrée pour la fonction glmnet_lasso_poisson_withoffset
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
# # Appel de la fonction glmnet_lasso_poisson_withoffset avec les données d'exemple
#resultats <- glmnet_lasso_poisson_withoffset(Abundance, Covariate)
# # Affichage des coefficients estimés
#print("Coefficients estimés :")
#print(resultats$coefficient)


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

# Fonction calculate_VIF : ----
calculate_VIF <- function(Covariate, seuil, pval_threshold = 0.05) {
  # Créer un vecteur pour stocker les résultats des VIF
  VIF_values <- rep(NA, ncol(Covariate))
  
  # Créer une matrice pour stocker les p-values du test du chi-carré pour les variables catégorielles
  cat_pvals <- matrix(NA, ncol(Covariate), ncol(Covariate))
  colnames(cat_pvals) <- colnames(Covariate)
  rownames(cat_pvals) <- colnames(Covariate)
  
  # Parcourir chaque variable
  for (i in 1:ncol(Covariate)) {
    # Sélectionner la variable dépendante (y)
    y <- Covariate[, i]
    
    if (is.factor(y) || is.character(y)) {
      # Tester la colinéarité entre variables catégorielles
      for (j in 1:ncol(Covariate)) {
        if (j != i && (is.factor(Covariate[, j]) || is.character(Covariate[, j]))) {
          # Effectuer le test du chi-carré
          chi2_test <- chisq.test(table(y, Covariate[, j]))
          cat_pvals[i, j] <- chi2_test$p.value
        }
      }
      # Ne pas calculer le VIF pour les variables catégorielles
      VIF_values[i] <- NA
    } else {
      # Sélectionner toutes les autres variables indépendantes (X_i) en excluant les variables catégorielles
      X_i <- Covariate[, -i]
      
      # Créer un nouveau dataframe avec y et X_i
      df <- data.frame(y = y, X_i)
      
      # Calculer le coefficient de détermination (R²) du modèle de régression linéaire
      R_squared <- summary(lm(y ~ ., data = df))$r.squared 
      
      # Calculer le VIF à partir de R²
      VIF_values[i] <- ifelse(R_squared == 0, NA, 1 / (1 - R_squared))
    }
  }
  
  # Attribution des noms aux valeurs
  names(VIF_values) <- colnames(Covariate)
  
  # Identifier les variables avec un VIF inférieur au seuil
  selected_variables <- colnames(Covariate)[VIF_values < seuil]
  
  # Identifier les variables catégorielles avec des p-values significatives
  cat_pairs <- which(cat_pvals < pval_threshold, arr.ind = TRUE)
  cat_vars_to_remove <- unique(c(colnames(Covariate)[cat_pairs[, 1]], colnames(Covariate)[cat_pairs[, 2]]))
  
  # Filtrer les valeurs NA de selected_variables
  selected_variables <- selected_variables[!is.na(selected_variables)]
  
  # Concaténer les variables catégorielles non colinéaires
  selected_variables <- unique(c(selected_variables, colnames(Covariate)[!colnames(Covariate) %in% cat_vars_to_remove & sapply(Covariate, is.factor)]))
  
  # Créer une liste pour stocker les résultats
  results <- list(VIF_values = VIF_values, selected_variables = selected_variables, cat_pvals = cat_pvals)
  
  # Si le max des VIF est supérieur au seuil
  if (max(VIF_values, na.rm = TRUE) > seuil) {
    # Identifier l'index de la variable avec le VIF le plus élevé
    index_to_remove <- which.max(VIF_values)
    
    # Enlever la variable correspondante de Covariate
    Covariate <- Covariate[, -index_to_remove]
    
    # Rappeler la fonction pour recalculer les VIF sans cette variable
    results <- calculate_VIF(Covariate, seuil, pval_threshold)
  }
  
  # Retourner les résultats
  return(results)
}
# ----


# Définition de la fonction pour tracer les coefficients : ----
plot_coefficient <- function(coef_matrix, axis_names = c("Genus", "Variables", "Coefficient Value", "Title"), grad_min, grad_max) {
  # Renommer les axes et le gradient de couleur à partir des paramètres d'entrée
  #grad_min <- -3
  #grad_max <- 3
  y_axis_name <- axis_names[2]
  x_axis_name <- axis_names[1]
  color_gradient_name <- axis_names[3]
  title <- axis_names[4]
  
  # Convertir la matrice de coefficients en un format long pour ggplot
  melted_coef_matrix <- melt(coef_matrix)
  
  # Création du graphique avec ggplot
  plot <- ggplot(melted_coef_matrix, aes(x = as.factor(Var2), y = as.factor(Var1), fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "darkred",high = "darkgreen",mid = "white", midpoint = 0, limit = c(grad_min, grad_max), name = color_gradient_name) +
    coord_equal() +
    labs(x = x_axis_name, y = y_axis_name, title = title, color = "black") +
    theme(legend.position = "right", plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 20), legend.text = element_text(size = 15), strip.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black"), axis.text.x = element_text(size = 20, colour = "black"))
  
  # Retourner le graphique
  return(plot)
}
#----


# Définition de la fonction pour tracer les coefficients :----
coef_genus_barre <- function(coef_matrix, axis_names = c("Genus", "Variables", "Coefficient Value", "Title"), grad_min, grad_max) {
  # Renommer les axes et le gradient de couleur à partir des paramètres d'entrée
  y_axis_name <- axis_names[2]
  x_axis_name <- axis_names[1]
  color_gradient_name <- axis_names[3]
  title <- axis_names[4]
  
  # Convertir la matrice de coefficients en un format long pour ggplot
  melted_coef_matrix <- melt(coef_matrix)
  melted_coef_matrix$sparsity <- (ifelse(melted_coef_matrix$value == 0, 'Zero', "Nonzero"))
  
  # Création du graphique avec ggplot
  plot <- ggplot(melted_coef_matrix, aes(x =as.factor(Var2), y = as.factor(Var1),pattern = sparsity, fill = value)) +
    geom_tile(aes(fill = value), colour = "white") + 
    geom_tile_pattern(pattern_color = NA, pattern_fill = "black", pattern_angle = 45,
                      pattern_density = 0.15, pattern_spacing = 0.065,pattern_key_scale_factor = 1) +
    scale_pattern_manual(values = c(Zero = "circle", Nonzero = "none"),name="") +
    scale_fill_gradient2( low = "darkred",high = "darkgreen",mid = "white", midpoint = 0, limit = c(grad_min, grad_max), name = color_gradient_name) +
    coord_equal() + labs(x = x_axis_name, y = y_axis_name, title = title, color = "black") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")))+
    theme(legend.position = "right", plot.title = element_text(size = 10, hjust = 0.5), axis.title = element_text(size = 10), legend.text = element_text(size = 10), strip.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 10, colour = "black"), axis.text.x = element_text(size = 10, colour = "black"))
  
  # Retourner le graphique
  return(plot)
}
#----
# MAtrice de precision
plot_precision_matrix <- function(coef_matrix, 
                                  axis_names = c("Genus","Variables","Coefficients value","SICPLN"),
                                  limit = c(-1,1)){
  # Renommer les axes et le gradient de couleur à partir des paramètres d'entrée
  y_axis_name <- axis_names[2]
  x_axis_name <- axis_names[1]
  color_gradient_name <- axis_names[3]
  title <- axis_names[4]
  
  # Convertir la matrice de coefficients en un format long pour ggplot
  melted_coef_matrix <- melt(coef_matrix)
  
  # Création du graphique avec ggplot
  plot <- ggplot(melted_coef_matrix, aes(Var1, Var2,value)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient2(low = "darkred", high = "darkgreen",mid="white", midpoint = 0, limit = limit,name=color_gradient_name)+ labs(x = x_axis_name, y = y_axis_name, title = title,color = "black")+theme(legend.position = "right",plot.title = element_text(size = 20,hjust = 0.5),axis.title = element_text(size = 15),legend.text = element_text(size = 12),strip.text.x = element_text(size = 15, colour = "black"),axis.text.y = element_text(size = 15,colour = "black"),axis.text.x = element_text(size = 15,colour = "black"))
  return(plot)
}
#----

# Graphique Regularisation : ----
# Fonction pour générer les données du graphique des coefficients :----
datat_coef_path <- function(stockage_coef, numero_espece, nombre_variable, matrice_variables){
  # Nombre de variables
  d <- nombre_variable
  
  # Coefficients à afficher
  gg <- stockage_coef[-1,]
  
  # Numéro de l'espèce
  sp <- numero_espece
  
  # Calcul des limites pour l'axe y
  ymin <- min(gg[,sp])
  ymax <- max(gg[,sp])
  
  # Séquence d'epsilon
  sequence_eps = seq(1, 100, by = 1)
  
  # Initialisation de la matrice pour les données du graphique
  graph <- matrix(NA, nrow = nrow(gg), ncol = d + 1)
  colnames(graph) <- c(name_var, "epsilon")
  graph[, d + 1] <- sequence_eps
  
  # Extraction des données pour chaque variable et epsilon
  gg <- as.data.frame(gg)
  for(i in 1:(d)){
    graph[,i] <- gg[seq(i, nrow(gg), by = (d)), sp]
  }
  graph <- as.data.frame(graph) 
  
  # Création du graphique avec ggplot
  graph_espece <- ggplot(graph, aes(x = epsilon)) +
    geom_line(aes(x = epsilon, y = vegetation_index, color = "vegetation_index"), linewidth = 2) +
    geom_line(aes(x = epsilon, y = wetness, color = "wetness"), linewidth = 2) +
    geom_line(aes(x = epsilon, y = center_x, color = "center_x"), linewidth = 2) +
    geom_line(aes(x = epsilon, y = center_y, color = "center_y"), linewidth = 2) +
    labs(x = "Epsilon indices", y = "Coefficients", title = "",
         color = "") +
    scale_color_manual(name = "variables", 
                       values = c("vegetation_index" = "darkred",
                                  "wetness" = "darkblue",
                                  "center_x" = "darkgreen",
                                  "center_y" = "black")) +
    labs(title = paste0("Regularization path of genus ", sp)) +
    theme(plot.title = element_text(size = 25, hjust = 0.5),
          axis.title = element_text(size = 25, colour = "black"),
          legend.text = element_text(size = 25, colour = "black"),
          legend.title = element_text(color = "black", size = 25),
          legend.key.size = unit(1.5, "cm"),
          axis.text.x = element_text(size = 25, colour = "black"),
          strip.text.x = element_text(size = 25, colour = "black"),
          axis.text.y = element_text(size = 25, colour = "black"))
  
  return(graph_espece)
}

boxplot_graph <- function(data_df, x, y, title){
  plt <- suppressMessages(
    ggplot(data_df, aes(x = !!rlang::sym(x), y = !!rlang::sym(y), fill = !!rlang::sym(x))) +
      geom_boxplot() +
      labs(x = x, y = y, title = title) +
      scale_fill_discrete(name = x)
  )
  return(plt)
}


plot_density <- function(data1, data2, col1 = "red", col2 = "blue", main = "Density Plot", xlab) {
  dx <- density(data1)
  dy <- density(data2)
  
  df1 <- data.frame(x = dx$x, y = dx$y, group = "Density of Abundance")
  df2 <- data.frame(x = dy$x, y = dy$y, group = "Density of Abundance_hat")
  df <- rbind(df1, df2)
  
  plt <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1, show.legend = TRUE) +  # Adjusted size to linewidth
    labs(title = main, x = xlab, y = "Density") +
    scale_color_manual(values = c(col1, col2)) +
    theme_minimal()
  
  return(plt)
}

plot_mse_barplot <- function(abundance, fitted_abundance,label_order=NULL) {
  # Calcul des MSE pour chaque espèce
  mse_tmp <- rep(NA, ncol(abundance))
  for(j in seq_along(mse_tmp)) {
    mse_tmp[j] <- mse_predict(abundance[, j], fitted_abundance[, j])  
  }
  
  # Créer un data frame pour les données du graphique
  mse_data <- data.frame(Species = colnames(abundance), MSE = mse_tmp)
  
  if(!is.null(label_order)){
    mse_data$Species <- factor(mse_data$Species,
                               levels = label_order)
  }
  
  # Créer le graphique avec ggplot2
  gg <- ggplot(mse_data, aes(x = Species, y = MSE)) +
    geom_bar(stat = "identity", fill = "skyblue", na.rm = TRUE) +
    labs(title = "MSE entre les valeurs réelles et prédites pour chaque espèce",
         x = "Espèces",
         y = "MSE") +
    theme_minimal()
}

plot_residus_density <- function(abundance, fitted_abundance, file_path, verbose = FALSE) {
  # Crée le répertoire s'il n'existe pas
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
  }
  
  # Initialiser la matrice des résidus
  residus <- matrix(NA, nrow = nrow(abundance), ncol = ncol(abundance))
  
  for (j in 1:ncol(abundance)) {
    # Calcul des résidus
    residus[, j] <- abundance[, j] - fitted_abundance[, j]
    
    # Déterminer la largeur de bande
    if (bw.nrd0(residus[, j]) < 0.3) {
      mon_bw <- 0.3
    } else {
      mon_bw <- bw.nrd0(residus[, j])
    }
    
    # Graphique de densité des résidus
    dens <- density(residus[, j], bw = mon_bw)
    file_name <- file.path(file_path, paste("residus_density_species_", j, ".png", sep = ""))
    png(file_name, width = 800, height = 600)
    plot(dens, main = paste("Espèce", j), xlab = "Résidus", ylab = "Densité", col = "blue")
    
    # Ajout des densités pour les autres espèces avec une couleur graduelle
    for (k in 1:ncol(abundance)) {
      if (k != j) {
        dens_tmp <- density(abundance[, k] - fitted_abundance[, k])
        lines(dens_tmp, col = rgb(0, 0, 0, alpha = 0.2))
      }
    }
    dev.off()
    
    # Affiche un message de confirmation si verbose est TRUE
    if (verbose) {
      message("Fichier sauvegardé: ", file_name)
    }
  }
  # Retourner un message de succès ou autre information si nécessaire
  return(invisible(NULL))
}

plot_abundance_vs_environment <- function(B_hat, data , output_dir) {
  selected_vars <- colnames(B_hat[-1,])[apply(B_hat[-1,], 2,
                                              function(x) any(x != 0 & x < max(B_hat[-1,])))]
  selected_vars_and_coeffs <- as.matrix(B_hat[-1, selected_vars])
  colnames(selected_vars_and_coeffs) <- selected_vars

  # Boucle pour parcourir chaque élément de la matrice selected_vars_and_coeffs
  for (i in 1:nrow(selected_vars_and_coeffs)) {
    for (j in 1:ncol(selected_vars_and_coeffs)) {
      row_name <- rownames(selected_vars_and_coeffs)[i]
      col_name <- colnames(selected_vars_and_coeffs)[j]

      # Vérifier si le coefficient est différent de zéro
      if (selected_vars_and_coeffs[i, j] != 0) {
        # Récupérer les données de l'abondance et de l'environnement
        environment <- data[[row_name]]
        abundance <- data$Abundance[, which(colnames(data$Abundance) == col_name)]

        # Créer un data.frame pour ggplot
        plot_data <- data.frame(environment = environment, abundance = abundance)

        # Créer le graphique
        plotcovar_Abund <- ggplot(plot_data, aes(x = environment, y = abundance)) +
          geom_point() +
          labs(x = col_name, y = row_name) +
          theme_minimal()

        # Sauvegarder le graphique
        ggsave(file.path(output_dir, paste("plotcovar_Abund_", row_name, "_", col_name, ".pdf", sep = "")), plotcovar_Abund, width = 10, height = 6, units = "in")      }
    }
  } 
}


