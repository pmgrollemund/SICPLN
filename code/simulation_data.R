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
sparsity_recognition <- function(vrai_beta, beta_estime){
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
  sp <- nombre_espece
  n <- taille_echant
  
  nb_times <- nbre_repetition       # Nombre de répétitions pour la simulation
  
  
  n <- n  # Nombre d'observations
  p <- p  # Nombre de covariables
  d <- d # Nombre d'espèces
  
  # Définition des paramètres pour le calcul des erreurs 
  beta <- matrix(1:(p*d),nrow = d,ncol = p)
  Sigma <- diag(p)
  
  # Génération des covariables (matrice aléatoire)
  Covariate <- matrix(runif(n * d), nrow = n, ncol = d)
  
  # Génération des données d'exemple
  gen_data <- simul_Abundance(B = beta, Sigma, Covariate)
  
  #================================================================================
  rownames(gen_data$Covariate) <- rownames(gen_data$Abundance)
  
  gen_data$Covariate <- cbind(gen_data$Covariate, o = gen_data$o)
  
  
  data_tmp <- list(Abundance = gen_data$Abundance, Covariate = gen_data$Covariate)
  data_to_PLN <- prepare_data(data_tmp$Abundance, data_tmp$Covariate)
  data_to_PLN <- data_to_PLN[, -ncol(data_to_PLN)]
  data_to_PLN <- data_to_PLN[, -ncol(data_to_PLN)]
  
  formule <- Abundance ~ .
  data <- data_to_PLN
  
  # Spécifiez les paramètres de contrôle
  control <- PLN_param(backend = "nlopt", covariance = "full")
  
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
  # Initialisation des matrices pour stocker les coefficients estimés
  list_mat_B_sp <- list(SIC_PLN = matrix(NA, 1, sp + 2),
                        PLN = matrix(NA, 1, sp + 2),
                        GLMNET = matrix(NA, 1, sp + 2))
  
  
  
  res_SICPLN <- SICPLN(formule ,offset_column_name = NULL, data = data_to_PLN, control = control)
  
  res_glmpois <- glmnet_lasso_poisson(data_tmp$Abundance,data_tmp$Covariate[, -ncol(data_tmp$Covariate)])
  
  # Boucle sur le nombre de répétitions
  for (t in 1:nb_times) {  
    # Code pour chaque répétition de simulation
    
    cat("Preprocissing for time(s) :", t , "\n")
    
    # Calcul des erreurs pour chaque méthode
    err_sic_sp[t] <- norm(((beta)-(res_SICPLN$res_fisher$model_par$B[-1,])),"F")/norm(beta,"F")
    err_pln_sp[t] <- norm(((beta)-(res_SICPLN$res_pln$model_par$B[-1,])),"F")/norm(beta,"F")
    err_glm_sp[t] <- norm(((beta)-(res_glmpois$hat_B[-1,])),"F")/norm(beta,"F")
    
    # Calcul des MSE
    # Prediction SIC_PLN
    # params <-params1
    data_forma_pred <- list(Covariate=gen_data$Covariate,Abundance=gen_data$Abundance
                            ,O=res_SICPLN$res_fisher$Offset)
    data <- data_forma_pred
    params1 <- list(B=res_SICPLN$res_fisher$model_par$B,
                    S=res_SICPLN$res_fisher$var_par$S, M=res_SICPLN$res_fisher$var_par$M)
    n <- nrow(res_SICPLN$res_fisher$var_par$M)
    pred_sic <- Error_var(data=data_forma_pred, params=params1)  #Couche latente Z avec SIC
    
    # Prediction PLN
    # data_forma_pred$Y
    params2 <- list(B=res_SICPLN$res_pln$model_par$B,S=res_SICPLN$res_pln$var_par$S
                    ,M=res_SICPLN$res_pln$var_par$M)
    pred_pln <- Error_var(data=data_forma_pred, params=params2)  #couche latente Z avec PLNmodel 
    MSE_SIC[t] <- mse_predict(gen_data$Abundance,pred_sic)  #Erreur genere en utilisant SIC
    MSE_PLN[t] <- mse_predict(gen_data$Abundance,pred_pln)  #Erreur genere en utilisant PLNmodel
    MSE_GLMNET[t] <- mse_predict(gen_data$Abundance,res_glmpois$prediction)  #Erreur genere en utilisant GLMNET (la fonction calcule ses valeurs predicts)
    
    
    # Calcul des taux de faux positifs et faux négatifs
    SIC_sparsity <- sparsity_recognition(beta,res_SICPLN$res_fisher$model_par$B[-1,])
    GLMNET_sparsity <- sparsity_recognition(beta,res_glmpois$hat_B[-1,])
    PLN_sparsity <- sparsity_recognition(beta,res_SICPLN$res_pln$model_par$B[-1,])
    TPR_SIC[t] <- SIC_sparsity$TPR
    TPR_PLN[t] <- PLN_sparsity$TPR
    TPR_GLMNET[t] <- GLMNET_sparsity$TPR
    TNR_SIC[t] <- SIC_sparsity$TNR
    TNR_PLN[t] <- PLN_sparsity$TNR
    TNR_GLMNET[t] <- GLMNET_sparsity$TNR
    
    # Stockage des résultats dans les matrices
    
    # Stockage des résultats dans les matrices
    list_mat_B_sp[["SIC_PLN"]] <- rbind(list_mat_B_sp[["SIC_PLN"]], cbind(abs(beta - res_SICPLN$res_fisher$model_par$B[-1,]), rep(n, d), rep(sp, d)))
    list_mat_B_sp[["PLN"]] <- rbind(list_mat_B_sp[["PLN"]], cbind(abs(beta - res_SICPLN$res_pln$model_par$B[-1,]), rep(n, d), rep(sp, d)))
    list_mat_B_sp[["GLMNET"]] <- rbind(list_mat_B_sp[["GLMNET"]], cbind(abs(beta - res_glmpois$hat_B[-1,]), rep(n, d), rep(sp, d)))
  }
  
  
  df_sp_size=data.frame(nom_fb=c(err_sic_sp,err_pln_sp, err_glm_sp),
                        lab_methode = as.factor(rep(c("SIC_PLN","PLN","GLMNET"),each=length(err_sic_sp))),
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
  
  # Création des structures de données pour les résultats
  SICPLN_matb <- as.matrix(list_mat_B_sp$SIC_PLN)
  PLN_matb <- as.matrix(list_mat_B_sp$PLN)
  GLM_matb <- as.matrix(list_mat_B_sp$GLMNET)
  head(list_mat_B_sp$SIC_PLN)
  GLM_bpt <- as.numeric(GLM_matb[,1:sp])
  SIC_bpt <- as.numeric(SICPLN_matb[,1:sp])
  PLN_bpt <- as.numeric(PLN_matb[,1:sp])
  lab_methode <- rep(c("GLMNET","SIC_PLN","PLN"),each=length(PLN_bpt))
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
