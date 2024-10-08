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


# Fonction pour calculer le log-vraisemblance
calculate_log_likelihood <- function(Abundance, Covariate, B, Sigma, M, S) {
  Z <- as.matrix(Covariate %*% B)
  loglik <- sum(dpois(Abundance, exp(Z + M + S^2 / 2), log = TRUE))
  return(loglik)
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
  #Covariate <- cbind(rep(1,n),Covariate)            # Ajout de l'intercept pour la matrices des covariables 
  
  O <- O #offset PMG
  #dim(O)
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
        # cat("**** QR decomposition***** \n ")
        # QR decomposition de pour calculer l'inverse de l'information de fisher
        qr_decomposition <- qr(fish_Sigma)
        # Matrice Q
        Q <- qr.Q(qr_decomposition)
        
        # Matrice R
        R <- qr.R(qr_decomposition)
        
        # Calcul de l'inverse de R
        inverse_R <- MASS::ginv(R)  #MASS::ginv
        
        # Transposée de Q
        transpose_Q <- t(Q)
        
        # Calcul de l'inverse de l'information de fisher
        inverse_fish_Sigma <-inverse_R %*% transpose_Q
        # cat("*************** ",sum(inverse_fish_Sigma), "\n")
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
        
        #--------------------------------- 
        
        
        
        diff_B <- sum(abs(B_newmat-B_oldf))
        if(diff_B <= B_tol)
        {
          #cat("Convergence atteint avant max_iter","\n")
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
          #cat(" pas de Convergence","\n")
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
      loglik <- calculate_log_likelihood(Abundance, Covariate, B, Sigma, M, S)
    }
    
    model_par <- list(B = B, Sigma = Sigma,Omega = Omega)
    var_par <- list(M= M, S = plngendata$var_par$S, S2 = plngendata$var_par$S2^2)
    
    
    res1 <- list(model_par = model_par, var_par = var_par,vcov_model = plngendata$vcov_model, Offset = O, Covariate = Covariate, dd = dd, loglik = loglik)
    return(res1)
    
  } , error = function(e){ return(res1) }) #Spécification de la classe d'erreur et gestion de l'erreur
  cat("\n\n\t Adjusting a SIC-PLN Model : \n")
  cat("\n\t Execution of SICPLN is DONE! \n\t")
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
SICPLN <- function(formule, offset_column_name = NULL , data , numb_eps = 100, maxIter = 300, control = control) {
  if (!is.null(offset_column_name)) {
    offsett <- calcul_Offset(data, column_name = offset_column_name)
    offset_values <- data[[offset_column_name]]
    
    data <- data[, colnames(data) != offset_column_name]
    
    Covariate <- model.matrix(formule, data)
    Abundance <- model.response(model.frame(formule, data))
    
    capture.output({plngendata <- PLN(Abundance ~ . + offset(offset_values), data = data, control = control)}, file = NULL)
    colnames(offsett) <- colnames(plngendata$model_par$B)
    
    res_fs <- compute_fisher_opt_pln(Covariate = Covariate, Abundance = Abundance, O = offsett, numb_eps = numb_eps, maxIter = maxIter, plngendata = plngendata)
    
    colnames(res_fs$model_par$B) <- colnames(plngendata$model_par$B)
    rownames(res_fs$model_par$B) <- rownames(plngendata$model_par$B)
    
    predictions <- predict_sicpln(Covariate = Covariate, Offset = offsett, params = list(B = res_fs$model_par$B, S = res_fs$var_par$S, M = res_fs$var_par$M))
    predict_pln <- predict_sicpln(Covariate = Covariate, Offset = offsett, params = list(B = plngendata$model_par$B, S = plngendata$var_par$S, M = plngendata$var_par$M))
    residuals <- sum(abs(Abundance - predictions))
    residuals_pln <- sum(abs(Abundance - predict_pln))
    Rsquared <- compute_r_squared(Abundance, predictions)
    
    loglik <- res_fs$loglik
    nb_param <- length(res_fs$model_par$B)
    BIC <- -2 * loglik + log(nrow(data)) * nb_param
    
    resf <- list(res_pln = plngendata, res_fisher = res_fs, fitted = predictions, residuals = residuals, residuals_pln = residuals_pln, R_sqared = Rsquared,
                 loglik = loglik, BIC = BIC, nb_param = nb_param)
  } else {
    offsett <- calcul_Offset(data, column_name = offset_column_name)
    
    Covariate <- model.matrix(formule, data)
    Abundance <- model.response(model.frame(formule, data))
    
    capture.output({plngendata <- PLN(Abundance ~ ., data = data, control = control)}, file = NULL)
    colnames(offsett) <- colnames(plngendata$model_par$B)
    
    res_fs <- compute_fisher_opt_pln(Covariate = Covariate, Abundance = Abundance, O = offsett, numb_eps = numb_eps, maxIter = maxIter, plngendata = plngendata)
    
    colnames(res_fs$model_par$B) <- colnames(plngendata$model_par$B)
    rownames(res_fs$model_par$B) <- rownames(plngendata$model_par$B)
    
    predictions <- predict_sicpln(Covariate = Covariate, Offset = offsett, params = list(B = res_fs$model_par$B,S = res_fs$var_par$S , M = res_fs$var_par$M))
    predict_pln <- predict_sicpln(Covariate = Covariate, Offset = offsett, params = list(B = plngendata$model_par$B, S = plngendata$var_par$S, M = plngendata$var_par$M))
    residuals <- sum(abs(Abundance - predictions))
    residuals_pln <- sum(abs(Abundance - predict_pln))
    Rsquared <- compute_r_squared(Abundance, predictions)
    
    loglik <- res_fs$loglik
    nb_param <- length(res_fs$model_par$B)
    BIC <- -2 * loglik + log(nrow(data)) * nb_param
    
    resf <- list(res_pln = plngendata, res_fisher = res_fs, fitted = predictions, residuals = residuals, residuals_pln = residuals_pln, R_sqared = Rsquared,
                 loglik = loglik, BIC = BIC, nb_param = nb_param)
  }
  
  class(resf) <- "SICPLN"
  return(resf)
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
# Fonction d'impression pour le modèle SICPLN
print.SICPLN <- function(x, ...) {
  cat("A multivariate Poisson Lognormal fit with full covariance model.\n")
  cat("==================================================================\n")
  cat(sprintf("  nb_param   loglik      BIC\n"))
  cat(sprintf("%-10d %-10.3f %-10.3f\n", x$nb_param, x$loglik, x$BIC))
  cat("==================================================================\n")
  cat("  * Useful fields\n")
  cat("  $model_par,  $var_par, res_pln, res_fisher\n")
  cat("  $loglik, $BIC, $nb_param, $criteria\n")
  cat("* Useful S3 methods\n")
  cat("  print(), coef(), sigma(), vcov(), fitted()\n")
  cat("  predict_sicpln() \n")
}

# Méthode coef pour la classe SICPLN
coef.SICPLN <- function(object, ...) {
  return(object$res_fisher$model_par$B)
}

# Méthode sigma pour la classe SICPLN
sigma.SICPLN <- function(object, ...){
  return(object$res_fisher$var_par$S)
}

# Méthode vcov pour la classe SICPLN
vcov.SICPLN <- function(object, ...){
  return(object$res_fisher$vcov_model)
}

# Méthode fitted pour la classe SICPLN
fitted.SICPLN <- function(object, ...){
  return(object$fitted)
}

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
          chi2_test <- suppressWarnings(chisq.test(table(y, Covariate[, j])))
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

