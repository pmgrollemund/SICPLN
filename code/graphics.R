 # Définition de la fonction pour tracer les coefficients : ----
plot_coefficient <- function(coef_matrix, axis_names = c("Genus", "Variables", "Coefficient Value", "Title")) {
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
    scale_fill_gradient2(low = "darkred",high = "darkgreen",mid = "white", midpoint = 0, name = color_gradient_name) +
    coord_equal() +
    labs(x = x_axis_name, y = y_axis_name, title = title, color = "black") +
    theme(legend.position = "right", plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 20), legend.text = element_text(size = 15), strip.text.x = element_text(size = 20, colour = "black"), axis.text.y = element_text(size = 20, colour = "black"), axis.text.x = element_text(size = 20, colour = "black"))
  
  # Retourner le graphique
  return(plot)
}
#----


# Détecter les lignes de séparation pour les variables catégorielles
factor_lines <- function(data, coef_matrix) {
  # Initialiser le vecteur index
  index <- vector()
  
  # Obtenir les noms des lignes de la matrice de coefficients, en excluant la première ligne
  coef_rownames <- rownames(coef_matrix)  #[-1]
  
  # Obtenir les noms des colonnes du data frame, en excluant la première colonne
  data_colnames <- colnames(data)  #[-1]
  
  # Parcourir les noms des lignes de la matrice de coefficients
  for (i in seq_along(coef_rownames)) {
    # Si le nom de la ligne n'apparaît pas dans les noms de colonnes des données
    if (!(coef_rownames[i] %in% data_colnames)) {
      # Ajouter l'index à la liste
      index <- c(index, i)
    }
  }
  
  # Retourner le vecteur d'index
  return(index)
}


# Définition de la fonction pour tracer les coefficients :----
coef_genus_barre <- function(coef_matrix, axis_names = c("Genus", "Variables", "Coefficient Value", "Title"), data) {
  # Renommer les axes et le gradient de couleur à partir des paramètres d'entrée
  y_axis_name <- axis_names[2]
  x_axis_name <- axis_names[1]
  color_gradient_name <- axis_names[3]
  title <- axis_names[4]
  
  # Convertir la matrice de coefficients en un format long pour ggplot
  melted_coef_matrix <- melt(coef_matrix)
  melted_coef_matrix$sparsity <- ifelse(melted_coef_matrix$value == 0, 'Zero', 'Nonzero')
  
  index <- factor_lines(data, coef_matrix)
  
  # Création du graphique avec ggplot
  plot <- ggplot(melted_coef_matrix, aes(x = as.factor(Var2), y = as.factor(Var1), pattern = sparsity, fill = value)) +
    geom_tile(aes(fill = value), colour = "white") + 
    geom_tile_pattern(pattern_color = NA, pattern_fill = "black", pattern_angle = 45,
                      pattern_density = 0.15, pattern_spacing = 0.065, pattern_key_scale_factor = 1) +
    scale_pattern_manual(values = c(Zero = "circle", Nonzero = "none"), name = "") +
    scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "white", midpoint = 0, name = color_gradient_name) +
    coord_equal() + 
    labs(x = x_axis_name, y = y_axis_name, title = title, color = "black") +
    guides(pattern = guide_legend(override.aes = list(fill = "white"))) +
    theme(legend.position = "right", plot.title = element_text(size = 10, hjust = 0.5), 
          axis.title = element_text(size = 10), legend.text = element_text(size = 10), 
          strip.text.x = element_text(size = 10, colour = "black"), 
          axis.text.y = element_text(size = 10, colour = "black"), 
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, colour = "black")) # Ajuster l'angle du texte sur l'axe x
  
  # Ajouter des lignes horizontales pour séparer les variables
  if (length(index) > 0) {
    plot <- plot + 
      geom_hline(yintercept = min(index) - 0.5, color = "black") +
      geom_hline(yintercept = max(index) + 0.5, color = "black")
  }
  
  # Retourner le graphique
  return(plot)
}
#----
# MAtrice de precision
plot_precision_matrix <- function(coef_matrix, 
                                  axis_names = c("Genus","Variables","Coefficients value","SICPLN"),
                                  limit = c(-1,1)) {
  # Renommer les axes et le gradient de couleur à partir des paramètres d'entrée
  y_axis_name <- axis_names[2]
  x_axis_name <- axis_names[1]
  color_gradient_name <- axis_names[3]
  title <- axis_names[4]
  
  # Convertir la matrice de coefficients en un format long pour ggplot
  melted_coef_matrix <- melt(coef_matrix)
  
  # Création du graphique avec ggplot
  plot <- ggplot(melted_coef_matrix, aes(Var1, Var2, value)) +
    geom_tile(aes(fill = value), colour = "white") +
    scale_fill_gradient2(low = "darkred", high = "darkgreen", mid = "white", midpoint = 0, 
                         limit = limit, name = color_gradient_name) +
    labs(x = x_axis_name, y = y_axis_name, title = title, color = "black") +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 20, hjust = 0.5),
      axis.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      strip.text.x = element_text(size = 15, colour = "black"),
      axis.text.y = element_text(size = 10, colour = "black"),  # Taille réduite pour l'axe y
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5, colour = "black")  # Texte incliné à 90°
    )
  
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

plot_abundance_vs_environment <- function(B_hat, data ) {
  selected_vars <- colnames(B_hat[-1,])[apply(B_hat[-1,], 2,
                                              function(x) any(x != 0 & x < max(B_hat[-1,])))]
  
  if(length(selected_vars) > 0){
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
            geom_point() + geom_smooth(method = lm ,formula = y ~ x)
          labs(x = col_name, y = row_name) +
            theme_minimal()
          
          # Sauvegarder le graphique
          #ggsave(file.path(output_dir, paste("plotcovar_Abund_", row_name, "_", col_name, ".pdf", sep = "")), plotcovar_Abund, width = 10, height = 6, units = "in")      }
      }
        
    }
  }
   return(plotcovar_Abund)
} 
}

plot_combined_density <- function(data1, data2, data3, data4, col1 = "red", col2 = "blue", col3 = "green", col4 = "black" , main = "Combined Density Plot", xlab) {
  dx1 <- density(data1)
  dy1 <- density(data2)
  dy2 <- density(data3)
  dy3 <- density(data4)
  
  
  df1 <- data.frame(x = dx1$x, y = dx1$y, group = "Density of Abundance (PLN)")
  df2 <- data.frame(x = dy1$x, y = dy1$y, group = "Density of Abundance_hat (PLN)")
  df3 <- data.frame(x = dy2$x, y = dy2$y, group = "Density of Abundance_hat (SICPLN)")
  df4 <- data.frame(x = dy3$x, y = dy3$y, group = "Density of Abundance_hat (GLMNET)")
  
  
  df <- rbind(df1, df2, df3, df4)
  
  plt <- ggplot(df, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1, show.legend = TRUE) +
    labs(title = main, x = xlab, y = "Density") +
    scale_color_manual(values = c(col1, col2, col3, col4)) +
    theme_minimal()
  
  return(plt)
}


# Fonction pour créer les graphiques de régression et les sauvegarder
plot_abundance_regressions <- function(data, file_path) {
  # Crée le répertoire s'il n'existe pas
  if (!dir.exists(file_path)) {
    dir.create(file_path, recursive = TRUE)
  }
  
  # Convertir la matrice d'abondance en format long
  abundance_long <- as.data.frame(as.table(as.matrix(data %>%
                                                       select(starts_with("Abundance")))))
  colnames(abundance_long) <- c("Observation", "Species", "Abundance")
  
  # Convertir les variables explicatives en format long
  explanatory_long <- data %>%
    select(-starts_with("Abundance")) %>%
    mutate(Observation = as.character(row.names(data))) %>%
    pivot_longer(cols = -Observation, names_to = "Variable", values_to = "Value")
  
  # Joindre les données d'abondance avec les variables explicatives
  suppressWarnings(combined_data <- abundance_long %>%
                     left_join(explanatory_long, by = "Observation"))
  
  # Créer une liste pour stocker les graphiques
  plots <- list()
  
  # Itérer sur chaque variable d'abondance
  for (species in unique(abundance_long$Species)) {
    # Filtrer les données pour l'espèce actuelle
    species_data <- combined_data %>%
      filter(Species == species)
    
    # Créer un graphique de régression pour l'espèce actuelle
    p <- suppressMessages(ggplot(species_data, aes(x = Value, y = Abundance)) +
                            geom_point() +
                            geom_smooth(method = "lm", se = FALSE, color = "blue") +
                            facet_wrap(~ Variable, scales = "free") +
                            labs(x = "Variables explicatives", y = species, 
                                 title = paste("Régression pour", species)) +
                            theme_minimal() + theme(panel.background = element_rect(fill = "white", color = NA),
                                                    plot.background = element_rect(fill = "white", color = NA),
                                                    panel.grid.major = element_line(color = "grey80"),
                                                    panel.grid.minor = element_line(color = "grey90"),
                                                    text = element_text(size = 12, color = "black"),
                                                    axis.text = element_text(color = "black"),
                                                    plot.title = element_text(hjust = 0.5, face = "bold")))
    
    # Ajouter le graphique à la liste
    plots[[species]] <- p
    
    # Sauvegarder le graphique
    file_name <- file.path(file_path, paste("regression_plot_", species, ".png", sep = ""))
    suppressMessages(ggsave(file_name, plot = p, width = 8, height = 6))
  }
  
  return(plots)
}

create_heatmap_with_groups <- function(data, num_groups_col = 3) {
  # Vérifier que les données sont une matrice
  if (!is.matrix(data$Abundance)) {
    stop("La donnée 'Abundance' doit être une matrice.")
  }
  
  # Palette de couleurs pour le heatmap
  colMain <- viridis::viridis(100)  # 100 couleurs pour une bonne gradation
  
  # Calculer les distances et créer le dendrogramme pour les colonnes uniquement
  col_dist <- dist(t(data$Abundance))
  col_hclust <- hclust(col_dist)
  
  # Convertir l'objet hclust en dendrogramme
  col_dendro <- as.dendrogram(col_hclust)
  
  # Découper le dendrogramme pour former des groupes
  col_groups <- cutree(col_hclust, k = num_groups_col)
  
  # Créer des couleurs pour les groupes
  col_colors <- viridis::viridis(num_groups_col)[col_groups]
  
  # Créer le heatmap sans dendrogramme pour les lignes
  graph <- suppressWarnings(heatmap.2(
    data$Abundance,                # Les données de la matrice de chaleur
    Colv = col_dendro,             # Ajouter le dendrogramme des colonnes
    Rowv = FALSE,                  # Retirer le dendrogramme des lignes
    scale = "row",                 # Normaliser les données par ligne
    col = colMain,                 # Palette de couleurs pour le heatmap
    trace = "none",                # Désactiver le tracé des lignes de trace
    margins = c(12, 12),           # Marges pour les étiquettes des axes
    key = TRUE,                    # Afficher la légende de la couleur
    keysize = 1.5,                 # Taille de la légende de la couleur
    key.title = "Intensité",       # Titre de la légende de la couleur
    key.xlab = "Valeurs",          # Étiquette de l'axe des x de la légende
    key.ylab = "",                 # Étiquette de l'axe des y de la légende (vide)
    cexRow = 1.0,                  # Taille du texte des lignes
    cexCol = 1.0,                  # Taille du texte des colonnes
    main = "Heatmap des Abondances",# Titre du graphique
    ColSideColors = col_colors,    # Couleurs des barres latérales des colonnes
    labRow = NULL,                 # Suppression des noms des lignes si trop nombreux
    labCol = NULL                  # Suppression des noms des colonnes si trop nombreux
  ))
  
  return(graph)
}


# Fonction pour générer un graphique ACP amélioré
create_pca_plot <- function(covarnumeric_data, covarfactor_data = NULL, title = "Analyse en Composantes Principales (ACP)") {
  # Réaliser l'ACP sur les covariates
  acp_result <- prcomp(covarnumeric_data, scale. = TRUE)
  
  # Visualiser l'ACP avec une présentation améliorée
  pca_plot <- fviz_pca_biplot(acp_result,
                              geom.ind = "point",       # Représenter les individus par des points
                              geom.var = c("arrow", "text"),       # Représenter les variables par des flèches
                              label = "var",            # Afficher les étiquettes des variables
                              addEllipses = TRUE,       # Ajouter des ellipses de confiance
                              ellipse.level = 0.95,     # Niveau de confiance à 95%
                              habillage = covarfactor_data,
                              palette = "jco",          # Palette de couleurs
                              arrowsize = 1.2,          # Taille des flèches
                              pointshape = 21,          # Forme des points
                              pointsize = 2,            # Taille des points
                              fill.ind = "black",       # Couleur de remplissage des points
                              col.var = "blue",         # Couleur des flèches des variables
                              col.ind = "black",        # Couleur des points des individus
                              repel = TRUE,             # Repousser les étiquettes pour éviter le chevauchement
                              legend.title = list(color = "Groupes", shape = "Individus"),
                              title = title) +          # Titre du graphique
    theme_minimal(base_size = 14) +         # Utiliser un thème minimaliste
    theme(axis.title = element_text(size = 14, face = "bold"), # Taille et style des titres d'axes
          axis.text = element_text(size = 12),  # Taille du texte des axes
          plot.title = element_text(hjust = 0.5, face = "bold"), # Centrer et mettre en gras le titre
          legend.position = "top",        # Position de la légende
          legend.text = element_text(size = 12))  # Taille du texte de la légende
  
  return(pca_plot)
}
