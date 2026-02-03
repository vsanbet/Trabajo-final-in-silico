# ---------------- Funciones auxiliares ----------------

get_protein_ids <- function(df) {
  if ("Majority.protein.IDs" %in% colnames(df)) {
    return(df$Majority.protein.IDs)
  } else {
    return(rownames(df))
  }
}

# -------------- Función principal ------------------

run_geyer_pipeline <- function(
    path,
    group_regex,
    group_names,
    analysis_title,
    mod_file = NULL
) {
  
  
  # ---------------- Libraries --------------------
  library(dplyr)
  library(ggplot2)
  library(tidyverse)

  
  # ---------------- Output directories ----------------
  plot_dir <- file.path(path, "plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  
  plot_prefix <- gsub("[^A-Za-z0-9]+", "_", analysis_title)
  
  
  message("Leyendo archivos...")
  
  # ---------------- QC ----------------
  summary <- read.delim(
  file.path(path, "summary.txt"),
  check.names = FALSE
)
  
  
  # ---------------- proteinGroups ----------------
  pg <- read.delim(
  file.path(path, "proteinGroups.txt"),
  check.names = FALSE
)
  
  pg_clean <- pg %>%
    dplyr::filter(
      Reverse != "+",
      `Potential contaminant` != "+",
      `Only identified by site` != "+"
    )
  
  # ---------------- Selección de columnas de cuantificación ----------------
  lfq_cols <- grep("^LFQ intensity", colnames(pg_clean), value = TRUE)
  intensity_cols <- grep("^Intensity", colnames(pg_clean), value = TRUE)
  
  if (length(lfq_cols) > 0) {
    expr <- pg_clean[, lfq_cols]
    colnames(expr) <- sub("^LFQ intensity ", "", colnames(expr))
  } else if (length(intensity_cols) > 0) {
    expr <- pg_clean[, intensity_cols]
    colnames(expr) <- sub("^Intensity ", "", colnames(expr))
  } else {
    stop("No se encontraron columnas LFQ ni Intensity para cuantificación")
  }
  
  rownames(expr) <- pg_clean$`Protein IDs`
  
  expr[expr == 0] <- NA
  
  keep <- rowSums(!is.na(expr)) >= 0.7 * ncol(expr)
  if (sum(keep) == 0) stop("No hay filas con suficientes datos para PCA.")
  
  expr_log <- log2(expr[keep, ])
  expr_for_bias <- expr_log   # SIN centrar
  expr_norm <- sweep(expr_log, 2, apply(expr_log, 2, median, na.rm = TRUE), FUN = "-")
  
  expr_noimp <- as.matrix(expr_norm)
  mode(expr_noimp) <- "numeric"
  
  # ---------------- PCA ----------------
  expr_pca <- expr_noimp[complete.cases(expr_noimp), ]
  if (nrow(expr_pca) < 2) stop("No hay suficientes filas completas para PCA.")
  
  pca <- prcomp(t(expr_pca), scale. = FALSE)
  
  sample_info <- data.frame(Sample = colnames(expr_pca), stringsAsFactors = FALSE)
  sample_info$Group <- "Other"
  
  for (i in seq_along(group_regex)) {
    sample_info$Group[
      grepl(group_regex[i], sample_info$Sample, ignore.case = TRUE)
    ] <- group_names[i]
  }
  
  message("Grupos asignados:")
  print(table(sample_info$Group))
  
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  ) %>% left_join(sample_info, by = "Sample")
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4, alpha = 0.85) +
    stat_ellipse(aes(group = Group), linewidth = 0.8, linetype = "dashed") +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = analysis_title,
      x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
      y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
    ) +
    scale_color_brewer(palette = "Set1")
  
  print(p_pca)
  
  ggsave(
    filename = file.path(plot_dir, paste0(plot_prefix, "_PCA.png")),
    plot = p_pca, width = 7, height = 6, dpi = 300
  )
  
  
  # ---------------- Oxidation (M) ----------------
  if (!is.null(mod_file)) {
    
    mods <- read.delim(
      file.path(path, mod_file),
      check.names = FALSE
    )

    
    intensity_cols <- grep("^Intensity ", colnames(mods), value = TRUE)
    if (length(intensity_cols) == 0)
      stop("No se encontraron columnas de intensidad")
    
    if ("Number of Oxidation (M)" %in% colnames(mods)) {
      ox_rows <- mods[["Number of Oxidation (M)"]] > 0
    } else if ("Oxidation (M) Probabilities" %in% colnames(mods)) {
      ox_rows <- !is.na(mods[["Oxidation (M) Probabilities"]])
    } else if ("Oxidation (M)" %in% colnames(mods)) {
      ox_rows <- mods[["Oxidation (M)"]] > 0
    } else {
      stop("No se encontró columna de Oxidation (M)")
    }
    
    ox_score <- colSums(mods[ox_rows, intensity_cols, drop = FALSE], na.rm = TRUE)
    names(ox_score) <- sub("^Intensity ", "", names(ox_score))
    
    ox_df <- data.frame(
      Sample = names(ox_score),
      Oxidation_M = ox_score
    ) %>% left_join(sample_info, by = "Sample")
    
    p_ox <- ggplot(ox_df, aes(Group, Oxidation_M, fill = Group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = paste("Oxidation (M) –", basename(mod_file)),
        y = "Sum of intensities"
      ) +
      scale_fill_brewer(palette = "Set1")
    
    print(p_ox)
    ggsave(
      filename = file.path(plot_dir, paste0(plot_prefix, "_Oxidation_M.png")),
      plot = p_ox, width = 7, height = 6, dpi = 300
    )
    
  }
  
  # ---------------- Missing values per sample ----------------
  
  missing_df <- data.frame(
    Sample = colnames(expr),
    Missing_fraction = colMeans(is.na(expr))
  ) %>% left_join(sample_info, by = "Sample")
  
  p_missing_density <- ggplot(missing_df,
                              aes(Missing_fraction, color = Group, fill = Group)) +
    geom_density(alpha = 0.3) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("Distribution of missing values –", analysis_title),
      x = "Fraction of missing proteins",
      y = "Density"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1")
  
  print(p_missing_density)
  
  
  ggsave(
    filename = file.path(plot_dir, paste0(plot_prefix, "_MissingValues.png")),
    plot = p_missing_density, width = 7, height = 6, dpi = 300
  )
  
  # ---------------- Coefficient of Variation (CV) ----------------
  
  cv_list <- lapply(unique(sample_info$Group), function(g) {
    
    samples_g <- sample_info$Sample[sample_info$Group == g]
    samples_g <- intersect(samples_g, colnames(expr_log))
    
    if (length(samples_g) < 2) return(NULL)
    
    submat <- expr_log[, samples_g, drop = FALSE]
    
    cv <- apply(submat, 1, function(x) {
      if (all(is.na(x))) return(NA)
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    })
    
    data.frame(
      Group = g,
      CV = cv
    )
  })
  
  cv_df <- bind_rows(cv_list)
  
  p_cv_ecdf <- ggplot(cv_df, aes(CV, color = Group)) +
    stat_ecdf(size = 1) +
    coord_cartesian(xlim = c(0, quantile(cv_df$CV, 0.95, na.rm = TRUE))) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste("ECDF of CV –", analysis_title),
      x = "CV (sd / mean)",
      y = "Cumulative fraction"
    ) +
    scale_color_brewer(palette = "Set1")
  
  
  print(p_cv_ecdf)
  ggsave(
    filename = file.path(plot_dir, paste0(plot_prefix, "_CV.png")),
    plot = p_cv_ecdf, width = 7, height = 6, dpi = 300
  )
  
}

# ===============================
# Ejemplo de uso de run_geyer_pipeline
# ===============================
print(
  'run_geyer_pipeline( 
      path = "ruta/a/tu/dataset", 
      group_regex = c("Grupo1", "Grupo2", ...), 
      group_names = c("Grupo1", "Grupo2", ...), 
      analysis_title = "Título del análisis", 
      mod_file = "Oxidation (M)Sites.txt"'
     )



