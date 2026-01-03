# Script completo: Analisi e Model Selection tra Normale, t-Student, Laplace
# Per OGNI colonna di dati MEMS/VN100:
# 1. Rimuove outlier (IQR method)
# 2. Centra i dati
# 3. Calcola statistiche (sd, skewness, kurtosis)
# 4. Fitta 3 distribuzioni: Normale, t-Student, Laplace
# 5. Model Selection via AIC/BIC + Anderson-Darling
# 6. Salva conclusioni oggettive (quale distribuzione è più appropriata)

# Installa/carica pacchetti necessari
packages <- c("moments", "MASS", "nortest")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Impostiamo la working directory
script_dir <- getSrcDirectory(function() {})[1]
if (script_dir != "") {
  setwd(script_dir)
}

cartella <- "Euroc12025_separati/"
files <- list.files(path = cartella, pattern = "*.csv", full.names = TRUE)

output_dir <- file.path(cartella, "Risultati")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

colonne <- c("accelerationX", "accelerationY", "accelerationZ",
             "angularSpeedX", "angularSpeedY", "angularSpeedZ")

# ===== FUNZIONI AUSILIARIE =====

# Fit Laplace (location + scale)
fitdistr_laplace <- function(x) {
  mu <- median(x)
  scale <- mean(abs(x - mu))
  return(list(mu = mu, scale = scale))
}

# Log-likelihood distribuzioni
loglik_normal <- function(x, mu, sigma) {
  sum(dnorm(x, mu, sigma, log = TRUE))
}

loglik_t <- function(x, mu, sigma, df) {
  sum(dt((x - mu) / sigma, df, log = TRUE) - log(sigma))
}

loglik_laplace <- function(x, mu, scale) {
  sum(log(1 / (2 * scale)) - abs(x - mu) / scale)
}

# Funzione principale di confronto distribuzioni
compare_distributions <- function(x_clean, name_col, output_file = NULL) {
  cat("\n===== CONFRONTO DISTRIBUZIONI:", name_col, "=====\n")
  
  x_clean <- x_clean[!is.na(x_clean)]
  n <- length(x_clean)
  
  # Data frame per risultati
  results <- data.frame(
    Distribution = c("Normale", "t-Student", "Laplace"),
    n_params = c(2, 3, 2),
    LogLik = NA,
    AIC = NA,
    BIC = NA,
    AD_pvalue = NA,
    stringsAsFactors = FALSE
  )
  
  # 1. FIT NORMALE
  tryCatch({
    fit_norm <- fitdistr(x_clean, "normal")
    mu_n <- fit_norm$estimate["mean"]
    sigma_n <- fit_norm$estimate["sd"]
    ll_norm <- loglik_normal(x_clean, mu_n, sigma_n)
    results[1, "LogLik"] <- ll_norm
    results[1, "AIC"] <- -2 * ll_norm + 2 * 2
    results[1, "BIC"] <- -2 * ll_norm + 2 * log(n)
    
    ad_norm <- ad.test(x_clean)
    results[1, "AD_pvalue"] <- ad_norm$p.value
  }, error = function(e) {
    cat("Errore in fit normale:", e$message, "\n")
  })
  
  # 2. FIT t-STUDENT
  tryCatch({
    fit_t <- fitdistr(x_clean, "t")
    mu_t <- fit_t$estimate["m"]
    sigma_t <- fit_t$estimate["s"]
    df_t <- fit_t$estimate["df"]
    ll_t <- loglik_t(x_clean, mu_t, sigma_t, df_t)
    results[2, "LogLik"] <- ll_t
    results[2, "AIC"] <- -2 * ll_t + 2 * 3
    results[2, "BIC"] <- -2 * ll_t + 3 * log(n)
    
    x_std <- (x_clean - mu_t) / sigma_t
    ad_t <- ad.test(x_std)
    results[2, "AD_pvalue"] <- ad_t$p.value
  }, error = function(e) {
    cat("Errore in fit t-Student:", e$message, "\n")
  })
  
  # 3. FIT LAPLACE
  tryCatch({
    fit_lap <- fitdistr_laplace(x_clean)
    mu_l <- fit_lap$mu
    scale_l <- fit_lap$scale
    ll_lap <- loglik_laplace(x_clean, mu_l, scale_l)
    results[3, "LogLik"] <- ll_lap
    results[3, "AIC"] <- -2 * ll_lap + 2 * 2
    results[3, "BIC"] <- -2 * ll_lap + 2 * log(n)
    
    x_lap_std <- (x_clean - mu_l) / scale_l
    ad_lap <- ad.test(x_lap_std)
    results[3, "AD_pvalue"] <- ad_lap$p.value
  }, error = function(e) {
    cat("Errore in fit Laplace:", e$message, "\n")
  })
  
  # Ranking
  results$AIC_rank <- rank(results$AIC)
  results$BIC_rank <- rank(results$BIC)
  
  # Stampa a video
  print(results)
  
  # Calcoli per interpretazione
  best_aic <- results[which.min(results$AIC), "Distribution"]
  delta_aic_norm_t <- results[1, "AIC"] - results[2, "AIC"]
  delta_aic_norm_lap <- results[1, "AIC"] - results[3, "AIC"]
  
  cat("\n--- RANKING AIC ---\n")
  cat("1.", results[results$AIC_rank == 1, "Distribution"], 
      "(AIC =", round(min(results$AIC, na.rm = T), 2), ")\n")
  cat("2.", results[results$AIC_rank == 2, "Distribution"], "\n")
  cat("3.", results[results$AIC_rank == 3, "Distribution"], "\n")
  
  cat("\n--- DIFFERENZE AIC ---\n")
  cat("ΔAIC (Normale - t-Student) =", round(delta_aic_norm_t, 2), "\n")
  cat("ΔAIC (Normale - Laplace) =", round(delta_aic_norm_lap, 2), "\n")
  cat("Regola: |ΔAIC| < 2 = ambiguo; |ΔAIC| > 10 = chiara preferenza\n")
  
  cat("\n--- ANDERSON-DARLING TEST (p-value > 0.05 = non respinto) ---\n")
  cat("Normale:", round(results[1, "AD_pvalue"], 4), "\n")
  cat("t-Student:", round(results[2, "AD_pvalue"], 4), "\n")
  cat("Laplace:", round(results[3, "AD_pvalue"], 4), "\n")
  
  # CONCLUSIONE FINALE
  cat("\n=== CONCLUSIONE FINALE ===\n")
  
  conclusion <- ""
  
  # Controlli NA per evitare errori
  norm_aic <- results[1, "AIC"]
  t_aic <- results[2, "AIC"]
  lap_aic <- results[3, "AIC"]
  norm_pval <- results[1, "AD_pvalue"]
  
  # Verifica che almeno uno dei modelli sia stato fittato
  valid_aics <- !is.na(c(norm_aic, t_aic, lap_aic))
  
  if (!any(valid_aics)) {
    conclusion <- "ERRORE: nessun modello è stato fittato correttamente"
  } else {
    # Rimuovi NA per il ranking
    best_idx <- which.min(results$AIC)
    best_aic <- results[best_idx, "Distribution"]
    
    # Calcola differenze solo se validi
    delta_aic_norm_t <- if (!is.na(norm_aic) & !is.na(t_aic)) norm_aic - t_aic else NA
    delta_aic_norm_lap <- if (!is.na(norm_aic) & !is.na(lap_aic)) norm_aic - lap_aic else NA
    
    if (!is.na(norm_pval) & norm_pval > 0.05 & best_aic == "Normale") {
      conclusion <- "GAUSSIANO (Normale compatibile, AD test p > 0.05)"
    } else if (!is.na(norm_aic) & !is.na(t_aic) & best_aic == "Normale" & 
               !is.na(delta_aic_norm_t) & abs(delta_aic_norm_t) < 2) {
      conclusion <- "GAUSSIANO QUANTIZZATO (AIC ambiguo, ma Normale è preferita)"
    } else if (best_aic == "t-Student" & !is.na(delta_aic_norm_t) & delta_aic_norm_t < -2) {
      conclusion <- "t-STUDENT (AIC significativamente migliore di Normale; code pesanti)"
    } else if (best_aic == "Laplace") {
      conclusion <- "LAPLACE (doppia esponenziale: picco alto, code pesanti)"
    } else {
      conclusion <- "AMBIGUO (quantizzazione dominante o modelli quasi equivalenti)"
    }
  }
  
  cat(conclusion, "\n")
  
  # Salva su file se specificato
  if (!is.null(output_file)) {
    sink(output_file, append = TRUE)
    cat("\n===== MODEL SELECTION:", name_col, "=====\n")
    print(results)
    cat("\nCONCLUSIONE:", conclusion, "\n")
    cat("AIC ranking: 1.", results[results$AIC_rank == 1, "Distribution"],
        "2.", results[results$AIC_rank == 2, "Distribution"],
        "3.", results[results$AIC_rank == 3, "Distribution"], "\n")
    sink()
  }
  
  return(list(results = results, conclusion = conclusion))
}

# ===== CICLO PRINCIPALE =====

for (file in files) {
  cat("\n\n##################################################################\n")
  cat("ELABORAZIONE FILE:", basename(file), "\n")
  cat("##################################################################\n")
  
  dati <- read.csv(file, header = TRUE)
  dati_selezionati <- dati[, colonne]
  dati_selezionati <- tail(dati_selezionati, min(5000, nrow(dati_selezionati)))
  
  for (colonna in colonne) {
    x_raw <- dati_selezionati[[colonna]]
    
    # STEP 1: Rimuovi outlier (IQR)
    Q1 <- quantile(x_raw, 0.25, na.rm = TRUE)
    Q3 <- quantile(x_raw, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower_bound <- Q1 - 3 * IQR_val
    upper_bound <- Q3 + 3 * IQR_val
    x_clean <- x_raw[x_raw >= lower_bound & x_raw <= upper_bound]
    n_outliers <- length(x_raw) - length(x_clean)
    
    # STEP 2: Centra
    mu_raw <- mean(x_clean, na.rm = TRUE)
    x_centered <- x_clean - mu_raw
    
    # STEP 3: Statistiche base
    sigma_noise <- sd(x_centered, na.rm = TRUE)
    skew <- skewness(x_centered)
    kurt <- kurtosis(x_centered)
    
    # STEP 4: Stampa preliminare
    cat("\n--- COLONNA:", colonna, "---\n")
    cat("Campioni:", length(x_raw), "| Outlier rimossi:", n_outliers, "\n")
    cat("Media:", round(mu_raw, 4), "| SD:", round(sigma_noise, 4), "\n")
    cat("Skewness:", round(skew, 3), "| Kurtosis:", round(kurt, 3), "\n")
    
    # STEP 5: CONFRONTO DISTRIBUZIONI
    model_result <- compare_distributions(x_centered, colonna)
    
    # STEP 6: Salva risultati su file
    stats_file <- file.path(output_dir,
                            paste0(sub(".csv", "", basename(file)),
                                   "_", colonna, "_model_selection.txt"))
    
    sink(stats_file)
    cat("=== MODEL SELECTION REPORT ===\n")
    cat("File:", basename(file), "\n")
    cat("Colonna:", colonna, "\n")
    cat("Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    
    cat("--- STATISTICHE DESCRITTIVE ---\n")
    cat("N totali:", length(x_raw), "\n")
    cat("Outlier rimossi (IQR method):", n_outliers, "\n")
    cat("N analizzati:", length(x_centered), "\n\n")
    
    cat("Media originale:", mu_raw, "\n")
    cat("SD rumore:", sigma_noise, "\n")
    cat("Skewness:", skew, "(0 = simmetrico)\n")
    cat("Kurtosis:", kurt, "(3 = gaussiana)\n\n")
    
    cat("--- CONFRONTO MODELLI ---\n")
    print(model_result$results)
    
    cat("\n--- CONCLUSIONE ---\n")
    cat(model_result$conclusion, "\n")
    
    sink()
    
    # STEP 7: Grafico comparativo (QQ-plot)
    qq_file <- file.path(output_dir,
                         paste0(sub(".csv", "", basename(file)),
                                "_", colonna, "_model_qqplot.png"))
    png(qq_file, width = 1000, height = 700)
    qqnorm(sample(x_centered, min(1000, length(x_centered))),
           main = paste("QQ-Plot:", colonna, "-", basename(file)))
    qqline(sample(x_centered, min(1000, length(x_centered))), col = "red", lwd = 2)
    dev.off()
  }
}

cat("\n\n=== ANALISI COMPLETATA ===\n")
cat("Risultati salvati in:", output_dir, "\n")
cat("Cerca i file *_model_selection.txt per le conclusioni dettagliate\n")
