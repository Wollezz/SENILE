# Questo script esegue un test statistico per verificare se un insieme di dati
# segue una distribuzione gaussiana quantizzata.
# Lo script:
# 1. Rimuove outlier (spike) usando IQR method
# 2. Centra i dati (sottrae la media per isolare il rumore)
# 3. Calcola statistiche: media, sd, skewness, kurtosis
# 4. Esegue Shapiro-Wilk test su dati centrati e puliti
# 5. Salva grafici e statistiche complete

# Installa/carica pacchetti necessari
if (!require("moments"))
  install.packages("moments")
if (!require("pracma"))
  install.packages("pracma")
library(moments)
library(pracma)

# Impostiamo la working directory alla cartella dello script (metodo universale)
script_dir <- getSrcDirectory(function() {})[1]
if (script_dir != "") {
  setwd(script_dir)
}

# Cartella di lavoro
cartella <- "Euroc12025_separati/"

# Elenchiamo tutti i file CSV nella cartella specificata
files <- list.files(path = cartella, pattern = "*.csv", full.names = TRUE)

# Cartella per i risultati
output_dir <- file.path(cartella, "Risultati")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Nomi delle colonne da analizzare
colonne <- c("accelerationX", "accelerationY", "accelerationZ",
             "angularSpeedX", "angularSpeedY", "angularSpeedZ")

for (file in files){
  # Carichiamo i dati dal workspace SENILE
  dati <- read.csv(file, header = TRUE)

  # Prendiamo gli ultimi min(5000, nrow) dati per gestire file più piccoli
  n_keep <- min(5000, nrow(dati))
  idx_tail <- (nrow(dati) - n_keep + 1):nrow(dati)
  dati_tail <- dati[idx_tail, ]
  # Estraiamo solo le colonne di interesse
  dati_selezionati <- dati_tail[, colonne]

  for (colonna in colonne){
    x_raw <- dati_selezionati[[colonna]]

    # Se esiste il timestamp associato, usiamolo per Allan (altrimenti usiamo un
    # indice uniforme)
    ts_col <- if (grepl("acceleration", colonna, fixed = TRUE))
      "accelerationTimestamp"
    else if (grepl("angularSpeed", colonna, fixed = TRUE))
      "angularSpeedTimestamp"
    else NA
    t_raw <- if (!is.na(ts_col) && ts_col %in% names(dati_tail))
      dati_tail[[ts_col]]
    else seq_along(x_raw)

    # STEP 1: Rimuovi outlier usando metodo IQR (spike casuali)
    q1 <- quantile(x_raw, 0.25, na.rm = TRUE)
    q3 <- quantile(x_raw, 0.75, na.rm = TRUE)
    iqr_val <- q3 - q1
    lower_bound <- q1 - 3 * iqr_val  # 3 IQR = outlier estremi
    upper_bound <- q3 + 3 * iqr_val
    keep_mask <- x_raw >= lower_bound & x_raw <= upper_bound
    x_clean <- x_raw[keep_mask]
    t_clean <- t_raw[keep_mask]
    n_outliers <- length(x_raw) - length(x_clean)

    # STEP 2: Centra i dati (rimuovi bias per isolare rumore)
    mu_raw <- mean(x_clean, na.rm = TRUE)
    x_centered <- x_clean - mu_raw

    # STEP 3: Calcola statistiche complete
    sigma_noise <- sd(x_centered, na.rm = TRUE)
    skew <- skewness(x_centered)
    kurt <- kurtosis(x_centered)

    # STEP 3b: Allan variance/deviation (robusto) su TUTTA la serie grezza
    x_allan <- as.numeric(stats::na.omit(dati[[colonna]]))
    t_allan <- if (!is.na(ts_col) && ts_col %in% names(dati)) dati[[ts_col]]
    else seq_along(dati[[colonna]])
    t_allan <- t_allan[is.finite(t_allan)]
    # riallinea lunghezza
    n_sync <- min(length(x_allan), length(t_allan))
    x_allan <- x_allan[seq_len(n_sync)]
    t_allan <- t_allan[seq_len(n_sync)]
    dt_vec <- diff(t_allan)
    dt_vec <- dt_vec[is.finite(dt_vec) & dt_vec > 0]
    # Se i timestamp sono troppo irregolari, usa l'indice come tempo uniforme
    use_index_time <- FALSE
    if (length(dt_vec) > 1) {
      cv_dt <- sd(dt_vec) / mean(dt_vec)
      if (!is.finite(cv_dt) || cv_dt > 0.1) {
        use_index_time <- TRUE
      }
    } else {
      use_index_time <- TRUE
    }
    if (use_index_time) {
      t_allan <- seq_along(x_allan)
      dt_vec <- rep(1, length(x_allan) - 1)
    }

    dt_med <- median(dt_vec, na.rm = TRUE)
    # Heuristica unità timestamp: se >1e6 => ns, >1e3 => us, >1 => ms,
    # altrimenti s
    if (is.na(dt_med) || dt_med <= 0) {
      dt_step <- 1
    } else if (dt_med > 1e6) {
      dt_step <- dt_med / 1e9
    } else if (dt_med > 1e3) {
      dt_step <- dt_med / 1e6
    } else if (dt_med > 1) {
      dt_step <- dt_med / 1e3
    } else {
      dt_step <- dt_med
    }

    fs <- if (is.finite(dt_step) && dt_step > 0) 1 / dt_step else 1

    # Tau: griglia robusta predefinita (potenze di 2) limitata a n/2 campioni
    if (length(x_allan) < 2) {
      av <- NULL
      av_tau <- numeric(0)
      adev <- numeric(0)
    } else {
      m_max <- max(1, floor(length(x_allan) / 2))
      m_candidates <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
      m_grid <- m_candidates[m_candidates <= m_max]
      if (length(m_grid) == 0) m_grid <- 1
      if (length(m_grid) == 1 && m_max > m_grid) m_grid <- c(m_grid, m_max)
      # pracma::allanvar non accetta fs; tau = m_grid / fs lo calcoliamo a parte
      av <- tryCatch(allanvar(x_allan, m = m_grid), error = function(e) NULL)
      if (is.null(av)) {
        # Fallback minimale: ADEV a tau=1 step
        av_tau <- if (length(x_allan) >= 2) 1 / fs else numeric(0)
        adev <- if (length(x_allan) >= 2) sqrt(0.5 * mean(diff(x_allan)^2)) else numeric(0)
      } else {
        av_tau <- m_grid / fs
        adev <- sqrt(av)
      }
    }

    # STEP 4: Test di normalità su dati centrati e puliti
    # (Campione limitato a 1000 per evitare iper-sensibilità)
    n_test <- min(1000, length(x_centered))
    x_test <- sample(x_centered, n_test)
    risultato_shapiro <- shapiro.test(x_test)

    # Stampa risultati
    cat("\n=== Risultati per", colonna, "in", basename(file), "===\n")
    cat("Campioni totali:", length(x_raw), "\n")
    cat("Outlier rimossi:", n_outliers, "\n")
    cat("Media originale:", mu_raw, "\n")
    cat("SD rumore:", sigma_noise, "\n")
    cat("Varianza rumore:", sigma_noise^2, "\n")
    cat("Skewness:", skew, "(gaussiana = 0)\n")
    cat("Kurtosis:", kurt, "(gaussiana = 3)\n")
    cat("Shapiro-Wilk p-value:", risultato_shapiro$p.value, "\n")
    cat("Sample rate (Hz) stimato:", fs, "\n")

    # STEP 5: Salva statistiche in file di testo
    stats_file <- file.path(output_dir, paste0(sub(".csv", "", basename(file)),
                                               "_", colonna, "_statistics.txt"))
    sink(stats_file)
    cat("=== STATISTICHE SENSORE ===\n")
    cat("File:", basename(file), "\n")
    cat("Colonna:", colonna, "\n")
    cat("Campioni totali:", length(x_raw), "\n")
    cat("Outlier rimossi:", n_outliers, "\n\n")
    cat("--- Dati originali ---\n")
    cat("Media:", mu_raw, "\n")
    cat("Range: [", min(x_clean), ",", max(x_clean), "]\n\n")
    cat("--- Rumore (dati centrati) ---\n")
    cat("Standard Deviation:", sigma_noise, "\n")
    cat("Varianza:", sigma_noise^2, "\n")
    cat("Skewness:", skew, "(0=simmetrico)\n")
    cat("Kurtosis:", kurt, "(3=gaussiana, >3=picco alto)\n\n")
    cat("--- Allan variance ---\n")
    cat("Sample rate stimato (Hz):", fs, "\n")
    cat("Tau usati (s):", if (length(av_tau) > 0) paste(round(av_tau, 6), collapse = ", ") else "n/a", "\n\n")
    cat("--- Test di normalità (Shapiro-Wilk) ---\n")
    cat("Campione testato:", n_test, "punti\n")
    print(risultato_shapiro)
    cat("\nInterpretazione:\n")
    if (risultato_shapiro$p.value > 0.05) {
      cat("p > 0.05: Dati compatibili con gaussiana\n")
    } else {
      cat("p < 0.05: Dati NON gaussiani\n")
      if (abs(kurt - 3) > 1) {
        cat("Kurtosis =", kurt, "-> distribuzione leptocurtica (picco alto)\n")
      }
    }
    sink()

    # STEP 6: Plot distribuzione dati centrati con gaussiana
    png_filename <- file.path(output_dir,
                              paste0(sub(".csv", "", basename(file)),
                                     "_", colonna, "_distribution.png"))
    png(png_filename, width = 800, height = 600)
    hist(x_centered, breaks = 50, probability = TRUE,
         main = paste("Distribuzione rumore -", colonna),
         xlab = "Valore centrato",
         col = "lightblue")
    curve(dnorm(x, mean = 0, sd = sigma_noise), col = "red", lwd = 2,
          add = TRUE)
    legend("topright", legend = c("Dati", "Gaussiana teorica"),
           col = c("lightblue", "red"), lwd = c(10, 2))
    text(x = par("usr")[1] + 0.05 * diff(par("usr")[1:2]),
         y = par("usr")[4] - 0.05 * diff(par("usr")[3:4]),
         labels = paste0("SD=", round(sigma_noise, 4), "\nKurt=",
                         round(kurt, 2)), pos = 4, cex = 0.9)
    dev.off()

    # STEP 7: QQ plot su dati centrati
    qqplot_filename <- file.path(output_dir,
                                 paste0(sub(".csv", "", basename(file)),
                                        "_", colonna, "_qqplot.png"))
    png(qqplot_filename, width = 800, height = 600)
    qqnorm(x_test, main = paste("QQ Plot -", colonna))
    qqline(x_test, col = "red", lwd = 2)
    dev.off()

    # STEP 8: Allan deviation plot (log-log) e CSV
    adev_plot <- file.path(output_dir,
                           paste0(sub(".csv", "", basename(file)),
                                  "_", colonna, "_allan.png"))
    png(adev_plot, width = 900, height = 700)
        if (length(av_tau) > 0) {
       plot(av_tau, adev, log = "xy", type = "b",
         xlab = "Tau [s]", ylab = "Allan deviation",
         main = paste("Allan deviation -", colonna))
       grid()
       if (length(av_tau) == 1) {
         mtext("Solo 1 tau disponibile: serie troppo corta o timestamps non uniformi",
            side = 3, line = 0.5, col = "red", cex = 0.9)
       }
        } else {
       plot(1, type = "n", xlab = "Tau", ylab = "ADEV",
         main = "Allan non disponibile")
       text(1, 1, "Allan non calcolabile: pochi punti o timestamp assenti",
         cex = 1.2)
        }
    dev.off()

    adev_csv <- file.path(output_dir,
                          paste0(sub(".csv", "", basename(file)),
                                 "_", colonna, "_allan.csv"))
    if (length(av_tau) > 0) {
      write.csv(data.frame(tau = av_tau, adev = adev), adev_csv,
                row.names = FALSE)
    }
  }
}

cat("\n\n=== ANALISI COMPLETATA ===\n")
cat("Risultati salvati in:", output_dir, "\n")