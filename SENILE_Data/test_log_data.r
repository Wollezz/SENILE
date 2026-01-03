# ============================================================================
# SCRIPT: Analisi statistica distribuzionale di dati sensori MEMS (IMU)
# ============================================================================
# Obiettivo: Classificare dati accelerometro/giroscopio in tre categorie:
#   1. GAUSSIANO: distribuzione Normale (processo casuale puro)
#   2. t-STUDENT: code pesanti (outlier/anomalie)
#   3. LAPLACE: picco elevato (quantizzazione dominante)
#
# Metodo: Model Selection via AIC (Akaike Information Criterion)
#   - Fitta tre distribuzioni ai dati
#   - Calcola log-likelihood di ognuna
#   - AIC = -2*LL + 2*k (k=parametri). AIC basso = modello migliore.
#   - Usa Anderson-Darling test per validare fit (p > 0.05 = compatibile)
#
# Input: CSV con colonne [accelerationX/Y/Z, angularSpeedX/Y/Z]
# Output: Report *_model_selection.txt per ogni log
# ============================================================================

# Installazione/caricamento pacchetti necessari
packages <- c("moments", "MASS", "nortest")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Impostazione della working directory
setwd(getSrcDirectory(function() {})[1])

# Cartella di lavoro
cartella <- "Euroc12025_separati/"

# Caricamento funzioni ausiliarie da file esterni
source("fit_distribution_functions.r")
source("compare_distributions.r")

# Elenco di tutti i file CSV nella cartella specificata
files <- list.files(path = cartella, pattern = "*.csv", full.names = TRUE)

# Cartella per i risultati
output_dir <- file.path(cartella, "Risultati")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Nomi delle colonne da analizzare
colonne <- c("accelerationX", "accelerationY", "accelerationZ",
             "angularSpeedX", "angularSpeedY", "angularSpeedZ")

# Inizializzazione contatore file
file_count <- 0

# Ciclo per processare ogni file CSV
for (file in files) {
  file_count <- file_count + 1 # aggiorna contatore

  # Output a console
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("FILE", file_count, "DI", length(files), ":", basename(file), "\n")
  cat(rep("=", 70), "\n", sep = "")

  # Lettura file
  dati <- tryCatch({
    read.csv(file, header = TRUE)
  }, error = function(e) {
    cat("ERRORE LETTURA FILE:", e$message, "\n")
    return(NULL)
  })

  # Salto al file successivo in caso di errore
  if (is.null(dati)) next

  # Selezione delle colonne di interesse e prendi ultimi 5000 campioni
  dati_selezionati <- dati[, colonne]

  # SALVATAGGIO RISULTATI SU FILE
  stats_file <- file.path(output_dir,
                          paste0(sub(".csv", "", basename(file)),
                                 "_model_selection.txt"))

  # Reindirizzazione output a file di testo
  sink(stats_file)
  cat("=== MODEL SELECTION REPORT ===\n")
  cat("File:", basename(file), "\n")
  cat("Timestamp:", format(Sys.time()), "\n\n")

  # Analisi di ogni colonna
  for (colonna in colonne) {
    x_raw <- dati_selezionati[[colonna]]

    # PREPROCESSING: Rimozione outlier via metodo IQR
    # Outlier = valori oltre ±3×IQR dai quartili
    q1 <- quantile(x_raw, 0.25, na.rm = TRUE)
    q3 <- quantile(x_raw, 0.75, na.rm = TRUE)
    iqr_val <- q3 - q1
    x_clean <<- x_raw[x_raw >= q1 - 3 * iqr_val & x_raw <= q3 + 3 * iqr_val]
    n_outliers <- length(x_raw) - length(x_clean)

    # CENTERING: Sottrai media per isolare il rumore dalla componente DC
    mu_raw <- mean(x_clean, na.rm = TRUE)
    x_centered <- x_clean - mu_raw

    # STATISTICHE DESCRITTIVE
    sigma_noise <- sd(x_centered, na.rm = TRUE)
    skew <- skewness(x_centered)  # 0 = simmetrica, >0 = coda destra
    kurt <- kurtosis(x_centered)  # 3 = normocurtica / mesocurtica,
    # >3 = leptocurtica ("appuntita")

    cat("\n", rep("=", 70), "\n", sep = "")
    cat("COLONNA:", colonna, "\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("N totali:", length(x_raw), "\n")
    cat("N outlier rimossi:", n_outliers, "\n")
    cat("N analizzati:", length(x_centered), "\n\n")

    cat("Media:", mu_raw, "\nSD:", sigma_noise, "\n")
    cat("Skewness:", skew, "(0=simmetrico)\n")
    cat("Kurtosis:", kurt, "(3=gaussiana)\n")

    # ANALISI PRINCIPALE: Fit delle tre distribuzioni
    model_result <- compare_distributions(x_centered, colonna)

    # Shapiro-Wilk test sui dati senza outlier per completezza
    # Il test di Shapiro-Wilk è molto sensibile alle dimensioni del campione.
    # Con campioni molto grandi, anche piccole deviazioni dalla normalità
    # possono risultare in un rifiuto dell'ipotesi nulla. Eseguendo il test
    # su un campione più piccolo (100 dati) si ottiene una valutazione più
    # permissiva della normalità dei dati.
    # Per alcuni esempi si fa riferimento al seguente link:
    # https://statorials.com/shapiro-wilk-test-in-r/
    # Si segnala che la fonte non è accademica e quindi non è ritenuta
    # completamente affidabile.
    shapiro_result <- shapiro.test(tail(x_clean,
                                        min(100, length(x_clean))))
    cat("\n--- TEST DI SHAPIRO-WILK SU 100 DATI SENZA OUTLIER ---\n")
    print(shapiro_result)
    shapiro_result <- shapiro.test(tail(x_clean,
                                        min(5000, length(x_clean))))
    cat("\n--- TEST DI SHAPIRO-WILK SU 5000 DATI SENZA OUTLIER ---\n")
    print(shapiro_result)

    # Creazione e salvataggio dei plot (istogramma + fit, QQ-plot)
    # Istogramma con curva gaussiana sovrapposta
    png_filename <- file.path(output_dir,
                              paste0(sub(".csv", "", basename(file)),
                                     "_", colonna, "_distribution.png"))
    png(png_filename)
    hist(x_clean, breaks = 30, probability = TRUE,
         main = paste("Distribuzione di", colonna),
         xlab = colonna)
    mu <- mean(x_clean, na.rm = TRUE)
    sigma <- sd(x_clean, na.rm = TRUE)
    curve(dnorm(x, mean = mu, sd = sigma), col = "blue", lwd = 2, add = TRUE)
    dev.off()

    # QQ plot
    qqplot_filename <- file.path(output_dir,
                                 paste0(sub(".csv", "", basename(file)),
                                        "_", colonna, "_qqplot.png"))
    png(qqplot_filename)
    qqnorm(x_clean, main = paste("QQ Plot di", colonna))
    qqline(x_clean, col = "red")
    dev.off()
  }

  # Chiusura del sink per tornare all'output console
  sink()

  # Output a console
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("FILE", file_count, "COMPLETATO:", basename(file), "\n")
  cat("Salvato:", stats_file, "\n")
  cat(rep("=", 70), "\n", sep = "")
}

# Stampa del messaggio finale
cat("\n\n", rep("=", 70), "\n", sep = "")
cat("ANALISI COMPLETATA\n")
cat("Risultati in:", output_dir, "\n")
cat("Totale file elaborati:", file_count, "\n")