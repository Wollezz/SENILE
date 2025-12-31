# Questo script esegue un test statistico per verificare se un insieme di dati
# segue una distribuzione gaussiana quantizzata.
# Lo script:
# 1. Rimuove outlier (spike) usando IQR method
# 2. Centra i dati (sottrae la media per isolare il rumore)
# 3. Calcola statistiche: media, sd, skewness, kurtosis
# 4. Esegue Shapiro-Wilk test su dati centrati e puliti
# 5. Salva grafici e statistiche complete

# Installa/carica pacchetti necessari
if (!require("moments")) install.packages("moments")
library(moments)

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

  # Estraiamo solo le colonne di interesse
  dati_selezionati <- dati[, colonne]

  # Prendiamo gli ultimi min(5000, nrow) dati per gestire file più piccoli
  dati_selezionati <- tail(dati_selezionati, min(5000, nrow(dati_selezionati)))

  for (colonna in colonne){
    x_raw <- dati_selezionati[[colonna]]
    
    # STEP 1: Rimuovi outlier usando metodo IQR (spike casuali)
    Q1 <- quantile(x_raw, 0.25, na.rm = TRUE)
    Q3 <- quantile(x_raw, 0.75, na.rm = TRUE)
    IQR_val <- Q3 - Q1
    lower_bound <- Q1 - 3 * IQR_val  # 3 IQR = outlier estremi
    upper_bound <- Q3 + 3 * IQR_val
    x_clean <- x_raw[x_raw >= lower_bound & x_raw <= upper_bound]
    n_outliers <- length(x_raw) - length(x_clean)
    
    # STEP 2: Centra i dati (rimuovi bias per isolare rumore)
    mu_raw <- mean(x_clean, na.rm = TRUE)
    x_centered <- x_clean - mu_raw
    
    # STEP 3: Calcola statistiche complete
    sigma_noise <- sd(x_centered, na.rm = TRUE)
    skew <- skewness(x_centered)
    kurt <- kurtosis(x_centered)
    
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
    png(png_filename, width=800, height=600)
    hist(x_centered, breaks = 50, probability = TRUE,
         main = paste("Distribuzione rumore -", colonna),
         xlab = "Valore centrato",
         col = "lightblue")
    curve(dnorm(x, mean = 0, sd = sigma_noise), col = "red", lwd = 2, add = TRUE)
    legend("topright", legend = c("Dati", "Gaussiana teorica"),
           col = c("lightblue", "red"), lwd = c(10, 2))
    text(x = par("usr")[1] + 0.05*diff(par("usr")[1:2]),
         y = par("usr")[4] - 0.05*diff(par("usr")[3:4]),
         labels = paste0("SD=", round(sigma_noise, 4), "\nKurt=", round(kurt, 2)),
         pos = 4, cex = 0.9)
    dev.off()
    
    # STEP 7: QQ plot su dati centrati
    qqplot_filename <- file.path(output_dir,
                                 paste0(sub(".csv", "", basename(file)),
                                        "_", colonna, "_qqplot.png"))
    png(qqplot_filename, width=800, height=600)
    qqnorm(x_test, main = paste("QQ Plot -", colonna))
    qqline(x_test, col = "red", lwd = 2)
    dev.off()
  }
}

cat("\n\n=== ANALISI COMPLETATA ===\n")
cat("Risultati salvati in:", output_dir, "\n")