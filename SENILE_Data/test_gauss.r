# Questo script esegue un test statistico per verificare se un insieme di dati
# segue una distribuzione gaussiana.
# Lo script carica i dati dal workspace SENILE, esegue e salva un boxplot;
# salva la differenza tra il numero di dati totali ed il numero di dati dopo
# aver rimosso gli outliers; esegue il test di normalità di Shapiro-Wilk, ne
# stampa i risultati; esegue e salva un grafico della distribuzione dei dati
# con la curva gaussiana sovrapposta; esegue e salva un grafico qqplot.

# Impostiamo la working directory alla cartella dello script
setwd(getSrcDirectory(function() {})[1])

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

  for (colonna in colonne){
    # Creiamo e salviamo il boxplot
    boxplot_filename <- file.path(output_dir,
                                  paste0(sub(".csv", "", basename(file)),
                                         "_", colonna, "_boxplot.png"))
    png(boxplot_filename)
    boxplot(dati_selezionati[[colonna]], main = paste("Boxplot di", colonna))
    dev.off()

    # Rimuoviamo gli outliers usando il metodo IQR
    q1 <- quantile(dati_selezionati[[colonna]], 0.25, na.rm = TRUE)
    q3 <- quantile(dati_selezionati[[colonna]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    dati_senza_outliers <- dati_selezionati[[colonna]][dati_selezionati[[colonna]] >= (q1 - 3 * iqr) &
                                                      dati_selezionati[[colonna]] <= (q3 + 3 * iqr)]
    num_outliers <- nrow(dati_selezionati) - length(dati_senza_outliers)
    print(paste("Numero di outliers rimossi per", colonna, "nel file", file, ":", num_outliers))

    # Eseguiamo il test di normalità di Shapiro-Wilk per ogni colonna
    # Prendiamo solo gli ultimi min(100, nrow) dati per non avere un test
    # iper sensibile
    risultato_test <- shapiro.test(tail(dati_senza_outliers,
      min(100, length(dati_senza_outliers))))
    print(paste("Risultati per", colonna, "nel file", file))
    print(risultato_test)

    # Eseguiamo e salviamo il plot della distribuzione dei dati con la curva
    # gaussiana sovrapposta
    png_filename <- file.path(output_dir,
                              paste0(sub(".csv", "", basename(file)),
                                     "_", colonna, "_distribution.png"))
    png(png_filename)
    hist(dati_senza_outliers, breaks = 30, probability = TRUE,
         main = paste("Distribuzione di", colonna),
         xlab = colonna)
    mu <- mean(dati_senza_outliers, na.rm = TRUE)
    sigma <- sd(dati_senza_outliers, na.rm = TRUE)
    curve(dnorm(x, mean = mu, sd = sigma), col = "blue", lwd = 2, add = TRUE)
    dev.off()

    # Eseguiamo e salviamo il QQ plot
    qqplot_filename <- file.path(output_dir,
                                 paste0(sub(".csv", "", basename(file)),
                                        "_", colonna, "_qqplot.png"))
    png(qqplot_filename)
    qqnorm(dati_senza_outliers, main = paste("QQ Plot di", colonna))
    qqline(dati_senza_outliers, col = "red")
    dev.off()
  }

  # Salviamo un file riassuntivo dei risultati per ogni file
  summary_file <- file.path(output_dir,
                            paste0(sub(".csv", "", basename(file)),
                                   "_summary.txt"))
  summary_lines <- c()
  for (colonna in colonne){
    # Rimuoviamo gli outliers usando il metodo IQR
    q1 <- quantile(dati_selezionati[[colonna]], 0.25, na.rm = TRUE)
    q3 <- quantile(dati_selezionati[[colonna]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    dati_senza_outliers <- dati_selezionati[[colonna]][dati_selezionati[[colonna]] >= (q1 - 1.5 * iqr) &
                                                      dati_selezionati[[colonna]] <= (q3 + 1.5 * iqr)]
    # Calcoliamo il numero di outliers rimossi
    num_outliers <- nrow(dati_selezionati) - length(dati_senza_outliers)
    # Eseguiamo il test di Shapiro-Wilk sui dati senza outliers
    # Per spiegare perché si fanno due test, uno con 100 e uno con 5000 dati:
    # il test di Shapiro-Wilk è molto sensibile alle dimensioni del campione.
    # Con campioni molto grandi, anche piccole deviazioni dalla normalità
    # possono risultare in un rifiuto dell'ipotesi nulla. Eseguendo il test
    # su un campione più piccolo (100 dati) si ottiene una valutazione più
    # permissiva della normalità dei dati.
    # Per alcuni esempi si fa riferimento al seguente link:
    # https://statorials.com/shapiro-wilk-test-in-r/
    # Si segnala che la fonte non è accademica e quindi non è ritenuta
    # completamente affidabile.
    risultato_test_completo <- shapiro.test(tail(dati_senza_outliers, 
      min(5000, length(dati_senza_outliers))))
    risultato_test <- shapiro.test(tail(dati_senza_outliers,
      min(100, length(dati_senza_outliers))))
    mean_val <- mean(dati_senza_outliers, na.rm = TRUE)
    sd_val <- sd(dati_senza_outliers, na.rm = TRUE)
    summary_lines <- c(summary_lines,
                       paste("Risultati per", colonna, ":"),
                       paste("Dati dentro il file originale:", nrow(dati_selezionati)),
                       paste("Dati dopo rimozione outliers:", length(dati_senza_outliers)),
                       paste("Numero di outliers rimossi:", num_outliers),
                       paste("Media:", mean_val),
                       paste("Deviazione Standard:", sd_val),
                       capture.output(risultato_test),
                       capture.output(risultato_test_completo),
                       "------------------------------------------------", "")
  }
  writeLines(summary_lines, summary_file)
}