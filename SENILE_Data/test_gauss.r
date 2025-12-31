# Questo script esegue un test statistico per verificare se un insieme di dati
# segue una distribuzione gaussiana.
# Lo script carica i dati dal workspace SENILE, esegue il test di normalità di
# Shapiro-Wilk, ne stampa i risultati; esegue, mostra e salva un grafico della
# distribuzione dei dati con la curva gaussiana sovrapposta; esegue, mostra e
# salva un grafico qqplot.

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
    # Eseguiamo il test di normalità di Shapiro-Wilk per ogni colonna
    risultato_test <- shapiro.test(dati_selezionati[[colonna]])
    print(paste("Risultati per", colonna, "nel file", file))
    print(risultato_test)

    # Salviamo il risultato del test in un file di testo
    output_file <- file.path(output_dir, paste0(sub(".csv", "", basename(file)),
                                                "_", colonna,
                                                "_shapiro_test.txt"))
    writeLines(capture.output(risultato_test), output_file)

    # Eseguiamo il plot della distribuzione dei dati con la curva gaussiana
    # sovrapposta
    png_filename <- file.path(output_dir,
                              paste0(sub(".csv", "", basename(file)),
                                     "_", colonna, "_distribution.png"))
    png(png_filename)
    hist(dati_selezionati[[colonna]], breaks = 30, probability = TRUE,
         main = paste("Distribuzione di", colonna),
         xlab = colonna)
    mu <- mean(dati_selezionati[[colonna]], na.rm = TRUE)
    sigma <- sd(dati_selezionati[[colonna]], na.rm = TRUE)
    curve(dnorm(x, mean = mu, sd = sigma), col = "blue", lwd = 2, add = TRUE)
    dev.off()

    # Eseguiamo il QQ plot
    qqplot_filename <- file.path(output_dir,
                                 paste0(sub(".csv", "", basename(file)),
                                        "_", colonna, "_qqplot.png"))
    png(qqplot_filename)
    qqnorm(dati_selezionati[[colonna]], main = paste("QQ Plot di", colonna))
    qqline(dati_selezionati[[colonna]], col = "red")
    dev.off()
  }
}