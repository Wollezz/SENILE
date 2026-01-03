# ============================================================================
# FUNZIONE PRINCIPALE: Confronto fit simultaneo: Normale, t-Student, Laplace
# ============================================================================

compare_distributions <- function(x_centered, name_col) {
  cat("\n===== CONFRONTO DISTRIBUZIONI:", name_col, "=====\n")

  # Rimuovi NA
  x_centered <- x_centered[!is.na(x_centered)]

  # Data frame per raccogliere risultati dei tre fit
  results <- data.frame(
    Distribution = c("normal", "t-Student", "Laplace", "Poisson"),
    LogLik = NA_real_,
    AIC = NA_real_,
    stringsAsFactors = FALSE
  )

  # df stimati della t-Student (servono per regolare la soglia di parsimonia)
  df_t <- NA_real_

  # Anderson-Darling test è valido SOLO per Normale; salviamo il p-value separatamente
  ad_normal_pvalue <- NA_real_

  # ==========================================================================
  # 1. FIT DATI
  # ==========================================================================

  # Inizializzazione counter
  counter <- 1

  for (distribution in results$Distribution) {
    if (distribution != "t-Student") {
      cat("Fitting", distribution, " ... ")

      if (distribution == "Laplace") {
        fit_dist <- tryCatch({
          fitdistr_laplace(x_centered)
        }, error = function(e) {
          return(NULL)
        })
      } else {
        fit_dist <- tryCatch({
          MASS::fitdistr(x_centered, distribution)
        }, error = function(e) {
          return(NULL)
        })
      }

      if (!is.null(fit_dist)) {
        if (distribution == "Laplace") {
          mu_dist <- fit_dist$mu
          scale_dist <- fit_dist$scale
          ll_dist <- loglik_laplace(x_centered, mu_dist, scale_dist)
          results[counter, "LogLik"] <- ll_dist
          results[counter, "AIC"] <- -2 * ll_dist + 2 * 2  # 2 parametri: mu, scale
          cat("OK\n")
          next
        } else if (distribution == "normal") {
          mu_dist <- fit_dist$estimate["mean"]
          sigma_dist <- fit_dist$estimate["sd"]
          ll_dist <- loglik_normal(x_centered, mu_dist, sigma_dist)
          results[counter, "LogLik"] <- ll_dist
          results[counter, "AIC"] <- -2 * ll_dist + 2 * 2  # 2 parametri: mu, sigma
        }

        if (distribution == "normal") {
          # Test Anderson-Darling: H0=dati da Normale; p>0.05 supporta H0
          ad_normal_pvalue <- suppressWarnings(tryCatch(nortest::ad.test(x_centered),
                                       error = function(e) list(p.value = NA)))$p.value

          # Shapiro-Wilk test sui dati senza outlier per completezza
          # Il test di Shapiro-Wilk è molto sensibile alle dimensioni del
          # campione. Con campioni molto grandi, anche piccole deviazioni dalla
          # normalità possono risultare in un rifiuto dell'ipotesi nulla.
          # Eseguendo il test su un campione più piccolo (100 dati) si ottiene
          # una valutazione più permissiva della normalità dei dati.
          # Per alcuni esempi si fa riferimento al seguente link:
          # https://statorials.com/shapiro-wilk-test-in-r/
          # Si segnala che la fonte non è accademica e quindi non è ritenuta
          # completamente affidabile.
          # Test di Shapiro-Wilk su un campione di 5000 dati (se possibile)
          shapiro_result_pvalue <- suppressWarnings(tryCatch(shapiro.test(
            tail(x_centered,
                 min(5000, length(x_centered)))
          ), error = function(e) list(p.value = NA)))$p.value

          # Test di Shapiro-Wilk su un campione di 100 dati (se possibile)
          shapiro_result_100 <- suppressWarnings(tryCatch(shapiro.test(
            tail(x_centered,
                 min(100, length(x_centered)))
          ), error = function(e) list(p.value = NA)))$p.value
        }

        cat("OK\n")
      } else {
        cat("SKIPPED\n")
      }
    } else {
      # ========================================================================
      # - FIT t-STUDENT (distribuzioni con code pesanti)
      # ========================================================================
      # Prova fitdistr() con densità custom; se fallisce, usa optim
      # in parametrizzazione log-space per garantire sigma > 0, df > 0
      cat("Fitting t-Student... ")

      # Tentativo 1: fitdistr con densità esplicita dt((x-m)/s)/s
      fit_t <- tryCatch({
        # Densità t con scala positiva: dt((x-m)/s)/s
        # fitdistr usa questa per calcolare MLE via ottimizzazione numerica
        t_density <- function(x, m, s, df) {
          dt((x - m) / s, df = df) / s
        }

        # Starting values estratti dai dati (nessun hardcoding):
        #   s_start = sd(x) ma protetto da underflow
        #   df_start = 4 + 6/(kurt_obs - 3) = formula inversa kurtosis t-Student
        #            (t-Student ha kurtosis = 3 + 6/(df-4) per df>4)
        s_start <- max(sd(x_centered), sqrt(.Machine$double.eps))
        kurt_obs <- suppressWarnings(moments::kurtosis(x_centered))
        df_start <- if (is.finite(kurt_obs) && kurt_obs > 3) {
          4 + 6 / (kurt_obs - 3)
        } else {
          10  # fallback se kurt <= 3
        }

        suppressWarnings(
          moments::fitdistr(x_centered, t_density,
          start = list(m = median(x_centered), s = s_start, df = df_start),
          lower = c(-Inf, sqrt(.Machine$double.eps), sqrt(.Machine$double.eps))))
      }, error = function(e) {
        return(NULL)
      })

      # Tentativo 2: Se fitdistr fallisce, usa ottimizzazione diretta (optim)
      # con parametrizzazione log-space per evitare dominio invalido
      if (is.null(fit_t)) {
        fit_t <- tryCatch({
          # Negative log-likelihood con parametri trasformati:
          #   theta[1] = m (illimitato)
          #   theta[2] = log(s) => s = exp(theta[2]) sempre > 0
          #   theta[3] = log(df) => df = exp(theta[3]) sempre > 0
          # Questo garantisce automaticamente validità dei parametri
          negLL_t <- function(theta) {
            m <- theta[1]
            s <- exp(theta[2])
            df <- exp(theta[3])

            z <- (x_centered - m) / s
            # Calcola log-likelihood con protezione da NaN
            ll <- suppressWarnings(sum(dt(z, df = df, log = TRUE))) - length(x_centered) * log(s)

            # Se non finito (overflow), ritorna penalità massima
            if (!is.finite(ll)) return(sqrt(.Machine$double.xmax))
            return(-ll)  # optim minimizza, quindi ritorna -LL
          }

          # Starting values in spazio trasformato
          theta0 <- c(
            median(x_centered),
            log(max(sd(x_centered), sqrt(.Machine$double.eps))),
            log(max(1, df_start))
          )

          # Ottimizza con BFGS (metodo del gradiente coniugato)
          opt <- suppressWarnings(optim(theta0, negLL_t, method = "BFGS"))

          # Se convergenza OK, estrai parametri
          if (opt$convergence == 0) {
            m_est <- opt$par[1]
            s_est <- exp(opt$par[2])
            df_est <- exp(opt$par[3])
            list(estimate = c(m = m_est, s = s_est, df = df_est))
          } else {
            NULL
          }
        }, error = function(e) {
          NULL
        })
      }

      # Calcola statistiche se il fit è riuscito
      if (!is.null(fit_t)) {
        mu_t <- fit_t$estimate["m"]
        sigma_t <- fit_t$estimate["s"]
        df_t <- fit_t$estimate["df"]

        ll_t <- loglik_t(x_centered, mu_t, sigma_t, df_t)
        results[counter, "LogLik"] <- ll_t
        results[counter, "AIC"] <- -2 * ll_t + 2 * 3  # 3 parametri: mu, sigma, df

        # Per t-Student usiamo solo AIC e log-likelihood per il confronto.
        cat("OK (df=", round(df_t, 2), ")\n", sep = "")
      } else {
        cat("SKIPPED\n")
      }
    }
    counter <- counter + 1
  }

  # ==========================================================================
  # RANKING E CONCLUSIONI
  # ==========================================================================

  cat("\n")

  # Stampa tabella risultati
  print(results)

  # Trova modello vincente (AIC minimo)
  valid_idx <- !is.na(results$AIC)
  best_idx <- which.min(results$AIC[valid_idx])
  if (length(best_idx) > 0) {
    best_dist <- results[valid_idx, ][best_idx, "Distribution"]
    best_aic <- min(results$AIC, na.rm = TRUE)
    cat("\nMigliore modello (AIC):", best_dist, "AIC =", round(best_aic, 2),
        "\n")
  } else {
    best_dist <- NA
    cat("\nNessun modello valido fittato\n")
  }

  # Differenze AIC tra modelli (interpretazione: ΔAIC < 2 ambiguo,
  # > 10 decisivo)
  cat("\n--- DIFFERENZE AIC ---\n")
  while (counter > 2) {
    if (!is.na(results[counter, "AIC"]) && !is.na(results[counter - 1, "AIC"])) {
      delta <- results[counter - 1, "AIC"] - results[counter, "AIC"]
      cat(results$Distribution[counter - 1], "-", results$Distribution[counter], ":", round(delta, 2))
      if (abs(delta) < 2) cat(" (AMBIGUO)\n")
      else if (delta > 0) cat(" (", results$Distribution[counter - 1], " migliore)\n")
      else cat(" (", results$Distribution[counter], " migliore)\n")
    }
    counter <- counter - 1
  }

  # Anderson-Darling test: p > 0.05 significa "compatibile"
  cat("\n--- ANDERSON-DARLING p-value (p > 0.05 = compatibile) ---\n")
  if (!is.na(ad_normal_pvalue)) {
    status <- if (ad_normal_pvalue > 0.05) "✓ Compatibile" else "✗ Respinto"
    cat("Normale:", round(ad_normal_pvalue, 4), status, "\n")
  } else {
    cat("Normale: N/A\n")
  }

  # Shapiro-Wilk test: p > 0.05 significa "compatibile"
  cat("\n--- SHAPIRO-WILK p-value (p > 0.05 = compatibile) ---\n")
  if (!is.na(shapiro_result_pvalue)) {
    status_5000 <- if (shapiro_result_pvalue > 0.05) "✓ Compatibile" else "✗ Respinto"
    cat("Normale (5000 dati):", round(shapiro_result_pvalue, 4), status_5000, "\n")
  } else {
    cat("Normale (5000 dati): N/A\n")
  }
  if (!is.na(shapiro_result_100)) {
    status_100 <- if (shapiro_result_100 > 0.05) "✓ Compatibile" else "✗ Respinto"
    cat("Normale (100 dati):", round(shapiro_result_100, 4), status_100, "\n")
  } else {
    cat("Normale (100 dati): N/A\n")
  }

  # Conclusione finale: basata su AIC + (per Normale) compatibilità AD
  cat("\n=== CONCLUSIONE ===\n")

  conclusion <- "INDETERMINATO (errore nel fit)"

  if (!is.na(best_dist)) {
    norm_aic <- results[results$Distribution == "normal", "AIC"]

    # Soglia di decisione: se la t ha df alti, richiediamo un vantaggio
    # AIC < df
    parsimony_tol <- function(df) {
      if (is.na(df) || df <= 0) return(2)
      df
    }

    # Se il best non è Normale ma la Normale è entro la soglia di decisione,
    # scegliamo la Normale (modello più semplice; Laplace ha la stessa
    # complessità di Normale, t ha 1 parametro in più)
    if (!is.na(norm_aic) && !is.na(best_aic) && best_dist != "normal") {
      if ((norm_aic - best_aic) < parsimony_tol(df_t)) {
        best_dist <- "normal"
      }
    }

    if (best_dist == "normal") {
      if (!is.na(ad_normal_pvalue) && ad_normal_pvalue > 0.05) {
        conclusion <- "✓ GAUSSIANO (Normale miglior AIC e AD compatibile)"
      } else {
        conclusion <- "~ Normale miglior AIC ma AD non compatibile (consideriamo comunque)"
      }
    } else if (best_dist == "t-Student") {
      conclusion <- "✓ t-STUDENT (AIC migliore; AD non applicabile)"
    } else if (best_dist == "Laplace") {
      conclusion <- "✓ LAPLACE (AIC migliore; AD non applicabile)"
    }
  }

  cat(conclusion, "\n")

  return(list(results = results, conclusion = conclusion))
}