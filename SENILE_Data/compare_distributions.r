# ============================================================================
# FUNZIONE PRINCIPALE: Confronto fit simultaneo: Normale, t-Student, Laplace
# ============================================================================

compare_distributions <- function(x_centered, name_col) {
  cat("\n===== CONFRONTO DISTRIBUZIONI:", name_col, "=====\n")

  # Rimuovi NA
  x_centered <- x_centered[!is.na(x_centered)]

  # Data frame per raccogliere risultati dei tre fit
  results <- data.frame(
    Distribution = c("Normale", "t-Student", "Laplace"),
    LogLik = NA_real_,
    AIC = NA_real_,
    AD_pvalue = NA_real_,
    stringsAsFactors = FALSE
  )

  # df stimati della t-Student (servono per regolare la soglia di parsimonia)
  df_t <- NA_real_

  # ==========================================================================
  # 1. FIT NORMALE (Gaussiano)
  # ==========================================================================
  cat("Fitting Normale... ")
  fit_norm <- tryCatch({
    MASS::fitdistr(x_centered, "normal")
  }, error = function(e) {
    return(NULL)
  })

  if (!is.null(fit_norm)) {
    mu_n <- fit_norm$estimate["mean"]
    sigma_n <- fit_norm$estimate["sd"]
    ll_norm <- loglik_normal(x_centered, mu_n, sigma_n)
    results[1, "LogLik"] <- ll_norm
    results[1, "AIC"] <- -2 * ll_norm + 2 * 2  # 2 parametri: mu, sigma

    # Test Anderson-Darling: H0=dati da Normale; p>0.05 supporta H0
    ad_norm <- tryCatch(nortest::ad.test(x_centered), error = function(e)
                          list(p.value = NA))
    results[1, "AD_pvalue"] <- ad_norm$p.value
    cat("OK\n")
  } else {
    cat("SKIPPED\n")
  }

  # ==========================================================================
  # 2. FIT t-STUDENT (distribuzioni con code pesanti)
  # ==========================================================================
  # Strategia: Prova fitdistr() con densità custom; se fallisce, usa optim
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
    kurt_obs <- MASS::kurtosis(x_centered)
    df_start <- if (is.finite(kurt_obs) && kurt_obs > 3) {
      4 + 6 / (kurt_obs - 3)
    } else {
      10  # fallback se kurt <= 3
    }

    suppressWarnings(
      MASS::fitdistr(
        x_centered,
        t_density,
        start = list(m = median(x_centered), s = s_start, df = df_start),
        lower = c(-Inf, sqrt(.Machine$double.eps), sqrt(.Machine$double.eps))
      )
    )
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
    results[2, "LogLik"] <- ll_t
    results[2, "AIC"] <- -2 * ll_t + 2 * 3  # 3 parametri: mu, sigma, df

    # Anderson-Darling: H0=dati da t-Student
    ad_t <- tryCatch(nortest::ad.test(x_centered), 
                     error = function(e) list(p.value = NA))
    results[2, "AD_pvalue"] <- ad_t$p.value
    cat("OK (df=", round(df_t, 2), ")\n", sep = "")
  } else {
    cat("SKIPPED\n")
  }

  # ==========================================================================
  # 3. FIT LAPLACE (doppia esponenziale, picco acuto)
  # ==========================================================================
  cat("Fitting Laplace... ")
  fit_lap <- tryCatch({
    fitdistr_laplace(x_centered)
  }, error = function(e) {
    return(NULL)
  })

  if (!is.null(fit_lap)) {
    mu_l <- fit_lap$mu
    scale_l <- fit_lap$scale
    ll_lap <- loglik_laplace(x_centered, mu_l, scale_l)
    results[3, "LogLik"] <- ll_lap
    results[3, "AIC"] <- -2 * ll_lap + 2 * 2  # 2 parametri: mu, scale

    ad_lap <- tryCatch(ad.test(x_centered), error = function(e) list(p.value = NA))
    results[3, "AD_pvalue"] <- ad_lap$p.value
    cat("OK\n\n")
  } else {
    cat("SKIPPED\n\n")
  }

  # ==========================================================================
  # RANKING E CONCLUSIONI
  # ==========================================================================

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
  if (!is.na(results[1, "AIC"]) && !is.na(results[2, "AIC"])) {
    delta_nt <- results[1, "AIC"] - results[2, "AIC"]
    cat("Normale - t-Student:", round(delta_nt, 2))
    if (abs(delta_nt) < 2) cat(" (AMBIGUO)\n")
    else if (delta_nt > 0) cat(" (t-Student migliore)\n")
    else cat(" (Normale migliore)\n")
  }

  if (!is.na(results[1, "AIC"]) && !is.na(results[3, "AIC"])) {
    delta_nl <- results[1, "AIC"] - results[3, "AIC"]
    cat("Normale - Laplace:", round(delta_nl, 2))
    if (abs(delta_nl) < 2) cat(" (AMBIGUO)\n")
    else if (delta_nl > 0) cat(" (Laplace migliore)\n")
    else cat(" (Normale migliore)\n")
  }

  # Anderson-Darling goodness-of-fit test: p > 0.05 significa "compatibile"
  cat("\n--- ANDERSON-DARLING p-values (p > 0.05 = compatibile) ---\n")
  for (i in 1:3) {
    pval <- results[i, "AD_pvalue"]
    dist <- results[i, "Distribution"]
    if (!is.na(pval)) {
      status <- if (pval > 0.05) "✓ Compatibile" else "✗ Respinto"
      cat(dist, ":", round(pval, 4), status, "\n")
    } else {
      cat(dist, ": N/A\n")
    }
  }

  # Conclusione finale: logica basata su AIC + compatibilita' AD, con
  # preferenza al modello piu' semplice se il vantaggio e' minimo
  cat("\n=== CONCLUSIONE ===\n")

  conclusion <- "INDETERMINATO (errore nel fit)"

  if (!is.na(best_dist)) {
    norm_pval <- results[1, "AD_pvalue"]
    t_pval <- results[2, "AD_pvalue"]
    lap_pval <- results[3, "AD_pvalue"]
    norm_aic <- results[1, "AIC"]

    # Soglia di parsimonia: se la t ha df alti, richiediamo un vantaggio
    # AIC < df
    parsimony_tol <- function(df) {
      if (is.na(df) || df <= 0) return(2)
      df
    }

    # Se il best non è Normale ma la Normale è entro la soglia di parsimonia,
    # scegliamo la Normale (modello più semplice; Laplace ha la stessa
    # complessità di Normale, t ha 1 parametro in più)
    if (!is.na(norm_aic) && !is.na(best_aic) && best_dist != "Normale") {
      if ((norm_aic - best_aic) < parsimony_tol(df_t)) {
        best_dist <- "Normale"
      }
    }

    if (best_dist == "Normale") {
      if (!is.na(norm_pval) && norm_pval > 0.05) {
        conclusion <- "✓ GAUSSIANO (Normale compatibile, AD p > 0.05)"
      } else {
        conclusion <- "~ Normale vince AIC ma AD respinge (p <= 0.05)"
      }
    } else if (best_dist == "t-Student") {
      if (!is.na(t_pval) && t_pval > 0.05) {
        conclusion <- "✓ t-STUDENT (AIC migliore e AD compatibile)"
      } else {
        conclusion <- "~ t-Student vince AIC ma AD respinge (p <= 0.05)"
      }
    } else if (best_dist == "Laplace") {
      if (!is.na(lap_pval) && lap_pval > 0.05) {
        conclusion <- "✓ LAPLACE (AIC migliore e AD compatibile)"
      } else {
        conclusion <- "~ Laplace vince AIC ma AD respinge (p <= 0.05)"
      }
    }
  }

  cat(conclusion, "\n")

  return(list(results = results, conclusion = conclusion))
}