# ============================================================================
# FUNZIONI AUSILIARIE per fit distribuzioni
# ============================================================================

# Fit Laplace manuale (mediana + deviazione assoluta media)
# La distribuzione di Laplace (doppia esponenziale) ha densità:
#   f(x) = 1/(2*scale) * exp(-|x - mu| / scale)
# Stimatori non-parametrici: mu = mediana, scale = MAD - Mean Absolute Deviation
# scale = E[|X - mu|]
fitdistr_laplace <- function(x) {
  mu <- median(x)
  scale <- mean(abs(x - mu))
  return(list(mu = mu, scale = scale))
}

# Log-likelihood per distribuzione Normale
# LL(x; mu, sigma) = sum(log(dnorm(x, mu, sigma)))
loglik_normal <- function(x, mu, sigma) {
  sum(dnorm(x, mu, sigma, log = TRUE))
}

# Log-likelihood per distribuzione t-Student
# LL(x; mu, sigma, df) = sum(log(dt((x-mu)/sigma, df))) - n*log(sigma)
# (fattore log(sigma) perché usiamo parametrizzazione con scala)
loglik_t <- function(x, mu, sigma, df) {
  sum(dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma))
}

# Log-likelihood per distribuzione Laplace (doppia esponenziale)
# LL(x; mu, scale) = sum(log(1/(2*scale)) - |x - mu| / scale)
loglik_laplace <- function(x, mu, scale) {
  sum(log(1 / (2 * scale)) - abs(x - mu) / scale)
}