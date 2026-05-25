library(mvtnorm)
library(MASS)
library(coda)
library(foreach)
library(doParallel)
library(parallel)
library(readr)
library(dplyr)

# ==============================================================================
# SECTION 1: Load and Prepare Data
# ==============================================================================
raw_data <- read_csv("test_cenluz_empi_2.csv")

raw_data <- raw_data %>%
  mutate(
    Year_num = as.integer(sub("-.*", "", week)),
    Week_num = as.integer(sub(".*-", "", week))
  ) %>%
  arrange(Year_num, Week_num) %>%
  filter(Year_num >= 2014 & Year_num <= 2025)

Y  <- as.integer(raw_data$lepto_count)
X1 <- raw_data$total_rain
X2 <- raw_data$mean_temp
n  <- length(Y)

cat("=== Descriptive Statistics ===\n")
cat("Total observations:", n, "\n")
cat("Number of zeros:   ", sum(Y == 0), "\n")
cat("Proportion zeros:  ", round(sum(Y == 0) / n, 4), "\n")
cat("Mean:              ", round(mean(Y), 4), "\n")
cat("Variance:          ", round(var(Y), 4), "\n")
cat("Var/Mean ratio:    ", round(var(Y) / mean(Y), 4), "\n")
cat("Min:", min(Y), " Max:", max(Y), "\n\n")

extreme_idx <- which(Y > 200)
if (length(extreme_idx) > 0) {
  cat("Observations above 200:\n")
  cat("Indices:", extreme_idx, "\n")
  cat("Values: ", Y[extreme_idx], "\n")
  cat("Weeks:  ", raw_data$week[extreme_idx], "\n\n")
}

X1_scaled <- as.numeric(scale(X1))
X2_scaled <- as.numeric(scale(X2))

X1_mat <- matrix(X1_scaled, 1, ncol = n)
X2_mat <- matrix(X2_scaled, 1, ncol = n)

# ==============================================================================
# SECTION 2: BNB Distribution Utilities
# ==============================================================================
rbnb <- function(n, r, gamma_t, phi) {
  p <- rbeta(n, phi, gamma_t)
  y <- rnbinom(n, size = r, prob = p)
  return(y)
}

rbnb_truncated <- function(n, r, gamma_t, phi) {
  result <- integer(n)
  for (i in 1:n) {
    repeat {
      val <- rbnb(1, r, gamma_t, phi)
      if (val > 0) { result[i] <- val; break }
    }
  }
  return(result)
}

# ==============================================================================
# SECTION 3: Parameter Constraints
# ==============================================================================
check_alpha_constraints <- function(alpha_prop) {
  return(abs(alpha_prop[2]) < 1 &&
           abs(alpha_prop[3]) < 1 &&
           abs(alpha_prop[2] + alpha_prop[3]) < 1)
}

check_beta_constraints <- function(beta_prop) {
  return(beta_prop[1] > 0 &&
           beta_prop[2] > 0 &&
           beta_prop[3] >= 0 &&
           (beta_prop[2] + beta_prop[3]) < 1)
}

# ==============================================================================
# SECTION 4: Log-Likelihood Function
# ==============================================================================
log_likelihood <- function(Y, alpha, beta, omega, r, phi, X_list, b, b0, k) {
  n        <- length(Y)
  log_like <- 0
  pi_t     <- rep(0.08, n)
  lambda_t <- rep(mean(Y[Y > 0]), n)
  
  for (t in (b0 + 1):n) {
    logit_pi_prev <- log(pi_t[t - 1] / (1 - pi_t[t - 1]))
    exponent_pi   <- alpha[1] + alpha[2] * Y[t - 1] + alpha[3] * logit_pi_prev
    pi_t[t]       <- exp(exponent_pi) / (1 + exp(exponent_pi))
    
    lambda_t[t] <- beta[1] + beta[2] * Y[t - 1] + beta[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega[i] * X_list[[i]][1, t - b[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    gamma_t     <- (phi - 1) / r * lambda_t[t]
    
    indicator_y_zero <- ifelse(Y[t] == 0, 1, 0)
    
    log_fbnb_y <- lgamma(Y[t] + r) - lgamma(Y[t] + 1) - lgamma(r) +
      lbeta(phi + r, gamma_t + Y[t]) - lbeta(phi, gamma_t)
    log_fbnb_0 <- lbeta(phi + r, gamma_t) - lbeta(phi, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    
    log_like <- log_like +
      indicator_y_zero * log(pi_t[t]) +
      (1 - indicator_y_zero) * (log(1 - pi_t[t]) +
                                  log_fbnb_y - log(1 - fbnb_0))
  }
  return(log_like)
}

# ==============================================================================
# SECTION 5: MCMC Function — Two-Phase Adaptive (Chen & So, 2006)
#
#  ALIGNED WITH BNB-INGARCHX BASELINE:
#    N = 20000, burn_in = 8000  (default)
#    Phase 2 warmup: burn_in_samples[-1:-1000, ] — drops first 1000 of burn-in
#                    (proportional to baseline's -1:-1000 with burn_in = 8000)
#
#  ARGUMENT: step_sizes — named list:
#    alpha  : length-3 vector  (alpha0, alpha1, alpha2)
#    beta   : length-3 vector  (beta0, beta1, beta2)
#    r      : scalar
#    phi    : scalar
#    omega  : scalar or length-k vector (used when k > 0)
# ==============================================================================

default_step_sizes <- list(
  alpha = c(0.10, 0.10, 0.10),
  beta  = c(0.05, 0.05, 0.50),
  r     = 0.80,
  phi   = 0.60,
  omega = 0.50
)

run_mcmc_empirical <- function(Y, X_list, k, N = 20000, burn_in = 8000,
                               b0_mcmc = 3,
                               step_sizes = default_step_sizes,
                               prior_hyp = list(c1 = 1, c2 = 1,
                                                a1 = 3, a2 = 1,
                                                d1 = 3, d2 = 1)) {
  c1 <- prior_hyp$c1; c2 <- prior_hyp$c2
  a1 <- prior_hyp$a1; a2 <- prior_hyp$a2
  d1 <- prior_hyp$d1; d2 <- prior_hyp$d2
  
  # --- Unpack step sizes from argument ---
  step_size_alpha <- step_sizes$alpha                          # length 3
  step_size_beta  <- step_sizes$beta                           # length 3
  step_size_r     <- step_sizes$r                              # scalar
  step_size_phi   <- step_sizes$phi                            # scalar
  step_size_omega <- if (k > 0) {
    ov <- step_sizes$omega
    if (length(ov) >= k) ov[1:k] else rep(ov, length.out = k)
  } else NULL
  
  # --- Initial values ---
  alpha_current <- c(-2.0, -0.4, 0.25)
  beta_current  <- c(0.5, 0.1, 0.1)
  omega_current <- if (k > 0) rep(0.1, k) else NULL
  r_current     <- 2
  phi_current   <- 3
  b             <- if (k > 0) rep(1, k) else NULL
  
  # --- Storage ---
  # Continuous: alpha(3) + beta(3) + omega(k) + r(1) + phi(1) = 8 + k
  # Discrete:   b(k)
  n_param         <- 3 + 3 + k + 1 + 1 + k
  samples         <- matrix(0, nrow = N - burn_in, ncol = n_param)
  burn_in_samples <- matrix(0, nrow = burn_in,     ncol = n_param)
  
  alpha_accept_count <- 0
  beta_accept_count  <- 0
  omega_accept_count <- 0
  r_accept_count     <- 0
  phi_accept_count   <- 0
  
  # --- Index maps ---
  idx_alpha <- 1:3
  idx_beta  <- 4:6
  idx_omega <- if (k > 0) (7):(6 + k)                else integer(0)
  idx_r     <- 6 + k + 1
  idx_phi   <- 6 + k + 2
  idx_b     <- if (k > 0) (6 + k + 3):(6 + 2*k + 2) else integer(0)
  
  pack_params <- function() {
    c(alpha_current, beta_current,
      if (k > 0) omega_current else numeric(0),
      r_current, phi_current,
      if (k > 0) b else numeric(0))
  }
  
  # ============================================================
  # PHASE 1: Random Walk Metropolis-Hastings (Burn-in)
  # ============================================================
  for (iter in 1:burn_in) {
    alpha_proposal <- numeric(3)
    beta_proposal  <- numeric(3)
    omega_proposal <- if (k > 0) numeric(k) else NULL
    
    if (iter %% 2000 == 0) cat("  [k=", k, "] Burn-in Iteration:", iter, "\n")
    
    # --- alpha ---
    repeat {
      alpha_proposal[1] <- alpha_current[1] + rnorm(1, 0, step_size_alpha[1])
      alpha_proposal[2] <- alpha_current[2] + rnorm(1, 0, step_size_alpha[2])
      alpha_proposal[3] <- alpha_current[3] + rnorm(1, 0, step_size_alpha[3])
      if (check_alpha_constraints(alpha_proposal)) break
    }
    log_accept_ratio <- log_likelihood(Y, alpha_proposal, beta_current, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k)
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        alpha_current <- alpha_proposal
        alpha_accept_count <- alpha_accept_count + 1
      }
    }
    
    # --- beta ---
    repeat {
      beta_proposal[1] <- beta_current[1] + rnorm(1, 0, step_size_beta[1])
      beta_proposal[2] <- beta_current[2] + rnorm(1, 0, step_size_beta[2])
      beta_proposal[3] <- beta_current[3] + rnorm(1, 0, step_size_beta[3])
      if (check_beta_constraints(beta_proposal)) break
    }
    log_accept_ratio <- log_likelihood(Y, alpha_current, beta_proposal, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k)
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        beta_current <- beta_proposal
        beta_accept_count <- beta_accept_count + 1
      }
    }
    
    # --- omega ---
    if (k > 0) {
      repeat {
        for (i in 1:k) omega_proposal[i] <- omega_current[i] + rnorm(1, 0, step_size_omega[i])
        if (all(omega_proposal > 0)) break
      }
      omega_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_proposal,
                                      r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE))
      omega_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                      r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE))
      log_accept_ratio <- omega_prop_lp - omega_curr_lp
      accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
      if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
        if (U < accept_prob) {
          omega_current <- omega_proposal
          omega_accept_count <- omega_accept_count + 1
        }
      }
    }
    
    # --- r ---
    r_proposal <- r_current
    repeat {
      r_proposal <- r_current + rnorm(1, 0, step_size_r)
      if (r_proposal > 0) break
    }
    r_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                r_proposal, phi_current, X_list, b, b0_mcmc, k) +
      dgamma(r_proposal, shape = a1, rate = a2, log = TRUE)
    r_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                r_current, phi_current, X_list, b, b0_mcmc, k) +
      dgamma(r_current, shape = a1, rate = a2, log = TRUE)
    log_accept_ratio <- r_prop_lp - r_curr_lp
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        r_current <- r_proposal
        r_accept_count <- r_accept_count + 1
      }
    }
    
    # --- phi ---
    phi_proposal <- phi_current
    repeat {
      phi_proposal <- phi_current + rnorm(1, 0, step_size_phi)
      if (phi_proposal > 2) break
    }
    phi_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                  r_current, phi_proposal, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_proposal) - d2 * phi_proposal
    phi_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                  r_current, phi_current, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_current) - d2 * phi_current
    log_accept_ratio <- phi_prop_lp - phi_curr_lp
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        phi_current <- phi_proposal
        phi_accept_count <- phi_accept_count + 1
      }
    }
    
    # --- b (discrete, Gibbs-style) ---
    if (k > 0) {
      for (bi_idx in 1:k) {
        b_temp <- b
        lik_b  <- sapply(1:b0_mcmc, function(j) {
          b_candidate <- b_temp; b_candidate[bi_idx] <- j
          log_likelihood(Y, alpha_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_candidate, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        lik_bm   <- exp(lik_b - max_lik)
        prob     <- lik_bm / sum(lik_bm)
        cum_prob <- cumsum(prob); U <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1
          } else {
            B <- FALSE; I <- 1
            while (!B) {
              if (U > cum_prob[I] && U < cum_prob[I + 1]) {
                b[bi_idx] <- I + 1; B <- TRUE
              } else { I <- I + 1 }
            }
          }
        }
      }
    }
    
    burn_in_samples[iter, ] <- pack_params()
  }
  
  cat("\n  [k=", k, "] Phase 1 Acceptance Rates:\n")
  cat("    alpha:", round(alpha_accept_count / burn_in, 4), "\n")
  cat("    beta: ", round(beta_accept_count  / burn_in, 4), "\n")
  if (k > 0) cat("    omega:", round(omega_accept_count / burn_in, 4), "\n")
  cat("    r:    ", round(r_accept_count   / burn_in, 4), "\n")
  cat("    phi:  ", round(phi_accept_count / burn_in, 4), "\n")
  
  # --- Compute Phase 2 proposal moments from last (burn_in - 1000) burn-in draws ---
  # Aligned with BNB-INGARCHX baseline: drop first 1000, use remaining 7000.
  warmup_smp <- burn_in_samples[-1:-1000, , drop = FALSE]
  
  mu_alpha  <- colMeans(warmup_smp[, idx_alpha, drop = FALSE])
  cov_alpha <- cov(warmup_smp[, idx_alpha, drop = FALSE])
  mu_beta   <- colMeans(warmup_smp[, idx_beta,  drop = FALSE])
  cov_beta  <- cov(warmup_smp[, idx_beta,  drop = FALSE])
  
  if (k > 0) {
    if (k == 1) {
      mu_omega  <- mean(warmup_smp[, idx_omega])
      cov_omega <- var(warmup_smp[,  idx_omega])
    } else {
      mu_omega  <- colMeans(warmup_smp[, idx_omega, drop = FALSE])
      cov_omega <- cov(warmup_smp[,      idx_omega, drop = FALSE])
    }
  }
  
  mu_r    <- mean(warmup_smp[, idx_r])
  cov_r   <- var(warmup_smp[,  idx_r])
  mu_phi  <- mean(warmup_smp[, idx_phi])
  cov_phi <- var(warmup_smp[,  idx_phi])
  
  alpha_accept_count <- 0; beta_accept_count  <- 0
  omega_accept_count <- 0; r_accept_count     <- 0; phi_accept_count <- 0
  
  # ============================================================
  # PHASE 2: Independent Kernel Metropolis-Hastings
  # ============================================================
  for (iter in 1:(N - burn_in)) {
    if (iter %% 2000 == 0) cat("  [k=", k, "] IK Iteration:", iter, "\n")
    
    # --- alpha ---
    alpha_proposal <- numeric(3)
    repeat {
      alpha_proposal <- as.numeric(mvrnorm(1, mu = mu_alpha, Sigma = cov_alpha))
      if (check_alpha_constraints(alpha_proposal)) break
    }
    log_gau_curr <- mvtnorm::dmvnorm(alpha_current,  mean = mu_alpha, sigma = cov_alpha, log = TRUE)
    log_gau_prop <- mvtnorm::dmvnorm(alpha_proposal, mean = mu_alpha, sigma = cov_alpha, log = TRUE)
    log_accept_ratio <- log_likelihood(Y, alpha_proposal, beta_current, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) +
      log_gau_curr -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k) - log_gau_prop
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        alpha_current <- alpha_proposal
        alpha_accept_count <- alpha_accept_count + 1
      }
    }
    
    # --- beta ---
    beta_proposal <- numeric(3)
    repeat {
      beta_proposal <- as.numeric(mvrnorm(1, mu = mu_beta, Sigma = cov_beta))
      if (check_beta_constraints(beta_proposal)) break
    }
    log_gau_curr <- mvtnorm::dmvnorm(beta_current,  mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_gau_prop <- mvtnorm::dmvnorm(beta_proposal, mean = mu_beta, sigma = cov_beta, log = TRUE)
    log_accept_ratio <- log_likelihood(Y, alpha_current, beta_proposal, omega_current,
                                       r_current, phi_current, X_list, b, b0_mcmc, k) +
      log_gau_curr -
      log_likelihood(Y, alpha_current, beta_current, omega_current,
                     r_current, phi_current, X_list, b, b0_mcmc, k) - log_gau_prop
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        beta_current <- beta_proposal
        beta_accept_count <- beta_accept_count + 1
      }
    }
    
    # --- omega ---
    if (k > 0) {
      if (k == 1) {
        omega_proposal <- numeric(1)
        repeat {
          omega_proposal <- rnorm(1, mu_omega, sqrt(cov_omega))
          if (omega_proposal > 0) break
        }
        log_gau_curr <- dnorm(omega_current,  mu_omega, sqrt(cov_omega), log = TRUE)
        log_gau_prop <- dnorm(omega_proposal, mu_omega, sqrt(cov_omega), log = TRUE)
      } else {
        omega_proposal <- numeric(k)
        repeat {
          omega_proposal <- as.numeric(mvrnorm(1, mu = mu_omega, Sigma = cov_omega))
          if (all(omega_proposal > 0)) break
        }
        log_gau_curr <- mvtnorm::dmvnorm(omega_current,  mean = mu_omega, sigma = cov_omega, log = TRUE)
        log_gau_prop <- mvtnorm::dmvnorm(omega_proposal, mean = mu_omega, sigma = cov_omega, log = TRUE)
      }
      omega_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_proposal,
                                      r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_proposal, shape = c1, rate = c2, log = TRUE))
      omega_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                      r_current, phi_current, X_list, b, b0_mcmc, k) +
        sum(dgamma(omega_current, shape = c1, rate = c2, log = TRUE))
      log_accept_ratio <- omega_prop_lp + log_gau_curr - omega_curr_lp - log_gau_prop
      accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
      if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
        if (U < accept_prob) {
          omega_current <- omega_proposal
          omega_accept_count <- omega_accept_count + 1
        }
      }
    }
    
    # --- r ---
    r_proposal <- r_current
    repeat {
      r_proposal <- rnorm(1, mu_r, sqrt(cov_r))
      if (r_proposal > 0) break
    }
    log_gau_curr <- dnorm(r_current,  mu_r, sqrt(cov_r), log = TRUE)
    log_gau_prop <- dnorm(r_proposal, mu_r, sqrt(cov_r), log = TRUE)
    r_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                r_proposal, phi_current, X_list, b, b0_mcmc, k) +
      dgamma(r_proposal, shape = a1, rate = a2, log = TRUE)
    r_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                r_current, phi_current, X_list, b, b0_mcmc, k) +
      dgamma(r_current, shape = a1, rate = a2, log = TRUE)
    log_accept_ratio <- r_prop_lp + log_gau_curr - r_curr_lp - log_gau_prop
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        r_current <- r_proposal
        r_accept_count <- r_accept_count + 1
      }
    }
    
    # --- phi ---
    phi_proposal <- phi_current
    repeat {
      phi_proposal <- rnorm(1, mu_phi, sqrt(cov_phi))
      if (phi_proposal > 2) break
    }
    log_gau_curr <- dnorm(phi_current,  mu_phi, sqrt(cov_phi), log = TRUE)
    log_gau_prop <- dnorm(phi_proposal, mu_phi, sqrt(cov_phi), log = TRUE)
    phi_prop_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                  r_current, phi_proposal, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_proposal) - d2 * phi_proposal
    phi_curr_lp <- log_likelihood(Y, alpha_current, beta_current, omega_current,
                                  r_current, phi_current, X_list, b, b0_mcmc, k) +
      (d1 - 1) * log(phi_current) - d2 * phi_current
    log_accept_ratio <- phi_prop_lp + log_gau_curr - phi_curr_lp - log_gau_prop
    accept_prob <- min(1, exp(log_accept_ratio)); U <- runif(1)
    if (!is.na(log_accept_ratio) && !is.nan(log_accept_ratio)) {
      if (U < accept_prob) {
        phi_current <- phi_proposal
        phi_accept_count <- phi_accept_count + 1
      }
    }
    
    # --- b (discrete) ---
    if (k > 0) {
      for (bi_idx in 1:k) {
        b_temp <- b
        lik_b  <- sapply(1:b0_mcmc, function(j) {
          b_candidate <- b_temp; b_candidate[bi_idx] <- j
          log_likelihood(Y, alpha_current, beta_current, omega_current,
                         r_current, phi_current, X_list, b_candidate, b0_mcmc, k)
        })
        max_lik  <- max(lik_b)
        lik_bm   <- exp(lik_b - max_lik)
        prob     <- lik_bm / sum(lik_bm)
        cum_prob <- cumsum(prob); U <- runif(1)
        if (!is.na(cum_prob[1]) && !is.nan(cum_prob[1])) {
          if (U < cum_prob[1]) {
            b[bi_idx] <- 1
          } else {
            B <- FALSE; I <- 1
            while (!B) {
              if (U > cum_prob[I] && U < cum_prob[I + 1]) {
                b[bi_idx] <- I + 1; B <- TRUE
              } else { I <- I + 1 }
            }
          }
        }
      }
    }
    
    samples[iter, ] <- pack_params()
  }
  
  n_phase2     <- N - burn_in
  n_blocks     <- if (k > 0) 5 else 4
  total_accept <- alpha_accept_count + beta_accept_count + r_accept_count + phi_accept_count
  if (k > 0) total_accept <- total_accept + omega_accept_count
  accept_rate  <- (total_accept / n_blocks) / n_phase2
  
  cat("\n  [k=", k, "] Phase 2 Acceptance Rates:\n")
  cat("    alpha:", round(alpha_accept_count / n_phase2, 4), "\n")
  cat("    beta: ", round(beta_accept_count  / n_phase2, 4), "\n")
  if (k > 0) cat("    omega:", round(omega_accept_count / n_phase2, 4), "\n")
  cat("    r:    ", round(r_accept_count   / n_phase2, 4), "\n")
  cat("    phi:  ", round(phi_accept_count / n_phase2, 4), "\n")
  cat("    Overall:", round(accept_rate, 4), "\n")
  
  return(list(samples   = samples,   n_param   = n_param,   k       = k,
              idx_alpha = idx_alpha, idx_beta  = idx_beta,  idx_omega = idx_omega,
              idx_r     = idx_r,     idx_phi   = idx_phi,   idx_b   = idx_b))
}

# ==============================================================================
# SECTION 6: Per-Model Step Sizes and Model Specs
#
#  Target acceptance rate in Phase 1 RWM: 23-44% per block.
#    - acceptance too high (>60%): step too small, increase it
#    - acceptance too low  (<15%): step too large, reduce it
# ==============================================================================

step_sizes_m1 <- list(
  # Model 1 — No exogenous.
  alpha = c(0.28, 0.08, 0.18),
  beta  = c(0.40, 0.08, 0.12),
  r     = 0.80,
  phi   = 0.60,
  omega = 0.50   # unused for k=0
)

step_sizes_m2 <- list(
  # Model 2 — Rainfall Only.
  alpha = c(0.28, 0.08, 0.18),
  beta  = c(0.40, 0.08, 0.12),
  r     = 0.80,
  phi   = 0.60,
  omega = 0.45   # rainfall, posterior range ~1.5
)

step_sizes_m3 <- list(
  # Model 3 — Mean Temp Only.
  alpha = c(0.28, 0.08, 0.18),
  beta  = c(0.40, 0.08, 0.12),
  r     = 0.80,
  phi   = 0.60,
  omega = 0.35   # temperature, posterior range ~1.2
)

step_sizes_m4 <- list(
  # Model 4 — Rain + Mean Temp.
  alpha = c(0.28, 0.08, 0.18),
  beta  = c(0.40, 0.08, 0.12),
  r     = 0.80,
  phi   = 0.60,
  omega = c(0.45, 0.35)   # rainfall, temperature
)

model_specs <- list(
  list(k = 0, X_list = NULL,                 seed = 42, label = "Model 1: No Exo",
       step_sizes = step_sizes_m1),
  list(k = 1, X_list = list(X1_mat),         seed = 42, label = "Model 2: Rainfall Only",
       step_sizes = step_sizes_m2),
  list(k = 1, X_list = list(X2_mat),         seed = 42, label = "Model 3: Mean Temp Only",
       step_sizes = step_sizes_m3),
  list(k = 2, X_list = list(X1_mat, X2_mat), seed = 42, label = "Model 4: Rain + Mean Temp",
       step_sizes = step_sizes_m4)
)

ncores <- min(4L, max(1L, parallel::detectCores() - 1L))
cl     <- makeCluster(ncores)
registerDoParallel(cl)

cat("Running 4 empirical models on", ncores, "cores...\n")
t0 <- Sys.time()

all_results <- foreach(
  spec      = model_specs,
  .packages = c("mvtnorm", "MASS"),
  .export   = c("run_mcmc_empirical", "log_likelihood",
                "check_alpha_constraints", "check_beta_constraints",
                "rbnb", "rbnb_truncated", "default_step_sizes",
                "Y", "X1_mat", "X2_mat", "n")
) %dopar% {
  set.seed(spec$seed)
  cat("Starting:", spec$label, "\n")
  res <- run_mcmc_empirical(Y = Y, X_list = spec$X_list, k = spec$k,
                            step_sizes = spec$step_sizes)
  cat("Completed:", spec$label, "\n")
  res
}

stopCluster(cl)
cat("All models completed in",
    round(difftime(Sys.time(), t0, units = "mins"), 2), "minutes.\n")

res_m1 <- all_results[[1]]
res_m2 <- all_results[[2]]
res_m3 <- all_results[[3]]
res_m4 <- all_results[[4]]

# ==============================================================================
# SECTION 7: Summary Function
# ==============================================================================
mode_fn <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

summarize_empirical <- function(res, param_names, b_names = NULL) {
  samples <- res$samples
  k       <- res$k
  n_cont  <- 3 + 3 + k + 1 + 1   # alpha(3) + beta(3) + omega(k) + r + phi
  
  mean_est   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, mean),   4)
  median_est <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, median), 4)
  sd_est     <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, sd),     4)
  p025_est   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, quantile, probs = 0.025), 4)
  p975_est   <- round(apply(samples[, 1:n_cont, drop = FALSE], 2, quantile, probs = 0.975), 4)
  
  result <- data.frame(
    Parameter = param_names,
    Mean      = mean_est,
    Median    = median_est,
    Std       = sd_est,
    P0.025    = p025_est,
    P0.975    = p975_est,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  if (k > 0 && !is.null(b_names)) {
    for (i in 1:k) {
      mode_b <- mode_fn(samples[, res$idx_b[i]])
      result <- rbind(result, data.frame(
        Parameter = b_names[i], Mean = mode_b,
        Median = "", Std = "", P0.025 = "", P0.975 = "",
        stringsAsFactors = FALSE))
    }
  }
  return(result)
}

# ==============================================================================
# SECTION 8: Geweke Z-Score and Inefficiency Factor
# ==============================================================================
compute_geweke_IF <- function(res, param_names, model_label) {
  samples  <- res$samples
  k        <- res$k
  n_cont   <- 3 + 3 + k + 1 + 1
  cont_smp <- samples[, 1:n_cont, drop = FALSE]
  
  geweke_z <- numeric(n_cont)
  IF_vals  <- numeric(n_cont)
  
  for (i in 1:n_cont) {
    chain_i     <- mcmc(cont_smp[, i])
    geweke_z[i] <- round(geweke.diag(chain_i, frac1 = 0.1, frac2 = 0.5)$z, 4)
    chain_var   <- var(cont_smp[, i])
    spec0       <- spectrum0.ar(chain_i)$spec
    IF_vals[i]  <- round(spec0 / chain_var, 4)
  }
  
  diag_table <- data.frame(
    Parameter = param_names,
    Geweke_Z  = geweke_z,
    IF        = IF_vals,
    Converged = ifelse(abs(geweke_z) <= 2, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  cat("\n===", model_label, "— Geweke Z-Score and Inefficiency Factor ===\n")
  print(diag_table)
  cat("Convergence:", sum(abs(geweke_z) <= 2), "of", n_cont,
      "parameters within [-2, 2].\n")
  cat("IF range: [", min(IF_vals), ",", max(IF_vals), "]\n")
  cat("IF mean: ",   round(mean(IF_vals), 4), "\n")
  
  return(diag_table)
}

# ==============================================================================
# SECTION 9: DIC Computation
# ==============================================================================
compute_DIC <- function(samples, Y, X_list, b0, k,
                        idx_alpha, idx_beta, idx_omega,
                        idx_r, idx_phi, idx_b) {
  n_samples <- nrow(samples)
  
  alpha_bar <- colMeans(samples[, idx_alpha, drop = FALSE])
  beta_bar  <- colMeans(samples[, idx_beta,  drop = FALSE])
  r_bar     <- mean(samples[, idx_r])
  phi_bar   <- mean(samples[, idx_phi])
  omega_bar <- if (k > 0) colMeans(samples[, idx_omega, drop = FALSE]) else NULL
  b_bar     <- if (k > 0) sapply(1:k, function(i) mode_fn(samples[, idx_b[i]])) else NULL
  
  llik_bar <- log_likelihood(Y, alpha_bar, beta_bar, omega_bar,
                             r_bar, phi_bar, X_list, b_bar, b0, k)
  D_bar    <- -2 * llik_bar
  
  llik_smp <- sapply(1:n_samples, function(s) {
    alpha_s <- samples[s, idx_alpha]
    beta_s  <- samples[s, idx_beta]
    r_s     <- samples[s, idx_r]
    phi_s   <- samples[s, idx_phi]
    omega_s <- if (k > 0) samples[s, idx_omega] else NULL
    b_s     <- if (k > 0) sapply(1:k, function(i) samples[s, idx_b[i]]) else NULL
    log_likelihood(Y, alpha_s, beta_s, omega_s, r_s, phi_s, X_list, b_s, b0, k)
  })
  
  E_D  <- mean(-2 * llik_smp)
  pD   <- E_D - D_bar
  DIC  <- D_bar + 2 * pD
  
  return(list(DIC   = round(DIC,   4),
              pD    = round(pD,    4),
              D_bar = round(D_bar, 4),
              E_D   = round(E_D,   4)))
}

# ==============================================================================
# SECTION 10: In-Sample Prediction, RMSE, Bias, and Residual Diagnostics
#
#  RMSE = sqrt( mean( (Y_t - mu_t)^2 ) )   over t = b0+1 .. n
#  Bias = mean( mu_t - Y_t )               over t = b0+1 .. n
# ==============================================================================
compute_predicted <- function(res, Y, X_list, b0_mcmc = 3,
                              model_label, year_start = 2014) {
  samples   <- res$samples
  k         <- res$k
  
  alpha_bar <- colMeans(samples[, res$idx_alpha, drop = FALSE])
  beta_bar  <- colMeans(samples[, res$idx_beta,  drop = FALSE])
  r_bar     <- mean(samples[, res$idx_r])
  phi_bar   <- mean(samples[, res$idx_phi])
  omega_bar <- if (k > 0) colMeans(samples[, res$idx_omega, drop = FALSE]) else NULL
  b_bar     <- if (k > 0) sapply(1:k, function(i) mode_fn(samples[, res$idx_b[i]])) else NULL
  
  n        <- length(Y)
  pi_t     <- rep(0.08, n)
  lambda_t <- rep(mean(Y[Y > 0]), n)
  mu_t     <- numeric(n)
  
  for (t in (b0_mcmc + 1):n) {
    logit_pi_prev <- log(pi_t[t - 1] / (1 - pi_t[t - 1]))
    exponent_pi   <- alpha_bar[1] + alpha_bar[2] * Y[t - 1] + alpha_bar[3] * logit_pi_prev
    pi_t[t]       <- exp(exponent_pi) / (1 + exp(exponent_pi))
    
    lambda_t[t] <- beta_bar[1] + beta_bar[2] * Y[t - 1] + beta_bar[3] * lambda_t[t - 1]
    if (k >= 1) {
      for (i in 1:k) {
        lambda_t[t] <- lambda_t[t] + omega_bar[i] * X_list[[i]][1, t - b_bar[i]]
      }
    }
    lambda_t[t] <- max(lambda_t[t], 1e-6)
    
    gamma_t    <- (phi_bar - 1) / r_bar * lambda_t[t]
    log_fbnb_0 <- lbeta(phi_bar + r_bar, gamma_t) - lbeta(phi_bar, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    mu_t[t]    <- (1 - pi_t[t]) * lambda_t[t] / (1 - fbnb_0)
  }
  
  # --- RMSE and Bias (computed over the modelled window only) ---
  valid_idx <- (b0_mcmc + 1):n
  rmse      <- sqrt(mean((Y[valid_idx] - mu_t[valid_idx])^2))
  bias      <- mean(mu_t[valid_idx] - Y[valid_idx])
  
  cat("  [", model_label, "] RMSE:", round(rmse, 4),
      "  Bias:", round(bias, 4), "\n")
  
  # --- Standardised residuals ---
  residuals <- Y - mu_t
  var_t     <- numeric(n)
  
  for (t in (b0_mcmc + 1):n) {
    gamma_t    <- (phi_bar - 1) / r_bar * lambda_t[t]
    log_fbnb_0 <- lbeta(phi_bar + r_bar, gamma_t) - lbeta(phi_bar, gamma_t)
    fbnb_0     <- exp(log_fbnb_0)
    denom      <- 1 - fbnb_0
    p_pos      <- 1 - pi_t[t]
    
    second_moment_pos <-
      ((phi_bar + r_bar - 1) / (phi_bar - 2)) *
      lambda_t[t] * (1 + lambda_t[t] / r_bar) +
      lambda_t[t]^2
    
    mean_hurdle <- p_pos * lambda_t[t] / denom
    var_t[t]    <- p_pos * second_moment_pos / denom - mean_hurdle^2
    var_t[t]    <- max(var_t[t], 1e-6)
  }
  
  std_resid              <- residuals / sqrt(var_t)
  std_resid[1:b0_mcmc]  <- NA
  
  time_seq <- seq(from = as.Date(paste0(year_start, "-01-01")),
                  by = "week", length.out = n)
  year_breaks <- seq(from = as.Date(paste0(year_start, "-01-01")),
                     to   = max(time_seq), by = "3 years")
  
  safe_label <- gsub("[: ]", "_", model_label)
  pdf(paste0("prediction_diagnostics_CenLuz_", safe_label, ".pdf"),
      width = 12, height = 8)
  
  layout(matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(2, 1.2))
  
  par(mar = c(4, 4, 3, 2))
  plot(time_seq, Y, type = "l", lty = 2, col = "blue",
       xlab = "", ylab = "Cases",
       main = paste0("Central Luzon Leptospirosis — BNB Hurdle-INGARCHX\n",
                     model_label,
                     "   RMSE = ", round(rmse, 2),
                     "   Bias = ", round(bias, 2)),
       cex.main = 0.95, xaxt = "n")
  axis.Date(1, at = year_breaks, labels = format(year_breaks, "%Y"),
            las = 1, cex.axis = 0.85)
  lines(time_seq, mu_t, col = "red", lwd = 1.5)
  legend("topleft",
         legend = c("Observed", "Predicted"),
         lty = c(2, 1), col = c("blue", "red"), bty = "n", cex = 0.85)
  
  par(mar = c(4, 4, 3, 1))
  plot(std_resid, type = "l", col = "black",
       main = "Standardized Residuals",
       xlab = "Time", ylab = "Residual", cex.main = 0.95)
  abline(h = 0, col = "red", lwd = 1.5)
  
  par(mar = c(4, 4, 3, 1))
  acf(na.omit(std_resid), main = "ACF of Residuals",
      col = "black", cex.main = 0.95)
  
  par(mar = c(4, 4, 3, 1))
  acf(na.omit(std_resid)^2, main = "ACF of Squared Residuals",
      col = "black", cex.main = 0.95)
  
  dev.off()
  cat("Prediction plot saved: CenLuz_", safe_label, "\n")
  
  return(list(Y_pred    = mu_t,
              std_resid = std_resid,
              lambda_t  = lambda_t,
              pi_t      = pi_t,
              rmse      = rmse,
              bias      = bias))
}

# ==============================================================================
# SECTION 11: Traceplots and ACF Plots — Symbol-form parameter labels
#
#  R plotmath expressions used for axis titles, keyed by param string name.
#  Passed directly as `main` via the expr_map lookup.
#
#  Supported names: alpha0, alpha1, alpha2, beta0, beta1, beta2,
#                   omega1, omega2, r, phi
# ==============================================================================

# Named lookup table of plotmath call objects.
# Each entry is the [[1]] element of an expression(), which is what
# base graphics accepts when passed as `main`.
param_expr_map <- list(
  alpha0 = expression(alpha[0])[[1]],
  alpha1 = expression(alpha[1])[[1]],
  alpha2 = expression(alpha[2])[[1]],
  beta0  = expression(beta[0])[[1]],
  beta1  = expression(beta[1])[[1]],
  beta2  = expression(beta[2])[[1]],
  omega1 = expression(omega[1])[[1]],
  omega2 = expression(omega[2])[[1]],
  r      = expression(r)[[1]],
  phi    = expression(phi)[[1]]
)

# Helper: return the expression object for a given param name string,
# falling back to a plain character label if not in the map.
get_param_expr <- function(name) {
  if (name %in% names(param_expr_map)) {
    param_expr_map[[name]]
  } else {
    as.name(name)
  }
}

plot_empirical <- function(res, param_names, model_label) {
  n_cont  <- length(param_names)
  samples <- res$samples
  nr      <- ceiling(n_cont / 3)
  
  # --- Traceplots ---
  pdf(paste0("traceplot_cenluz_", model_label, ".pdf"))
  par(mfrow = c(nr, 3))
  for (i in 1:n_cont) {
    lbl <- get_param_expr(param_names[i])
    plot(samples[, i], type = "l",
         main = lbl, ylab = "Value", xlab = "Iteration",
         cex.main = 1.2)
  }
  dev.off()
  
  # --- ACF plots ---
  pdf(paste0("acf_cenluz_", model_label, ".pdf"))
  par(mfrow = c(nr, 3))
  for (i in 1:n_cont) {
    lbl <- get_param_expr(param_names[i])
    acf(samples[, i], main = lbl, cex.main = 1.2)
  }
  dev.off()
  
  cat("Plots saved for CenLuz", model_label, "\n")
}

# ==============================================================================
# SECTION 12: Run All Summaries, Diagnostics, Predictions
# ==============================================================================
param_m1 <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2", "r", "phi")
param_m2 <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2", "omega1", "r", "phi")
param_m3 <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2", "omega1", "r", "phi")
param_m4 <- c("alpha0", "alpha1", "alpha2", "beta0", "beta1", "beta2",
              "omega1", "omega2", "r", "phi")

summary_m1 <- summarize_empirical(res_m1, param_m1)
summary_m2 <- summarize_empirical(res_m2, param_m2, b_names = "b1")
summary_m3 <- summarize_empirical(res_m3, param_m3, b_names = "b1")
summary_m4 <- summarize_empirical(res_m4, param_m4, b_names = c("b1", "b2"))

cat("\n========== Model 1 Summary ==========\n"); print(summary_m1)
cat("\n========== Model 2 Summary ==========\n"); print(summary_m2)
cat("\n========== Model 3 Summary ==========\n"); print(summary_m3)
cat("\n========== Model 4 Summary ==========\n"); print(summary_m4)

diag_m1 <- compute_geweke_IF(res_m1, param_m1, "Model 1: No Exo")
diag_m2 <- compute_geweke_IF(res_m2, param_m2, "Model 2: Rainfall Only")
diag_m3 <- compute_geweke_IF(res_m3, param_m3, "Model 3: Mean Temp Only")
diag_m4 <- compute_geweke_IF(res_m4, param_m4, "Model 4: Rain + Mean Temp")

cat("\n========== DIC Comparison ==========\n")

dic_m1 <- compute_DIC(res_m1$samples, Y, NULL,        3, 0,
                      res_m1$idx_alpha, res_m1$idx_beta, res_m1$idx_omega,
                      res_m1$idx_r, res_m1$idx_phi, res_m1$idx_b)

dic_m2 <- compute_DIC(res_m2$samples, Y, list(X1_mat), 3, 1,
                      res_m2$idx_alpha, res_m2$idx_beta, res_m2$idx_omega,
                      res_m2$idx_r, res_m2$idx_phi, res_m2$idx_b)

dic_m3 <- compute_DIC(res_m3$samples, Y, list(X2_mat), 3, 1,
                      res_m3$idx_alpha, res_m3$idx_beta, res_m3$idx_omega,
                      res_m3$idx_r, res_m3$idx_phi, res_m3$idx_b)

dic_m4 <- compute_DIC(res_m4$samples, Y, list(X1_mat, X2_mat), 3, 2,
                      res_m4$idx_alpha, res_m4$idx_beta, res_m4$idx_omega,
                      res_m4$idx_r, res_m4$idx_phi, res_m4$idx_b)

dic_table <- data.frame(
  Model = c("Model 1: No Exo",
            "Model 2: Rainfall Only",
            "Model 3: Mean Temp Only",
            "Model 4: Rain + Mean Temp"),
  D_bar = c(dic_m1$D_bar, dic_m2$D_bar, dic_m3$D_bar, dic_m4$D_bar),
  E_D   = c(dic_m1$E_D,   dic_m2$E_D,   dic_m3$E_D,   dic_m4$E_D),
  pD    = c(dic_m1$pD,    dic_m2$pD,    dic_m3$pD,    dic_m4$pD),
  DIC   = c(dic_m1$DIC,   dic_m2$DIC,   dic_m3$DIC,   dic_m4$DIC)
)

print(dic_table)
cat("\nBest model (lowest DIC):",
    dic_table$Model[which.min(dic_table$DIC)], "\n")

# --- Predictions with RMSE and Bias ---
cat("\n========== In-Sample RMSE and Bias ==========\n")

pred_m1 <- compute_predicted(res_m1, Y, NULL,
                             model_label = "Model1_No_Exo",
                             year_start  = 2014)
pred_m2 <- compute_predicted(res_m2, Y, list(X1_mat),
                             model_label = "Model2_Rainfall_Only",
                             year_start  = 2014)
pred_m3 <- compute_predicted(res_m3, Y, list(X2_mat),
                             model_label = "Model3_MeanTemp_Only",
                             year_start  = 2014)
pred_m4 <- compute_predicted(res_m4, Y, list(X1_mat, X2_mat),
                             model_label = "Model4_Rain_MeanTemp",
                             year_start  = 2014)

# Consolidated RMSE / Bias table
pred_summary <- data.frame(
  Model = c("Model 1: No Exo",
            "Model 2: Rainfall Only",
            "Model 3: Mean Temp Only",
            "Model 4: Rain + Mean Temp"),
  RMSE  = round(c(pred_m1$rmse, pred_m2$rmse, pred_m3$rmse, pred_m4$rmse), 4),
  Bias  = round(c(pred_m1$bias, pred_m2$bias, pred_m3$bias, pred_m4$bias), 4)
)
cat("\n")
print(pred_summary)
cat("\nBest model (lowest RMSE):",
    pred_summary$Model[which.min(pred_summary$RMSE)], "\n")

# --- Traceplots and ACF (symbol labels) ---
plot_empirical(res_m1, param_m1, "model1")
plot_empirical(res_m2, param_m2, "model2")
plot_empirical(res_m3, param_m3, "model3")
plot_empirical(res_m4, param_m4, "model4")

# ==============================================================================
# SECTION 13: Save All Results
# ==============================================================================
write.csv(summary_m1,   "cenluz_result_model1.csv",      row.names = FALSE)
write.csv(summary_m2,   "cenluz_result_model2.csv",      row.names = FALSE)
write.csv(summary_m3,   "cenluz_result_model3.csv",      row.names = FALSE)
write.csv(summary_m4,   "cenluz_result_model4.csv",      row.names = FALSE)
write.csv(dic_table,    "cenluz_DIC_comparison.csv",     row.names = FALSE)
write.csv(diag_m1,      "cenluz_diagnostics_model1.csv", row.names = FALSE)
write.csv(diag_m2,      "cenluz_diagnostics_model2.csv", row.names = FALSE)
write.csv(diag_m3,      "cenluz_diagnostics_model3.csv", row.names = FALSE)
write.csv(diag_m4,      "cenluz_diagnostics_model4.csv", row.names = FALSE)
write.csv(pred_summary, "cenluz_RMSE_Bias.csv",          row.names = FALSE)

cat("\nAll CenLuz Hurdle-INGARCHX results saved.\n")