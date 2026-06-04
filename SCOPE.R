############################################################
# SCOPE: Symmetric COvariance Population Estimator
# A hierarchical Bayesian regression model implemented in Stan
# Version 1.0.0
# Author: Alister W. Graham (2026)
# Repository: https://github.com/A-Graham/SCOPE
#
# PURPOSE:
#   Fits a linear scaling relation between two observed quantities,
#   fully accounting for measurement uncertainties in both variables.
#   The intrinsic population is modelled as a bivariate distribution,
#   ensuring symmetry in (X, Y) during inference.
#   The regression slope is derived from the intrinsic covariance:
#       beta = Cov(X,Y) / Var(X) = rho * sigma_y / sigma_x
#   This corresponds to the population conditional slope E[Y|X], but is 
#   obtained from a symmetric covariance model rather than from 
#   directional regression (which avoids standard attenuation bias).
#
# HOW TO RUN:
#   1. Place this script (SCOPE.R), the Stan file (Scoped.stan), and
#      the data file (Scin.csv) in the same directory.
#   2. Edit 'input_file' in this script if using a different CSV name.
#   3. Run the script from the terminal: Rscript SCOPE.R
#
# OUTPUT:
#   - Posterior summaries printed to the console.
#   - Figure saved to: SCOPE_fit.pdf
#   - Posterior samples saved to: Scout.dat (loadable in R via load())
############################################################


############################################################
# License: MIT License
#
# Copyright (c) 2026 Alister W. Graham
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
############################################################


############################################################
# CITATION
#
# If you use SCOPE in published work, please cite both 
# the following paper that describes the code:
#
# Graham, A. W. (2026),
# "Galaxy morphology dependent (black hole mass)-(velocity
#  dispersion) relations: implications for gravitational
#  wave forecasts and cosmological simulations", 
# Monthly Notices of the Royal Astronomical Society, submitted
# https://doi.org/10.1093/mnras/(to be assigned)
#
# and the SCOPE software release:
#
# Graham, A. W. (2026),
# SCOPE: Symmetric COvariance Population Estimator,
# Zenodo DOI: (to be assigned)
#
# The latest citation information is provided in CITATION.cff
# in the GitHub repository.
############################################################

############################################################
# INITIALISATION & STARTUP MESSAGE
############################################################
cat("\n================================================================================\n")
cat(" SCOPE: Symmetric COvariance Population Estimator\n")
cat(" Version 1.0.0\n")
cat(" Author : Alister W. Graham (2026)\n")
cat(" GitHub : https://github.com/A-Graham/SCOPE\n")
cat("================================================================================\n")
cat("\nInitialising environment and loading packages...\n\n")
cat(" NOTE: Stan must compile the model to C++ on first use in each session.\n")
cat("       This typically takes 1-3 minutes. Subsequent runs in the same\n")
cat("       session will start sampling almost immediately.\n")
cat("--------------------------------------------------------------------------------\n\n")
flush.console()  # Only effective in the Windows RGui interactive console.
                 # No-op on macOS, Linux, or when running via Rscript on Windows.
############################################################

# Plot Customisation
plot_x_label <- expression(Log[10] ~ plain("Velocity Dispersion") ~ (km ~ s^{-1}))
plot_y_label <- expression(Log[10] ~ plain("Black Hole Mass") ~ (M[solar]))

suppressPackageStartupMessages({
  library(rstan)
  library(sn)
  library(MASS)
  library(ellipse)
})

num_chains <- 4  # Number of independent MCMC chains
# Safely detect cores: leave 1 for the OS, cap at num_chains, never go below 1.
# Handles rare cases where detectCores() returns NA on restricted computing clusters.
n_cores <- parallel::detectCores()
n_cores <- if (is.na(n_cores)) 1L else as.integer(n_cores)
options(mc.cores = max(1L, min(n_cores - 1L, num_chains)))

# Note: Increasing the number of cores can speed things up but will increase memory usage.
# If memory errors occur, reduce 'chains' or 'mc.cores'. For example
#   options(mc.cores = min(2, parallel::detectCores()))
rstan_options(auto_write = TRUE)

############################################################
# INITIALISE FILE PATHS AND VERIFY DIRECTORY
############################################################

# If input_file is not already set in the R environment (e.g., via RStudio), 
# check for command-line arguments. Default to "Scin.csv" if neither exist.
if (!exists("input_file", inherits = FALSE)) {
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- if (length(args) > 0) args[1] else "Scin.csv"
}

# Verify required input files are present in the working directory.
for (required_file in c("Scoped.stan", input_file)) {
  if (!file.exists(required_file)) {
    stop(sprintf(
      "Required file '%s' not found in working directory '%s'.\n  Set your working directory to the folder containing SCOPE.R.",
      required_file, getwd()
    ))
  }
}

############################################################
# REQUIRED INPUT FORMAT (CSV)
#
# Required columns:
#   label      : object name (string)
#   Mbh        : log10 black hole mass (Y variable)
#   dMbh       : symmetric 1-sigma uncertainty in Mbh
#   Sig        : log10 velocity dispersion, or other X variable
#   dSig       : 1-sigma uncertainty in Sig (symmetric case)
#
# Optional asymmetric errors (preferred when available):
#   replace the single column dMbh with the following two columns: 
#     dMbh_pos   : +1-sigma uncertainty in Mbh
#     dMbh_neg   : -1-sigma uncertainty in Mbh
#     If dMbh_pos/dMbh_neg are present, dMbh is ignored.
#
#   replace the single column dSig with the following two columns: 
#     dSig_pos   : +1-sigma uncertainty in Sig
#     dSig_neg   : -1-sigma uncertainty in Sig
#     If dSig_pos/dSig_neg are present, dSig is ignored.
#
# All quantities should be in log10 units.
############################################################


## Read in data
population.galaxies <- read.csv(input_file)
# Clean data: remove any completely empty trailing rows (often caused by Excel)
population.galaxies <- population.galaxies[
  rowSums(is.na(population.galaxies)) < ncol(population.galaxies), ]

# Safety check: drop rows with missing values in critical columns
error_cols <- intersect(names(population.galaxies),
                        c("dMbh", "dMbh_pos", "dMbh_neg",
                          "dSig", "dSig_pos", "dSig_neg"))
valid_rows <- complete.cases(population.galaxies[, c("Mbh", "Sig", error_cols)])
if (sum(!valid_rows) > 0) {
  warning(sprintf(
    "Dropped %d row(s) due to missing (NA) values in required columns.",
    sum(!valid_rows)))
  population.galaxies <- population.galaxies[valid_rows, ]
}

# Count the remaining valid objects
Npopulation <- nrow(population.galaxies)

# Validate numeric columns
if (!is.numeric(population.galaxies$Sig)) stop("Column 'Sig' is not numeric. Check your CSV for formatting issues or stray characters.")
if (!is.numeric(population.galaxies$Mbh)) stop("Column 'Mbh' is not numeric. Check your CSV for formatting issues or stray characters.")

X_est <- population.galaxies$Sig
if ("dSig_pos" %in% names(population.galaxies)) {
  if (!"dSig_neg" %in% names(population.galaxies))
    stop("Column 'dSig_pos' found but 'dSig_neg' is missing from the input CSV.")
  X_err_pos <- abs(population.galaxies$dSig_pos)
  X_err_neg <- abs(population.galaxies$dSig_neg)
} else {
  if (!"dSig" %in% names(population.galaxies))
    stop("Neither 'dSig' nor 'dSig_pos'/'dSig_neg' columns found in the input CSV.")
  X_err_pos <- abs(population.galaxies$dSig)
  X_err_neg <- abs(population.galaxies$dSig)
}

Y_est <- population.galaxies$Mbh
if ("dMbh_pos" %in% names(population.galaxies)) {
  if (!"dMbh_neg" %in% names(population.galaxies))
    stop("Column 'dMbh_pos' found but 'dMbh_neg' is missing from the input CSV.")
  Y_err_pos <- abs(population.galaxies$dMbh_pos)
  Y_err_neg <- abs(population.galaxies$dMbh_neg)
} else {
  if (!"dMbh" %in% names(population.galaxies))
    stop("Neither 'dMbh' nor 'dMbh_pos'/'dMbh_neg' columns found in the input CSV.")
  Y_err_pos <- abs(population.galaxies$dMbh)
  Y_err_neg <- abs(population.galaxies$dMbh)
}

# Holdout vectors: set individual entries to 1L to exclude specific
# objects from the likelihood (e.g., for cross-validation).
# Currently, all objects are included (all zeros).
inholdout_X <- rep(0L, Npopulation)
inholdout_Y <- rep(0L, Npopulation)

# Print a quick data summary to the console before Stan starts
n_symmetric_X <- sum(abs(X_err_pos - X_err_neg) < 1e-6)
n_symmetric_Y <- sum(abs(Y_err_pos - Y_err_neg) < 1e-6)
cat(sprintf("Data loaded successfully : N = %d objects\n", Npopulation))
cat(sprintf("Symmetric error bars detected: X = %d/%d, Y = %d/%d\n",
            n_symmetric_X, Npopulation,
            n_symmetric_Y, Npopulation))
cat("Preparing statistical model...\n")

############################################################
# Data-driven weakly informative priors
# Note: priors are derived from the observed X_est values directly.
# X_est is not passed to Stan but is retained in R for prior computation
# and for residual calculations in the observed scatter section below.

x_mean_data  <- mean(X_est)
x_sd_data    <- stats::sd(X_est)

y_mean_data  <- mean(Y_est)
y_sd_data    <- stats::sd(Y_est)

prior_mean_X_mu  <- x_mean_data
prior_mean_X_sd  <- 2 * x_sd_data 

prior_mean_Y_mu  <- y_mean_data
prior_mean_Y_sd  <- 2 * y_sd_data

prior_sd_X_rate  <- 1 / x_sd_data
prior_sd_Y_rate  <- 1 / y_sd_data
############################################################
# test
# prior_mean_X_sd  <- 5 * 2 * x_sd_data   # 5x broader population mean prior
# prior_mean_Y_sd  <- 5 * 2 * y_sd_data   # 5x broader population mean prior
# prior_sd_X_rate  <- 1 / (5 * x_sd_data) # 5x broader intrinsic scatter prior
# prior_sd_Y_rate  <- 1 / (5 * y_sd_data) # 5x broader intrinsic scatter prior
############################################################


############################################################
## Skew-normal approximation for X uncertainties
X_mu <- X_sigma <- X_alpha <- numeric(Npopulation)

# Initialize warning collectors
fallback_X <- character(0)
failure_X  <- character(0)
boundary_X <- character(0)

for (i in seq_len(Npopulation)) {
  # If errors are symmetric, analytically assign normal parameters and skip optimisation
  tol <- 1e-6 * max((X_err_pos[i] + X_err_neg[i]) / 2, 1e-3)
  if (abs(X_err_pos[i] - X_err_neg[i]) < tol) {
    X_mu[i]    <- X_est[i]
    X_sigma[i] <- X_err_pos[i]
    X_alpha[i] <- 0
  } else {
    qs <- c(X_est[i]-3*X_err_neg[i], X_est[i]-2*X_err_neg[i], X_est[i]-1*X_err_neg[i],
            X_est[i], X_est[i]+1*X_err_pos[i], X_est[i]+2*X_err_pos[i], X_est[i]+3*X_err_pos[i])

    # sqdiff_x is defined inside the loop because it forms a closure over 'qs'
    # Use squared differences (^2) as the function name implies
    sqdiff_x <- function(theta) {
      theta[2] <- max(theta[2], 1e-6)
      # suppressWarnings silences internal precision complaints from qsn() 
      suppressWarnings(
        sum((qsn(pnorm(c(-3,-2,-1,0,1,2,3)),
             xi = theta[1], omega = theta[2], alpha = theta[3]) - qs)^2)
      )
    }

    # Give the optimizer a starting nudge so it doesn't get trapped at zero
    start_alpha_X <- if (X_err_pos[i] > X_err_neg[i]) 2 else -2
    xfit <- nlminb(
      start     = c(X_est[i], max((X_err_neg[i]+X_err_pos[i])/2, 1e-6), start_alpha_X),
      objective = sqdiff_x,
      lower     = c(-Inf, 1e-6, -20),
      upper     = c( Inf,  Inf,  20), 
      control   = list(iter.max = 1000, eval.max = 2000)
    )

    if (xfit$convergence != 0 && abs(xfit$par[3]) >= 0.1) {
      failure_X <- c(failure_X, sprintf("  Object %d (%s): code %d", i, population.galaxies$label[i], xfit$convergence))
      X_mu[i]    <- X_est[i]
      X_sigma[i] <- (X_err_pos[i] + X_err_neg[i]) / 2
      X_alpha[i] <- 0
    } else if (abs(xfit$par[3]) < 0.1) {
      fallback_X <- c(fallback_X, sprintf("  Object %d (%s)", i, population.galaxies$label[i]))
      X_mu[i]    <- X_est[i]
      X_sigma[i] <- (X_err_pos[i] + X_err_neg[i]) / 2
      X_alpha[i] <- 0
    } else {
      xpar <- xfit$par
      if (abs(xpar[3]) > 15) {
        boundary_X <- c(boundary_X, sprintf("  Object %d (%s): alpha = %.2f", i, population.galaxies$label[i], xpar[3]))
      }
      X_mu[i]    <- xpar[1]
      X_sigma[i] <- xpar[2]
      X_alpha[i] <- xpar[3]
    }
  }
}

## Skew-normal approximation for Y uncertainties
Y_mu <- Y_sigma <- Y_alpha <- numeric(Npopulation)

# Initialize warning collectors
fallback_Y <- character(0)
failure_Y  <- character(0)
boundary_Y <- character(0)

for (i in seq_len(Npopulation)) {
  # If errors are symmetric, analytically assign normal parameters and skip optimisation
  tol <- 1e-6 * max((Y_err_pos[i] + Y_err_neg[i]) / 2, 1e-3)
  if (abs(Y_err_pos[i] - Y_err_neg[i]) < tol) {
    Y_mu[i]    <- Y_est[i]
    Y_sigma[i] <- Y_err_pos[i]
    Y_alpha[i] <- 0
  } else {
    qs <- c(Y_est[i]-3*Y_err_neg[i], Y_est[i]-2*Y_err_neg[i], Y_est[i]-1*Y_err_neg[i],
            Y_est[i], Y_est[i]+1*Y_err_pos[i], Y_est[i]+2*Y_err_pos[i], Y_est[i]+3*Y_err_pos[i])

    # sqdiff_y is defined inside the loop because it forms a closure over 'qs'
    sqdiff_y <- function(theta) {
      theta[2] <- max(theta[2], 1e-6)
      # suppressWarnings silences internal precision complaints from qsn() 
      suppressWarnings(
        sum((qsn(pnorm(c(-3,-2,-1,0,1,2,3)),
                    xi = theta[1], omega = theta[2], alpha = theta[3]) - qs)^2)
      )
    }

    start_alpha_Y <- if (Y_err_pos[i] > Y_err_neg[i]) 2 else -2

    yfit <- nlminb(
      start     = c(Y_est[i], max((Y_err_neg[i]+Y_err_pos[i])/2, 1e-6), start_alpha_Y),
      objective = sqdiff_y,
      lower     = c(-Inf, 1e-6, -20),
      upper     = c( Inf,  Inf,  20), 
      control   = list(iter.max = 1000, eval.max = 2000)
    )

    if (yfit$convergence != 0 && abs(yfit$par[3]) >= 0.1) {
      failure_Y <- c(failure_Y, sprintf("  Object %d (%s): code %d", i, population.galaxies$label[i], yfit$convergence))
      Y_mu[i]    <- Y_est[i]
      Y_sigma[i] <- (Y_err_pos[i] + Y_err_neg[i]) / 2
      Y_alpha[i] <- 0
    } else if (abs(yfit$par[3]) < 0.1) {
      fallback_Y <- c(fallback_Y, sprintf("  Object %d (%s)", i, population.galaxies$label[i]))
      Y_mu[i]    <- Y_est[i]
      Y_sigma[i] <- (Y_err_pos[i] + Y_err_neg[i]) / 2
      Y_alpha[i] <- 0
    } else {
      ypar <- yfit$par
      if (abs(ypar[3]) > 15) {
        boundary_Y <- c(boundary_Y, sprintf("  Object %d (%s): alpha = %.2f", i, population.galaxies$label[i], ypar[3]))
      }
      Y_mu[i]    <- ypar[1]
      Y_sigma[i] <- ypar[2]
      Y_alpha[i] <- ypar[3]
    }
  }
}

# Print consolidated notes and warnings if any exist
if (length(fallback_X) > 0) message(sprintf("Note: %d object(s) yielded a near-zero fitted skewness parameter (|alpha| < 0.1) — approximating with a symmetric normal distribution.", length(fallback_X)))
if (length(failure_X) > 0) warning(sprintf("%d object(s) had genuinely asymmetric X errors but the skew-normal fit failed to converge.\n  Falling back to symmetric normal (mean of pos/neg errors used).\n  Check these objects manually:\n%s", length(failure_X), paste(failure_X, collapse = "\n")))
if (length(boundary_X) > 0) warning(sprintf("%d object(s) had X alpha near boundary:\n%s", length(boundary_X), paste(boundary_X, collapse = "\n")))
if (length(fallback_Y) > 0) message(sprintf("Note: %d object(s) yielded a near-zero fitted skewness parameter (|alpha| < 0.1) — approximating with a symmetric normal distribution.", length(fallback_Y)))
if (length(failure_Y) > 0) warning(sprintf("%d object(s) had genuinely asymmetric Y errors but the skew-normal fit failed to converge.\n  Falling back to symmetric normal (mean of pos/neg errors used).\n  Check these objects manually:\n%s", length(failure_Y), paste(failure_Y, collapse = "\n")))
if (length(boundary_Y) > 0) warning(sprintf("%d object(s) had Y alpha near boundary:\n%s", length(boundary_Y), paste(boundary_Y, collapse = "\n")))
############################################################


############################################################
stan_data <- list(
  Npopulation = Npopulation,
  X_mu        = X_mu,
  X_sigma     = X_sigma,
  X_alpha     = X_alpha,
  Y_mu        = Y_mu,
  Y_sigma     = Y_sigma,
  Y_alpha     = Y_alpha,
  inholdout_X = inholdout_X,
  inholdout_Y = inholdout_Y,

  prior_mean_X_mu = prior_mean_X_mu,
  prior_mean_X_sd = prior_mean_X_sd,
  prior_mean_Y_mu = prior_mean_Y_mu,
  prior_mean_Y_sd = prior_mean_Y_sd,
  prior_sd_X_rate = prior_sd_X_rate,
  prior_sd_Y_rate = prior_sd_Y_rate
)

# Note: X_est is intentionally absent from stan_data. The observed X values
# are used in R only, for computing priors and observed residuals. Stan
# receives the skew-normal parameters (X_mu, X_sigma, X_alpha) which fully
# describe the X measurement uncertainty distribution for each object.
###################################################


###################################################
# Updated inits() for Non-Centered Parameterisation

inits <- function() {
  list(
    # true_masses_raw is standard normal, so initializing at 0 is optimal
    true_masses_raw = matrix(0, nrow = Npopulation, ncol = 2),
    sigma = c(stats::sd(X_est), stats::sd(Y_est)),
    Lcorr = matrix(c(1, 0, 0, 1), nrow = 2),  # Cholesky factor of identity (rho = 0)
    mean_X = mean(X_est),
    mean_Y = mean(Y_est)
  )
}
###############################################


##################################################
# Read Stan code and auto-upgrade syntax if needed
##################################################
stan_file <- "Scoped.stan"
stan_code <- readLines(stan_file)

# If the user has a modern version of rstan (>= 2.26.0),
# translate the old array syntax to the new array syntax on the fly.
if (requireNamespace("rstan", quietly = TRUE) &&
    packageVersion("rstan") >= package_version("2.26.0")) {

  stan_code_original <- stan_code
  stan_code <- tryCatch({
    # Apply a generic regex translator to convert 1D array declarations:
    #   TYPE name[SIZE];   -->   array[SIZE] TYPE name;
    # Automatically handles real, int, vector[2], constraints like <lower=0>,
    # and dynamically captures any valid Stan identifier, underscore, or
    # integer literal as the array size (e.g., Npopulation, N, K, 50).
    #
    # NOTE: This translator does not support multi-dimensional arrays
    # (e.g., real foo[N, K]; or real foo[N][K];).
    # Users adding 2D arrays to the Stan file should write them
    # directly in new syntax: array[N, K] real foo;
    sapply(stan_code, function(line) {
      if (!grepl("^\\s*//", line)) {
        line <- gsub(
          pattern = paste0('^(\\s*)([A-Za-z_][A-Za-z0-9_]*',
                           '(?:<[^>]+>)?(?:\\[[^\\]]*\\])?)\\s+',
                           '([A-Za-z_][A-Za-z0-9_]*)\\s*',
                           '\\[\\s*([A-Za-z0-9_]+)\\s*\\]\\s*;'),
          replacement = '\\1array[\\4] \\2 \\3;',
          x = line,
          perl = TRUE
        )
      }
      line
    }, USE.NAMES = FALSE)

  }, error = function(e) {
    warning("Stan syntax auto-translation failed: ", conditionMessage(e),
            "\nProceeding with original syntax. ",
            "Compilation may fail on rstan >= 2.26.0.")
    stan_code_original
  })
}

# Collapse back into a single string
stan_code_string <- paste(stan_code, collapse = "\n")
###############################################


###############################################
# MCMC settings and model fitting
###############################################
cat("\n--------------------------------------------------------------------------------\n")
cat(" Starting MCMC sampling. Stan will print progress every 500 iterations.\n")
cat(" If sampling is too slow, consider:\n")
cat("   - Reducing 'num_iter' and 'num_warmup' in the MCMC settings block\n")
cat("   - Increasing 'mc.cores' if your hardware allows\n")
cat("   - Reducing the dataset size for exploratory runs\n")
cat(" The PDF output will be saved regardless of how long sampling takes.\n")
cat("--------------------------------------------------------------------------------\n")
# Dynamically scale iterations based on dataset size
# Small datasets need more draws to get good Effective Sample Size (ESS)
# Large datasets constrain the posterior well and take longer per step.
if (Npopulation < 100) {
  num_iter <- 4000
  num_warmup <- 2000
} else {
  num_iter <- 2000
  num_warmup <- 1000
}

# Note: We now use `model_code = stan_code_string` instead of `file = ...`
scope_fit <- stan(
  model_code = stan_code_string, 
  data = stan_data,
  init = inits,
  iter = num_iter,           
  warmup = num_warmup,         
  chains = num_chains,
  seed = 1234,
  refresh = 500,  # Prints progress every 500 iterations; reduces scrambled console output
  control = list(
    # Slightly elevated from default 0.8; increase to 0.99 or 0.999 if divergences occur
    adapt_delta   = 0.95,
    # Increase to 15 or 20 if max_treedepth warnings appear
    max_treedepth = 12
  )
)
###############################################


############################################################
# Extract samples and print results
############################################################

# 1. Extract the samples into the variable 'parameters'
parameters <- extract(scope_fit)
n_draws    <- length(parameters$slope)

# Verify true_masses array dimensions: expect [n_draws, Npopulation, 2]
stopifnot(
  length(dim(parameters$true_masses)) == 3,
  dim(parameters$true_masses)[2] == Npopulation,
  dim(parameters$true_masses)[3] == 2
)

# 2. Define a helper function for clean table printing
print_param <- function(name, values) {
  # Calculate statistics
  val_mean   <- mean(values)
  val_sd     <- stats::sd(values)
  # Calculate median and exact 2-sigma Credible Interval (~2.28% - 97.72%)
  val_quant  <- quantile(values, probs = c(pnorm(-2), 0.5, pnorm(2)))
  
  # Print formatted output
  # Format: Name | Mean +/- SD | Median [Lower, Upper]
  cat(sprintf("%-25s : Mean: %6.3f +/- %5.3f | Median: %6.3f [%6.3f, %6.3f]\n", 
              name, val_mean, val_sd, val_quant[2], val_quant[1], val_quant[3]))
}

###############################################
cat("\n\n------------------------------------------------------------------------------------\n")
cat("  *** SCOPE: Symmetric COvariance Population Estimator *** \n\n")
# Run summary: dataset size and posterior sample count.
# Reported here for reproducibility — these numbers could 
# be quoted alongside the posterior estimates in publications.
# ------------------------------------------------------------------
total_post_warmup_draws <- (num_iter - num_warmup) * num_chains
cat(sprintf("Input data   : N = %d objects\n", Npopulation))
cat(sprintf("MCMC samples : %d chains x %d post-warmup draws = %d total\n",
            num_chains, num_iter - num_warmup, total_post_warmup_draws))

# SAMPLER DIAGNOSTICS
#
# Three complementary convergence diagnostics are reported:
#
#   Rhat (potential scale reduction factor):
#     Compares within-chain and between-chain variance across the 4
#     chains. Values close to 1.00 indicate the chains have converged
#     to the same distribution. Rhat > 1.01 for any parameter is a
#     warning sign that the sampler has not fully converged and results
#     should be treated with caution.
#
#   n_eff (effective sample size):
#     Accounts for autocorrelation within chains. Because successive
#     MCMC draws are not independent, n_eff is typically less than the
#     total number of post-warmup draws (here: 4 chains x
#     (num_iter - num_warmup) draws each). A rough rule of thumb is
#     that n_eff > 400 is adequate for reliable posterior summaries;
#     values below ~100 suggest poor mixing and unreliable estimates.
#
#   Divergent transitions:
#     A divergence occurs when the numerical integrator makes an error
#     large enough to leave the typical set of the posterior. Even a
#     small number of divergences can bias results. The target is zero.
#     If divergences are detected, increase adapt_delta toward 1
#     (e.g. 0.99 or 0.999) in the stan() control list.
# ------------------------------------------------------------------
summary_fit <- summary(scope_fit)$summary
max_rhat    <- max(summary_fit[, "Rhat"], na.rm = TRUE)
bad_rhat    <- sum(summary_fit[, "Rhat"] > 1.01, na.rm = TRUE)
min_neff    <- min(summary_fit[, "n_eff"], na.rm = TRUE)

sampler_params <- get_sampler_params(scope_fit, inc_warmup = FALSE)
divergences    <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

cat("\nSampler Diagnostics\n")

# --- Rhat ---
cat(sprintf("  Max Rhat                  : %.3f", max_rhat))
if (max_rhat <= 1.01) {
  cat("  (Good: all chains converged)\n")
} else {
  cat("  WARNING: Rhat > 1.01 — chains may not have converged\n")
  bad_params <- rownames(summary_fit)[summary_fit[, "Rhat"] > 1.01]
  # Exclude the hundreds of latent 'true_masses_raw' from flooding the console
  bad_params_clean <- bad_params[!grepl("true_masses_raw", bad_params)]
  if (length(bad_params_clean) > 0) {
    cat(sprintf("    -> Problematic parameters: %s\n", paste(bad_params_clean, collapse=", ")))
  }
}

cat(sprintf("  Parameters with Rhat>1.01 : %d", bad_rhat))
if (bad_rhat == 0) {
  cat("  (Good: no parameters flagged)\n")
} else {
  cat(sprintf("  WARNING: %d parameter(s) show poor convergence\n", bad_rhat))
}

# --- n_eff ---
cat(sprintf("  Min n_eff                 : %.0f", min_neff))
if (min_neff >= 400) {
  cat("  (Good: adequate effective sample size)\n")
} else if (min_neff >= 100) {
  cat("  (Acceptable, but consider more iterations)\n")
} else {
  cat("  WARNING: n_eff < 100 — estimates may be unreliable\n")
}

# --- Divergences ---
if (divergences > 0) {
  cat(sprintf("  Divergent transitions     : %d  (Target: 0)\n", divergences))
  cat("\n  WARNING: Divergences detected. The posterior may be biased.\n")
  cat("  ACTION : Increase 'adapt_delta' (e.g. 0.99 or 0.999) in the stan() control list.\n")
} else {
  cat(sprintf("  Divergent transitions     : %d  (Good: none detected)\n", divergences))
}
cat("------------------------------------------------------------------------------------\n")

# Save posterior samples to disk only after diagnostics have been processed/flagged
save(parameters, file="Scout.dat")

cat("\n====================================================================================\n")
cat(" Model form: \n")
cat("   Y = mu_y + beta * (X - mu_x) \n")
cat("   beta = rho * sigma_y / sigma_x  \n\n")
cat(" Posterior summaries:\n")
cat("   Mean +/- 1-sigma (posterior standard deviation) \n")
#    (the posterior is the probability distribution of each parameter
#     given the data and the model; analogous to a likelihood surface)
cat("   Median with 2-sigma (~95.4%) credible interval \n")
#     (slope derived from intrinsic population covariance;
#      equivalent to the conditional slope Y|X, corrected for errors in X)
#      A 95% credible interval: the true parameter lies within this range
#     with 95% probability, given the data (Bayesian analog of a
#     frequentist confidence interval, but with a more direct interpretation).
cat("------------------------------------------------------------------------------------\n")

cat("Regression relation: \n")
print_param("Pivot Y (mu_y)", parameters$pivot_y)
print_param("Pivot X (mu_x)", parameters$pivot_x)
print_param("Slope (beta)", parameters$slope)
#   (relation is Y = mu_y + beta*(X - mu_x), anchored at the sample means)
#  print_param("Intercept at X=0 (alpha)", parameters$intercept)
#  cat("(note: highly covariant with slope; use pivot-centred form for prediction) \n")

cat("------------------------------------------------------------------------------------\n")
cat("Population covariance structure: \n")
#   Lcorr is an array [iterations, 2, 2]. 
#   We want the off-diagonal element [row 2, col 1] which represents rho
print_param("Intrinsic population correlation, rho", parameters$rho) 
#   (rho is the correlation between true X and true Y in the population,
#    after removing the contribution of measurement noise)
#    Sigma is a matrix [iterations, 2]. 
#    Column 1 = sigma or Bulge mass (X), Column 2 = BH mass (Y)
print_param("Intrinsic dispersion, sig_x   ", parameters$sigma[,1]) 
print_param("Intrinsic dispersion, sig_y   ", parameters$sigma[,2])
#   (intrinsic 1-sigma spread of the true X and Y population distributions,
#    distinct from the per-object measurement uncertainties in the input data)
print_param("Intrinsic scatter, (sig_y|x)", parameters$intrinsic_scatter) 
#   (irreducible physical scatter of true Y values around the mean relation,
#    after accounting for all measurement uncertainties; analogous to the
#    intrinsic scatter epsilon in FITEXY)
##       Note: beta = rho*sigma_y/sigma_x holds exactly per posterior draw,
##       but not when combining the marginal means/medians reported above,
##       due to posterior correlation between rho, sigma_x, and sigma_y.
#   (note: the reported slope is the primary result; rho, sigma_x, and sigma_y
#   describe the population structure and cannot be combined arithmetically
#   to reproduce the slope from their individual reported means or medians). 


cat("------------------------------------------------------------------------------------\n")
cat("Observed vertical scatter about fitted relation:\n")
# cat("(Observed scatter includes measurement uncertainties; compare with Intrinsic Scatter above) \n")
# Posterior means
mu_x  <- mean(parameters$pivot_x)
mu_y  <- mean(parameters$pivot_y)
beta  <- mean(parameters$slope)
# Observed data
x_obs <- X_est
y_obs <- Y_est
# Fitted conditional mean
y_hat <- mu_y + beta * (x_obs - mu_x)
# Residuals
residuals <- y_obs - y_hat
# RMS vertical scatter (population definition)
rms_vertical <- sqrt(mean(residuals^2))
# Sample SD version (division by N-1)
rms_vertical_sd <- stats::sd(residuals)
cat(sprintf("Classical RMS (division by N)        = %6.3f (point estimate only)\n", rms_vertical))
cat(sprintf("Sample Standard Deviation (division by N-1) = %6.3f (point estimate only)\n", rms_vertical_sd))

rms_draws <- numeric(length(parameters$slope))
for (i in seq_len(n_draws)) {
  y_hat_i <- parameters$pivot_y[i] +
             parameters$slope[i] *
             (X_est - parameters$pivot_x[i])
  rms_draws[i] <- sqrt(mean((Y_est - y_hat_i)^2))
}
cat(sprintf("Observed RMS (posterior mean), Delta_rms = %6.3f +/- %6.3f\n",
            mean(rms_draws), stats::sd(rms_draws)))

cat("  The posterior RMS is the most complete summary as it propagates\n")
cat("  uncertainty in the fitted relation; report this value in publications.\n")
cat("  (note: the posterior RMS marginalises over uncertainty in the fitted \n")
cat("   relation parameters, and is typically slightly larger than the simple RMS) \n")
cat("====================================================================================\n\n")


############################################################
# PLOTTING
#
# Produces two identical figures:
#   (1) An interactive screen window (held open until user dismisses it)
#   (2) SCOPE_fit.pdf saved to the working directory
#
# Figure contents:
#   - Per-galaxy 1-sigma posterior ellipses from the joint posterior of
#     (true X, true Y), reflecting both observational uncertainty and
#     hierarchical shrinkage
#   - 68% and 95% pointwise credible bands for the posterior mean relation
#   - Posterior median mean relation line
#   - Observed data points with asymmetric error bars on both X and Y
#   - Equation annotation in pivot-centred form with 1-sigma intervals
#
# Note: the screen device does not support semi-transparent colours, so
# draw_plot(use_alpha = FALSE) uses solid opaque approximations for the
# screen render, while draw_plot(use_alpha = TRUE) uses full transparency
# for the PDF.
############################################################

# ------------------------------------------------------------------
# Posterior summaries for annotation
# Mean values are used throughout to match the console output above.
# ~68% credible intervals (quantiles 0.1586553, 0.8413447) correspond to 1-sigma
# for a symmetric posterior.
# ------------------------------------------------------------------
slope_draws <- parameters$slope
pivx_draws  <- parameters$pivot_x
pivy_draws  <- parameters$pivot_y

slope_mean <- mean(slope_draws)
slope_lo   <- quantile(slope_draws, pnorm(-1))
slope_hi   <- quantile(slope_draws, pnorm(1))

pivy_mean  <- mean(pivy_draws)
pivy_lo    <- quantile(pivy_draws, pnorm(-1))
pivy_hi    <- quantile(pivy_draws, pnorm(1))

pivx_mean  <- mean(pivx_draws)
pivx_lo    <- quantile(pivx_draws, pnorm(-1))
pivx_hi    <- quantile(pivx_draws, pnorm(1))
# Pivot-centred relation per draw, re-expressed as intercept at X=0
# for efficient evaluation over the x_grid below:
#   Y = pivot_y + slope*(X - pivot_x)  =>  Y = (pivot_y - slope*pivot_x) + slope*X
int_draws <- pivy_draws - slope_draws * pivx_draws

# ------------------------------------------------------------------
# Data aliases — keeps the draw_plot() function readable and ensures
# the plotting block is self-contained regardless of earlier variable names
# ------------------------------------------------------------------
x_data    <- X_est
y_data    <- Y_est
x_err_pos <- X_err_pos
x_err_neg <- X_err_neg
y_err_pos <- Y_err_pos
y_err_neg <- Y_err_neg

# ------------------------------------------------------------------
# Axis limits: driven tightly by the data and their error bars,
# with a small fixed pad so no point touches the plot border
# ------------------------------------------------------------------
x_pad <- 0.10
y_pad <- 0.20
x_lim <- c(min(x_data - x_err_neg) - x_pad,
            max(x_data + x_err_pos) + x_pad)
y_lim <- c(min(y_data - y_err_neg) - y_pad,
            max(y_data + y_err_pos) + y_pad)
# ------------------------------------------------------------------
# Credible bands for the posterior mean relation
# Evaluated pointwise over a dense x_grid spanning the plot window
# ------------------------------------------------------------------
x_grid     <- seq(x_lim[1], x_lim[2], length.out = 200)
y_mean_mat <- matrix(NA, nrow = n_draws, ncol = length(x_grid))
for (d in seq_len(n_draws)) {
  y_mean_mat[d, ] <- int_draws[d] + slope_draws[d] * x_grid
}
# Exact 1-sigma bounds (approx. 15.87% and 84.13%)
band_lo68 <- apply(y_mean_mat, 2, quantile, pnorm(-1))
band_hi68 <- apply(y_mean_mat, 2, quantile, pnorm(1))

# Exact 2-sigma bounds (approx. 2.28% and 97.72%)
band_lo95 <- apply(y_mean_mat, 2, quantile, pnorm(-2))
band_hi95 <- apply(y_mean_mat, 2, quantile, pnorm(2))

band_mean <- apply(y_mean_mat, 2, mean)
band_med  <- apply(y_mean_mat, 2, median)

# ------------------------------------------------------------------
# Per-galaxy 1-sigma posterior ellipses
#
# For each galaxy i, the joint posterior of (true X_i, true Y_i) is
# summarised by its mean vector and 2x2 covariance matrix, from which
# a 68.3% (1-sigma) bivariate ellipse is drawn.  This is much more
# compact in file size than plotting individual posterior scatter points,
# and cleanly visualises both the size and orientation of each galaxy's
# positional uncertainty in the regression plane.
# ------------------------------------------------------------------
ellipse_list <- vector("list", Npopulation)
for (i in seq_len(Npopulation)) {
  px <- parameters$true_masses[, i, 1]   # posterior draws for true X_i
  py <- parameters$true_masses[, i, 2]   # posterior draws for true Y_i
  ellipse_list[[i]] <- list(
    centre = c(mean(px), mean(py)),
    cov    = cov(cbind(px, py))
  )
}

# ------------------------------------------------------------------
# Equation annotation (built outside draw_plot so it is computed once)
#
# Main line:  log M_bh = <mean pivy> + <mean slope> * (log sigma_0 - <mean pivx>)
# Second line: 1-sigma credible intervals on each parameter
#
# NOTE: The main equation uses posterior means, while the CIs use quantile-based 
# 1-sigma intervals. (See the detailed note at the end of the file regarding 
# this deliberate distinction for skewed posteriors).
#
# plotmath is used for both lines to ensure proper rendering of
# subscripts (sigma_0) and the multiplication symbol (%*%).
# ------------------------------------------------------------------
eq_expr <- bquote(
  log~M[bh] == .(sprintf("%.2f", pivy_mean)) +
    .(sprintf("%.2f", slope_mean)) *
    (log~sigma[0] - .(sprintf("%.2f", pivx_mean)))
)
str_y <- sprintf("[%.2f, %.2f]", pivy_lo, pivy_hi)
str_b <- sprintf("[%.2f, %.2f]", slope_lo, slope_hi)
str_x <- sprintf("[%.2f, %.2f]", pivx_lo, pivx_hi)

ci_expr <- bquote(.(str_y) ~ "+" ~ .(str_b) %*% (log~sigma[0] - .(str_x)))

# ------------------------------------------------------------------
# draw_plot(): renders the complete figure to the current graphics device
#
# Arguments:
#   use_alpha  logical; TRUE for PDF (supports transparency),
#              FALSE for screen (uses opaque colour approximations)
# ------------------------------------------------------------------
draw_plot <- function(use_alpha = TRUE) {

  par(mai = c(0.75, 1.10, 0.15, 0.15), family = "serif",
    cex.axis = 1.4, col.axis = "black")
  plot(NA, xlim = x_lim, ylim = y_lim,
    xlab = "", ylab = "", xaxt = "n", yaxt = "n",
    xaxs = "i", yaxs = "i")


  # --- 1. Per-galaxy 1-sigma posterior ellipses ---
  # Ellipse fill and border colours differ between devices
  ell_fill <- if (use_alpha) rgb(0.60, 0.60, 0.60, 0.22) else rgb(0.84, 0.84, 0.84)
  ell_bord <- if (use_alpha) rgb(0.35, 0.35, 0.35, 0.40) else rgb(0.58, 0.58, 0.58)

  # Exact probability mass for a 1-sigma normal interval (~0.6826895)
  mass_1sigma <- pnorm(1) - pnorm(-1) 

  for (i in seq_len(Npopulation)) {
    ell_pts <- ellipse(ellipse_list[[i]]$cov,
                     centre = ellipse_list[[i]]$centre,
                     level  = mass_1sigma)
    polygon(ell_pts, border = ell_bord, col = ell_fill, lwd = 1.2)
  }


  # --- 2. Credible bands for the posterior mean relation ---
  # 95% band (lighter) drawn first, then 68% band (darker) on top
  if (use_alpha) {
    col95 <- hsv(0.85, s = 0.35, v = 0.95, alpha = 0.30)
    col68 <- hsv(0.85, s = 0.40, v = 0.85, alpha = 0.50)
  } else {
    col95 <- hsv(0.85, s = 0.12, v = 0.93)
    col68 <- hsv(0.85, s = 0.22, v = 0.82)
  }
  polygon(c(x_grid, rev(x_grid)), c(band_lo95, rev(band_hi95)),
        border = NA, col = col95)
  polygon(c(x_grid, rev(x_grid)), c(band_lo68, rev(band_hi68)),
        border = NA, col = col68)


  # --- 3a. Posterior mean relation (solid) ---
  lines(x_grid, band_mean,
      col = hsv(0.85, s = 0.8, v = 0.50), lwd = 3.5)
      
  # --- 3b. Posterior median relation (dashed) ---
  lines(x_grid, band_med,
      col = hsv(0.85, s = 0.8, v = 0.50), lwd = 1.5, lty = 2)


  # --- 4. Observed data points with error bars ---
  # Define transparent colors: rgb(red, green, blue, alpha)
  col_err <- rgb(0.3, 0.3, 0.3, 0.5) # Semi-transparent dark grey for error bars
  col_pt  <- rgb(0.1, 0.1, 0.1, 0.8) # Barely transparent near-black for points
  # Drawn last so they appear on top of the ellipses and bands  
  for (i in seq_len(Npopulation)) {
   # Vertical error bar (asymmetric BH mass uncertainties)
    lines(c(x_data[i], x_data[i]),
          c(y_data[i] - y_err_neg[i], y_data[i] + y_err_pos[i]),
          col = col_err, lwd = 2.0)
    # Horizontal error bar (asymmetric X uncertainties)
    lines(c(x_data[i] - x_err_neg[i], x_data[i] + x_err_pos[i]),
          c(y_data[i], y_data[i]),
          col = col_err, lwd = 2.0)
    points(x_data[i], y_data[i], col = col_pt, pch = 19, cex = 1.4)
  }

  # ---  Subtle Background Grid ---
  #  abline(v = x_maj, h = y_maj, col = "grey90", lty = "dotted", lwd = 1.2)


  # --- 5. Equation annotation ---
  # Positioned in the upper-left corner of the plot area.
  # Main line uses plotmath for proper M_bh and sigma_0 rendering.
  # CI line is plain text positioned just below.
  x_ann <- x_lim[1] + 0.03 * diff(x_lim)
  y_ann <- y_lim[2] - 0.05 * diff(y_lim)
  y_ci  <- y_lim[2] - 0.13 * diff(y_lim)

  # White backing rectangle sized to actual annotation text width
  rect(x_lim[1],
       y_lim[2] - 0.166 * diff(y_lim),
       x_lim[1] + 0.55 * diff(x_lim),
       y_lim[2],
       col = rgb(1, 1, 1, 0.75), border = NA)
       
  text(x_ann, y_ann, eq_expr,
       adj = c(0, 1), cex = 1.4)
  text(x_ann, y_ci, ci_expr,
       adj = c(0, 1), cex = 0.95, col = "grey35")


  # --- 6. Axis tick marks and labels ---
  
  # Generate exact sequence for major X ticks at 0.1 intervals
  x_maj <- seq(floor(x_lim[1] * 10) / 10, ceiling(x_lim[2] * 10) / 10, by = 0.1)
  y_maj <- axTicks(2)  # Keep R's default 'pretty' algorithm for Y major ticks
  
  # Major axes with ticks inside the frame
  axis(1, at = x_maj, tck = 0.02, padj = -0.5, lwd = 2.5)               # Bottom
  axis(2, at = y_maj, tck = 0.02, las = 2, hadj = 0.6, lwd = 2.5)      # Left
  axis(3, at = x_maj, tck = 0.02, lwd = 2.5, labels = FALSE)           # Top
  axis(4, at = y_maj, tck = 0.02, lwd = 2.5, labels = FALSE, las = 2)  # Right
  
  # Minor ticks (sub-ticks) inside the frame
  x_min <- seq(floor(x_lim[1] * 100) / 100, ceiling(x_lim[2] * 100) / 100, by = 0.01)
  y_min <- seq(floor(y_lim[1] * 10) / 10, ceiling(y_lim[2] * 10) / 10, by = 0.1)
  
  axis(1, at = x_min, labels = FALSE, tck = 0.01, lwd = 2.5)  # Bottom minors
  axis(3, at = x_min, labels = FALSE, tck = 0.01, lwd = 2.5)  # Top minors
  axis(2, at = y_min, labels = FALSE, tck = 0.01, lwd = 2.5)  # Left minors
  axis(4, at = y_min, labels = FALSE, tck = 0.01, lwd = 2.5)  # Right minors

  mtext(plot_x_label, side = 1, line = 2.5, cex = 1.6)
  mtext(plot_y_label, side = 2, line = 2.5, cex = 1.6)
}

# ------------------------------------------------------------------
# Render to PDF first — ensures output is saved regardless of what
# happens with the interactive screen window.
# SCOPE_fit.pdf is overwritten on each run; no accumulating files.
# ------------------------------------------------------------------
pdf("SCOPE_fit.pdf", width = 8, height = 6)   # Default 
draw_plot(use_alpha = TRUE)
invisible(dev.off())
cat("Plot saved to SCOPE_fit.pdf\n")
cat("NOTE: The figure's first equation line uses posterior means,\n")
cat("      matching the 'Mean' column in the console output above.\n")
cat("      The figure's second line shows 1-sigma quantile intervals\n")
cat("      [15.9%, 84.1%] centred on the posterior median. These are\n")
cat("      narrower than the 2-sigma intervals [2.3%, 97.7%] reported\n")
cat("      in the console's 'Median' column.\n")

# ------------------------------------------------------------------
# Attempt to render to an interactive screen window.
#
# The appropriate device varies by OS:
#   macOS   — quartz()
#   Linux   — X11()  (requires an X display; may be unavailable on
#                      headless servers or remote sessions)
#   Windows — windows()
#
# If the screen device cannot be opened (e.g. no display available,
# or Rscript was compiled without GUI support), a warning is issued
# and the script continues — the PDF is your figure in that case.
#
# Keeping the window open:
#   A while-loop periodically checks dev.cur() to see if the plot window 
#   is still active. Sys.sleep(0.2) pauses execution briefly to allow 
#   the OS window server to process UI events (like resizing or rendering) 
#   without hogging the CPU or causing OS-level freezes.
#   The script exits gracefully once the user closes the plot window.
# ------------------------------------------------------------------
screen_ok <- tryCatch({
  sysname <- Sys.info()[["sysname"]]
  if (sysname == "Darwin") {
    quartz(width = 8, height = 6)
  } else if (sysname == "Linux") {
    X11(width = 8, height = 6)       # Linux MUST use X11
  } else {
    windows(width = 8, height = 6)   # Windows uses windows
  }
  TRUE
}, error = function(e) {
  message("Note: could not open a screen graphics window on this system.")
  message("      The figure has been saved to SCOPE_fit.pdf instead.")
  FALSE
})

# ------------------------------------------------------------------
if (screen_ok) {
  draw_plot(use_alpha = FALSE)
  cat("\nQuicklook figure displayed!\n")
  cat("--> CLOSE THE FIGURE WINDOW to end the script...\n")
  
  tryCatch({
    # Keep the script alive as long as a graphics device is actively open.
    # Sys.sleep(0.2) allows the OS window manager to process drawing/resizing
    # without hogging the CPU or causing the macOS "spinning beachball".
    while (dev.cur() > 1) {    
      Sys.sleep(0.2)
    }
  }, error = function(e) {
    invisible()
  })
}


# ------------------------------------------------------------------
# NOTE ON FIGURE vs CONSOLE UNCERTAINTY SUMMARIES
#
# The console reports parameter uncertainties as the posterior standard
# deviation (SD): a symmetric +/- value centred on the posterior mean.
#
# The figure shows two central lines:
#   Solid  — posterior mean relation (matches the first annotation line)
#   Dashed — posterior median relation (matches the second annotation line
#             and the shaded credible bands)
#
# The shaded bands and the second annotation line both use the quantile-
# based 1-sigma interval [Phi(-1), Phi(1)] = [15.87%, 84.13%], which is
# centred on the posterior median and makes no symmetry assumption.
#
# For symmetric posteriors the mean and median coincide, the two lines
# are indistinguishable, and the quantile interval equals the SD interval.
# For skewed posteriors they will differ slightly, and the figure makes
# this explicit by showing both lines.
# ------------------------------------------------------------------
