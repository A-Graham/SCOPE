/*
===============================================================================
SCOPE: Symmetric COvariance Population Estimator
Stan model file: Scoped.stan

Version    : 1.0.0
Author     : Alister W. Graham (2026)
Institution: Centre for Astrophysics and Supercomputing,
             Swinburne University of Technology
Repository : https://github.com/A-Graham/SCOPE
License    : MIT

Part of the SCOPE package. This file should be used alongside SCOPE.R
and Scin.csv. See the repository for full documentation and usage instructions.

-------------------------------------------------------------------------------
Hierarchical Bayesian model for fitting a linear relation between two
latent variables (X, Y) with:

  - Measurement uncertainties in both variables (skew-normal likelihood
    for asymmetric errors, reducing to normal when errors are symmetric)
  - Intrinsic covariance between X and Y
  - Slope derived from intrinsic covariance:
        beta = rho * sigma_y / sigma_x
  - LKJ prior on correlation matrix (eta = 2 by default)
  - Non-centred parameterisation for improved sampler geometry

Model form in population space:

    Y = mu_y + beta * (X - mu_x)

This formulation is symmetric in (X, Y) and derives the slope from
the intrinsic covariance structure rather than assuming conditional
regression in one direction.

Designed for astrophysical scaling relations (e.g., M_BH-sigma, M_BH-M_bulge)
but applicable to any two-variable dataset where both variables carry
measurement uncertainty and the true underlying relation is linear.
===============================================================================
*/


//=============================================================================
data {
  int<lower=0> Npopulation;

// Arrays (written in R version 2.21 syntax; should auto-update to v.2.26+ if needed)
  real Y_mu[Npopulation];          
  real<lower=0> Y_sigma[Npopulation];
  real Y_alpha[Npopulation];

  real X_mu[Npopulation];
  real<lower=0> X_sigma[Npopulation];
  real X_alpha[Npopulation];
  
  int<lower=0,upper=1> inholdout_X[Npopulation];
  int<lower=0,upper=1> inholdout_Y[Npopulation];

// Scalars (these do not change between versions)
  real prior_mean_X_mu;
  real<lower=0> prior_mean_X_sd;

  real prior_mean_Y_mu;
  real<lower=0> prior_mean_Y_sd;

  real<lower=0> prior_sd_X_rate;
  real<lower=0> prior_sd_Y_rate;
}
//=============================================================================


//=============================================================================
parameters {
  // RAW masses (sampled from a standard normal distribution)
  vector[2] true_masses_raw[Npopulation];

  vector<lower=0>[2] sigma;        //  intrinsic dispersions
  cholesky_factor_corr[2] Lcorr;   //  population correlation

  real mean_X;
  real mean_Y;
}
//=============================================================================


//=============================================================================
transformed parameters {
  vector[2] xmean;
  matrix[2,2] L_cov;
  
  // ACTUAL true masses (calculated deterministically here)
  vector[2] true_masses[Npopulation];

  xmean[1] = mean_X;
  xmean[2] = mean_Y;

  // Cholesky factor of covariance matrix
  L_cov = diag_pre_multiply(sigma, Lcorr);
  
  // Non-Centered Parameterization: 
  // Scale and shift the raw standard normal values into the actual distribution
  for (i in 1:Npopulation) {
    true_masses[i] = xmean + L_cov * true_masses_raw[i];
  }
}
//=============================================================================


//=============================================================================
model {

// --------------------------------------------------------------------------
// LKJ prior on the intrinsic correlation matrix
// (Lcorr is the Cholesky factor of a 2x2 correlation matrix.)
//
// lkj_corr_cholesky(eta):
//   eta = 1  -> uniform over all valid correlation matrices
//   eta = 2  -> mild regularisation toward zero correlation (default)
//   eta > 2  -> stronger shrinkage toward zero correlation
//
// The choice of eta = 2 is a conservative default that provides mild
// stabilisation without strongly biasing the posterior. It is most
// beneficial for small samples (N < 50) where the data alone may not
// constrain rho well. For large samples (N > 200) the likelihood
// dominates and the choice of eta makes little practical difference.
//
// To use a uniform prior instead, replace eta = 2 with eta = 1 below.
// --------------------------------------------------------------------------
//  Lcorr ~ lkj_corr_cholesky(1);  // Uniform over correlations  
  Lcorr ~ lkj_corr_cholesky(2);   // eta > 1, e.g. 2, shrinks toward zero correlation (rho = 0) 

  sigma[1] ~ exponential(prior_sd_X_rate);
  sigma[2] ~ exponential(prior_sd_Y_rate);

  mean_X ~ normal(prior_mean_X_mu, prior_mean_X_sd);
  mean_Y ~ normal(prior_mean_Y_mu, prior_mean_Y_sd);

  // Sample the raw parameters from a flat, easy-to-navigate standard normal space
  for (i in 1:Npopulation) {
    true_masses_raw[i] ~ std_normal();
  }

  // Evaluate the likelihood of the data given the derived true_masses
  for (i in 1:Npopulation) {

if (inholdout_X[i]==0) {
  if (fabs(X_alpha[i]) < 1e-12)
    target += normal_lpdf(true_masses[i,1] |
                          X_mu[i],
                          X_sigma[i]);
  else
    target += skew_normal_lpdf(true_masses[i,1] |
                               X_mu[i],
                               X_sigma[i],
                               X_alpha[i]);
}
// X measurement uncertainty is modelled as a skew-normal, fitted to the
// reported asymmetric errors. When X errors are symmetric, X_alpha = 0 
// and this reduces exactly to a normal_lpdf, giving identical results 
// to the previous symmetric-only treatment.

if (inholdout_Y[i]==0) {
  if (fabs(Y_alpha[i]) < 1e-12)
    target += normal_lpdf(true_masses[i,2] |
                          Y_mu[i],
                          Y_sigma[i]);
  else
    target += skew_normal_lpdf(true_masses[i,2] |
                               Y_mu[i],
                               Y_sigma[i],
                               Y_alpha[i]);
}
// For the Skew-Normal (Black Hole mass), symmetry does not hold, but your pre-processing 
// steps in R (fitting the quantiles) correctly set this up to represent the probability 
// of the True mass given the Data, consistent with the "inverted likelihood" approach.
  }
}
//=============================================================================


//=============================================================================
generated quantities {

  real slope;
  real intercept;
  real intrinsic_scatter;
  real pivot_x;
  real pivot_y;
  real rho; 

// For a 2x2 Cholesky correlation factor, the [2,1] element is exactly rho
  rho = Lcorr[2,1]; 

// Slope derived from intrinsic covariance:
// beta = rho * sigma_y / sigma_x
  slope = rho * sigma[2] / sigma[1];

  pivot_x = mean_X;
  pivot_y = mean_Y;

  intercept = mean_Y - slope * mean_X;

// Conditional intrinsic scatter in Y at fixed X:
// sigma_{y|x} = sigma_y * sqrt(1 - rho^2)
  intrinsic_scatter = sigma[2] * sqrt(1 - square(rho));
}
//=============================================================================



