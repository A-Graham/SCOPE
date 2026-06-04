# Changelog

All notable changes to SCOPE will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
SCOPE uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html):
MAJOR.MINOR.PATCH, where MAJOR = breaking change, MINOR = new feature,
PATCH = bug fix or documentation update.

---

## [1.0.0] — 2026-03-08

### Initial release

**Core model**
- Hierarchical Bayesian linear regression via a bivariate normal population model
- Regression slope derived from the intrinsic covariance matrix: β = ρ σ_Y / σ_X
- Exact symmetry guaranteed: swapping X and Y yields 1/β to within numerical precision
- LKJ prior on the intrinsic correlation matrix for stable, weakly informative regularisation
- Non-centred parameterisation for efficient sampling in Stan

**Measurement error handling**
- Full propagation of measurement uncertainties in both X and Y
- Native support for asymmetric (skewed) error bars via a skew-normal likelihood
  approximation, without reducing them to a symmetric average
- Symmetric errors supported as a special case (skew-normal reduces to normal)

**Implementation**
- Stan model (`Scoped.stan`) compiled via `rstan` with automatic syntax upgrade
  for compatibility with both rstan < 2.26 and rstan >= 2.26
- Adaptive MCMC iteration count: 4000 iterations for N < 100, 2000 for N >= 100
- Sampler diagnostics reported: Rhat, n_eff, divergent transitions
- Posterior summaries: mean +/- 1-sigma (SD) and median with exact 2-sigma (~95.4%) credible interval

**Output**
- Publication-quality PDF figure (`SCOPE_fit.pdf`) with posterior mean and median
  regression lines, 1-sigma credible bands, error ellipses, and data points
- Full posterior samples saved to `Scout.dat` for downstream analysis
- Interactive cross-platform screen device fallback for Quicklook figure viewing

**Documentation**
- `README.md`: model overview, file descriptions, input format, usage instructions
- `INSTALL.md`: platform-specific installation guide for macOS, Linux, and Windows,
  with verification steps and troubleshooting
- `CITATION.cff`: machine-readable citation metadata
- `LICENSE`: MIT licence