# SCOPE: Symmetric COvariance Population Estimator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20536917.svg)](https://doi.org/10.5281/zenodo.20536917)
[![arXiv](https://img.shields.io/badge/arXiv-2606.05808-b31b1b.svg)](https://arxiv.org/abs/2606.05808)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SCOPE** is a hierarchical Bayesian regression framework implemented in R and Stan for fitting a linear scaling relation between two observed quantities while fully accounting for:

* Measurement uncertainties in **both** variables
* Asymmetric (skewed) measurement errors
* Intrinsic population covariance
* Intrinsic conditional scatter

The regression slope is derived from the **intrinsic population covariance**, avoiding the bias that arises in ordinary least squares when measurement errors in the independent variable are ignored.

Companion software for:
> **Alister W. Graham (2026)**  
> *"Galaxy morphology dependent (black hole mass)-(velocity dispersion) relations: implications for gravitational wave forecasts and cosmological simulations"*  
> **arXiv:2606.05808 [astro-ph.GA]** -- [https://arxiv.org/abs/2606.05808](https://arxiv.org/abs/2606.05808)

---

## Model Overview

SCOPE models the intrinsic population as a bivariate normal distribution:

$$(X, Y) \sim \mathcal{N}\left(\begin{bmatrix} \mu_X \\ \mu_Y \end{bmatrix}, \Sigma\right)$$

with covariance matrix:

$$\Sigma = \begin{bmatrix} \sigma_X^2 & \rho \sigma_X \sigma_Y \\ \rho \sigma_X \sigma_Y & \sigma_Y^2 \end{bmatrix}$$

The reported regression slope is derived from the intrinsic covariance:

$$\beta = \rho \frac{\sigma_Y}{\sigma_X}$$

The linear relation is expressed in pivoted form:

$$Y = \mu_Y + \beta (X - \mu_X)$$

The conditional intrinsic scatter is:

$$\sigma_{Y|X} = \sigma_Y \sqrt{1 - \rho^2}$$

---

## Symmetry of the Model

A known failure of one-sided conditional estimators (e.g. ordinary least squares or standard conditional Bayesian fits) is **directional bias**: fitting $Y$ on $X$ yields a slope $\beta_1$, whereas fitting $X$ on $Y$ yields a slope $\beta_2$ whose reciprocal does not match ($1/\beta_2 \neq \beta_1$). This discrepancy ($\Delta\beta$) is often as large as the $1\sigma$ uncertainty on the slope itself.

**SCOPE resolves this by fitting the joint 2D intrinsic population distribution.**

Because the regression parameters are derived directly from the underlying intrinsic covariance matrix $\Sigma$:

$$\beta_{Y|X} = \rho \frac{\sigma_Y}{\sigma_X} \qquad \text{and} \qquad \beta_{X|Y} = \rho \frac{\sigma_X}{\sigma_Y}$$

the model mathematically guarantees exact directional invariance:

$$\beta_{Y|X} \equiv \frac{1}{\beta_{X|Y}}$$

across **every single MCMC posterior draw**. 

Whether a user assigns a physical quantity to the $X$-axis or the $Y$-axis, the fitted intrinsic relation describes the exact same physical line in parameter space without regression dilution or axis-dependent bias.

---

## Files in This Repository

| File | Description |
|---|---|
| `SCOPE.R` | Main driver script:<br>— reads data from `Scin.csv` (or a custom-named file; see Running SCOPE)<br>— prepares skew-normal likelihood parameters<br>— compiles and runs Stan<br>— produces diagnostic output<br>— saves posterior summaries |
| `Scoped.stan` | Stan model implementing:<br>— non-centred parameterisation<br>— LKJ prior on the correlation matrix<br>— intrinsic covariance modelling<br>— skew-normal likelihood for asymmetric errors<br>— derived slope and intrinsic scatter in `generated quantities` |
| `Scin.csv` | Example input data file (default input; rename your data file to `Scin.csv` or pass a custom filename at runtime) |

---

## Input Data Format (`Scin.csv`)

A CSV file named Scin.csv with columns:
```
label     : object name
Mbh       : log10 black hole mass (Y variable)
dMbh_pos  : +1-sigma uncertainty in Mbh
dMbh_neg  : -1-sigma uncertainty in Mbh
Sig       : log10 velocity dispersion, or other X variable
dSig_pos  : +1-sigma uncertainty in Sig
dSig_neg  : -1-sigma uncertainty in Sig
```

Symmetric fallbacks: supply dMbh instead of dMbh_pos/dMbh_neg,
and/or dSig instead of dSig_pos/dSig_neg.

The format in `Scin.csv` can therefore be any of the following:
```
label,Mbh,dMbh_pos,dMbh_neg,Sig,dSig_pos,dSig_neg
label,Mbh,dMbh_pos,dMbh_neg,Sig,dSig
label,Mbh,dMbh,Sig,dSig_pos,dSig_neg
label,Mbh,dMbh,Sig,dSig
```

Asymmetric errors are converted in `SCOPE.R` into skew-normal parameters so that the likelihood correctly represents:

p(true value | reported measurement)

When errors are symmetric, the skew-normal reduces exactly to a normal distribution.

---

## System Requirements

**R version:** 4.0 or later — download from https://www.r-project.org

**R packages** (install once from within R):
```r
install.packages(c("rstan", "sn", "MASS", "ellipse"))
```

**C++ compiler** (required by rstan to compile the Stan model):
- macOS: Xcode Command Line Tools — install by running `xcode-select --install` in Terminal
- Windows: RTools — download from https://cran.r-project.org/bin/windows/Rtools
- Linux: GCC — typically pre-installed; if not, run `sudo apt install build-essential` (Ubuntu/Debian)

For detailed rstan installation instructions see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

**Tested on:**
- macOS 10.15 (Catalina) and later
- Linux (Ubuntu 22.04)
- Windows 10

---

## Sampling Configuration

SCOPE runs 4 independent MCMC chains. The number of iterations scales
automatically with sample size:

| Sample size | Iterations | Warmup |
|---|---|---|
| N < 100 | 4000 | 2000 |
| N ≥ 100 | 2000 | 1000 |

Other settings:
- `adapt_delta = 0.95` — controls step size; increase toward 0.99 if divergences are reported
- `max_treedepth = 12` — controls sampler depth; increase to 15 or 20 if warnings appear
- `seed = 1234` — set for reproducibility

**Compatibility:** the Stan model syntax is automatically upgraded at
runtime for `rstan` ≥ 2.26, so the same `Scoped.stan` file works with
both old and modern versions of `rstan` without any manual editing.

---

## Running SCOPE

Place `SCOPE.R`, `Scoped.stan`, and `Scin.csv` in the same folder, then
run using either of the following methods — both work identically on
macOS, Linux, and Windows:

**From the terminal** (Command Prompt on Windows, Terminal on macOS/Linux):
```
Rscript SCOPE.R                       # uses Scin.csv by default
Rscript SCOPE.R my_data.csv           # uses a custom-named input file
```

**From within R or RStudio:**
```r
setwd("/path/to/your/SCOPE/folder")   # set to wherever your files are
source("SCOPE.R")                      # uses Scin.csv by default
```

To use a custom-named input file from within R or RStudio, set
`input_file` before sourcing:
```r
setwd("/path/to/your/SCOPE/folder")   # set to wherever your files are
input_file <- "my_data.csv"
source("SCOPE.R")
```

---

## Output Files

Results are printed to the console, and an interactive figure is displayed.
Each run of SCOPE creates the following files in the same folder:

| File | Description | Safe to delete? |
|---|---|---|
| `SCOPE_fit.pdf` | Publication-quality figure of the fitted relation | Yes, recreated each run |
| `Scout.dat` | Saved posterior samples (see below) | Yes, recreated each run |
| `Scoped.rds` | Cached compiled Stan model (speeds up repeat runs) | Yes, recreated automatically |

`SCOPE_fit.pdf` is the primary output for most users. `Scout.dat`
stores the full posterior samples in R's native binary format and can
be reloaded in a later session using `load("Scout.dat")`, allowing
further analysis without re-running the model. Both `SCOPE_fit.pdf`
and `Scout.dat` are regenerated on each run and may be safely deleted
at any time. `Scoped.rds` caches the compiled Stan model so subsequent
runs start immediately without recompilation — it is safe to delete,
but keeping it avoids the 1–3 minute compile step. `Scoped.rds` is 
tied to your installed version of `rstan` — if you upgrade `rstan` 
and encounter errors, deleting `Scoped.rds` and re-running will 
resolve them.

---

## Intended Applications

SCOPE was developed for astrophysical scaling relations, including:

- Black hole mass — velocity dispersion (M_BH — sigma)
- Black hole mass — bulge mass (M_BH — M_bulge)
- Galaxy structural scaling relations

but is applicable to any two-variable dataset where both variables carry
measurement uncertainty and the true underlying relation is linear.

---

## Author

Alister W. Graham

---

## Citation

If you use SCOPE in published work, please cite:

```bibtex
@MISC{Graham2026_SCOPE_software,
       author = {{Graham}, Alister W.},
        title = "{SCOPE: Symmetric COvariance Population Estimator}",
         year = 2026,
        month = jun,
          doi = {10.5281/zenodo.20536917},
    publisher = {Zenodo},
          url = {https://doi.org/10.5281/zenodo.20536917}
}

@ARTICLE{Graham2026_SCOPE_paper,
       author = {{Graham}, Alister W.},
        title = "{Galaxy morphology dependent (black hole mass)-(velocity dispersion) relations: implications for gravitational wave forecasts and cosmological simulations}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Astrophysics of Galaxies},
         year = 2026,
        month = jun,
          eid = {arXiv:2606.05808},
        pages = {arXiv:2606.05808},
archivePrefix = {arXiv},
       eprint = {2606.05808},
 primaryClass = {astro-ph.GA},
          url = {https://arxiv.org/abs/2606.05808}
}
```

The latest citation information is provided in CITATION.cff 
in the GitHub repository.

---

## License

MIT License

Copyright (c) 2026 Alister W. Graham

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

---

## Bug Reports and Feedback

Please use the [GitHub Issues](https://github.com/A-Graham/SCOPE/issues)
page to report bugs or suggest improvements.