# SCOPE Installation Guide

This guide walks through installing SCOPE on **macOS**, **Linux**, and **Windows**.
SCOPE requires R, the `rstan` package, and a C++ compiler. Follow the section
for your operating system, then complete the [Common Steps](#common-steps-all-platforms)
section and the [Verification](#verification) section at the end.

Estimated time: 15–30 minutes on a fresh system (most of this is download time).

---

## Table of Contents

- [macOS](#macos)
- [Linux](#linux)
- [Windows](#windows)
- [Common Steps (all platforms)](#common-steps-all-platforms)
- [Verification](#verification)
- [Troubleshooting](#troubleshooting)

---

## macOS

### Step 1 — Install R

Download and install the latest R release for macOS from:

```
https://cran.r-project.org/bin/macosx/
```

Choose the correct installer for your Mac:
- **Apple Silicon (M1/M2/M3/M4):** download the `arm64` `.pkg` file
- **Intel Mac:** download the `x86_64` `.pkg` file

If you are unsure which chip your Mac has, go to **Apple menu → About This Mac**.
Open the downloaded `.pkg` file and follow the installer prompts.

### Step 2 — Install Xcode Command Line Tools

`rstan` requires a C++ compiler. Open the **Terminal** application
(found in Applications → Utilities) and run:

```bash
xcode-select --install
```

A dialog box will appear asking you to install the tools. Click **Install**
and wait for the download to complete (this may take several minutes).
If the command returns `xcode-select: error: command line tools are already installed`,
you can skip this step.

### Step 3 — Configure R for rstan (Apple Silicon Macs only)

If you are on an Apple Silicon Mac (M1/M2/M3/M4), you need to tell R where
the compiler lives. In Terminal, run:

```bash
mkdir -p ~/.R
echo "CC=clang" >> ~/.R/Makevars
echo "CXX=clang++" >> ~/.R/Makevars
echo "CXX14=clang++" >> ~/.R/Makevars
```

Intel Mac users can skip this step.

### Step 4 — Proceed to [Common Steps](#common-steps-all-platforms)

---

## Linux

These instructions are written for Ubuntu/Debian. For Fedora/RHEL, replace
`apt` commands with the `dnf` equivalents.

### Step 1 — Install R

```bash
sudo apt update
sudo apt install r-base r-base-dev
```

This installs R and the development headers needed to compile R packages.
To verify the installation:

```bash
R --version
```

You should see R version 4.0 or later. If your distribution's package
manager provides an older version, follow the instructions at
`https://cran.r-project.org/bin/linux/ubuntu/` to add the CRAN repository
and install a current version.

### Step 2 — Install C++ build tools and dependencies

```bash
sudo apt install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev
```

These libraries are required to compile `rstan` and its dependencies.
`libcairo2-dev` is included to ensure the Cairo graphics engine is available,
which SCOPE uses to correctly render the solar mass symbol (⊙) in the output
figure. Ubuntu includes Cairo by default, but this command guarantees it.

### Step 3 — Proceed to [Common Steps](#common-steps-all-platforms)

---

## Windows

### Step 1 — Install R

Download and install the latest R release for Windows from:

```
https://cran.r-project.org/bin/windows/base/
```

Run the downloaded `.exe` installer and follow the prompts, accepting
all default options.

### Step 2 — Install RTools

`rstan` requires a C++ compiler on Windows, provided by RTools.
Download RTools from:

```
https://cran.r-project.org/bin/windows/Rtools/
```

Select the version of RTools that matches your installed R version
(the page lists which version to use). Run the downloaded installer
and follow the prompts, accepting all default options.

> **Important:** After installing RTools, open R and run the following
> to confirm that R can find the compiler:
>
> ```r
> Sys.which("make")
> ```
>
> If this returns an empty string, restart R (or your computer) and try again.

### Step 3 — Proceed to [Common Steps](#common-steps-all-platforms)

---

## Common Steps (all platforms)

These steps are the same on macOS, Linux, and Windows.
Open R (or RStudio) and run the following commands.

### Step 1 — Install rstan

```r
install.packages("rstan", repos = "https://cloud.r-project.org")
```

This will download and compile `rstan` along with its dependencies.
It may take 5–15 minutes. You will see a lot of compiler output — this
is normal. The install is complete when you see the R prompt (`>`) again.

If prompted to install from source or to update existing packages,
answer `n` (no) to update existing packages — this avoids inadvertently
breaking other installed packages.

### Step 2 — Install remaining R packages

```r
install.packages(c("sn", "MASS", "ellipse"), repos = "https://cloud.r-project.org")
```

These are smaller packages and should install in under a minute each.

### Step 3 — Download SCOPE

Download the SCOPE repository from GitHub. You can do this in two ways:

**Option A — using git (recommended):**
```bash
git clone https://github.com/A-Graham/SCOPE.git
cd SCOPE
```

**Option B — as a ZIP file:**  
Go to `https://github.com/A-Graham/SCOPE`, click the green **Code** button,
select **Download ZIP**, then unzip the downloaded file.

The folder should contain:
```
SCOPE.R
Scoped.stan
Scin.csv
```

### Step 4 — Run SCOPE

Place your data file (named `Scin.csv`) in the same folder as `SCOPE.R`
and `Scoped.stan`. Then run SCOPE using either of the following methods:

**From the terminal / Command Prompt:**
```bash
Rscript SCOPE.R
```

**From within R or RStudio:**
```r
setwd("/path/to/your/SCOPE/folder")
source("SCOPE.R")
```

On the first run in a given R session, `rstan` will compile the Stan model
into C++ code in the background. This takes 1–3 minutes. Within the same
R session, subsequent runs will use the cached model in memory and start
sampling almost immediately.

---

## Verification

To confirm that everything is installed correctly before running SCOPE
on your own data, run the included example using the provided `Scin.csv`
file. From within R:

```r
setwd("/path/to/your/SCOPE/folder")
source("SCOPE.R")
```

A successful run produces output similar to the following:

```
================================================================================
 SCOPE: Symmetric COvariance Population Estimator
 Version 1.0
 Author : Alister W. Graham (2026)
 GitHub : https://github.com/A-Graham/SCOPE
================================================================================

Initialising environment and loading packages...

 NOTE: Stan must compile the model to C++ on first use in each session.
       This typically takes 1-3 minutes. Subsequent runs in the same
       session will start sampling almost immediately.
--------------------------------------------------------------------------------

Input data   : N = 38 objects
MCMC samples : 4 chains x 1000 post-warmup draws = 4000 total

Sampler Diagnostics
  Max Rhat                  : 1.00X  (Good: all chains converged)
  Parameters with Rhat>1.01 : 0  (Good: no parameters flagged)
  Min n_eff                 : XXXX  (Good: adequate effective sample size)
  Divergent transitions     : 0  (Good: none detected)
------------------------------------------------------------------------------------

====================================================================================
 Model form:
   Y = mu_y + beta * (X - mu_x)
   beta = rho * sigma_y / sigma_x

 Posterior summaries:
   Mean +/- 1-sigma (posterior standard deviation)
   Median with 2-sigma (~95.4%) credible interval
------------------------------------------------------------------------------------
Regression relation:
Pivot Y (mu_y)            : Mean:  9.087 +/- 0.102 | Median:  9.088 [ 8.889,  9.289]
Pivot X (mu_x)            : Mean:  2.428 +/- 0.013 | Median:  2.428 [ 2.402,  2.454]
Slope (beta)              : Mean:  7.848 +/- 1.386 | Median:  7.788 [ 5.294, 10.675]
```

followed by further posterior summaries and scatter diagnostics, and ending with:

```
Plot saved to SCOPE_fit.pdf
```

If you see this output the installation is complete and working correctly.

---

## Troubleshooting

### rstan fails to install or compile

The most common cause is a missing or misconfigured C++ compiler.

- **macOS:** confirm Xcode Command Line Tools are installed by running
  `xcode-select -p` in Terminal — it should print a path. If it prints nothing,
  re-run `xcode-select --install`.

- **Windows:** confirm RTools is on the PATH by running `Sys.which("make")` in R.
  If it returns an empty string, reinstall RTools and restart R.

- **Linux:** confirm `build-essential` is installed by running
  `g++ --version` in a terminal. If the command is not found, run
  `sudo apt install build-essential`.

If the issue persists, consult the rstan installation wiki at:
```
https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
```

### "object not found" or "package not found" errors when running SCOPE.R

Run the following in R to check all required packages are installed:

```r
required <- c("rstan", "sn", "MASS", "ellipse")
missing  <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
    cat("All required packages are installed.\n")
} else {
    cat("Missing packages:", paste(missing, collapse=", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org")
}
```

### Stan model fails with "syntax error" or "parser error"

Because SCOPE automatically upgrades the Stan syntax to match your installed
version of `rstan`, parser errors are rare. If you encounter one (often caused
by modifying `Scoped.stan` manually), simply restart your R session to clear
the temporary compiler cache, and run the script again.

### SCOPE runs but produces divergence warnings

Divergences indicate that the sampler is having difficulty exploring
the posterior. Try increasing `adapt_delta` in `SCOPE.R`:

```r
adapt_delta = 0.99   # increase from default 0.95
```

If warnings about maximum treedepth appear, also increase `max_treedepth`:

```r
max_treedepth = 15   # increase from default 12
```

### The PDF figure is blank or does not open

On Linux, confirm that a PDF viewer is installed. On macOS, `SCOPE_fit.pdf`
can be opened with Preview. On Windows, use any PDF reader such as
Adobe Acrobat or SumatraPDF.

If the figure file is not created at all, check that R has write permission
to the folder containing `SCOPE.R`.

### Getting further help

If you encounter a problem not covered here, please open an issue at:

```
https://github.com/A-Graham/SCOPE/issues
```

Include the output of `sessionInfo()` (run from within R) and the full
error message — this provides the information needed to diagnose most problems quickly.
