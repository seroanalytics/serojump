# serojump

![Schematic of the methods behind the `serojump` package](./man/figures/schematic_A.png)

## Overview

The `serojump` package provides tools for fitting serological models to antibody kinetics data using reversible-jump Markov Chain Monte Carlo (RJ-MCMC). It enables researchers to model the dynamics of antibody levels in response to infections, incorporating both observational and antibody kinetics models. The package supports the inclusion of priors and various exposure scenarios, making it highly flexible for serological data analysis.

This package is ideal for researchers looking to:

- Model antibody kinetics over time
- Estimate infection rates and antibody waning using serological data
- Perform MCMC inference on serological models

The best place to get started with serojump is our preprint (coming soon). In this repo, ceveral model templates are provided to simplify usage, while also allowing users to customize models for specific research questions.

## Installation

To install the `serojump` package, follow these steps:

### Step 1: Install R

Make sure you have R installed on your system. You can download R from [https://cran.r-project.org/](https://cran.r-project.org/).

### Step 2: Install `serojump` from GitHub

You can install the development version of `serojump` from GitHub using the `devtools` package. If you don't already have `devtools` installed, you can install it with:


```

install.packages("devtools")
devtools::install_github("seroanalytics/serojump")

```

### Step 3: Load `serojump` from GitHub
After installation, you can load the serojump package into your R session with:


```

library(serojump)

```

## Usage

For detailed usage instructions, please refer to the package vignettes and examples. These vignettes explain:

- The format of the data require for serojump to work, click [here](https://seroanalytics.org/serojump/articles/data_format.html).
- How to define observational models and antibody kinetics models and run the serojump sampler, click [here](https://seroanalytics.org/serojump/articles/model_define.html).
- Specify priors for your models and integrate empirical data [coming soon]

In addition we have several worked examples of `serojump` on:
- Simulated data click [here](https://seroanalytics.org/serojump/articles/sim_recovery.html).

---

## Contributing

We welcome contributions and suggestions! If you'd like to contribute to the `serojump` package or report issues, please feel free to:

- Submit a pull request on GitHub.
- Open an issue on the repository.

If you have any questions or feedback, or would like more informative vignettes, you can contact the package maintainer at:

**David Hodgson**  
Email: [david.hodgson@lshtm.ac.uk](mailto:david.hodgson@lshtm.ac.uk)

---

## Project Status

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/seroanalytics/serojump/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/seroanalytics/serojump/actions/workflows/R-CMD-check.yaml)

This package is actively maintained and in a stable, usable state. New features and improvements are continually being developed.
