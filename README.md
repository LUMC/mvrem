# JSM: Joint Scaled Models for INLA

An R package for constructing and fitting Joint Scaled Models (JSM) using the Integrated Nested Laplace Approximation (INLA) approach.

## Overview

This package provides tools for specifying complex hierarchical models with multiple longitudinal and outcome processes. It facilitates handling of both fixed and random effects in a unified framework, with support for multiple distributions and flexible random effects structures.

## Installation

### Prerequisites

This package requires **R-INLA**. Please install it following the official instructions at: https://www.r-inla.org/download

### Installing jsm

You can install the jsm package from the GitLab repository:

```r
# Install from GitHub
remotes::install_github("LUMC/mvrem", ref = "devel")

# Or during development
devtools::install()
devtools::load_all()
```

## Quick Start

```r
library(mvrem)

# Load example data
d <- read.table("data/simulated_t_x1x2_y.tsv", header = TRUE, sep = "\t")

# Create a JSM model
model <- jsm(
  formLong = list(x1 ~ 1 + t + (1 + t | id)),
  dataLong = d,
  formOut = list(y ~ 1 + t + (1 | id)),
  dataOut = d,
  fit = FALSE  # Set to TRUE to auto-fit
)

# Inspect the model formula
cat(model$getFormulaText())

# Get the model data
data <- model$getData()

# Fit using INLA
result <- model$callInla()
result
```

## Main Functions

- **`jsm()`** - Create a Joint Scaled Model object
- **`JSMMdl`** - The main R6 class for Joint Scaled Models (JSM) models

## JSMMdl Class Methods

- `getFormula()` - Get the INLA formula
- `getFormulaText()` - Get formula as text with comments
- `getData()` - Get the prepared data for INLA
- `getFamilies()` - Get outcome distributions
- `getControlFamilies()` - Get family control parameters
- `getLevels()` - Get factor levels for random effects
- `callInla()` - Fit the model using INLA

## Features

- Multiple outcome variables
- Flexible random effects structures (intercepts, slopes, correlated effects)
- Support for multiple distributions (gaussian, binomial, poisson, etc.)
- Automatic handling of fixed and random effects
- Copy mechanism for linking longitudinal and outcome processes
- Comprehensive model diagnostics

## Documentation

- Example code in `doc/demo.Rmd`

## References

- R-INLA: https://www.r-inla.org/
- Joint modeling with INLA: https://becarioprecario.bitbucket.io/inla-gitbook/
