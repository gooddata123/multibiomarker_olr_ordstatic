# multibiomarker_olr_ordstatic

Multivariate Ordinal Logistic Regression (OLR) for pro-arrhythmic risk prediction of 28 CiPA reference drugs, using biomarkers derived from *in silico* single-cell electrophysiology (EP) simulations under a static CiPAORd framework (i.e., standard conductance-scaling replaces the Markov-based dynamic IKr formulation).

---

## Overview

This repository implements an OLR-based risk classification pipeline that accepts 15 candidate biomarkers as inputs, evaluated either independently (single-input, *n* = 1) or in combination (*n* = 2–6). The upper limit of six biomarkers per combination is intentional: it preserves biological interpretability and statistical robustness, guarding against the overfitting and "black-box" opacity that arise in high-dimensional models.

---

## Repository Structure

## Data

Input biomarker values are pre-processed and split into two files under `/data/`:

| File | Description |
|------|-------------|
| `merged_biomarkers_zscore_train.csv` | Training set — 12 drugs, biomarker values z-score normalized |
| `merged_biomarkers_zscore_val.csv` | Validation set — normalized using training set mean and SD to prevent data leakage |

> **Note:** Standardization of the validation set using training statistics prevents information leakage and improves numerical stability during optimization.

---

## Usage

Set the `dimension` variable in the main script to control the number of biomarker inputs:

```r
dimension <- 1   # Single-input mode (each biomarker evaluated independently)
dimension <- 2   # All pairwise combinations (e.g., qNet + APD90)
dimension <- 3   # All 3-biomarker combinations  (e.g., qNet + APD90 + CaD50)
# ... up to
dimension <- 6   # Maximum supported combination size
```

Then run the script in R:

```r
source("main.R")
```

---

## Requirements

- R (≥ 4.x recommended)
- R packages: `MASS`, `pROC`, `ggplot2` , `readr` , `foreach`, `doParallel`

---
