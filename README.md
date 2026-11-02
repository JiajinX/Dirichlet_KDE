## Overview

Companion package for the article *Bandwidth rules for Dirichlet kernel density estimation on the simplex*, including:

-   simulation
-   real data application
-   tables / figures

## Directory structure

```         
project/
├── data/       # Data for applications
├── R/          # General-purpose functions
├── scripts/    # Analysis pipelines (call functions from R/)
├── results/    # Data sources for tables and images
├── figures/    # Figures for the paper
```

## Requirements

-   R version: 4.4.1 (or later)
-   R packages:\
    `compositions`, `dplyr`, `ggplot2`, `gtools`, `ks`, `lsei`, `misc3d`, `mvtnorm`, `paletteer`, `rgl`, `rnaturalearth`, `sf`, `stats`, `viridis`, `WDI`

## Data

The folder `/data` stores the datasets for applications.

-   *Wholesle Dataset*: the dataset is used in the script `wholesale.R`.
-   *Employment Dataset*: the dataset is processed and used in the script `employ.R`. The raw data comes from the `WDI` package.

## Replicability

The folder `/scripts` contains the scripts to replicate the simulation and real data application of the paper.

-   The script `fig1.R` generates (projected) gradient fields the Dirichlet distribution (Figures `5.1-projected_gradient.pdf` and `5.1-gradient.png`).
-   The script `MISE_MC.R`, `fig2.R` and `fig3.R` are related to two Monte Carlo evaluation methods for MISE:
    - `MISE_MC.R` calculates the MISE values (Results in folder `MISE_MC`).
    - `fig2.R` shows the MISE curves (Figure `6.1-MISE_curve.pdf`).
    - `MISE_time_SE.R` records the computation times and the standard errors (Results in folder `MISE_time_SE`).  
-   The script `LSCV_RSLS.R` and `fig3.R` are related to LSCV and RSLS selectors:
    -  `LSCV_RSLS.R` calculates the LSCV and RSLS values (Results in folder `LSCV_RSLS`)
    -  `fig3.R` shows the mean curves and the corresponding distributions of the selected bandwidths (Figure `6.2-LSCV_RSLS.pdf`).

For simulation:

-   The script `fig4.R` generates contour plots of Gaussian and Dirichlet mixture simulations (Figure `7-simumix.pdf`).
-   The script `simu.R` reproduces the simulation experiments (Results in folder `Simulation`).

For real data application:

-   The script `wholesale.R` reproduces the selector comparison experiment (Results in folder `Wholesale`) and generates the density contour from the EMI selector (Figure `8.1-wholesale.pdf`).
-   The script `employ.R` reproduces the clustering experiment and generates density-ascent trajectories (Figure `8.2-clustering.pdf`), basin of attraction for clusters (Figure `8.2-clustering.pdf`), and clusters assignment on a world map (Figure `8.2-worldmap.pdf`).

## Figures and Results

The folder `figures/` and `/results` contain all outputs generated from simulations and real data applications.
