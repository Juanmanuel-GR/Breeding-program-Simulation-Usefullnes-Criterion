

# Optimal Mate Allocation in a Multi-Trait Wheat Breeding Program



Cloning, selection, and cross-optimization stochastic simulation framework for wheat breeding using `AlphaSimR`.



[!\[R-Version](https://img.shields.io/badge/R-%2B4.0-blue.svg)](https://www.r-project.org/)

[!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



---



##  Project Overview



This repository contains the R code and simulation framework developed for my Master's Thesis (TFM). The project implements a comprehensive stochastic simulation to evaluate long-term genetic gain and genetic diversity in wheat (\*Triticum aestivum\*) over multiple breeding cycles.



The core objective is to evaluate the impact of the \*\*Usefulness Criterion (UC)\*\* within a crossing optimization framework. Unlike traditional methods that only consider the family mean, the UC integrates both the family mean and the within-family genetic standard deviation. This allows the program to identify specific crosses that are highly likely to produce elite, transgressive segregants, thereby maximizing long-term genetic gain.



Additionally, this project assesses the performance of the \*\*Desired Gain Index\*\* in mitigating the well-documented negative correlation between Grain Yield and Protein Content.



---



##  Scenarios Evaluated



The framework evaluates \*\*4 distinct scenarios\*\*, each simulated under both \*\*Single-Trait (Univariate)\*\* and \*\*Multi-Trait (Multivariate)\*\* selection schemes:



1. **`PS-RM` (Phenotypic Selection + Random Mating):** Selection is based strictly on observed phenotypic values, and crosses are performed at random.

2. **`GS-RM` (Genomic Selection + Random Mating):** Selection is guided by Genomic Estimated Breeding Values (GEBVs), with random parental mating.

3. **`MtR-2f` (Genomic Mating - 2 Factors):** Selection on GEBVs, with crossing optimization balancing the **Family Mean** against **Inbreeding Control** (penalizing crosses between closely related individuals).

4. **`MtR-3f` (Genomic Mating - 3 Factors):** Selection on GEBVs, with crossing optimization considering three dimensions: **Family Mean**, **Inbreeding Control**, and **Within-Family Variance** (The Usefulness Criterion framework).



---

## How to use this script

The repository is structured into two main tracks based on dimensionality: \*\*Univariate\*\* (single-trait) and \*\*Multivariate\*\* (multi-trait) selection. Both tracks follow the exact same operational workflow.







---

\## Software \& Dependencies



The simulation framework is built entirely in **R (version $\\ge$ 4.0)** and relies on a specialized ecosystem of quantitative genetics and breeding optimization packages:



* **`AlphaSimR`**: Used as the core engine to simulate the wheat genome architecture (14 chromosomes, SNPs, and pleiotropic QTLs), individual founder lines, and the operational breeding pipeline (replications, stages, and phenotypic noise).

* **`sommer`**: Used in the multi-trait phenotypic scenarios to fit the linear mixed models and estimate the genetic and phenotypic variance-covariance matrices ($\\Gamma$ and $P$).

* **`AGHmatrix`**: Utilized via the `Amatrix` function to compute the pedigree-based relationship matrix ($A$), providing the explicit covariance structure required for variance partitioning under phenotypic selection.

* **`StageWise`**: Implements the Pešek and Baker (1969) desired gain index framework. It calculates the index coefficients ($b$) to scale and target specific genetic gains ($d = (1, 0)'$).

* **`MateR`**: Executes the genomic mating optimization via the `GenomicMatingMT` function. It drives the greedy search algorithm (2,000 iterations) to find near-optimal mating plans maximizing the Usefulness Criterion.

&#x20;

