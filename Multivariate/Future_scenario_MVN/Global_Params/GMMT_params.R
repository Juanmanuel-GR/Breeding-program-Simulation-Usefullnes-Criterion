# ---- GMMT Parameters ----
control <- TrainSel::SetControlDefault(
  size       = "large",
  complexity = "high_complexity",
  verbose    = FALSE
)
control$niterations <- 2000   # Increase if no convergence achieved

Proportion_selected     <- 50 / 100   # Selection intensity limit
chromosome_length_cM    <- 143        # Genetic length of a chromosome (cM)
n_families              <- nPYT / famMax   # Number of families
selfing_cycles          <- 0
mult_trait_coefficients <- NULL       # For multi-trait selection
phi                     <- 2          # Ploidy level
LD                      <- "Approx"
replication             <- FALSE
DHs                     <- TRUE

# ---- License Information ----
## CONTACT j.isidro@upm.es or javier.fgonzalez@upm.es for license information

# MateR license:
#username
#password

# TrainSel license:
#username_TSel
#password_TSel