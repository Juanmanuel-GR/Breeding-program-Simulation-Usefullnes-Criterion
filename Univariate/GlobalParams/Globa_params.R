# Global parameters i need through the simulation

#number of simulations replications and breeding cycles 
nReps <- 40
nBurnin <- 20
nfuture <- 20
nCycles <- nBurnin+ nfuture

startTP <- 16
start_Pedigree <- 15

# Genone Simulation 
nChr <- 14
nQtl <- 500
nSnp <- 200

#Initial parents mean and variance
initmeanG <- 1
initVarG <- 2
initVarEnv <- 1e-5
initVarGE <- 2
VarE <- 4 # Yield trial error variance, bushels per acre
# Relates to error variance for an entry mean

# Breeding program details
# ---- Breeding program details ----
nParents = 50  # Number of parents to start a breeding cycle
nCrosses = 100 # Number of crosses per year
nDH      = 100 # DH lines produced per cross
famMax   = 2#5  # The maximum number of DH lines per cross to enter PYT
nPYT     = nDH*famMax#500 # Entries per preliminary yield trial
nAYT     = 50  # Entries per advanced yield trial
nEYT     = 10  # Entries per elite yield trial

#Effective replication of yield trials 
repHDRW <- 4/9
repPYT <- 1
repAYT <- 4
repEYT <- 8
