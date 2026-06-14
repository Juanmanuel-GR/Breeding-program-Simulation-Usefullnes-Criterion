# Global parameters i need through the simulation

#number of simulations replications and breeding cycles 
nReps <- 40
nBurnin <- 20
nfuture <- 20
nCycles <- nBurnin + nfuture
startTP <- 16


# Genone Simulation 
nChr <- 14
nQtl <- 500
nSnp <- 200

#Initial parents mean and variance
initmeanG <- c(1,2)
initVarG <- c(2,4)
initVarEnv <- c(1e-5, 1e-5)
initVarGE <- c(2, 2)
VarE <- c(4,4) # Yield trial error variance, bushels per acre
# Relates to error variance for an entry mean

# Breeding program details
# ---- Breeding program details ----
nParents = 50  # Number of parents to start a breeding cycle
nCrosses = 100 # Number of crosses per year
nDH      = 100 # DH lines produced per cross
famMax   = 2  # The maximum number of DH lines per cross to enter PYT
nPYT     = nDH*famMax # Entries per preliminary yield trial
nAYT     = 50  # Entries per advanced yield trial
nEYT     = 10  # Entries per elite yield trial

#Effective replication of yield trials 
repHDRW <- 4/9
repPYT <- 1
repAYT <- 4
repEYT <- 8


#params for population evolution 

#nEvolution <- 10 #number of years of the population to be develop





