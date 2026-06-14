# Future_GS_Yield_RM.R

# ---- Clean Environment and Load Packages ----
rm(list = ls())
library(AlphaSimR)
library(lme4breeding)
library(tidyr)
library(dplyr)

#It was run in a SLURN HPC
REP      <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "GS_UNV_RM"



# ---- Helper Functions ---- 
source(file = "GlobalParams/sim_functions_mult.R")
source(file = "GlobalParams/Globa_params.R")



#-----Input parameters------
burnin_dir <- "../../Multivariate/Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
output_csv_gs <- file.path(results_dir, paste0("GS_UNV_RM_Results_", REP, ".csv"))

#---- Create Results Directory ----
if (dir.exists(results_dir)) {
  dir.create(results_dir)
}

# ---- Initializa CSV File -----
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#Initialize CSV file
if(file.exists(output_csv_gs)) file.remove(output_csv_gs)

#create the header for the csv file
header <- c("rep",
            "year",
            "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein", "heritability")

write.table(
  t(header),
  file = output_csv_gs,
  sep = ",",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

#---Future PS Phase 
cat("Processing the future scenario of REP", REP, "SCENARIO", SCENARIO, "\n")

burnin_file <- file.path(burnin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)



# ---- Extend 'output' for PS phase ----
output <- bind_rows(
  output,
  data.frame(
    rep         = REP,
    year        = (nBurnin + 1):(nBurnin + nfuture),
    meanG_yield       = numeric(nfuture),
    meanG_protein = numeric(nfuture),# si incluyera ambos datos aunq fuera genetic selection solo para este 
    varG_yield        = numeric(nfuture),
    varG_protein  = numeric(nfuture),
    genicVarG_yield   = numeric(nfuture),
    genicVarG_protein = numeric(nfuture),
    heritability = numeric(nfuture)
  ))



#start the year loop
# Future GS Simulation Loop
for (year in (nBurnin+1):(nBurnin+nfuture)) {
  cat("Future year:", year, "\n")
  
  #genomic Selection
  M_all <- pullSnpGeno(TrainPop)
  ids_all <- TrainPop@id
  id_unique <- as.integer(unique(ids_all))
  unique_positions <- match(id_unique, ids_all)
  M_train <- M_all[unique_positions, ]
  rm(M_all)
  #Filter 0 variance Markers
  M_train <- filter_zero_variance(M_train)
  
  #Compute the VanRaden additive GRM
  GRM <- getVR_GRM(M_train)
  GRM <- GRM +diag(1e-6, nrow(GRM))
  #Create the dataframe for the phenotype
  phenotype <- pheno(TrainPop)
  fixed <- TrainPop@fixEff # year fixed effect
  
  y_train <- data.frame(GID=ids_all,
                        Yield = as.numeric(phenotype[,1]))
  #add fixed effect as factor
  y_train$YearEff <- as.factor(fixed)
  
  #genomic selection using lmbe4breeding
  #building models
  starting_time <- Sys.time()
  cat("Runnin lme4breeding models...\n")
  model <- lmeb(Yield~ YearEff +
                  (1|GID),
                relmat = list(GID=GRM),
                verbose = 0L, trace = 0L,
                data = y_train)
  
  end_time <- Sys.time()
  time_took <- end_time - starting_time
  cat("Model Converged, took" , time_took, "\n")  
  ################################################################################
  # Extract estimates from model 
  ################################################################################
  
  #GEBVs - the BLUPS
  GEBVs <- ranef(model)$GID
  
  #recover the same order
  GEBVs <- GEBVs[rownames(M_train),]
  hist(as.numeric(GEBVs))
  #Compute heritability
  vc <- VarCorr(model) ; print(vc, comp=c("Variance"))
  ve <- ve <- attr(VarCorr(model), "sc")^2
  #heritability
  h2_mod <- vc$GID[1] / ( vc$GID[1] + ve[1])
  
  h2_1 <- var(PYT@gv)[1,1] / var(PYT@pheno)[1,1]
  
  print(h2_mod)
  print(h2_1)
  
  output$heritability[year] <- h2_mod
 
  
  ################################################################################
  # Extract additive marker effects
  ################################################################################
  
  marker_eff <- extract_meff(M_train, GEBVs)
  trait_name   <- SP$traitNames[[1]]
  colnames(marker_eff) <- trait_name
  
  #Predict EBVs and select new Parents, we predict on the EBV object of the population
  
  PYT@ebv <- setGEAVs(PYT, marker_eff)
  
  modAcc <- sqrt(h2_mod)
  
  # Store general model accuracies
  #output$accSel[year] <- cor(PYT@gv, PYT@ebv)
  #output$modAcc[year] <- modAcc
  
  PYT@misc$Family <- factor(paste(PYT@father, PYT@mother, sep = "/"))
  
  # Create dataframe for within-family accuracy
  phenoDF_TS <- data.frame(
    GID      = PYT@id,
    Pheno    = c(PYT@pheno),
    EBVs     = c(PYT@ebv),
    Family   = PYT@misc$Family,
    True_BVs = c(PYT@gv)
  )
  
  # Center values by family mean for PYT population
  for (family in unique(phenoDF_TS$Family)) {
    rows                 <- which(phenoDF_TS$Family == family)
    meanPheno            <- mean(phenoDF_TS$Pheno[rows])
    meanEBVs             <- mean(phenoDF_TS$EBVs[rows])
    meanBVs              <- mean(phenoDF_TS$True_BVs[rows])
    phenoDF_TS$Pheno[rows]    <- phenoDF_TS$Pheno[rows]    - meanPheno
    phenoDF_TS$EBVs[rows]     <- phenoDF_TS$EBVs[rows]     - meanEBVs
    phenoDF_TS$True_BVs[rows] <- phenoDF_TS$True_BVs[rows] - meanBVs
  }
  
  # True within-family h2
  h2_2 <- var(phenoDF_TS$True_BVs) / var(phenoDF_TS$Pheno)
  
  accSel_fam <- cor(phenoDF_TS$True_BVs, phenoDF_TS$EBVs)
  #modAcc_fam <- sqrt(h2_mod2)
  
  # Store within-family accuracies
  #output$accSel_fam[year] <- accSel_fam
  #output$modAcc_fam[year] <- modAcc_fam
  
  # Store within-family heritabilities
  #output$h2_2[year]     <- h2_2
  #output$h2_mod2[year]  <- h2_mod2  
  
  # ---- Recycling 10 parents from PYT (based on EBV) ----
  newParents <- selectInd(PYT, 10, use = "ebv")
  Parents    <- c(Parents[11:nParents], newParents)
  
  # ---- Advance Year ----
  # Stage 6: Release variety (no code here)
  
  # Stage 5 (Elite yield trials (EYT))
  AYT@ebv <- setGEAVs(AYT, marker_eff)
  EYT <- selectInd(AYT, nEYT, use= "ebv")
  EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
  
  # Stage 4 (Advanced yield trials (AYT))
  
  AYT <- selectInd(PYT, nAYT, use = "ebv")
  AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
  
  # Stage 3 (Preliminary yield trials (PYT))
  DH@ebv <- setGEAVs(
    DH, marker_eff
  )
  
  #output$accSel_True[year] <- cor(DH@gv, DH@ebv)
  
  PYT <- selectWithinFam(DH, famMax, use = "ebv")
  PYT <- setPheno(PYT, varE = VarE, reps = repPYT)
  
  # Stage 2 (Generating DHs)
  DH <- makeDH(F1, nDH)
  
  # ---- Update Global Parameters ----
  #source(file = "GlobalParameters/GlobalParameters_Control.R")
  
  # Stage 1 (Crossing block)
  F1 <- randCross(Parents, nCrosses)
  
  # ---- Store Training Population ----
  PYT@fixEff <- as.integer(rep(year, nInd(PYT)))
  AYT@fixEff <- as.integer(rep(year, nInd(AYT)))
  EYT@fixEff <- as.integer(rep(year, nInd(EYT)))
  
  if (year > nBurnin) {
    cat("  Maintaining training population \n")
    TrainPop <- c(TrainPop[-c(1:c(PYT, EYT, AYT)@nInd)], PYT, EYT, AYT)
  }
  
  # ---- Store DH Metrics ----
  output$meanG_yield[year]     <- meanG(DH)[1]
  output$meanG_protein[year]     <- meanG(DH)[2]
  
  output$varG_yield[year]      <- varG(DH)[1]
  output$varG_protein[year]      <- varG(DH)[2,2]
  
  output$genicVarG_yield[year] <- genicVarG(DH)[1]
  output$genicVarG_protein[year] <- genicVarG(DH)[2]
  
  
}

write.table(
  output[output$rep == REP & output$year > nBurnin, ],
  file      = output_csv_gs,
  sep       = ",",
  col.names = FALSE,
  row.names = FALSE,
  append    = TRUE
)

cat("PS_RM phase completed and results saved in '", output_csv_gs, "'.\n")