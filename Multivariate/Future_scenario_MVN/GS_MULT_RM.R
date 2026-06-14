# Multitrait Genomic Selection, I'll be using sommer 

# ---- Clean Environment and Load Packages ----
rm(list = ls())
library(AlphaSimR)
library(sommer)
library(tidyr)
library(dplyr)



REP      <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "GS_MVN_RM"

y 

# ---- Helper Functions ---- 
source(file = "Global_Params/sim_functions_mult.R")
source(file = "Global_Params/global_param_mult.R")

#-----Input parameters------
burnin_dir <- "../Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
gs_csv <- file.path(results_dir, paste0("GS_MVN_RM_Results_", REP, ".csv"))

#---- Create Results Directory ----
if (dir.exists(results_dir)) {
  dir.create(results_dir)
}

# ---- Initializa CSV File -----
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

header <- c(
  "rep", "year", "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein", "weighted_heritability"
)

write.table(
  t(header),
  file = gs_csv,
  sep = ",",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

burnin_file <- file.path(burnin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)



output <- rbind(
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
    weighted_heritability = numeric(nfuture)
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
  
  #Compute the VanRaden additive GRM, using sommer
  GRM <- getVR_GRM(M_train)
  GRM <- GRM +diag(1e-6, nrow(GRM))
  #Create the dataframe for the phenotype
  phenotype <- pheno(TrainPop)
  fixed <- TrainPop@fixEff # year fixed effect
  
  y_train <- data.frame(GID=ids_all,
                        Yield = as.numeric(phenotype[,1]),
                        Protein = as.numeric(phenotype[,2]))
  #add fixed effect as factor
  y_train$YearEff <- as.factor(fixed)
  

  
  
  # I need to transform from wide format to long format
  #id_samp <- sample(id_unique, size = length(id_unique)*0.3, replace=FALSE)
  
  #GRM_samp <- GRM[rownames(GRM)%in% id_samp, colnames(GRM)%in% id_samp ]
  
  #y_train_samp <- y_train[y_train$GID %in% id_samp,]
  
  #M_train_samp <- M_train[rownames(M_train)%in% id_samp,]
  cat(" Fitting the mmer model \n")
  starting_time <- Sys.time()
  model <- mmer(cbind(Yield, Protein)~ YearEff,
                random = ~vsr(GID, Gu=GRM, Gtc=unsm(2)),
                rcov = ~vsr(units, Gtc=unsm(2)),
                data = y_train, verbose = TRUE)
  
  end_time <- Sys.time()
  time_took <- end_time - starting_time
  cat("Model Converged, took" , time_took, "\n")
  
  
  #first transform the dataset we used in sommer `mmer` and transform it in longformat
  pheno_long <- y_train %>% pivot_longer(cols = c("Yield", "Protein"), names_to = "Trait", values_to = "value") %>%
    arrange(desc(Trait)) %>% mutate(across(c(Trait, GID), as.factor))
  #have name of the traits 
  traits <- unique(pheno_long$Trait)
  
  n <- n_distinct(pheno_long$GID) #number of disting genotypes
  t <- n_distinct(pheno_long$Trait) #number of traits
  ids <- unique(pheno_long$GID) #ids of genotypes
  
  names_G <- as.vector(outer(traits, ids, paste, sep = ":")) #maybe not need
  
  sigma_g <- model$sigma$`u:GID`
  print(sigma_g)
  sigma_e <- model$sigma$`u:units`
  
  G <- GRM %x% sigma_g #kronecker product
  R <- diag(nrow(y_train)) %x% sigma_e #
  Z <- model.matrix(~ -1 + GID:Trait, pheno_long) #Z design matrix
  X <- model.matrix(~ Trait, pheno_long) #X design matrix
  
  
  var_uhat <- get_marginal_variance(Z, G, R )
  
  b <- compute_coefficients_blups(var_uhat, y_train)
  
  print("Coefficients retrive")
  print(b)
  #genomic selection using sommer
  #starting_time <- Sys.time()
  #cat("Runnin sommer (mmes) models...\n")

  #compute the index heritability 
  
  sigma_p <- sigma_g + sigma_e
  
  # trait order in theta is: Protein, Yield
  w <- matrix(b, ncol = 1)
  
  h2_index <- as.numeric(t(w) %*% sigma_g %*% w / (t(w) %*% sigma_p %*% w))
  h2_index
  output$weighted_heritability[year] <- h2_index
  Acc_sel <- sqrt(h2_index)
  
  
  
  
  ################################################################################
  # Extract estimates from model 
  ################################################################################
  
  #GEBVs - the BLUPS
  GEBVs <- data.frame(model$U$`u:GID`)
  
  #recover the same order
  GEBVs <- GEBVs[rownames(GRM),]
  #hist(as.numeric(GEBVs[,1]))
  #hist(as.numeric(GEBVs[,2]))
  
  #head(GEBVs[,1])
  
  
  
  
  ################################################################################
  # Extract additive marker effects
  ################################################################################
  
  marker_eff <- extract_meff(M_train, GEBVs)
  dim(M_train)
  dim(GEBVs)
  #trait_name   <- SP$traitNames[[2]]
  #colnames(marker_eff) <- trait_name

  #Predict EBVs and select new Parents, we predict on the EBV object of the population
  
  PYT@ebv <- setGEAVs(PYT, marker_eff)
  head(PYT@ebv)
  #modAcc <- sqrt(h2_mod)
  
  # Store general model accuracies
  #output$accSel[year] <- cor(PYT@gv, PYT@ebv)
  #output$modAcc[year] <- modAcc
  
  PYT@misc$Family <- factor(paste(PYT@father, PYT@mother, sep = "/"))
  
  #pheno(PYT)[,2]
  # Create dataframe for within-family accuracy
  phenoDF_TS <- data.frame(
    GID      = PYT@id,
    Pheno    = round(c(PYT@pheno[,2]), digits = 5),
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
  
  # ---- Recycling 10 parents from PYT (based on EBV) ----
  #newParents <- selectInd(PYT, 10, use = "ebv")
  newParents <- geno_selection_index(PYT, b, 10)
  Parents    <- c(Parents[11:nParents], newParents)
  
  # ---- Advance Year ----
  # Stage 6: Release variety (no code here)
  
  # Stage 5 (Elite yield trials (EYT))
  AYT@ebv <- setGEAVs(AYT, marker_eff)
  EYT <- geno_selection_index(AYT, b, nEYT)
  EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
  
  # Stage 4 (Advanced yield trials (AYT))
  PYT@ebv <- setGEAVs(PYT, marker_eff)
  AYT <- geno_selection_index(PYT, b, nAYT)
  AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
  
  # Stage 3 (Preliminary yield trials (PYT))
  DH@ebv <- setGEAVs(
    DH, marker_eff
  )
  
  #output$accSel_True[year] <- cor(DH@gv, DH@ebv)
  PYT <- select_index_within(DH, b, famMax, pattern = "ebv")
  #PYT <- selectWithinFam(DH, famMax, use = "ebv")
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

# ---- Append Results to CSV ----
write.table(
  output[output$rep == REP & output$year > nBurnin, ],
  file      = gs_csv,
  sep       = ",",
  col.names = FALSE,
  row.names = FALSE,
  append    = TRUE
)

cat("GS_MVN_RM phase completed and results saved in '", gs_csv, "'.\n")
