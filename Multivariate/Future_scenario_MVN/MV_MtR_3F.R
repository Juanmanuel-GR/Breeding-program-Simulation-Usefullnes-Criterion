# Multitrait Genomic Selection, I'll be using sommer 

# ---- Clean Environment and Load Packages ----
rm(list = ls())
library(AlphaSimR)
library(sommer)
library(tidyr)
library(dplyr)
library(MateR)


REP      <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "MVN_MtR_3F"


#try 

# ---- Helper Functions ---- poner aqui las mias
source(file = "Global_Params/sim_functions_mult.R")
source(file = "Global_Params/global_param_mult.R")

#-----Input parameters------
burnin_dir <- "../Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
MtR_3F_csv <- file.path(results_dir, paste0("MVN_MtR_3F_Results_", REP, ".csv"))

#---- Create Results Directory ----
if (dir.exists(results_dir)) {
  dir.create(results_dir)
}

# ---- Initializa CSV File -----
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

mating_dir <- paste0("mating_dir_", SCENARIO)

mating_file <- file.path( mating_dir, paste0("Mating_Plan_Total", REP, ".csv"))

# ---- Create Results Directory ----
if (!dir.exists(mating_dir)) dir.create(mating_dir)

# ---- Initialize CSV File ----
if (file.exists(mating_file)) file.remove(mating_file) # Remove existing file to avoid appending


correlation_dir <- paste0("correlation_dir_", SCENARIO)

correlation_file <- file.path(correlation_dir, paste0("Correlation_true_Total", REP, ".csv"))

# ---- Create Results Directory ----
if (!dir.exists(correlation_dir)) dir.create(correlation_dir)

# ---- Initialize CSV File ----
if (file.exists(correlation_file)) file.remove(correlation_file)


header <- c(
  "rep", "year", "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein",
  "weighted_heritability", "var_FamMean", "var_FamDev", "var_ratio", "UI", "genicSD_lost_yield", "genic_SD_lost_protein"
)

write.table(
  t(header),
  file = MtR_3F_csv,
  sep = ",",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)




burnin_file <- file.path(burnin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)



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
    weighted_heritability = numeric(nfuture),
    var_FM = numeric(nfuture),
    var_Fdev = numeric(nfuture),
    var_ratio = numeric(nfuture),
    UI = numeric(nfuture),
    genic_SDlost_yield = numeric(nfuture),
    genic_SDlost_protein = numeric(nfuture)
  ))


mating_plan_out <- data.frame(Rep = REP,
                              year = numeric(nCrosses*nfuture),
                              Parent1 = numeric(nCrosses*nfuture),
                              Parent2 = numeric(nCrosses*nfuture),
                              NumberOfCrosses = numeric(nCrosses*nfuture),
                              FM = numeric(nCrosses*nfuture),
                              genetic_deviation = numeric(nCrosses*nfuture),
                              usefulness = numeric(nCrosses*nfuture))


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
  
  
  
  
  
  cat(" Fitting the mmer model \n")
  
  model <- mmer(cbind(Yield, Protein)~ YearEff,
                random = ~vsr(GID, Gu=GRM, Gtc=unsm(2)),
                rcov = ~vsr(units, Gtc=unsm(2)),
                data = y_train, verbose = TRUE)
  
  
  
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
  
  # MateR runner part --------------------------------------------------------
  cat("  Finding an optimal mating plan...\n")
  
  source("Global_Params/GMMT_params.R")
  
  ParentalPool <- c(Parents, PYT, AYT, EYT) #crear una pool d padres nueva, donde meteremos : Parents, PYT, AYT, EYT
  print(length(ParentalPool@id))
  ParentalPool <- ParentalPool[unique(ParentalPool@id)] # comprobar qu solo hay padres unique
  
  Parental_ids <- ParentalPool@id
  newP_Markers <- pullSnpGeno(ParentalPool)
  filt_newA_Markers <- newP_Markers[, colnames(newP_Markers) %in% colnames(M_train), drop = FALSE]
  marker_names <- colnames(filt_newA_Markers)
  
  
  # Data prep for computing the c_matrix
  physical_positions <- as.integer(sapply(strsplit(marker_names, "_"), `[`, 2))
  
  # Compute recombination frequency matrix with improved genetic distance calculation
  c_matrix <- compute_recombination_matrix(
    physical_positions = physical_positions,
    marker_names = marker_names,
    chromosome_length_cM = chromosome_length_cM
  )
  
  # ---------- Creating c_list ---------:
  
  # Extract chromosome numbers from marker names
  chromosomes <- as.integer(sub("_.*", "", rownames(c_matrix)))
  
  # Identify unique chromosomes
  unique_chromosomes <- sort(unique(chromosomes))
  
  # Split c_matrix into a list of matrices by chromosome
  c_list <- lapply(unique_chromosomes, function(chr) {
    markers_on_chr <- rownames(c_matrix)[chromosomes == chr]
    c_matrix[markers_on_chr, markers_on_chr, drop = FALSE]
  })
  
  # Assign names to the list elements
  names(c_list) <- paste0("Chr_", unique_chromosomes)
  
  rm(c_matrix)
  
  h2 <- Acc_sel

  #h2 <- 0.5
  if (h2 < 0) {
  h2 <- 0
  }
  
  #h2 <- 0  #esto es para q usefulness sea solo el family mean, mirare a ver si lo tengo q cambiar para q usefulness sea mas cosas
  
  crossingBlockT <- 0.5/100
  
  #head(pheno(ParentalPool))
  
  #head(marker_eff)
  #marker_eff <- list(marker_eff)
  #marker_eff[,1]
  markereffects <- list(Yield = marker_eff[,1],
                        Protein = marker_eff[,2] )
  markereffects$Yield[1:4]
  
  
  
  # Supongamos que tu lista tiene dos elementos
  names(markereffects) <- c("Yield", "Protein")
  markereffects$Yield
  
  
  coefficients_vector <- as.numeric(b)
  
  # Le ponemos exactamente los mismos nombres
  names(coefficients_vector) <- c("Yield", "Protein")
  # Running MateR:
  
  coefficients_vector
  names(markereffects)
  names(coefficients_vector)
  out <- GenomicMatingMT(Parents1 = Parental_ids,
                         Parents2 = Parental_ids,
                         parametrization = "Genotypic",
                         opt_type = "Greedy",
                         Markers = filt_newA_Markers,
                         phi = phi,                           
                         markereffects = markereffects,
                         n = selfing_cycles,                  
                         PropSD = crossingBlockT,
                         size = nCrosses,
                         c_list = c_list,
                         coefficients = coefficients_vector,#mult_trait_coefficients,
                         offspring_per_cross = nDH,
                         within_family_accuracy = h2,   
                         control = control,
                         n_selected_per_family= famMax,
                         n_families = n_families,             
                         Username= username, 
                         Password= password, # To obtain license keys, please contact javier.fgonzalez@upm.es or j.isidro@upm.es
                         replication = replication,
                         DHs = DHs,
                         LD = LD,
                         Username_TrainSel = username_trainsel, # To obtain license keys, please contact javier.fgonzalez@upm.es or j.isidro@upm.es
                         Password_TrainSel = password_trainsel
  )
  
  # Extract the opt_matingPlan from GMMT
  opt_matingPlan <- out$OptimalMatingScheme
  # Expand rows based on Number_Of_Crosses
  expandedPlan <- opt_matingPlan[rep(1:nrow(opt_matingPlan), opt_matingPlan$Number_Of_Crosses), ]
  rownames(expandedPlan) <- NULL
  # Create crossPlan matrix for AlphaSimR
  crossPlan <- as.matrix(expandedPlan[, c("Parent1", "Parent2")])
  # 0 means no convergence
  
  #Retrieve Params for the UI
  Family_names <- names(out$FamilyValues$mean)
  parents <- do.call(rbind, strsplit(Family_names, "/"))
  means <- as.numeric(out$FamilyValues$mean[Family_names])
  sds <- as.numeric(out$FamilyValues$sd[Family_names])
  
  UC <- means + compute_i(famMax/nDH) *sds * h2
  
  # Create dataframe
  UC_df <- data.frame(
    Parent1 = parents[, 1],
    Parent2 = parents[, 2],
    Mean    = means,
    SD      = sds,
    UC      = UC,
    stringsAsFactors = FALSE
  )
  
  UC_df <- UC_df[order(UC_df$UC, decreasing = TRUE), ]
  rownames(UC_df) <- NULL
  
  varM      <- var(UC_df$Mean)
  output$var_FM[year] <- varM
  
  varSD     <- var(UC_df$SD)
  output$var_Fdev[year] <- varSD
  
  RatioM_SD <- varM / varSD
  output$var_ratio[year] <- RatioM_SD
  
  # varRatio = var(Sigma) / var(Mu)
  rho_mu_UC <- function(varRatio, i, r) {
    1 / sqrt(1 + (i * r)^2 * (1 / varRatio))
  }
  
  rhoValue <- rho_mu_UC(
    varRatio = RatioM_SD,
    i        = compute_i(famMax / nDH),
    r        = h2
  )
  
  UI <- 1 / rhoValue
  output$UI[year] <- UI
  
  
  mating_plan_anio <- data.frame(
    rep               = REP,
    Year              = year,
    Parent1           = out$OptimalMatingScheme$Parent1,
    Parent2           = out$OptimalMatingScheme$Parent2,
    NumberOfCrosses   = out$OptimalMatingScheme$Number_Of_Crosses,
    FM                = out$OptimalMatingScheme$Family_average,
    genetic_deviation = out$OptimalMatingScheme$Deviation_From_Average_Selected_Best,
    usefulness        = out$OptimalMatingScheme$Usefulness
  )
  
  # 5. GUARDAR AL CSV (Modo Append)
  # col.names solo será TRUE si el archivo NO existe aún (para poner la cabecera una sola vez)
  write.table(mating_plan_anio, 
              file = mating_file, 
              sep = ",", 
              append = TRUE, 
              quote = FALSE,
              col.names = !file.exists(mating_file), 
              row.names = FALSE)
  
  
  
  
  
  # ---- Update Global Parameters ----
  #source(file = "GlobalParameters/GlobalParameters_Control.R")
  
  # all unique parent IDs used in the mating plan
  uniqueParents <- unique(as.vector(crossPlan))
  Parents_w <- ParentalPool[match(as.vector(crossPlan), ParentalPool@id)]
  
  genicSD_C <- sqrt(as.numeric(genicVarG(ParentalPool)))
  #geneticSD_C <- sqrt(as.numeric(varG(ParentalPool)))
  
  genicSD_Pw <- sqrt(as.numeric(genicVarG(Parents_w)))
  #geneticSD_Pw <- sqrt(as.numeric(varG(Parents_w)))
  
  genicSD_lost_CPw <- 1 - (genicSD_Pw/genicSD_C)
  #geneticSD_lost_CPw <- 1 - (geneticSD_Pw/geneticSD_C)
  
  output$genic_SDlost_yield[year] <- genicSD_lost_CPw[1]
  output$genic_SDlost_protein[year] <- genicSD_lost_CPw[2]
  #output$genetic_SDlost[year] <- geneticSD_lost_CPw
  
  #correlation part
  
  # I only simulate the progeny for all possible crosses on year 21, due to computational time, but it could be done untill the end of the breeding program
  if (year == 21) {
    
    results <- data.frame(
      family    = family_vect,
      est_mean  = unlist(out$FamilyValues$mean[family_vect]),
      est_sd    = unlist(out$FamilyValues$sd[family_vect]),
      true_mean = NA_real_,
      true_sd   = NA_real_,
      ebv_mean = NA_real_
    )
    
    find_stage <- function(id, stages, stage_names) {
      for (i in seq_along(stages)) {
        if (id %in% stages[[i]]@id) {
          return(stage_names[i])
        }
      }
      return("not found")
    } 
    
    get_ebv_from_stages <- function(id, stages, stage_names) {
      stage_name <- find_stage(id, stages, stage_names)
      if (stage_name == "not found") {
        warning(paste("ID", id, "not found in any stage"))
        return(NA)
      }
      stage <- stages[[which(stage_names == stage_name)]]
      return(stage@ebv[stage@id == id, , drop = FALSE])
    }
    
  
    stages        <- list(EYT, AYT, PYT, Parents)
    stage_names   <- c("EYT", "AYT", "PYT", "Parents")
    
    
    for (i in seq_along(family_vect)) {
      
      ids <- strsplit(family_vect[i], "/")[[1]]
      crossing_pop <- ParentalPool[ParentalPool@id %in% ids]
      crossing_pop@ebv <- setGEAVs(crossing_pop, marker_eff)
      
      crossPlan <-matrix(ids, nrow = 1, dimnames = list(NULL, c("Parent1", "Parent2")))
      F1 <- makeCross(crossing_pop,
                      crossPlan = crossPlan,
                      nProgeny  = 10,
                      simParam  = SP)
      
      DH <- makeDH(F1, nDH = 50, simParam = SP)
      
      
      
      results$true_mean[i] <- meanG(DH)[1]
      results$true_sd[i]   <- sqrt(varG(DH)[1,1])
      
      
      
      cat("Cross", i, "/", length(family_vect), ":", family_vect[i], "\n")
      
    }
    results$RatioM_SD = var(results$true_mean) / var(results$true_sd)
    results$rep <- REP
    
    write.csv(results, file = correlation_file, row.names = FALSE)
    
  }
  
  
  
  
  
  # Optimal crosses ---------------------------------------------------------
  F1 <- makeCross(ParentalPool, crossPlan = crossPlan)
  
  # all unique parent IDs used in the mating plan
  uniqueParents <- unique(as.vector(crossPlan))
  
  Parents <- ParentalPool[ParentalPool@id %in% uniqueParents]
  
  # ---- Store Training Population ----
  PYT@fixEff <- as.integer(rep(year,nInd(PYT)))
  AYT@fixEff <- as.integer(rep(year,nInd(AYT)))
  EYT@fixEff <- as.integer(rep(year,nInd(EYT)))
  
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

# Append results to CSV file
write.table(output[output$rep == REP & output$year > nBurnin, ],
            file = MtR_3F_csv, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

cat("MtR phase completed and results saved in '", MtR_3F_csv, "'.\n")
