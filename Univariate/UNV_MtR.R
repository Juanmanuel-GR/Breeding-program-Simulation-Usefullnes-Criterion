# Future Scenario with MateR

#clean environment and load packages 
rm(list = ls())
library(AlphaSimR)
library(MateR)
library(lme4breeding)
library(tidyr)
library(dplyr)


REP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "UNV_MtR_2F"
#REP <- 1

#load the globar parameters, 

source(file = "GlobalParams/Globa_params.R")
#source(file = "FutureParams/GMMT_params.R")

#functions 
source(file = "GlobalParams/sim_functions_mult.R")

# --- Input Parameters ---
burin_dir <- "../../Multivariate/Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
output_csv_2f <- file.path(results_dir, paste0("UNV_MtR_Results_", REP, ".csv"))

# ---- Create Results Directory ----
if (!dir.exists(results_dir)) dir.create(results_dir)

# ---- Initialize CSV File ----
if (file.exists(output_csv_2f)) file.remove(output_csv_2f) # Remove existing file to avoid appending


mating_dir <- paste0("mating_dir_", SCENARIO)

mating_file <- file.path( mating_dir, paste0("Mating_Plan_Total", REP, ".csv"))

# ---- Create Results Directory ----
if (!dir.exists(mating_dir)) dir.create(mating_dir)

# ---- Initialize CSV File ----
if (file.exists(mating_file)) file.remove(mating_file) # Remove existing file to avoid appending

header <- c(
  "rep", "year", "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein",
  "heritability", "var_FamMean", "var_FamDev", "var_ratio", "UI", "genicSD_lost"
)


write.table(t(header), file = output_csv_2f, sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE)


#getwd()
#load(file = "../BurnInGen/burning_states/BurnIn_REP_1.RData")

#---- Future GS Phase ----
cat("Processing the future scenario of REP", REP, "Scenario", SCENARIO, "\n")

#Load the burn-in state
burnin_file <- file.path(burin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)

#output$accSel <- NULL

# Ensure 'output'NULL# Ensure 'output' data frame is extended for GS phase

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
    heritability = numeric(nfuture),
    var_FM = numeric(nfuture),
    var_Fdev = numeric(nfuture),
    var_ratio = numeric(nfuture),
    UI = numeric(nfuture),
    genic_SDlost = numeric(nfuture)
  ))



mating_plan_out <- data.frame(Rep = REP,
                              year = numeric(nCrosses*nfuture),
                              Parent1 = numeric(nCrosses*nfuture),
                              Parent2 = numeric(nCrosses*nfuture),
                              NumberOfCrosses = numeric(nCrosses*nfuture),
                              FM = numeric(nCrosses*nfuture),
                              genetic_deviation = numeric(nCrosses*nfuture),
                              usefulness = numeric(nCrosses*nfuture))

# Future GS simulation loop
for (year in (nBurnin+1):(nBurnin+nfuture)) {
  cat("  Future year:", year, "\n")
  
  
  #genomic Selection
  M_all <- pullSnpGeno(TrainPop) # extract the marker matrix
  ids_all <- TrainPop@id 
  id_unique <- as.integer(unique(ids_all))
  posiciones_unicas <- match(id_unique, ids_all)
  M_train <- M_all[posiciones_unicas, ] #filter to have only once each genotype, maybe thats what filter_zero_variance does
  rm(M_all)
  #Filter 0 variance Markers
  M_train <- filter_zero_variance(M_train)
  
  #Compute the VanRaden additive GRM
  GRM <- getVR_GRM(M_train)
  GRM <- GRM +diag(1e-6, nrow(GRM)) #have the relationship matrix invertible
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
  
  h2_1 <- var(PYT@gv) / var(PYT@pheno)
  
  
  # Store general heritabilities
  output$heritability[year] <- h2_mod
  #output$h2_1[year]    <- h2_1
  
  ################################################################################
  # Extract additive marker effects
  ################################################################################
  
  marker_eff <- extract_meff(M_train, GEBVs)
  trait_name   <- SP$traitNames[[1]]
  colnames(marker_eff) <- trait_name
  
  #Predict EBVs and select new Parents, we predict on the EBV object of the population
  
  PYT@ebv <- setGEAVs(PYT, marker_eff)
  
  modAcc <- sqrt(h2_mod)
  
  # Storing general model accuracies
  #output$accSel[year] <- cor(PYT@gv, PYT@ebv)
  #utput$modAcc[year] <- modAcc
  
  PYT@misc$Family <- factor(paste(PYT@father, PYT@mother, sep = "/"))
  
  # Creating dataframe for within family accuracy:
  phenoDF_TS <- data.frame(
    GID = PYT@id,
    Pheno = c(PYT@pheno),
    EBVs = c(PYT@ebv),
    Family = PYT@misc$Family,
    True_BVs = c(PYT@gv)
  )
  
  # Centering different values by the family mean for PYT population
  for (family in unique(phenoDF_TS$Family)) {
    rows <- which(phenoDF_TS$Family == family)
    meanPheno <- mean(phenoDF_TS$Pheno[rows])
    meanEBVs <- mean(phenoDF_TS$EBVs[rows])
    meanBVs <- mean(phenoDF_TS$True_BVs[rows])
    phenoDF_TS$Pheno[rows] <- phenoDF_TS$Pheno[rows]-meanPheno
    phenoDF_TS$EBVs[rows] <- phenoDF_TS$EBVs[rows]-meanEBVs
    phenoDF_TS$True_BVs[rows] <- phenoDF_TS$True_BVs[rows]-meanBVs
  }
  
  h2_2 <- var(phenoDF_TS$True_BVs)/var(phenoDF_TS$Pheno) # True within family h2
  
  accSel_fam <- cor(phenoDF_TS$True_BVs, phenoDF_TS$EBVs)
  #modAcc_fam <- sqrt(h2_mod2)
  
  # Storing within family accuracies
  #output$accSel_fam[year] <- accSel_fam
  #output$modAcc_fam[year] <- modAcc_fam
  
  # Storing within family heritabilities:
  #output$h2_2[year] <- h2_2
  #output$h2_mod2[year] <- h2_mod2
  
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
  DH@ebv <- setGEAVs(DH,
                     marker_eff)
  
  #output$accSel_True[year] <- cor(DH@gv, DH@ebv)
  
  PYT <- selectWithinFam(DH, famMax, use = "ebv")
  PYT <- setPheno(PYT, varE = VarE, reps = repPYT)
  
  # Stage 2 (Generating DHs)
  DH <- makeDH(F1, nDH)
  
  # MateR runner part --------------------------------------------------------
  cat("  Finding an optimal mating plan...\n")
  
  source("GlobalParams/GMMT_params.R")
  
  ParentalPool <- c(Parents, PYT, AYT, EYT) #crear una pool d padres nueva, donde meteremos : Parents, PYT, AYT, EYT
  print(length(ParentalPool@id))
  ParentalPool <- ParentalPool[unique(ParentalPool@id)] # comprobar qu solo hay padres unique
  
  Parental_ids <- ParentalPool@id
  newP_Markers <- pullSnpGeno(ParentalPool)
  filt_newA_Markers <- newP_Markers[, colnames(newP_Markers) %in% colnames(M_train), drop = FALSE]
  marker_names <- colnames(filt_newA_Markers)
  
  # Compute a recombination frequency matrix
  #para predecir como sera una fam, se necesita saber q tan probable es q dos genes se hereden juntos o haya recombinacion entre ellos 
  compute_recombination_matrix <- function(physical_positions, marker_names, chromosome_length_cM) {
    # Extract chromosome numbers from marker names
    chromosomes <- as.integer(sub("_.*", "", marker_names))  # Extract portion before "_", suele ser el marker name 1_12434, siendo el 1 el numero del croms 
    
    # Scale physical positions dynamically
    genetic_distances_cM <- physical_positions / max(physical_positions) * chromosome_length_cM
    
    n_markers <- length(genetic_distances_cM)
    rec_matrix <- matrix(0, nrow = n_markers, ncol = n_markers)
    colnames(rec_matrix) <- marker_names
    rownames(rec_matrix) <- marker_names
    
    # Haldane mapping function
    haldane <- function(d) 0.5 * (1 - exp(-2 * d))
    
    # Compute pairwise recombination fractions
    for (i in 1:(n_markers - 1)) {
      for (j in (i + 1):n_markers) {
        if (chromosomes[i] != chromosomes[j]) {
          # Markers on different chromosomes segregate independently
          rec_matrix[i, j] <- rec_matrix[j, i] <- 0.5
        } else {
          # Markers on the same chromosome
          d <- abs(genetic_distances_cM[i] - genetic_distances_cM[j]) / 100  # Convert to Morgans
          rec_matrix[i, j] <- rec_matrix[j, i] <- haldane(d)
        }
      }
    }
    
    # Ensure the diagonal remains 0 (no recombination within the same marker)
    diag(rec_matrix) <- 0
    
    return(rec_matrix)
  }
  
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
  
  #h2 <- modAcc_fam
  
  #if (h2 < 0) {
    #h2 <- 0
  #}
  
  h2 <- 0  #esto es para q usefulness sea solo el family mean, mirare a ver si lo tengo q cambiar para q usefulness sea mas cosas
  
  crossingBlockT <- 0.5/100

  #head(pheno(ParentalPool))
  
  #head(marker_eff)
  #marker_eff <- list(marker_eff)
  marker_eff_list <- list(c(marker_eff))
  
  names(marker_eff_list[[1]]) <- rownames(marker_eff)
  names(marker_eff_list) <- colnames(marker_eff)
  # Running MateR:
  out <- GenomicMatingMT(Parents1 = Parental_ids,
                         Parents2 = Parental_ids,
                         parametrization = "Genotypic",
                         opt_type = "Greedy",
                         Markers = filt_newA_Markers,
                         phi = phi,                           
                         markereffects = marker_eff_list,
                         n = selfing_cycles,                  
                         PropSD = crossingBlockT,
                         size = nCrosses,
                         c_list = c_list,
                         coefficients = NULL,#mult_trait_coefficients,
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
  con_val <- out$TrainSelout$convergence
  cat("  An optimal plan was found with", con_val, "convergence.\n")
  
  
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
    
  

  
  # all unique parent IDs used in the mating plan
  uniqueParents <- unique(as.vector(crossPlan))
  Parents_w <- ParentalPool[match(as.vector(crossPlan), ParentalPool@id)]
  
  genicSD_C <- sqrt(as.numeric(genicVarG(ParentalPool)))
  #geneticSD_C <- sqrt(as.numeric(varG(ParentalPool)))
  
  genicSD_Pw <- sqrt(as.numeric(genicVarG(Parents_w)))
  #geneticSD_Pw <- sqrt(as.numeric(varG(Parents_w)))
  
  genicSD_lost_CPw <- 1 - (genicSD_Pw/genicSD_C)
  #geneticSD_lost_CPw <- 1 - (geneticSD_Pw/geneticSD_C)
  
  output$genic_SDlost[year] <- genicSD_lost_CPw[1]
  
  

  
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
            file = output_csv_2f, sep = ",", col.names = FALSE, row.names = FALSE, append = TRUE)

cat("MtR phase completed and results saved in '", output_csv_2f, "'.\n")




