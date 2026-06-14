#simulation where i simulate a population with multitrait

library(AlphaSimR)
library(sommer)
library(tidyr)
library(dplyr)
library(AGHmatrix)


rm(list = ls())




#load the global Parameters 
source(file = "Global_Params/global_param_mult.R")

#load the functions needed for the index selection
source(file = "Global_Params/sim_functions_mult.R")

#Retrieve Parameters -----
REP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "PS_MULT_RM"
#REP <- 1

burnin_dir <- "../Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
pheno_csv <- file.path(results_dir, paste0("PS_MULT_RM_Results", REP, ".csv"))


# Create Results Directory
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#Initialize CSV file
if(file.exists(pheno_csv)) file.remove(pheno_csv)

#create the header for the csv file
header <- c("rep",
            "year",
            "meanG_yield", 
            "meanG_protein",
            "varG_yield",
            "varG_protein",
            "genicVarG_yield",
            "genicVarG_protein")


write.table(t(header),
            file = pheno_csv,
            sep = ",",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

#---Future PS Phase 
cat("Processing the future scenario of REP", REP, "SCENARIO", SCENARIO, "\n")
getwd()
#load the burnIn state data
burnin_file <- file.path(burnin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)



# ---- Extend 'output' for PS phase ----
output <- bind_rows(
  output,
  data.frame(
    rep               = REP,
    year              = (nBurnin + 1):(nBurnin + nfuture),
    meanG_yield       = numeric(nfuture),
    meanG_protein     = numeric(nfuture),
    varG_yield        = numeric(nfuture),
    varG_protein      = numeric(nfuture),
    genicVarG_yield   = numeric(nfuture),
    genicVarG_protein = numeric(nfuture)
  ))

# Future PS Simulation Loop ----
for (year in (nBurnin+1):(nBurnin + nfuture)) {
  cat( " Future year:", year, "\n")
  
  #current_row <- year - nBurnin
  # ---- Update Parents ----
  # Replace 10 oldest inbred parents with 10 new parents from EYT stage
  Parents <- c(Parents[11:nParents], EYT)
  
  #genomic Selection
  #M_all <- pullSnpGeno(TrainPop)
  ids_all <- TrainPop@id
  id_unique <- as.integer(unique(ids_all))
 
  
  
  print(length(TrainPop@id))
  phenotype <- pheno(TrainPop)
  fixed <- TrainPop@fixEff # year fixed effect
  
  y_train <- data.frame(GID=ids_all,
                        Yield = as.numeric(phenotype[,1]),
                        Protein = as.numeric(phenotype[,2]))
  #add fixed effect as factor
  y_train$YearEff <- as.factor(fixed)
  
  
  # I need to transform from wide format to long format
  #id_samp <- sample(id_unique, size = length(id_unique)*0.5, replace=FALSE)
  
  #y_train_samp <- y_train[y_train$GID %in% id_samp,]
  # 2. Calcular la Matriz A
  # ploidy=2 para diploides (trigo/cebada)
  ped_sommer <- master_ped[, c("id", "mother", "father")]

 
  if (year <= 30) { #from year 30 bc pedigree matrix consumes lot of ram memory, we start replacing individuals from this year, typical in real life
    ped_sommer$mother[is.na(ped_sommer$mother)] <- 0
    ped_sommer$father[is.na(ped_sommer$father)] <- 0
    
    A_matrix <- Amatrix(data = ped_sommer, 
                        ploidy = 2, 
                        verify = TRUE)
    
  } else{
    cont_recorte <- year - 30
    n_delete <- 2860 * cont_recorte
    
    ped_sommer <- ped_sommer[-c(1:n_delete),]
    
    valid_ids <- ped_sommer$id
    ped_sommer$mother <- ifelse(ped_sommer$mother %in% valid_ids, 
                                ped_sommer$mother, 
                                NA_character_)
    
    
    ped_sommer$father <- ifelse(ped_sommer$father %in% valid_ids, 
                                ped_sommer$father, 
                                NA_character_)
    
    ped_sommer$mother[is.na(ped_sommer$mother)] <- 0
    ped_sommer$father[is.na(ped_sommer$father)] <- 0
    
    A_matrix <- Amatrix(data = ped_sommer, 
                        ploidy = 2, 
                        verify = TRUE)
    
    
 
    
  }

  
  A_matrix_fin <- A_matrix[rownames(A_matrix) %in% TrainPop@id, colnames(A_matrix) %in% TrainPop@id]
  
  
  model <- mmer(cbind(Yield, Protein)~ YearEff,
                random = ~vs(GID,Gu=A_matrix_fin, Gtc=unsm(2)),
                rcov = ~vs(units, Gtc=unsm(2)),
                data = y_train, verbose = TRUE)
  
  model$sigma$`u:GID`
  #extract coeffictients 
  
  coefficients <- obtain_coefficients(model, desired)
  
  # ---- Advance Year ----
  # Stage 7: Release variety (no code here, placeholder for additional logic if needed)
  
  # Stage 6
  #stat
  #EYT <- selectInd(AYT, nEYT) previous funtion, just take into account one trait
  
  EYT <- pheno_selection_index(AYT, coefficients, nEYT) # phenotypic selection based on index selection
  
  #EYT <- AYT[AYT@id %in% EYT_id,]
  
  EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
  
  
  
  # Stage 5
  #AYT <- selectInd(PYT, nAYT)
  
  AYT <- pheno_selection_index(PYT, coefficients, nAYT)
  #AYT <- PYT[PYT@id %in% AYT_id,]
  AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
  
  # Stage 4
  #output$accSel[year] <- cor(PYT@gv, PYT@pheno)
  
  PYT <- select_index_within(HDRW, coefficients,famMax, pattern = "phenotype") # aqui ira una funcion para seleccion withinfam en multitrait
  
  
  PYT <- pheno_selection_index(PYT,coefficients, nPYT)
  #PYT <- PYT[PYT@id %in% PYT_id]
  PYT <- setPheno(PYT, varE = VarE, reps = repPYT)
  
  # Stage 3
  HDRW <- setPheno(DH, varE = VarE, reps = repHDRW)
  
  #output$accSel_True[year] <- cor(HDRW@gv, HDRW@pheno)
  
  # Stage 2
  DH <- makeDH(F1, nDH)
  
  # ---- Load Global Parameters ----
  #source(file = "GlobalParameters/GlobalParameters_Control.R")
  
  # Stage 1
  F1 <- randCross(Parents, nCrosses)
  
  # ---- Store Training Population ----
  PYT@fixEff <- as.integer(rep(year, nInd(PYT)))
  AYT@fixEff <- as.integer(rep(year, nInd(AYT)))
  EYT@fixEff <- as.integer(rep(year, nInd(EYT)))
  
  if (year > nBurnin) {
    cat("  Maintaining training population \n")
    TrainPop <- c(TrainPop[-c(1:c(PYT, EYT, AYT)@nInd)], PYT, EYT, AYT)
  }
  
  # ---------------- Pedigree & Inbreeding metrics ----------------
  if (year == Start_Pedigree) {
    cat("  Start recording pedigree and set founders \n")
    
    # Initialize master_ped with a year column
    master_ped <- data.frame(
      id     = character(),
      mother = character(),
      father = character(),
      year   = integer(),
      stringsAsFactors = FALSE
    )
    
    # Create the initial pedigree with founders
    master_ped <- getMaster_ped(
      master_ped,
      Parents = TrainPop,
      year    = year,
      cycles  = 0:10
    )
  }
  
  if (year > Start_Pedigree) {
    cat("  Update pedigree and record inbreeding \n")
    
    # Update the existing master pedigree dataframe
    master_ped <- updateMaster_ped(
      master_ped,
      Parents = TrainPop,
      year    = year,
      cycles  = 0:10
    )
    
    # Get the pedigree inbreeding value
    #Fvalue_A <- getFvalue_A(master_ped, Parents)
    #output$Fvalue_A[year] <- Fvalue_A
  }
  
  if (year == Start_Pedigree) {
    # Store the allelic frequencies of the reference population
    #markers_ref <- pullSnpGeno(TrainPop)
    #pvec_ref    <- colMeans(markers_ref) / 2
  }
  
  if (year > Start_Pedigree) {
    # Get the population inbreeding from vanRaden GRM using reference freqs
    #Fvalue_G <- getFvalue_G(Parents, pvec_ref)
    #output$Fvalue_G[year] <- Fvalue_G
  }
  
  # ---- Store DH Metrics ----
  output$meanG_yield[year]     <- meanG(DH)[1]
  output$meanG_protein[year]     <- meanG(DH)[2]
  
  output$varG_yield[year]      <- varG(DH)[1]
  output$varG_protein[year]      <- varG(DH)[2]
  output$genicVarG_yield[year] <- genicVarG(DH)[1]
  output$genicVarG_protein[year] <- genicVarG(DH)[2]
}

write.table(
  output[output$rep == REP & output$year > nBurnin, ],
  file      = pheno_csv,
  sep       = ",",
  col.names = FALSE,
  row.names = FALSE,
  append    = TRUE
)




cat("PS_RM phase completed and results saved in '", pheno_csv, "'.\n")
