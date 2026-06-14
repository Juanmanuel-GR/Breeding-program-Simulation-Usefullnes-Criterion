#simulation where i simulate a population with multitrait

library(AlphaSimR)
#library(tidyverse)
#library(sommer)
library(AGHmatrix)
rm(list = ls())



#load the global Parameters 
source(file = "Global_Params/global_param_mult.R")
source(file="Global_Params/sim_functions_mult.R")

output_dir <- "burnin_states"
output_csv <- "burnIn_Results_mv.csv"

# Create Results Directory
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

#Initialize CSV file
if(file.exists(output_csv)) file.remove(output_csv)

#create the header for the csv file
header <- c(
  "rep", "year", "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein"
)
write.table(t(header),#entiendo q hay q poner la traspuesta pq se pondria en las filas y lo quieres como columna
            file = output_csv,
            sep = ",",
            col.names = FALSE,
            row.names = FALSE)


for (REP in 1:40) { #40 reps 
  cat("Working on Burn-in REP:", REP, "\n")
  
  # Initialize the simulation
  source(file = "Global_Params/create_parents_mult.R")
  source(file = "Global_Params/fillPipeline.R")
  
  #create a data frame to track key params
  output <- data.frame(
    rep         = REP,
    year        = 1:nBurnin, # un vector q cada linea le pone un year
    meanG_yield       = numeric(nBurnin),
    meanG_protein = numeric(nBurnin),
    varG_yield = numeric(nBurnin),
    varG_protein        = numeric(nBurnin), #funciona pq todos tienen el mismo numero d filas
    genicVarG_yield = numeric(nBurnin),
    genicVarG_protein = numeric(nBurnin),
    weighted_heritability = numeric(nBurnin)
  )
 
  
  # Run the burn-in for the defined number of years
  for (year in 1:nBurnin) {
    cat(" Burn-in year", year, "\n")
    
    #Replace 10 parents
    Parents <- c(Parents[11:nParents], EYT)
    # -------------------- Stage 7 --------------------
    # Release variety
    
    # -------------------- Stage 6 --------------------
    EYT <- selectInd(AYT, nEYT, trait = 1)
    EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
    
    # -------------------- Stage 5 --------------------
    AYT <- selectInd(PYT, nAYT, trait = 1)
    AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
    
    # -------------------- Stage 4 --------------------
    #output$accSel[year] <- cor(HDRW@gv, HDRW@pheno)
    PYT <- selectWithinFam(HDRW, famMax)
    PYT <- selectInd(PYT, nPYT, trait = 1)
    PYT <- setPheno(PYT, varE = VarE, reps = repPYT)
    
    # -------------------- Stage 3 --------------------
    HDRW <- setPheno(DH, varE = VarE, reps = repHDRW)
    
    # -------------------- Stage 2 --------------------
    DH <- makeDH(F1, nDH)
    
    # -------------------- Stage 1 --------------------
    F1 <- randCross(Parents, nCrosses)
    
    # ---------------- Training population ----------------
    PYT@fixEff <- as.integer(rep(year, nInd(PYT)))
    AYT@fixEff <- as.integer(rep(year, nInd(AYT)))
    EYT@fixEff <- as.integer(rep(year, nInd(EYT)))
    
    
    #train Pop
    source("Global_Params/store_trainpop.R")
    
    # ---------------- Pedigree & Inbreeding metrics ----------------
    # peedigree needed for Phenotypic selection Multivariate
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
        cycles  = 0:10 #10 generations for fully homozygous
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
    
    # ---------------- Genetic metrics ----------------
    output$meanG_yield[year]     <- meanG(DH)[1]
    output$meanG_protein[year] <- meanG(DH)[2]
    output$varG_yield[year]      <- varG(DH)[1]
    output$varG_protein[year]      <- varG(DH)[2,2]
    output$genicVarG_yield[year] <- genicVarG(DH)[1]
    output$genicVarG_protein[year] <- genicVarG(DH)[2]
    
  }
  # ---- Save the state for the current replicate ----
  save.file <- paste0("BurnIn_REP_", REP, ".RData")
  save.image(file = file.path(output_dir, save.file))
  
  # ---- Append results to CSV file ----
  write.table(
    output,
    file      = output_csv,
    sep       = ",",
    col.names = FALSE,
    row.names = FALSE,
    append    = TRUE
  )
  
}
cat(
  "Burn-in phase completed and results saved in '",
  output_dir, "' and '", output_csv, "'.\n",
  sep = "")






