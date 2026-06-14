# Future_PS_RM.R

# ---- Clean Environment and Load Packages ----
rm(list = ls())
library(AlphaSimR)

#setwd("C:/Users/jmgrb/OneDrive/MASTER/Practicas/Tutoriales/AlphasimR/Simulation/future_scenario")
#getwd()
#load(file = "../BurnInGen/burning_states/BurnIn_REP_1.RData")
#Retrieve Parameters -----
REP <- Sys.getenv("SLURM_ARRAY_TASK_ID")
SCENARIO <- "PS_RM"
#REP <- 1
source(file = "GlobalParams/Globa_params.R")
#source(file = "GlobalParams/Globa_params.R")
#--------- Helper Functions ----- Idont need them rn, but...
#poner desps

burnin_dir <- "../../Multivariate/Burnin_Gen_MVN/burnin_states"
results_dir <- paste0("results_", SCENARIO)
output_csv_ps <- file.path(results_dir, paste0("PS_RM_Results_", REP, ".csv"))


# Create Results Directory
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

#Initialize CSV file
if(file.exists(output_csv_ps)) file.remove(output_csv_ps)

#create the header for the csv file
header <- c("rep",
            "year",
            "meanG_yield", "meanG_protein","varG_yield" ,"varG_protein", "genicVarG_yield", "genicVarG_protein")

write.table(t(header),#entiendo q hay q poner la traspuesta pq se pondria en las filas y lo quieres como columna
            file = output_csv_ps,
            sep = ",",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)




  
#---Future PS Phase 
cat("Processing the future scenario of REP", REP, "SCENARIO", SCENARIO, "\n")
  
#load the burnIn state data
burnin_file <- file.path(burnin_dir, paste0("BurnIn_REP_", REP, ".RData"))
load(burnin_file)
  
# ---- Extend 'output' for PS phase ----
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
  
  
# Future PS Simulation Loop ----
for (year in (nBurnin+1):(nBurnin+ nfuture)) {
  cat( " Future year:", year, "\n")
    
  #current_row <- year - nBurnin
  # ---- Update Parents ----
  # Replace 10 oldest inbred parents with 10 new parents from EYT stage
  Parents <- c(Parents[11:nParents], EYT)
    
  # ---- Advance Year ----
  # Stage 7: Release variety (no code here, placeholder for additional logic if needed)
    
  # Stage 6
  EYT <- selectInd(AYT, nEYT, trait = "Yield", use = "pheno")
  EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
    
  # Stage 5
  AYT <- selectInd(PYT, nAYT, trait = "Yield", use = "pheno")
  AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
    
  # Stage 4
  #output$accSel[year] <- cor(PYT@gv, PYT@pheno)
    
  PYT <- selectWithinFam(HDRW, famMax, trait = "Yield", use = "pheno")
  PYT <- selectInd(PYT, nPYT, trait = "Yield", use = "pheno")
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
  file      = output_csv_ps,
  sep       = ",",
  col.names = FALSE,
  row.names = FALSE,
  append    = TRUE
)
  

cat("PS_RM phase completed and results saved in '", output_csv_ps, "'.\n")
