# this scripts aims to create the initial population for filling the breeding pipeline
#set inital yield trials with unique individuals, so on each stage we dont find the same individuals bc the come from dif crosses

for (cohort in 1:7) {
  cat (" FillPipeline Stage:", cohort, "of 7\n")
  if (cohort < 7) {
    F1 <- randCross(Parents, nCrosses = nCrosses)
  }
  if (cohort < 6) {
    #stage 2
    DH <- makeDH(F1, nDH)
  }
  if (cohort < 5) {
   #Stage 3
    HDRW <- setPheno(DH, varE = VarE, reps = repHDRW) 
  }
  if (cohort < 4) {
    #stage 4
    PYT <- selectWithinFam(HDRW, famMax)
    PYT <- selectInd(PYT, nPYT)
    PYT <- setPheno(PYT, varE = VarE, reps = repPYT)
  }
  if (cohort < 3) {
    #stage 5
    AYT <- selectInd(PYT, nAYT)
    AYT <- setPheno(AYT, varE = VarE, reps = repAYT)
  }
  if (cohort<2) {
    #stage 6
    EYT <- selectInd(AYT, nEYT)
    EYT <- setPheno(EYT, varE = VarE, reps = repEYT)
  }
  if (cohort<1) {
    #stage 7
    
  }
}
