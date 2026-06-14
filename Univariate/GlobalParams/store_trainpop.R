# store train population 

if (year == startTP){
  cat("  Start collecting training population \n")
  TrainPop = c(PYT, EYT, AYT)
}

if (year > startTP & year < nBurnin+1){
  cat("  Collecting training population \n")
  TrainPop = c(TrainPop,
               PYT, EYT, AYT)
} #si entiendo bien guardaria 2 poblaciones hasta donde tengo hecho

if (year > nBurnin){
  cat("  Maintaining training population \n")
  TrainPop = c(TrainPop[-c(1:c(PYT, EYT, AYT)@nInd)], #con esto borras de la poblacion inicial el numero de individuos que vas a meter en tu pob nueva
               PYT, EYT, AYT)
}

