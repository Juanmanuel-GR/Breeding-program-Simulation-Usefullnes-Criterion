# Advance year
#este es ek script q servira para hacer el burning y seleccionar individ pa pasar a la siguiente poblacion

# Advance breeding program by 1 year
# Works backwards through pipeline to avoid copying data

# Stage 7
# Release variety

#stage 6
EYT <- selectInd(AYT, nEYT)
EYT <- setPheno(EYT, varE = VarE, reps = repEYT)

#stage 5
AYT <- selectInd(PYT, nAYT)
AYT <- setPheno(AYT, varE = varE, reps = repAYT)

#stage 4
#aqui viene lo del output para guardarlo en una tabla, x ahora no lo hago
PYT <- selectWithinFam(HDRW, famMax) #ver cuando indiv hay aqui
PYT <- selectInd(PYT, nInd = nPYT)
PYT <- setPheno(PYT, varE = varE, reps = repPYT)

#stage 3
HDRW <- setPheno(DH, varE = VarE, reps = repHDRW)

#Stage2 
DH <- makeDH(F1, nDH)

#Stage 1
F1 <- randCross(Parents, nCrosses = n_Crosses) #number progeny = nCrosses, default 1

PYT@nInd
