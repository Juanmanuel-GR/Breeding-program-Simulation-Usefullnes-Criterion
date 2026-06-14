# Create founders

#generate initial haplotypes
founderPop <- runMacs(nInd = nParents,
                      nChr = nChr,
                      segSites = nSnp + nQtl,
                      inbred = TRUE,
                      species = "WHEAT")

SP <- SimParam$new(founderPop)

#Add SNP chip 
SP$restrSegSites(nQtl, nSnp) #esto le dice al progrm q elija a la azar el numero d qtl q van 
#a marcar donde esta la info q da el fenotipo y tb q marque como imp los el n d snp asociados a esos qtl
#añadir param de snp 
if (nSnp > 0) {
  SP$addSnpChip(nSnp)
}

SP$addTraitAG(nQtlPerChr = nQtl,
              mean = initmeanG,
              var = initVarG,
              varEnv = initVarEnv,
              varGxE = initVarGE)

SP$traitNames <- "Yield"

# collect pedigree
SP$setTrackPed(TRUE)

#create parents population 
Parents <- newPop(founderPop)

# Add phenotype reflecting evalatuation in EYT
#The reps parameter is for convenient representation of replicated data. It is intended to represent
#replicated yield trials in plant breeding programs. In this case, varE is set to the plot error and reps is
#set to the number of plots per entry. The resulting phenotype represents the entry-means.
Parents <- setPheno(Parents, varE = VarE, reps = repEYT)
# la varianza se reducira pq lo has probado en muchos sitios, piensas q inicialmente la pob viene de eyt, los primeros padres
rm(founderPop)