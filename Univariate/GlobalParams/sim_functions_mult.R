
# modificacion de esta funcion para que trabaje en multitrait y univariate
extract_meff <- function(mmatrix, blups){
  # 1. CENTRADO DE LA MATRIZ DE MARCADORES
  # Es fundamental que sea la misma transformación que en setGEAVs
  mmatrix <- scale(mmatrix, center = TRUE, scale = FALSE)
  
  n_ind <- nrow(mmatrix)
  
  #add little penalization to invert the marker matrix
  mm_inv <- solve(tcrossprod(mmatrix)+ diag(1e-6, n_ind))
  
  #proyection from the blups to the marker effects
  marker_effects <- t(mmatrix)%*%mm_inv%*%as.matrix(blups)
  
  #i expect a m(number of markers) x n(number of traits), that show the marker effect of each marker for the 2 traits separately
  return(marker_effects)
}



setGEAVs <- function(geno_obj, addEffs_list, trait_name=NULL) {
  # 1. pull new individuals' SNP matrix (rows = individuals, cols = SNPs)
  new_markers <- pullSnpGeno(geno_obj)
  
  # --- OPERACIÓN DE CENTRADO POR MEDIA ---
  # 'scale' con center=TRUE y scale=FALSE resta la media de cada columna
  # Esto hace que la media de cada SNP sea exactamente 0
  new_markers <- scale(new_markers, center = TRUE, scale = FALSE)
  
  #Exrtact the marker names, used in the model
  ref_markers <- rownames(addEffs_list)
  
  # 2. filter to the SNPs present in the reference marker set
  common_snps <- intersect(colnames(new_markers), ref_markers)
  if (length(common_snps) == 0) {
    stop("No overlapping SNP columns between new_markers and ref_markers.")
  }
  filt_new <- new_markers[, common_snps, drop = FALSE]
  
  # 3. extract the effect vector for the chosen trait, default all the traits
  if (!is.null(trait_name)) {
    if (!trait_name %in% names(addEffs_list)) {
      stop("'", trait_name, "' not found in addEffs_list.")
    }
    eff_vec <- addEffs_list[[trait_name]]
    if (length(eff_vec) != length(common_snps)) {
      stop("Length of effect vector (", length(eff_vec), 
           ") does not match number of filtered SNPs (", length(common_snps), ").")
    }
    # 4. compute GEBVs as matrix multiplication
    GEBVs <- filt_new %*% eff_vec
    
    # 5. name the output
    colnames(GEBVs) <- trait_name
    rownames(GEBVs) <- NULL
  }
    
  GEBVs <- filt_new %*% addEffs_list
  return(GEBVs)
}




#' Remove columns with zero variance
#'
#' @param mat A numeric matrix or data.frame
#' @return A matrix (or data.frame) with only those columns whose variance is non-zero
#' @export
filter_zero_variance <- function(mat) {
  # compute variance of each column
  variances <- apply(mat, 2, var, na.rm = TRUE)
  # keep only columns with non-zero variance
  filtered <- mat[, variances != 0, drop = FALSE]
  return(filtered)
}



#' Compute VanRaden additive genomic relationship matrix (GRM)
#'
#' @param Markers A numeric matrix or data.frame of genotypes (individuals × SNPs),
#'        encoded as allele counts (e.g., 0, 1, 2).
#' @return A numeric matrix: the VanRaden additive GRM.
#' @importFrom stats scale
#' @export
getVR_GRM <- function(Markers) {
  # ensure a matrix
  M <- as.matrix(Markers)
  
  # 1. estimate allele frequencies p_j = mean genotype / 2
  pvec <- colMeans(M) / 2
  
  # 2. center each column: subtract 2 * p_j
  #    (scale=FALSE so only centering)
  Z <- scale(M, center = TRUE, scale = FALSE)
  
  # 3. compute denominator: 2 * sum[p_j * (1 - p_j)]
  denom <- 2 * sum(pvec * (1 - pvec))
  
  # 4. build GRM
  Ga <- tcrossprod(Z) / denom
  
  # 5. clean up temporary object
  rm(Z); gc()
  
  return(Ga)
}

#function to compute the index value of the individuals, and then order them to select the best individuals to pass to the next generation
# ' 
#' @param Population A population , who must have already the set pheno function loaded
#' @param WeigthVector Vector of economic weigths to create the Smith Hazel Index
#' @param NumberIndiv Number of individuals that are passed to the next generation
#' 
#weight <- c(1, 2)
pheno_selection_index <- function(Population, Weigthvector, Numberindiv){
  #extract the pheno data and create the data frame
  
  pheno_data <- data.frame(
    GID = Population@id,
    Yield = pheno(Population)[1],
    Protein_content = pheno(Population)[2])
  
  df_temp <- sweep(pheno_data[, c("Yield", "Protein_content")], 2, Weigthvector, "*")
  pheno_data$Index <- rowSums(df_temp)
  
  #order the individuals based on their index selection
  pheno_data <- pheno_data[order(pheno_data$Index, decreasing = TRUE),]
  
  pheno_data_select <- pheno_data[1:Numberindiv,]
  
  finalpop <- Population[Population@id %in% pheno_data_select$GID]
  return(finalpop)
}


#AYT_id <- pheno_selection_index(Population = PYT, Weigthverctor = weight, Numberindiv =nAYT)

#AYT <- PYT[PYT@id %in% AYT_id,]


select_index_within <- function(Population, coefficients, nFamMax, pattern = "phenotype") {
  
  # 1. Extraer la matriz de datos según el patrón ("phenotype" o "genotype")
  # Si es genotype, usamos @ebv; si no, usamos pheno()
  data_matrix <- if (pattern == "phenotype") pheno(Population) else Population@ebv
  
  #data_matrix <- pheno(Population)
  df_fam <- data.frame(data_matrix)
  # 2. Calcular el Índice (Suma ponderada)
  # Multiplicamos cada columna por su peso y sumamos por filas
  #index_values <- rowSums(t(t(data_matrix) * WeightVector))
  rownames(df_fam) <- Population@id
  # 3. Crear población temporal para la selección
  # Copiamos la población original para no perder sus datos reales
  temp_pop <- Population
  w <- as.numeric(coefficients)
  
  df_weigthed <- sweep(df_fam, 2, w, "*")
  
  df_fam$Index <- rowSums(df_weigthed)
  # Reemplazamos sus fenotipos por el índice calculado (en formato matriz)
  # Esto obliga a AlphaSimR a seleccionar basándose en nuestro Índice
  temp_pop@pheno <- as.matrix(df_fam$Index)
  
  # 4. Ejecutar la selección dentro de familias
  # Selecciona los 'nFamMax' mejores individuos de cada familia según el Índice
  selected_temp <- selectWithinFam(temp_pop, nFamMax)
  
  # 5. Retornar los individuos originales
  # Usamos los IDs seleccionados para filtrar la población inicial
  final_pop <- Population[Population@id %in% selected_temp@id,]
  
  return(final_pop)
}





extract_heritability <- function(model, trait_value){
  #trait value must be 1, 2 ....
  param <- trait_value
  
  h2 <- model$theta[[1]][param,param] / (model$theta[[1]][param,param] + model$theta[[2]][param, param])
  return(h2)

}

#economic_weight <- c(1,10)
geno_selection_index <- function(population ,economicWeight, number_selected) {
  #the population must have already compute it its EBVs
  
  population_df <- as.data.frame(population@ebv)
  # 3. Calcular el Índice usando sweep (multiplica cada fila por el peso)
  # El '2' indica que la operación es sobre las columnas
  df_weighted <- sweep(population_df, 2, economicWeight, "*")
  
  # 4. Sumar los valores pesados para obtener el Índice final
  population_df$Index <- rowSums(df_weighted)
  
  population_df$GID <- rownames(population_df)
  
  population_df <- population_df[order(population_df$Index, decreasing = TRUE),]
  
  #Select the new members of the next population 
  geno_data_select <- population_df[1:number_selected,]
  
  finalpop <- population[population@id %in% geno_data_select$GID,]
  
  return(finalpop)
}

#newParents <- geno_selection_index(PYT, economic_weight, 10)




#extract coefficients
obtain_coefficients <- function(model, desired_index ) {
  gamma <-  model$sigma$`u:GID`
  
  H <- diag(sqrt(diag(gamma)))
  
  P <- cov(y_train[,c("Yield", "Protein")])
  
  Q <- H %*% solve(gamma) %*% P %*% solve(gamma) %*% H
  
  d <- c(1,0.5)
  i <- 1
  l <- i / sqrt(t(desired_index)%*%Q%*%desired_index)
  x <- l*desired_index
  #x
  b <- solve(gamma)%*%H %*% x
  
  return(as.vector(b))
}

#coefficients <- obtain_coefficients(model, desired)



get_marginal_variance <- function(Z,G,R) {
  # we need many different things to compute the marginal variance of the blups 

  
  
  V <- Z%*%G%*%t(Z) + R #variance of the Y in MVN model
  V_inv <- solve(V)
  P <- V_inv - V_inv %*%X%*%solve(t(X)%*%V_inv%*%X)%*%t(X)%*%V_inv 
  
  PEV <- G - G%*%t(Z)%*%P%*%Z%*%G
  
  var_uhat <- G - PEV #marginal variance of the blups
  #para muchos datos es mas eficiuante quiza hacer el ec de henderrson
  return(var_uhat)
}


compute_coefficients_blups <- function(var_uhat, dataset ) {
  #first transform the dataset we used in sommer `mmer` and transform it in longformat
  pheno_long <- dataset %>% pivot_longer(cols = c("Yield", "Protein"), names_to = "Trait", values_to = "value") %>%
    arrange(desc(Trait)) %>% mutate(across(c(Trait, GID), as.factor))
  #have name of the traits 
  traits <- unique(pheno_long$Trait)
  
  n <- n_distinct(pheno_long$GID) #number of disting genotypes
  t <- n_distinct(pheno_long$Trait) #number of traits
  ids <- unique(pheno_long$GID) #ids of genotypes
  
  names <- as.vector(outer(traits, unique(dataset$GID), paste, sep = ":"))
  #dimnames(G) <- 
    
  dimnames(var_uhat) <- list(names, names)
  #var_uhat[1:4, 1:4]
  
  nt <- nrow(var_uhat)
  idx <- split(1:nt, f = rep(1:t, times = n))  #contains the position for each traits variances 
  B <- matrix(0, ncol = t, nrow = t)
  for (i in 1:t) {
    for (j in i:t) {
      idxi <- idx[[i]]
      idxj <- idx[[j]]
      L <- var_uhat[idxi, idxj]
      B[i, j] <- B[j, i] <- mean(diag(L)) - mean(L)
    }
  }
  dimnames(B) <- list(traits, traits)
  P <- Gamma <- B # P equals Gamma since we are working under BLUP
  Gamma
  
  H <- diag(sqrt(diag(sigma_g)))
  dimnames(H) <- list(traits, traits)
  H
  Q <- H%*%solve(Gamma)%*%H
  desired <- c(1, 0)
  intensity <- 1
  x.opt <- desired * intensity/as.numeric(sqrt(crossprod(desired, Q%*%desired))) # apply dimensionality
  
  b <- (solve(Gamma) %*% H %*% x.opt)/intensity #substitute response to solve for coefficients
  b_desired <- b/sqrt(sum(b^2))
  
  
  return(b_desired)
}


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
