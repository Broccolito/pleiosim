#' Simulate Genomes
#'
#' This function simulates genotype data for a cohort based on the specified
#' effect allele frequency (EAF). It generates genotypes with probabilities
#' corresponding to homozygous (2), heterozygous (1), and wild-type (0) states.
#'
#' @param n_participant Integer specifying the number of participants.
#' @param eaf Numeric vector representing the effect allele frequencies for each SNP.
#'
#' @return A matrix of simulated genotypes with individuals as rows and SNPs as columns.
#'
#' @export
simulate_genomes = function(n_participant, eaf){
  f_hom = eaf^2
  f_het = 2*eaf*(1-eaf)
  f_wt = (1-eaf)^2
  prob_matrix = cbind(f_hom, f_het, f_wt)
  genomes = apply(prob_matrix, 1, function(p){
    sample(c(2, 1, 0), size = n_participant, replace = TRUE, prob = p)
  })
  return(genomes)
}

#' Simulate Cohort Genomes
#'
#' This function simulates genotype data for multiple cohorts based on cohort makeup
#' and effect allele frequency (EAF) matrices. It handles overlapping cohorts and
#' assigns participant identifiers.
#'
#' @param cohort_makeup_matrix A matrix specifying the number of participants in each cohort.
#' @param eaf_matrix A matrix representing effect allele frequencies for each variant across cohorts.
#'
#' @return A list of matrices where each matrix represents the simulated genotypes for a cohort.
#'
#' @export
simulate_cohort_genomes = function(cohort_makeup_matrix, eaf_matrix){
  n_phenotype = sum(!grepl(pattern = "_", colnames(cohort_makeup_matrix)))
  cohort_names = colnames(cohort_makeup_matrix)[!grepl(pattern = "_", colnames(cohort_makeup_matrix))]
  variant_names = rownames(eaf_matrix)
  intersection_names = unlist(lapply(2:n_phenotype, function(k){
    combn(cohort_names, k, FUN = function(x) paste(x, collapse = "_"))
  }))

  simulated_genomes = lapply(1:(2^n_phenotype-1), function(i){
    simulate_genomes(n = cohort_makeup_matrix[1,i], eaf = eaf_matrix[,i])
  })
  names(simulated_genomes) = c(cohort_names, intersection_names)

  for(i in 1:length(simulated_genomes)){
    if(!identical(simulated_genomes[[i]], numeric(0))){
      participant_names = paste0(names(simulated_genomes)[i], "_participant",
                                 1:nrow(simulated_genomes[[i]]))
      rownames(simulated_genomes[[i]]) = participant_names
      colnames(simulated_genomes[[i]]) = variant_names
    }else{
      simulated_genomes[[i]] = matrix(ncol = length(variant_names), nrow = 0)
      colnames(simulated_genomes[[i]]) = variant_names
    }
  }

  cohort_genomes = lapply(as.list(cohort_names), function(x){
    return(simulated_genomes[[x]])
  })
  names(cohort_genomes) = cohort_names

  for(intersection_name in intersection_names){
    for(cohort_name in unlist(strsplit(intersection_name, split = "_"))){
      cohort_genomes[[cohort_name]] = rbind(cohort_genomes[[cohort_name]],
                                            simulated_genomes[[intersection_name]])
    }
  }

  return(cohort_genomes)
}

#' Simulate Non-Heritable Phenotype Matrix
#'
#' This function simulates non-heritable phenotypes based on the specified
#' non-heritable correlation matrix, sample sizes, and heritability values. It accounts
#' for environmental variance to generate realistic non-genetic influences on phenotypes.
#'
#' @param heritable_phenotype_matrix A list of matrices representing the heritable component of phenotypes.
#' @param nonheritable_correlation_matrix A matrix specifying the correlation between non-heritable phenotypes.
#' @param defacto_sample_size_matrix A matrix indicating the sample size for each cohort.
#' @param heritability A numeric vector representing the heritability for each phenotype.
#'
#' @return A list of matrices where each matrix represents non-heritable phenotypes for a cohort.
#'
#' @export
simulate_nonheritable_phenotype_matrix = function(heritable_phenotype_matrix,
                                                  nonheritable_correlation_matrix,
                                                  defacto_sample_size_matrix,
                                                  heritability){
  n_phenotype = nrow(nonheritable_correlation_matrix)
  phenotype_names = paste0("phenotype", 1:n_phenotype)

  eigen_decomp = eigen(nonheritable_correlation_matrix)
  eigen_vectors = eigen_decomp$vectors
  eigen_values = diag(pmax(eigen_decomp$values, 0))

  simulate_environmental_noise = function(n, eigen_vectors, eigen_values){
    uncorrelated_noise = matrix(rnorm(n * n_phenotype), nrow = n, ncol = n_phenotype)
    correlated_noise = uncorrelated_noise %*% eigen_vectors %*% sqrt(eigen_values)
    return(correlated_noise)
  }

  nonheritable_phenotype_matrix = lapply(defacto_sample_size_matrix, function(n){
    simulate_environmental_noise(n, eigen_vectors, eigen_values)
  })
  names(nonheritable_phenotype_matrix) = phenotype_names

  environmental_variance = mapply(function(genetic_component, h2){
    genetic_variance = var(genetic_component)
    env_variance = (genetic_variance * (1 - h2)) / h2
    return(env_variance)
  }, heritable_phenotype_matrix, heritability)

  for (i in 1:n_phenotype){
    nonheritable_phenotype_matrix[[i]][, i] = scale(nonheritable_phenotype_matrix[[i]][, i]) * sqrt(environmental_variance[i])
  }
  return(nonheritable_phenotype_matrix)
}


