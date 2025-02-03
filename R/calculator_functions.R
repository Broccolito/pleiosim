#' Calculate De Facto EAF Matrix
#'
#' This function calculates the de facto Effect Allele Frequency (EAF) matrix
#' for a given list of cohort genomes.
#'
#' @param cohort_genomes A list of matrices where each matrix represents the genomes of a cohort.
#' Each matrix has individuals as rows and SNPs as columns.
#'
#' @return A matrix where each column represents the EAF for each cohort.
#'
#' @export
calculate_defacto_eaf_matrix = function(cohort_genomes){
  defacto_eaf_matrix = lapply(cohort_genomes, function(x){
    colSums(x)/(dim(x)[1]*2)
  })
  defacto_eaf_matrix = do.call(cbind, defacto_eaf_matrix)
  return(defacto_eaf_matrix)
}

#' Calculate De Facto Sample Size Matrix
#'
#' This function computes the sample size matrix for a given list of cohort genomes.
#'
#' @param cohort_genomes A list of matrices representing the genomes of different cohorts.
#'
#' @return A matrix with the number of participants for each cohort.
#'
#' @export
calculate_defacto_sample_size_matrix = function(cohort_genomes){
  defacto_sample_size_matrix = lapply(cohort_genomes, function(x){
    nrow(x)
  })
  defacto_sample_size_matrix = do.call(cbind, defacto_sample_size_matrix)
  rownames(defacto_sample_size_matrix) = "n_participant"
  return(defacto_sample_size_matrix)
}

#' Calculate Genetic Correlation Matrix
#'
#' This function calculates the genetic correlation matrix using the effect size matrix
#' and the de facto EAF matrix.
#'
#' @param efs_matrix_template A matrix of effect sizes with phenotypes as rows and SNPs as columns.
#' @param defacto_eaf_matrix A matrix representing the effect allele frequencies across cohorts.
#'
#' @return A genetic correlation matrix.
#'
#' @export
calculate_genetic_correlation = function(efs_matrix_template, defacto_eaf_matrix){
  n_phenotype = nrow(efs_matrix_template)
  genotype_variance = 2*defacto_eaf_matrix*(1-defacto_eaf_matrix)
  defacto_genetic_covariance_matrix = matrix(0, nrow = n_phenotype, ncol = n_phenotype)
  for(i in 1:n_phenotype){
    for(j in 1:n_phenotype){
      defacto_genetic_covariance_matrix[i, j] = sum(efs_matrix_template[i, ] * efs_matrix_template[j, ] * genotype_variance[, i])
    }
  }

  defacto_genetic_correlation_matrix = matrix(0, nrow = n_phenotype, ncol = n_phenotype)
  for(i in 1:n_phenotype){
    for(j in 1:n_phenotype){
      defacto_genetic_correlation_matrix[i, j] = defacto_genetic_covariance_matrix[i, j] /
        sqrt(defacto_genetic_covariance_matrix[i, i] * defacto_genetic_covariance_matrix[j, j])
    }
  }
  return(defacto_genetic_correlation_matrix)
}

#' Calculate Heritable Phenotype Matrix
#'
#' This function computes the heritable phenotype matrix based on cohort genomes
#' and the effect size (EFS) matrix.
#'
#' @param cohort_genomes A list of matrices where each matrix represents the genomes of a cohort.
#' @param efs_matrix A matrix of effect sizes with phenotypes as rows and SNPs as columns.
#'
#' @return A list of matrices, each representing heritable phenotypes for a cohort.
#'
#' @export
calculate_heritable_phenotype_matrix = function(cohort_genomes, efs_matrix){
  n_phenotype = nrow(efs_matrix)
  phenotype_names = paste0("phenotype", 1:n_phenotype)
  heritable_phenotype_matrix = lapply(1:n_phenotype, function(i){
    heritable_phenotype = cohort_genomes[[i]] %*% efs_matrix[i,]
    colnames(heritable_phenotype) = phenotype_names[i]
    return(heritable_phenotype)
  })
  names(heritable_phenotype_matrix) = phenotype_names
  return(heritable_phenotype_matrix)
}

#' Calculate Combined Phenotype Matrix
#'
#' This function combines heritable and non-heritable phenotype matrices
#' to generate the overall phenotype matrix for each cohort.
#'
#' @param heritable_phenotype_matrix A list of matrices representing heritable phenotypes.
#' @param nonheritable_phenotype_matrix A list of matrices representing non-heritable phenotypes.
#'
#' @return A list of matrices with combined phenotypes for each cohort.
#'
#' @export
calculate_phenotype_matrix = function(heritable_phenotype_matrix,
                                      nonheritable_phenotype_matrix){
  n_phenotype = length(heritable_phenotype_matrix)
  phenotype_names = paste0("phenotype", 1:n_phenotype)
  cohort_phenotypes = lapply(1:n_phenotype, function(i){
    combined_phenotype_matrix = heritable_phenotype_matrix[[i]] + nonheritable_phenotype_matrix[[i]][, i]
    return(combined_phenotype_matrix)
  })
  names(cohort_phenotypes) = phenotype_names
  return(cohort_phenotypes)
}

#' Calculate Summary Statistics Matrix
#'
#' This function calculates the summary statistics (beta, standard error, p-value, EAF)
#' for each SNP based on the cohort genomes and phenotypes.
#'
#' @param cohort_genomes A list of matrices representing the genomes of different cohorts.
#' @param cohort_phenotypes A list of matrices representing the phenotypes of different cohorts.
#'
#' @return A list of matrices containing summary statistics for each phenotype.
#'
#' @export
calculate_summary_stats_matrix = function(cohort_genomes,
                                          cohort_phenotypes){
  n_phenotype = length(cohort_phenotypes)
  phenotype_names = paste0("phenotype", 1:n_phenotype)
  summary_stats_matrix = lapply(1:n_phenotype, function(i){
    genotypes = cohort_genomes[[i]]
    phenotypes = cohort_phenotypes[[i]]

    summary_stats = t(apply(genotypes, 2, function(snp){
      l = summary(lm(phenotypes[,1] ~ snp))$coefficients
      beta = l[2,1]
      se = l[2,2]
      p_value = l[2,4]
      eaf = sum(snp)/(length(snp)*2)
      summary_stat = matrix(c(beta, se, p_value, eaf), ncol = 4)
      return(summary_stat)
    }))
    colnames(summary_stats) = c("beta", "se", "p_value", "eaf")
    return(summary_stats)
  })
  names(summary_stats_matrix) = phenotype_names
  return(summary_stats_matrix)
}

#' Calculate Estimated Heritable Correlation Matrix
#'
#' This function estimates the heritable correlation matrix based on summary statistics.
#'
#' @param summary_stats_matrix A list of matrices containing summary statistics for each phenotype.
#'
#' @return A matrix representing the estimated heritable correlation between phenotypes.
#'
#' @export
calculate_estimated_heritable_correlation_matrix = function(summary_stats_matrix){
  n_phenotype = length(summary_stats_matrix)

  defacto_eaf_matrix = lapply(summary_stats_matrix, function(x){x[,4]})
  defacto_eaf_matrix = do.call(cbind, defacto_eaf_matrix)

  defacto_efs_matrix = lapply(summary_stats_matrix, function(x){x[,1]})
  defacto_efs_matrix = do.call(rbind, defacto_efs_matrix)

  genotype_variance = 2 * defacto_eaf_matrix * (1 - defacto_eaf_matrix)

  genetic_covariance = matrix(0, nrow = n_phenotype, ncol = n_phenotype)
  for (i in 1:n_phenotype){
    for (j in 1:n_phenotype){
      genetic_covariance[i, j] = sum(defacto_efs_matrix[i, ] * defacto_efs_matrix[j, ] * genotype_variance[, i])
    }
  }

  estimated_heritable_correlation_matrix = matrix(0, nrow = n_phenotype, ncol = n_phenotype)
  for(i in 1:n_phenotype){
    for(j in 1:n_phenotype){
      estimated_heritable_correlation_matrix[i, j] = genetic_covariance[i, j] /
        sqrt(genetic_covariance[i, i] * genetic_covariance[j, j])
    }
  }

  return(estimated_heritable_correlation_matrix)
}
