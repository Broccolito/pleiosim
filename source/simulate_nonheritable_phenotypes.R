simulate_nonheritable_phenotypes = function(heritable_phenotype_matrix, heritability,
                                            n_phenotype, n_sample,
                                            nonheritable_correlation_matrix,
                                            variable_names){
  
  nonheritable_phenotype_variance = apply(heritable_phenotype_matrix, MARGIN = 2, var) * ((1/heritability)-1)
  nonheritable_phenotype_std = sqrt(nonheritable_phenotype_variance)
  
  sigma = (nonheritable_phenotype_std %*% t(nonheritable_phenotype_std)) * nonheritable_correlation_matrix
  
  eigen_object = eigen(sigma, symmetric = TRUE)
  eigen_vector = eigen_object$vectors
  eigen_value = eigen_object$values
  
  nonheritable_phenotype_matrix = matrix(rnorm(n_phenotype*n_sample), n_sample)
  nonheritable_phenotype_matrix = eigen_vector %*% diag(sqrt(pmax(eigen_value,0)), 
                                                        n_phenotype) %*% t(nonheritable_phenotype_matrix)
  nonheritable_phenotype_matrix = t(nonheritable_phenotype_matrix)
  
  colnames(nonheritable_phenotype_matrix) = variable_names[["phenotype_names"]]
  rownames(nonheritable_phenotype_matrix) = variable_names[["subject_names"]]
  
  return(nonheritable_phenotype_matrix)
  
}