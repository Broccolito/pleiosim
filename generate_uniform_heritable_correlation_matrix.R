generate_uniform_heritable_correlation_matrix = function(heritable_correlation, 
                                                         n_phenotype){
  
  heritable_correlation_matrix = matrix(heritable_correlation, nrow = n_phenotype, ncol = n_phenotype)
  diag(heritable_correlation_matrix) = 1
  
  return(heritable_correlation_matrix)
}