#' @export
generate_uniform_nonheritable_correlation_matrix = function(nonheritable_correlation,
                                                          n_phenotype){

  nonheritable_correlation_matrix = matrix(nonheritable_correlation, nrow = n_phenotype, ncol = n_phenotype)
  diag(nonheritable_correlation_matrix) = 1

  return(nonheritable_correlation_matrix)
}
