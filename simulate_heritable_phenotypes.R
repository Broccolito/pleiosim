simulate_heritable_phenotypes = function(efs_matrix, genotype_matrix,
                                         variable_names){
  
  heritable_phenotype_matrix = t(efs_matrix %*% t(genotype_matrix))
  
  colnames(heritable_phenotype_matrix) = variable_names[["phenotype_names"]]
  rownames(heritable_phenotype_matrix) = variable_names[["subject_names"]]
  
  return(heritable_phenotype_matrix)
  
}