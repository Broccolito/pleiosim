simulate_phenotypes = function(heritable_phenotype_matrix, 
                               nonheritable_phenotype_matrix,
                               variable_names){
  
  phenotype_matrix = heritable_phenotype_matrix + nonheritable_phenotype_matrix
  
  colnames(phenotype_matrix) = variable_names[["phenotype_names"]]
  rownames(phenotype_matrix) = variable_names[["subject_names"]]
  
  return(phenotype_matrix)
  
}