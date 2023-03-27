generate_variable_names = function(eaf, n_phenotype, n_sample){
  
  n_variant = length(eaf)
  variant_names = paste0("snp", 1:n_variant)
  phenotype_names = paste0("pheno", 1:n_phenotype)
  subject_names = paste0("subject",1:n_sample)
  
  variable_names = list(
    variant_names = variant_names,
    phenotype_names = phenotype_names,
    subject_names = subject_names
  )
  
  return(variable_names)
  
}