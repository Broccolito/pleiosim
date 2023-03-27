generate_efs_matrix_template = function(n_variant_pleiotropic,
                                        n_variant_nonpleiotropic,
                                        n_variant_null,
                                        n_phenotype,
                                        variable_names){
  
  efs_matrix_template_pleiotropic = matrix(1, nrow = n_phenotype, ncol = n_variant_pleiotropic)
  
  efs_matrix_template_nonpleiotropic = vector()
  for(i in 1:length(n_variant_nonpleiotropic)){
    efs_i = matrix(0, nrow = n_phenotype, ncol = n_variant_nonpleiotropic[i])
    efs_i[i,] = 1
    efs_matrix_template_nonpleiotropic = cbind(
      efs_matrix_template_nonpleiotropic,
      efs_i
    )
  }
  
  efs_matrix_template_null = matrix(0, nrow = n_phenotype, ncol = n_variant_null)
  
  efs_matrix_template = cbind(
    efs_matrix_template_pleiotropic,
    efs_matrix_template_nonpleiotropic,
    efs_matrix_template_null
  )
  
  colnames(efs_matrix_template) = variable_names[["variant_names"]]
  rownames(efs_matrix_template) = variable_names[["phenotype_names"]]
  
  return(efs_matrix_template)
}

## NOT RUN
# generate_efs_matrix_template(n_variant_pleiotropic = 10,
#                              n_variant_nonpleiotropic = c(10, 10, 10),
#                              n_variant_null = 960,
#                              n_phenotype = 3)