#' @export
heterogeneity_transform = function(efs_matrix_template,
                                   n_phenotype,
                                   crosstrait_heterogeneity,
                                   withintrait_heterogeneity,
                                   random_crosstrait_heterogeneity,
                                   random_withintrait_heterogeneity,
                                   variable_names){

  n_variant = dim(efs_matrix_template)[2]

  crosstrait_heterogeneity_matrix = rep(1, n_phenotype)
  if(crosstrait_heterogeneity){
    crosstrait_heterogeneity_matrix = (-1)^(1:n_phenotype)
  }

  withintrait_heterogeneity_matrix = rep(1, n_variant)
  if(withintrait_heterogeneity){
    withintrait_heterogeneity_matrix = (-1)^(1:n_variant)
  }

  if(random_crosstrait_heterogeneity){
    crosstrait_heterogeneity_matrix = sample(
      crosstrait_heterogeneity_matrix,
      n_phenotype,
      replace = TRUE
    )
  }

  if(random_withintrait_heterogeneity){
    withintrait_heterogeneity_matrix = sample(
      withintrait_heterogeneity_matrix,
      n_variant,
      replace = TRUE
    )
  }

  efs_matrix_template = apply(efs_matrix_template, MARGIN = 2, function(x){
    x * crosstrait_heterogeneity_matrix
  })

  efs_matrix_template = t(apply(efs_matrix_template, MARGIN = 1, function(x){
    x * withintrait_heterogeneity_matrix
  }))

  colnames(efs_matrix_template) = variable_names[["variant_names"]]
  rownames(efs_matrix_template) = variable_names[["phenotype_names"]]

  return(efs_matrix_template)
}

# NOT RUN
# heterogeneity_transform(efs_matrix_template = efs_matrix,
#                         crosstrait_heterogeneity = TRUE,
#                         withintrait_heterogeneity = TRUE,
#                         random_crosstrait_heterogeneity = FALSE,
#                         random_withintrait_heterogeneity = TRUE)

