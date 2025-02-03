#' Create Heritable Correlation Matrix
#'
#' This function generates a heritable correlation matrix with uniform correlation values
#' for off-diagonal elements and 1 for diagonal elements.
#'
#' @param n_phenotype Integer specifying the number of phenotypes. Default is 3.
#' @param uniform_heritable_correlation Numeric value representing the uniform correlation between phenotypes. Default is 0.4.
#'
#' @return A square matrix representing the heritable correlation between phenotypes.
#'
#' @export
create_heritable_correlation_matrix = function(n_phenotype = 3,
                                               uniform_heritable_correlation = 0.4){
  heritable_correlation_matrix = matrix(uniform_heritable_correlation, nrow = n_phenotype, ncol = n_phenotype)
  diag(heritable_correlation_matrix) = 1
  return(heritable_correlation_matrix)
}

#' Create Non-Heritable Correlation Matrix
#'
#' This function generates a non-heritable correlation matrix with uniform correlation values
#' for off-diagonal elements and 1 for diagonal elements.
#'
#' @param n_phenotype Integer specifying the number of phenotypes.
#' @param uniform_nonheritable_correlation Numeric value representing the uniform correlation between non-heritable phenotypes. Default is 0.2.
#'
#' @return A square matrix representing the non-heritable correlation between phenotypes.
#'
#' @export
create_nonheritable_correlation_matrix = function(n_phenotype,
                                                  uniform_nonheritable_correlation = 0.2){
  nonheritable_correlation_matrix = matrix(uniform_nonheritable_correlation, nrow = n_phenotype, ncol = n_phenotype)
  diag(nonheritable_correlation_matrix) = 1
  return(nonheritable_correlation_matrix)
}

#' Create Cohort Makeup Matrix
#'
#' This function generates a cohort makeup matrix defining the number of participants
#' for each phenotype and their intersections.
#'
#' @param n_phenotype Integer specifying the number of phenotypes. Default is 3.
#' @param n_participant Numeric value or vector indicating the number of participants. Default is 1000.
#' @param unique_cohorts Logical indicating whether participants are unique to each cohort. Default is TRUE.
#'
#' @return A matrix representing the cohort makeup with the number of participants for each cohort and their intersections.
#'
#' @export
create_cohort_makeup_matrix = function(n_phenotype = 3,
                                       n_participant = 1000,
                                       unique_cohorts = TRUE){

  cohort_names = paste0("cohort", 1:n_phenotype)
  cohort_makeup_matrix = matrix(0, ncol = 2^n_phenotype-1, nrow = 1)

  intersection_names = unlist(lapply(2:n_phenotype, function(k){
    combn(cohort_names, k, FUN = function(x) paste(x, collapse = "_"))
  }))

  colnames(cohort_makeup_matrix) = c(cohort_names, intersection_names)
  rownames(cohort_makeup_matrix) = c("n_participant")

  if(length(n_participant) == 1){
    if(unique_cohorts){
      cohort_makeup_matrix[1,1:n_phenotype] = n_participant
    }else{
      cohort_makeup_matrix[1,2^n_phenotype-1] = n_participant
    }
  }else{
    if(length(n_participant) != n_phenotype){
      stop("Wrong participant number dimensions specified...\n")
    }else{
      cohort_makeup_matrix[1,1:n_phenotype] = n_participant
    }
  }
  return(cohort_makeup_matrix)
}

#' Create Effect Allele Frequency (EAF) Matrix
#'
#' This function generates an EAF matrix for different variants across multiple cohorts.
#'
#' @param n_phenotype Integer specifying the number of phenotypes. Default is 3.
#' @param uniform_eaf Numeric value representing the uniform effect allele frequency. Default is 0.4.
#' @param n_variant_pleiotropic Number of pleiotropic variants. Default is 10.
#' @param n_variant_nonpleiotropic Vector specifying the number of non-pleiotropic variants per phenotype. Default is c(10, 10, 10).
#' @param n_variant_null Number of null variants. Default is 960.
#'
#' @return A matrix representing effect allele frequencies for variants across cohorts.
#'
#' @export
create_eaf_matrix = function(n_phenotype = 3,
                             uniform_eaf = 0.4,
                             n_variant_pleiotropic = 10,
                             n_variant_nonpleiotropic = c(10, 10, 10),
                             n_variant_null = 960){

  n_variant = n_variant_pleiotropic + sum(n_variant_nonpleiotropic) + n_variant_null

  variant_names = c(paste0("pleio_variant", 1:n_variant_pleiotropic),
                    paste0("nonpleio_variant", 1:sum(n_variant_nonpleiotropic)),
                    paste0("null_variant", 1:n_variant_null))

  cohort_names = paste0("cohort", 1:n_phenotype)
  intersection_names = unlist(lapply(2:n_phenotype, function(k){
    combn(cohort_names, k, FUN = function(x) paste(x, collapse = "_"))
  }))

  eaf_matrix = matrix(uniform_eaf, ncol = 2^n_phenotype-1, nrow = n_variant)
  colnames(eaf_matrix) = c(cohort_names, intersection_names)
  rownames(eaf_matrix) = variant_names

  return(eaf_matrix)
}

#' Create Effect Size (EFS) Matrix Template
#'
#' This function generates a template for the effect size (EFS) matrix,
#' specifying the effect of each variant on different phenotypes.
#'
#' @param n_phenotype Integer specifying the number of phenotypes. Default is 3.
#' @param n_variant_pleiotropic Number of pleiotropic variants. Default is 10.
#' @param n_variant_nonpleiotropic Vector specifying the number of non-pleiotropic variants per phenotype. Default is c(10, 10, 10).
#' @param n_variant_null Number of null variants. Default is 960.
#' @param crosstrait_heterogeneity Logical indicating if cross-trait heterogeneity should be applied. Default is TRUE.
#' @param withintrait_heterogeneity Logical indicating if within-trait heterogeneity should be applied. Default is TRUE.
#' @param random_crosstrait_heterogeneity Logical indicating if cross-trait heterogeneity should be random. Default is FALSE.
#' @param random_withintrait_heterogeneity Logical indicating if within-trait heterogeneity should be random. Default is FALSE.
#'
#' @return A matrix template representing the effect sizes for each variant and phenotype.
#'
#' @export
create_efs_matrix_template = function(n_phenotype = 3,
                                      n_variant_pleiotropic = 10,
                                      n_variant_nonpleiotropic = c(10, 10, 10),
                                      n_variant_null = 960,
                                      crosstrait_heterogeneity = TRUE,
                                      withintrait_heterogeneity = TRUE,
                                      random_crosstrait_heterogeneity = FALSE,
                                      random_withintrait_heterogeneity = FALSE){

  phenotype_names = paste0("phenotype", 1:n_phenotype)
  variant_names = c(paste0("pleio_variant", 1:n_variant_pleiotropic),
                    paste0("nonpleio_variant", 1:sum(n_variant_nonpleiotropic)),
                    paste0("null_variant", 1:n_variant_null))

  n_variant = n_variant_pleiotropic + sum(n_variant_nonpleiotropic) + n_variant_null

  efs_matrix_template_pleiotropic = matrix(1, nrow = n_phenotype, ncol = n_variant_pleiotropic)

  efs_matrix_template_nonpleiotropic = vector()
  for(v in 1:length(n_variant_nonpleiotropic)){
    efs_v = matrix(0, nrow = n_phenotype, ncol = n_variant_nonpleiotropic[v])
    efs_v[v,] = 1
    efs_matrix_template_nonpleiotropic = cbind(
      efs_matrix_template_nonpleiotropic,
      efs_v
    )
  }

  efs_matrix_template_null = matrix(0, nrow = n_phenotype, ncol = n_variant_null)

  efs_matrix_template = cbind(
    efs_matrix_template_pleiotropic,
    efs_matrix_template_nonpleiotropic,
    efs_matrix_template_null
  )

  rownames(efs_matrix_template) = phenotype_names
  colnames(efs_matrix_template) = variant_names

  if(random_crosstrait_heterogeneity){
    crosstrait_heterogeneity = FALSE
    sampled_traits = sample(1:n_phenotype, floor(n_phenotype/2))
    efs_matrix_template[sampled_traits,] = -1 * efs_matrix_template[sampled_traits,]
  }

  if(random_withintrait_heterogeneity){
    withintrait_heterogeneity = FALSE
    sampled_variants = sample(1:n_variant, floor(n_variant/2))
    efs_matrix_template[,sampled_variants] = -1 * efs_matrix_template[,sampled_variants]
  }

  if(crosstrait_heterogeneity){
    selected_traits = seq(1, n_phenotype, 2)
    efs_matrix_template[selected_traits,] = -1 * efs_matrix_template[selected_traits,]
  }

  if(withintrait_heterogeneity){
    selected_variants = seq(1, n_variant, 2)
    efs_matrix_template[,selected_variants] = -1 * efs_matrix_template[,selected_variants]
  }

  return(efs_matrix_template)
}

#' Create Effect Size (EFS) Matrix
#'
#' This function generates the effect size (EFS) matrix by optimizing parameters
#' to minimize the loss function based on the genetic correlation. The optimization
#' aims to align the simulated genetic correlation with the provided heritable correlation matrix.
#'
#' @param efs_matrix_template A matrix template representing initial effect sizes.
#' @param defacto_eaf_matrix A matrix representing the effect allele frequencies.
#' @param heritable_correlation_matrix A matrix specifying the target heritable correlation
#'   between phenotypes, used to guide the optimization process.
#'
#' @return A matrix representing optimized effect sizes for each variant and phenotype.
#'
#' @export
create_efs_matrix = function(efs_matrix_template,
                             defacto_eaf_matrix,
                             heritable_correlation_matrix){

  n_variant_pleiotropic = sum(apply(efs_matrix_template, 2, prod) != 0)
  n_phenotype = nrow(efs_matrix_template)

  loss_function = function(params){
    adjusted_matrix = efs_matrix_template
    adjusted_matrix[, 1:n_variant_pleiotropic] = adjusted_matrix[, 1:n_variant_pleiotropic] * params

    adjusted_genetic_correlation_matrix = calculate_genetic_correlation(adjusted_matrix,
                                                                        defacto_eaf_matrix = defacto_eaf_matrix)

    loss = sum((adjusted_genetic_correlation_matrix - heritable_correlation_matrix)^2)
    return(loss)
  }

  initial_params = rep(1, n_phenotype)
  optimized_multiplier = optim(par = initial_params, fn = loss_function,
                               method = "L-BFGS-B", lower = -100, upper = 100)$par

  efs_matrix = efs_matrix_template
  efs_matrix[, 1:n_variant_pleiotropic] = efs_matrix[, 1:n_variant_pleiotropic] * optimized_multiplier

  return(efs_matrix)
}

