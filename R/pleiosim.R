#' Pleiotropic Simulation (pleiosim)
#'
#' This function performs a comprehensive simulation of pleiotropic genetic data,
#' including genotype generation, phenotype calculation, and summary statistics.
#'
#' @param n_phenotype Number of phenotypes. Default is 3.
#' @param heritable_correlation Uniform heritable correlation value. Default is 0.4.
#' @param nonheritable_correlation Uniform non-heritable correlation value. Default is 0.2.
#' @param n_participant Number of participants. Default is 10000.
#' @param eaf Effect allele frequency. Default is 0.4.
#' @param n_variant_pleiotropic Number of pleiotropic variants. Default is 10.
#' @param n_variant_nonpleiotropic Number of non-pleiotropic variants per phenotype. Default is c(10, 10, 10).
#' @param n_variant_null Number of null variants. Default is 960.
#' @param heritability Vector of heritability values for each phenotype. Default is c(0.1, 0.2, 0.3).
#' @param unique_cohorts Logical indicating whether participants are unique across cohorts. Default is TRUE.
#' @param crosstrait_heterogeneity Logical for cross-trait heterogeneity. Default is TRUE.
#' @param withintrait_heterogeneity Logical for within-trait heterogeneity. Default is TRUE.
#' @param random_crosstrait_heterogeneity Logical for random cross-trait heterogeneity. Default is FALSE.
#' @param random_withintrait_heterogeneity Logical for random within-trait heterogeneity. Default is FALSE.
#' @param customized_heritable_correlation_matrix Optional custom heritable correlation matrix.
#' @param customized_nonheritable_correlation_matrix Optional custom non-heritable correlation matrix.
#' @param customized_cohort_makeup_matrix Optional custom cohort makeup matrix. Ensure correct row/column names.
#' @param customized_eaf_matrix Optional custom EAF matrix. Ensure correct structure with appropriate row/column names.
#' @param customized_efs_matrix_template Optional custom EFS matrix template. Follow given data structure.
#' @param customized_efs_matrix Optional custom EFS matrix. Ensure consistency with matrix structure.
#'
#' @return A "pleio" class object containing simulated genetic, phenotype, and summary statistics data.
#'
#' @details
#' For `customized_cohort_makeup_matrix`, `customized_eaf_matrix`, `customized_efs_matrix_template`,
#' and `customized_efs_matrix`, users can provide custom matrices by specifying the data directly.
#' Ensure the provided matrices follow the given data structure, including proper row and column names.
#'
#' @examples
#' pleio_object = pleiosim(
#'   n_phenotype = 3,
#'   heritable_correlation = 0.4,
#'   nonheritable_correlation = 0.2,
#'   n_participant = 10000,
#'   eaf = 0.4,
#'   n_variant_pleiotropic = 10,
#'   n_variant_nonpleiotropic = c(10, 10, 10),
#'   n_variant_null = 960,
#'   heritability = c(0.1, 0.2, 0.3),
#'   unique_cohorts = TRUE,
#'   crosstrait_heterogeneity = TRUE,
#'   withintrait_heterogeneity = TRUE,
#'   random_crosstrait_heterogeneity = FALSE,
#'   random_withintrait_heterogeneity = FALSE,
#'   customized_heritable_correlation_matrix = NULL,
#'   customized_nonheritable_correlation_matrix = NULL,
#'   customized_cohort_makeup_matrix = NULL,
#'   customized_eaf_matrix = NULL,
#'   customized_efs_matrix_template = NULL,
#'   customized_efs_matrix = NULL
#' )
#'
#' summary(pleio_object)
#'
#' @export
pleiosim = function(n_phenotype = 3,
                    heritable_correlation = 0.4,
                    nonheritable_correlation = 0.2,
                    n_participant = 10000,
                    eaf = 0.4,
                    n_variant_pleiotropic = 10,
                    n_variant_nonpleiotropic = c(10, 10, 10),
                    n_variant_null = 960,
                    heritability = c(0.1, 0.2, 0.3),
                    unique_cohorts = TRUE,
                    crosstrait_heterogeneity = TRUE,
                    withintrait_heterogeneity = TRUE,
                    random_crosstrait_heterogeneity = FALSE,
                    random_withintrait_heterogeneity = FALSE,
                    customized_heritable_correlation_matrix = NULL,
                    customized_nonheritable_correlation_matrix = NULL,
                    customized_cohort_makeup_matrix = NULL,
                    customized_eaf_matrix = NULL,
                    customized_efs_matrix_template = NULL,
                    customized_efs_matrix = NULL){

  if(is.null(customized_heritable_correlation_matrix)){
    heritable_correlation_matrix = create_heritable_correlation_matrix(n_phenotype = n_phenotype,
                                                                       uniform_heritable_correlation = heritable_correlation)
  }else{
    heritable_correlation_matrix = customized_heritable_correlation_matrix
  }

  if(is.null(customized_nonheritable_correlation_matrix)){
    nonheritable_correlation_matrix = create_nonheritable_correlation_matrix(n_phenotype = n_phenotype,
                                                                             uniform_nonheritable_correlation = nonheritable_correlation)
  }else{
    nonheritable_correlation_matrix = customized_nonheritable_correlation_matrix
  }

  if(is.null(customized_cohort_makeup_matrix)){
    cohort_makeup_matrix = create_cohort_makeup_matrix(n_phenotype = n_phenotype,
                                                       n_participant = n_participant,
                                                       unique_cohorts = unique_cohorts)
  }else{
    cohort_makeup_matrix = customized_cohort_makeup_matrix
  }

  if(is.null(customized_eaf_matrix)){
    eaf_matrix = create_eaf_matrix(n_phenotype = n_phenotype,
                                   uniform_eaf = eaf,
                                   n_variant_pleiotropic = n_variant_pleiotropic,
                                   n_variant_nonpleiotropic = n_variant_nonpleiotropic,
                                   n_variant_null = n_variant_null)
  }else{
    eaf_matrix = customized_eaf_matrix
  }

  if(is.null(customized_efs_matrix_template)){
    efs_matrix_template = create_efs_matrix_template(n_phenotype = n_phenotype,
                                                     n_variant_pleiotropic = n_variant_pleiotropic,
                                                     n_variant_nonpleiotropic = n_variant_nonpleiotropic,
                                                     n_variant_null = n_variant_null,
                                                     crosstrait_heterogeneity = crosstrait_heterogeneity,
                                                     withintrait_heterogeneity = withintrait_heterogeneity,
                                                     random_crosstrait_heterogeneity = random_crosstrait_heterogeneity,
                                                     random_withintrait_heterogeneity = random_withintrait_heterogeneity)
  }else{
    efs_matrix_template = customized_efs_matrix_template
  }

  cohort_genomes = simulate_cohort_genomes(cohort_makeup_matrix = cohort_makeup_matrix,
                                           eaf_matrix = eaf_matrix)

  defacto_eaf_matrix = calculate_defacto_eaf_matrix(cohort_genomes = cohort_genomes)
  defacto_sample_size_matrix = calculate_defacto_sample_size_matrix(cohort_genomes = cohort_genomes)

  if(is.null(customized_efs_matrix)){
    efs_matrix = create_efs_matrix(efs_matrix_template = efs_matrix_template,
                                   defacto_eaf_matrix = defacto_eaf_matrix,
                                   heritable_correlation_matrix = heritable_correlation_matrix)
  }else{
    efs_matrix = customized_efs_matrix
  }

  heritable_phenotype_matrix = calculate_heritable_phenotype_matrix(cohort_genomes = cohort_genomes,
                                                                    efs_matrix = efs_matrix)

  nonheritable_phenotype_matrix = simulate_nonheritable_phenotype_matrix(heritable_phenotype_matrix = heritable_phenotype_matrix,
                                                                         nonheritable_correlation_matrix = nonheritable_correlation_matrix,
                                                                         defacto_sample_size_matrix = defacto_sample_size_matrix,
                                                                         heritability = heritability)

  cohort_phenotypes = calculate_phenotype_matrix(heritable_phenotype_matrix = heritable_phenotype_matrix,
                                                 nonheritable_phenotype_matrix = nonheritable_phenotype_matrix)

  summary_stats_matrix = calculate_summary_stats_matrix(cohort_genomes = cohort_genomes,
                                                        cohort_phenotypes = cohort_phenotypes)

  estimated_heritable_correlation_matrix = calculate_estimated_heritable_correlation_matrix(summary_stats_matrix)

  pleio = new("pleio",
              n_phenotype = n_phenotype,
              n_variant_pleiotropic = n_variant_pleiotropic,
              n_variant_nonpleiotropic = n_variant_nonpleiotropic,
              n_variant_null = n_variant_null,
              heritability = heritability,
              unique_cohorts = unique_cohorts,
              crosstrait_heterogeneity = crosstrait_heterogeneity,
              withintrait_heterogeneity = withintrait_heterogeneity,
              random_crosstrait_heterogeneity = random_crosstrait_heterogeneity,
              random_withintrait_heterogeneity = random_withintrait_heterogeneity,
              heritable_correlation_matrix = heritable_correlation_matrix,
              nonheritable_correlation_matrix = nonheritable_correlation_matrix,
              cohort_makeup_matrix = cohort_makeup_matrix,
              eaf_matrix = eaf_matrix,
              efs_matrix_template = efs_matrix_template,
              efs_matrix = efs_matrix,
              heritable_phenotype_matrix = heritable_phenotype_matrix,
              nonheritable_phenotype_matrix = nonheritable_phenotype_matrix,
              defacto_eaf_matrix = defacto_eaf_matrix,
              defacto_sample_size_matrix = defacto_sample_size_matrix,
              estimated_heritable_correlation_matrix = estimated_heritable_correlation_matrix,
              cohort_genomes = cohort_genomes,
              cohort_phenotypes = cohort_phenotypes,
              summary_stats_matrix = summary_stats_matrix
  )

  return(pleio)
}

