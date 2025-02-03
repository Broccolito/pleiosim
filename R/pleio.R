#' Pleiotropic Simulation Class (pleio)
#'
#' Defines the S4 class 'pleio' for storing pleiotropic simulation results,
#' including genotype data, phenotype data, and summary statistics.
#'
#' @slot n_phenotype Number of phenotypes simulated.
#' @slot n_variant_pleiotropic Number of pleiotropic variants.
#' @slot n_variant_nonpleiotropic Number of non-pleiotropic variants per phenotype.
#' @slot n_variant_null Number of null variants.
#' @slot heritability Vector of heritability values for each phenotype.
#' @slot unique_cohorts Logical indicating if participants are unique across cohorts.
#' @slot crosstrait_heterogeneity Logical indicating cross-trait heterogeneity.
#' @slot withintrait_heterogeneity Logical indicating within-trait heterogeneity.
#' @slot random_crosstrait_heterogeneity Logical for random cross-trait heterogeneity.
#' @slot random_withintrait_heterogeneity Logical for random within-trait heterogeneity.
#' @slot heritable_correlation_matrix Matrix of heritable correlation values.
#' @slot nonheritable_correlation_matrix Matrix of non-heritable correlation values.
#' @slot cohort_makeup_matrix Matrix defining the cohort structure.
#' @slot eaf_matrix Matrix of effect allele frequencies.
#' @slot efs_matrix_template Template matrix for effect size simulation.
#' @slot efs_matrix Matrix of effect sizes.
#' @slot heritable_phenotype_matrix List of heritable phenotype data.
#' @slot nonheritable_phenotype_matrix List of non-heritable phenotype data.
#' @slot defacto_eaf_matrix Matrix of de facto effect allele frequencies.
#' @slot defacto_sample_size_matrix Matrix of sample sizes for each cohort.
#' @slot estimated_heritable_correlation_matrix Estimated matrix of heritable correlations.
#' @slot cohort_genomes List of genotype data for each cohort.
#' @slot cohort_phenotypes List of phenotype data for each cohort.
#' @slot summary_stats_matrix List of summary statistics for each phenotype.
#'
#' @export
setClass(
  Class = "pleio",
  slots = list(
    n_phenotype = "numeric",
    n_variant_pleiotropic = "numeric",
    n_variant_nonpleiotropic = "numeric",
    n_variant_null = "numeric",
    heritability = "numeric",
    unique_cohorts = "logical",
    crosstrait_heterogeneity = "logical",
    withintrait_heterogeneity = "logical",
    random_crosstrait_heterogeneity = "logical",
    random_withintrait_heterogeneity = "logical",
    heritable_correlation_matrix = "matrix",
    nonheritable_correlation_matrix = "matrix",
    cohort_makeup_matrix = "matrix",
    eaf_matrix = "matrix",
    efs_matrix_template = "matrix",
    efs_matrix = "matrix",
    heritable_phenotype_matrix = "list",
    nonheritable_phenotype_matrix = "list",
    defacto_eaf_matrix = "matrix",
    defacto_sample_size_matrix = "matrix",
    estimated_heritable_correlation_matrix = "matrix",
    cohort_genomes = "list",
    cohort_phenotypes = "list",
    summary_stats_matrix = "list"
  ),
  prototype = list(
    n_phenotype = NULL,
    n_variant_pleiotropic = NULL,
    n_variant_nonpleiotropic = NULL,
    n_variant_null = NULL,
    heritability = NULL,
    unique_cohorts = NULL,
    crosstrait_heterogeneity = NULL,
    withintrait_heterogeneity = NULL,
    random_crosstrait_heterogeneity = NULL,
    random_withintrait_heterogeneity = NULL,
    heritable_correlation_matrix = NULL,
    nonheritable_correlation_matrix = NULL,
    cohort_makeup_matrix = NULL,
    eaf_matrix = NULL,
    efs_matrix_template = NULL,
    efs_matrix = NULL,
    heritable_phenotype_matrix = NULL,
    nonheritable_phenotype_matrix = NULL,
    defacto_eaf_matrix = NULL,
    defacto_sample_size_matrix = NULL,
    estimated_heritable_correlation_matrix = NULL,
    cohort_genomes = NULL,
    cohort_phenotypes = NULL,
    summary_stats_matrix = NULL
  )
)

#' Display Method for pleio Class
#'
#' Custom display method for objects of class 'pleio' to summarize key simulation parameters.
#'
#' @param object An object of class 'pleio'.
#'
#' @export
setMethod(
  f = "show",
  signature = "pleio",
  definition = function(object){
    cat("[pleio object]\n")
    cat("A simulation of: \n")
    cat(paste0("- ", object@n_phenotype, " phenotypes\n"))
    cat(paste0("- ", (object@n_variant_pleiotropic + sum(object@n_variant_nonpleiotropic) + object@n_variant_null), " SNPs\n"))
    cat(paste0(" - ", object@n_variant_pleiotropic, " pleiotropic SNPs\n"))
    cat(paste0(" - ", sum(object@n_variant_nonpleiotropic), " non-pleiotropic SNPs\n"))
    cat(paste0(" - ", object@n_variant_null, " null SNPs\n"))
    cat(paste0("- Total of ", sum(object@defacto_sample_size_matrix), " participants\n"))
  }
)
