#' @export
check_input = function(
    n_sample,
    n_variant_pleiotropic,
    n_variant_nonpleiotropic,
    n_variant_null,
    eaf,
    n_phenotype,
    heritable_correlation_matrix,
    nonheritable_correlation_matrix,
    heritability,
    crosstrait_heterogeneity,
    withintrait_heterogeneity,
    random_crosstrait_heterogeneity,
    random_withintrait_heterogeneity
){

  if(length(eaf) != n_variant_pleiotropic + sum(n_variant_nonpleiotropic) + n_variant_null){
    stop("The number of variants specified is not consistant with the effect size matrix...")
  }

  if(!all(n_phenotype == dim(heritable_correlation_matrix))){
    stop("The number of phenotypes specified is not consistant with the heritable correlation matrix...")
  }

  if(!all(n_phenotype == dim(nonheritable_correlation_matrix))){
    stop("The number of phenotypes specified is not consistant with the non-heritable correlation matrix...")
  }

  if(sum(!(t(heritable_correlation_matrix) == heritable_correlation_matrix))!=0){
    stop("The heritable correlation matrix is not a symmetrical matrix...")
  }

  if(sum(!(t(nonheritable_correlation_matrix) == nonheritable_correlation_matrix))!=0){
    stop("The heritable correlation matrix is not a symmetrical matrix...")
  }

  if(n_phenotype != length(heritability)){
    stop("The number of phenotypes specified is not consistant with the heritabilities specified...")
  }

  if(max(heritable_correlation_matrix[!diag(n_phenotype)]) >= 0.8){
    warning("High genotypic correlations may compromise the validity of some pleiotropic methods...")
  }

  if(max(nonheritable_correlation_matrix[!diag(n_phenotype)]) >= 0.8){
    warning("High environmental correlations may compromise the validity of some pleiotropic methods...")
  }

  if(n_phenotype < 3 && crosstrait_heterogeneity){
    warning("The introduction of crosstrait heterogeneity cannot be randomized...")
  }

  if(length(eaf) < 3 && withintrait_heterogeneity){
    warning("The introduction of withintrait heterogeneity cannot be randomized; consider simulating more variants...")
  }

}


