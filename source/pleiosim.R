pleiosim = function(
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
  
  ## Add error check
  
  cat("\n\n----------------------------------------\n")
  cat(paste0("Simulating ", n_sample, " subjects...\n"))
  cat(paste0("Simulating ", n_variant_pleiotropic, " pleiotropic variants...\n"))
  cat(paste0("Simulating ", sum(n_variant_nonpleiotropic), " non-pleiotropic functional variants...\n"))
  cat(paste0("Simulating ", n_variant_null, " null variants...\n"))
  cat(paste0("Simulating ", n_phenotype, " phenotypes...\n"))
  cat("----------------------------------------\n\n")
  cat(paste0("Simulation starts...\n"))
  
  cat(paste0("Generating variable names for genotypes, phenotypes, and subjects by default...\n"))
  variable_names = generate_variable_names(
    eaf = eaf, 
    n_phenotype = n_phenotype,
    n_sample = n_sample
  )
  
  cat(paste0("Generating genotype matrix...\n"))
  genotype_matrix = simulate_genotypes(
    n_sample = n_sample, 
    eaf = eaf, 
    variable_names = variable_names
  )
  
  cat(paste0("Generating effect size matrix template...\n"))
  efs_matrix_template = generate_efs_matrix_template(
    n_variant_pleiotropic = n_variant_pleiotropic,
    n_variant_nonpleiotropic = n_variant_nonpleiotropic,
    n_variant_null = n_variant_null,
    n_phenotype = n_phenotype,
    variable_names = variable_names
  )
  
  cat(paste0("Transforming effect size matrix based on heterogeneity...\n"))
  efs_matrix_template = heterogeneity_transform(
    efs_matrix_template = efs_matrix_template, 
    n_phenotype = n_phenotype,
    crosstrait_heterogeneity = TRUE,
    withintrait_heterogeneity = TRUE,
    random_crosstrait_heterogeneity = FALSE,
    random_withintrait_heterogeneity = TRUE,
    variable_names = variable_names
  )
  
  cat(paste0("Transforming effect size matrix based on pre-defined correlations...\n"))
  efs_matrix = correlation_transform(
    efs_matrix_template = efs_matrix_template, 
    genotype_matrix = genotype_matrix,
    n_phenotype = n_phenotype,
    n_variant_pleiotropic = n_variant_pleiotropic,
    heritable_correlation_matrix = heritable_correlation_matrix,
    variable_names = variable_names
  )
  
  cat(paste0("Generating heritable portions of the phenotype...\n"))
  heritable_phenotype_matrix = simulate_heritable_phenotypes(
    efs_matrix = efs_matrix, 
    genotype_matrix = genotype_matrix, 
    variable_names = variable_names
  )
  defacto_heritable_correlation_matrix = cor(heritable_phenotype_matrix)
  
  cat(paste0("Generating non-heritable portions of the phenotype...\n"))
  nonheritable_phenotype_matrix = simulate_nonheritable_phenotypes(
    heritable_phenotype_matrix = heritable_phenotype_matrix,
    heritability = heritability,
    n_phenotype = n_phenotype,
    n_sample = n_sample,
    nonheritable_correlation_matrix = nonheritable_correlation_matrix,
    variable_names = variable_names
  )
  defacto_nonheritable_correlation_matrix = cor(nonheritable_phenotype_matrix)
  
  cat(paste0("Combining the heritable and non-heritable portions of the phenotype ...\n"))
  phenotype_matrix = simulate_phenotypes(
    heritable_phenotype_matrix = heritable_phenotype_matrix,
    nonheritable_phenotype_matrix = nonheritable_phenotype_matrix,
    variable_names = variable_names
  )
  
  cat(paste0("Running Genome-Wide Association analysis...\n"))
  gwas_result = run_gwas(
    genotype_matrix = genotype_matrix,
    phenotype_matrix = phenotype_matrix
  )
  
  cat(paste0("Simulation completed!\n"))
  pleio = list(
    n_sample = n_sample,
    n_variant_pleiotropic = n_variant_pleiotropic,
    n_variant_nonpleiotropic = n_variant_nonpleiotropic,
    n_variant_null = n_variant_null,
    eaf = eaf,
    n_phenotype = n_phenotype,
    desired_heritable_correlation_matrix = heritable_correlation_matrix,
    desired_nonheritable_correlation_matrix = nonheritable_correlation_matrix,
    heritability = heritability,
    crosstrait_heterogeneity = crosstrait_heterogeneity,
    withintrait_heterogeneity = withintrait_heterogeneity,
    random_crosstrait_heterogeneity = random_crosstrait_heterogeneity,
    random_withintrait_heterogeneity = random_withintrait_heterogeneity,
    genotype_matrix = genotype_matrix,
    efs_matrix = efs_matrix,
    heritable_phenotype_matrix = heritable_phenotype_matrix,
    defacto_heritable_correlation_matrix = defacto_heritable_correlation_matrix,
    nonheritable_phenotype_matrix = nonheritable_phenotype_matrix,
    defacto_nonheritable_correlation_matrix = defacto_nonheritable_correlation_matrix,
    phenotype_matrix = phenotype_matrix,
    gwas_result = gwas_result
  )
  
  return(pleio)
  
}

