simulate_genotypes = function(n_sample,
                              eaf,
                              variable_names){
  
  f_hom = eaf^2
  f_het = 2*eaf*(1-eaf)
  f_wt = (1-eaf)^2
  
  allele_counts = rbind(
    round(n_sample * f_hom),
    round(n_sample * f_het),
    round(n_sample * f_wt)
  )
  
  allele_counts[3,] = allele_counts[3,] - (apply(allele_counts, 2, sum) - n_sample)
  
  genotype_matrix = mapply(function(n_hom, n_het, n_wt){
    sample(c(rep(2, n_hom), rep(1, n_het), rep(0, n_wt)),
           size = n_sample, replace = FALSE)
  }, allele_counts[1,], allele_counts[2,], allele_counts[3,])
  
  colnames(genotype_matrix) = variable_names[["variant_names"]]
  rownames(genotype_matrix) = variable_names[["subject_names"]]
  
  return(genotype_matrix)
  
}

### NOT RUN
# simulate_genotypes(n_sample = 10000, eaf = rep(0.4, 1000))