run_gwas = function(genotype_matrix, phenotype_matrix){
  
  gwas_coef = apply(phenotype_matrix, MARGIN = 2, function(ph){
    apply(genotype_matrix, MARGIN = 2, function(gt){
      l = lm(ph ~ gt)
      beta = coef(l)[2]
      p_value = summary(l)[["coefficients"]][2,4]
      paste(beta, p_value, sep = "|")
    })
  })
  
  beta_matrix = apply(gwas_coef, MARGIN = c(1,2), function(x){
    as.numeric(unlist(strsplit(x, split = "[|]"))[1])
  })
  
  pvalue_matrix = apply(gwas_coef, MARGIN = c(1,2), function(x){
    as.numeric(unlist(strsplit(x, split = "[|]"))[2])
  })
  
  result_list = list(
    beta_matrix = beta_matrix,
    pvalue_matrix = pvalue_matrix
  )
  
  return(result_list)
  
}