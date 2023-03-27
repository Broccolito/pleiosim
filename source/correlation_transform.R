correlation_transform = function(efs_matrix_template, genotype_matrix, n_phenotype,
                                 n_variant_pleiotropic, heritable_correlation_matrix,
                                 variable_names){
  
  determine_transform_factor = function(m, efs_mat, gt_mat, 
                                        npx, gc_mat){
    
    px_efs_mat = efs_mat[,1:npx]
    px_efs_mat = apply(px_efs_mat, MARGIN = 2, function(x){x*as.matrix(m)})
    
    efs_mat[,1:npx] =  px_efs_mat
    
    pheno_mat = t(efs_mat %*% t(gt_mat))
    pred_corr_mat = cor(pheno_mat)
    loss = sum((pred_corr_mat - gc_mat)^2)
    
    return(loss)
  }
  
  m = optim(par = rep(1, n_phenotype), 
            fn = determine_transform_factor, 
            efs_mat = efs_matrix_template,
            gt_mat = genotype_matrix, 
            npx = n_variant_pleiotropic,
            gc_mat = heritable_correlation_matrix)$par
  
  efs_matrix_template[,1:n_variant_pleiotropic] = apply(efs_matrix_template[,1:n_variant_pleiotropic], 
                                                        MARGIN = 2, function(x){x*as.matrix(m)})
  
  colnames(efs_matrix_template) = variable_names[["variant_names"]]
  rownames(efs_matrix_template) = variable_names[["phenotype_names"]]
  
  return(efs_matrix_template)
  
}
