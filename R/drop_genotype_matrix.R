#' @export
drop_genotype_matrix = function(pleio){
  pleio[["genotype_matrix"]] = "Genotype matrix has been dropped from the results to reduce object size..."
  return(pleio)
}
