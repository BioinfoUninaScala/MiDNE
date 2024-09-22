#' Generate co-methylation network
#' 
#' @name gen_coDNAmethNet
#' @param omics_matrix A matrix of dimensions genes X samples representing DNA methylation data (beta-values).
#' @param correction_method The method used to correct the p-value of the statistical test (either "bonferroni" or "fdr").
#' @param th A positive integer used to binarize the input DNA methylation matrix.
#' @param cpu The number of cores to use for parallel processing.
#' @return A list containing two undirected networks (co-amplification and co-deletion), each defined as a 3-column table (source, dest, weight).
#' @export


gen_coDNAmethNet <- function( 
                                omics_matrix, 
                                correction_method = NULL,
                                th = 0.3,
                                cpu = 1
                                )
{
  methMat <- ifelse(omics_matrix > th, 1, 0)
  meth_Mat <- base::structure(factor(methMat), dim = dim(methMat), class = c('matrix', 'factor'))
  rownames(meth_Mat) <- rownames(methMat)
  colnames(meth_Mat) <- colnames(methMat)
  
  cometh <- fisher_test_post_hoc(matrix = meth_Mat, cpu = cpu, rows = 4, correction_method = correction_method )
  cometh1 <- as.matrix(cometh)
  cometh2 <- ifelse(cometh1 == -Inf, 0, cometh1)
  colnames(cometh2) <- rownames(cometh2) <- rownames(omics_matrix)
  network <- get_adjList(cometh2)
  
  result <- network %>% dplyr::mutate_if(is.factor, as.character)
  return(result)
}
