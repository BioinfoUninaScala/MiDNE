#' Generate co-CNV network
#' 
#' @name gen_coCNVnet
#' @param omics_matrix A matrix of dimensions genes X samples representing Copy Number Variation (CNV) data.
#' @param correction_method The method used to correct the p-value of the statistical test (either "bonferroni" or "fdr").
#' @param cpu The number of cores to use for parallel processing.
#' @return A list containing two undirected networks (co-amplification and co-deletion), each represented as a 3-column table (source, dest, weight).
#' @export


gen_coCNVnet.R <- function( 
                            omics_matrix, 
                            correction_method = NULL,
                            cpu = 1
                            )
{
  cnv_mat <- as.data.frame(omics_matrix)
  amp <- ifelse(cnv_mat >= 1, 1, 0)
  del <- ifelse(cnv_mat <= -1, 1, 0)
  
  amp_mat <- structure(factor(amp), dim = dim(amp), class = c('matrix', 'factor'))
  del_mat <- structure(factor(del), dim = dim(del), class = c('matrix', 'factor'))
  rownames(amp_mat) <-  rownames(amp)
  rownames(del_mat) <-  rownames(del)
  colnames(amp_mat) <-  colnames(amp)
  colnames(del_mat) <-  colnames(del)
  
  coamp <- fisher_test_post_hoc(matrix = amp_mat, cpu = cpu, rows = 4, correction_method)
  codel <- fisher_test_post_hoc(matrix = del_mat, cpu = cpu, rows = 4, correction_method)
  net_matrix <- list('coamp' = coamp, 'codel' = codel)
  
  result <- list()
  for (mat_name in names(net_matrix)) {
    res_mat <- as.matrix(net_matrix[[mat_name]])
    network <- get_adjList(res_mat)
    network <- network %>% dplyr::mutate_if(is.factor, as.character)
    if (nrow(network) > 0){
      result[[mat_name]] <- network
    }
  }
 return(result) 
}
