#' Generate co-expression network
#'
#' @name gen_coExpressionNet
#' @param omics_matrix A matrix of dimensions genes X samples representing transcriptomics data.
#' @param correction_method The method used to correct the p-value of the statistical test (either "bonferroni" or "fdr").
#' @param cpu The number of cores to use for parallel processing.
#' @param pth A numeric value, ranging from 0 and 1, that will be applied to the p-value of Pearson Correlation test.
#' @return An undirected co-expression network, represented as a 3-column table (source, dest, weight).
#' @export


gen_coExpressionNet <- function( 
                                  omics_matrix, 
                                  correction_method = NULL,
                                  cpu = 1, pth = 0.05
                                  )
{
  net_matrix <- get_filtered_corMat_by_adj_pval(input_matrix = omics_matrix, 
                                                cpu = cpu,
                                                cor_method = 'pearson', 
                                                adj_method = correction_method,
                                                pth = pth)
  
  network <- get_adjList(net_matrix)
  result <- network %>% dplyr::mutate_if(is.factor, as.character)
  return(result)
}
