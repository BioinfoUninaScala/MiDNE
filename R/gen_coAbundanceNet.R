#' Generate co-abundance network
#' 
#' @name gen_coAbundanceNet
#' @param omics_matrix A matrix of dimensions genes X samples representing proteomics data.
#' @param correction_method The method used to correct the p-value of the statistical test (either "bonferroni" or "fdr").
#' @param cpu The number of cores to use for parallel processing.
#' @param pth A numeric value, ranging from 0 and 1, that will be applied to the p-value of Spearman Correlation test.
#' @return An undirected co-expression network, represented as a 3-column table (source, dest, weight).
#' @export


gen_coAbundanceNet <- function( 
                                  omics_matrix, 
                                  correction_method = NULL,
                                  cpu = 1,
                                  pth = 0.05
                                  )
{
  network <- get_filtered_corMat_by_adj_pval(input_matrix = omics_matrix, 
                                                cpu = cpu,
                                                cor_method = 'spearman', 
                                                adj_method = correction_method,
                                                pth = pth
                                                )
  result <- network %>% tibble::as_tibble()
  return(result)
}
