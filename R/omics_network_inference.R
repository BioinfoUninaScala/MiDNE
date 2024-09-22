#' Generate co-expression network
#' 
#' @name omics_network_inference
#' @param omics_matrix a genesXsamples matrix.
#' @param omics_type character that indicates the molecular domain of the omics matrix (Transcriptomics, Proteomics, Genetics, Epigenetics).
#' @param correction_method method for correcting the p-value of the statistical test (bonferroni or fdr).
#' @param th a numeric value to binarize a DNA methylation matrix (expressed as beta-values).
#' @param cpu number of cores to work in parallel.
#' @return an undirected biological network, defined as a 3-columns table (source, dest, weight).
#' @export 


omics_network_inference <- function(omics_matrix, 
                                    omics_type = NULL,
                                    #inference_method = NULL,
                                    correction_method = NULL, 
                                    th = NA, 
                                    cpu = 1)
  {
  
  if (omics_type == 'Transcriptomics'){
    result <- gen_coExpressionNet(omics_matrix = omics_matrix,
                                  correction_method = correction_method,
                                  cpu = cpu)
  } else if (omics_type == 'Proteomics'){
    result <- gen_coAbundanceNet(omics_matrix = omics_matrix,
                                 correction_method = correction_method,
                                 cpu = cpu)
  } else if (omics_type == 'Epigenomics'){
    result <- gen_coDNAmethNet(omics_matrix = omics_matrix, 
                               correction_method = correction_method, 
                               th = th, cpu = cpu)
  } else {
    result <- gen_coCNVnet(omics_matrix = omics_matrix, 
                           correction_method = correction_method, 
                           cpu = cpu)
  }
  
  return(result)
  }

