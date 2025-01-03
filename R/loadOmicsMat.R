#' Load TCGA-BRCA omics matrices.
#'
#' @description This function loads specific omics matrices from the TCGA-BRCA dataset included in the `MiDNE` package.
#'
#' @param mat_type A character vector containing the names of omics matrices to upload. Accepted values are: 
#' 'expression', 'proteomics', 'methylation', 'scnv', and 'all'. You can provide a combination of the first four values 
#' (e.g., `c('expression', 'proteomics')`).
#' @return A named list of TCGA-BRCA matrices with genes on rows and samples on columns.
#' The names of the list correspond to the omics types provided in `mat_type`.
#' @examples
#' # Load expression and proteomics data
#' omics_data <- loadOmicsMat(mat_type = c('expression', 'proteomics'))
#'
#' # Load all data types
#' all_data <- loadOmicsMat(mat_type = 'all')
#' @export

loadOmicsMat <- function(mat_type = c('expression', 'proteomics', 'methylation', 'scnv', 'all')) {
  valid_types <- c('expression', 'proteomics', 'methylation', 'scnv', 'all')
  mat_type <- unique(mat_type)
  
  if (!all(mat_type %in% valid_types)) {
    stop(sprintf("The 'mat_type' parameter is invalid. Use one or more of the following values: %s",
                 paste(valid_types, collapse = ", ")))
  }
  
  omics_path_list <- list()
  
  if (length(mat_type) == 1 && mat_type == 'all') {
    omics_path_list <- list(
      expr = system.file("extdata", "data/biological/BRCA_expr_HiSeq.RDS", package = "MiDNE"),
      meth = system.file("extdata", "data/biological/BRCA_Methylation_Meth450.RDS", package = "MiDNE"),
      cnv  = system.file("extdata", "data/biological/BRCA_SCNA.RDS", package = "MiDNE"),
      prot = system.file("extdata", "data/biological/BRCA_proteome_CDAP.RDS", package = "MiDNE")
    )
  } else {
    # Caso: mat_type contiene uno o più valori specifici
    if ('expression' %in% mat_type) {
      omics_path_list[['expr']] <- system.file("extdata", "data/biological/BRCA_expr_HiSeq.RDS", package = "MiDNE")
    }
    if ('proteomics' %in% mat_type) {
      omics_path_list[['prot']] <- system.file("extdata", "data/biological/BRCA_proteome_CDAP.RDS", package = "MiDNE")
    }
    if ('methylation' %in% mat_type) {
      omics_path_list[['meth']] <- system.file("extdata", "data/biological/BRCA_Methylation_Meth450.RDS", package = "MiDNE")
    }
    if ('scnv' %in% mat_type) {
      omics_path_list[['scnv']] <- system.file("extdata", "data/biological/BRCA_SCNA.RDS", package = "MiDNE")
    }
  }
  
  OMICS_LIST <- lapply(omics_path_list, readRDS)
  return(OMICS_LIST)
}