#' Load TCGA-BRCA omics networks
#'
#' @description This function loads specific omics networks included in the `MiDNE` package.
#'
#' @param net_type A character vector containing the names of omics networks to upload. Accepted values are: 
#' 'expression', 'proteomics', 'methylation', 'amp', 'del' and 'all'. You can provide a combination of the first four values 
#' (e.g., `c('expression', 'proteomics')`).
#' @return A named list of TCGA-BRCA networks.
#' The names of the list correspond to the omics types provided in `net_type`.
#' @export

loadOmicsNet<- function(net_type = c('expression', 'proteomics', 'methylation', 'del', 'amp', 'all')) {
  valid_types <- c('expression', 'proteomics', 'methylation', 'del', 'amp', 'all')
  net_type <- unique(net_type)
  
  if (!all(net_type %in% valid_types)) {
    stop(sprintf("The 'net_type' parameter is invalid. Use one or more of the following values: %s",
                 paste(valid_types, collapse = ", ")))
  }
  
  omics_path_list <- list()
  
  if (length(net_type) == 1 && net_type == 'all') {
    omics_path_list <- list(
      expr = system.file("extdata", "networks/biological/BRCA_filt_coexpr_network.RDS", package = "MiDNE"),
      meth = system.file("extdata", "networks/biological/BRCA_filt_cometh_network.RDS", package = "MiDNE"),
      amp  = system.file("extdata", "networks/biological/BRCA_filt_coamp_network.RDS", package = "MiDNE"),
      del  = system.file("extdata", "networks/biological/BRCA_filt_codel_network.RDS", package = "MiDNE"),
      prot = system.file("extdata", "networks/biological/BRCA_filt_prot_net.RDS", package = "MiDNE")
    )
  } else {
    # Caso: net_type contiene uno o più valori specifici
    if ('expression' %in% net_type) {
      omics_path_list[['expr']] <- system.file("extdata", "networks/biological/BRCA_filt_coexpr_network.RDS", package = "MiDNE")
    }
    if ('proteomics' %in% net_type) {
      omics_path_list[['prot']] <- system.file("extdata", "networks/biological/BRCA_filt_prot_net.RDS", package = "MiDNE")
    }
    if ('methylation' %in% net_type) {
      omics_path_list[['meth']] <- system.file("extdata", "networks/biological/BRCA_filt_cometh_network.RDS", package = "MiDNE")
    }
    if ('amp' %in% net_type) {
      omics_path_list[['amp']] <- system.file("extdata", "networks/biological/BRCA_filt_coamp_network.RDS", package = "MiDNE")
    }
    if ('del' %in% net_type) {
      omics_path_list[['del']] <- system.file("extdata", "networks/biological/BRCA_filt_codel_network.RDS", package = "MiDNE")
    }
  }
  
  OMICS_LIST <- lapply(omics_path_list, readRDS)
  return(OMICS_LIST)
}