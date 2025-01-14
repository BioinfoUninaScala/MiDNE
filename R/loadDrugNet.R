#' Load table of drug IDs
#'
#' @description This function loads a table of drug IDs included in the `MiDNE` package.
#'
#' @param net_type A character defining the drug network to upload. Accepted values are:
#' 'vn' to upload the drug network generated with the "virtual nodes" startegy. 
#' @return A dataframe with three columns: the source column reports the drug IDs,
#'  the second one the virtual nodes, and the third one the weights of the edges.
#' @export

loadDrugNet <- function(net_type = 'vn') {
  valid_types <- c('vn')
  
  if (!(net_type %in% valid_types)) {
    stop(sprintf("The 'net_type' parameter is invalid. Use one of the following values: %s",
                 paste(valid_types, collapse = ", ")))
  }
  
  # Load the appropriate drug table
  if (net_type == 'vn') {
    drug_path <- system.file("extdata", "networks/pharmacological/FDAdrugs_net.csv", package = "MiDNE")
    drug_net <- readr::read_delim(file = drug_path, delim = ',' )
  }
  
  return(drug_net)
}