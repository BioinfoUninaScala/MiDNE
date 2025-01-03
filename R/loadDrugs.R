#' Load table of drug IDs
#'
#' @description This function loads a table of drug IDs included in the `MiDNE` package.
#'
#' @param drug_type A character defining the drug table to upload. Accepted values are:
#' 'all' to upload all drugs interacting with at least one protein in the DrugBank database, and 'FDA'
#' to upload only the pharmacologically active FDA-approved drugs. 
#' @return A dataframe with one column reporting drug IDs.
#' @examples
#' # Load all drugs interacting with proteins
#' all_drugs <- loadDrugs(drug_type = 'all')
#'
#' # Load FDA-approved drugs
#' fda_drugs <- loadDrugs(drug_type = 'FDA')
#' @export

loadDrugs <- function(drug_type = 'FDA') {
  # Valid values for drug_type
  valid_types <- c('all', 'FDA')
  
  # Ensure drug_type is valid
  if (!(drug_type %in% valid_types)) {
    stop(sprintf("The 'drug_type' parameter is invalid. Use one of the following values: %s",
                 paste(valid_types, collapse = ", ")))
  }
  
  # Load the appropriate drug table
  if (drug_type == 'all') {
    drug_path <- system.file("extdata", "data/pharmacological/ALLdrugs.RDS", package = "MiDNE")
    drugs <- readRDS(file = drug_path)
  } else if (drug_type == 'FDA') {
    drug_path <- system.file("extdata", "data/pharmacological/FDAdrugs.RDS", package = "MiDNE")
    drugs <- readRDS(file = drug_path)
  }
  
  return(drugs)
}