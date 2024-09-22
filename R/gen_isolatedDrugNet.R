#' Generate isolated drug network
#'
#' @name drug_network_inference
#' @param drug_input A single-column dataframe containing the drug names.
#' @param drug_type The method used to construct the drug network. Using the "Virtual Nodes" strategy, each drug is linked to a virtual node.
#' @param cpu The number of cores to use for parallel processing.
#' @return An unweighted isolated drug network, represented as a 3-column table (source = drug, dest = virtual node, weight = 1).
#' @export


drug_network_inference <- function(
                                    drug_input,
                                    drug_type,
                                    #dist, 
                                    cpu = 1
    )
  {
  if (drug_type == "Virtual Nodes"){
    colnames(drug_input)[1] <- 'source'
    f_drug_table <- drug_input %>% tibble::as_tibble(.) %>% dplyr::distinct() %>% dplyr::arrange(source) %>%
      dplyr::mutate(dest = paste0('v_', 1:nrow(drug_input)),
             weight = rep(1, nrow(drug_input)))
    return(f_drug_table)
  }
}
