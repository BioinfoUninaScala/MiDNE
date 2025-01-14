#' Load gene-drug network
#'
#' @description This function loads gene-drug network included in the `MiDNE` package.
#' @import readr
#' @return A dataframe with three columns reporting genes as source nodes, drugs as target nodes and the weights of the edges.
#' @export

loadBipartiteNet <- function() {
  
  file_path <- system.file("extdata", "networks/bipartite/FDA_active_DRUGBANK_bnet.RDS", package = "MiDNE")
  file <- readr::read_delim(file = file_path, delim = ',')

  return(file)
}