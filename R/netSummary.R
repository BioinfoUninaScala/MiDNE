#' Network summary
#'
#' @description Provides basic information about a network object, including number of edges, unique nodes, and max node degree.
#' 
#' @importFrom dplyr as_tibble arrange
#' @param network A two-column data frame representing the edge list of a network.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{dim}: dimensions of the network (rows, columns)
#'   \item \code{edges}: total number of edges
#'   \item \code{nodes}: total number of unique nodes
#'   \item \code{maxConnections}: maximum number of connections per node
#' }
#'
#' @export


netSummary <- function(network){
  
  nodes <- c(network[[1]], network[[2]])
  degree <- table(nodes) %>%
    dplyr::as_tibble(.) %>%
    dplyr::arrange(dplyr::desc(n))
  netSummary_list <- list( 'dim'= dim(network),
                           'edges'= dim(network)[1],
                           'nodes'= length(unique(nodes)),
                           'maxConnections'= degree[[1, 'n']]
  )
  return(netSummary_list)
}

