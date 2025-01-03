#' Get intersection between omics matrices
#'
#' @param omics_list a list of matrices with features on columns
#' @return a list of matrices sharing the same column ids
#' @export

get_intersection_matrices <- function(omics_list) {
  col_inter <- base::Reduce(intersect, base::lapply(omics_list, colnames))
  base::lapply(omics_list, function(x) x[,col_inter])
}