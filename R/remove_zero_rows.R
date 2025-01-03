#' Remove rows containing only zeros from a matrix
#'
#' @description This function removes rows that contain only zeros from a given numeric matrix.
#'
#' @param matrix A numeric matrix (features x samples).
#' @return A matrix of the same type as \code{matrix}, with rows containing only zeros removed.
#' @examples
#' # Example matrix
#' mat <- matrix(c(1, 0, 3, 0, 0, 0, 7, 8, 9), nrow = 3, byrow = TRUE)
#' print(mat)
#' # Remove rows with all zeros
#' mat_no_zeros <- remove_zeros_rows(mat)
#' print(mat_no_zeros)
#' @export

remove_zeros_rows <- function(matrix) {
  if (!is.matrix(matrix)) {
    stop("The input must be a matrix.")
  }
  matrix[rowSums(abs(matrix)) > 0, , drop = FALSE]
}