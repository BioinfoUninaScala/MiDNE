#' Function to infer correlation network filtered by corrected p-value
#'
#' @name get_filtered_corMat_by_adj_pval
#' @import stats
#' @import utils
#' @param input_matrix A matrix of dimensions genes X samples.
#' @param cor_method The correlation method to use (e.g., "pearson", "spearman").
#' @param adj_method The method used to correct the p-value (either "bonferroni" or "fdr").
#' @param cpu The number of cores to use for parallel processing.
#' @return A gene X gene correlation matrix, filtered by a corrected p-value lower than 0.05.
#' @export


get_filtered_corMat_by_adj_pval <- function(
    input_matrix,
    cpu = 10, 
    cor_method = 'spearman', 
    adj_method = 'fdr'){
  
  if (!requireNamespace("snow", quietly = TRUE)) {
    stop("The 'snow' package is required but not installed.")
  }
  if (!requireNamespace("doSNOW", quietly = TRUE)) {
    stop("The 'doSNOW' package is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but not installed.")
  }
  
  # input_matrix --> features x samples
  
  matrix = input_matrix
  pval_mat <- matrix(0, nrow(matrix), nrow(matrix))
  rownames(pval_mat) <- colnames(pval_mat) <- rownames(matrix)
  
  message(paste('Constructing p-value matrix of correlation (method: ', cor_method, ') ...'))
  
  cl <- snow::makeCluster(cpu)
  doSNOW::registerDoSNOW(cl)
  iterations <- nrow(matrix)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  output_list <- list()
  result <- foreach::foreach(i = 1:nrow(matrix), 
                    .options.snow=opts) %dopar% {
                      vector <- c(rep(NA, nrow(matrix)))
                      for (j in i:nrow(matrix)){
                        x <- matrix[i,]
                        y <- matrix[j,]
                        res <- cor.test(x,y, method = cor_method)
                        vector[j] <- res$p.value
                      }
                      output_list[[i]] <- vector
                    }
  base::close(pb)
  snow::stopCluster(cl)
  
  pval_mat <- base::do.call(rbind, result)
  
  message(paste('Adjusting p-value for multiple testing (method: ', adj_method, ') ...'))
  pval_vec <- pval_mat[upper.tri(pval_mat)]
  adj_pval_vec <- stats::p.adjust(pval_vec, 
                           adj_method, 
                           length(pval_vec) + nrow(matrix))
  
  padj_mat <- matrix(0, nrow(matrix), nrow(matrix))
  padj_mat[upper.tri(padj_mat, diag=FALSE)] <- adj_pval_vec
  
  
  final_padj_mat <- Matrix::forceSymmetric(padj_mat)
  colnames(final_padj_mat) <- rownames(final_padj_mat) <- rownames(matrix)
  
  MASK <- ifelse(final_padj_mat < 0.05, TRUE, FALSE)
  
  message(paste('Constructing correlation matrix (method: ', cor_method, ') ...'))
  cor_mat <- stats::cor(t(input_matrix), method = cor_method)
  
  message(paste('Filtering correlation matrix by adjusted p-value matrix ...'))
  filtered_corMat_by_padj <- base::replace(cor_mat, !MASK, 0)
  return(filtered_corMat_by_padj)
}

