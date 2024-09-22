#' Fisher's Exact Test and post-hoc analysis
#' 
#' @name fisher_test_post_hoc
#' @import snow
#' @import utils
#' @import foreach
#' @import stats
#' @param matrix A matrix of dimensions genes X samples.
#' @param correction_method The method used to correct the p-value (either "bonferroni" or "fdr").
#' @param rows The number of rows associated with the pair of events of interest (e.g., row 4 indicates the co-occurrence of a molecular event, including co-methylation, co-amplification, and co-deletion).
#' @param cpu The number of cores to use for parallel processing.
#' @return A gene X gene log(expected/observed) matrix.
#' @export



fisher_test_post_hoc <- function(matrix, 
                                 cpu,
                                 rows, 
                                 correction_method){
  
  if (!requireNamespace("RVAideMemoire", quietly = TRUE)) {
    stop("The 'RVAideMemoire' package is required but not installed.")
  }
  if (!requireNamespace("doSNOW", quietly = TRUE)) {
    stop("The 'doSNOW' package is required but not installed.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The 'Matrix' package is required but not installed.")
  }
  
  cl <- snow::makeCluster(cpu)
  doSNOW::registerDoSNOW(cl)
  iterations <- base::nrow(matrix)
  pb <- utils::txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  output_list <- list()
  result <- foreach::foreach(i = 1:nrow(matrix),
                     .packages = c('tidyverse', 'RVAideMemoire','doSNOW', 'Matrix'),
                     .export = c('post_hoc_analysis_2'),
                     .options.snow=opts)  %dopar% {
                       vector <- c(rep(NA, nrow(matrix)))
                       for (j in i:nrow(matrix)) {
                         cont_table <- base::table(matrix[i, , drop = FALSE], matrix[j, ,drop = FALSE])
                         # Perform Fisher's exact test
                         fisher_p  <- stats::fisher.test(cont_table)$p.value
                         vector[j] <- ifelse(fisher_p < 0.05, post_hoc_analysis_2(cont_table, fisher_p, rows, correction_method), 0)
                       }
                       output_list[[i]] <- vector
                       
                     }
  
  base::close(pb)
  snow::stopCluster(cl)
  
  res_mat <- base::do.call(rbind, result)
  colnames(res_mat) <- rownames(res_mat) <- rownames(matrix)
  final_res_mat <- Matrix::forceSymmetric(res_mat)
  return(final_res_mat)
}


post_hoc_analysis_2 <- function(cont_table, fisher_p, rows, correction_method){
  
  if (!requireNamespace("RVAideMemoire", quietly = TRUE)) {
    stop("The 'RVAideMemoire' package is required but not installed.")
  }
  
  post_hoc <- RVAideMemoire::chisq.theo.multcomp(cont_table, p.method = correction_method)
  pval <- post_hoc$p.value[rows, 6]
  obs <- post_hoc$p.value$observed.Freq
  exp <- post_hoc$p.value$expected
  logFC <- ifelse(all(pval < 0.05), log2(base::mean(obs[rows])/exp[1]), 0)
  return(logFC)
}
