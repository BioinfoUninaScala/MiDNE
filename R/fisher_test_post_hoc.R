#' Fisher's Exact Test and post-hoc analysis
#' 
#' @name fisher_test_post_hoc
#' @importFrom doSNOW registerDoSNOW
#' @param matrix A matrix of dimensions genes X samples.
#' @param correction_method The method used to correct the p-value (either "bonferroni" or "fdr").
#' @param cpu The number of cores to use for parallel processing.
#' @param pth A numeric value, ranging from 0 and 1, that will be applied to the p-value of the Fisher's Exact Test and of post-hoc analysis.
#' @return A gene X gene matrix.
#' @export



fisher_test_post_hoc <- function(matrix, 
                                 cpu,
                                 correction_method,
                                 pth = 0.05
                                 ){
  
  # if (!requireNamespace("RVAideMemoire", quietly = TRUE)) {
  #   stop("The 'RVAideMemoire' package is required but not installed.")
  # }
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
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  result <- foreach::foreach(i = 1:nrow(matrix),
                             .packages = c('tidyverse', 'Matrix'),
                             .export = c('post_hoc_analysis_2', 'chisq_theo_multcomp', 'psignif'),
                             .options.snow = opts,
                             .combine = 'rbind') %dopar% {
                               sub_results <- list()  # raccoglie tutti i risultati per i-th riga
                               
                               for (j in i:nrow(matrix)) {
                                 gene_i <- rownames(matrix)[i]
                                 gene_j <- rownames(matrix)[j]
                                 cont_table <- base::table(matrix[i, , drop = FALSE], matrix[j, , drop = FALSE])
                                 
                                 fisher_p <- tryCatch(stats::fisher.test(cont_table)$p.value, error = function(e) NA)
                                 
                                 if (!is.na(fisher_p) && fisher_p < pth) {
                                   df <- post_hoc_analysis_2(cont_table, correction_method, fisher_p, pth, gene_i, gene_j)
                                   sub_results[[length(sub_results) + 1]] <- df
                                 }
                               }
                               
                               if (length(sub_results) > 0) {
                                 do.call(rbind, sub_results)
                               } else {
                                 NULL
                               }
                             }
  
  base::close(pb)
  snow::stopCluster(cl)
  
  final_network <- result %>% filter(source != dest)
  return(final_network)
}


post_hoc_analysis_2 <- function(cont_table, correction_method, fisher_p, pth, gene_i, gene_j){
  
  # if (!requireNamespace("RVAideMemoire", quietly = TRUE)) {
  #   stop("The 'RVAideMemoire' package is required but not installed.")
  # }
  
  expected = outer(rowSums(cont_table), colSums(cont_table), "*")/sum(cont_table)
  norm_expected <- expected/sum(expected)

  post_hoc <- chisq_theo_multcomp(x = cont_table, 
                                  p = norm_expected, 
                                  p.method = correction_method)
  adj_pval <- post_hoc$p.value2
    
  # if (any(diff_obs_exp[-1] > 0)) {
  if (any(adj_pval < pth)) {
    chisq <- stats::chisq.test(cont_table)
    stat <- unname(chisq$statistic)
    
    network <- data.frame(
      source = gene_i, 
      dest   = gene_j,
      weight = stat,
      adj_pval = fisher_p,
      
      # p-values
      adj_pval_00 = post_hoc$p.value$`Pr(>Chi)`[1],
      adj_pval_10 = post_hoc$p.value$`Pr(>Chi)`[2],
      adj_pval_01 = post_hoc$p.value$`Pr(>Chi)`[3],
      adj_pval_11 = post_hoc$p.value$`Pr(>Chi)`[4],
      
      # osservati
      obs_00 = post_hoc$p.value$observed.Freq[1],
      obs_10 = post_hoc$p.value$observed.Freq[2],
      obs_01 = post_hoc$p.value$observed.Freq[3],
      obs_11 = post_hoc$p.value$observed.Freq[4],
      
      # attesi
      exp_00 = post_hoc$p.value$expected[1],
      exp_10 = post_hoc$p.value$expected[2],
      exp_01 = post_hoc$p.value$expected[3],
      exp_11 = post_hoc$p.value$expected[4],
      
      # differenze
      diff_00 = post_hoc$p.value$observed.Freq[1] - post_hoc$p.value$expected[1],
      diff_10 = post_hoc$p.value$observed.Freq[2] - post_hoc$p.value$expected[2],
      diff_01 = post_hoc$p.value$observed.Freq[3] - post_hoc$p.value$expected[3],
      diff_11 = post_hoc$p.value$observed.Freq[4] - post_hoc$p.value$expected[4]
    )
    
  } else {
    network <- NULL
  }
  
  return(network)
}


#' Pairwise comparisons after a chi-squared test for given probabilities (adapted from RVAideMemoire::chisq.theo.multcomp)
#' 
#' @description Performs pairwise comparisons after a global chi-squared test for given probabilities.
#' 
#' @param x a contigency table
#' @param p theoretical proportions
#' @param p.method	method for p-values correction. See help of p.adjust.
#'
#' @return a table with the results of the pairwise comparisons
#' @export

chisq_theo_multcomp <- function (x, p = rep(1/length(x), length(x)), p.method = "fdr") 
{
  if (!all.equal(sum(p), 1)) {
    stop("sum of probabilities must be 1")
  }
  theo <- integer(length(x))
  chi2 <- integer(length(x))
  pval <- integer(length(x))
  for (i in 1:length(x)) {
    test <- suppressWarnings(stats::chisq.test(c(x[i], sum(x) - x[i]), p = c(p[i], 1 - p[i])))
    theo[i] <- as.numeric(test$expected[1])
    chi2[i] <- as.numeric(test$statistic)
    pval[i] <- as.numeric(test$p.value)
  }
  p.adj <- stats::p.adjust(pval, method = p.method)
  comp <- data.frame(observed = x, expected = theo, Chi = chi2, 
                     `Pr(>Chi)` = p.adj, ` ` = psignif(p.adj), stringsAsFactors = FALSE, 
                     check.names = FALSE)
  call <- match.call()
  dname.x <- if (length(call$x) == 1) {
    call$x
  }
  else {
    paste(call$x[1], "(", paste(call$x[-1], collapse = ","), 
          ")", sep = "")
  }
  dname.p <- if (length(call$p) == 1) {
    call$p
  }
  else {
    paste(call$p[1], "(", paste(call$p[-1], collapse = ","), 
          ")", sep = "")
  }
  dname <- paste(dname.x, " and ", dname.p, sep = "")
  result <- list(method = "chi-squared tests", data.name = dname, 
                 observed = x, expected = theo, p.adjust.method = p.method, 
                 statistic = chi2, p.value2 = p.adj, p.value = comp)
  class(result) <- "RV.multcomp"
  return(result)
}



#' Transform p-value significance in characters (adapted from RVAideMemoire)
#' 
#' @description transform p-value significance in characters.
#' 
#' @param p p-value
#'
#' @return a string
#' @export

psignif <- function(p) 
{
  result <- character(length(p))
  for (i in 1:length(p)) {
    if (p[i] != "NA") {
      if (as.numeric(p[i]) >= 0.1) {
        result[i] <- " "
      }
      else if (as.numeric(p[i]) < 0.1 & as.numeric(p[i]) >= 
               0.05) {
        result[i] <- "."
      }
      else if (as.numeric(p[i]) < 0.05 & as.numeric(p[i]) >= 
               0.01) {
        result[i] <- "*"
      }
      else if (as.numeric(p[i]) < 0.01 & as.numeric(p[i]) >= 
               0.001) {
        result[i] <- "**"
      }
      else if (as.numeric(p[i]) < 0.001) {
        result[i] <- "***"
      }
    }
    else {
      result[i] <- " "
    }
  }
  return(result)
}



# post_hoc_analysis_1 <- function(cont_table, fisher_p, rows, correction_method){
#   
#   if (!requireNamespace("RVAideMemoire", quietly = TRUE)) {
#     stop("The 'RVAideMemoire' package is required but not installed.")
#   }
#   
#   post_hoc <- RVAideMemoire::chisq.theo.multcomp(cont_table, p.method = correction_method)
#   
#   pval <- post_hoc$p.value[rows, 6]
#   obs <- post_hoc$p.value$observed.Freq
#   exp <- post_hoc$p.value$expected
#   logODD <- log2(base::mean(obs[rows])/exp[1])
#   
#   res <- ifelse(all(pval < pth) & logODD > 0, logODD, 0)
#   
#   #other_pval <- post_hoc$p.value[-rows, 6]
#   # if (all(pval < other_pval)) {
#   #   obs <- post_hoc$p.value$observed.Freq
#   #   exp <- post_hoc$p.value$expected
#   #   logFC <- log2(base::mean(obs[rows])/exp[1])
#   # } else {
#   #   logFC <- 0
#   # }
#   
#   return(res)
# }
