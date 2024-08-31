#############################################################################################################
#                                              NETWORK INFERENCE                                            #
#############################################################################################################

## PACKAGES
library(tidyverse)
library(reshape2)
source('MiDNE/R/myfunctions.R')
conflicts_prefer(dplyr::filter)
conflicts_prefer(snow::makeCluster)


omics_network_inference <- function(omics_matrix, omics_type = NULL, inference_method = NULL, correction_method = NULL, th = NA, cpu = 1){
  # omics_matrix = (samples x features)
  
  if (inference_method == 'Pearson Correlation Coefficient'){
    net_matrix <- get_filtered_corMat_by_adj_pval(input_matrix = omics_matrix, 
                                                                 cpu = cpu,
                                                                 cor_method = 'pearson', 
                                                                 adj_method = correction_method)
    
    network <- get_adjList(net_matrix)
    result <- network %>% mutate_if(is.factor, as.character)
    
  } else if (inference_method == 'Spearman Correlation Coefficient'){
    net_matrix <- get_filtered_corMat_by_adj_pval(input_matrix = omics_matrix, 
                                                                 cpu = cpu,
                                                                 cor_method = 'spearman', 
                                                                 adj_method = correction_method)
    network <- get_adjList(net_matrix)
    result <- network %>% mutate_if(is.factor, as.character)
    
  } else {
    if (omics_type == 'Epigenomics'){
      methMat <- ifelse(omics_matrix > th, 1, 0)
      meth_Mat <- structure(factor(methMat), dim = dim(methMat), class = c('matrix', 'factor'))
      rownames(meth_Mat) <- rownames(methMat)
      colnames(meth_Mat) <- colnames(methMat)
      
      cometh <- fisher_test_post_hoc(matrix = meth_Mat, cpu = cpu, rows = 4, correction_method = correction_method )
      cometh1 <- as.matrix(cometh)
      cometh2 <- ifelse(cometh1 == -Inf, 0, cometh1)
      colnames(cometh2) <- rownames(cometh2) <- rownames(omics_matrix)
      network <- get_adjList(cometh2)
       
      result <- network %>% mutate_if(is.factor, as.character)
      
    } else if (omics_type == 'Genomics'){
      cnv_mat <- as.data.frame(omics_matrix)
      amp <- ifelse(cnv_mat >= 1, 1, 0)
      del <- ifelse(cnv_mat <= -1, 1, 0)
      
      amp_mat <- structure(factor(amp), dim = dim(amp), class = c('matrix', 'factor'))
      del_mat <- structure(factor(del), dim = dim(del), class = c('matrix', 'factor'))
      rownames(amp_mat) <-  rownames(amp)
      rownames(del_mat) <-  rownames(del)
      colnames(amp_mat) <-  colnames(amp)
      colnames(del_mat) <-  colnames(del)
      
      coamp <- fisher_test_post_hoc(matrix = amp_mat, cpu = cpu, rows = 4, correction_method)
      codel <- fisher_test_post_hoc(matrix = del_mat, cpu = cpu, rows = 4, correction_method)
      net_matrix <- list('coamp' = coamp, 'codel' = codel)
      
      result <- list()
      for (mat_name in names(net_matrix)) {
        res_mat <- as.matrix(net_matrix[[mat_name]])
        network <- get_adjList(res_mat)
        network <- network %>% mutate_if(is.factor, as.character)
        if (nrow(network) > 0){
          result[[mat_name]] <- network
        }
      }
    }
  }
  
  return(result)
}


drug_network_inference <- function(drug_input, drug_type, dist, cpu){
  if (drug_type == "Virtual Nodes"){
      colnames(drug_input)[1] <- 'source'
      f_drug_table <- drug_input %>% as_tibble(.) %>% unique() %>% arrange(source) %>%
                                     mutate(dest = paste0('v_', 1:nrow(drug_input)),
                                            weight = rep(1, nrow(drug_input)))
      return(f_drug_table)
  }
}
