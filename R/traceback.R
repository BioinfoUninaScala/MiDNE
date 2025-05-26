#' Traceback a node-node association
#'
#' @name traceback_link
#' @param extended_RWR A proximity matrix generated with gen_sim_mat_M or gen_sim_mat_MH functions by setting 'get_completeRWRmat = TRUE'. 
#' @param links_list A list of pairs of nodes (e.g., list(c('TP53', 'SOX2'), c('TP53', 'KLF5')))
#' @param reverse Logical, whether to include reverse direction (target --> source)
#' @param cpu The number of cores to use for parallel processing.
#' @return A tibble with columns: pair, source, target, proximity, layer
#' @export

traceback_link <- function(extended_RWR, links_list, reverse = FALSE, cpu = 1){
  if (!is.list(links_list)) {
    stop("Error: links_list must be a list!")
  }
  
  if (!is.matrix(extended_RWR)) {
    stop("Error: extended_RWR must be a numeric matrix!")
  }
  
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required!")
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required!")
  
  cl <- parallel::makeCluster(cpu)
  doParallel::registerDoParallel(cl)
  
  res_list <- foreach::foreach(i = seq_along(links_list), .combine = rbind,
                               .packages = c("dplyr", "stringr", "tibble")) %dopar% {
                                 pair <- links_list[[i]]
                                 
                                 if (length(pair) == 2 && is.vector(pair)) {
                                   source_node <- pair[1]
                                   target_node <- pair[2]
                                   
                                   res <- list()
                                   
                                   if (source_node %in% colnames(extended_RWR)) {
                                     source_target <- extended_RWR[grepl(paste0(target_node, "_"), rownames(extended_RWR)), source_node, drop = FALSE]
                                     res[[1]] <- tibble::tibble(
                                       pair = paste(pair, collapse = "_"),
                                       source = source_node,
                                       target = rownames(source_target),
                                       proximity = as.numeric(source_target)
                                     )
                                   }
                                   
                                   if (reverse && target_node %in% colnames(extended_RWR)) {
                                     target_source <- extended_RWR[grepl(paste0(source_node, "_"), rownames(extended_RWR)), target_node, drop = FALSE]
                                     res[[2]] <- tibble::tibble(
                                       pair = paste(pair, collapse = "_"),
                                       source = target_node,
                                       target = rownames(target_source),
                                       proximity = as.numeric(target_source)
                                     )
                                   }
                                   
                                   dplyr::bind_rows(res) %>%
                                     mutate(layer = stringr::str_extract(target, "(?<=_)[^_]+$"))
                                   
                                 } else {
                                   stop("Each element of the links_list must be a vector of length 2!")
                                 }
                               }
  
  parallel::stopCluster(cl)
  return(res_list)
}


#' Plot traceback output
#'
#' @name plot_traceback
#' @param traceback_RES Output of 'traceback_link()'
#' @return A ggplot object
#' @export

plot_traceback <- function(traceback_RES) {
  f_res <- traceback_RES %>%
    dplyr::mutate(
      link_base = paste(source, stringr::str_remove(target, "_[^_]+$"), sep = " -> "),
      pair_source = stringr::str_extract(pair, "^[^_]+"),
      direction = ifelse(source == pair_source, "forward", "reverse"),
      proximity_signed = ifelse(direction == "forward", proximity, -proximity)
    )
  
  ggplot2::ggplot(f_res, aes(x = layer, y = proximity_signed, fill = layer)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~ pair, scales = "free") +
    ggplot2::labs(x = "Layer", y = "Signed proximity", fill = "Layer",
                  title = "RWR proximity: direct vs. reverse transitions") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggplot2::geom_hline(yintercept = 0, color = "black")
}

