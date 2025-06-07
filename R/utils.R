# ################################################################################
# ############################      MY FUNCTIONS      ############################
# ################################################################################


# get_adj_list <- function(adj_mat, th){
#
#    # -------------------------------------------------------------------------
#    #  This function outputs an adjacency list given an adjacency matrix and
#    #  a double as cut-off.
#    # -------------------------------------------------------------------------
#     net_0 <- ifelse(adj_mat <= -th | adj_mat >= th, adj_mat, 0)
#     net_0[upper.tri(net_0)] <- 999
#     net_1 <- reshape2::melt(net_0)
#     net_2 <- net_1 %>%
#                         dplyr::filter(value != 999) %>%
#                         dplyr::filter(value != 0) %>%
#                         dplyr::filter(Var1 != Var2) %>%
#                         as_tibble(.)
#     names(net_2) <- c('source', 'dest', 'weight')
#     #net <- net_2 %>% as_tibble() %>%
#     #                      mutate_if(., is.factor, as.character) %>%
#     #                      separate_rows(source) %>%
#     #                      separate_rows(dest)
#     return(net_2)}
# 
# 
# 
# netSummary_1 <- function(network){
#   # -------------------------------------------------------------------------
#   #  This function outputs a summary for a given network, like the network
#   #  dimension, the number of genes and edges, the nodes with the highest degree
#   #  and their degree.
#   # -------------------------------------------------------------------------
#   nodes <- c(network$source, network$dest)
#   degree <- table(nodes) %>%
#                          as_tibble(.) %>%
#                          arrange(desc(n))
#   top10 <- c(degree[c(1:10), 'nodes'])
#   netSummary_list <- list( 'dim'= dim(network),
#                            'edges'= dim(network)[1],
#                            'nodes'= length(unique(nodes)),
#                            'maxConnections'= degree[[1, 'n']],
#                            'maxDegreeNodes'= top10$nodes
#                           )
#   return(netSummary_list)
# }
# 
# 
# 
# from_adjMat_to_Ugraph <- function(adj_mat, th){
#   # -------------------------------------------------------------------------
#   #  This function outputs an undirected graph from an adjacency matrix.
#   #     INPUT: an adjacency matrix (features x features) and a threshold
#   #
#   #     OUTPUT: an undirected igraph object
#   # -------------------------------------------------------------------------
# 
#      boolNet <- ifelse(adj_mat >= th | adj_mat <= -th, 1, 0)
#      graph <- graph.adjacency(boolNet, mode = "undirected",
#                               add.colnames = NULL, diag = FALSE)
#      return(graph)
# }
# 
# 
# 
# get_hubs <- function(graph, clusters){
#   # -------------------------------------------------------------------------
#   #  This function finds the hub genes for every cluster of a graph.
#   #     INPUT: a graph object created with igraph and
#   #            a cluster object created with clustering algorithm, such as
#   #            louvain_cluster
#   #     OUTPUT: a list of hubs gene for each cluster
#   # -------------------------------------------------------------------------
#   hub <- list()
#   for (i in 1:length(clusters)) {
#     subgraph_nodes <- V(graph)[membership(clusters) == i]
#     subgraph <- induced_subgraph(graph, subgraph_nodes)
#     degrees <- degree(subgraph)
#     sort_deg <- sort(degrees, decreasing=TRUE)
#     top_genes <- names(sort_deg[sort_deg >= max(sort_deg)])
#     hub[[i]] <- top_genes
#   }
#   return(hub)
# }
# 
# 
# 
# clusterAnno <- function(clusters, hub){
#   # -------------------------------------------------------------------------
#   #  This function creates an annotation table keeping information on the gene
#   #  belonging to a cluster and the gene identity as an hub.
#   #     INPUT: a cluster object and a list of lists containing the hub genes
#   #            for each cluster in clusters
# 
#   #     OUTPUT: a tibble (genes x 3) with three columns: gene name, cluster id
#   #             and a boolean column in order to identify the hub genes
#   # -------------------------------------------------------------------------
# 
#   clust_tab <- clusters %>% Reduce(c, .) %>% as_tibble(.) %>%
#                             mutate('cluster' = membership(clusters)) %>%
#                             separate_rows(value)
#   colnames(clust_tab)[1] <- 'Gene'
#   pre_clust_anno <- clust_tab[!duplicated(clust_tab$Gene), ]
# 
#   hubs <- Reduce(c, hub)
#   clust_anno <- pre_clust_anno %>%
#                           mutate('hub'= if_else(Gene %in% hubs, 'TRUE', 'FALSE'))
#   return(clust_anno)
# }
# 
# 
# 
# get_padj_mat <- function(matrix){
#   ## matrix has samples on rows and genes on columns
# 
#   pval_mat <- matrix(0, ncol(matrix), ncol(matrix))
#   rownames(pval_mat) <- colnames(pval_mat) <- colnames(matrix)
# 
#   #options(width =  ncol(matrix))
#   #n <- ncol(matrix)
#   for (i in 1:ncol(matrix)){
#     #extra <- nchar('||100%')
#     #width <- options()$width
#     #step <- round(i / n * (width - extra))
#     #text <- sprintf('|%s%s|% 3s%%', strrep('=', step),
#     #                strrep(' ', width - step - extra), round((i / n )* 100))
#     #cat(text)
#     #Sys.sleep(0.05)
#     #cat(if (i == n) '\n' else '\014')
# 
#     for (j in i:ncol(matrix)){
#       x <- matrix[,i]
#       y <- matrix[,j]
#       res <- cor.test(x,y, method = 'spearman')
#       pval_mat[i,j] <- pval_mat[j,i] <- res$p.value
#     }
#   }
# 
#   pval_vec <- pval_mat[upper.tri(pval_mat)]
#   adj_pval_vec <- p.adjust(pval_vec, 'fdr', length(pval_vec)+ ncol(matrix))
# 
#   padj_mat <- matrix(0, ncol(matrix), ncol(matrix))
#   padj_mat[upper.tri(padj_mat, diag=FALSE)] <- adj_pval_vec
# 
#   p_adj_mat <- padj_mat + t(padj_mat)
#   colnames(p_adj_mat) <- rownames(p_adj_mat) <- colnames(matrix)
#   return(p_adj_mat)
# }
# 
# 
# 
# from_cpg_to_gene_promoter <- function(annotation.table, matrix, filter_list){
# 
#     geneAnno <- as_tibble(annotation.table) %>%
#         dplyr::filter(Name %in% rownames(matrix)) %>%
#         mutate_at('UCSC_RefGene_Group', ~gsub("'", "", .)) %>%
#         dplyr::select(Name, UCSC_RefGene_Name,
#                  UCSC_RefGene_Group, Regulatory_Feature_Group)
# 
#     geneAnnotation <- geneAnno %>%
#         mutate(UCSC_RefGene_Name =  strsplit(UCSC_RefGene_Name, split=';', fixed=TRUE),
#                UCSC_RefGene_Group =  strsplit(UCSC_RefGene_Group, split=';', fixed=TRUE)) %>%
#         unnest(c(UCSC_RefGene_Name, UCSC_RefGene_Group))
# 
# 
#     Annotation <- geneAnnotation %>%
#                              dplyr::filter( UCSC_RefGene_Group %in% filter_list) %>%
#                              dplyr::select(Name, UCSC_RefGene_Name)  %>%
#                              distinct(.) %>%
#                              arrange(Name)
#     
#     return(Annotation)
#  }


#' Convert adjacency matrix to adjacency list
#'
#' This function takes an adjacency matrix and returns an adjacency list.
#' It removes the upper triangle, zeros, and self-loops.
#'
#' @param adj_mat A square adjacency matrix where rows and columns correspond to nodes.
#'
#' @return A tibble with three columns: `source`, `dest`, and `weight`.
#' 
#' @importFrom reshape2 melt
#' @importFrom dplyr filter as_tibble
#' @export

get_adjList <- function(adj_mat){
  # -------------------------------------------------------------------------
  #  This function outputs an adjacency list given an adjacency matrix
  # -------------------------------------------------------------------------
  net_0 <- adj_mat
  net_0[upper.tri(net_0)] <- 999
  net_1 <- reshape2::melt(net_0)
  net_2 <- net_1 %>%
    dplyr::filter(value != 999) %>%
    dplyr::filter(value != 0) %>%
    dplyr::filter(Var1 != Var2) %>%
    dplyr::as_tibble()
  names(net_2) <- c('source', 'dest', 'weight')
  return(net_2)}



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
    dplyr::as_tibble %>%
    dplyr::arrange(dplyr::desc(n))
  netSummary_list <- list( 'dim'= dim(network),
                           'edges'= dim(network)[1],
                           'nodes'= length(unique(nodes)),
                           'maxConnections'= degree[[1, 'n']]
  )
  return(netSummary_list)
}

