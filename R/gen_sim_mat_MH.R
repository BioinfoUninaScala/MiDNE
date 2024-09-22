#################################################################################################################
#                   Random Walk with Restart on Multiplex and Heterogeneous biological Networks                 #
#################################################################################################################

#' Generate similarity matrix
#' Adapted from:
#' https://github.com/LPioL/MultiVERSE/blob/master/RWR/Functions_RWRMH.R
#' https://github.com/LPioL/MultiVERSE/blob/master/RWR/GenerateSimMatrix_MH.R
#'
#' Research group:  Léo Pio-Lopez, Alberto Valdeolivas, Laurent Tichit, Élisabeth Remy, Anaïs Baudot
#' Publication:     https://arxiv.org/abs/2008.10085
#'
#' @name generate_sim_mt_MH
#' @import utils
#' @import methods
#' @import igraph
#' @import foreach
#' @param network1 a biological network
#' @param network2 a pharmacological network
#' @param network_B a bipartite network
#' @param tau1 tau1
#' @param tau2 tau2
#' @param restart A real in the range of 0-1, it is the probability to restart the algorithm in the starting point
#' @param delta1 delta
#' @param delta2 delta
#' @param lambda lambda
#' @param layer_transition_1 cond_jump
#' @param layer_transition_2 cond_jump
#' @param jump_neighborhood_1 A boolean, if false, the nodes between omics will be connected if they are the same, if true, the nodes will be connected with itself and its neighborhood but in the other omics
#' @param jump_neighborhood_2 A boolean, if false, the nodes between omics will be connected if they are the same, if true, the nodes will be connected with itself and its neighborhood but in the other omics
#' @param weighted_multiplex_1 A boolean, if true, the edge between omics will be weighted. It is considered only if jump_neighborhood is true
#' @param weighted_multiplex_2 A boolean, if true, the edge between omics will be weighted. It is considered only if jump_neighborhood is true
#' @param aggregation_method aggregation method
#' @param get_completeRWRmat boolean 
#' @param no_seed_nodes vector of characterss
#' @param cores Number of threads for Parallelization. It has to be positive integer. If it is equal to 1, no parallelization is not performed
#' @return RWRMH_similarity
#' @export
#' 


gen_sim_mat_MH <- function(network1, network2,                                  # outputs of MoNETA::create_multiplex()
                           network_B,                                           # Bipartite network 
                           restart = 0.7, 
                           tau1 = NA, tau2 = NA,                                # restarting probabilities per multiplex 
                           delta1 = 0.5, delta2 = 0.5, lambda = 0.5,            # inter-layers and inter-multiplex probabilities
                           layer_transition_1 = NULL, layer_transition_2 = NULL,
                           jump_neighborhood_1 = FALSE, weighted_multiplex_1 = FALSE,
                           jump_neighborhood_2 = FALSE, weighted_multiplex_2 = FALSE,
                           aggregation_method = 'Sum', get_completeRWRmat = FALSE, no_seed_nodes = NULL,
                           cores = 1){
  
  cond_jump1 = layer_transition_1
  cond_jump2 = layer_transition_2
  
  InputNetwork1 <- network1
  InputNetwork2 <- network2
  
  if (ncol(InputNetwork1) != 4 | ncol(InputNetwork2) != 4){
    stop("Both of the input networks should have 4 columns. See help for format details",
         call. = FALSE)
  } else {
    base::colnames(InputNetwork1) <- c("EdgeType", "source", "target", "weight")
    base::colnames(InputNetwork2) <- c("EdgeType", "source", "target", "weight")
    Number_Layers1 <- length(unique(InputNetwork1$EdgeType))
    Number_Layers2 <- length(unique(InputNetwork2$EdgeType))
    
    print(paste0("Input multiplex networks with ", Number_Layers1, " and ", Number_Layers2, ' Layers respectively.'))
    
    if (!is.numeric(InputNetwork1$weight) | !is.numeric(InputNetwork2$weight)){
      stop("The weights in the input networks should be numeric")
    }
  }
  
  if (restart > 1 || restart < 0){
    stop("Restart parameter should range between 0 and 1")
  }
  
  if (length(tau1) == 1 && is.na(tau1)){
    tau1 <- rep(1, Number_Layers1)/Number_Layers1
  }
  if (length(tau2) == 1 && is.na(tau2)){
    tau2 <- rep(1, Number_Layers2)/Number_Layers2
  }
  print(paste0('Restarting probability per layer in gene Network: ',  paste(tau1, collapse = ', ')))
  print(paste0('Restarting probability per layer in drug Network: ',  paste(tau2, collapse = ', ')))
  
  
  
  ##############################################################################
  #        Data Transformation and associated calculations to apply RWR        #
  ##############################################################################  
  ## We transform the input multiplex format to a L-length list of igraphs objects.
  ## We also scale the weigths of every layer between 1 and the minimun
  ## divided by the maximun weight.
  
  
  Layers_list <- list()
  i = 1
  for (InputNetwork in list(InputNetwork1, InputNetwork2)){
    LayersNames <- unique(InputNetwork$EdgeType)
    Layers <- lapply(LayersNames, function(x) {
      Original_weight <- InputNetwork$weight[InputNetwork$EdgeType == x]
      if (all(!is.na(Original_weight)) && any(Original_weight > 0) && !length(unique(Original_weight)) == 1){
        b <- 1
        a <- min(Original_weight) / max(Original_weight)
        range01 <- (b - a) * (Original_weight - min(Original_weight)) / (max(Original_weight) - min(Original_weight)) + a
        InputNetwork$weight[InputNetwork$EdgeType == x] <- range01
      }else {
        InputNetwork$weight[InputNetwork$EdgeType == x] <- 1
      }
      
      igraph::simplify(igraph::graph_from_data_frame(InputNetwork[InputNetwork$EdgeType == x, 2:4],
                                                     directed = FALSE), edge.attr.comb=mean)
    })
    names(Layers) <- LayersNames
    Layers_list[[i]] <- Layers
    i = i + 1
  }
  
  ## In case the network is monoplex, we have to be aware of isolated nodes.
  IsolatedVertex1 = NULL
  if (Number_Layers1 == 1){
    message("Dealing with isolated nodes in gene Network ...")
    IsolatedVertex1 <- igraph::V(Layers_list[[1]][[1]])$name[which(igraph::degree(Layers_list[[1]][[1]]) == 0)]
    mynetPrev1 <- Layers_list[[1]][[1]]
    Layers_list[[1]][[1]] <- igraph::delete.vertices(igraph::simplify(Layers_list[[1]][[1]]),
                                           igraph::degree(Layers_list[[1]][[1]]) == 0)
  }
  
  IsolatedVertex2 = NULL
  if (Number_Layers2 == 1){
    message("Dealing with isolated nodes in drug Network ...")
    IsolatedVertex2 <- igraph::V(Layers_list[[2]][[1]])$name[which(igraph::degree(Layers_list[[2]][[1]]) == 0)]
    mynetPrev2 <- Layers_list[[2]][[1]]
    Layers_list[[2]][[1]] <- igraph::delete.vertices(igraph::simplify(Layers_list[[2]][[1]]),
                                           igraph::degree(Layers_list[[2]][[1]]) == 0)
  }
  
  
  ## We prepare the data to compute RWR-MH 
  # (Compute the adjacency matrices of both multiplex networks and their normalization)
  
  Multiplex_Object1 <- create.multiplex(Layers_list[[1]])
  Multiplex_Object2 <- create.multiplex(Layers_list[[2]])
  
  MultiplexHet_Object <- create.multiplexHet(Multiplex_Object1, Multiplex_Object2, network_B)
  
  message('Check jumping parameters for omics multiplex network...')
  
  jump_mat_nodes_1 = NULL
  if (jump_neighborhood_1) {
    jump_mat_nodes_1 = list()
    attr = NULL
    if (weighted_multiplex_1) {
      attr = "weight"
    }
    for (mo_name in names(Multiplex_Object1)[seq(Multiplex_Object1$Number_of_Layers)]) {
      mo = Multiplex_Object1[[mo_name]]
      adjacency = as_adjacency_matrix(mo, attr = attr, sparse = T)
      rsum = rowSums(adjacency)
      rsum[rsum == 0] = 1
      wadj = adjacency/(rsum * 2)
      diag(wadj) = 0.5
      wadj[rowSums(wadj) == 0.5 & wadj == 0.5] = 1
      jump_mat_nodes_1[[mo_name]] = wadj
    }
  }
  
  message('Check jumping parameters for drug multiplex network...')
  
  jump_mat_nodes_2 = NULL
  if (jump_neighborhood_2) {
    jump_mat_nodes_2 = list()
    attr = NULL
    if (weighted_multiplex_2) {
      attr = "weight"
    }
    for (mo_name in names(Multiplex_Object2)[seq(Multiplex_Object2$Number_of_Layers)]) {
      mo = Multiplex_Object2[[mo_name]]
      adjacency = as_adjacency_matrix(mo, attr = attr, sparse = T)
      rsum = rowSums(adjacency)
      rsum[rsum == 0] = 1
      wadj = adjacency/(rsum * 2)
      diag(wadj) = 0.5
      wadj[rowSums(wadj) == 0.5 & wadj == 0.5] = 1
      jump_mat_nodes_2[[mo_name]] = wadj
    }
  }
  
  message('Computing the transitiona matrix...')
  trans_matrix <- compute_transition_matrix(MultiplexHet_Object, 
                                            lambda = 0.5, delta1=0.5, delta2=0.5, 
                                            cond_jump1, cond_jump2,
                                            jump_mat_nodes_1, jump_mat_nodes_2 )
                             
  
   
  ##############################################################################
  #   Apply RWR on every node of the multiplex het network (Parallel version)
  ##############################################################################
  
  ## Apply Random Walk with Restart for each node of the Multiplex and Heterogeneous network.

  Allnodes <- c(Multiplex_Object1$Pool_of_Nodes, Multiplex_Object2$Pool_of_Nodes)     
  f_Allnodes <- Allnodes[!(Allnodes %in% no_seed_nodes)]
  numberNodes <- length(Allnodes)
  f_numberNodes <- length(f_Allnodes)
  
  message("Computing RWR for every network node ...")
  
  cl <- snow::makeCluster(cores)
  doSNOW::registerDoSNOW(cl)
  iterations <- length(f_Allnodes)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  Results <- foreach::foreach(i = 1:length(f_Allnodes),
                              .packages = c("Matrix"),
                              .export = c('Random.Walk.Restart.MultiplexHet', 'isMultiplexHet',
                                          'get.seed.scoresMultiplex', "sumValues"),
                              .options.snow=opts) %dopar% {
                                Random.Walk.Restart.MultiplexHet(x = trans_matrix, 
                                                                 MultiplexHet_Object = MultiplexHet_Object,
                                                                 Multiplex1_Multiplex2_Seeds = f_Allnodes[i],
                                                                 r = restart, 
                                                                 tau1, tau2, 
                                                                 DispResults = "Alphabetic", 
                                                                 MeanType = aggregation_method,
                                                                 get_completeRWRmat = get_completeRWRmat)
                              }
  
  close(pb)
  snow::stopCluster(cl)

  
  
  ##############################################################################
  #               Generation of the output RWR_MH Similarity Matrix            #
  ##############################################################################
  
  if(get_completeRWRmat) {
    message("Building pre-aggregation Similarity Matrix ...")
    preRWRMH_similarity <- matrix(data = 0, nrow = nrow(Results[[1]]$preRWRMH_Results), ncol = f_numberNodes)
    for (j in seq(f_numberNodes)) {
      preRWRMH_similarity[,j] <- Results[[j]]$preRWRMH_Results$Score
    }
    base::rownames(preRWRMH_similarity) <- Results[[1]]$preRWRMH_Results$NodeNames
    base::colnames(preRWRMH_similarity) <- f_Allnodes
  }else{ 
    preRWRMH_similarity <- NULL 
  }
  
  message("Building Similarity Matrix ...")
  RWRMH_similarity <- matrix(data = 0, nrow = numberNodes, ncol = f_numberNodes)
  for (j in seq(f_numberNodes)) {
    RWRMH_similarity[,j] <- Results[[j]]$RWRMH_Results$Score
  }
  base::rownames(RWRMH_similarity) <- Allnodes
  base::colnames(RWRMH_similarity) <- f_Allnodes
  
  
  if (Number_Layers1 == 1){
    numberIso <- length(IsolatedVertex1)
    if (numberIso > 0){
      print("Including isolated nodes from the first Network in the Similarity Matrix ...")
      TotalSizeMatrix <-   numberNodes + numberIso
      Index_IsoNodes <- IsolatedVertex1
      
      NewMatrix <- matrix(data = 0, nrow = TotalSizeMatrix, ncol = TotalSizeMatrix)
      base::rownames(NewMatrix) <- c(base::rownames(RWRMH_similarity), Index_IsoNodes)
      base::colnames(NewMatrix) <- c(base::colnames(RWRMH_similarity), Index_IsoNodes)
      
      NewMatrix[1:numberNodes, 1:numberNodes] <- RWRMH_similarity
      NewMatrix[(numberNodes + 1):TotalSizeMatrix, (numberNodes + 1):TotalSizeMatrix] <- diag(1, numberIso, numberIso)
      
      RWRMH_similarity <- NewMatrix
    }
  }
    
  if (Number_Layers2 == 1){
    numberIso <- length(IsolatedVertex2)
    if (numberIso > 0){
      print("Including isolated nodes from the second Network in the Similarity Matrix ...")
      TotalSizeMatrix <-  numberNodes + numberIso
      Index_IsoNodes <- IsolatedVertex2
      
      NewMatrix <- matrix(data = 0, nrow = TotalSizeMatrix, ncol = TotalSizeMatrix)
      base::rownames(NewMatrix) <- c(base::rownames(RWRMH_similarity), Index_IsoNodes)
      base::colnames(NewMatrix) <- c(base::colnames(RWRMH_similarity), Index_IsoNodes)
      
      NewMatrix[1:numberNodes, 1:numberNodes] <- RWRMH_similarity
      NewMatrix[(numberNodes + 1):TotalSizeMatrix, (numberNodes + 1):TotalSizeMatrix] <- diag(1, numberIso, numberIso)
      
      RWRMH_similarity <- NewMatrix
    }
  }
  
  
  if (is.numeric(base::rownames(RWRMH_similarity))){
    RWRMH_similarity <- RWRMH_similarity[order(as.numeric(base::rownames(RWRMH_similarity))),
                                       order(as.numeric(base::colnames(RWRMH_similarity)))]
  } else {
    sorted_genes <- sort(c(Multiplex_Object1$Pool_of_Nodes, IsolatedVertex1))
    sorted_drugs <- sort(c(Multiplex_Object2$Pool_of_Nodes, IsolatedVertex2))
    
    f_sorted_genes <- sorted_genes[!sorted_genes %in% no_seed_nodes]
    f_sorted_drugs <- sorted_drugs[!sorted_drugs %in% no_seed_nodes]
    
    RWRMH_similarity <- RWRMH_similarity[c(sorted_genes, sorted_drugs),
                                         c(f_sorted_genes, f_sorted_drugs)]
  }

  return(list(RWRMH_sim_mat = RWRMH_similarity, 
              wholeRWRMH_sim_mat = preRWRMH_similarity))
}


################################################################################



isMultiplex <- function (x)
{
  is(x, "Multiplex")
}

isMultiplexHet <- function (x)
{
  is(x, "MultiplexHet")
}




## Add missing nodes in some of the layers.
add.missing.nodes <- function (Layers,Nr_Layers,NodeNames) {
  
  add_vertices(Layers,
               length(NodeNames[which(!NodeNames %in% V(Layers)$name)]),
               name=NodeNames[which(!NodeNames %in%  V(Layers)$name)])
}


# create.multiplex <- function(...){
#   UseMethod("create.multiplex")
# }

#create.multiplex.default <- function(LayersList,...)
create.multiplex <- function(LayersList,...)
{
  
  Number_of_Layers <- length(LayersList)
  SeqLayers <- seq(Number_of_Layers)
  
  if (!all(sapply(SeqLayers, function(x) is.igraph(LayersList[[x]])))){
    stop("Not igraph objects")
  }
  
  ## We get a pool of nodes (Nodes in any of the layers.)
  Pool_of_Nodes <-
    sort(unique(unlist(lapply(SeqLayers, function(x) V(LayersList[[x]])$name))))
  
  if (is.numeric(Pool_of_Nodes)){
    Pool_of_Nodes <- sort(as.numeric(Pool_of_Nodes))
  } else {
    Pool_of_Nodes <- sort(Pool_of_Nodes)
  }
  
  Number_of_Nodes <- length(Pool_of_Nodes)
  
  Layer_List <-
    lapply(LayersList, add.missing.nodes,Number_of_Layers,Pool_of_Nodes)
  
  
  MultiplexObject <- c(Layer_List,list(Pool_of_Nodes=Pool_of_Nodes,
                                       Number_of_Nodes_Multiplex=Number_of_Nodes,
                                       Number_of_Layers=Number_of_Layers))
  
  class(MultiplexObject) <- "Multiplex"
  
  return(MultiplexObject)
}


# create.multiplexHet <- function(...) {
#   UseMethod("create.multiplexHet")
# }


#create.multiplexHet.default  <- function(MultiObject1, MultiObject2,
create.multiplexHet  <- function(MultiObject1, MultiObject2,
                                 BipartiteNetwork,...)
{
  
  Allnodes1 <- MultiObject1$Pool_of_Nodes
  Allnodes2 <- MultiObject2$Pool_of_Nodes
  
  bipartiteNodesNetwork1 <- unique(c(as.character(BipartiteNetwork$source)))
  bipartiteNodesNetwork2 <- unique(c(as.character(BipartiteNetwork$target)))
  
  if (!all(bipartiteNodesNetwork1 %in% Allnodes1)){
    stop("Some of the source nodes of the bipartite are not present in
         the first network")
  }
  
  if (!all(bipartiteNodesNetwork2 %in% Allnodes2)){
    stop("Some of the target nodes of the bipartite are not present in
         the second network")
  }
  
  ## Multiplex graph features
  NumberNodes1 <- MultiObject1$Number_of_Nodes
  NumberLayer1 <- MultiObject1$Number_of_Layers
  
  NumberNodes2 <- MultiObject2$Number_of_Nodes
  NumberLayer2 <- MultiObject2$Number_of_Layers
  
  message("Generating bipartite matrix ...")
  Bipartite_Matrix <-
    get.bipartite.graph(Allnodes1,Allnodes2,BipartiteNetwork,NumberNodes1,
                        NumberNodes2)
  
  message("Expanding bipartite matrix to fit the multiplex network ...")
  Supra_Bipartite_Matrix <- expand.bipartite.graph(NumberNodes1,NumberLayer1,
                                                   NumberNodes2,NumberLayer2,Bipartite_Matrix)
  
  Multiplex_HetObject <- list(Multiplex1 = MultiObject1,
                              Multiplex2 = MultiObject2,
                              BipartiteNetwork = Supra_Bipartite_Matrix)
  
  class(Multiplex_HetObject) <- "MultiplexHet"
  return(Multiplex_HetObject)
}



## Bipartite graph construction.
get.bipartite.graph <- function(Names_Mul1, Names_Mul2, BipartiteNetwork,
                                Number_Nodes_1,Number_Nodes_2){
  
  Bipartite_matrix <- Matrix(data=0, nrow=Number_Nodes_1, ncol=Number_Nodes_2)
  Names_Mul1_order <- sort(Names_Mul1)
  Names_Mul2_order <- sort(Names_Mul2)
  rownames(Bipartite_matrix) <- Names_Mul1_order
  colnames(Bipartite_matrix) <- Names_Mul2_order
  
  for (i in seq(nrow(BipartiteNetwork))){
    Bipartite_matrix[BipartiteNetwork$source[i], BipartiteNetwork$target[i]] <-
      BipartiteNetwork$weight[i]
  }
  
  return(Bipartite_matrix)
}


## Fitting the bipartite graph to the multiplex networks.
expand.bipartite.graph <-
  function(Number_Nodes_1,Number_Layers_1,Number_Nodes_2,
           Number_Layers_2,Bipartite_matrix){
    
    Supra_Bipartite_Matrix <-
      do.call(rbind, replicate(Number_Layers_1,Bipartite_matrix,simplify=FALSE))
    
    rownames(Supra_Bipartite_Matrix) <-
      paste0(rownames(Bipartite_matrix), sep="_",rep(seq(Number_Layers_1),
                                                     each=Number_Nodes_1))
    
    
    Supra_Bipartite_Matrix <-
      do.call(cbind, replicate(Number_Layers_2,Supra_Bipartite_Matrix,
                               simplify=FALSE))
    
    colnames(Supra_Bipartite_Matrix) <-
      paste0(colnames(Bipartite_matrix), sep="_",rep(seq(Number_Layers_2),
                                                     each=Number_Nodes_2))
    
    return(Supra_Bipartite_Matrix)
  }



### Transition functions
get.transition.multiplex1.multiplex2 <-
  function(Number_Nodes_Multiplex1, Number_Layers1,Number_Nodes_Multiplex2,
           Number_Layers2, SupraBipartiteMatrix,lambda){
    
    TransitionMat_Multiplex1_Multiplex2 <-
      Matrix(0, nrow=Number_Nodes_Multiplex1*Number_Layers1,
             ncol=Number_Nodes_Multiplex2*Number_Layers2,sparse = TRUE)
    
    colnames(TransitionMat_Multiplex1_Multiplex2) <-
      colnames(SupraBipartiteMatrix)
    rownames(TransitionMat_Multiplex1_Multiplex2) <-
      rownames(SupraBipartiteMatrix)
    
    Col_Sum_Bipartite <-
      Matrix::colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
                       sparseResult = FALSE)
    
    m <- lambda * t(t(SupraBipartiteMatrix) / Col_Sum_Bipartite)
    idx <- Col_Sum_Bipartite != 0
    TransitionMat_Multiplex1_Multiplex2[,idx] = m[,idx]
    
    return(TransitionMat_Multiplex1_Multiplex2)
  }


get.transition.multiplex2.multiplex1 <-
  function(Number_Nodes_Multiplex1, Number_Layers1,Number_Nodes_Multiplex2,
           Number_Layers2,SupraBipartiteMatrix,lambda){
    
    TransitionMat_Multiplex2_Multiplex1 <-
      Matrix(0,nrow=Number_Nodes_Multiplex2*Number_Layers2,
             ncol=Number_Nodes_Multiplex1*Number_Layers1,sparse = TRUE)
    
    colnames(TransitionMat_Multiplex2_Multiplex1) <-
      rownames(SupraBipartiteMatrix)
    rownames(TransitionMat_Multiplex2_Multiplex1) <-
      colnames(SupraBipartiteMatrix)
    
    Row_Sum_Bipartite <-
      Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
                       sparseResult = FALSE)
    
    m <- lambda * t((SupraBipartiteMatrix) / Row_Sum_Bipartite)
    idx <- Row_Sum_Bipartite != 0
    TransitionMat_Multiplex2_Multiplex1[,idx] = m[,idx]
    
    return(TransitionMat_Multiplex2_Multiplex1)
  }


get.transition.multiplex <-
  function(Number_Nodes,Number_Layers, lambda,SupraAdjacencyMatrix,
           SupraBipartiteMatrix) {
    
    Transition_Multiplex_Network <-
      Matrix(0, nrow=Number_Nodes*Number_Layers,
             ncol=Number_Nodes*Number_Layers,sparse = TRUE)
    
    rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
    colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
    
    Col_Sum_Multiplex <-
      Matrix::colSums(SupraAdjacencyMatrix,na.rm=FALSE, dims=1,
                      sparseResult=FALSE)
    Row_Sum_Bipartite <-
      Matrix::rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,
                       sparseResult = FALSE)
    
    idx <- Row_Sum_Bipartite != 0
    Transition_Multiplex_Network[,idx] <-
      ((1-lambda)*t(t(SupraAdjacencyMatrix[,idx])/Col_Sum_Multiplex[idx]))
    
    Transition_Multiplex_Network[,!idx] <-
      t(t(SupraAdjacencyMatrix[,!idx]) / Col_Sum_Multiplex[!idx])
    
    return(Transition_Multiplex_Network)
  }


compute_transition_matrix <- function(x, lambda = 0.5, delta1=0.5, delta2=0.5, 
                                      cond_jump1 = NULL, cond_jump2 = NULL,
                                      jump_mat_nodes_1 = NULL, jump_mat_nodes_2 = NULL)
{
  if (!isMultiplexHet(x)) {
    stop("Not a Multiplex Heterogeneous object")
  }
  
  if (delta1 > 1 || delta1 <= 0) {
    stop("Delta should be between 0 and 1")
  }
  
  if (delta2 > 1 || delta2 <= 0) {
    stop("Delta should be between 0 and 1")
  }
  
  if (lambda > 1 || lambda <= 0) {
    stop("Lambda should be between 0 and 1")
  }
  
  NumberNodes1 <- x$Multiplex1$Number_of_Nodes_Multiplex
  NumberLayers1 <- x$Multiplex1$Number_of_Layers
  
  NumberNodes2 <- x$Multiplex2$Number_of_Nodes_Multiplex
  NumberLayers2 <- x$Multiplex2$Number_of_Layers
  
  SupraBipartiteMatrix <- x$BipartiteNetwork
  
  message("Computing adjacency matrix of the first input network ...")
  AdjMatrix_Multiplex1 <- compute_adjacency_matrix2(x$Multiplex1, delta1, cond_jump1, jump_mat_nodes_1)
  norm_AdjMatrix_Multiplex1 <- normalize.multiplex.adjacency(AdjMatrix_Multiplex1)
  
  message("Computing adjacency matrix of the second input network ...")
  AdjMatrix_Multiplex2 <- compute_adjacency_matrix2(x$Multiplex2, delta2, cond_jump2, jump_mat_nodes_2)
  norm_AdjMatrix_Multiplex2 <- normalize.multiplex.adjacency(AdjMatrix_Multiplex2)
  
  ## Transition Matrix for the inter-subnetworks links
  message("Computing inter-subnetworks transitions ...")
  Transition_Multiplex1_Multiplex2 <-
    get.transition.multiplex1.multiplex2(NumberNodes1,NumberLayers1,
                                         NumberNodes2,NumberLayers2,SupraBipartiteMatrix,lambda)
  
  Transition_Multiplex2_Multiplex1 <-
    get.transition.multiplex2.multiplex1(NumberNodes1,NumberLayers1,
                                         NumberNodes2, NumberLayers2, SupraBipartiteMatrix,lambda)
  
  ## Transition Matrix for the intra-subnetworks links
  message("Computing intra-subnetworks transitions ...")
  Transition_Multiplex_Network1 <-
    get.transition.multiplex(NumberNodes1, NumberLayers1, lambda,
                             norm_AdjMatrix_Multiplex1,
                             SupraBipartiteMatrix)
  Transition_Multiplex_Network2 <-
    get.transition.multiplex(NumberNodes2, NumberLayers2, lambda,
                             t(norm_AdjMatrix_Multiplex2),t(SupraBipartiteMatrix))
  
  ## We generate the global transition matrix and we return it.
  message("Combining inter e intra layer probabilities into the global Transition Matrix")
  Transition_Multiplex_Heterogeneous_Matrix_1 <-
    cbind(Transition_Multiplex_Network1, Transition_Multiplex1_Multiplex2)
  Transition_Multiplex_Heterogeneous_Matrix_2 <-
    cbind(Transition_Multiplex2_Multiplex1, Transition_Multiplex_Network2)
  Transition_Multiplex_Heterogeneous_Matrix <-
    rbind(Transition_Multiplex_Heterogeneous_Matrix_1,
          Transition_Multiplex_Heterogeneous_Matrix_2)
  
  return(Transition_Multiplex_Heterogeneous_Matrix)
}



##### Compute adjacency matrix (MoNETA version) + Normalization

compute_adjacency_matrix2 <- function(x, delta = 0.5, cond_jump = NULL, jump_mat_nodes = NULL)
{
  if (!isMultiplex(x)) {
    stop("Not a Multiplex or Multiplex Heterogeneous object")
  }
  if (delta > 1 || delta < 0) {
    stop("Delta should be between 0 and 1")
  }
  
  
  N <- x$Number_of_Nodes_Multiplex
  L <- x$Number_of_Layers
  
  ###############################################################################
  if (!is.null(cond_jump) && (any(base::colnames(cond_jump) != names(x)[seq(L)]) || any(base::rownames(cond_jump) != names(x)[seq(L)]))) {
    stop("cond_jump elements order differs from x")
  }
  ###############################################################################
  ## We impose delta=0 in the monoplex case.
  if (L  ==1){
    delta = 0
  }
  
  Layers_Names <- names(x)[seq(L)]
  
  ###############################################################################
  if (is.null(cond_jump))
    cond_jump <- matrix(1 / (L - 1), nrow = L, ncol = L, dimnames = list(Layers_Names, Layers_Names))
  ###############################################################################
  
  diag(cond_jump) <- 0
  ## IDEM_MATRIX.
  if (is.null(jump_mat_nodes)) {
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)
  }
  
  counter <- 0
  Layers_List <- lapply(x[Layers_Names], function(x){
    
    counter <<- counter + 1;
    if (igraph::is_weighted(x)) {
      Adjacency_Layer <- igraph::as_adjacency_matrix(x, sparse = TRUE,
                                                     attr = "weight")
    } else {
      Adjacency_Layer <- igraph::as_adjacency_matrix(x, sparse = TRUE)
    }
    
    if (is.numeric(base::rownames(Adjacency_Layer))){
      Adjacency_Layer <- Adjacency_Layer[order(as.numeric(base::rownames(Adjacency_Layer))),
                                         order(as.numeric(base::colnames(Adjacency_Layer)))]
    } else {
      Adjacency_Layer <- Adjacency_Layer[order(base::rownames(Adjacency_Layer)),
                                         order(base::colnames(Adjacency_Layer))]
    }
    
    base::colnames(Adjacency_Layer) <-
      paste0(base::colnames(Adjacency_Layer), "_", counter)
    base::rownames(Adjacency_Layer) <-
      paste0(base::rownames(Adjacency_Layer), "_", counter)
    Adjacency_Layer
  })
  
  MyColNames <- unlist(lapply(Layers_List, function(x) unlist(base::colnames(x))))
  MyRowNames <- unlist(lapply(Layers_List, function(x) unlist(base::rownames(x))))
  names(MyColNames) <- c()
  names(MyRowNames) <- c()
  SupraAdjacencyMatrix <- (1 - delta) * (Matrix::bdiag(unlist(Layers_List)))
  base::colnames(SupraAdjacencyMatrix) <- MyColNames
  base::rownames(SupraAdjacencyMatrix) <- MyRowNames
  
  #offdiag <- (delta/(L-1))*Idem_Matrix
  
  i <- seq_len(L)
  Position_ini_row <- 1 + (i - 1) * N
  Position_end_row <- N + (i - 1) * N
  j <- seq_len(L)
  Position_ini_col <- 1 + (j - 1) * N
  Position_end_col <- N + (j - 1) * N
  
  for (i in seq_len(L)){
    if (is.null(jump_mat_nodes)) {
      mat_tmp = Idem_Matrix
    } else {
      mat_tmp = jump_mat_nodes[[i]]
    }
    for (j in seq_len(L)){
      if (j != i){
        
        SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                             (Position_ini_col[j]:Position_end_col[j])] <- delta * cond_jump[i,j] * mat_tmp
      }
    }
  }
  
  SupraAdjacencyMatrix <- methods::as(SupraAdjacencyMatrix, "dgCMatrix")
  return(SupraAdjacencyMatrix)
}



normalize.multiplex.adjacency <- function(x)
{
  if (!is(x,"dgCMatrix")){
    stop("Not a dgCMatrix object of Matrix package")
  }
  
  Adj_Matrix_Norm <- t(t(x)/(Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                             sparseResult = FALSE)))
  
  return(Adj_Matrix_Norm)
}



#### Compute seed scores

get.seed.scoresMultiplex <- function(Seeds, Number_Layers, tau) {
  
  Nr_Seeds <- length(Seeds)
  
  Seeds_Seeds_Scores <- rep(tau / Nr_Seeds, Nr_Seeds)
  Seed_Seeds_Layer_Labeled <-
    paste0(rep(Seeds, Number_Layers), sep = "_", rep(seq(Number_Layers),
                                                     length.out = Nr_Seeds * Number_Layers, each = Nr_Seeds))
  
  Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                            Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)
  
  return(Seeds_Score)
}



###### Get final seed scores

geometric.mean <- function(Scores, L, N) {
  FinalScore <- numeric(length = N)
  for (i in seq_len(N)){
    FinalScore[i] <- prod(Scores[seq(from = i, to = N * L, by = N)]) ^ (1 / L)
  }
  return(FinalScore)
}


regular.mean <- function(Scores, L, N) {
  FinalScore <- numeric(length = N)
  for (i in seq_len(N)){
    FinalScore[i] <- mean(Scores[seq(from = i, to = N * L, by = N)])
  }
  return(FinalScore)
}


sumValues <- function(Scores, L, N) {
  FinalScore <- numeric(length = N)
  for (i in seq_len(N)){
    FinalScore[i] <- sum(Scores[seq(from = i, to = N * L, by = N)])
  }
  return(FinalScore)
}




###### RWR-MH

Random.Walk.Restart.MultiplexHet <-
  function(x, MultiplexHet_Object, 
           Multiplex1_Multiplex2_Seeds,
           r = 0.7, tau1, tau2, 
           MeanType= aggregation_method,
           get_completeRWRmat = get_completeRWRmat, 
           DispResults="Alphabetic", ...){
    
    ## We control the different values.
    if (!"dgCMatrix" %in% class(x)){
      stop("Not a dgCMatrix object of Matrix package")
    }
    
    if (!isMultiplexHet(MultiplexHet_Object)) {
      stop("Not a Multiplex Heterogeneous object")
    }
    
    NumberLayers1 <- MultiplexHet_Object$Multiplex1$Number_of_Layers
    NumberNodes1 <- MultiplexHet_Object$Multiplex1$Number_of_Nodes_Multiplex
    NumberLayers2 <- MultiplexHet_Object$Multiplex2$Number_of_Layers
    NumberNodes2 <- MultiplexHet_Object$Multiplex2$Number_of_Nodes_Multiplex
    
    All_nodes_Multiplex1 <- MultiplexHet_Object$Multiplex1$Pool_of_Nodes
    All_nodes_Multiplex2 <- MultiplexHet_Object$Multiplex2$Pool_of_Nodes
    
    
    if (length(Multiplex1_Multiplex2_Seeds) < 1) {
      stop("You did not provided any seeds")
    } else {
      if (length(Multiplex1_Multiplex2_Seeds) >= NumberNodes1 | length(Multiplex1_Multiplex2_Seeds) >= NumberNodes2){
        stop("The length of some of the vectors containing the seed nodes
             is not correct")
      }  else {
        if (!all(Multiplex1_Multiplex2_Seeds %in% c(All_nodes_Multiplex1, All_nodes_Multiplex2))){
          stop("Some of the  input seeds are not nodes of the first
               input network")
        }
      }
    }
    
    if (r >= 1 || r <= 0) {
      stop("Restart parameter should be between 0 and 1")
    }
    
    
    if(!(as.character(MeanType) %in% c("Geometric","Arithmetic","Sum"))){
      stop("The type mean should be Geometric, Arithmetic or Sum")
    }
    
    if(!(DispResults %in% c("TopScores","Alphabetic"))){
      stop("The way to display RWRM results should be TopScores or
           Alphabetic")
    }
    
    ## We define the threshold and the number maximum of iterations
    ## for the random walker.
    Threeshold <- 1e-10
    NetworkSize <- ncol(x)
    
    ## We initialize the variables to control the flux in the RW algo.
    residue <- 1
    iter <- 1
    
    ## We compute the scores for the different seeds.
    if (Multiplex1_Multiplex2_Seeds %in% All_nodes_Multiplex1){
      NumberLayers = NumberLayers1
      tau = tau1
    } else {
      NumberLayers = NumberLayers2
      tau = tau2
    }
    
    Seeds_Score <-
      get.seed.scoresMultiplex(Multiplex1_Multiplex2_Seeds, NumberLayers, 
                               tau)
    
    ## We define the prox_vector(The vector we will move after the first
    ## RWR iteration. We start from The seed. We have to take in account
    ## that the walker with restart in some of the Seed genes,
    ## depending on the score we gave in that file).
    prox_vector <- matrix(0, nrow = NetworkSize, ncol=1)
    
    prox_vector[which(colnames(x) %in% Seeds_Score[,1])] <- (Seeds_Score[,2])
    
    prox_vector  <- prox_vector/sum(prox_vector)
    restart_vector <-  prox_vector
    
    while(residue >= Threeshold){
      
      old_prox_vector <- prox_vector
      prox_vector <- (1-r)*(x %*% prox_vector) + r*restart_vector
      residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
      iter <- iter + 1;
    }
    
    IndexSep <- NumberNodes1*NumberLayers1
    prox_vector_1 <- prox_vector[1:IndexSep,]
    prox_vector_2 <- prox_vector[(IndexSep+1):nrow(prox_vector),]
    
    NodeNames1 <-
      gsub("_1", "",names(prox_vector_1)[seq_len(NumberNodes1)])
    NodeNames2 <-
      gsub("_1", "",names(prox_vector_2)[seq_len(NumberNodes2)])
    
    
    if (get_completeRWRmat){
      preAggregated_df <- data.frame(NodeNames = c(names(prox_vector_1), names(prox_vector_2)),
                                     Score = c(prox_vector_1, prox_vector_2))
      rownames(preAggregated_df) <- c()
    } else{
      preAggregated_df <- NULL
    }
    
    
    if (MeanType=="Geometric"){
      rank_global1 <- geometric.mean(prox_vector_1,NumberLayers1,NumberNodes1)
      rank_global2 <- geometric.mean(prox_vector_2,NumberLayers2,NumberNodes2)
    } else {
      if (MeanType=="Arithmetic") {
        rank_global1 <- regular.mean(prox_vector_1,NumberLayers1,NumberNodes1)
        rank_global2 <- regular.mean(prox_vector_2,NumberLayers2,NumberNodes2)
      } else {
        rank_global1 <- sumValues(prox_vector_1,NumberLayers1,NumberNodes1)
        rank_global2 <- sumValues(prox_vector_2,NumberLayers2,NumberNodes2)
      }
    }
    
    Global_results <- data.frame(NodeNames = c(NodeNames1,NodeNames2),
                                 Score = c(rank_global1,rank_global2))
    
    if (DispResults=="TopScores"){
      ## We sort the nodes according to their score.
      Global_results <-
        Global_results[with(Global_results, order(-Score, NodeNames)), ]
      
      ### We remove the seed nodes from the Ranking and we write the results.
      Global_results <-
        Global_results[which(!Global_results$NodeNames %in% Seeds),]
    } else {
      Global_results <- Global_results
    }
    
    rownames(Global_results) <- c()
    
    RWRMH_ranking <- list(Seed_Nodes = c(Multiplex1_Multiplex2_Seeds),
                          RWRMH_Results = Global_results,
                          preRWRMH_Results = preAggregated_df)
    
    return(RWRMH_ranking)
  }

