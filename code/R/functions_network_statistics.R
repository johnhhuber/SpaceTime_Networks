# mc3TransmissionNetworkStats
networkStatistic <- function(chain.file, burnIn=0, numLines=-1L, node.file=NULL, verbose=F){
  require(tools)
  require(igraph)
  
  if(verbose){ message("File: ", chain.file) }
  # set up stream connetion
  extFn <- file_ext(chain.file)
  if(extFn=="bz2")
    chain <- bzfile(chain.file, open="rb")
  else if(extFn=="gz")
    chain <- gzfile(chain.file, open="rb") 
  else if(extFn=="zip")
    chain <- unz(chain.file, open="rb")  
  else
    chain <- file(chain.file, open="rt")
  
  # read in burnin samples
  if(burnIn>1){ startChunk <- readLines(chain, n=burnIn) }
  
  # read in posterior network samples
  numYield <- 1
  resLst <- list()
  while(length(chunk <- readLines(chain, n=numLines))>0){
    if(verbose){ message("Yield: ", numYield) }
    edgesList <- strsplit(chunk, ";")
    numNodes = as.numeric(strsplit(tail(edgesList[[1]],1),'-')[[1]][3])
    
    # set up array to store networks
    networks = matrix(0, nrow=numNodes+1, ncol=numNodes+1)
    ancestor = descendant = character()
    k <- numeric()
    k.distribution = matrix(0, length(chunk), numNodes)
    R0.distribution = matrix(0, nrow=length(chunk), ncol=numNodes)
    
    
    # populate network array
    for(ii in seq_along(chunk)){ # over all samples
      edges <- do.call(rbind, strsplit(as.character(edgesList[[ii]]), '-'))
      edges[edges[,1] == 's',1] <- 0
      ancestor <- as.integer(edges[,1]) + 1
      descendant <- as.integer(edges[,3]) + 1
      k <- as.integer(edges[,2])
      
      networksNew = matrix(0, nrow=numNodes+1, ncol=numNodes+1)
      networksNew[cbind(ancestor, descendant)] = 1
      
      # calculate R0
      R0 <- table(ancestor)[descendant]
      R0[is.na(R0)] <- 0
      names(R0) <- descendant
      
      networks <- networks + networksNew
      k.distribution[ii,seq_along(descendant)] = k
      R0.distribution[ii,] = R0 # TODO take mean of R0
    }
    
    # return a list with posterior samples of network metrics
    resLst[[numYield]] <- list(edgeProbs = networks/length(chunk),
                               k.distribution = k.distribution,
                               R0.distribution = R0.distribution)
    numYield <- numYield+1
  }
  close(chain)
  if(verbose){ message("Done.") }
  
  if(numLines=="-1L")
    return(resLst[[1]])
  else
    return(resLst)
  
}

# calculates the network accuracy
# returns a list with accuracy measure for each network
networkAccuracy <- function(chain.file, trueNetworkFile, burnIn=0, numLines=-1L,  nodeFile=NULL, fullNetworkFile=NULL, fullNodeFile=NULL, verbose=F){
  require(tools)
  
  if(verbose){ message("File: ", chain.file) }
  # set up stream connetion
  extFn <- file_ext(chain.file)
  if(extFn=="bz2")
    chain <- bzfile(chain.file, open="rb")
  else if(extFn=="gz")
    chain <- gzfile(chain.file, open="rb") 
  else if(extFn=="zip")
    chain <- unz(chain.file, open="rb")  
  else
    chain <- file(chain.file, open="rt")
  
  # load true network
  chunk <- readLines(trueNetworkFile, n=numLines)
  edgesList <- strsplit(chunk, ";")
  numNodes = as.numeric(strsplit(tail(edgesList[[1]],1),'-')[[1]][3])
  mat.true = matrix(0, nrow=numNodes+1, ncol=numNodes+1)
  edges <- do.call(rbind, strsplit(as.character(edgesList[[1]]), '-'))
  edges[edges[,1] == 's',1] <- 0
  ancestor <- as.integer(edges[,1]) + 1
  k.true <- as.integer(edges[,2])
  descendant <- as.integer(edges[,3]) + 1
  mat.true[cbind(ancestor, descendant)] = 1
  
  # calculate true R0
  r0.true <- table(ancestor)[descendant]
  r0.true[is.na(r0.true)] <- 0
  names(r0.true) <- descendant
  
  # calculate the true outbreak membership of each node
  adj.true = graph.adjacency(mat.true[-1,-1], mode = 'directed')
  membership.true = components(adj.true)$membership
  
  # identify the true founders and true locally acquired nodes
  founders.true = which(mat.true[1,] == 1) - 1
  locals.true = which(mat.true[1,] == 0)[-1] - 1
  parents.true = sapply(locals.true, function(x){ return(which(mat.true[,x+1] == 1) - 1) })
  which.founder = sapply(locals.true, function(x){ return(which(membership.true == membership.true[x])[1]) })
  ancestor.true = lapply(1:length(locals.true), function(x){
    path <- all_simple_paths(adj.true, which.founder[x], locals.true[x])[[1]]
    return(head(as.numeric(path), n=-1))
  })
  
  # read in burnin network samples
  if(burnIn>1){ startChunk <- readLines(chain, n=burnIn) }
  
  # read in posterior network samples
  numYield <- 1
  resLst <- list()
  while(length(chunk <- readLines(chain, n=numLines))>0){
    if(verbose){ message("Yield: ", numYield) }
    edgesList <- strsplit(chunk, ";")
    numNodes = as.numeric(strsplit(tail(edgesList[[1]],1),'-')[[1]][3])
    
    # set up array to store networks
    founder.acc = numeric(0)
    edge.acc = numeric(0)
    ancestor.acc = numeric(0)
    outbreak.acc = numeric(0)
    k.acc = numeric(0)
    r0.acc = numeric(0)
    
    # populate network array
    for(ii in seq_along(chunk)){ # over all samples
      # creat adjacency matrix
      edges <- do.call(rbind, strsplit(as.character(edgesList[[ii]]), '-'))
      edges[edges[,1] == 's',1] <- 0
      ancestor <- as.integer(edges[,1]) + 1
      k.inferred <- as.integer(edges[,2])
      descendant <- as.integer(edges[,3]) + 1
      networks = matrix(0, nrow=numNodes+1, ncol=numNodes+1)
      networks[cbind(ancestor, descendant)] = 1
      
      # calculate inferred R0
      r0.inferred <- table(ancestor)[descendant]
      r0.inferred[is.na(r0.inferred)] <- 0
      names(r0.inferred) <- descendant
      
      adj.inferred = graph.adjacency(networks[-1,-1], mode = 'directed')
      membership.inferred = components(adj.inferred)$membership
      parents.inferred = sapply(locals.true, function(x){return(which(networks[,x+1] == 1) - 1)})
      
      founder.acc = c(founder.acc, 
                      sum(networks[1,] == mat.true[1,]) / length(networks[1,]))
      edge.acc = c(edge.acc, 
                   sum(parents.true == parents.inferred) / length(parents.true))
      ancestor.acc = c(ancestor.acc,
                       sum(sapply(1:length(parents.inferred), function(x){ is.element(parents.inferred[x], ancestor.true[[x]])})) / length(parents.inferred))
      outbreak.acc = c(outbreak.acc,
                       sum(membership.true[parents.true[which(parents.inferred > 0)]] == membership.true[parents.inferred[which(parents.inferred > 0)]]) / length(parents.inferred))
      k.acc <- c(k.acc, sum(k.true == k.inferred)/length(k.true))
      r0.acc <- c(r0.acc, sum(r0.true == r0.inferred)/length(r0.true))
    }
    #browser()
    # return a list with posterior samples of network metrics
    resLst[[numYield]] <- list(founder.acc = founder.acc,
                               edge.acc = edge.acc,
                               ancestor.acc = ancestor.acc,
                               outbreak.acc = outbreak.acc,
                               k.acc = k.acc,
                               r0.acc = r0.acc)
    numYield <- numYield+1
  }
  close(chain)
  if(verbose){ message("Done.") }
  
  if(numLines=="-1L")
    return(resLst[[1]])
  else{
    return(list(founder.acc = do.call(c, lapply(resLst,"[[", 1)), 
                parent.acc = do.call(c, lapply(resLst,"[[", 2)),
                ancestor.acc = do.call(c, lapply(resLst,"[[", 3)),
                outbreak.acc = do.call(c, lapply(resLst,"[[", 4)),
                k.acc = do.call(c, lapply(resLst,"[[", 5)),
                r0.acc = do.call(c, lapply(resLst,"[[", 6))))
  }
  
}
