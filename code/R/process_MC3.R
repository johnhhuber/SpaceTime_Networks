# if(!require(igraph)){install.packages('igraph'); library(igraph)}
#library(igraph)
#if(!require(igraph)){install.packages('igraph', lib = '~/myRlibs', repos='http://cran.us.r-project.org'); library('igraph', lib = '~/myRlibs')}
require(igraph)

networkstats = function(start, path.in, chain.file)
{
  # read in posterior network samples
  chain = readLines(paste(path.in, chain.file, sep = ''))
  chain = strsplit(chain,';')
  numNodes = as.numeric(strsplit(tail(chain[[1]],1),'-')[[1]][3])

  # set up array to store networks
  networks = matrix(0, nrow = numNodes + 1, ncol = numNodes + 1)
  ancestor = descendant = character()
  k = ii.in = jj.in = numeric()
  k.distribution = matrix(0,length(chain)-start+1,length(chain[[1]]))

  # populate network array
  for(ii in start : length(chain)){
    for(jj in 1 : length(chain[[ii]])){
      ancestor = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][1]
      k = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][2]
      descendant = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][3]

      if(ancestor == 's'){ii.in = 1}
      if(ancestor != 's'){ii.in = as.integer(ancestor) + 1}
      jj.in = as.integer(descendant) + 1

      networks[ii.in, jj.in] = networks[ii.in, jj.in] + 1
      k.distribution[ii-start+1,jj] = as.numeric(k)
    }
  }

  # return a list with posterior samples of network metrics
  return(list(
    edgeProbs = networks / (length(chain)-start+1),
    k.distribution = k.distribution))
}

# generate path in
path.in = commandArgs(trailingOnly = TRUE)
print(path.in)
path.out = path.in

# compute network statistics on posterior networks
files = list.files(path.in)
files = files[which(grepl('networks',files))]
network.data = list()
for(ff in 1:length(files)){
  chain.length =
    as.numeric(
      strsplit(
        system(paste(
          'wc -l ', path.in, files[ff], sep=''), intern=T), split = ' ')[[1]][1])
  if(chain.length > 10000){
    network.data[[ff]] = networkstats(10001, path.in, files[ff])
  }
  save(network.data,file=paste(path.in, 'networks_processed.RData', sep = ''))
}
