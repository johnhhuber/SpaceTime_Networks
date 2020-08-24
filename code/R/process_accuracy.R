# load necessary packages
require(igraph)

# load necessary functions
source('functions_network_statistics.R')

# generate path in
path.mc3  = commandArgs(trailingOnly = TRUE)[1]
file.network.true = commandArgs(trailingOnly = TRUE)[2]

# identify network mc3 files
files.network.mc3 = list.files(path.mc3)[grepl('networks', list.files(path.mc3))]

# generate accuracy
accuracy = networkAccuracy(file.path(path.mc3, files.network.mc3), file.network.true, burnIn = 10000, numLines = 10, verbose = F)

# save output
save(accuracy, file = file.path(path.mc3, 'validation_metrics.RData'))
