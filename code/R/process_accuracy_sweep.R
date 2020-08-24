# load necessary packages
require(igraph)

# load necessary functions
source('processNetworkStatistics.R')

# generate path in
path.mc3 = commandArgs(trailingOnly = TRUE)[1]
path.data = commandArgs(trailingOnly = TRUE)[2]
path.out = commandArgs(trailingOnly = TRUE)[3]
index = as.integer(commandArgs(trailingOnly = TRUE)[4])
scenario = commandArgs(trailingOnly = TRUE)[5]

# identify network mc3 files
file.network.mc3 = paste(path.mc3, 'networks_', scenario, '_', index, '.csv.bz2', sep = '')

# specify true network file 
file.network.true = paste(path.data, 'network_true_', index, '_FULL.csv', sep = '')

# generate accuracy
accuracy = networkAccuracy(file.network.mc3, file.network.true, burnIn = 1000, numLines = 10, verbose = F)

# generate statistics
statistics = networkStatistic(file.network.mc3, burnIn = 1000, numLines = 10, verbose = F)

# save output
save(accuracy, statistics, file = paste(path.out, 'validation_', scenario, '_', index, '.RData', sep = ''))

