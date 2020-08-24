# load necessary packages
if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(pomp)){install.packages('pomp'); library(pomp)}
if(!require(spatstat)){install.packages('spatstat');library(spatstat)}

# load necessary functions
source('functions_simulate_data.R')
load('../../data/network_sim/Huber_MalariaJournal/Fig2.RData')

# specify paths in and out 
path.in = commandArgs(trailingOnly = T)[1]
path.out = path.in

# load parameters
load(file.path(path = path.in, pattern = 'validation_parameters.RData', full.names = t))

# load eswatini nodes file
nodes = read.csv('../../data/eswatini_inference/nodes.csv')
numNodes = nrow(nodes)

# specify number of each case type
num.ST = length(which(nodes$type == 'ST'))
num.SU = length(which(nodes$type == 'SU'))
num.AT = length(which(nodes$type == 'AT'))
num.AU = length(which(nodes$type == 'AU'))

prob.symp = (num.ST + num.SU) / numNodes
prob.trtd.symp = (num.ST) / (num.ST + num.SU)
prob.trtd.asym = (num.AT) / (num.AT + num.AU)

# write files to necessary output folder
max.cases = numNodes
num.files = 1
for(ii in 1:num.files)
{
  print(ii)
  write.sim.data(
    file = ii - 1,
    num.mcmcs = 1,
    diffusion.coef = median.diffusion.coef,
    tau.s = median.tau.s,
    tau.l = median.tau.l,
    founder.prob = median.founder.prob,
    max_founder_date = max.founder.date,
    radius = founder.radius,
    imported.prob = 1.0,
    max_cases = max.cases - 1,
    max_tries = max.cases,
    prob.symp = prob.symp,
    prob.trtd.symp = prob.trtd.symp,
    prob.trtd.asym = prob.trtd.asym,
    data.dir = path.out)
}
