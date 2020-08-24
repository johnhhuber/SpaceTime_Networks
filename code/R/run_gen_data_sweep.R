# set seed 
set.seed(123)

# command line argument to 
args <- commandArgs(trailingOnly = T)

# get the index for the simulation to run 
index <- as.integer(args[1])

# load necessary functions 
source('functions_simulate_data.R')

# specify path to data and output
path.data = '../../data/sim_sweep/' 
path.out = '../../data/sim_sweep/data/'

# load parameter simulation sweep 
params.sweep <- read.csv(paste(path.data, 'parameters_sim_sweep.csv', sep = ''))

# specify parameters that are invariant across simulation sweep 
n_nodes <- 200
clustered <- TRUE
prob_symp <- 1
prob_trtd_symp <- 1
prob_trtd_asym <- 0

# simulate data set 
write.sim.data(file = index - 1,
               diffusion.coef = params.sweep$diffusion_coef[index],
               tau.s = params.sweep$tau_s[index],
               tau.l = params.sweep$tau_l[index],
               founder.prob = params.sweep$founder_prob[index],
               max_founder_date = params.sweep$max_founder_date[index],
               clustered = clustered,
               frac.cluster = params.sweep$frac_cluster[index],
               max_cases = n_nodes,
               max_tries = n_nodes,
               prob.symp = prob_symp,
               prob.trtd.asym = prob_trtd_asym,
               prob.trtd.symp = prob_trtd_symp,
               data.dir = path.out)
