# set seed 
set.seed(123)

# install necessary packages
if(!require(pomp)){install.packages('pomp'); library(pomp)}

# specify number of simulations for simulation sweep
n.sim <- 2000

# specify parameter ranges for simulation sweep
max_founder_date <- c(1, 25 * 365)
diffusion_coef <- c(0, 30)
tau_s <- c(0,1)
tau_l <- c(0,1)
frac_cluster <- c(0,1)
founder_prob <- c(0,1)

# construct sobol design
sweep <- sobolDesign(lower = c(max_founder_date = max_founder_date[1],
                               diffusion_coef = diffusion_coef[1],
                               tau_s = tau_s[1],
                               tau_l = tau_l[1],
                               frac_cluster = frac_cluster[1],
                               founder_prob = founder_prob[1]),
                     upper = c(max_founder_date = max_founder_date[2],
                               diffusion_coef = diffusion_coef[2],
                               tau_s = tau_s[2],
                               tau_l = tau_l[2],
                               frac_cluster = frac_cluster[2],
                               founder_prob = founder_prob[2]),
                     nseq = n.sim)

# write to file 
write.csv(x = sweep, file = '../../data/sim_sweep/parameters_sim_sweep.csv', row.names = F)