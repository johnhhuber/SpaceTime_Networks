# install necesary packages
if(!require(boot)){install.packages('boot'); library(boot)}
if(!require(mgcv)){install.packages('mgcv'); library(mgcv)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# list the accuracy files 
files.acc <- list.files('../../output/sim_sweep/accuracy/', full.names = T)

# load the parameters used in the simulation sweep 
params <- read.csv('../../data/sim_sweep/parameters_sim_sweep.csv')

# get the number of observations 
n.obs <- nrow(params)

# create data frame to store results 
df.stb <- data.frame(sim_id = rep(NA, n.obs),
                     max_founder_date = rep(NA, n.obs),
                     diffusion_coef = rep(NA, n.obs),
                     tau_s = rep(NA, n.obs),
                     tau_l = rep(NA, n.obs),
                     frac_cluster = rep(NA, n.obs), 
                     founder_prob = rep(NA, n.obs), 
                     founder_acc = rep(NA, n.obs),
                     local_acc = rep(NA, n.obs),
                     classification_acc = rep(NA, n.obs),
                     parent_acc = rep(NA, n.obs),
                     outbreak_acc = rep(NA, n.obs), 
                     Rc = rep(NA, n.obs))

df.stn <- data.frame(sim_id = rep(NA, n.obs),
                     max_founder_date = rep(NA, n.obs),
                     diffusion_coef = rep(NA, n.obs),
                     tau_s = rep(NA, n.obs),
                     tau_l = rep(NA, n.obs),
                     frac_cluster = rep(NA, n.obs), 
                     founder_prob = rep(NA, n.obs), 
                     founder_acc = rep(NA, n.obs),
                     local_acc = rep(NA, n.obs),
                     classification_acc = rep(NA, n.obs),
                     parent_acc = rep(NA, n.obs),
                     outbreak_acc = rep(NA, n.obs), 
                     Rc = rep(NA, n.obs))

df.stt <- data.frame(sim_id = rep(NA, n.obs),
                     max_founder_date = rep(NA, n.obs),
                     diffusion_coef = rep(NA, n.obs),
                     tau_s = rep(NA, n.obs),
                     tau_l = rep(NA, n.obs),
                     frac_cluster = rep(NA, n.obs), 
                     founder_prob = rep(NA, n.obs), 
                     founder_acc = rep(NA, n.obs),
                     local_acc = rep(NA, n.obs),
                     classification_acc = rep(NA, n.obs),
                     parent_acc = rep(NA, n.obs),
                     outbreak_acc = rep(NA, n.obs), 
                     Rc = rep(NA, n.obs))

# loop through and fill the data frame
for(ii in 1:length(files.acc))
{
  # track progress
  print(ii)
  
  # load the file 
  load(files.acc[ii])
  
  # get the number associated with the simulation 
  sim_id <- as.integer(strsplit(strsplit(files.acc[ii], split = '_')[[1]][4], '.R')[[1]][1])
  
  # get the scenario associated with the simulation 
  scenario <- strsplit(files.acc[ii], split = '_')[[1]][3]
  
  # add in the parameters and accuracy metrics 
  if(scenario == 'stb')
  {
    df.stb$sim_id[sim_id] <- sim_id
    
    # add the parameters 
    df.stb$max_founder_date[sim_id] <- params$max_founder_date[sim_id+1]
    df.stb$diffusion_coef[sim_id] <- params$diffusion_coef[sim_id+1]
    df.stb$tau_s[sim_id] <- params$tau_s[sim_id+1]
    df.stb$tau_l[sim_id] <- params$tau_l[sim_id+1]
    df.stb$frac_cluster[sim_id] <- params$frac_cluster[sim_id+1]
    df.stb$founder_prob[sim_id] <- params$founder_prob[sim_id+1]
    
    # add in the accuracy metrics - we are considering the median of the posterior distribution for simplicity
    df.stb$founder_acc[sim_id] <- quantile(accuracy$classif.import.acc, probs = 0.50)
    df.stb$local_acc[sim_id] <- quantile(accuracy$classif.local.acc, probs = 0.50)
    df.stb$classification_acc[sim_id] <- quantile(accuracy$classif.acc, probs = 0.50)
    df.stb$parent_acc[sim_id] <- quantile(accuracy$edge.acc, probs = 0.50)
    df.stb$outbreak_acc[sim_id] <- quantile(accuracy$outbreak.acc, probs = 0.50)
    
    # compute the median percent error of Rc estimate 
    Rc <- as.numeric(unlist(statistics))
    df.stb$Rc[sim_id] <- quantile(Rc, probs = 0.50)
  }
  if(scenario == 'stn')
  {
    df.stn$sim_id[sim_id] <- sim_id
    
    # add the parameters 
    df.stn$max_founder_date[sim_id] <- params$max_founder_date[sim_id+1]
    df.stn$diffusion_coef[sim_id] <- params$diffusion_coef[sim_id+1]
    df.stn$tau_s[sim_id] <- params$tau_s[sim_id+1]
    df.stn$tau_l[sim_id] <- params$tau_l[sim_id+1]
    df.stn$frac_cluster[sim_id] <- params$frac_cluster[sim_id+1]
    df.stn$founder_prob[sim_id] <- params$founder_prob[sim_id+1]
    
    # add in the accuracy metrics - we are considering the median of the posterior distribution for simplicity
    df.stn$founder_acc[sim_id] <- quantile(accuracy$classif.import.acc, probs = 0.50)
    df.stn$local_acc[sim_id] <- quantile(accuracy$classif.local.acc, probs = 0.50)
    df.stn$classification_acc[sim_id] <- quantile(accuracy$classif.acc, probs = 0.50)
    df.stn$parent_acc[sim_id] <- quantile(accuracy$edge.acc, probs = 0.50)
    df.stn$outbreak_acc[sim_id] <- quantile(accuracy$outbreak.acc, probs = 0.50)
    
    # compute the median percent error of Rc estimate 
    Rc <- as.numeric(unlist(statistics))
    df.stn$Rc[sim_id] <- quantile(Rc, probs = 0.50)
  }
  if(scenario == 'stt')
  {
    df.stt$sim_id[sim_id] <- sim_id
    
    # add the parameters 
    df.stt$max_founder_date[sim_id] <- params$max_founder_date[sim_id+1]
    df.stt$diffusion_coef[sim_id] <- params$diffusion_coef[sim_id+1]
    df.stt$tau_s[sim_id] <- params$tau_s[sim_id+1]
    df.stt$tau_l[sim_id] <- params$tau_l[sim_id+1]
    df.stt$frac_cluster[sim_id] <- params$frac_cluster[sim_id+1]
    df.stt$founder_prob[sim_id] <- params$founder_prob[sim_id+1]
    
    # add in the accuracy metrics - we are considering the median of the posterior distribution for simplicity
    df.stt$founder_acc[sim_id] <- quantile(accuracy$classif.import.acc, probs = 0.50)
    df.stt$local_acc[sim_id] <- quantile(accuracy$classif.local.acc, probs = 0.50)
    df.stt$classification_acc[sim_id] <- quantile(accuracy$classif.acc, probs = 0.50)
    df.stt$parent_acc[sim_id] <- quantile(accuracy$edge.acc, probs = 0.50)
    df.stt$outbreak_acc[sim_id] <- quantile(accuracy$outbreak.acc, probs = 0.50)
    
    # compute the median percent error of Rc estimate 
    Rc <- as.numeric(unlist(statistics))
    df.stt$Rc[sim_id] <- quantile(Rc, probs = 0.50)
  }
}

# subset to remove any NAs 
df.stb <- subset(df.stb, !(is.na(founder_acc) | is.na(local_acc) | is.na(classification_acc) | 
                             is.na(parent_acc) | is.na(outbreak_acc)))
df.stn <- subset(df.stn, !(is.na(founder_acc) | is.na(local_acc) | is.na(classification_acc) | 
                             is.na(parent_acc) | is.na(outbreak_acc)))
df.stt <- subset(df.stt, !(is.na(founder_acc) | is.na(local_acc) | is.na(classification_acc) | 
                             is.na(parent_acc) | is.na(outbreak_acc)))

# specify the indices to accent in Rc plots 
indices.stb <- which(abs((df.stb$tau_l / (1 + df.stb$tau_l - df.stb$tau_s)) - df.stb$founder_prob) < 0.05)
indices.stt <- which(df.stt$max_founder_date / (200 * df.stt$founder_prob) > 100)
indices.stn <- which(df.stn$max_founder_date / (200 * df.stn$founder_prob) > 100)

# save output 
save(df.stb, 
     df.stn,
     df.stt,
     indices.stb,
     indices.stt,
     indices.stn,
     file = '../../data/figs/fig_7/fig7.RData')
