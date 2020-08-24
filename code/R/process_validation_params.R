# specify path in 
path.in = commandArgs(trailingOnly = T)[1]
path.setting = commandArgs(trailingOnly = T)[2]

# specifyp path out 
path.out = path.in

# read in nodes 
nodes = read.csv('../../data/eswatini_inference/nodes.csv')

# read in scalar estimates
scalars = read.csv(list.files(path = path.in, pattern = 'scalars_pooled.csv', full.names = T))

# read in settings files
inference.setting <- strsplit(path.in, '/')[[1]][5]
setting.file = list.files(path.setting, pattern = paste(inference_setting, '.csv', sep = ''), full.names = T)
setting = read.csv(paste(path.setting, settings.file, sep = ''), header = FALSE)

# compute posterior distribution of founder probabilities
length.chain = nrow(scalars)

chain = readLines(paste(path.in, '/networks_pooled.csv.bz2', sep = ''))
chain = strsplit(chain, ';')
numNodes = as.numeric(strsplit(chain[[1]][length(chain[[1]])],'-')[[1]][3])

founder.prob = rep(NA, length.chain)

for(ii in 1:length.chain)
{
  if(ii %% 1000 == 0){print(ii)}
  networks = matrix(0, numNodes + 1, numNodes + 1)
  for(jj in 1 : length(chain[[ii]])){
    ancestor = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][1]
    k = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][2]
    descendant = strsplit(as.character(chain[[ii]][jj]), '-')[[1]][3]
    
    if(ancestor == 's'){ii.in = 1}
    if(ancestor != 's'){ii.in = as.integer(ancestor) + 1}
    jj.in = as.integer(descendant) + 1
    
    networks[ii.in, jj.in] = as.integer(k)
  }
  founder.prob[ii] = sum(networks[1,-1]) / numNodes
}

# establish median parameters 
median.founder.prob = median(founder.prob)
median.diffusion.coef = ifelse(setting[which(setting[,1] == "PROP_DIFFUSION_COEF"), 2] > 0, 
                               median(scalars$diffusion_coef), setting[which(setting[,1] == "DIFFUSION_min"), 2])
median.tau.s = ifelse(setting[which(setting[,1] == "PROP_TAU_S"), 2] > 0,
                      median(scalar$tau_s), 1.0)
median.tau.l = ifelse(setting[which(setting[,1] == "PROP_TAU_L"), 2] > 0,
                      median(scalar$tau_l), 0.0)

# establish other necessary parameters
num.files = 1
max.founder.date = max(nodes$time, na.rm = T) - min(nodes$time, na.rm = T)
area.swazi = 17364.0
founder.radius = sqrt(area.swazi / pi)

# save output to RData structure
save(median.founder.prob, 
     median.diffusion.coef, 
     median.tau.s, 
     median.tau.l,
     num.files, 
     max.founder.date, 
     founder.radius,
     founder.prob,
     file = paste(path.out, '/validation_parameters.RData', sep = '')) 


