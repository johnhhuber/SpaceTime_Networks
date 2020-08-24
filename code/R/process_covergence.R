# load necessary packages
if(!require(coda)){install.packages('coda'); library(coda)}

# specify path
path.in = commandArgs(trailingOnly = TRUE)[1]

# load all scalars
folders = list.files(path.in, pattern = 'replicate', full.names = T)
scalars = list()
mats = list()
num.replicates = length(folders)
Rc = list()

for(ii in 1:num.replicates)
{
  scalars[[ii]] = read.csv(paste(folders[ii], '/scalars.csv', sep  = ''))
  load(paste(folders[ii], '/networks_processed.RData', sep = ''))
  mats[[ii]] = network.data[[length(network.data)]]$edgeProbs
  Rc[[ii]] = rowSums(mats[[ii]][-1,-1])
}

# get the chains for convergence 
min.length.chain = min(sapply(1:length(scalars), function(x){return(nrow(scalars[[x]]))}))

combined.chains = mcmc.list(mcmc(scalars[[1]]$diffusion_coef[1:min.length.chain]), mcmc(scalars[[2]]$diffusion_coef[1:min.length.chain]), mcmc(scalars[[3]]$diffusion_coef[1:min.length.chain]), mcmc(scalars[[4]]$diffusion_coef[1:min.length.chain]), mcmc(scalars[[5]]$diffusion_coef[1:min.length.chain]))
gelman.diag(combined.chains)
gelman.plot(combined.chains)

combined.chains = mcmc.list(mcmc(scalars[[1]]$tau_s[1:min.length.chain]), mcmc(scalars[[2]]$tau_s[1:min.length.chain]), mcmc(scalars[[3]]$tau_s[1:min.length.chain]), mcmc(scalars[[4]]$tau_s[1:min.length.chain]), mcmc(scalars[[5]]$tau_s[1:min.length.chain]))
gelman.diag(combined.chains)
gelman.plot(combined.chains)

combined.chains = mcmc.list(mcmc(scalars[[1]]$tau_l[1:min.length.chain]), mcmc(scalars[[2]]$tau_l[1:min.length.chain]), mcmc(scalars[[3]]$tau_l[1:min.length.chain]), mcmc(scalars[[4]]$tau_l[1:min.length.chain]), mcmc(scalars[[5]]$tau_l[1:min.length.chain]))
gelman.diag(combined.chains)
gelman.plot(combined.chains)

# combine the posteriors of the chains after burn-in
burn.in = 1:5000
scalars.pooled = rbind(scalars[[1]][-burn.in, ], scalars[[2]][-burn.in, ],
                       scalars[[3]][-burn.in, ], scalars[[4]][-burn.in, ],
                       scalars[[5]][-burn.in, ])

write.csv(scalars.pooled, file = paste(path.in, '/scalars_pooled.csv', sep = ''))

# combine the posterior distributions of networks for the chains after burn-in
networks = list()
for(ii in 1:num.replicates)
{
  network.file <- list.files(path = folders[ii], pattern = 'networks.*.bz2', full.names = T)
  networks[[ii]] = read.csv(network.file, header = F)
}

networks.pooled = rbind(as.matrix((networks[[1]][-burn.in,1]), nrow = 5000, ncol = 1), 
                        as.matrix((networks[[2]][-burn.in,1]), nrow = 5000, ncol = 1),
                        as.matrix((networks[[3]][-burn.in,1]), nrow = 5000, ncol = 1),
                        as.matrix((networks[[4]][-burn.in,1]), nrow = 5000, ncol = 1),
                        as.matrix((networks[[5]][-burn.in,1]), nrow = 5000, ncol = 1))

write.table(networks.pooled, file = paste(path.in, '/networks_pooled.csv', sep = ''), sep = ',', row.names = F, col.names = F, quote = F)

# plot pairwise network plots
plot.replicates = 5
jpeg(file = paste(path.in, '/founder_correlation_plots.jpg', sep = ''), width = 13, height = 6.8, units = 'in', res = 200)
par(mfrow=c(2,5))
for(ii in 1:(plot.replicates-1))
{
  for(jj in (ii+1):plot.replicates)
  {
    plot(mats[[ii]][1,-1], mats[[jj]][1,-1], pch = 16, xlim = c(0,1), ylim = c(0,1), cex = 0.5, 
         bty = 'n', xaxs = 'i', yaxs = 'i', las = 1, xlab = paste('Network ', ii, ' Posterior Probabilities', sep = ''),
         ylab = paste('Network ', jj, ' Posterior Probabilities', sep = ''))
    mtext(paste('R = ', round(cor(c(mats[[ii]][1,-1]), c(mats[[jj]][1,-1])), digits = 2), sep = ''), 
          side = 1, at = 0.9, line = -1, col = 'red')
  }
}
dev.off()

jpeg(file = paste(path.in, '/edge_correlation_plots.jpg', sep = ''), width = 13, height = 6.8, units = 'in', res = 200)
par(mfrow=c(2,5))
for(ii in 1:(plot.replicates-1))
{
  for(jj in (ii+1):plot.replicates)
  {
    plot(mats[[ii]][-1,-1], mats[[jj]][-1,-1], pch = 16, xlim = c(0,1), ylim = c(0,1), cex = 0.5, 
         bty = 'n', xaxs = 'i', yaxs = 'i', las = 1, xlab = paste('Network ', ii, ' Posterior Probabilities', sep = ''),
         ylab = paste('Network ', jj, ' Posterior Probabilities', sep = ''))
    mtext(paste('R = ', round(cor(c(mats[[ii]][-1,-1]), c(mats[[jj]][-1,-1])), digits = 2), sep = ''), 
          side = 1, at = 0.9, line = -1, col = 'red')
  }
}
dev.off()

plot.replicates = 5
jpeg(file = paste(path.in, '/Rc_correlation_plots.jpg', sep = ''), width = 13, height = 6.8, units = 'in', res = 200)
par(mfrow=c(2,5))
max.Rc = max(unlist(Rc))
print(max.Rc)
for(ii in 1:(plot.replicates-1))
{
  for(jj in (ii+1):plot.replicates)
  {
    plot(Rc[[ii]], Rc[[jj]], pch = 16, xlim = c(0,max.Rc+1), ylim = c(0,max.Rc+1), cex = 0.5, 
         bty = 'n', xaxs = 'i', yaxs = 'i', las = 1, xlab = paste('Network ', ii, ' Rc', sep = ''),
         ylab = paste('Network ', jj, ' Rc', sep = ''))
    mtext(paste('R = ', round(cor(Rc[[ii]], Rc[[jj]]), digits = 2), sep = ''), 
          side = 1, at = max.Rc, line = -1, col = 'red')
  }
}
dev.off()
