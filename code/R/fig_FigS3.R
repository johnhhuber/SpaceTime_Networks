# specify path 
path.in = '../../data/eswatini_inference/space_time_travel'
path.out = '../../output/figs/'

# load all scalars
folders = list.files(path.in, pattern = 'replicate', full.names = T)
mats = list()
num.replicates = length(folders)

for(ii in 1:num.replicates)
{
  load(paste(folders[ii], '/networks_processed.RData', sep = ''))
  mats[[ii]] = network.data[[length(network.data)]]$edgeProbs
}

# plot 
plot.replicates = 5
tiff(file = paste(path.out, '/S3_Fig.tiff', sep = ''), width = 7.5, height = 7.5, units = 'in', res = 300, compression = 'lzw')
par(mar = c(4.3,4.1,1.1,1.3))
par(mfrow=c(4,3))
for(ii in 1:(plot.replicates-1))
{
  for(jj in (ii+1):plot.replicates)
  {
    plot(mats[[ii]][-1,-1], mats[[jj]][-1,-1], pch = 16, xlim = c(0,1), ylim = c(0,1), cex = 0.5, 
         bty = 'n', xaxs = 'i', yaxs = 'i', las = 1, xlab = paste('Chain ', ii, ' Posterior Probs', sep = ''),
         ylab = paste('Chain ', jj, ' Posterior Probs', sep = ''))
    mtext(paste('R = ', round(cor(c(mats[[ii]][-1,-1]), c(mats[[jj]][-1,-1])), digits = 2), sep = ''), 
          side = 1, at = 0.9, line = -1, col = 'red')
  }
}
dev.off()