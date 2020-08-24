# specify path 
path.in = '../../data/eswatini_inference/space_time_travel'
path.out = '../../output/figs/'

# load all scalars
folders = list.files(path.in, pattern = 'replicate', full.names = T)
Rc = list()
num.replicates = length(folders)

for(ii in 1:num.replicates)
{
  load(paste(folders[ii], '/networks_processed.RData', sep = ''))
  Rc[[ii]] = rowSums(mats[[ii]][-1,-1])
}

# plot 
plot.replicates = 5
tiff(file = paste(path.out, '/S4_Fig.tiff', sep = ''), width = 7.5, height = 7.5, units = 'in', res = 300, compression = 'lzw')
par(mar = c(4.3,4.1,1.1,1.3))
par(mfrow=c(4,3))
max.Rc = max(unlist(Rc))
print(max.Rc)
for(ii in 1:(plot.replicates-1))
{
  for(jj in (ii+1):plot.replicates)
  {
    plot(Rc[[ii]], Rc[[jj]], pch = 16, xlim = c(0,1.1 * max.Rc), ylim = c(0,1.1 * max.Rc), cex = 0.5, 
         bty = 'n', xaxs = 'i', yaxs = 'i', las = 1, xlab = paste('Chain ', ii, ' Rc', sep = ''),
         ylab = paste('Chain ', jj, ' Rc', sep = ''))
    mtext(paste('R = ', round(cor(Rc[[ii]], Rc[[jj]]), digits = 2), sep = ''), 
          side = 1, at = max.Rc, line = -1, col = 'red')
  }
}
dev.off()