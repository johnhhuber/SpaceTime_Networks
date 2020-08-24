# install necessary packages
if(!require(MASS)){install.packages('MASS'); library(MASS)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}

# function to sample from null distribution of spatial and temporal distances
sample.null = function(num.samples)
{
  indices.sampled = arrayInd(sample(length(mat.SI), num.samples, replace = TRUE), dim(mat.SI))
  SI.sampled = mat.SI[indices.sampled]
  dist.sampled = sqrt(mat.dist[cbind(indices.sampled, rep(1, num.samples))]^2 + 
                        mat.dist[cbind(indices.sampled, rep(2, num.samples))]^2)
  return(cbind(SI.sampled, dist.sampled))
}

# function to sample from null distribution of spatial and temporal distances when you believe the travel history
sample.null.bth = function(num.samples)
{
  nodes = read.csv(file.path(path.in, 'nodes.csv'))
  indices.founders = which(nodes$History == 1)
  
  mat.SI.bth = mat.SI[,-indices.founders]
  mat.dist.bth = mat.dist[,-indices.founders,]
  
  indices.sampled = arrayInd(sample(length(mat.SI.bth), num.samples, replace = TRUE), dim(mat.SI.bth))
  SI.sampled = mat.SI.bth[indices.sampled]
  dist.sampled = sqrt(mat.dist.bth[cbind(indices.sampled, rep(1, num.samples))]^2 + 
                        mat.dist.bth[cbind(indices.sampled, rep(2, num.samples))]^2)
  return(cbind(SI.sampled, dist.sampled))
}


# specify paths in and out 
path.in = '../../data/eswatini_inference'
path.out = '../../output/figs/'

# load spatial and temporal distance matrices
load(list.files(path = path.in, pattern = 'distances.RData', full.names = T))

# list all directories
files <- list.files(path = path.in, pattern = '*distances.RData', 
                    recursive = T, full.names = T)
files <- files[grep(files, pattern = 'time')]

# aggregate posterior distances into single vector for each scenario
SI = list()
spat.dist = list()
for(ii in 1:length(files))
{
  load(files[ii])
  SI[[ii]] = unlist(posterior.SI)
  spat.dist[[ii]] = unlist(posterior.dist)
}

# generate null distribution
null.distribution = sample.null(1000000)
density.spatial = density(null.distribution[,2], na.rm = T)
indices.density = which(density.spatial$x < 150)

density.temporal = density(null.distribution[,1], na.rm = T)

null.distribution.bth = sample.null.bth(1000000)
density.spatial.bth = density(null.distribution.bth[,2], na.rm = T)
indices.density.bth = which(density.spatial.bth$x < 150)

# generate SI distribution from temporal likelihood (per Huber et al., 2016)
shift.STS = 44
k.STS = 737.7606011
p.STS = 0.8902353
days = seq(from = -20, to = 100, by = 1)

SI.dist = dnbinom(x = days + shift.STS, size = k.STS, prob = p.STS)

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot
tiff(filename = paste(path.out, 'Fig3.tiff', sep = ''), width = 6, height = 1.25 * 6, units = 'in', res = 300,
     compression = 'lzw')

layout(mat = matrix(1:15, nrow = 5, ncol = 3, byrow = F))
par(mar=c(3.1,4.5,1.1,0.1), oma = c(2.1,0,2.1,0))

density.SI = density(SI[[2]])
plot(density.SI, col = col2alpha(palette[2], alpha = 0.8), type = 'n',
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(-100,200), ylim = c(0,0.06), yaxs = 'i', xaxs = 'i', axes = FALSE, lwd = 1)
polygon(c(days, rev(days)), c(SI.dist, rep(0, length(days))), col = adjustcolor('#222222', 0.4), border = NA)
polygon(c(density.SI$x, rev(density.SI$x)), c(density.SI$y, rep(0, length(density.SI$x))), col = col2alpha(palette[2], alpha = 0.6),
        border = col2alpha(palette[2], alpha = 0.8))
lines(density.SI, col = col2alpha(palette[2], alpha = 0.8))
lines(density.temporal, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'Temporal (days)', side = 3, at = 75, line = 0.75, cex = 1, col = '#222222')
mtext(text = 'A', side = 3, line = 0.75, cex = 1, col = '#222222', at = -100)

par(mar=c(3.1,4.5,2.1,0.1))
density.SI = density(SI[[1]])
plot(density.SI, col = col2alpha(palette[1], alpha = 0.8), type = 'n',
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(-100,200), ylim = c(0,0.06), yaxs = 'i', xaxs = 'i', axes = FALSE, lwd = 1)
polygon(c(days, rev(days)), c(SI.dist, rep(0, length(days))), col = adjustcolor('#222222', 0.4), border = NA)
polygon(c(density.SI$x, rev(density.SI$x)), c(density.SI$y, rep(0, length(density.SI$x))), col = col2alpha(palette[1], alpha = 0.6),
        border = col2alpha(palette[1], alpha = 0.8))
lines(density.SI, col = col2alpha(palette[1], alpha = 0.8))
lines(density.temporal, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'C', side = 3, line = 0.75, cex = 1, col = '#222222', at = -100)

density.SI = density(SI[[3]])
plot(density.SI, col = col2alpha(palette[3], alpha = 0.8), type = 'n',
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(-100,200), ylim = c(0,0.06), yaxs = 'i', xaxs = 'i', axes = FALSE, lwd = 1)
polygon(c(days, rev(days)), c(SI.dist, rep(0, length(days))), col = adjustcolor('#222222', 0.4), border = NA)
polygon(c(density.SI$x, rev(density.SI$x)), c(density.SI$y, rep(0, length(density.SI$x))), col = col2alpha(palette[3], alpha = 0.6),
        border = col2alpha(palette[3], alpha = 0.8))
lines(density.SI, col = col2alpha(palette[3], alpha = 0.8))
lines(density.temporal, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'E', side = 3, line = 0.75, cex = 1, col = '#222222', at = -100)

density.SI = density(SI[[5]])
plot(density.SI, col = col2alpha(palette[5], alpha = 0.8), type = 'n',
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(-100,200), ylim = c(0,0.06), yaxs = 'i', xaxs = 'i', axes = FALSE, lwd = 1)
polygon(c(days, rev(days)), c(SI.dist, rep(0, length(days))), col = adjustcolor('#222222', 0.4), border = NA)
polygon(c(density.SI$x, rev(density.SI$x)), c(density.SI$y, rep(0, length(density.SI$x))), col = col2alpha(palette[5], alpha = 0.6),
        border = col2alpha(palette[5], alpha = 0.8))
lines(density.SI, col = col2alpha(palette[5], alpha = 0.8))
lines(density.temporal, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'G', side = 3, line = 0.75, cex = 1, col = '#222222', at = -100)

density.SI = density(SI[[4]])
plot(density.SI, col = col2alpha(palette[4], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(-100,200), ylim = c(0,0.06), yaxs = 'i', xaxs = 'i', axes = FALSE, lwd = 1)
polygon(c(days, rev(days)), c(SI.dist, rep(0, length(days))), col = adjustcolor('#222222', 0.4), border = NA)
polygon(c(density.SI$x, rev(density.SI$x)), c(density.SI$y, rep(0, length(density.SI$x))), col = col2alpha(palette[4], alpha = 0.6),
        border = col2alpha(palette[4], alpha = 0.8))
lines(density.SI, col = col2alpha(palette[4], alpha = 0.8))
lines(density.temporal, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'I', side = 3, line = 0.75, cex = 1, col = '#222222', at = -100)

par(mar=c(3.1,4.5,1.1,0.1), xpd = FALSE)
density.spat = density(spat.dist[[2]], na.rm = T)
plot(density.spat, col = col2alpha(palette[2], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(0,175), ylim = c(0,0.08), axes = FALSE, yaxs = 'i', xaxs = 'i', lwd = 1)
polygon(c(density.spat$x, rev(density.spat$x)), c(density.spat$y, rep(0, length(density.spat$x))), col = col2alpha(palette[2], alpha = 0.6),
        border = col2alpha(palette[2], alpha = 0.8))
lines(density.spatial$x, density.spatial$y, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'Spatial (km)', side = 3, line = 0.75, cex = 1, col = '#222222')
mtext(text = expression(underline('Scales of Transmission')), side = 3, line = 2, at = -5.5, cex = 1, font = 2, col = '#222222')
mtext(text = 'B', side = 3, line = 0.75, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,4.5,2.1,0.1))
density.spat = density(spat.dist[[1]], na.rm = T)
plot(density.spat, col = col2alpha(palette[1], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(0,175), ylim = c(0,0.08), axes = FALSE, yaxs = 'i', xaxs = 'i', lwd = 1)
polygon(c(density.spat$x, rev(density.spat$x)), c(density.spat$y, rep(0, length(density.spat$x))), col = col2alpha(palette[1], alpha = 0.6),
        border = col2alpha(palette[1], alpha = 0.8))
lines(density.spatial.bth$x, density.spatial.bth$y, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'D', side = 3, line = 0.75, cex = 1, col = '#222222', at = 0)

density.spat = density(spat.dist[[3]], na.rm = T)
plot(density.spat, col = col2alpha(palette[3], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(0,175), ylim = c(0,0.08), axes = FALSE, yaxs = 'i', xaxs = 'i', lwd = 1)
polygon(c(density.spat$x, rev(density.spat$x)), c(density.spat$y, rep(0, length(density.spat$x))), col = col2alpha(palette[3], alpha = 0.6),
        border = col2alpha(palette[3], alpha = 0.8))
lines(density.spatial$x, density.spatial$y, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'F', side = 3, line = 0.75, cex = 1, col = '#222222', at = 0)

density.spat = density(spat.dist[[5]], na.rm = T)
plot(density.spat, col = col2alpha(palette[5], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(0,175), ylim = c(0,0.08), axes = FALSE, yaxs = 'i', xaxs = 'i', lwd = 1)
polygon(c(density.spat$x, rev(density.spat$x)), c(density.spat$y, rep(0, length(density.spat$x))), col = col2alpha(palette[5], alpha = 0.6),
        border = col2alpha(palette[5], alpha = 0.8))
lines(density.spatial$x, density.spatial$y, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'H', side = 3, line = 0.75, cex = 1, col = '#222222', at = 0)

density.spat = density(spat.dist[[4]], na.rm = T)
plot(density.spat, col = col2alpha(palette[4], alpha = 0.8),
     bty = 'n', xlab = '', ylab = '', main = '', las = 1, xlim = c(0,175), ylim = c(0,0.08), axes = FALSE, yaxs = 'i', xaxs = 'i', lwd = 1)
polygon(c(density.spat$x, rev(density.spat$x)), c(density.spat$y, rep(0, length(density.spat$x))), col = col2alpha(palette[4], alpha = 0.6),
        border = col2alpha(palette[4], alpha = 0.8))
lines(density.spatial.bth$x, density.spatial.bth$y, col = adjustcolor('#222222', 0.8), lwd = 1.5, lty = 2)
axis(side = 1, lwd = 1.5, las = 1, col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col.axis = '#222222')
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'J', side = 3, line = 0.75, cex = 1, col = '#222222', at = 0)

par(mar=c(2.1,1.1,1.5,0.1))

# data sources included
plot.new()
mtext(text = expression(underline('Inference Settings')), side = 3, line = 2, cex = 1, font = 1, adj = 0)
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -2.3, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, line = -2.3, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -2.3, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, line = -2.3, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Ignore')), col = '#222222', side = 1, line = -2.3, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Ignore'), col = palette[3], side = 1, line = -2.3, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[5], side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[5], side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -2.3, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[5], side = 1, line = -2.3, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[4], side = 1, line = -6.3, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[4], side = 1, line = -4.3, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -2.3, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[4], side = 1, line = -2.3, cex = 1, adj = 0)


par(fig = c(0,1,0,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend('bottom', legend = c('Null Distribution', 'Serial Interval Distribution'),
       lwd = c(1.5, NA), col = c(adjustcolor('#222222', 0.8), adjustcolor('#222222', 0.4)), lty = c(2,NA), pch = c(NA, 15), pt.cex = c(NA, 1.5), bty = 'n', xpd = TRUE, horiz = TRUE, cex = 1, inset = c(-20,0))

dev.off()