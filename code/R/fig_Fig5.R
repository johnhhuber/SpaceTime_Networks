# load necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# specify paths in and out 
path.in = '../../data/eswatini_inference'
path.out = '../../output/figs/'

# load the posterior scalars files 
files<- list.files(path = path.in, pattern = 'scalars_pooled.csv', recursive = T, full.names = T)
files <- files[grepl('validation', files)]

scalars = list()
for(ii in 1:length(files))
{
  scalars[[ii]] = read.csv(file = files[ii], header = TRUE)
}

# load the true scalars
files<- list.files(path = path.in, pattern = 'validation_parameters.RData', recursive = T, full.names = T)
diffusion.true = rep(NA, length(files))
tau.s.true = rep(NA, length(files))
tau.l.true = rep(NA, length(files))
for(ii in 1:length(files))
{
  load(files[ii])
  diffusion.true[ii] = median.diffusion.coef
  tau.s.true[ii] = median.tau.s
  tau.l.true[ii] = median.tau.l
}

# create color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot
tiff(filename = paste(path.out, 'Fig5.tiff', sep = ''), width = 1.25 * 6, height = 6, units = 'in', res = 300, compression = 'lzw')

layout(mat = matrix(1:20, nrow = 5, ncol = 4, byrow = F))
par(mar=c(3.1,4.5,1.1,0.1), oma = c(0,0,3.1,0))

hist(scalars[[2]]$diffusion_coef, col = col2alpha(palette[2], alpha = 0.8), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,10), ylim = c(0,1), yaxs ='i', xaxs = 'i', axes = FALSE, freq = FALSE)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = diffusion.true[2], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'D', side = 3, line = 1, cex = 1, col = '#222222', font = 2)
mtext(text = 'A', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,4.5,2.1,0.1))
hist(scalars[[1]]$diffusion_coef, col = col2alpha(palette[1], alpha = 0.8), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,10), ylim = c(0,1), yaxs ='i', xaxs = 'i', axes = FALSE, freq = FALSE)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = diffusion.true[1], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'D', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

hist(scalars[[3]]$diffusion_coef, col = col2alpha(palette[3], alpha = 0.8), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,10), ylim = c(0,1), yaxs ='i', xaxs = 'i', axes = FALSE, freq = FALSE)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = diffusion.true[3], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'G', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
legend('center', legend = c('NA'), text.col = palette[5], bty = 'n', cex = 1.5)
mtext(text = 'J', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
title(xlab = '', ylab = 'Density', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
legend('center', legend = c('NA'), text.col = palette[4], bty = 'n', cex = 1.5)
mtext(text = 'M', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,3.5,1.1,0.1))
hist(scalars[[2]]$tau_s, col = col2alpha(palette[2], alpha = 0.0), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,5), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE)
polygon(x = c(seq(from = 0, to = 1, by = 0.01), rev(seq(from = 0, to = 1, by = 0.01))), c(dbeta(seq(from = 0, to = 1, by = 0.01), 12, 3), rep(0, length(seq(from = 0, to = 1, by = 0.01)))),
        col = adjustcolor('#222222', 0.4), border = NA)
hist(scalars[[2]]$tau_s, col = col2alpha(palette[2], alpha = 0.8), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,5), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE, add = T)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = tau.s.true[2], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = expression(tau['s']), side = 3, line = 1, cex = 1, col = '#222222', font = 2)
mtext(text = expression(underline('Parameters')), side = 3, line = 3, cex = 1, col = '#222222')
mtext(text = 'B', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,3.5,2.1,0.1))
plot.new()
legend('center', legend = c('NA'), text.col = palette[1], bty = 'n', cex = 1.5)
mtext(text = 'E', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
legend('center', legend = c('NA'), text.col = palette[3], bty = 'n', cex = 1.5)
mtext(text = 'H', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

hist(scalars[[5]]$tau_s, col = col2alpha(palette[5], alpha = 0.0), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,5), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE)
polygon(x = c(seq(from = 0, to = 1, by = 0.01), rev(seq(from = 0, to = 1, by = 0.01))), c(dbeta(seq(from = 0, to = 1, by = 0.01), 12, 3), rep(0, length(seq(from = 0, to = 1, by = 0.01)))),
        col = adjustcolor('#222222', 0.4), border = NA)
hist(scalars[[5]]$tau_s, col = col2alpha(palette[5], alpha = 0.8), border = 'white',
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,5), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE, add = T)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = tau.s.true[5], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'K', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
legend('center', legend = c('NA'), text.col = palette[4], bty = 'n', cex = 1.5)
mtext(text = 'N', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,3.5,1.1,1.1))
hist(scalars[[2]]$tau_l, col = col2alpha(palette[2], alpha = 0.0), border = palette[2],
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,30), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE)
polygon(x = c(seq(from = 0, to = 1, by = 0.01), rev(seq(from = 0, to = 1, by = 0.01))), c(dbeta(seq(from = 0, to = 1, by = 0.01), 3, 12), rep(0, length(seq(from = 0, to = 1, by = 0.01)))),
        col = adjustcolor('#222222', 0.4), border = NA)
hist(scalars[[2]]$tau_l, col = col2alpha(palette[2], alpha = 1.0), border = palette[2],
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0,1), ylim = c(0,30), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE, add = T)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline(v = tau.l.true[2], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = expression(tau['l']), side = 3, line = 1, cex  = 1, col = '#222222', font = 2)
mtext(text = 'C', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(3.1,3.5,2.1,1.1))
plot.new()
legend('center', legend = c('NA'), text.col = palette[1], bty = 'n', cex = 1.5)
mtext(text = 'F', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
legend('center', legend = c('NA'), text.col = palette[3], bty = 'n', cex = 1.5)
mtext(text = 'I', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

hist(scalars[[5]]$tau_l, col = col2alpha(palette[5], alpha = 0.0), border = palette[5],
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0, 1), ylim = c(0,30), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE)
polygon(x = c(seq(from = 0, to = 1, by = 0.01), rev(seq(from = 0, to = 1, by = 0.01))), c(dbeta(seq(from = 0, to = 1, by = 0.01), 3, 12), rep(0, length(seq(from = 0, to = 1, by = 0.01)))),
        col = adjustcolor('#222222', 0.4), border = NA)
hist(scalars[[5]]$tau_l, col = col2alpha(palette[5], alpha = 0.8), border = palette[5],
     xlab = '', ylab = '', main = '', las = 1, xlim = c(0, 1), ylim = c(0,30), axes = FALSE, yaxs = 'i', xaxs = 'i', freq = FALSE, add = T)
axis(side = 1, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
axis(side = 2, lwd = 1.5, las = 1, col = '#222222', col.axis = '#222222')
abline( v= tau.l.true[5], lwd = 1.5, lty = 1, col = adjustcolor('#222222', 0.8))
title(xlab = '', ylab = '', font.lab = 1, cex.lab = 1.5, col.lab = '#222222')
mtext(text = 'L', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

plot.new()
legend('center', legend = c('NA'), text.col = palette[4], bty = 'n', cex = 1.5)
mtext(text = 'O', side = 3, line = 1, cex = 1, col = '#222222', at = 0)

par(mar=c(2.1,1.1,1.5,0.1))

# data sources included
plot.new()
mtext(text = expression(underline('Inference Settings')), side = 3, line = 3, cex = 1, adj = 0)
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, line = -1.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, line = -1.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' No')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' No'), col = palette[3], side = 1, line = -1.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[5], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[5], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[5], side = 1, line = -1.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[4], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[4], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[4], side = 1, line = -1.25, cex = 1, adj = 0)

par(fig = c(0,1,0,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
legend('bottom', legend = c('True Value','Prior Distribution'),
       lwd = c(1.5, NA), col = c(adjustcolor('#222222', 0.8), adjustcolor('#222222', 0.4)), lty = c(1,NA), pch = c(NA,15), bty = 'n', xpd = TRUE, horiz = TRUE, cex = 1, pt.cex = 2, inset = c(-20,0))

dev.off()
