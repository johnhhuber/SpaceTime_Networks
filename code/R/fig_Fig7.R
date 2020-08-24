# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load necessary data
load('../../data/figs/fig_7/fig7.RData')

# specify color palette 
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = '../../output/figs/Fig7.tiff', width = 7.5, height = 3, units = 'in', res = 300, compression = 'lzw')

layout(mat = matrix(c(1,2,3,1,2,3,4,4,4), nrow = 3, ncol = 3, byrow = T))
par(mar = c(3.1,3.6,1.1,1.1))

plot(1- df.stt$founder_prob, df.stt$Rc, type = 'n', bty = 'n', las = 1, xaxs = 'i', yaxs = 'i',
     xlim = c(0,1.025), ylim = c(0,1.025), xlab = '', ylab = '')
lines(0:1, 0:1)
points(1 - df.stt$founder_prob[-indices.stt], df.stt$Rc[-indices.stt], pch = 16, cex = 0.5,
       col = col2alpha(palette[2], alpha = 0.5))

points(1 - df.stt$founder_prob[indices.stt], df.stt$Rc[indices.stt], pch = 21, cex = 0.6,
       bg = palette[2], col = '#222222', lwd = 0.5)
mtext(expression(paste('True ', 'R'['c'], sep = '')), side = 1, line = 2.3, cex = 0.8)
mtext(expression(paste('Estimated ', 'R'['c'], sep = '')), side = 2, line = 2.3, cex = 0.8)
mtext(side = 3, at = 0, 'A', cex = 0.8)

plot(1- df.stb$founder_prob, df.stb$Rc, type = 'n', bty = 'n', las = 1, xaxs = 'i', yaxs = 'i',
     xlim = c(0,1.025), ylim = c(0,1.025), xlab = '', ylab = '')
lines(0:1, 0:1)
points(1 - df.stb$founder_prob[-indices.stb], df.stb$Rc[-indices.stb], pch = 16, cex = 0.6,
       col = col2alpha(palette[1], alpha = 0.5))

points(1 - df.stb$founder_prob[indices.stb], df.stb$Rc[indices.stb], pch = 21, cex = 0.6,
       bg = palette[1], col = '#222222', lwd = 0.5)
mtext(expression(paste('True ', 'R'['c'], sep = '')), side = 1, line = 2.3, cex = 0.8)
mtext(expression(paste('Estimated ', 'R'['c'], sep = '')), side = 2, line = 2.3, cex = 0.8)
mtext(side = 3, at = 0, 'B', cex = 0.8)

plot(1- df.stn$founder_prob, df.stn$Rc, type = 'n', bty = 'n', las = 1, xaxs = 'i', yaxs = 'i',
     xlim = c(0,1.025), ylim = c(0,1.025), xlab = '', ylab = '')
lines(0:1, 0:1)
points(1 - df.stn$founder_prob[-indices.stn], df.stn$Rc[-indices.stn], pch = 16, cex = 0.5,
       col = col2alpha(palette[3], alpha = 0.5))
points(1 - df.stn$founder_prob[indices.stn], df.stn$Rc[indices.stn], pch = 21, cex = 0.6,
       bg = palette[3], col = '#222222', lwd = 0.5)
mtext(expression(paste('True ', 'R'['c'], sep = '')), side = 1, line = 2.3, cex = 0.8)
mtext(expression(paste('Estimated ', 'R'['c'], sep = '')), side = 2, line = 2.3, cex = 0.8)
mtext(side = 3, at = 0, 'C', cex = 0.8)

plot(0:6, 0:6, xlim = c(0.5, 5), ylim = c(0,1.02), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
mtext(text = expression(underline('Inference Settings')), side = 1, at = 0.5, line = -3.75, cex = 0.8, font = 1, adj = 0, col = '#222222')
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 0.5, line = -2, cex = 0.8, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, at = 0.5, line = -2, cex = 0.8, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 0.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, at = 0.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, at = 0.5, line = 1, cex = 0.8, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, at = 0.5, line = 1, cex = 0.8, adj = 0)

mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 1.5, line = -2, cex = 0.8, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, at = 1.5, line = -2, cex = 0.8, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 1.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, at = 1.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, at = 1.5, line = 1, cex = 0.8, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, at = 1.5, line = 1, cex = 0.8, adj = 0)

mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 2.5, line = -2, cex = 0.8, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, at = 2.5, line = -2, cex = 0.8, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 2.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, at = 2.5, line = -0.5, cex = 0.8, adj = 0)
mtext(expression('Travel:' * phantom(' Exclude')), col = '#222222', side = 1, at = 2.5, line = 1, cex = 0.8, adj = 0)
mtext(expression(phantom('Travel:') * ' Exclude'), col = palette[3], side = 1, at = 2.5, line = 1, cex = 0.8, adj = 0)

dev.off()
