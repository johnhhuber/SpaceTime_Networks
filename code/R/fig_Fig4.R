# install necessary packages
if(!require(raster)){install.packages('raster'); library(raster)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# specify paths in and out 
path.in = '../../data/figs/fig_4'
path.out = '../../output/figs/'

# load necessary data
load(list.files(path = path.in, 
                pattern = '*.RData', 
                full.names = T))
# load shp file 
load('../../data/eswatini_inference/rasters.RData')

# create color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot the maps of proportion imported and Rc 
tiff(filename = paste(path.out, 'Fig4.tiff', sep = ''), width = 6.5, height = 8, units = 'in', res = 300,
     compression = 'lzw')
layout(mat = matrix(1:15, nrow = 5, ncol = 3, byrow = FALSE))
par(mar = c(0,0,1.1,0), oma = c(2.1,0,2.1,0))

colors <- colorRampPalette(c('white', palette[2]))(11)

plot(map.swazi.reproject, border = '#222222')
plot(imported.coord.pred.2, zlim = c(0,1), add = T,legend.mar = 6.1,
     col = colors, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = expression(underline('Importation Proportion')), side = 3, line = 1, cex = 1, font = 2, col = '#222222')
mtext(text = 'A', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

colors <- colorRampPalette(c('white', palette[2]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(Rc.coord.pred.2, add = T, zlim = c(0,1.1), legend.mar = 6.1,
     col = colors, legend = T, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = expression(underline('R'['c'])), side = 3, line = 1, cex = 1, font = 2, col = '#222222')
mtext(text = 'B', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

par(mar=c(2.1,1.1,1.1,0.1))
plot.new()
mtext(text = expression(underline('Inference Settings')), side = 3, line = 1, cex = 1, font = 1, adj = 0, col = '#222222')
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, line = -5.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, line = -1.25, cex = 1, adj = 0)

par(mar=c(0,0,0,0))
colors <- colorRampPalette(c('white', palette[1]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(imported.coord.pred.1, zlim = c(0,1), add = T, legend.mar = 6.1,
     col = colors, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'C', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

colors <- colorRampPalette(c('white', palette[1]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(Rc.coord.pred.1, add = T, zlim = c(0,1.1), legend.mar = 6.1,
     col = colors, legend = T, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'D', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

par(mar=c(2.1,1.1,1.1,0.1))
plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -2.125, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, line = -2.125, cex = 1, adj = 0)

par(mar=c(0,0,0,0))
colors <- colorRampPalette(c('white', palette[3]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(imported.coord.pred.3, zlim = c(0,1), add = T, legend.mar = 6.1,
     col = colors, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'E', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

colors <- colorRampPalette(c('white', palette[3]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(Rc.coord.pred.3, add = T, zlim = c(0,1.1), legend.mar = 6.1,
     col = colors, legend = T, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'F', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

par(mar=c(2.1,1.1,1.1,0.1))
plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Ignore')), col = '#222222', side = 1, line = -2.125, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Ignore'), col = palette[3], side = 1, line = -2.125, cex = 1, adj = 0)

par(mar=c(0,0,0,0))
colors <- colorRampPalette(c('white', palette[5]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(imported.coord.pred.5, zlim = c(0,1), add = T, legend.mar = 6.1,
     col = colors, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'G', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

colors <- colorRampPalette(c('white', palette[5]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(Rc.coord.pred.5, add = T, zlim = c(0,1.1), legend.mar = 6.1,
     col = colors, legend = T, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'H', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

par(mar=c(2.1,1.1,1.1,0.1))
plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[5], side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[5], side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = -2.125, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[5], side = 1, line = -2.125, cex = 1, adj = 0)

par(mar=c(0,0,0,0))
colors <- colorRampPalette(c('white', palette[4]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(imported.coord.pred.4, zlim = c(0,1), add = T, legend.mar = 6.1,
     col = colors, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'I', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

colors <- colorRampPalette(c('white', palette[4]))(11)
plot(map.swazi.reproject, border = '#222222')
plot(Rc.coord.pred.4, add = T, zlim = c(0,1.1), legend.mar = 6.1,
     col = colors, legend = T, legend.width = 1.4)
plot(map.swazi.reproject, border = '#222222', add = T)
mtext(text = 'J', side = 3, line = -0.5, cex = 1, font = 1, col = '#222222', at = extent(map.swazi.reproject)@xmin)

par(mar=c(2.1,1.1,1.1,0.1))
plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[4], side = 1, line = -6.125, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[4], side = 1, line = -4.125, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = -2.125, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[4], side = 1, line = -2.125, cex = 1, adj = 0)
dev.off()

