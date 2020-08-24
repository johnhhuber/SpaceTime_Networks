# install necessary packages
if(!require(grDevices)){install.packages('grDevices'); library(grDevices)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# specify paths in and out 
path.in = '../../data/eswatini_inference'
path.out = '../../output/figs/'

# load validation metrics 
files <- list.files(path = path.in, pattern = 'validation_metrics.RData', full.names = T, recursive = T)

assump.edge = list()
assump.founder = list()
assump.outbreak = list()
assump.classif.total = list()
assump.classif.local = list()
assump.classif.import = list()
Rc = list()

for(ii in 1:length(files))
{
  print(ii)
  load(file.path(files[ii]))
  assump.edge[[ii]] <- quantile(accuracy$edge.acc, probs = c(0.025, 0.50, 0.975))
  assump.founder[[ii]] <- quantile(accuracy$classif.acc, probs = c(0.025, 0.50, 0.975))
  assump.classif.total[[ii]] <- quantile(accuracy$classif.acc, probs = c(0.025, 0.50, 0.975))
  assump.classif.local[[ii]] <- quantile(accuracy$classif.local.acc, prob = c(0.025, 0.50, 0.975))
  assump.classif.import[[ii]] <- quantile(accuracy$classif.import.acc, prob = c(0.025, 0.50, 0.975))
  assump.outbreak[[ii]] <- quantile(accuracy$outbreak.acc, probs = c(0.025, 0.50, 0.975))
  Rc[[ii]] <- quantile(unlist(lapply(1:length(statistics), function(x){(statistics[[x]]$R0.distribution)})),
                       probs = c(0.025, 0.50, 0.975))
}

# load the true networks and get the correct Rc values 
files.networks = list.files(path.in, pattern = 'network_true_0_FULL.csv', recursive = T, full.names = T)
files.nodes = list.files(path.in, pattern = 'nodes_0_FULL.csv', recursive = T, full.names = T)

nodes = lapply(1:length(files.nodes), function(x){n = read.csv(file.path(path.in, files.nodes[x]), header = T); return(n)})
Rc.true = rep(NA, length(files.networks))

for(ii in 1:length(files.networks))
{
  chain = readLines(files.networks[ii])
  chain = strsplit(chain,';')
  numNodes = as.numeric(strsplit(tail(chain[[1]],1),'-')[[1]][3])
  networks = matrix(0, numNodes+1, numNodes+1)
  ancestor = descendant = character()
  k = ii.in = jj.in = numeric()
  
  for(jj in 1 : length(chain[[1]])){
    ancestor = strsplit(as.character(chain[[1]][jj]), '-')[[1]][1]
    k = strsplit(as.character(chain[[1]][jj]), '-')[[1]][2]
    descendant = strsplit(as.character(chain[[1]][jj]), '-')[[1]][3]
    
    if(ancestor == 's'){ii.in = 1}
    if(ancestor != 's'){ii.in = as.integer(ancestor) + 1}
    jj.in = as.integer(descendant) + 1
    
    networks[ii.in, jj.in] = as.integer(k)
  }
  
  Rc.true[ii] = mean(rowSums(networks[-1,-1]))
}

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot
tiff(filename = file.path(path.out, 'Fig6.tiff'), width = 12, height = 6 , units = 'in', res = 300, compression = 'lzw')
layout(mat = matrix(c(rep(1,5), 2, 2, 2, 2, 2), nrow = 2, ncol = 5, byrow = T))
par(mar=c(2.1,4.2,0.1,0.0))
plot(0:6, 0:6, xlim = c(0.5, 9), ylim = c(0,1.02), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
axis(side = 2, col.axis = '#222222', las = 1, cex.axis = 1.25)
mtext('Accuracy', side = 2, at = 0.5, line = 2.75, cex = 1.25)

rect(xleft = 0.5, xright = 1.75, ybottom = 0, ytop = 1, col = 'white', border = '#222222')
segments(x0 = 0.75, y0 = assump.founder[[2]][1], y1 = assump.founder[[2]][3], lwd = 2, col = palette[2])
points(0.75, assump.founder[[2]][2], pch = 15, col = palette[2], cex = 2)
segments(x0 = 1.0, y0 = assump.edge[[2]][1], y1 = assump.edge[[2]][3], lwd = 2, col = palette[2])
points(1.0, assump.edge[[2]][2], pch = 16, col = palette[2], cex = 2)
segments(x0 = 1.25, y0 = assump.outbreak[[2]][1], y1 = assump.outbreak[[2]][3], lwd = 2, col = palette[2])
points(1.25, assump.outbreak[[2]][2], pch = 17, col = palette[2], cex = 2)
segments(x0 = 1.375, x1 = 1.625, Rc.true[2], lwd = 2, col = 'gray60')
segments(x0 = 1.5, y0 = Rc[[2]][1], y1 = Rc[[2]][3], lwd = 2, col = palette[2])
points(1.5, Rc[[2]][2], pch = 23, col = palette[2], bg = palette[2], cex = 2)
mtext(text = 'A', side = 3, cex = 1, col = '#222222', at = 0.6, line = -2.5)

rect(xleft = 1.75, xright = 3, ybottom = 0, ytop = 1, col = 'white', border = '#222222')
segments(x0 = 2, y0 = assump.founder[[1]][1], y1 = assump.founder[[1]][3], lwd = 2, col = palette[1])
points(2, assump.founder[[1]][2], pch = 15, col = palette[1], cex = 2)
segments(x0 = 2.25, y0 = assump.edge[[1]][1], y1 = assump.edge[[1]][3], lwd = 2, col = palette[1])
points(2.25, assump.edge[[1]][2], pch = 16, col = palette[1], cex = 2)
segments(x0 = 2.5, y0 = assump.outbreak[[1]][1], y1 = assump.outbreak[[1]][3], lwd = 2, col = palette[1])
points(2.5, assump.outbreak[[1]][2], pch = 17, col = palette[1], cex = 2)
segments(x0 = 2.625, x1 = 2.875, y0 = Rc.true[1], lwd = 2, col = 'gray60')
segments(x0 = 2.75, y0 = Rc[[1]][1], y1 = Rc[[1]][3], lwd = 2, col = palette[1])
points(2.75, Rc[[1]][2], pch = 23, col = palette[1], bg = palette[1], cex = 2)
mtext(text = 'B', side = 3, cex = 1, col = '#222222', at = 1.85, line = -2.5)

rect(xleft = 3, xright = 4.25, ybottom = 0, ytop = 1, col = 'white', border = '#222222')
segments(x0 = 3.25, y0 = assump.founder[[3]][1], y1 = assump.founder[[3]][3], lwd = 2, col = palette[3])
points(3.25, assump.founder[[3]][2], pch = 15, col = palette[3], cex = 2)
segments(x0 = 3.5, y0 = assump.edge[[3]][1], y1 = assump.edge[[3]][3], lwd = 2, col = palette[3])
points(3.5, assump.edge[[3]][2], pch = 16, col = palette[3], cex = 2)
segments(x0 = 3.75, y0 = assump.outbreak[[3]][1], y1 = assump.outbreak[[3]][3], lwd = 2, col = palette[3])
points(3.75, assump.outbreak[[3]][2], pch = 17, col = palette[3], cex = 2)
segments(x0 = 3.875, x1 = 4.125, y0 = Rc.true[3], lwd = 2, col = 'gray60')
segments(x0 = 4, y0 = Rc[[3]][1], y1 = Rc[[3]][3], lwd = 2, col = palette[3])
points(4, Rc[[3]][2], pch = 23, col = palette[3], bg = palette[3], cex = 2)
mtext(text = 'C', side = 3, cex = 1, col = '#222222', at = 3.1, line = -2.5)

rect(xleft = 4.25, xright = 5.5, ybottom = 0, ytop = 1, col = 'white', border = '#222222')
segments(x0 = 4.5, y0 = assump.founder[[5]][1], y1 = assump.founder[[5]][3], lwd = 2, col = palette[5])
points(4.5, assump.founder[[5]][2], pch = 15, col = palette[5], cex = 2)
segments(x0 = 4.75, y0 = assump.edge[[5]][1], y1 = assump.edge[[5]][3], lwd = 2, col = palette[5])
points(4.75, assump.edge[[5]][2], pch = 16, col = palette[5], cex = 2)
segments(x0 = 5.0, y0 = assump.outbreak[[5]][1], y1 = assump.outbreak[[5]][3], lwd = 2, col = palette[5])
points(5.0, assump.outbreak[[5]][2], pch = 17, col = palette[5], cex = 2)
segments(x0 = 5.125, x1 = 5.375, y0 = Rc.true[5], lwd = 2, col = 'gray60')
segments(x0 = 5.25, y0 = Rc[[5]][1], y1 = Rc[[5]][3], lwd = 2, col = palette[5])
points(5.25, Rc[[5]][2], pch = 23, col = palette[5], bg = palette[5], cex = 2)
mtext(text = 'D', side = 3, cex = 1, col = '#222222', at = 4.35, line = -2.5)

rect(xleft = 5.5, xright = 6.75, ybottom = 0, ytop = 1, col = 'white', border = '#222222')
segments(x0 = 5.75, y0 = assump.founder[[4]][1], y1 = assump.founder[[4]][3], lwd = 2, col = palette[4])
points(5.75, assump.founder[[4]][2], pch = 15, col = palette[4], cex = 2)
segments(x0 = 6, y0 = assump.edge[[4]][1], y1 = assump.edge[[4]][3], lwd = 2, col = palette[4])
points(6, assump.edge[[4]][2], pch = 16, col = palette[4], cex = 2)
segments(x0 = 6.25, y0 = assump.outbreak[[4]][1], y1 = assump.outbreak[[4]][3], lwd = 2, col = palette[4])
points(6.25, assump.outbreak[[4]][2], pch = 17, col = palette[4], cex = 2)
segments(x0 = 6.375, x1 = 6.624, y0 = Rc.true[4], lwd = 2, col = 'gray60')
segments(x0 = 6.5, y0 = Rc[[4]][1], y1 = Rc[[4]][3], lwd = 2, col = palette[4])
points(6.5, Rc[[4]][2], pch = 23, col = palette[4], bg = palette[4], cex = 2)
mtext(text = 'E', side = 3, cex = 1, col = '#222222', at = 5.6, line = -2.5)

# add in axis for Rc 
axis(side = 4, col.axis = '#222222', las = 1, pos = 6.75, cex.axis = 1.25)
mtext(expression('R'['c']), side = 4, line = -19, cex = 1.25)

# add in legend
points(7.5, 0.75, pch = 15, col = '#222222', cex = 2)
text(x = 7.6, y = 0.75, 'Case Classification', cex = 1.5, pos = 4)
points(7.5, 0.65, pch = 16, col = '#222222', cex = 2)
text(x = 7.6, y = 0.65, 'Transmission Linkage', cex = 1.5, pos = 4)
points(7.5, 0.55, pch = 17, col = '#222222', cex = 2)
text(x = 7.6, y = 0.55, 'Outbreak', cex = 1.5, pos = 4)
points(7.5, 0.45, pch = 23, col = '#222222', bg = '#222222', cex = 2)
text(x = 7.6, y = 0.45, expression('R'['c']), cex = 1.5, pos = 4)

segments(x0 = 7.5, y0 = 0.30, y1 = 0.40, col = '#222222', lwd = 2)
text(x = 7.6, y = 0.35, '95% CI', cex = 1.5, pos = 4)

segments(x0 = 7.425, x1 = 7.575, y0 = 0.25, col = 'gray60', lwd = 2)
text(x = 7.6, y = 0.25, 'True', cex = 1.5, pos = 4)

par(mar=c(2.1,4.2,0.5,0.0))
plot(0:6, 0:6, xlim = c(0.5, 7), ylim = c(0,1.02), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
mtext(text = expression(underline('Inference Settings')), side = 1, at = 0.5, line = -22.75, cex = 1, font = 1, adj = 0, col = '#222222')
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 0.5, line = -21, cex = 1.0, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, at = 0.5, line = -21, cex = 1.0, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 0.5, line = -19, cex = 1.0, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, at = 0.5, line = -19, cex = 1.0, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, at = 0.5, line = -17, cex = 1.0, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, at = 0.5, line = -17, cex = 1.0, adj = 0)

mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 1.5, line = -21, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, at = 1.5, line = -21, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 1.5, line = -19, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, at = 1.5, line = -19, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, at = 1.5, line = -17, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, at = 1.5, line = -17, cex = 1, adj = 0)

mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, at = 2.5, line = -21, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, at = 2.5, line = -21, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 2.5, line = -19, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, at = 2.5, line = -19, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Exclude')), col = '#222222', side = 1, at = 2.5, line = -17, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Exclude'), col = palette[3], side = 1, at = 2.5, line = -17, cex = 1, adj = 0)

mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, at = 3.5, line = -21, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[5], side = 1, at = 3.5, line = -21, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 3.5, line = -19, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[5], side = 1, at = 3.5, line = -19, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, at = 3.5, line = -17, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[5], side = 1, at = 3.5, line = -17, cex = 1, adj = 0)

mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, at = 4.5, line = -21, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[4], side = 1, at = 4.5, line = -21, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, at = 4.5, line = -19, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[4], side = 1, at = 4.5, line = -19, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, at = 4.5, line = -17, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[4], side = 1, at = 4.5, line = -17, cex = 1, adj = 0)

dev.off()
