# install necessary packages
if(!require(igraph)){install.packages('igraph'); library(igraph)}
if(!require(grDevices)){install.packages('grDevices'); library(grDevices)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# specify paths in and out 
path.in = '../../data/figs/fig_1/'
path.out = '../../output/figs/'

# load all of the validation metrics and store in list
files = list.files(path.in, pattern = 'validation_metrics.RData', recursive = T)
files.10local = files[grep('10local', files)]
files.50local = files[grep('50local', files)]
files.95local = files[grep('95local', files)]

assump.edge = list()
assump.edge[[1]] = lapply(1:length(files.10local), function(x){load(file.path(path.in, files.10local[x])); return(quantile(accuracy$parent.acc, probs = c(0.025, 0.50, 0.975)))})
assump.edge[[2]] = lapply(1:length(files.50local), function(x){load(file.path(path.in, files.50local[x])); return(quantile(accuracy$parent.acc, probs = c(0.025, 0.50, 0.975)))})
assump.edge[[3]] = lapply(1:length(files.95local), function(x){load(file.path(path.in, files.95local[x])); return(quantile(accuracy$parent.acc, probs = c(0.025, 0.50, 0.975)))})

assump.founder = list()
assump.founder[[1]] = lapply(1:length(files.10local), function(x){load(file.path(path.in, files.10local[x])); return(quantile(accuracy$founder.acc, probs = c(0.025, 0.50, 0.975)))})
assump.founder[[2]] = lapply(1:length(files.50local), function(x){load(file.path(path.in, files.50local[x])); return(quantile(accuracy$founder.acc, probs = c(0.025, 0.50, 0.975)))})
assump.founder[[3]] = lapply(1:length(files.95local), function(x){load(file.path(path.in, files.95local[x])); return(quantile(accuracy$founder.acc, probs = c(0.025, 0.50, 0.975)))})

assump.outbreak = list()
assump.outbreak[[1]] = lapply(1:length(files.10local), function(x){load(file.path(path.in, files.10local[x])); return(quantile(accuracy$outbreak.acc, probs = c(0.025, 0.50, 0.975)))})
assump.outbreak[[2]] = lapply(1:length(files.50local), function(x){load(file.path(path.in, files.50local[x])); return(quantile(accuracy$outbreak.acc, probs = c(0.025, 0.50, 0.975)))})
assump.outbreak[[3]] = lapply(1:length(files.95local), function(x){load(file.path(path.in, files.95local[x])); return(quantile(accuracy$outbreak.acc, probs = c(0.025, 0.50, 0.975)))})

Rc = list()
Rc[[1]] = lapply(1:length(files.10local), function(x){load(file.path(path.in, files.10local[x])); quantile(unlist(lapply(1:length(statistics), function(x){rowMeans(statistics[[x]]$R0.distribution)})),
                                                                                                           probs = c(0.025, 0.50, 0.975))})
Rc[[2]] = lapply(1:length(files.50local), function(x){load(file.path(path.in, files.50local[x])); quantile(unlist(lapply(1:length(statistics), function(x){rowMeans(statistics[[x]]$R0.distribution)})),
                                                                                                           probs = c(0.025, 0.50, 0.975))})
Rc[[3]] = lapply(1:length(files.95local), function(x){load(file.path(path.in, files.95local[x])); quantile(unlist(lapply(1:length(statistics), function(x){rowMeans(statistics[[x]]$R0.distribution)})),
                                                                                                           probs = c(0.025, 0.50, 0.975))})

# load all networks and compute adjacency matrices
files.networks = list.files(path.in, pattern = 'network_true_1_FULL.csv', recursive = T)
files.nodes = list.files(path.in, pattern = 'nodes_1_FULL_adjusted.csv', recursive = T)

nodes = lapply(1:length(files.nodes), function(x){n = read.csv(file.path(path.in, files.nodes[x]), header = T); return(n)})
graphs = list()
Rc.true = rep(NA, length(files.networks))

for(ii in 1:length(files.networks))
{
  chain = readLines(file.path(path.in, files.networks[ii]))
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
  
  g = graph.adjacency(networks[-1,-1]>0,mode='directed')
  graphs[[ii]] = g
}

# specify color palette
palette <- rcartocolor::carto_pal(n = 10, name = 'Bold')

tiff(filename = file.path(path.out, 'Fig1.tiff'), width = 7.5, height = 5, units = 'in', res = 300, compression = 'lzw')
layout(mat=matrix(c(rep(4,5), rep(5,5), rep(6,5),
                    rep(1,5), rep(2,5), rep(3,5),
                    rep(1,5), rep(2,5), rep(3,5),
                    rep(1,5), rep(2,5), rep(3,5),
                    rep(7,3), rep(8,3), rep(9,3),
                    rep(10,3), rep(11,3)),
                    nrow = 15, ncol = 5, byrow = F))
par(mar=c(3.1,4.2,0.3,4.1))
plot(0:4, 0:4, xlim = c(0.5, 4.5), ylim = c(0,1.1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
segments(x0 = 0.8, y0 = assump.founder[[1]][[2]][1], y1 = assump.founder[[1]][[2]][3], lwd = 2, col = palette[2])
points(0.8, assump.founder[[1]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 0.9, y0 = assump.founder[[1]][[1]][1], y1 = assump.founder[[1]][[1]][3], lwd = 2, col = palette[1])
points(0.9, assump.founder[[1]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 1.0, y0 = assump.founder[[1]][[3]][1], y1 = assump.founder[[1]][[3]][3], lwd = 2, col = palette[3])
points(1.0, assump.founder[[1]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 1.1, y0 = assump.founder[[1]][[5]][1], y1 = assump.founder[[1]][[5]][3], lwd = 2, col = palette[5])
points(1.1, assump.founder[[1]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 1.2, y0 = assump.founder[[1]][[4]][1], y1 = assump.founder[[1]][[4]][3], lwd = 2, col = palette[4])
points(1.2, assump.founder[[1]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 1.8, y0 = assump.edge[[1]][[2]][1], y1 = assump.edge[[1]][[2]][3], lwd = 2, col = palette[2])
points(1.8, assump.edge[[1]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 1.9, y0 = assump.edge[[1]][[1]][1], y1 = assump.edge[[1]][[1]][3], lwd = 2, col = palette[1])
points(1.9, assump.edge[[1]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 2.0, y0 = assump.edge[[1]][[3]][1], y1 = assump.edge[[1]][[3]][3], lwd = 2, col = palette[3])
points(2.0, assump.edge[[1]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 2.1, y0 = assump.edge[[1]][[5]][1], y1 = assump.edge[[1]][[5]][3], lwd = 2, col = palette[5])
points(2.1, assump.edge[[1]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 2.2, y0 = assump.edge[[1]][[4]][1], y1 = assump.edge[[1]][[4]][3], lwd = 2, col = palette[4])
points(2.2, assump.edge[[1]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 2.8, y0 = assump.outbreak[[1]][[2]][1], y1 = assump.outbreak[[1]][[2]][3], lwd = 2, col = palette[2])
points(2.8, assump.outbreak[[1]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 2.9, y0 = assump.outbreak[[1]][[1]][1], y1 = assump.outbreak[[1]][[1]][3], lwd = 2, col = palette[1])
points(2.9, assump.outbreak[[1]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 3.0, y0 = assump.outbreak[[1]][[3]][1], y1 = assump.outbreak[[1]][[3]][3], lwd = 2, col = palette[3])
points(3.0, assump.outbreak[[1]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 3.1, y0 = assump.outbreak[[1]][[5]][1], y1 = assump.outbreak[[1]][[5]][3], lwd = 2, col = palette[5])
points(3.1, assump.outbreak[[1]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 3.2, y0 = assump.outbreak[[1]][[4]][1], y1 = assump.outbreak[[1]][[4]][3], lwd = 2, col = palette[4])
points(3.2, assump.outbreak[[1]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 3.7, x1 = 4.3, y0 = Rc.true[1], lwd = 2, lty = 1, col = 'gray60')
segments(x0 = 3.8, y0 = Rc[[1]][[2]][1], y1 = Rc[[1]][[2]][3], lwd = 2, col = palette[2])
points(3.8, Rc[[1]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 3.9, y0 = Rc[[1]][[1]][1], y1 = Rc[[1]][[1]][3], lwd = 2, col = palette[1])
points(3.9, Rc[[1]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 4.0, y0 = Rc[[1]][[3]][1], y1 = Rc[[1]][[3]][3], lwd = 2, col = palette[3])
points(4.0, Rc[[1]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 4.1, y0 = Rc[[1]][[5]][1], y1 = Rc[[1]][[5]][3], lwd = 2, col = palette[5])
points(4.1, Rc[[1]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 4.2, y0 = Rc[[1]][[4]][1], y1 = Rc[[1]][[4]][3], lwd = 2, col = palette[4])
points(4.2, Rc[[1]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

axis(side = 4, col.axis = '#222222', las = 1)
axis(side = 2, col.axis = '#222222', las = 1)
axis(side = 1, col.axis = '#222222', labels = FALSE, lwd.ticks = 0, at = c(0.5, 4.5))
mtext('Accuracy', side = 2, line = 3)
mtext(expression('R'['c']), side = 4, line = 3)
mtext('A', side = 3, line = -1, at = 0.5, font = 1)

par(mar=c(3.1,4.2,0.3,4.1))
plot(0:4, 0:4, xlim = c(0.5, 4.5), ylim = c(0,1.1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
segments(x0 = 0.8, y0 = assump.founder[[2]][[2]][1], y1 = assump.founder[[2]][[2]][3], lwd = 2, col = palette[2])
points(0.8, assump.founder[[2]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 0.9, y0 = assump.founder[[2]][[1]][1], y1 = assump.founder[[2]][[1]][3], lwd = 2, col = palette[1])
points(0.9, assump.founder[[2]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 1.0, y0 = assump.founder[[2]][[3]][1], y1 = assump.founder[[2]][[3]][3], lwd = 2, col = palette[3])
points(1.0, assump.founder[[2]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 1.1, y0 = assump.founder[[2]][[5]][1], y1 = assump.founder[[2]][[5]][3], lwd = 2, col = palette[5])
points(1.1, assump.founder[[2]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 1.2, y0 = assump.founder[[2]][[4]][1], y1 = assump.founder[[2]][[4]][3], lwd = 2, col = palette[4])
points(1.2, assump.founder[[2]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 1.8, y0 = assump.edge[[2]][[2]][1], y1 = assump.edge[[2]][[2]][3], lwd = 2, col = palette[2])
points(1.8, assump.edge[[2]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 1.9, y0 = assump.edge[[2]][[1]][1], y1 = assump.edge[[2]][[1]][3], lwd = 2, col = palette[1])
points(1.9, assump.edge[[2]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 2.0, y0 = assump.edge[[2]][[3]][1], y1 = assump.edge[[2]][[3]][3], lwd = 2, col = palette[3])
points(2.0, assump.edge[[2]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 2.1, y0 = assump.edge[[2]][[5]][1], y1 = assump.edge[[2]][[5]][3], lwd = 2, col = palette[5])
points(2.1, assump.edge[[2]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 2.2, y0 = assump.edge[[2]][[4]][1], y1 = assump.edge[[2]][[4]][3], lwd = 2, col = palette[4])
points(2.2, assump.edge[[2]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 2.8, y0 = assump.outbreak[[2]][[2]][1], y1 = assump.outbreak[[2]][[2]][3], lwd = 2, col = palette[2])
points(2.8, assump.outbreak[[2]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 2.9, y0 = assump.outbreak[[2]][[1]][1], y1 = assump.outbreak[[2]][[1]][3], lwd = 2, col = palette[1])
points(2.9, assump.outbreak[[2]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 3.0, y0 = assump.outbreak[[2]][[3]][1], y1 = assump.outbreak[[2]][[3]][3], lwd = 2, col = palette[3])
points(3.0, assump.outbreak[[2]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 3.1, y0 = assump.outbreak[[2]][[5]][1], y1 = assump.outbreak[[2]][[5]][3], lwd = 2, col = palette[5])
points(3.1, assump.outbreak[[2]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 3.2, y0 = assump.outbreak[[2]][[4]][1], y1 = assump.outbreak[[2]][[4]][3], lwd = 2, col = palette[4])
points(3.2, assump.outbreak[[2]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 3.7, x1 = 4.3, y0 = Rc.true[2], lwd = 2, lty = 1, col = 'gray60')
segments(x0 = 3.8, y0 = Rc[[2]][[2]][1], y1 = Rc[[2]][[2]][3], lwd = 2, col = palette[2])
points(3.8, Rc[[2]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 3.9, y0 = Rc[[2]][[1]][1], y1 = Rc[[2]][[1]][3], lwd = 2, col = palette[1])
points(3.9, Rc[[2]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 4.0, y0 = Rc[[2]][[3]][1], y1 = Rc[[2]][[3]][3], lwd = 2, col = palette[3])
points(4.0, Rc[[2]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 4.1, y0 = Rc[[2]][[5]][1], y1 = Rc[[2]][[5]][3], lwd = 2, col = palette[5])
points(4.1, Rc[[2]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 4.2, y0 = Rc[[2]][[4]][1], y1 = Rc[[2]][[4]][3], lwd = 2, col = palette[4])
points(4.2, Rc[[2]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

axis(side = 4, col.axis = '#222222', las = 1)
axis(side = 2, col.axis = '#222222', las = 1)
axis(side = 1, col.axis = '#222222', labels = FALSE, lwd.ticks = 0, at = c(0.5,4.5))
mtext('Accuracy', side = 2, line = 3)
mtext(expression('R'['c']), side = 4, line = 3)
mtext('B', side = 3, line = -1, at = 0.5, font = 1)

par(mar=c(3.3,4.2,0.3,4.1))
plot(0:4, 0:4, xlim = c(0.5,4.5), ylim = c(0,1.1), type = 'n', axes = FALSE, xlab = '', ylab = '', main = '', xaxt = 'n', xaxs = 'i')
segments(x0 = 0.8, y0 = assump.founder[[3]][[2]][1], y1 = assump.founder[[3]][[2]][3], lwd = 2, col = palette[2])
points(0.8, assump.founder[[3]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 0.9, y0 = assump.founder[[3]][[1]][1], y1 = assump.founder[[3]][[1]][3], lwd = 2, col = palette[1])
points(0.9, assump.founder[[3]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 1.0, y0 = assump.founder[[3]][[3]][1], y1 = assump.founder[[3]][[3]][3], lwd = 2, col = palette[3])
points(1.0, assump.founder[[3]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 1.1, y0 = assump.founder[[3]][[5]][1], y1 = assump.founder[[3]][[5]][3], lwd = 2, col = palette[5])
points(1.1, assump.founder[[3]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 1.2, y0 = assump.founder[[3]][[4]][1], y1 = assump.founder[[3]][[4]][3], lwd = 2, col = palette[4])
points(1.2, assump.founder[[3]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 1.8, y0 = assump.edge[[3]][[2]][1], y1 = assump.edge[[3]][[2]][3], lwd = 2, col = palette[2])
points(1.8, assump.edge[[3]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 1.9, y0 = assump.edge[[3]][[1]][1], y1 = assump.edge[[3]][[1]][3], lwd = 2, col = palette[1])
points(1.9, assump.edge[[3]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 2.0, y0 = assump.edge[[3]][[3]][1], y1 = assump.edge[[3]][[3]][3], lwd = 2, col = palette[3])
points(2.0, assump.edge[[3]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 2.1, y0 = assump.edge[[3]][[5]][1], y1 = assump.edge[[3]][[5]][3], lwd = 2, col = palette[5])
points(2.1, assump.edge[[3]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 2.2, y0 = assump.edge[[3]][[4]][1], y1 = assump.edge[[3]][[4]][3], lwd = 2, col = palette[4])
points(2.2, assump.edge[[3]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 2.8, y0 = assump.outbreak[[3]][[2]][1], y1 = assump.outbreak[[3]][[2]][3], lwd = 2, col = palette[2])
points(2.8, assump.outbreak[[3]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 2.9, y0 = assump.outbreak[[3]][[1]][1], y1 = assump.outbreak[[3]][[1]][3], lwd = 2, col = palette[1])
points(2.9, assump.outbreak[[3]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 3.0, y0 = assump.outbreak[[3]][[3]][1], y1 = assump.outbreak[[3]][[3]][3], lwd = 2, col = palette[3])
points(3.0, assump.outbreak[[3]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 3.1, y0 = assump.outbreak[[3]][[5]][1], y1 = assump.outbreak[[3]][[5]][3], lwd = 2, col = palette[5])
points(3.1, assump.outbreak[[3]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 3.2, y0 = assump.outbreak[[3]][[4]][1], y1 = assump.outbreak[[3]][[4]][3], lwd = 2, col = palette[4])
points(3.2, assump.outbreak[[3]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 3.7, x1 = 4.3, y0 = Rc.true[3], lwd = 2, lty = 1, col = 'gray60')
segments(x0 = 3.8, y0 = Rc[[3]][[2]][1], y1 = Rc[[3]][[2]][3], lwd = 2, col = palette[2])
points(3.8, Rc[[3]][[2]][2], pch = 15, col = palette[2], cex = 1.5)
segments(x0 = 3.9, y0 = Rc[[3]][[1]][1], y1 = Rc[[3]][[1]][3], lwd = 2, col = palette[1])
points(3.9, Rc[[3]][[1]][2], pch = 15, col = palette[1], cex = 1.5)
segments(x0 = 4.0, y0 = Rc[[3]][[3]][1], y1 = Rc[[3]][[3]][3], lwd = 2, col = palette[3])
points(4.0, Rc[[3]][[3]][2], pch = 15, col = palette[3], cex = 1.5)
segments(x0 = 4.1, y0 = Rc[[3]][[5]][1], y1 = Rc[[3]][[5]][3], lwd = 2, col = palette[5])
points(4.1, Rc[[3]][[5]][2], pch = 15, col = palette[5], cex = 1.5)
segments(x0 = 4.2, y0 = Rc[[3]][[4]][1], y1 = Rc[[3]][[4]][3], lwd = 2, col = palette[4])
points(4.2, Rc[[3]][[4]][2], pch = 15, col = palette[4], cex = 1.5)

segments(x0 = 3.70, x1 = 3.90, y0 = 0.4, lwd = 2, lty = 1, col = 'gray60')
text(x = 4.15, y = 0.4, 'True', cex = 1.25)

segments(x0 = 3.8, y0 = 0.1, y1 = 0.3, lwd = 2, col = adjustcolor('#222222',1.0))
points(3.8, 0.2, pch = 15, col = adjustcolor('#222222', 1.0), cex = 1.5)
text(x = 4.15, y = 0.2, '95% CI', cex = 1.25)

axis(side = 4, col.axis = '#222222', las = 1)
axis(side = 2, col.axis = '#222222', las = 1)
axis(side = 1, col.axis = '#222222', labels = FALSE, lwd.ticks = 0, at = c(0.5, 4.5))
mtext('Accuracy', side = 2, line = 3)
mtext(expression('R'['c']), side = 4, line = 3, col = '#222222')
mtext('C', side = 3, line = -1, at = 0.5, font = 1)

mtext('Case', side = 1, at = 1, line = 1, cex = 0.8)
mtext('Classification', side = 1, at = 1, line = 2.25, cex = 0.8)
mtext('Transmission', side = 1, at = 2, line = 1, cex = 0.8)
mtext('Linkage', side = 1, at = 2, line = 2.25, cex = 0.8)
mtext('Outbreak', side = 1, at = 3, line = 1, cex = 0.8)
mtext(expression('R'['c']), side = 1, at = 4, line = 1, col = '#222222',cex = 0.8)

par(mar=c(3,0,0,1))
plot(
  graphs[[1]],
  vertex.size=6,vertex.label=NA,
  vertex.color = '#222222',
  vertex.frame.color = '#222222',
  vertex.shape = 'circle',
  edge.arrow.size=0.2,
  edge.arrow.width = 1.5,
  edge.width = 1,
  edge.color = 'darkgrey',
  layout=layout.auto(graphs[[1]])
)

par(mar=c(3,0,0,1))
plot(
  graphs[[2]],
  vertex.size=6,vertex.label=NA,
  vertex.color = '#222222',
  vertex.frame.color = '#222222',
  vertex.shape = 'circle',
  edge.arrow.size=0.2,
  edge.arrow.width = 1.5,
  edge.width = 1,
  edge.color = 'darkgrey',
  layout=layout.auto(graphs[[2]])
)

par(mar=c(3,0,0,1))
plot(
  graphs[[3]],
  vertex.size=6,vertex.label=NA,
  vertex.color = '#222222',
  vertex.frame.color = '#222222',
  vertex.shape = 'circle',
  edge.arrow.size=0.2,
  edge.arrow.width = 1.5,
  edge.width = 1,
  edge.color = 'darkgrey',
  layout=layout.auto(graphs[[3]])
)

par(mar=c(2.1,1.1,2.5,0.1))

# data sources included
plot.new()
mtext(text = expression(underline('Inference Settings')), side = 3, line = 1, cex = 1, font = 1, adj = 0, col = '#222222')
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[2], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[2], side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = 0.75, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[2], side = 1, line = 0.75, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[1], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[1], side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = 0.75, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[1], side = 1, line = 0.75, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' Yes')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' Yes'), col = palette[3], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[3], side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Ignore')), col = '#222222', side = 1, line = 0.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Ignore'), col = palette[3], side = 1, line = 0.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[5], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[5], side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Estimate')), col = '#222222', side = 1, line = 0.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Estimate'), col = palette[5], side = 1, line = 0.25, cex = 1, adj = 0)

plot.new()
mtext(expression('Spatial:' * phantom(' No')), col = '#222222', side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression(phantom('Spatial:') * ' No'), col = palette[4], side = 1, line = -3.25, cex = 1, adj = 0)
mtext(expression('Temporal:' * phantom(' Yes')), col = '#222222', side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression(phantom('Temporal:') * ' Yes'), col = palette[4], side = 1, line = -1.25, cex = 1, adj = 0)
mtext(expression('Travel:' * phantom(' Believe')), col = '#222222', side = 1, line = 0.25, cex = 1, adj = 0)
mtext(expression(phantom('Travel:') * ' Believe'), col = palette[4], side = 1, line = 0.25, cex = 1, adj = 0)

dev.off()

