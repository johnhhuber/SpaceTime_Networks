# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load nodes file
#nodes <- read.csv('../../data/eswatini_inference/nodes.csv')

# calculate number of nodes
#numNodes = nrow(nodes)
numNodes = 775

# specify prior hyperparameters
tau.s.alpha.prior = 12
tau.s.beta.prior = 3

tau.l.alpha.prior = 3
tau.l.beta.prior = 12

# load posterior scalars
scalars <- read.csv('../../data/eswatini_inference/space_time_travel/scalars_pooled.csv')
load('../../data/eswatini_inference/space_time_travel/networks_processed.RData')

# get posterior imported probability 
imported.prob <- network.data[[1]]$edgeProbs[1,-1]

# load median posterior estimates 
load('../../data/eswatini_inference/space_time_travel/validation/validation_parameters.RData')

# estimate median number of imported 
num.imported <- round(numNodes * median.founder.prob)

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = '../../output/figs/S5_Fig.tiff', width = 6.5, height = 5, units = 'in', res = 300, compression = 'lzw')

x <- seq(from = 0, to = 1, by = 0.001)
par(mfrow = c(2,1), mar = c(4.1, 4.1, 1.1, 1.1))

plot(x, dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior), type = 'n',
     bty = 'n', las = 1, xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', ylim = c(0, 1.75 * max(dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior))),
     axes = F)
polygon(c(x, rev(x)), c(rep(0, length(x)), rev(dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior))),
        col = col2alpha('darkgrey'), border = NA)
hist(scalars$tau_s, add = T, freq = F, col = col2alpha(palette[2], alpha = 0.75), border = 'white')
founders <- which(imported.prob > 0.25 & !is.na(nodes$History))
# calculate the posterior distribution for tau_s and tau_l 
tau.s.alpha.posterior <- tau.s.alpha.prior + sum(nodes$History[founders])
tau.s.beta.posterior <- tau.s.beta.prior + length(founders) - sum(nodes$History[founders])

# overlay plot 
lines(x, dbeta(x, shape1 = tau.s.alpha.posterior, shape2 = tau.s.beta.posterior),
      col = '#222222', lwd = 2)
axis(side = 1)
axis(side = 2, las = 1)
title(xlab = expression(tau['s']), ylab = 'Density')

# loop through and plot posterior for tau_l
x <- seq(from = 0, to = 1, by = 0.001)
plot(x, dbeta(x, shape1 = tau.l.alpha.prior, shape2 = tau.l.beta.prior), type = 'n',
     bty = 'n', las = 1, xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', ylim = c(0, 25), axes = F)
polygon(c(x, rev(x)), c(rep(0, length(x)), rev(dbeta(x, shape1 = tau.l.alpha.prior, shape2 = tau.l.beta.prior))),
        col = col2alpha('darkgrey'), border = NA)
hist(scalars$tau_l, add = T, freq = F, col = palette[2], border = 'white')
# sample founders 
locals <- which(imported.prob < 0.25 & !is.na(nodes$History))

# calculate the posterior distribution for tau_s and tau_l 
tau.l.alpha.posterior <- tau.l.alpha.prior + sum(nodes$History[locals], na.rm = T)
tau.l.beta.posterior <- tau.l.beta.prior + length(locals) - sum(nodes$History[locals], na.rm = T)

# overlay plot 
lines(x, dbeta(x, shape1 = tau.l.alpha.posterior, shape2 = tau.l.beta.posterior),
      col = '#222222', lwd = 2)
axis(side = 1)
axis(side = 2, las = 1)
title(xlab = expression(tau['l']), ylab = 'Density')
legend('topleft', pch = c(15, 15, NA), lwd = c(NA, NA, 2), lty = c(NA, NA, 1), col = c('darkgrey', palette[2], '#222222'),
       legend = c('Prior', 'MC3 Posterior', 'Conjugate-Prior Posterior'),
       bty = 'n')

dev.off()
