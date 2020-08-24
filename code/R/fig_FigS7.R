# install necessary packages
if(!require(seqinr)){install.packages('seqinr'); library(seqinr)}
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load nodes
nodes <- read.csv('../../data/eswatini_inference/space_time_travel/validation/nodes_0_FULL.csv')

# specify prior hyperparameters
tau.s.alpha.prior = 12
tau.s.beta.prior = 3

tau.l.alpha.prior = 3
tau.l.beta.prior = 12

# load posterior 
scalars <- read.csv('../../data/eswatini_inference/space_time_travel/validation/scalars_pooled.csv')
load('../../data/eswatini_inference/space_time_travel/validation/networks_processed.RData')

# get the imported probs 
imported.prob <- network.data[[1]]$edgeProbs[1,-1]

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot the tau_s 
tiff(filename = '../../output/figs/S7_Fig.tiff', width = 6.5, height = 5, units = 'in', res = 300, compresson = 'lzw')
x <- seq(from = 0, to = 1, by = 0.001)
plot(x, dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior), type = 'n',
     bty = 'n', las = 1, xaxs = 'i', yaxs = 'i', xlab = '', ylab = '', ylim = c(0, 1.75 * max(dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior))),
     axes = F)
polygon(c(x, rev(x)), c(rep(0, length(x)), rev(dbeta(x, shape1 = tau.s.alpha.prior, shape2 = tau.s.beta.prior))),
        col = col2alpha('darkgrey'), border = NA)
hist(scalars$tau_s, add = T, freq = F, col = col2alpha(palette[2], alpha = 0.75), border = 'white')
founders <- which(imported.prob > 0.25)

# calculate the posterior distribution for tau_s and tau_l 
tau.s.alpha.posterior <- tau.s.alpha.prior + sum(nodes$History[founders])
tau.s.beta.posterior <- tau.s.beta.prior + length(founders) - sum(nodes$History[founders])

# overlay plot 
lines(x, dbeta(x, shape1 = tau.s.alpha.posterior, shape2 = tau.s.beta.posterior),
      col = '#222222', lwd = 2)
abline(v = 0.5, lty = 2, lwd = 2)
axis(side = 1)
axis(side = 2, las = 1)
title(xlab = expression(tau['s']), ylab = 'Density')
legend('topleft', col = c('darkgrey', palette[2], '#222222', '#222222'), pch = c(15, 15, NA, NA),
       lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 2), legend = c('Prior', 'MC3 Posterior', 'Conjugate-Prior Posterior', 'True Value'),
       bty = 'n')
dev.off()
