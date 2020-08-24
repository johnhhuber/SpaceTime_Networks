# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# specify paths in and out 
path.in = '../../data/eswatini_inference/'
path.out = '../../output/figs/'

# load likelihood profile data
load(list.files(path = path.in, pattern = 'validation_likelihood_profile.RData', full.names = T))

# load the true diffusion coefficients 
files.validation = list.files(path.in, pattern = 'validation_parameters', recursive = T, full.names = T)

diffusion.coef.true = rep(NA, length(files.validation))
for(ii in 1:length(files.validation))
{
  load(files.validation[ii])
  diffusion.coef.true[ii] = median.diffusion.coef
}

# load the estimated diffusion coefficient 
files.posterior = list.files(path.in, pattern = 'scalars_pooled', recursive = T, full.names = T)
files.posterior = files.posterior[grepl('validation', files.posterior)]

scalars = list()
for(ii in 1:length(files.posterior))
{
  scalars[[ii]] <- read.csv(files.posterior[ii], header = T)
}

# get the mle estimates of the diffusion coefficients
diffusion.default.mle <- diffusion.coef[which(-likelihood.default.true.tau == max(-likelihood.default.true.tau[-c(1:5)]))]
diffusion.bth.mle <- diffusion.coef[which(-likelihood.bth == max(-likelihood.bth[-c(1:5)]))]
diffusion.ignore.mle <- diffusion.coef[which(-likelihood.ignore == max(-likelihood.ignore[-c(1:5)]))]

diffusion.default.bounds <- diffusion.coef[which(-likelihood.default.true.tau >= (max(-likelihood.default.true.tau[-c(1:5)]) - 10))][-1]
diffusion.bth.bounds <- diffusion.coef[which(-likelihood.bth >= (max(-likelihood.bth[-c(1:5)]) - 10))][-1]
diffusion.ignore.bounds <- diffusion.coef[which(-likelihood.ignore >= (max(-likelihood.ignore[-c(1:5)]) - 10))][-1]

# get diffusion coefficient median estimates
diffusion.coef.estim <- lapply(1:length(scalars), function(x){return(quantile(scalars[[x]]$diffusion_coef, probs = c(0.025, 0.5, 0.975)))})

# create the log likelihood ratio 
ll.ratio.default <- min(likelihood.default.true.tau[-1]) - likelihood.default.true.tau

# specify color palette 
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = paste(path.out, 'S1_Fig.tiff', sep = ''), width = 6.5, height = 5, units = 'in', res = 300, compression = 'lzw')

par(mar=c(3.1,6.3,1.1,2.1))
layout(mat = matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = T))
plot((diffusion.coef[-c(1:5)]), -likelihood.default.true.tau[-c(1:5)], type = 'n', bty = 'n', axes = F, xlab = '', ylab = '', yaxs = 'i',
     lwd = 2, col = palette[2], xlim = c(0,10), ylim = c(-max(likelihood.default.true.tau[-c(1:5)]), -0.9999 * min(likelihood.default.true.tau[-c(1:5)])), xaxs = 'i')

rect(xleft = head(diffusion.default.bounds, n = 1), xright = tail(diffusion.default.bounds, n = 1),
     ybottom = -max(likelihood.default.true.tau[-c(1:5)]), ytop = -0.9999 * min(likelihood.default.true.tau[-c(1:5)]),
     col = 'lightgrey', border = NA)
rect(xleft = head(diffusion.default.mle, n = 1), xright = tail(diffusion.default.mle, n = 1),
     ybottom = -max(likelihood.default.true.tau[-c(1:5)]), ytop = -0.9999 * min(likelihood.default.true.tau[-c(1:5)]),
     col = 'darkgrey', border = NA)

lines((diffusion.coef[-c(1:5)]), -likelihood.default.true.tau[-c(1:5)], col = palette[2], lwd = 2)
abline(v = (diffusion.coef.true[2]), col = '#222222', lwd = 2, lty = 1)

axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 2, 'Log Likelihood', line = 5, cex = 0.75)

plot((diffusion.coef[-c(1:5)]), -likelihood.bth[-c(1:5)], type = 'n', bty = 'n', axes = F, xlab = '', ylab = '',
     lwd = 2, col = palette[1], xlim = c(0,10), yaxs = 'i', xaxs = 'i', ylim = c(-max(likelihood.bth[-c(1:5)]), -0.99975 * min(likelihood.bth[-c(1:5)])))

rect(xleft = head(diffusion.bth.bounds, n = 1), xright = tail(diffusion.bth.bounds, n = 1),
     ybottom = -max(likelihood.bth[-c(1:5)]), ytop = -0.99975 * min(likelihood.bth[-c(1:5)]),
     col = 'lightgrey', border = NA)
rect(xleft = head(diffusion.bth.mle, n = 1), xright = tail(diffusion.bth.mle, n = 1),
     ybottom = -max(likelihood.bth[-c(1:5)]), ytop = -0.99975 * min(likelihood.bth[-c(1:5)]),
     col = 'darkgrey', border = NA)
segments(x0 = head(diffusion.bth.mle, n = 1), x1 = tail(diffusion.bth.mle, n = 1),
         y0 = -max(likelihood.bth[-c(1:5)]), y1 = -0.99975 * min(likelihood.bth[-c(1:5)]),
         col = 'darkgrey', lwd = 2)

lines((diffusion.coef[-c(1:5)]), -likelihood.bth[-c(1:5)], col = palette[1], lwd = 2)
abline(v = (diffusion.coef.true[1]), col = '#222222', lwd = 2, lty = 1)
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 2, 'Log Likelihood', line = 5, cex = 0.75)

par(mar=c(4.1,6.3,1.1,2.1))
plot((diffusion.coef[-c(1:5)]), -likelihood.ignore[-c(1:5)], type = 'n', bty = 'n', axes = F, xlab = '', ylab = '', yaxs = 'i',
     lwd = 2, col = palette[3], xlim = c(0,10), ylim = c(-max(likelihood.ignore[-c(1:5)]), -0.9999 * min(likelihood.ignore[-c(1:5)])), xaxs = 'i')

rect(xleft = head(diffusion.ignore.bounds, n = 1), xright = tail(diffusion.ignore.bounds, n = 1),
     ybottom = -max(likelihood.ignore[-c(1:5)]), ytop = -0.9999 * min(likelihood.ignore[-c(1:5)]),
     col = 'lightgrey', border = NA)
rect(xleft = head(diffusion.ignore.mle, n = 1), xright = tail(diffusion.ignore.mle, n = 1),
     ybottom = -max(likelihood.ignore[-c(1:5)]), ytop = -0.9999 * min(likelihood.ignore[-c(1:5)]),
     col = 'darkgrey', border = NA)

lines((diffusion.coef[-c(1:5)]), -likelihood.ignore[-c(1:5)], col = palette[3], lwd = 2)
abline(v = (diffusion.coef.true[3]), col = '#222222', lwd = 2, lty = 1)
axis(side = 1)
axis(side = 2, las = 1)
mtext(side = 2, 'Log Likelihood', line = 5, cex = 0.75)
mtext(side = 1, expression(paste('Diffusion Coefficient (D)', sep = '')), line = 3, cex = 0.75)

dev.off()
