# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load necessary data
load('../../data/figs/fig_7/fig7.RData')

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = '../../output/figs/S10_Fig.tiff', width = 9, height = 9, units = 'in', res = 300, compression = 'lzw')

layout(mat = matrix(1:30, nrow = 5, ncol = 6, byrow = T))
par(mar = c(3.1,3.5,2.3,1.1))

plot(df.stn$max_founder_date, df.stn$founder_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Imported Classification', side = 2, line = 2.3, cex = 1)
mtext('Max Founder Date', side = 3, line = 0, cex = 1)

plot(df.stn$diffusion_coef, df.stn$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('D', side = 3, line = 0, cex = 1)

plot(df.stn$tau_s, df.stn$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')
mtext(expression(tau['s']), side = 3, line = 0, cex = 1)

plot(df.stn$tau_l, df.stn$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')
mtext(expression(tau['l']), side = 3, line = 0, cex = 1)

plot(df.stn$frac_cluster, df.stn$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Clusteredness', side = 3, line = 0, cex = 1)

plot(df.stn$founder_prob, df.stn$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Prop Imported', side = 3, line = 0, cex = 1)


plot(df.stn$max_founder_date, df.stn$local_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Local Classification', side = 2, line = 2.3, cex = 1)

plot(df.stn$diffusion_coef, df.stn$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_s, df.stn$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_l, df.stn$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$frac_cluster, df.stn$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$founder_prob, df.stn$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')


plot(df.stn$max_founder_date, df.stn$classification_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Case Classification', side = 2, line = 2.3, cex = 1)

plot(df.stn$diffusion_coef, df.stn$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_s, df.stn$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_l, df.stn$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$frac_cluster, df.stn$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$founder_prob, df.stn$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')


plot(df.stn$max_founder_date, df.stn$parent_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Transmission Linkage', side = 2, line = 2.3, cex = 1)

plot(df.stn$diffusion_coef, df.stn$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_s, df.stn$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_l, df.stn$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$frac_cluster, df.stn$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$founder_prob, df.stn$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')


plot(df.stn$max_founder_date, df.stn$outbreak_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[3], xlab = '', ylab = '')
mtext('Outbreak', side = 2, line = 2.3, cex = 1)

plot(df.stn$diffusion_coef, df.stn$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_s, df.stn$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$tau_l, df.stn$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$frac_cluster, df.stn$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

plot(df.stn$founder_prob, df.stn$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stn$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[3], xlab = '', ylab = '')

dev.off()