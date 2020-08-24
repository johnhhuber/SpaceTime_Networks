# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load necessary data
load('../../data/figs/fig_7/fig7.RData')

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = '../../output/figs/S9_Fig.tiff', width = 9, height = 9, units = 'in', res = 300, compression = 'lzw')

layout(mat = matrix(1:30, nrow = 5, ncol = 6, byrow = T))
par(mar = c(3.1,3.5,2.3,1.1))

plot(df.stb$max_founder_date, df.stb$founder_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Imported Classification', side = 2, line = 2.3, cex = 1)
mtext('Max Founder Date', side = 3, line = 0, cex = 1)

plot(df.stb$diffusion_coef, df.stb$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('D', side = 3, line = 0, cex = 1)

plot(df.stb$tau_s, df.stb$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')
mtext(expression(tau['s']), side = 3, line = 0, cex = 1)

plot(df.stb$tau_l, df.stb$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')
mtext(expression(tau['l']), side = 3, line = 0, cex = 1)

plot(df.stb$frac_cluster, df.stb$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Clusteredness', side = 3, line = 0, cex = 1)

plot(df.stb$founder_prob, df.stb$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Prop Imported', side = 3, line = 0, cex = 1)


plot(df.stb$max_founder_date, df.stb$local_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Local Classification', side = 2, line = 2.3, cex = 1)

plot(df.stb$diffusion_coef, df.stb$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_s, df.stb$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_l, df.stb$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$frac_cluster, df.stb$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$founder_prob, df.stb$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')


plot(df.stb$max_founder_date, df.stb$classification_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Case Classification', side = 2, line = 2.3, cex = 1)

plot(df.stb$diffusion_coef, df.stb$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_s, df.stb$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_l, df.stb$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$frac_cluster, df.stb$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$founder_prob, df.stb$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')


plot(df.stb$max_founder_date, df.stb$parent_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Transmission Linkage', side = 2, line = 2.3, cex = 1)

plot(df.stb$diffusion_coef, df.stb$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_s, df.stb$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_l, df.stb$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$frac_cluster, df.stb$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$founder_prob, df.stb$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')


plot(df.stb$max_founder_date, df.stb$outbreak_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[1], xlab = '', ylab = '')
mtext('Outbreak', side = 2, line = 2.3, cex = 1)

plot(df.stb$diffusion_coef, df.stb$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_s, df.stb$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$tau_l, df.stb$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$frac_cluster, df.stb$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

plot(df.stb$founder_prob, df.stb$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stb$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[1], xlab = '', ylab = '')

dev.off()
