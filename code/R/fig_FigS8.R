# install necessary packages
if(!require(rcartocolor)){install.packages('rcartocolor'); library(rcartocolor)}

# load necessary data
load('../../data/figs/fig_7/fig7.RData')

# specify color palette
palette <- carto_pal(n = 10, name = 'Bold')

# plot 
tiff(filename = '../../output/figs/S8_Fig.tiff', width = 9, height = 9, units = 'in', res = 300, compression = 'lzw')

layout(mat = matrix(1:30, nrow = 5, ncol = 6, byrow = T))
par(mar = c(3.1,3.5,2.3,1.1))

plot(df.stt$max_founder_date, df.stt$founder_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Imported Classification', side = 2, line = 2.3, cex = 1)
mtext('Max Founder Date', side = 3, line = 0, cex = 1)

plot(df.stt$diffusion_coef, df.stt$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('D', side = 3, line = 0, cex = 1)

plot(df.stt$tau_s, df.stt$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')
mtext(expression(tau['s']), side = 3, line = 0, cex = 1)

plot(df.stt$tau_l, df.stt$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')
mtext(expression(tau['l']), side = 3, line = 0, cex = 1)

plot(df.stt$frac_cluster, df.stt$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Clusteredness', side = 3, line = 0, cex = 1)

plot(df.stt$founder_prob, df.stt$founder_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Prop Imported', side = 3, line = 0, cex = 1)


plot(df.stt$max_founder_date, df.stt$local_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Local Classification', side = 2, line = 2.3, cex = 1)

plot(df.stt$diffusion_coef, df.stt$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_s, df.stt$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_l, df.stt$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$frac_cluster, df.stt$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$founder_prob, df.stt$local_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')


plot(df.stt$max_founder_date, df.stt$classification_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Case Classification', side = 2, line = 2.3, cex = 1)

plot(df.stt$diffusion_coef, df.stt$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_s, df.stt$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_l, df.stt$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$frac_cluster, df.stt$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$founder_prob, df.stt$classification_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')


plot(df.stt$max_founder_date, df.stt$parent_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Transmission Linkage', side = 2, line = 2.3, cex = 1)

plot(df.stt$diffusion_coef, df.stt$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_s, df.stt$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_l, df.stt$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$frac_cluster, df.stt$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$founder_prob, df.stt$parent_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')


plot(df.stt$max_founder_date, df.stt$outbreak_acc, bty = 'n', cex  = 0.5, pch = 16, 
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$max_founder_date)), ylim = c(-0.01,1.01),
     col = palette[2], xlab = '', ylab = '')
mtext('Outbreak', side = 2, line = 2.3, cex = 1)

plot(df.stt$diffusion_coef, df.stt$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$diffusion_coef)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_s, df.stt$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_s)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$tau_l, df.stt$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$tau_l)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$frac_cluster, df.stt$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$frac_cluster)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

plot(df.stt$founder_prob, df.stt$outbreak_acc, bty = 'n', cex = 0.5, pch = 16,
     las = 1, xaxs = 'i', yaxs = 'i', xlim = c(0, max(df.stt$founder_prob)), ylim = c(-0.01, 1.01),
     col = palette[2], xlab = '', ylab = '')

dev.off()
