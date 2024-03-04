d <- read.table(paste0("input/data_density_dimensions=1_multi_True"), sep=',')
# hist(d[,1], 100, col = "yellow", border = rgb(0.8, 0.7, 0, 1), xlab = '', 
#      main = 'Type: Density - Regions: k=3', cex.main = 1)

pdf("dens1d.pdf", width = 3.3, height = 2.2)
par(mar = c(2, 3, 1, 0.2), mgp = c(2, 0.5, 0))
hist(d[,1], 100, col = "yellow", border = rgb(0.8, 0.7, 0, 1), xlab = '', 
     main = 'Type: Density - Regions: k=3', freq = FALSE,
     cex = 1.3, cex.main = 1, font.main = 1)
dev.off()

d <- read.table(paste0("input/data_aggr_dimensions=1_multi_True"), sep=',')
# plot(d[,1], d[,2], bty = 'l', xlab = '', ylab = 'target',
#      pch = 21, bg = rgb(0, 0, 1, 0.2, maxColorValue = 1), col = NA,
#      main = 'Type: Aggregate - Regions: k=3', cex = 1, cex.main = 1)
#  

d <- d[sample(nrow(d),500),]

pdf("aggr1d.pdf", width = 3.3, height = 2.2)
par(mar = c(2, 3, 1, 0.2), mgp = c(2, 0.5, 0))
plot(d[,1], d[,2], bty = 'l', xlab = '', ylab = 'Target variable', 
     pch = 8, bg = rgb(0, 0, 1, 0.2, maxColorValue = 1), col = rgb(0, 0, 1),
     main = 'Type: Aggregate - Regions: k=3', cex = 1.3, cex.main = 1, font.main = 1)

dev.off()
