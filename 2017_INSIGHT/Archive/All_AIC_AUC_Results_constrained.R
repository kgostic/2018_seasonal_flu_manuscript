rm(list = ls())
library(ComplexHeatmap)
library(circlize)


pdf('AIC_summary_constrained.pdf')
par(mfrow = c(1,2))
load('master.AIC.table.constrained.RData')
write.csv(master.AIC.table, file = 'Master_AIC_constrained.csv')
## All the normal baseline models win except simple baseline for symdur
pmat = master.AIC.table[c(1, 2, 5, 6, 3, 4),] # RE-order rows
pmat[5, 1:3] = pmat[5, 4:6]
pmat[6, 1:3] = pmat[6, 4:6]
pmat = pmat[, -c(4:6)]
Heatmap(pmat, col = colorRamp2(c(0, 2, 5, 10), c(colors()[636], colors()[27], colors()[99], "red")), cluster_rows = FALSE, cluster_columns = FALSE, heatmap_legend_param = list(title = "Delta AIC"), column_title = "Constrained", 
        column_title_gp = gpar(fontsize = 20, fontface = "bold"))
rm(master.AIC.table)
dev.off()

