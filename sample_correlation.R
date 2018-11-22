data <- read.csv('all.csv')
data <- data[,2:11]

library(Hmisc)
res <- rcorr(as.matrix(data))
pdf('Correlation_plot.pdf',onefile = T)
library(pheatmap)
pheatmap(res$r,
         col = colorRampPalette(c("#3399FF","white","red"))(100),
         cluster_rows = T,cluster_cols=T,
         clustering_distance_rows = 'correlation',clustering_distance_cols = 'correlation',
         treeheight_row = 80,treeheight_col = 80,
         display_numbers = T,number_color = 'black',fontsize_number = 10,
         border_color = NA,legend=TRUE,
         show_colnames=T,show_rownames=T,
         fontsize_col = 10, 
         fontsize_row = 10,fontsize = 10,
         main = 'Correlation')

library(corrplot)
corrplot(res$r,col = colorRampPalette(c("blue","white","red"))(100),
         type='full',order='hclust',tl.col = "black", 
         diag=T,outline = T,bg = 'white',addgrid.col = 'black',
         addCoef.col = 'black',tl.cex = 1,
         p.mat = res$P,sig.level = 0.01,insig = 'blank',
         number.digits = 2)

symnum(res$r)

dev.off()
