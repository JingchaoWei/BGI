#"2018-12-06 21:20:56 PST"
rm(list=ls())
data <- read.csv('all.csv')
data <- data[,2:11]#提取表达矩阵
colnames(data)
colnames(data) <- c("FK2","FK4","FK5","FK6","FK7",
                    "HN1","HN2","HN3","HN4","HN5")


library(Hmisc)
library(pheatmap)
res <- rcorr(as.matrix(data))
pdf('Correlation_plot.pdf',onefile = T)
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



#还可以通过层次聚类来查看样本之间的关系：
names_cluster <- data.frame(sample=colnames(data),group=c(rep('FK',5),rep('HN',5)))
#如果group是factor,要转化成character
names_cluster$group <- as.character(names_cluster$group)
#自己构建新名字
newnames <- unlist(strsplit(colnames(data),split = '.',fixed = T))
newnames = newnames[-which(newnames=='FPKM')]
names_cluster$new_names <- newnames
#查看对应关系是否正确，应该是正确的。
names_cluster
colnames(data)
#对矩阵的名字进行替换
colnames(data) <- names_cluster$new_names
colnames(data)
#然后进行聚类
distance <- dist(t(data),method="euclidean")
clusters <- hclust(distance,method = "complete")
pdf("sample_cluster.pdf")
par(mar = c(5, 5, 5, 3))#图的下，左，上，右到页面的边距
par(mgp = c(1.5, 0.5, 0))
plot(clusters)
dev.off()
