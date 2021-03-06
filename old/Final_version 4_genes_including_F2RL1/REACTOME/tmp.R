#这里的interactions是通路里面的基因，term_name是通路。
#把通路里面的基因转化为character类型
class(REACTOME$intersections)
REACTOME$intersections <- as.character(REACTOME$intersections)
#把所有的基因集合到一起
DEG <- c()
i=1
for (i in 1:20) {
  tmp <- unlist(strsplit(REACTOME$intersections[i],split = ','))
  DEG <- c(DEG,tmp)
  i=i+1
}
DEG
length(unique(DEG))
nrow(REACTOME)
cell_count <- length(unique(DEG))*nrow(REACTOME)
#创建矩阵，行列分别对应基因和通路
data <- matrix(0,36,20)
data <- as.data.frame(data)
colnames(data) <- REACTOME$term_name
rownames(data) <- unique(DEG)
head(data)
#新创建的矩阵，其行列分别为基因和通路名称，去原始矩阵里面检索该通路里面是否有该基因
j=1
for (j in 1:20) {
  i=1
  for (i in 1:36) {
    data[i,j]
    index <- which(as.character(REACTOME$term_name) == colnames(data)[j])
    tmp_gene <- unlist(strsplit(REACTOME$intersections,split = ',')[index])
    data[i,j] <- ifelse(rownames(data)[i] %in% tmp_gene,1,0)
    i=i+1
  }
  j=j+1
}
View(data)
data <- t(data)
#可视化
library(pheatmap)
library(RColorBrewer)
pdf("heatmap.pdf",width = 20,height = 5,onefile = T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
map <-pheatmap(data,col = c("white","#FF9999"),
               cluster_rows = F,cluster_cols=F,legend = F,
               border_color = 'grey',
               show_colnames=T,show_rownames=T,
               fontsize_col = 8, 
               fontsize_row = 8,fontsize = 8,
               main = 'Mapping Genes to Pathways')
dev.off()





