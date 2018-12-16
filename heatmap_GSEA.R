b <- read.delim('GSEA_ranked_gene_list.txt',header = T,sep = '\t')
gene <- b[18338:18537,]


exprs <- read.csv('all.csv')
exprs <- exprs[,c(2:10)]

data_heatmap <- merge(exprs,gene,by.x = 'Other.Gene.ID',by.y = 'NAME',all = F)
data_heatmap <- data_heatmap[,1:9]
rownames(data_heatmap)<-data_heatmap$Other.Gene.ID
data_heatmap <- data_heatmap[,-1]
colnames(data_heatmap)<-c('FK2','FK4','FK5',paste0(rep('HN',5),1:5))

library(pheatmap)
library(RColorBrewer)
pdf("heatmap.pdf",width = 7,height = 10,onefile = T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
map <-pheatmap(data_heatmap,col = colorRampPalette(c("blue","white","red"))(100),
               scale = "row",
               cluster_rows = T,cluster_cols=T,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               border_color = NA,legend=TRUE,
               show_colnames=T,show_rownames=T,
               fontsize_col = 10, 
               fontsize_row = 3,fontsize = 8,
               main = 'TOP200(UP in HN)')
dev.off()
