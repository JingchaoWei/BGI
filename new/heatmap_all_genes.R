#准备矩阵
load('Data.Rdata')
exp_gene <- background_all[,2:11]

#修改名字，去掉多余的字符
library(stringr)
names <- unlist(str_split(colnames(exp_gene),pattern = '[.,]'))
names <- names[seq(from=1,to=19,by=2)]
colnames(exp_gene) <- names

#用pheatmap,这个最好用
library(pheatmap)
library(RColorBrewer)
pdf("heatmap.pdf",width = 8,height = 15,onefile = T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
map <-pheatmap(exp_gene,col = colorRampPalette(c("blue","white","red"))(100),
               scale = "row",
               cluster_rows = F,cluster_cols=T,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               border_color = NA,legend=TRUE,
               show_colnames=T,show_rownames=F,
               fontsize_col = 8, 
               fontsize = 8,
               main = 'Heatmap')
dev.off()


#用ggplot
library(ggplot2)
#进行聚类:先用dist计算距离（组与组之间的距离），默认是euclidean 距离，然后用hclust聚类，
#默认是最长距离法：complete 
dist1 <- dist(exp_gene)
hc_gene <- hclust(dist1,method = "complete")
rowInd <- hc_gene$order#将聚类后行的顺序存为rowInd
hc_sample <-hclust(dist(t(exp_gene)),method = "complete")  #对矩阵进行转置，对原本的列进行聚类
colInd<-hc_sample$order  #将聚类后列的顺序存为colInd
exp_gene <- exp_gene[rowInd,colInd] #将数据按照聚类结果重排行和列
#进行数据格式转化，用melt
library(reshape2)
exp_gene$gene <- rownames(exp_gene)
colnames(exp_gene)
heatmap_data <- melt(exp_gene,id.vars=541,measure.vars=1:540,
                     variable.name="patient",value.name="expression")
head(heatmap_data)
#也可以手动进行数据转化
dim(exp_gene)
name_genes <- rep(rownames(exp_gene),dim(exp_gene)[2])
name_patients <- vector()
value_expression <- vector()
for (i in 1:dim(exp_gene)[2]){
  name_patients <- c(name_patients, rep(colnames(exp_gene)[i],dim(exp_gene)[1]))
  value_expression <- c(value_expression,exp_gene[,i])
}
df_heatmap <- data.frame(genes = name_genes, patient=name_patients,
                         expression_level=value_expression)
#画图
ggplot(heatmap_data, aes(gene,patient)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradient2(low = "#436EEE", mid="white",high = "#EE4000",midpoint = 13.62823,aesthetics = "fill") +
  xlab("List of genes ") +
  ylab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")

#用heatmap.2
exp_gene <- as.matrix(exp_gene)
library(gplots)
pdf("heatmap.pdf",width = 5,height = 10)
heatmap.2(t(exp_gene),col = colorRampPalette(c("blue","white","red"))(100),
          density.info="none",
          hclustfun = function(c)hclust(c,method="complete"),
          distfun = function(c)dist(c,method = "euclidean"),
          keysize = 1,cexRow=0.5,cexCol = 1, 
          srtCol=0,adjCol = c(0.5,1),
          trace = "none",dendrogram="both",
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean))
dev.off()
#the following produce the same image
#col = colorRampPalette(c("blue","white","red"))(100) can be replaced by col = bluered(100)

