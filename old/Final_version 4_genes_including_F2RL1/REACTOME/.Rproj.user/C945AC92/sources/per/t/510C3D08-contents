#柱状图：
rm(list=ls()) 
library(ggplot2) 
library(Cairo) 
REACTOME <- read.csv('reactome_file.csv')
REACTOME <- REACTOME[1:20,]
colnames(REACTOME)
png_path="./REACTOME.png" 
CairoPNG(png_path, width = 18, height = 12, units='in', dpi=600) 
ggplot(data=REACTOME) + geom_bar(aes(x=reorder(term_name,intersection_size),y=intersection_size, fill=-log10(adjusted_p_value)), 
                             stat='identity') + 
  coord_flip() + scale_fill_gradient(expression(-log["10"](adjusted_p_value)),low="blue", high = "red") + 
  xlab("") + ylab("Gene count") + scale_y_continuous(expand=c(0, 0)) + 
  theme(axis.text.x=element_text(color="black",size=rel(1.5)), 
        axis.text.y=element_text(color="black", size=rel(1.6)), 
        axis.title.x = element_text(color="black", size=rel(1.6)), 
        legend.text=element_text(color="black",size=rel(1.0)), 
        legend.title = element_text(color="black",size=rel(1.1)))
# legend.position=c(0,1),legend.justification=c(-1,0) # legend.position="top", ) 
dev.off()

#气泡图：
rm(list=ls()) 
library(Cairo)
library(ggplot2) 
REACTOME <- read.csv('reactome_file.csv')
REACTOME <- REACTOME[1:20,]
REACTOME$rich_factor <- as.numeric(REACTOME$intersection_size/REACTOME$term_size)
png_path="./bubble_plot.png" 
CairoPNG(png_path, width = 10, height = 6, units='in', dpi=600) 
ggplot(REACTOME,aes(x=rich_factor,y=term_name))+
  geom_point(aes(size=intersection_size,color=-1*log10(adjusted_p_value)))+ 
  scale_colour_gradient(low="blue",high="red")+ 
  labs(color=expression(-log[10](P.value)), size="Gene number", 
       x="Rich Ratio", y="Pathway name", title="Term enrichment")+ 
  theme_bw()+ 
  theme( axis.text.y = element_text(size = rel(1.3)), 
         axis.title.x = element_text(size=rel(1.3)), axis.title.y = element_blank()) 
dev.off()
                                                                                                                                                                                                                        
#画基因-通路方格图（就是显示每条通路里面富集了哪些基因）
#这里的interactions是通路里面的基因，term_name是通路。
#把通路里面的基因转化为character类型
class(REACTOME$intersections)
REACTOME$intersections <- as.character(REACTOME$intersections)
#把所有的基因集合到一起，看看有多少个
DEG <- c()
i=1
for (i in 1:20) {
  tmp <- unlist(strsplit(REACTOME$intersections[i],split = ','))
  DEG <- c(DEG,tmp)
  i=i+1
}
DEG
length(unique(DEG))#基因的数目
nrow(REACTOME)#通路的数目
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







