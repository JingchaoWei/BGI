#柱状图：
rm(list=ls()) 
library(ggplot2) 
library(Cairo) 
GO <- read.csv('GO_enrichment.csv')
colnames(GO)
png_path="./GO.png" 
CairoPNG(png_path, width = 18, height = 12, units='in', dpi=600) 
ggplot(data=GO) + geom_bar(aes(x=reorder(GO.Term,Term.Candidate.Gene.Num),y=Term.Candidate.Gene.Num, fill=-log10(Q.value)), 
                           stat='identity') + 
  coord_flip() + scale_fill_gradient(expression(-log["10"](Q.value)),low="blue", high = "red") + 
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
GO <- read.csv('kegg_HN_up.csv')
colnames(GO)
#png_path="./bubble_plot.png" 
#CairoPNG(png_path, width = 10, height = 6, units='in', dpi=600) 
pdf(file = 'bubble_plot.pdf',width = 8,height = 5)
ggplot(GO,aes(x=Rich.Ratio,y=Pathway.Name))+
  geom_point(aes(size=Term.Candidate.Gene.Num,color=-1*log10(Q.value)))+ 
  scale_colour_gradient(low="blue",high="red")+ 
  labs(color=expression(-log[10](Q.value)), size="Gene number", 
       x="Rich Ratio", y="Pathway name", title="Term enrichment (KEGG_up_HN)")+ 
  theme_bw()+ 
  theme( axis.text.y = element_text(size = rel(1.3)), 
         axis.title.x = element_text(size=rel(1.3)), axis.title.y = element_blank()) 
dev.off()

#画基因-通路方格图（就是显示每条通路里面富集了哪些基因）/mapping genes to pathways
#这里的interactions是通路里面的基因，term_name是通路。
#把通路里面的基因转化为character类型
go <- read.csv('GO2.csv')
go$GO.Term <- as.character(go$GO.Term)
go$GO.Term[1]
unlist(strsplit(go$GO.Term[1],split = c(';|//')))
#把所有的基因集合到一起，看看有多少个
nrow(go)#基因的数目
20#通路的数目
#创建矩阵，行列分别对应基因和通路
data <- matrix(0,35,20)
data <- as.data.frame(data)
colnames(data) <- GO$GO.Term
rownames(data) <- go$Other.Gene.ID
head(data)
#新创建的矩阵，其行列分别为基因和通路名称，去原始矩阵里面检索该通路里面是否有该基因
j=1
for (j in 1:20) {
  i=1
  for (i in 1:35) {
    data[i,j]
    index <- which(as.character(go$Other.Gene.ID) == rownames(data)[i])
    tmp_pathway <- unlist(strsplit(go$GO.Term[i],split = c(';|//')))
    data[i,j] <- ifelse(colnames(data)[j] %in% tmp_pathway,1,0)
  }
  j=j+1
}
View(data)
data <- t(data)
#可视化
library(pheatmap)
library(RColorBrewer)
pdf("heatmap.pdf",width = 10,height = 5,onefile = T)
par(mar=c(5.1, 4.1, 4.1, 2.1))
map <-pheatmap(data,col = c("white","#FF9999"),
               cluster_rows = F,cluster_cols=F,legend = F,
               border_color = 'grey',
               show_colnames=T,show_rownames=T,
               fontsize_col = 8, 
               fontsize_row = 8,fontsize = 8,
               main = 'Mapping Genes to Pathways')
dev.off()







