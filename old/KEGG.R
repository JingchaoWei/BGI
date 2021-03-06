# Sat Nov 17 18:15:24 2018 ------------------------------
rm(list=ls())
load('Data.Rdata')

index <- grep('BGI',upgenes$GeneID)#Filter our BGI novel genes.
upgenes <- upgenes[-index,] 

index <- grep('BGI',downgenes$GeneID)#Filter our BGI novel genes.
downgenes <- downgenes[-index,] 

index <- grep('BGI',background_all$GeneID)
background_all <- background_all[-index,]#filter out BGI_novel genes.

index <- grep('BGI',background_DEG$GeneID)
background_DEG <- background_DEG[-index,]

entr_up <- as.vector(upgenes$GeneID)#这些上调的基因可以导入clueGO直接进行分析了。
entr_down <- as.vector(downgenes$GeneID)
entr_bg <- as.vector(background_all$GeneID)
entr_up_down <- c(entr_up,entr_down)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pathview)
library(stringr)
library(dplyr)
library(DOSE)
library(ggplot2)
library(magrittr)

#kegg over-representation test
kk <- enrichKEGG(gene         = entr_down,universe = entr_bg,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,minGSSize = 10, maxGSSize = 2000,
                 use_internal_data = F)
head(kk)
kk <- setReadable(kk,keytype = 'ENTREZID',OrgDb = org.Hs.eg.db)#转化ID
#plot
pdf('Kegg_plot_down.pdf',width = 14,onefile = T)
barplot(kk,drop=T,showCategory = 20,x = "GeneRatio")
dotplot(kk,showCategory = 20)
emapplot(kk, showCategory=20,color = "p.adjust",layout = 'kk')
cnetplot(kk, showCategory = 20)
dev.off()
#or plot using ggplot2
tmp <- kk@result[kk@result$qvalue<0.01,]#提取数据框
#将generatio由字符转化为数值
str(tmp)
Char2Num <- function(i) {
  tmp1 <- tmp[i,3]%>%str_split(pattern = '/')%>%unlist()%>%as.numeric()
  tmp2 <- divide_by(tmp1[1],tmp1[2])
  return(tmp2)
}
tmp3 <- vector()
for (i in 1:nrow(tmp)) {
  tmp4 <- Char2Num(i)
  tmp3 <- c(tmp3,tmp4)
}
tmp3
tmp$GeneRatio <- tmp3
str(tmp)
#plot
library(ggplot2)
ggplot(tmp,aes(x=GeneRatio,y=Description))+
  geom_point(aes(size=Count,color=qvalue))+
  scale_colour_gradient(low="blue",high="red")+ 
  scale_x_continuous()+
  labs( color=expression(qvalue), x="Gene Ratio",
        y="Pathway name", title="Pathway enrichment")+ 
  theme_bw()+
  theme(axis.text.y = element_text(size = rel(1)), 
                    axis.title.x = element_text(size=rel(1)), 
                    axis.title.y = element_blank() ) 



#browseKEGG
name_kk <- kk@result$ID
for (i in 1:20) {
  browseKEGG(kk, name_kk[i])
}
save(kk,name_kk,file = "kegg_result.Rdata")






#get genes in significant pathways or specific pathways
load('kegg_result.Rdata')
result <- kk@result
nrow(result)
sig_genes_tmp <- result[result$p.adjust<0.05,]
#or specify pathways that you are interested in:
result$Description
sig_genes_tmp <- result[c(4,8),]
result[c(2,3,4,8,20),]
MySplit <- function(i){
  tmp <- sig_genes_tmp[i,'geneID'] %>% str_split(pattern = '/') %>% unlist(recursive = T)
  return(tmp)
}
sig_genes <- vector()
for (i in 1:nrow(sig_genes_tmp)) {
  tmp <- MySplit(i)
  sig_genes <- c(sig_genes,tmp)
}
sig_genes <- unique(sig_genes)
length(sig_genes)
cat(sig_genes,sep = '\n')#按行打印出来，可以直接复制到stringDB里面做蛋白相互作用
write.table(sig_genes,'KEGG_sig_genes.txt',sep = '\n',quote = F,
            col.names = F,row.names = F)



#比较基因集的生物学功能，比如比较上下调基因
geneset <- list(entr_up,entr_down)
names(geneset) <- c('Up','Down')#一定要给list命名，不然后面会报错
compare <- compareCluster(geneClusters = geneset,fun = 'enrichKEGG')
dotplot(compare)
#clusterProfiler出的图，我们可以用ggplot2随意改，比如改颜色
dotplot(compare)+ scale_color_continuous(low='purple', high='green')


#KEGG Gene Set Enrichment Analysis
geneList <- background_all$log2.HN.FK.
names(geneList) <- as.character(background_all$GeneID)
geneList <- sort(geneList,decreasing = T)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 10,maxGSSize = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE,
               use_internal_data = F)
kk2 <- setReadable(kk2,keytype = 'ENTREZID',OrgDb = org.Hs.eg.db)

dotplot(kk2,showCategory=20)

#GSEA plot
pdf(file = 'Kegg_GSEA.pdf',width = 10,onefile = T)
dotplot(kk2,showCategory=20)
for (i in 1:20) {
  tmp <- gseaplot(kk2, geneSetID = kk2@result$ID[i],
                  title=paste(kk2@result$ID[i],kk2@result$Description[i],sep = '_'))
  print(tmp)
}
dev.off()

#use pathview to see gene changes in pathway
name_kk2 <- kk2@result$ID
for (i in 1:20) {
  hsa04110 <- pathview(gene.data  = geneList,
                       pathway.id = name_kk2[i],
                       species    = "hsa",
                       limit      = list(gene=max(abs(geneList)), cpd=1))
}





#KEGG Module over-representation test
mkk <- enrichMKEGG(gene = entr_up,
                   organism = 'hsa',pvalueCutoff = 0.05,
                   universe = entr_bg,minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.05)
dotplot(mkk,showCategory=20)

#KEGG Module Gene Set Enrichment Analysis(GSEA)
geneList <- background_all$log2.HN.FK.
names(geneList) <- as.character(background_all$ï..Gene.ID)
geneList <- sort(geneList,decreasing = T)
mkk2 <- gseMKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 nPerm        = 1000,
                 minGSSize    = 10,maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
gseaplot(ego3, geneSetID = "hsa04145")
