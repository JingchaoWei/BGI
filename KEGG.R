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

#kegg over-representation test
kk <- enrichKEGG(gene         = entr_up_down,
                 organism     = 'hsa',
                 pvalueCutoff = 0.01,minGSSize = 10, maxGSSize = 500,
                 use_internal_data = F)
head(kk)
kk <- setReadable(kk,keytype = 'ENTREZID',OrgDb = org.Hs.eg.db)#转化ID
#plot
pdf('Kegg_plot.pdf',width = 14,onefile = T)
barplot(kk,drop=T,showCategory = 20,x = "GeneRatio")
dotplot(kk,showCategory = 20)
emapplot(kk, showCategory=20,color = "p.adjust",layout = 'kk')
cnetplot(kk, showCategory = 20)
dev.off()
#browseKEGG
name_kk <- kk@result$ID
browseKEGG(kk, name_kk[2])
#get genes in significant pathways
sig_genes_tmp <- kk@result
nrow(sig_genes_tmp)
sig_genes_tmp <- sig_genes_tmp[sig_genes_tmp$p.adjust<0.05,]
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
