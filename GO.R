# Fri Nov 16 15:44:08 2018 ------------------------------
rm(list=ls())
load('Data.Rdata')

index <- grep('BGI',upgenes$ï..Gene.ID)#Filter our BGI novel genes.
upgenes <- upgenes[-index,] 

index <- grep('BGI',background_all$ï..Gene.ID)
background_all <- background_all[-index,]#filter out BGI_novel genes.

index <- grep('BGI',background_DEG$ï..Gene.ID)
background_DEG <- background_DEG[-index,]

entr_up <- as.vector(upgenes$ï..Gene.ID)#这些上调的基因可以导入clueGO直接进行分析了。
entr_bg <- as.vector(background_all$ï..Gene.ID)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#convert gene symbols to ENTREZID ID
#Attention: BGI report includes entrezID. No need to run the following 2 commands. Just filter out
#the BGI_novel genes.

#keytypes(org.Hs.eg.db)
#entr_up <- bitr(upgenes$Other.Gene.ID,fromType = 'SYMBOL',toType = c('ENTREZID','GENENAME'),OrgDb = 'org.Hs.eg.db')


#The ID type (both fromType & toType) should be one of ‘kegg’, ‘ncbi-geneid’, ‘ncbi-proteinid’ or 
#‘uniprot’. The ‘kegg’ is the primary ID used in KEGG database. The data source of KEGG was from NCBI. 
#A rule of thumb for the ‘kegg’ ID is entrez ID for eukaryote species and Locus ID for prokaryotes.

#keg2np <- bitr_kegg(entr_up$ENTREZID,fromType = 'kegg',toType = 'ncbi-proteinid',organism = 'hsa')
#head(keg2np)

#GO analysis
#GO classification/GoupGO 分析

ggo <- groupGO(gene = entr_up,OrgDb = org.Hs.eg.db,ont = 'BP',level = 3,keyType = 'ENTREZID',readable = T)
head(ggo)
barplot(ggo,drop=T,x = "GeneRatio",showCategory=20,title='GO_GeneRatio')

#GO over-representation test/enrichGO 分析
head(entr_bg)

ego <- enrichGO(gene = entr_up,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = 'BP',
                pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                readable = T,universe = entr_bg)#background genes: use all gene. Don't use DEGs only.
head(ego)

#plot
barplot(ego,drop=T,showCategory = 20,title=paste0("enrichGo_",'BP'))
dotplot(ego,showCategory = 20,title=paste0("enrichGo_",'BP'))
emapplot(ego, showCategory=20,color = "p.adjust",layout = 'kk',title=paste0("enrichGo_",'BP'))
cnetplot(ego, showCategory = 20,title=paste0("enrichGo_",'BP'))
goplot(ego,showCategory =20)

ego2 <- gofilter(ego,level = 3)#remove specific GO terms or GO level (in case of redundance)
barplot(ego2,drop=T,showCategory = 20,title=paste0("enrichGo_",'BP'))
dotplot(ego2,showCategory = 20,title=paste0("enrichGo_",'BP'))
emapplot(ego2, showCategory=20,color = "p.adjust",layout = 'kk',title=paste0("enrichGo_",'BP'))
cnetplot(ego2, showCategory = 20,title=paste0("enrichGo_",'BP'))
#get genes in significant pathways
sig_genes_tmp <- ego@result
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




#plot for BP/CC/MF individually
pdf('GO_plot.PDF',width = 14,onefile = T)
for (i in c('CC','BP','MF')) {
  print(paste0('Now processing:',i))
  ego <- enrichGO(gene = entr_up,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = i,
                  pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                  readable = T,universe = entr_bg)
  a <- barplot(ego,drop=T,showCategory = 20,title=paste0("enrichGo_",i))
  b <- dotplot(ego,showCategory = 20,title=paste0("enrichGo_",i))
  c <- emapplot(ego, showCategory=20,color = "p.adjust",layout = 'kk',title=paste0("enrichGo_",i))
  d <- cnetplot(ego, showCategory = 20,title=paste0("enrichGo_",i))
  print(a)
  print(b)
  print(c)
  print(d)
}
dev.off()

#GO GSEA: 
#1st:
#prepare your own geneList
#https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList

geneList <- background_all$log2.HN.FK.
names(geneList) <- as.character(background_all$ï..Gene.ID)
geneList <- sort(geneList,decreasing = T)
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

dotplot(ego3,showCategory = 20)

ego3 <- setReadable(ego3,keytype = 'ENTREZID',OrgDb = org.Hs.eg.db)

#GSEA plot
pdf(file = 'GO_GSEA.pdf',width = 10,onefile = T)
dotplot(ego3,showCategory = 20)
for (i in 1:20) {
  tmp <- gseaplot(ego3, geneSetID = ego3@result$ID[i],
                  title=paste(ego3@result$ID[i],ego3@result$Description[i],sep = '_'))
  print(tmp)
}
dev.off()



