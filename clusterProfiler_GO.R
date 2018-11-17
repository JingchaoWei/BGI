# Fri Nov 16 15:44:08 2018 ------------------------------
load('Data.Rdata')

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#convert gene symbols to ENTREZID ID
#Attention: BGI report includes entrezID. No need to run the following 2 commands. Just filter out
#the BGI_novel genes.

#keytypes(org.Hs.eg.db)
#entr_up <- bitr(upgenes$Other.Gene.ID,fromType = 'SYMBOL',toType = c('ENTREZID','GENENAME'),OrgDb = 'org.Hs.eg.db')

index <- grep('BGI',upgenes$ï..Gene.ID)#Filter our BGI novel genes.
upgenes <- upgenes[-index,] 


#The ID type (both fromType & toType) should be one of ‘kegg’, ‘ncbi-geneid’, ‘ncbi-proteinid’ or 
#‘uniprot’. The ‘kegg’ is the primary ID used in KEGG database. The data source of KEGG was from NCBI. 
#A rule of thumb for the ‘kegg’ ID is entrez ID for eukaryote species and Locus ID for prokaryotes.

#keg2np <- bitr_kegg(entr_up$ENTREZID,fromType = 'kegg',toType = 'ncbi-proteinid',organism = 'hsa')
#head(keg2np)

#GO analysis
#GO classification/GoupGO 分析
entr_up <- as.vector(upgenes$ï..Gene.ID)

ggo <- groupGO(gene = entr_up,OrgDb = org.Hs.eg.db,ont = 'BP',level = 3,keyType = 'ENTREZID',readable = T)
head(ggo)
barplot(ggo,drop=T,x = "GeneRatio",showCategory=20,title='GO_GeneRatio')

#GO over-representation test/enrichGO 分析
index <- grep('BGI',background_all$ï..Gene.ID)
background_all <- background_all[-index,]#filter out BGI_novel genes.

entr_bg <- as.vector(background_all$ï..Gene.ID)
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

ego2 <- gofilter(ego,level = 3)#remove specific GO terms or GO level (in case of redundance)
barplot(ego2,drop=T,showCategory = 20,title=paste0("enrichGo_",'BP'))
dotplot(ego2,showCategory = 20,title=paste0("enrichGo_",'BP'))
emapplot(ego2, showCategory=20,color = "p.adjust",layout = 'kk',title=paste0("enrichGo_",'BP'))
cnetplot(ego2, showCategory = 20,title=paste0("enrichGo_",'BP'))

#plot for BP/CC/MF individually
pdf('enrichGO_plot.PDF',width = 14,onefile = T)
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

#GO GSEA
#1st:
#prepare your own geneList
#https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList

geneList <- upgenes$log2.HN.FK.
names(geneList) <- as.character(upgenes$ï..Gene.ID)
geneList <- sort(geneList,decreasing = T)
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 1000,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
dotplot(ego3)






