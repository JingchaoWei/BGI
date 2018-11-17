# Fri Nov 16 15:44:08 2018 ------------------------------
load('Data.Rdata')

library(clusterProfiler)
library(org.Hs.eg.db)
#convert gene symbols to ENTREZID ID
keytypes(org.Hs.eg.db)
entr_up <- bitr(upgenes$Other.Gene.ID,fromType = 'SYMBOL',toType = c('ENTREZID','GENENAME'),
            OrgDb = 'org.Hs.eg.db')
head(entr_up)

#The ID type (both fromType & toType) should be one of ‘kegg’, ‘ncbi-geneid’, ‘ncbi-proteinid’ or 
#‘uniprot’. The ‘kegg’ is the primary ID used in KEGG database. The data source of KEGG was from NCBI. 
#A rule of thumb for the ‘kegg’ ID is entrez ID for eukaryote species and Locus ID for prokaryotes.
keg2np <- bitr_kegg(entr_up$ENTREZID,fromType = 'kegg',toType = 'ncbi-proteinid',organism = 'hsa')
head(keg2np)

#GO analysis
#GO classification/GoupGO 分析
gene <- entr_up$ENTREZID
ggo <- groupGO(gene = gene,OrgDb = org.Hs.eg.db,ont = 'BP',level = 3,keyType = 'ENTREZID',readable = T)
head(ggo)
barplot(ggo,drop=T,showCategory=20)

#GO over-representation test/enrichGO 分析
entr_bg <- bitr(background_all$Other.Gene.ID,fromType = 'SYMBOL',toType = 'ENTREZID',
                OrgDb = 'org.Hs.eg.db')
head(entr_bg)

pdf('enrichGO_plot.PDF',width = 14,onefile = T)
for (i in c('CC','BP','MF')) {
  ego <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = i,
                  pAdjustMethod = 'BH',pvalueCutoff = 0.01,qvalueCutoff = 0.05,
                  readable = T,universe = entr_bg$ENTREZID)
  head(as.data.frame(ego))
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

