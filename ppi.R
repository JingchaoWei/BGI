upgene_in_FK <- read.csv('up_genes_in_FK.csv')
write.table(upgene_in_FK$Gene.ID,'up_genes_in_FK.csv',
            quote = F,sep = '\t',col.names = F,row.names = F)
