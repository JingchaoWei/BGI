data <- read.csv('DEG.csv',check.names = T)
data <- data[data$FK.FPKM>1&data$HN.FPKM>1&data$Qvalue.FK.vs.HN.<0.01,]
data <- data[data$log2.HN.FK.>2|data$log2.HN.FK.<(-2),]
index <- grep(pattern = 'BGI',x = data$Other.Gene.ID)
data <- data[-index,]
write.table(x = data,file = 'DEG_for_pathwayanalysis_kickout_BGINovelGenes.txt',quote = F,sep = '\t',row.names = F,col.names = T)
ID <- data$ï..Gene.ID
write.table(ID,'ID_for_pathway_analysis.txt',quote = F,sep = '\n',row.names = F,col.names = F)
