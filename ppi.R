data <- read.csv('up.csv')
write.table(data$ï..Gene.ID,'PPI_input_down_HN.csv',
            quote = F,sep = '\t',col.names = F,row.names = F)
