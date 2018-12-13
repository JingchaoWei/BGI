#标准化（均一化）的矩阵不能作为DESeq2的输入数据

library(DESeq2)

data <- read.csv('all.csv')
colnames(data)
rownames(data) <- data$Other.Gene.ID
data <- data[,2:11]
data <- round(data)#数据readcount是RSEM定量的结果，并不是常见的htseq-count的结果，
#所以count值会有小数点，而DESeq2包不支持count数有小数点，所以这里需要round取整


#设置分组信息以及构建dds对象
condition<-factor(c(rep("FK",5),rep("HN",5)))
coldata<-data.frame(condition,row.names = colnames(data))
dds<-DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)

#使用DESeq函数进行估计离散度，然后进行标准的差异表达分析，得到res对象结果
dds <- DESeq(dds)
res <- results(dds)


#最后设定阈值，筛选差异基因，导出数据
table(res$padj <0.05)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),
                 by="row.names",sort=FALSE)
write.table(resdata,file = "FK_vs_HN_DESeq2.csv",quote = F,sep = "\t",
            col.names = T,row.names = F)



resdata <- resdata[abs(resdata$log2FoldChange)>=2 & resdata$padj<=0.05,]
resdata <- subset(x = resdata,subset=abs(resdata$log2FoldChange)>=2 & resdata$padj<=0.05)


BGI_DEG <- read.csv('DEG.csv')
BGI_DEG <- BGI_DEG[abs(BGI_DEG$log2.HN.FK.)>=2 & BGI_DEG$Qvalue.FK.vs.HN.<=0.05,]
head(BGI_DEG)
tmp <- merge(resdata,test,by.x = 'Row.names', by.y = 'Other.Gene.ID',all=F)


my_res <- as.character(resdata$Row.names)
BGI_res <- as.character(BGI_DEG$Other.Gene.ID)
setdiff(my_res,BGI_res)
setdiff(BGI_res,my_res)

BGI_DEG <- BGI_DEG[order(BGI_DEG$log2.HN.FK.),]
resdata <- resdata[order(resdata$log2FoldChange),]
head(BGI_DEG)
head(resdata)

