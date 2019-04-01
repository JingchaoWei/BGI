# Fri Dec 07 19:52:55 2018 ------------------------------
rm(list=ls())
load('Data.Rdata')
head(upgenes)

#制作GSEA分析文件

##1.制作表达文件gct格式
GSEA_exp_file <- background_all
colnames(GSEA_exp_file)
GSEA_exp_file <- GSEA_exp_file[,c(10,2:9)]

colnames(GSEA_exp_file)[1] <- 'Symbol'

index <- grep('BGI',GSEA_exp_file$Symbol)#Filter our BGI novel genes.
GSEA_exp_file <- GSEA_exp_file[-index,] 


#按照GSEA格式要求，给GSEA需要的矩阵添加第二列description,虽然会被GSEA忽略，
#但是需要添加这一列。先判断探针注释文件行数和gsea矩阵行数是否相同：
library(tibble)
GSEA_exp_file <- add_column(GSEA_exp_file,DESCRIPTION=rep('NA',nrow(GSEA_exp_file)),
                            .before = 2)


#手动添加前两行，如下：
#    #1.2     前面的“#”不能省略,这里总是#1.2
#    54675	50    分别为gene number，样本数. Tab-delimited
dim(GSEA_exp_file)#看看有多少个样本（列数-2），多少个基因（行数）
write.table(GSEA_exp_file,'GSEA_expression.gct',row.names = F,col.names = T,
            quote = F,sep = '\t')



##2.制作分组文件，cls格式
#格式如下所示：第一行50代表样本数目，2代表分2组，空格间隔，1照抄；
#第二行井号注释说明分组信息，空格间隔；第三行为每个样本对应的组名，空格分隔
#  50 2 1   
#  # tumor benign
#  1	1	1	1	1	1	1	1	1	1	1	1	1	1	0	0	1	1	1	1	1	1	1	1	1	1	1	0	1	1	0	1	1	1	1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0
colnames(GSEA_exp_file)

#or: 手动添加前两行，总体格式如下(use space or tab, both are fine)：
#Line format:(number of samples) (space) (number of classes) (space) 1
#Line format:# (space) (class 0 name) (space) (class 1 name)
#Line format:(sample 1 class) (space) (sample 2 class) (space) ... (sample N class)
#然后就可以用GSEA软件进行分析了

