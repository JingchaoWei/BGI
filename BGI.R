# Fri Nov 16 15:16:57 2018 ------------------------------


background_DEG <- read.csv('DEG.csv')
background_all <- read.csv('all.csv')
attach(background_DEG)
upgenes <- background_DEG[FK.FPKM>=1&HN.FPKM>=1&log2.HN.FK.>=2&Qvalue.FK.vs.HN.<0.01,]
downgenes <- background_DEG[FK.FPKM>=1&HN.FPKM>=1&log2.HN.FK.<=-2&Qvalue.FK.vs.HN.<0.01,]
detach(background_DEG)
save(background,upgenes,downgenes,background_all,file = 'Data.Rdata')


