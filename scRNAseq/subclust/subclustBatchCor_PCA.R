#This script should be used for each initial cluster dataframe and proceeds through batch correction and PCA.
source('https://bioconductor.org/biocLite.R')
library(sva)
library(igraph)
library(ggplot2)
library(Matrix)
library(gmodels)
library(RANN)
source("class042518.R")

SNIShamdge = read.table("samples.c1.txt",header=T,sep="\t",row.names=1)
#for 0 padding
SNIShamdge[is.na(SNIShamdge)] <- 0
print(dim(SNIShamdge))
#remove all SNI04 cells since they're dogshit
#again <- SNIShamdge[, -grep("SNI04", colnames(SNIShamdge))]
#SNIShamdge <- again


mt.genes = grep("mt-", rownames(SNIShamdge), value = TRUE)
cells.use = colnames(SNIShamdge)[colSums(SNIShamdge[mt.genes, ])/colSums(SNIShamdge) < 0.1]
SNIShamdge = SNIShamdge[, cells.use]
dim(SNIShamdge)
dsq.SNISham=scDrop(count.data=SNIShamdge)
dsq.SNISham=initialize(dsq.SNISham, min.genes = 400, min.cells = ncol(SNIShamdge)/100, min.counts=1)
table(dsq.SNISham@meta$sample)
pdf('violin.pdf')
violinplot(dsq.SNISham,c("num.genes","num.trans"))
dev.off()

batchname = as.character(dsq.SNISham@meta$sample)
batchid = rep(1,length(batchname))
#batch cor by cohort
batchid[batchname=="Sham01"] = 1
batchid[batchname=="Sham02"] = 1
batchid[batchname=="Sham03"] = 1
batchid[batchname=="Sham04"] = 1
batchid[batchname=="Sham05"] = 1
batchid[batchname=="SNI01"] = 1
batchid[batchname=="SNI02"] = 1
batchid[batchname=="SNI03"] = 1
#batchid[batchname=="SNI04"] = 1

batchid[batchname=="Sham06"] = 2
batchid[batchname=="Sham07"] = 2
batchid[batchname=="Sham08"] = 2
batchid[batchname=="Sham09"] = 2
batchid[batchname=="Sham10"] = 2
batchid[batchname=="SNI05"] = 2
batchid[batchname=="SNI06"] = 2
batchid[batchname=="SNI07"] = 2
batchid[batchname=="SNI08"] = 2
batchid[batchname=="SNI09"] = 2
batchid[batchname=="SNI10"] = 2

dsq.SNISham=doBatchCorrection(dsq.SNISham, batch.cov=batchid)
batchname = as.character(dsq.SNISham@meta$sample)
batchid = rep(1,length(batchname))

#batch cor by dropseq day
#batch by dropseq day
batchid[batchname=="Sham01"] = 3
batchid[batchname=="Sham02"] = 4
batchid[batchname=="Sham03"] = 5
batchid[batchname=="Sham04"] = 5
batchid[batchname=="Sham05"] = 5
batchid[batchname=="SNI01"] = 3
batchid[batchname=="SNI02"] = 4
batchid[batchname=="SNI03"] = 5
#batchid[batchname=="SNI04"] = 5

batchid[batchname=="Sham06"] = 3
batchid[batchname=="Sham07"] = 4
batchid[batchname=="Sham08"] = 5
batchid[batchname=="Sham09"] = 5
batchid[batchname=="Sham10"] = 5
batchid[batchname=="SNI05"] = 3
batchid[batchname=="SNI06"] = 4
batchid[batchname=="SNI07"] = 5
batchid[batchname=="SNI08"] = 5
batchid[batchname=="SNI09"] = 5
batchid[batchname=="SNI10"] = 5
dsq.SNISham= secondBatchCorrection(dsq.SNISham, batch.cov = batchid)


dsq.SNISham=doPCA(dsq.SNISham,pcs.store=200)
data.plot = dsq.SNISham@pca.scores
data.plot$group = dsq.SNISham@group
write.table(dsq.SISNI@scale.data,'subclust1.dsq.SISNI.scale.data.txt',quote=F,sep='\t',col.names=NA)
save.image(file='batch_200pca_3daysbc.RData')

#for each subclustering iteration, it is necessary to determine the number of optimal PCs (using permuteEigen.R; permutations = 500)
