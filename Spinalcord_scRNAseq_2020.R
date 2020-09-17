#To run an interactive session, e.g. R or matlab, with a time limit of your choosing (here is 1 day and 8 hours):
#for this do 100g, several days
#make sure display-variable is set in environment
ssh -Y kylius0@longleaf.unc.edu

srun -n 1 --mem=100g --x11=first --time=8-12 --pty R




install.packages(c("igraph","ggplot2","grid","Matrix","gmodels","RANN","extrafont"))
#needs to download from mirror
58
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
library(igraph)
library(ggplot2)
library(Matrix)
library(gmodels)
library(RANN)
library(extrafont)
source("class042518.R")
font_import()
y
loadfonts()
SISNIdge = read.table("SISNIdge.txt",header=T,sep="\t",row.names=1)
#for 0 padding
SISNIdge[is.na(SISNIdge)] <- 0
print(dim(SISNIdge))


mt.genes = grep("mt-", rownames(SISNIdge), value = TRUE)
cells.use = colnames(SISNIdge)[colSums(SISNIdge[mt.genes, ])/colSums(SISNIdge) < 0.1]
SISNIdge = SISNIdge[, cells.use]
dim(SISNIdge)
dsq.SISNI=scDrop(count.data=SISNIdge)
dsq.SISNI=initialize(dsq.SISNI, min.genes = 400, min.cells = ncol(SISNIdge)/100, min.counts=3)
table(dsq.SISNI@meta$sample)
violinplot(dsq.SISNI,c("num.genes","num.trans"))

batchname = as.character(dsq.SISNI@meta$sample)
batchid = rep(1,length(batchname))
#batch cor by cohort
batchid[batchname=="SI01"] = 1
batchid[batchname=="SI02"] = 1
batchid[batchname=="SI03"] = 1
batchid[batchname=="SI04"] = 1
batchid[batchname=="SI05"] = 1
batchid[batchname=="SNI01"] = 1
batchid[batchname=="SNI02"] = 1
batchid[batchname=="SNI03"] = 1
#batchid[batchname=="SNI04"] = 1

batchid[batchname=="SI06"] = 2
batchid[batchname=="SI07"] = 2
batchid[batchname=="SI08"] = 2
batchid[batchname=="SI09"] = 2
batchid[batchname=="SI10"] = 2
batchid[batchname=="SNI05"] = 2
batchid[batchname=="SNI06"] = 2
batchid[batchname=="SNI07"] = 2
batchid[batchname=="SNI08"] = 2
batchid[batchname=="SNI09"] = 2
batchid[batchname=="SNI10"] = 2

dsq.SISNI=doBatchCorrection(dsq.SISNI, batch.cov=batchid)
batchname = as.character(dsq.SISNI@meta$sample)
batchid = rep(1,length(batchname))

#batch cor by dropseq day
####do 3 batches based on dropseq day rather than 4
#e.g. do 3s, 4s, and lump 5s and 6s together into 5s
batchid[batchname=="SI01"] = 3
batchid[batchname=="SI02"] = 4
batchid[batchname=="SI03"] = 5
batchid[batchname=="SI04"] = 5
batchid[batchname=="SI05"] = 5
batchid[batchname=="SNI01"] = 3
batchid[batchname=="SNI02"] = 4
batchid[batchname=="SNI03"] = 5
#batchid[batchname=="SNI04"] = 5

batchid[batchname=="SI06"] = 3
batchid[batchname=="SI07"] = 4
batchid[batchname=="SI08"] = 5
batchid[batchname=="SI09"] = 5
batchid[batchname=="SI10"] = 5
batchid[batchname=="SNI05"] = 3
batchid[batchname=="SNI06"] = 4
batchid[batchname=="SNI07"] = 5
batchid[batchname=="SNI08"] = 5
batchid[batchname=="SNI09"] = 5
batchid[batchname=="SNI10"] = 5
#do second batch correction on batch.data
#batch corrected data will be batch.data
dsq.SISNI= secondBatchCorrection(dsq.SISNI, batch.cov = batchid)


dsq.SISNI=doPCA(dsq.SISNI,pcs.store=120)
#good time to save workspace

write.table(dsq.SISNI@scale.data,'dsq.SISNI.scale.data.txt',quote=F,sep='\t',col.names=NA)
#####################################
# Cell-wise permutation to determine number of significant eigenvalues and therefore significant PCs to look at
# Run on cluster, took just over a minute per shuffling for this dataset
#####################################
#this scale.data file needs to be in the directory directory 
#the following code will be copied made into an R file using 'nano'. e.g. 
# <- )
#nano permuteEignen.R
#paste the code
#cntrl +x to save and exit
# note the following was slightly changed based on what the file was named as (in this case "dsq.SISNI.scale.data.txt")
dsq.SISNI.scale.data <- read.table(dsq.SISNI.scale.data.txt, header=T, sep="\t",row.names=1)
library(gmodels)
nperm=1000 #can be adjusted or run in smaller permutation blocks
randEV=rep(NA,nperm)

for (i in 1:nperm) {
  perm.mat=matrix(NA,nrow=length(rownames(dsq.SISNI.scale.data)),ncol=length(colnames(dsq.SISNI.scale.data)))
  
  for (j in (1:length(rownames(dsq.SISNI.scale.data)))) {
    s = sample(1:length(dsq.SISNI.scale.data[j,]))
    perm.mat[j,]=rbind(dsq.SISNI.scale.data[j,s])	
  }
  
  rownames(perm.mat) = rownames(dsq.SISNI.scale.data)
  colnames(perm.mat) = c(1:length(colnames(dsq.SISNI.scale.data)))
  
  data.use=perm.mat
  pc.genes = rownames(perm.mat)
  
  #Remove genes with zero variation
  
  pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
  genes.use = pc.genes[pc.genes.var>0]
  pc.data = data.use[genes.use,]
  
  pca.obj = fast.prcomp(t(pc.data),center=FALSE, scale=FALSE)
  ev=pca.obj$sdev^2	
  randEV[i] = max(ev)
  print(max(ev))
}
print(randEV)


print(max(randEV))

maxpca=length(dsq.SISNI@pca.eigenvalues$ev[dsq.SISNI@pca.eigenvalues$ev>max(randEV)])


#tsne in python

write.table(dsq.SISNI@pca.scores[,1:maxpca],"SISNIcsc_top11pcs.txt",quote=F,sep="\t",col.names=NA)

####tsne visualization parameters tested on significant PCs in python using different perp and learning rate parameters
####.txt file with tsne coordinates written with tsne_write.py

tsne.y = read.table("Xperp5100lr.txt",sep="\t",header=T,row.names=1)#

tSNE1=tsne.y$X0
tSNE2=tsne.y$X1
data.plot = dsq.SISNI@tsne.y
data.plot$group = dsq.SISNI@group
ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()+geom_point(size=1, aes(color=factor(group))) + theme(axis.text=element_text(size=13), axis.title=element_text(size=14,face="bold"))+theme(legend.text=element_text(size=16)) +theme(legend.title = element_text(size=18)) + labs(color="Condition") + ggtitle("SC cells by Cluster") +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(plot.title=element_text(size=33))



ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
#cluster; note; try different number nearest neighbros
dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:maxpca,iterations=10:120,corMethod='pearson')
dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:maxpca, optimize.NN = FALSE, opt.NN=81,corMethod='pearson')#choose optimal NN parameter
length(table(dsq.SISNI@group))



TPM.mat = exp(dsq.SISNI@data) - 1
Count.mat = dsq.SISNI@count.data[rownames(dsq.SISNI@data), colnames(dsq.SISNI@data)]

ident.louvain = dsq.SISNI@group

#make condition column in meta file
dsq.SISNI@meta$condition <- ifelse(grepl("^SI",row.names(dsq.SISNI@meta))==TRUE,"SI","SNI")


#make marker gene dfs
for (i in 1:length(table(ident.louvain))) {
        markers = markers.binom(dsq.SISNI, clust.1 = i, effect.size=log(2), TPM.mat = TPM.mat, Count.mat = Count.mat)
        markers	= markers[order(markers$pval),]
        assign(paste0('markers.',i),markers)
        siggenes<-rownames(markers[markers$log.effect>1,])
}

################################################################
################################################################
#End of initial clustering. Next steps were to perform a second round of clustering on each initial cluster in order to identify populations with better acuity.
#Aim is to fill dsq.SISNI@subclust and dsq.SISNI@subgroup objects with new IDs.
#dataframe .txt files of each cluster (e.g. dsq.SISNI@count.data[,which(dsq.SISNI@group==1)]) were generated to run through the above clustering pipeline again.
#New running subgroup/subclust IDs were generated from each initial cluster. For example, initial cluster 1 became four clusters and were named 1, 2, 3, 4. Initial cluster 2 became 
################################################################


