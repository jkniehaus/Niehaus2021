load('thruPreprocessing.RData')
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
#randEV will be the maximum eigenvalue printed from the permuteEigen.R code rounded to the next highest integer

#read tSNE coordinates from tsne_write.py script
tsne.y = read.table("Xperp5100lr.txt",sep="\t",header=T,row.names=1)#

tSNE1=tsne.y$X0
tSNE2=tsne.y$X1
data.plot = dsq.SISNI@tsne.y
data.plot$group = dsq.SISNI@group
ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()+geom_point(size=1, aes(color=factor(group))) + theme(axis.text=element_text(size=13), axis.title=element_text(size=14,face="bold"))+theme(legend.text=element_text(size=16)) +theme(legend.title = element_text(size=18)) + labs(color="Condition") + ggtitle("SC cells by Cluster") +
  theme(plot.title=element_text(hjust=0.5)) + 
  theme(plot.title=element_text(size=33))



ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()
#cluster; note; try different number nearest neighbors; 10-120
dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:maxpca,iterations=10:120,corMethod='pearson') #computes silhouette scores for each NN 
dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:maxpca, optimize.NN = FALSE, opt.NN=81,corMethod='pearson')#Clusters by optimal NN, 81NN produced highest average silhouette score per cluster
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

#at this point we can explore marker genes for each of the initial 20 clusters; heterogeneity within each cluster suggests subclusters exist. 
#the preprocessing, batch correcting, and iterative clustering was then redone on each of the initial 20 clusters to identify subcluster populations
#to do this, we break out each initial cluster into a genes x cells dataframe
for (i in 1:length(table(ident.louvain))) {
    subcount=dsq.SISNI@count.data[,which(dsq.SISNI@meta.data$clust==i)]
    write.table(subcount,paste0('samples.c',i,'.txt'),quote=F,sep='\t',col.names=NA)
}

#For each samples.c#.txt file, run subclust_batch_and_pca.R
