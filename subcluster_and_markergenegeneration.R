#source('https://bioconductor.org/biocLite.R')
#biocLite('sva')
#install.packages(c("igraph","ggplot2","grid","Matrix","gmodels","RANN","purrr",'factoextra'), repos='http://archive.linux.duke.edu/cran/')
library(sva)
#library(factoextra)
library(igraph)
library(ggplot2)
library(Matrix)
library(gmodels)
library(RANN)
library(purrr)
source("class042518.R")


load('batch_200pca_3daysbc.RData')
eigens<-read.table('eigenvalues.txt',header=F)
eigens<-eigens[,2]
sigeig<-length(dsq.SNISham@pca.eigenvalues$ev[dsq.SNISham@pca.eigenvalues$ev>1+max(eigens)])

dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:sigeig,iterations=10:120,corMethod='pearson')
dsq.SISNI = perform_refined_clustering(dsq.SISNI, pcs.use = 1:sigeig, optimize.NN = FALSE, opt.NN=81,corMethod='pearson')#choose optimal NN parameter based on plots generated from previous line; highest avg silhouette score

dsq.SISNI<-merge.clusters.DE(dsq.SISNIm,min.de.genes=1,pcs.use=1:sigeig)
TPM.mat = exp(dsq.SISNI@data) - 1
Count.mat = dsq.SISNI@count.data[rownames(dsq.SISNI@data), colnames(dsq.SISNI@data)]

ident.louvain = dsq.SISNI@group

#at this point, marker genes were compared against a list of oligodendrocyte-enriched genes to account for myelin contamination.
#if any subcluster was enriched for myelin related genes that cluster was removed
#the exception to this being if the initial cluster was comprised of oligodendrocytes or OPCs (initial clusters 1, 3, and 20)
oligogenes=read.csv('Oligogenelist.csv',header=T)
oligogenes=oligogenes$x
badcells=c()
for (i in 1:length(table(ident.louvain))) {
markers= markers.binom(dsq.SISNI,clust.1=i)
assign(paste0('markers.',i),markers)
markplusoligo=c(rownames(markers),oligogenes)
myelincontam=markplusoligo[which(duplicated(markplusoligo))]
print(length(myelincontam))
if (length(myelincontam)>10){
  badcells=c(badcells,rownames(dsq.SISNI@meta.data[which(dsq.SISNI@meta.data$clust==i),]))
}
#
cleanmeta=dsq.SISNI@meta.data[-badcells,]

write.table(cleanmeta,'Submeta.c1.txt',quote=F,sep='\t',col.names=NA)
save.image('subclust1_thru_markers.RData')




  
