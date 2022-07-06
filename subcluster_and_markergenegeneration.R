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

dsq.SNISham<-merge.clusters.DE(dsq.SNISham,min.de.genes=1,pcs.use=1:sigeig)
TPM.mat = exp(dsq.SNISham@data) - 1
Count.mat = dsq.SNISham@count.data[rownames(dsq.SNISham@data), colnames(dsq.SNISham@data)]

ident.louvain = dsq.SNISham@group
for (i in 1:length(table(ident.louvain))) {
markers= markers.binom(dsq.SNISham,clust.1=i)
assign(paste0('markers.',i),markers)
}
save.image('subclust1_thru_markers.RData')

