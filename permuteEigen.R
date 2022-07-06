#####################################
# Cell-wise permutation to determine number of significant eigenvalues and therefore significant PCs to look at
# Run on cluster, took just over a minute per shuffling for this dataset
#####################################
#this scale.data file needs to be in the directory directory 
#the following code will be copied made into an R file using 'nano'. e.g. 
# <- )
#nano permuteEigen.R
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


print(max(randEV)) #this number will be input for next 
