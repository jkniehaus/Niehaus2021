#####################################
# Cell-wise permutation to determine number of significant eigenvalues and therefore significant PCs to look at
#Eigenvalues for each permutation will be stored in 'eigenvalues.txt' file
#####################################
#this scale.data file needs to be in the directory directory 

library(gmodels)

data=read.table("dsq.SNISham.scale.data.batchcor_noSNI4.txt",header=T,sep="\t",row.names=1)
data=as.matrix(data)

nperm=2
randEV=rep(NA,nperm)

for (i in 1:nperm) {
        perm.mat=matrix(NA,nrow=length(rownames(data)),ncol=length(colnames(data)))

        for (j in (1:length(rownames(data)))) {
                s = sample(1:length(data[j,]))
                perm.mat[j,]=rbind(data[j,s])
        }

        rownames(perm.mat) = rownames(data)
        colnames(perm.mat) = c(1:length(colnames(data)))

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
        sink('eigenvalues.txt',append=T)
        print(max(ev))
        sink()
}
print(randEV)

print(max(randEV)

load('thruPreprocessing.RData')



print(max(randEV)) #this number will be input for next 
maxpca=length(dsq.SISNI@pca.eigenvalues$ev[dsq.SISNI@pca.eigenvalues$ev>max(randEV)]) 
#for tsne in python, write table of first X PCs
                       
write.table(dsq.SISNI@pca.scores[,1:maxpca],"SISNIcsc_toppcs.txt",quote=F,sep="\t",col.names=NA)
