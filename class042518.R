# S4 class file for Shekhar et al., "Comprehensive classification of retinal bipolar cells using single-cell transcriptomics", Cell, 2016

# Required packages
library(genefilter)
library(sva,lib.loc="/proj/jmsimon/Rlibs")
library(igraph)
library(ggplot2)
library(Matrix)
library(gmodels)
library(RANN)
library(reshape)
library(cluster)


# Define slots
scDrop <- setClass("scDrop", slots = 
                     c(count.data="data.frame", data="data.frame",batch.data="data.frame", scale.data="matrix", batch.data.fin="data.frame",
                       group="vector", pca.load="data.frame",pca.scores="data.frame",pca.eigenvalues="data.frame",hdbgroup='vector',
                       meta="data.frame",tsne.y="data.frame", cell.names="vector",sex="vector",sils="vector",numclust="vector",subgroup="vector",submarker.scale="data.frame",
		       sub.avgs="data.frame",submarkscale.avgs="data.frame"))


# Initializes the S4 object
setGeneric("initialize", function(object,  min.cells=3, min.genes=2500, min.counts=10, scale=TRUE, center=TRUE, maxexpr=5000,...) standardGeneric("initialize"))
setMethod("initialize","scDrop",
          function(object,  min.cells=3, min.genes=2500, min.counts=10, scale=TRUE, center=TRUE, maxexpr=5000,...) {
            print("Initializing S4 object")
            
            # Cell filtering
            num.genes = colSums(object@count.data > 0); names(num.genes) = colnames(object@count.data)
            cells.use = names(num.genes[which(num.genes>min.genes)])
            temp.data=object@count.data[,cells.use]
            
            # Gene filtering
            num.cells= rowSums(temp.data > 0)            
            genes.use=names(num.cells[which(num.cells>min.cells)])
            genes.use1 = names(which(rowSums(object@count.data[,cells.use]) > min.counts))
            genes.use = intersect(genes.use, genes.use1)
            maxgene = apply(object@count.data[, cells.use],1, max)
            genes.use2 = rownames(object@count.data)[maxgene < maxexpr]
            genes.use=intersect(genes.use, genes.use2)
            temp.data=temp.data[genes.use,]
            
            # Normalize each library to the median of the transcript counts across all cells 
            # Then, log transform expression values   
            print("Median normalizing counts and log-transforming")
            col.sums=apply(temp.data,2, sum)
            med_trans = median(col.sums)
            norm_counts = med_trans* scale(temp.data, center=FALSE, scale=col.sums)
            rm(temp.data)
            object@data=as.data.frame(log(norm_counts+ 1))

            # Group IDs for each cell
            object@group=factor(unlist(lapply(colnames(object@data),function(x) strsplit(x,"_")[[1]][1] )))
            names(object@group)=colnames(object@data)
            
            print("z-scoring each gene")
            object@cell.names=names(object@group)
            object@scale.data=t(scale(t(object@data),center=center,scale=scale))
            object@scale.data=object@scale.data[complete.cases(object@scale.data),]
            object@meta=data.frame(num.genes[object@cell.names]); colnames(object@meta)[1]="num.genes"
            object@meta[cells.use,"num.trans"] = colSums(object@count.data[,cells.use])
            
            
            object@meta[names(object@group),"sample"]=object@group
            return(object)
          }         
)

setGeneric("doBatchCorrection", function(object,  batch.cov=NULL, max.val=6 ) standardGeneric("doBatchCorrection"))
# Batch Correction using ComBat
# Memory heavy - consider running on a cluster
setMethod("doBatchCorrection","scDrop",
          function(object,   batch.cov=NULL, max.val=6) {
            correct.data = ComBat(object@data,batch.cov, prior.plots=FALSE, par.prior=TRUE)
            correct.data[correct.data > max.val] = max.val
            
            object@batch.data = as.data.frame(correct.data)
            rm(correct.data)
            object@scale.data = t(scale(t(object@batch.data), center=TRUE, scale=TRUE))
            return(object)
          }
)

setGeneric("secondBatchCorrection", function(object, batch.cov =NULL, max.val=6) standardGeneric("secondBatchCorrection"))
setMethod("secondBatchCorrection","scDrop",
          function(object, batch.cov=NULL, max.val=6) {
            sec.correct.data = ComBat(object@batch.data, batch.cov, prior.plots=FALSE, par.prior=TRUE)
            sec.correct.data[sec.correct.data > max.val] = max.val
            
            object@batch.data.fin = as.data.frame(sec.correct.data)
            rm(sec.correct.data)
            object@scale.data = t(scale(t(object@batch.data.fin), center = TRUE, scale = TRUE))
            return(object)
          }
)

setGeneric("doPCA", function(object,pcs.store=100) standardGeneric("doPCA"))
setMethod("doPCA", "scDrop", 
          function(object,pcs.store=100) {
            data.use=object@scale.data
            pc.genes = rownames(object@scale.data)
            
            #Remove genes with zero variation
            pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
            genes.use = pc.genes[pc.genes.var>0]
            pc.data = data.use[genes.use,]
            
            pca.obj = fast.prcomp(t(pc.data),center=FALSE, scale=FALSE)
            object@pca.scores=data.frame(pca.obj$x[,1:pcs.store])
            object@pca.load=data.frame(pca.obj$rotation[,1:pcs.store])
            
            ev=pca.obj$sdev^2
            object@pca.eigenvalues=data.frame(ev)
            return(object)
          }
)
setGeneric("violinplot", function(object,genes.plot,nCol=NULL,ylab.max=12, ...)  standardGeneric("violinplot"))
setMethod("violinplot","scDrop",
          function(object,genes.plot,nCol=NULL,ylab.max=12, ...) {
            if (is.null(nCol)) {
              nCol=1
              if (length(genes.plot)>6) nCol=3
              if (length(genes.plot)>9) nCol=4
            }
            genes.plot1 =intersect(genes.plot, rownames(object@data))
            plot.data = data.frame()
            if (length(genes.plot1) > 0) {
              plot.data = object@data[genes.plot1,]
            }
            
            genes.plot2 = intersect(genes.plot, c("num.genes", "num.trans"))
            if (length(genes.plot2) > 0) {
              plot.data = rbind(plot.data, t(data.frame(object@meta[,genes.plot2])))
            }
            genes.plot = intersect(genes.plot, rownames(plot.data))
            group.use=object@group
            
            #Make individual violins
            pList=lapply(genes.plot,function(x) makeviolin(x,plot.data[x,],group.use))
            
            multiplotList(pList,cols = nCol)
            rp()
            
          }
)  

makeviolin=function(x,data,group) {
  data$x=as.character(rownames(data))
  data.use=data.frame(data[x,])
  group = factor(group, levels=sort(unique(group)))
  if (length(x)==1) {
    data.melt=data.frame(rep(x,length(group))); colnames(data.melt)[1]="x"
    data.melt$value=as.numeric(data[1,1:length(group)])
    data.melt$id=names(data)[1:length(group)]
    }
  
  if (length(x)>1) data.melt=melt(data.use,id="x")
  data.melt$group=group
  
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  p=ggplot(data.melt,aes(factor(group),value),useDingbats=FALSE)
  p2=p + geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(group))) + ylab(x) + xlab("Group/Cluster")
  
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0,size=0.7)
  p4=p3+theme_bw() + ylim(0.8*min(data.melt$value), quantile(data.melt$value,0.99))
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="black", size=10), axis.text.y  = element_text(angle=90, size=8),
               axis.title.x = element_text(face="bold", colour="black", size=10), axis.text.x  = element_text(size=8)))

  return(p5)
}

multiplotList <- function(plots, file, cols=1, layout=NULL) {
  require(grid)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

rp=function() {par(mfrow=c(1,1))}


# Graph clustering
setGeneric("doGraph_clustering", function(object,cells.use=NULL,pcs.use=1:10, num.nn=30, do.jaccard=FALSE, method="Louvain") standardGeneric("doGraph_clustering"))
setMethod("doGraph_clustering", "scDrop", function(object,cells.use=NULL,pcs.use=1:10,num.nn=30, do.jaccard=FALSE, method="Louvain") {
  
  
  if (do.jaccard){
    weights=TRUE;
    method_print = paste0(method,"-","Jaccard")
  } else {
    weights=NULL;
    method_print = method
  }
  
  print(paste0("Performing ", method_print, " clustering. Using ", num.nn, " nearest neighbors, and ", max(pcs.use), " PCs"))
  
  if (is.null(cells.use)){
    data.use=object@pca.scores[,pcs.use]
  } else {
    data.use=object@pca.scores[cells.use,pcs.use]
  } 
  
  Adj = get_edges(data.use,nn=num.nn,do.jaccard=do.jaccard)
  
  
  g=graph.adjacency(Adj, mode = "undirected", weighted=weights)
  if (method=="Louvain") graph.out = cluster_louvain(g)
  if (method=="Infomap") graph.out = cluster_infomap(g)
  
  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
  print("Outputting clusters ..")
  object@meta$clust = NULL
  object@meta[names(clust.assign),"clust"]=clust.assign
  object@group=clust.assign; names(object@group)=names(clust.assign);               
  
  
  return(object) 
}
)

# Build a nearest neighbor graph with or without edge weights, and return an adjacency matrix
get_edges=function(X,nn=30,do.jaccard=TRUE) {
  nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
  

  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1
  
  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  if (do.jaccard){
    
    NN = nearest$nn.idx
    jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges$C = jaccard_dist
    edges = subset(edges, C != 0)
    edges$C = edges$C/max(edges$C)
  }
  
  
  
    
    Adj = matrix(0, nrow=nrow(X), ncol=nrow(X))
    rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
    Adj[cbind(edges$A,edges$B)] = edges$C
    Adj[cbind(edges$B,edges$A)] = edges$C
    return(Adj)
  
}


# Graph clustering with auto-optimization and refinement
setGeneric("perform_refined_clustering", function(object,pcs.use=1:10, corMethod="spearman",iterations=10:120, filePrefix="DropSeqAnalysis", min.newcluster.size=10, min.proportion.outliers=0.10, optimize.NN=TRUE, opt.NN=30) standardGeneric("perform_refined_clustering"))
setMethod("perform_refined_clustering", "scDrop", function(object,pcs.use=1:10, corMethod="spearman",iterations=10:120, filePrefix="DropSeqAnalysis", min.newcluster.size=10, min.proportion.outliers=0.10, optimize.NN=TRUE, opt.NN=30) {

	# Iterate through NNs to find optimal based on average silhouette width. Plot silhouette widths and number of clusters reported
	# This section can be toggled if the number of optimum NN is already known
	if(optimize.NN==TRUE) {
		print("Optimizing number of nearest neighbors (NN) from 10 through 100")
		sils=rep(NA,120)
		numclust=rep(NA,120)
		for (i in iterations) {		# Range of NNs to use is from 10 through 120. First 9 values will be set to NA in final sils vector
			object = doGraph_clustering(object, pcs.use = pcs.use, num.nn = i, do.jaccard = TRUE, method = "Louvain")
			numclust[i] = length(table(object@group))
			data.use=object@pca.scores[,pcs.use]
			data.use=t(data.use)
			sub<-as.numeric(object@group)			
			dm<-as.dist(1-cor(data.use,method=corMethod))
			si<-silhouette(sub,dm,data.use)
			sils[i] = mean(si[,3])
		}
		object@sils = sils
		object@numclust = numclust

		min.sils = round(min(object@sils,na.rm=T),2) - 0.01
		max.sils = round(max(object@sils,na.rm=T),2) + 0.01
		min.clust = round(min(object@numclust,na.rm=T),2) - 5
		max.clust = round(max(object@numclust,na.rm=T),2) + 5
		opt.NN = which.max(object@sils)

		pdfName = paste0(filePrefix,"_sils_numclust_iter_",iterations,corMethod,".pdf")
		pdf(file=pdfName)
		print(ggplot(as.data.frame(object@sils), aes(x=seq(1:120),y=object@sils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,120,10)),limits=c(0,120),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
		print(ggplot(as.data.frame(object@numclust), aes(x=seq(1:120),y=object@numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,120,10)),limits=c(0,120),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
		dev.off()
	}
	
	# Re-compute clusters using optimized number of NN
	object = doGraph_clustering(object, pcs.use = pcs.use, num.nn = opt.NN, do.jaccard = TRUE, method = "Louvain")
	nclust=length(table(object@group))

	# Compute silhouette widths for initial cluster assignments
	data.use=object@pca.scores[,pcs.use]
	data.use=t(data.use)
	sub<-as.numeric(object@group)
	dm<-as.dist(1-cor(data.use,method=corMethod))
	si<-silhouette(sub,dm,data.use)
	pdfName = paste0(filePrefix,"_",opt.NN,"NN_initial_silhouettes.pdf")
	pdf(file=pdfName,height=12)
	plot(si,col=sort(as.integer(object@group)))
	dev.off()
	fileName = paste0(filePrefix,"_",opt.NN,"NN_initial_silhouettes.txt")
	write.table(si,fileName,sep="\t",col.names=NA)

	
	# Iteratively refine cluster assignments using silhouette widths for each cluster
	# Allow for cells to be reassigned to their next closest neighbor if they are a cluster outlier
	# Allow for new clusters to be formed given sufficiently large group of outlier cells and that they themselves form a robust cluster, otherwise retain initial cluster assignment until next round
	
	print("Refining clusters")
	iter=5
	newclust = rep(0,length(as.numeric(object@group)))
	names(newclust) = names(object@group)
	nclust=seq(1:length(levels(object@group)))
	totalclust = length(nclust)
	origclust = length(levels(object@group))
	addedclusters=numeric(0)
	flagged = rep("FALSE",length(as.numeric(object@group)))
	names(flagged) = names(object@group)

	for (n in 1:iter) {
		print(paste0("Cluster refinement iteration ",n, " of ",iter))
		data.use=object@pca.scores[,pcs.use]
		data.use=t(data.use)
		if(n==1) {
			sub<-as.numeric(object@group)
		} else {
			sub<-as.numeric(newclust)
		}
		dm<-as.dist(1-cor(data.use,method=corMethod))
		si<-silhouette(sub,dm,data.use)
		rownames(si) = names(object@group)
		if(length(addedclusters)>0) {
			for (k in addedclusters) {
				if (length(rownames(si)[si[,1]==k]) > 1) {
					subsil = si[si[,1]==k,]
					if( (mean(subsil[,3]) < 0) || (length(subsil[,3]) < min.newcluster.size) ) {			
						for (i in 1:length(subsil[,3])) {
							newclust[rownames(subsil)[i]] = as.numeric(object@group[rownames(subsil)[i]])
							flagged[rownames(subsil)[i]] = "TRUE"
						}
						addedclusters = addedclusters[-which(addedclusters==k)]
						nclust = nclust[-which(nclust==k)]
						sub<-as.numeric(newclust)
						dm<-as.dist(1-cor(data.use,method=corMethod))
						si<-silhouette(sub,dm,data.use)
						rownames(si) = names(object@group)
					}
				} else {
					newclust[rownames(si)[si[,1]==k]] = as.numeric(object@group[rownames(si)[si[,1]==k]])
					flagged[rownames(si)[si[,1]==k]] = "TRUE"
				}
			}
		}
		pdfName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_round",n-1,".pdf")
		pdf(file=pdfName,height=12)
		if(n==1) {
			plot(si,col=sort(as.integer(object@group)))
		} else {
			plot(si,col=(1+sort(newclust)))
		}
		dev.off()
		fileName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_round",n-1,".txt")
		write.table(si,fileName,sep="\t",col.names=NA)

		for (a in nclust) {
			addclust = "FALSE"
			if (length(rownames(si)[si[,1]==a]) > 1) {
				subsil = si[si[,1]==a,]
				propneg = length(subsil[,3][subsil[,3]<0]) / length(subsil[,3])
				if ( (propneg > min.proportion.outliers) && (length(subsil[,3][subsil[,3]<0]) >= min.newcluster.size) && (n==(iter-1))) {
					for (i in 1:length(subsil[,3])) {
						if ( (subsil[i,3] < 0) && (flagged[rownames(subsil)[i]] == "FALSE") ) {
							newclust[rownames(subsil)[i]] = totalclust + 1
							addclust = "TRUE"
						} else {
							newclust[rownames(subsil)[i]] = subsil[i,1]
						}
					}
				} else {
					for (i in 1:length(subsil[,3])) {
						if (subsil[i,3] < -0.1) {
							newclust[rownames(subsil)[i]] = subsil[i,2]
						} else {
							newclust[rownames(subsil)[i]] = subsil[i,1]
						}
					}
				}
			} else {
				newclust[rownames(si)[si[,1]==a]] = as.numeric(object@group[rownames(si)[si[,1]==a]])
				flagged[rownames(si)[si[,1]==a]] = "TRUE"
			}
			if (addclust=="TRUE") {
				totalclust = totalclust + 1
				addedclusters = c(addedclusters,totalclust)
				nclust = c(nclust,totalclust)
			}
		}
	}
	
	# Write and plot out last round of silhouette widths
	data.use=object@pca.scores[,pcs.use]
	data.use=t(data.use)
	sub<-as.numeric(newclust)
	dm<-as.dist(1-cor(data.use,method=corMethod))
	si<-silhouette(sub,dm,data.use)
	rownames(si) = names(object@group)
	pdfName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_round",n,".pdf")
	pdf(file=pdfName,height=12)
	plot(si,col=(1+sort(newclust)))
	dev.off()
	fileName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_round",n,".txt")
	write.table(si,fileName,sep="\t",col.names=NA)


	# If new clusters were formed, their numbering may not be consecutive with original cluster labels. Re-number the newly created ones to be consecutive
	print("Finalizing cluster refinement")
	names(newclust) = names(object@group)
	newclust=as.factor(newclust)
	counter=0
	for (z in addedclusters){
		levels(newclust)[levels(newclust)==z] <- origclust + 1 + counter
		counter = counter + 1
	}

	object@group = newclust

	# Make one last round of silhouette plots using final cluster assignments
	data.use=object@pca.scores[,pcs.use]
	data.use=t(data.use)
	sub<-as.numeric(object@group)
	dm<-as.dist(1-cor(data.use,method=corMethod))
	si<-silhouette(sub,dm,data.use)
	rownames(si) = names(object@group)
	pdfName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_final.pdf")
	pdf(file=pdfName,height=12)
	plot(si,col=sort(as.integer(object@group)))
	dev.off()
	fileName = paste0(filePrefix,"_",opt.NN,"NN_silhouettes_final.txt")
	write.table(si,fileName,sep="\t",col.names=NA)

	print(paste0("Cluster refinement process completed. ",length(levels(object@group))," clusters remain")) 
	return(object)
}
)





#Visualize single cells as clusters in a tSNE plot
plot.tsne=function(object,legend.position="topright",xlimits=NULL,ylimits=NULL,scaling=0.5) {
  
  cols=rainbow(length(levels(object@group))); cols[1]="lightgrey"
  group=as.numeric(object@group)
  
  order = sample(c(1:dim(object@tsne.y)[1]), replace=FALSE)
  
  if(is.null(xlimits) & is.null(ylimits)) {
    xmin=floor(min(object@tsne.y[,1]))
    xmax=ceiling(max(object@tsne.y[,1]))
    ymin=floor(min(object@tsne.y[,2]))
    ymax=ceiling(max(object@tsne.y[,2]))
    xlimits=c(xmin,xmax)
    ylimits=c(ymin,ymax)
  }
  
  plot(object@tsne.y[,1],object@tsne.y[,2],col=cols[as.integer(object@group)],pch=16,xlab="tSNE1",ylab="tSNE2",cex=scaling,xlim=xlimits,ylim=ylimits)
  k.centers=t(sapply(levels(object@group),function(x) apply(object@tsne.y[names(object@group[which(object@group %in% x)]),],2,mean)))
  points(k.centers[,1],k.centers[,2],cex=1.3,col="white",pch=16); text(k.centers[,1],k.centers[,2],levels(object@group),cex=1.25)
  legend(legend.position,pch=rep(16,length(as.numeric(levels(object@group)))),col=cols[unique(sort(as.integer(object@group)))],paste0("cluster",seq(1,length(as.numeric(levels(object@group))))))
}

# scatter plot coloring genes by their expression levels
setGeneric("gene.expression.scatter", function(object, genes,cells.use=NULL,cols.use=Greys(10),pch.use=16,nCol=NULL, xlim=NULL, ylim=NULL) standardGeneric("gene.expression.scatter"))
setMethod("gene.expression.scatter", "scDrop", 
          function(object, genes,cells.use=NULL,cols.use=terrain.colors( 10),pch.use=16,nCol=NULL, xlim=NULL, ylim=NULL) {
            
            if (is.null(cells.use)){
              cells.use = colnames(object@data)
            } else {
              cells.use = cells.use[cells.use %in% colnames(object@data)]
            }
          
            
            if (is.null(nCol)) {
              nCol=2
              if (length(genes)>6) nCol=3
              if (length(genes)>9) nCol=4
            }         
            num.row=floor(length(genes)/nCol-1e-5)+1
            
            par(mfrow=c(num.row,nCol))
            data.plot=object@tsne.y[cells.use,]
           
            group.id=as.factor(object@group[cells.use])
            data.plot$group=group.id
            x1="tSNE1"; x2="tSNE2"
            data.plot$x=data.plot[,1]; data.plot$y=data.plot[,2]
            data.plot$pt.size=1
            for(i in genes) {
              data.gene=as.numeric(object@data[i,cells.use])
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = length(cols.use))))
              data.col=rev(cols.use)[data.cut]
              plot(data.plot$x,data.plot$y,col=data.col,cex=0.7,pch=16,main=i,xlab=x1,ylab=x2, xlim=xlim, ylim=ylim)
            }
            rp()
          }
)

# NEW

# Binomial test to evaluate differentially expressed genes between two clusters
# if only one cluster is provided, then it will be compared against the rest of the cells
#setGeneric("markers.binom", function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom"))
#setMethod("markers.binom", "scDrop",
#          function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
#            genes.use=rownames(object@data)
#            clust.use=object@group
#            cells.1=names(clust.use[which(clust.use%in%clust.1)])
#            
#            if (is.vector(clust.2)) {
#              cells.2=names(clust.use[which(clust.use%in%clust.2)])
#            }
#            else if (is.null(clust.2)) {
#              clust.2="rest"
#              cells.2=names(clust.use)
#              cells.2=cells.2[!(cells.2%in%cells.1)]
#            } else {
#              cells.2=names(clust.use[which(clust.use%in%clust.2)])
#            }
#		            
#            Count.mat = object@count.data
#            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
#            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
#            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
#            
#            
#            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
#            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
#            
#            if (is.integer(clust.2) | length(clust.2)>1) {
#              genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
#			}
#			else {
#              genes.include = posFrac.1 >= 0.1			
#			}            
#           
#            result = result[genes.include,]
#            result = result[order(abs(result$log.effect), decreasing=TRUE),]
#            
#            #Mean number of transcripts per cell
#            if (!is.null(attr(object,"count.data"))){
#              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
#              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
#              result[,paste0("nTrans_", clust.1)] = nTrans.1
#              result[, paste0("nTrans_", clust.2)] = nTrans.2
#            }
#            
#            return(result)
#          } 
#)
#

# ORIGINAL

setGeneric("markers.binom", function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom"))
setMethod("markers.binom", "scDrop",
          function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@group
            cells.1=names(clust.use[which(clust.use%in%clust.1)])
            
            if (is.null(clust.2)) {
              clust.2="rest"
              cells.2=names(clust.use)
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }
            
            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
            
            
            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
            
            if (clust.2=="rest"){
              genes.include = posFrac.1 >= 0.1
            } else{
              genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
            }
            
            result = result[genes.include,]
            result = result[order(abs(result$log.effect), decreasing=TRUE),]
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              result[,paste0("nTrans_", clust.1)] = nTrans.1
              result[, paste0("nTrans_", clust.2)] = nTrans.2
            }
            
            return(result)
          } 
)


# Binomial test to evaluate differentially expressed genes between two clusters
# if only one cluster is provided, then it will be compared against the rest of the cells
setGeneric("markers.binom.multi", function(object, clust.1,clust.2,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom.multi"))
setMethod("markers.binom.multi", "scDrop",
          function(object, clust.1,clust.2,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@group

            cells.1=names(clust.use[which(clust.use%in%clust.1)])
            cells.2=names(clust.use[which(clust.use%in%clust.2)])
		            
            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
        
            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
            
            genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
			
			result = result[genes.include,]
            result = result[order(abs(result$log.effect), decreasing=TRUE),]
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              result[,paste0("nTrans_", clust.1)] = nTrans.1
              result[, paste0("nTrans_", clust.2)] = nTrans.2
            }
            
            return(result)
          } 
)





# Binomial test to evaluate differentially expressed genes between two genotypes within the same cluster

setGeneric("markers.binom.geno", function(object, groupName.1, fullReport=FALSE, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom.geno"))
setMethod("markers.binom.geno", "scDrop",
          function(object, groupName.1, fullReport=FALSE, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@group
            cells.1=grep(groupName.1,names(clust.use[which(clust.use%in%clust.1)]),value=TRUE)
            
            if (is.null(clust.2)) {
              clust.2=clust.1
              cells.2=names(clust.use)
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }
            
            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
            
            if(fullReport==FALSE) {
            	posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            	posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
            
            	genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
            
            	result = result[genes.include,]
            	result = result[order(abs(result$log.effect), decreasing=TRUE),]
            
            	#Mean number of transcripts per cell
            	if (!is.null(attr(object,"count.data"))){
            	  nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              	  nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              	  result[,paste0("nTrans_", groupName.1)] = nTrans.1
              	  result[, paste0("nTrans_base")] = nTrans.2
            	}
            	return(result)
            }
            else{
            	result = result[order(abs(result$log.effect), decreasing=TRUE),]
            	
            	if (!is.null(attr(object,"count.data"))){
            	  nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              	  nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              	  result[,paste0("nTrans_", groupName.1)] = nTrans.1
              	  result[, paste0("nTrans_base")] = nTrans.2
            	}
            	return(result)
            }
          } 
)

###TTEST#####
setGeneric("ttest.geno.sub", function(object, groupName.1,clust.1,clust.2=NULL,data=NULL,effect.size=log(2)) standardGeneric("ttest.geno.sub"))
setMethod("ttest.geno.sub", "scDrop",
          function(object, groupName.1, clust.1,clust.2=NULL,data=NULL,effect.size=log(2)) {
            data=data[,which(object@subgroup==clust.1)]
            keep<-rowSums(data>0) >= ncol(data)*.25 #expressed in 25% of subclust cells
            data=data[keep,]

	    genes.use=rownames(data)
            clust.use=object@subgroup
            cells.1=grep(groupName.1,names(clust.use[which(clust.use%in%clust.1)]),value=TRUE)

            if (is.null(clust.2)) {
              clust.2=clust.1
              cells.2=names(clust.use)
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }
	    genes<-rownames(TPM.mat[rownames(data),])
            if (is.null(data)) data = object@scale.data[, c(cells.1, cells.2)]
	    df<-data.frame(row.names=genes)
	    df$pval=0
	    df$qval=0
	    df$LFC=0
            for(j in 1:length(genes.use)) {
            	a<-t.test(data[j,cells.1],data[j,cells.2])
            	df[j,1]=a$p.value
            	cells.1.mean=mean(as.numeric(data[j,cells.1]))
		cells.2.mean=mean(as.numeric(data[j,cells.2]))
		logfc<-log2((cells.1.mean+.1)/(cells.2.mean+.1))
            	df[j,3]=logfc
            }
	p<-df$pval
	q<-p.adjust(p,method='fdr')
	df[,2]=q
			
return(df)
}
)



# Binomial test to evaluate differentially expressed genes between two genotypes within the same cluster

setGeneric("markers.binom.geno.sub", function(object, groupName.1, fullReport=TRUE, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom.geno.sub"))
setMethod("markers.binom.geno.sub", "scDrop",
          function(object, groupName.1, fullReport=TRUE, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@subgroup
            cells.1=grep(groupName.1,names(clust.use[which(clust.use%in%clust.1)]),value=TRUE)

            if (is.null(clust.2)) {
              clust.2=clust.1
              cells.2=names(clust.use)
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }

            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)

            if(fullReport==FALSE) {
                posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
                posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))

                genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)

                result = result[genes.include,]
                result = result[order(abs(result$log.effect), decreasing=TRUE),]

		#Mean number of transcripts per cell
                if (!is.null(attr(object,"count.data"))){
                  nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
                  nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
                  result[,paste0("nTrans_", groupName.1)] = nTrans.1
                  result[, paste0("nTrans_base")] = nTrans.2
                }
                return(result)
            }
            else{
                result = result[order(abs(result$log.effect), decreasing=TRUE),]

                if (!is.null(attr(object,"count.data"))){
                  nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
                  nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
                  result[,paste0("nTrans_", groupName.1)] = nTrans.1
                  result[, paste0("nTrans_base")] = nTrans.2
                }
                return(result)
            }
          }
)


#new for subclusters
setGeneric("markers.binom.subclust", function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom.subclust"))
setMethod("markers.binom.subclust", "scDrop",
          function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@subgroup
            cells.1=names(clust.use[which(clust.use%in%clust.1)])

            if (is.null(clust.2)) {
              clust.2="rest"
              cells.2=names(clust.use)
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }

            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)


            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))

            if (clust.2=="rest"){
              genes.include = posFrac.1 >= 0.1
            } else{
              genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
            }

            result = result[genes.include,]
            result = result[order(abs(result$log.effect), decreasing=TRUE),]

            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              result[,paste0("nTrans_", clust.1)] = nTrans.1
              result[, paste0("nTrans_", clust.2)] = nTrans.2
            }

            return(result)
          }
)

# Binomial test to evaluate differentially expressed genes between two clusters
# if only one cluster is provided, then it will be compared against the rest of the cells
setGeneric("markers.binom.multi.sub", function(object, clust.1,clust.2,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom.multi.sub"))
setMethod("markers.binom.multi.sub", "scDrop",
          function(object, clust.1,clust.2,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@subgroup

            cells.1=names(clust.use[which(clust.use%in%clust.1)])
            cells.2=names(clust.use[which(clust.use%in%clust.2)])
		            
            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
        
            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
            
            genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
			
			result = result[genes.include,]
            result = result[order(abs(result$log.effect), decreasing=TRUE),]
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              result[,paste0("nTrans_", clust.1)] = nTrans.1
              result[, paste0("nTrans_", clust.2)] = nTrans.2
            }
            
            return(result)
          } 
)



setGeneric("binomcount.test", function(object, cells.1,cells.2, effect.size, TPM.mat, Count.mat) standardGeneric("binomcount.test"))
setMethod("binomcount.test", "scDrop",
          function(object, cells.1,cells.2, effect.size, TPM.mat, Count.mat) {
            
            x=TPM.mat
            y=Count.mat
            
            #Test for enrichments in cluster #1
            m = apply(x[, cells.2], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #2
            m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
            n = apply(x[, cells.1], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #1
            #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
            pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            
           log_fold_express = log((n+1)*length(cells.2)/((m+1)*length(cells.1))) #log proportion of expressing cells

            d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
            if(!is.na(effect.size)) {
            	d1 <- subset(d1, log.effect >= effect.size)
            }
            d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
            #Enrichments in cells.2
            n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
            #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
            pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
            d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
            if(!is.na(effect.size)) {
            	d2 <- subset(d2, log.effect <= -effect.size)
            }            
            d2 <- d2[order(d2$pval,decreasing=FALSE),]
            
            d = rbind(d1, d2);
            d = d[order(d$pval, decreasing=FALSE),]
            return(d)
          } 
)



setGeneric("wrs.test", function(object, nclust, name.1) standardGeneric("wrs.test"))
setMethod("wrs.test", "scDrop",
       function(object, nclust, name.1) {

		wil.mat=matrix(0,nrow=length(rownames(object@data)),ncol=(2*nclust))
		colnames(wil.mat)=rep(c("Wilcoxon.P","log2fc"),nclust)
		rownames(wil.mat)=rownames(object@data)

		genes.use=rownames(object@data)
		clust.use=object@group

		for(i in 1:length(rownames(object@data))) {
			for (j in 1:nclust) {
				cells.1=grep(name.1,names(clust.use[which(clust.use%in%j)]),value=TRUE)
				cells.2=names(clust.use[which(clust.use%in%j)])
				cells.2=cells.2[!(cells.2%in%cells.1)]
				TPM.mat = object@scale.data

				wil=wilcox.test(as.numeric(TPM.mat[i,cells.1]),as.numeric(TPM.mat[i,cells.2]),paired=F,alternative="two.sided",exact=FALSE)
				p=wil$p.value

				cells.1.mean=mean(as.numeric(TPM.mat[i,cells.1]))
				cells.2.mean=mean(as.numeric(TPM.mat[i,cells.2]))
	
				logfc=log2((cells.1.mean+0.1)/(cells.2.mean+0.1))
			
				wil.mat[i,((j*2)-1)] = p
				wil.mat[i,(j*2)] = logfc
			}
		}
		return(wil.mat)
   }
)

setGeneric("wrs.test.sub", function(object, clust, name.1,data) standardGeneric("wrs.test.sub"))
setMethod("wrs.test.sub", "scDrop",
       function(object, clust, name.1,data) {

		data=data[,which(object@subgroup==clust)]
		keep<-rowSums(data>0) >= ncol(data)*.25 #expressed in 25% of subclust cells
		data=data[keep,]
                wil.mat=matrix(0,nrow=length(rownames(data)),ncol=3)
                colnames(wil.mat)=c("Wilcoxon.P",'Wilcoxon.Q',"log2fc")
                rownames(wil.mat)=rownames(data)

                genes.use=rownames(data)
		cells.1=grep(name.1,colnames(data))
                cells.2=colnames(data)                             
                cells.2=cells.2[!(cells.2%in%cells.1)]
	
                for(i in 1:length(rownames(data))) {
                                wil=wilcox.test(as.numeric(data[i,cells.1]),as.numeric(data[i,cells.2]),paired=F,alternative="two.sided",exact=FALSE)
                                p=wil$p.value
                                cells.1.mean=mean(as.numeric(data[i,cells.1]))
                                cells.2.mean=mean(as.numeric(data[i,cells.2]))

                                logfc=log2((cells.1.mean+0.1)/(cells.2.mean+0.1))
				wil.mat[i,1]=p
				wil.mat[i,3]=logfc
                }
		p<-wil.mat[,1]
		q<-p.adjust(p,method='fdr')
		wil.mat[,2]<-q
                return(wil.mat)
   }
)



setGeneric("merge.clusters.DE", function(object, min.de.genes = 25, effect.size=log(2), pval.cutoff=0.01, pcs.use=1:10,TPM.mat=NULL, Count.mat=NULL) standardGeneric("merge.clusters.DE"))
setMethod("merge.clusters.DE", "scDrop", 
          function(object,min.de.genes=25, effect.size=log(2),pval.cutoff=0.01, pcs.use=1:10, TPM.mat=NULL, Count.mat=NULL) {
          
            genes.use = rownames(object@data)
            clust.test = as.numeric(levels(object@group))
           # if (is.null(tag)){
           #   filename = "CLUSTER_PAIRWISE_MARKERS.txt"
           # } else {
           #   filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
           # }
           # zz = file(filename,open="wt")
            
            
            num.clust=length(clust.test) 
            print(paste0("Starting with ", num.clust, " clusters"))
            
            pass.thresh=1e6*data.frame(diag(length(levels(object@group)))); 
            
            for (i in setdiff(as.numeric(levels(object@group)), clust.test)){
              pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
            } 
            
            dist.clust = pass.thresh
            
            #Find the number of differentially expressed genes between every pair of clusters
            for(k in 1:num.clust) {
              i=clust.test[k]
              print(paste0("Testing Cluster ", i))
              for(m in ((k+1):num.clust)) {
                j=clust.test[m]
                #print(j)
                if (m>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=markers.binom(object,i,j,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  
                  marker.pass=subset(marker,pval<pval.cutoff)
                  #print(paste("Test b/w Clusters ",i,j, "-# DE = ",nrow(marker.pass)))
                  #print(head(subset(marker.pass, log.effect > 0),5))
                  #print(head(subset(marker.pass, log.effect < 0),5))
                  
                  num.de.genes = 2*min(nrow(subset(marker.pass, log.effect > 0)), nrow(subset(marker.pass, log.effect < 0)))
                  pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
              print(pass.thresh[i,])
              
            }
            
            colnames(pass.thresh) = levels(object@group)
            rownames(pass.thresh) = levels(object@group)
            
            write.table(pass.thresh, file=paste0("DE_genes_matrix_2.txt"), sep="\t", quote=FALSE)
            
            #iteratively merge clusters
            min.val = min(pass.thresh)
            min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            merge.ind=-1
            while(min.val <= min.de.genes) {
              merge.ind=merge.ind+1
              
              #In case of ties, merge clusters that are closest in PC space
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
              
              if (pass.thresh[test.1,test.2]<= min.de.genes) {
                object@group[which(object@group==test.2)]=test.1
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
                old.group.levels = as.numeric(levels(object@group))
                old.group.levels = setdiff(old.group.levels, test.2)
                clust.test = setdiff(clust.test, test.2)
                
                object@group = droplevels(object@group)
                levels(object@group) = c(1:length(levels(object@group)))
                object@meta[,"clust"] = object@group
                
                new.group.levels = as.numeric(levels(object@group))
                names(new.group.levels) = as.character(old.group.levels)
                clust.test = new.group.levels[as.character(clust.test)]
                
                
                
                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(clust.test, test.1)){
                  print(i)
                  marker= markers.binom(object,test.1,i,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  marker.pass=subset(marker,pval<pval.cutoff)
                  pass.thresh[test.1,i]=2*min(nrow(subset(marker.pass, log.effect>0)),nrow(subset(marker.pass, log.effect<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  #pass.thresh[test.1,i]=nrow(marker.pass); 
                  pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:length(levels(object@group))
              rownames(pass.thresh) = colnames(pass.thresh)
              
              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)
              
            }
            return(object)
          }
)

setGeneric("ComputeClusterDistances", function(object, reduction.use="pca", dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ComputeClusterDistances"))
setMethod("ComputeClusterDistances", "scDrop", 
          function(object,reduction.use="pca",dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use =  colnames(object@data)
            group.use=object@group[cells.use]
            if (reduction.use == "pca"){
              data.use = object@pca.scores[cells.use,pcs.use]
              centroids = ClusterCentroids(object, reduction.use="pca", pcs.use=pcs.use, cells.use=cells.use)
            }
            
            if (dist.type=="centroid"){
              clust.dists = as.matrix(dist(centroids, upper=TRUE))
              diag(clust.dists) = 1e6
            }
            
            num.clust = length(levels(group.use))
            
            
            if (dist.type == "nn"){
              clust.dists = matrix(0, nrow=num.clust, ncol=num.clust)
              diag(clust.dists) = 1e6
              rownames(clust.dists) = levels(group.use)
              colnames(clust.dists) = rownames(clust.dists)
              for (i in 1:nrow(clust.dists)){
                for(j in ((i+1):ncol(clust.dists))){
                  if (j>nrow(clust.dists)) break
                  cells.in.cluster_i = names(object@group)[object@group %in% i]
                  cells.in.cluster_i = cells.in.cluster_i[cells.in.cluster_i %in% cells.use]
                  cells.in.cluster_j = names(object@group)[object@group %in% j]
                  cells.in.cluster_j = cells.in.cluster_j[cells.in.cluster_j %in% cells.use]
                  
                  nnA = nn2(data.use[cells.in.cluster_i,], query = centroids[j,], k=1)
                  nnB = nn2(data.use[cells.in.cluster_j,], query = centroids[i,],k=1)
                  clust.dists[i,j] = min(c(nnA$nn.dists, nnB$nn.dists))
                  
                  clust.dists[j,i] = clust.dists[i,j]
                }
              }
            }
            
            colnames(clust.dists) = c(1:ncol(clust.dists))
            rownames(clust.dists) = colnames(clust.dists)
            return(clust.dists)
          }
          
          
)

setGeneric("ClusterCentroids", function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ClusterCentroids"))
setMethod("ClusterCentroids", "scDrop", 
          function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use = colnames(object@data)
            group.use=object@group[cells.use]
            if (reduction.use == "pca"){
              data.use = object@pca.scores[cells.use,pcs.use]
            }
            
            centroids = c()
            for (i in levels(group.use)){
              cells.in.cluster = names(object@group)[object@group %in% i]
              cells.in.cluster = cells.in.cluster[cells.in.cluster %in% cells.use]
              centroids = rbind(centroids, colMeans(data.use[cells.in.cluster,]))
            }
            centroids = as.data.frame(centroids)
            colnames(centroids) = colnames(data.use)
            rownames(centroids) = as.numeric(levels(object@group))
            
            return(centroids)
          }
          
          
)

setGeneric("dot.plot", function(object,features.use=NULL, group.use=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL,family=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot"))
setMethod("dot.plot", "scDrop", 
          function(object,features.use=NULL,group.use=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL,family=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            if (is.null(group.use)) group.use = levels(object@group)
            if (is.null(group.names)) {
            	group.names = group.use
            }
            else {
            	group.names = paste0(group.names," [",group.use,"]")
            }
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in group.use){
              cells.in.cluster = names(object@group)[which(object@group== i)]
              vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names

            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
              
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) + theme(text=element_text(family=family))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) + theme(text=element_text(family=family))
              print(p)
              
              
            }
              
          }
)

setGeneric("dot.plot.gt", function(object,features.use=NULL, group.use=NULL,geno.name=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot.gt"))
setMethod("dot.plot.gt", "scDrop", 
          function(object,features.use=NULL,group.use=NULL,geno.name=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            if (is.null(group.use)) group.use = levels(object@group)
            if (is.null(group.names)) group.names = group.use
            
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in group.use){
			  if(is.null(geno.name)) {
                cells.in.cluster = names(object@group)[which(object@group== i)]
              }
              else{
                cells.in.cluster = grep(geno.name,names(object@group)[which(object@group== i)],value=TRUE)
			  }
              vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names

            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
              
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic"))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) 
              print(p)
              
              
            }
              
          }
)

setGeneric("dot.plot.sub", function(object,features.use=NULL, group.use=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,family=NULL,min.perc=0,...) standardGeneric("dot.plot.sub"))
setMethod("dot.plot.sub", "scDrop", 
          function(object,features.use=NULL,group.use=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,family=NULL,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            if (is.null(group.use)) group.use = levels(object@subgroup)
            if (is.null(group.names)) {
            	group.names = group.use
            }
            else {
            	group.names = paste0(group.names," [",group.use,"]")
            }
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in group.use){
              cells.in.cluster = names(object@subgroup)[which(object@subgroup== i)]
              vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names

            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
              
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) +  theme(text=element_text(family=family))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) +  theme(text=element_text(family=family))
              print(p)
              
              
            }
              
          }
)

######dotplo sham vs sni subclusts###########; 2x rows; 1 control, 1 experimental. multiple genes
### for condition, use a list of grepable names; I.E. ...condition=c('WT','Top1'),
####Dotplot function for MULTIPLE genes, results in 2x rows; control and experimental; geno should be grepable. i.e. geno.name=c('control','experiment')
setGeneric("dot.plot.subgt", function(object,features.use=NULL,group.use=NULL,geno.name=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot.subgt"))
setMethod("dot.plot.subgt", "scDrop", function(object,features.use=NULL,group.use=NULL,geno.name=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {     
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            print(features.use)
            if (is.null(group.use)) group.use = levels(object@subgroup)
            #if (is.null(group.names)) group.names = group.use
            
            
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            for (j in group.use) {
            	for (i in geno.name) {
            		cells.in.cluster=grep(i,names(object@subgroup)[which(object@subgroup==j)],value=TRUE)
            		vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x))
            		PercMat=cbind(PercMat,vec.exp)
            		vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
            		ExpMat = cbind(ExpMat, vec.exp)
            		colnames(ExpMat)[length(colnames(ExpMat))]=paste0(i,j)
            		colnames(PercMat)[length(colnames(PercMat))]=paste0(i,j)
            		#print('This is ExpMat')
            		#print(ExpMat)
            	}
            }
            print(colnames(PercMat))
            print(dim(PercMat))
            print(dim(ExpMat))
            group.names=colnames(ExpMat)
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
            print(ExpVal)  
            print(do.transpose)
            if (!do.transpose==TRUE){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic"))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) 
              print(p)
              
              
            }
              
          }
)



####conditional dotplot for a single gene. Rows = clusters, 2 columns = conditions
###DEFAULT USES SUBCLUSTS
setGeneric("dot.plot.subgt.singlegene", function(object,features.use=NULL,group.use=NULL,geno.name=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot.subgt.singlegene"))
setMethod("dot.plot.subgt.singlegene", "scDrop", function(object,features.use=NULL,group.use=NULL,geno.name=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {         
         
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            print(features.use)
            if (is.null(group.use)) group.use = levels(object@subgroup)
            #if (is.null(group.names)) group.names = group.use
            
            
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            for (j in group.use) {
            	for (i in geno.name) {
            		cells.in.cluster=grep(i,names(object@subgroup)[which(object@subgroup==j)],value=TRUE)
            		vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x))
            		PercMat=cbind(PercMat,vec.exp)
            		vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
            		ExpMat = cbind(ExpMat, vec.exp)
            		colnames(ExpMat)[length(colnames(ExpMat))]=paste0(j)
            		colnames(PercMat)[length(colnames(PercMat))]=paste0(j)
            		#print('This is ExpMat')
            		#print(ExpMat)
            	}
            }
            print(colnames(PercMat))
            print(dim(PercMat))
            print(dim(ExpMat))
            print(ExpMat)
            #group.names=colnames(ExpMat)
            #rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            #PercMat = PercMat[rows.use,]
            #ExpMat = ExpMat[rows.use,]
            #features.use = rows.use
            #if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            #if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
            print(do.transpose)
            ExpVal$Condition=rep(c('Sham','SNI'),nrow(ExpVal)/2)
            ExpVal$cluster=as.character(ExpVal$cluster)
            group.names=as.character(ExpVal$cluster)
            print(ExpVal)
        
            
            if (!do.transpose==TRUE){
              ExpVal$Condition = factor(ExpVal$Condition, levels=geno.name)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ExpVal$cluster))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(Condition)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic"))
              print(p)
            } else {
              ExpVal$Condition = factor(ExpVal$Condition, levels=rev(geno.name))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(Condition),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) 
              print(p)
              
              
            }
              
          }
)



# Graph clustering with auto-optimization and refinement
setGeneric("silopt", function(object,pcs.use=1:10, corMethod="pearson",iterations=10:120, filePrefix="DropSeqAnalysis", min.newcluster.size=10, min.proportion.outliers=0.10, optimize.NN=TRUE, opt.NN=30) standardGeneric("silopt"))
setMethod("silopt", "scDrop", function(object,pcs.use=1:10, corMethod="pearson",iterations=10:120, filePrefix="DropSeqAnalysis", min.newcluster.size=10, min.proportion.outliers=0.10, optimize.NN=TRUE, opt.NN=30) {

	# Iterate through NNs to find optimal based on average silhouette width. Plot silhouette widths and number of clusters reported
	# This section can be toggled if the number of optimum NN is already known
	if(optimize.NN==TRUE) {
		print("Optimizing number of nearest neighbors (NN) from 10 through 120")
		sils=rep(NA,1000)
		numclust=rep(NA,1000)
		for (i in iterations) {		# Range of NNs to use is from 10 through 100. First 9 values will be set to NA in final sils vector
			object = doGraph_clustering(object, pcs.use = pcs.use, num.nn = i, do.jaccard = TRUE, method = "Louvain")
			numclust[i] = length(table(object@group))
			data.use=object@pca.scores[,pcs.use]
			data.use=t(data.use)
			sub<-as.numeric(object@group)
			dm<-as.dist(1-cor(data.use,method=corMethod))
			si<-silhouette(sub,dm,data.use)
			sils[i] = mean(si[,3])
			print(sils[i])
		}
		
		object@sils = sils
		object@numclust = numclust

		min.sils = round(min(object@sils,na.rm=T),2) - 0.01
		max.sils = round(max(object@sils,na.rm=T),2) + 0.01
		min.clust = round(min(object@numclust,na.rm=T),2) - 5
		max.clust = round(max(object@numclust,na.rm=T),2) + 5
		opt.NN = which.max(object@sils)
		print("The highest silhouette score is")
		print(max.sils)
		print("opt NN is ")
		print(opt.NN)
		pdfName = paste0(filePrefix,"_sils_numclust_iter",min(iterations),'to',max(iterations),'_',corMethod,".pdf")
		pdf(file=pdfName)
		print(ggplot(as.data.frame(object@sils), aes(x=seq(1:1000),y=object@sils),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.5) + scale_x_continuous(breaks=c(seq(0,1000,50)),limits=c(0,1000),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.sils,max.sils,0.01)),limits=c(min.sils,max.sils),minor_breaks=NULL))
		print(ggplot(as.data.frame(object@numclust), aes(x=seq(1:1000),y=object@numclust),useDingbats=FALSE) + geom_point() + geom_smooth(span=0.2) + scale_x_continuous(breaks=c(seq(0,1000,50)),limits=c(0,1000),minor_breaks=NULL) + scale_y_continuous(breaks=c(seq(min.clust,max.clust,5)),limits=c(min.clust,max.clust),minor_breaks=NULL))
		dev.off()
		write.table(as.data.frame(object@sils),paste0('sils.',iterations,'.txt'),quote=F,sep='\t',col.names=NA)
		write.table(as.data.frame(object@numclust),paste0('numclust.',iterations,'.txt'),quote=F,sep='\t',col.names=NA)
	}
}
)




setGeneric("merge.subclusters.DE", function(object, min.de.genes = 25, effect.size=log(2), pval.cutoff=0.01, pcs.use=1:10,TPM.mat=NULL, Count.mat=NULL) standardGeneric("merge.subclusters.DE"))
setMethod("merge.subclusters.DE", "scDrop", 
          function(object,min.de.genes=25, effect.size=log(2),pval.cutoff=0.01, pcs.use=1:10, TPM.mat=NULL, Count.mat=NULL) {
          
            genes.use = rownames(object@data)
            clust.test = as.numeric(levels(object@subgroup))
           # if (is.null(tag)){
           #   filename = "CLUSTER_PAIRWISE_MARKERS.txt"
           # } else {
           #   filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
           # }
           # zz = file(filename,open="wt")
            
            
            num.clust=length(clust.test) 
            print(paste0("Starting with ", num.clust, " clusters"))
            
            pass.thresh=1e6*data.frame(diag(length(levels(object@subgroup)))); 
            
            for (i in setdiff(as.numeric(levels(object@subgroup)), clust.test)){
              pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
            } 
            
            dist.clust = pass.thresh
            
            #Find the number of differentially expressed genes between every pair of clusters
            for(k in 1:num.clust) {
              i=clust.test[k]
              print(paste0("Testing Cluster ", i))
              for(m in ((k+1):num.clust)) {
                j=clust.test[m]
                #print(j)
                if (m>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=markers.binom.subclust(object,i,j,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  
                  marker.pass=subset(marker,pval<pval.cutoff)
                  #print(paste("Test b/w Clusters ",i,j, "-# DE = ",nrow(marker.pass)))
                  #print(head(subset(marker.pass, log.effect > 0),5))
                  #print(head(subset(marker.pass, log.effect < 0),5))
                  
                  num.de.genes = 2*min(nrow(subset(marker.pass, log.effect > 0)), nrow(subset(marker.pass, log.effect < 0)))
                  pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
              print(pass.thresh[i,])
              
            }
            
            colnames(pass.thresh) = levels(object@subgroup)
            rownames(pass.thresh) = levels(object@subgroup)
            
            write.table(pass.thresh, file=paste0("DE_genes_matrix_2.txt"), sep="\t", quote=FALSE)
            
            #iteratively merge clusters
            min.val = min(pass.thresh)
            min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            merge.ind=-1
            while(min.val <= min.de.genes) {
              merge.ind=merge.ind+1
              
              #In case of ties, merge clusters that are closest in PC space
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
              
              if (pass.thresh[test.1,test.2]<= min.de.genes) {
                object@subgroup[which(object@subgroup==test.2)]=test.1
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
                old.group.levels = as.numeric(levels(object@subgroup))
                old.group.levels = setdiff(old.group.levels, test.2)
                clust.test = setdiff(clust.test, test.2)
                
                object@subgroup = droplevels(object@subgroup)
                levels(object@subgroup) = c(1:length(levels(object@subgroup)))
                object@meta[,"clust"] = object@subgroup
                
                new.group.levels = as.numeric(levels(object@subgroup))
                names(new.group.levels) = as.character(old.group.levels)
                clust.test = new.group.levels[as.character(clust.test)]
                
                
                
                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(clust.test, test.1)){
                  print(i)
                  marker= markers.binom.subclust(object,test.1,i,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  marker.pass=subset(marker,pval<pval.cutoff)
                  pass.thresh[test.1,i]=2*min(nrow(subset(marker.pass, log.effect>0)),nrow(subset(marker.pass, log.effect<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  #pass.thresh[test.1,i]=nrow(marker.pass); 
                  pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:length(levels(object@subgroup))
              rownames(pass.thresh) = colnames(pass.thresh)
              
              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)
              
            }
            return(object)
          }
)


setGeneric("dot.plot.sub.count", function(object,features.use=NULL, group.use=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,family=NULL,min.perc=0,...) standardGeneric("dot.plot.sub.count"))
setMethod("dot.plot.sub.count", "scDrop", 
          function(object,features.use=NULL,group.use=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,family=NULL,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@count.data)]
            if (is.null(group.use)) group.use = levels(object@subgroup)
            if (is.null(group.names)) {
            	group.names = group.use
            }
            else {
            	group.names = paste0(group.names," [",group.use,"]")
            }
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@count.data)]
            
            
            for (i in group.use){
              cells.in.cluster = names(object@subgroup)[which(object@subgroup== i)]
              vec.exp = apply(object@count.data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names

            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
              
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
            print(head(ExpVal))
            print(head(PercVal))
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              print('Final Dataframe')
              print(head(ExpVal))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) +  theme(text=element_text(family=family))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              print('Final Dataframe')
              print(head(ExpVal))
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster)),useDingbats=FALSE) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) + theme(text=element_text(family=family))
              print(p)
              
              
            }
              
          }
)
