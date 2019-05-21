avgbar_lattice<-function(gene,sham,sni) {
  shamgene=sham[gene,]
  snigene=sni[gene,]
  value=c(shamgene,snigene)
  value=unlist(value)
  s=rep('Sham',69)
  r=rep('SNI',69)
  condition=c(s,r)###
  cluster=rep(names(shamgene),2)
  data=data.frame(cluster,condition,value)
  datamin=min(data$value)
  datamax=max(data$value)
  OL=c(1,2,3,4,8,9)
  OL2=c(OL,OL+69)
  VA=c(5,6,7,30,31,32,40,41,42,43,51,52,53,54)
  VA2=c(VA,VA+69)
  NE=c(10,11,21,22,23,24,25,26,33,34,35,36,37,38,39,46,47,57,58,59,60,61,62,63,64,65,66)
  NE2=c(NE,NE+69)
  IM=c(12,13,14,15,16,55,56)
  IM2=c(IM,IM+69)
  SC=c(17,18,19,20,48,49,50)
  SC2=c(SC,SC+69)
  AS=c(27,28,29,44,45)
  AS2=c(AS,AS+69)
  OP=c(67,68,69)
  OP2=c(OP,OP+69)
  olig=data[OL2,]
  vasc=data[VA2,]
  neur=data[NE2,]
  immu=data[IM2,]
  schw=data[SC2,]
  astr=data[AS2,]
  opc=data[OP2,]
  
  #trellis.par.set(list(axis.text))
  p1.2=barchart(value~cluster,data=olig,groups=condition,newpage=FALSE,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','magenta2'),
                ylab='Avg Expr',main='Oligodendrocytes', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','magenta2'),pch=15,cex=2)),par.settings=noMargins())
  p2.2=barchart(value~cluster,data=vasc,groups=condition,newpage=FALSE,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','firebrick3'),
                ylab='Avg Expr',main='Vascular cells', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','firebrick3'),pch=15,cex=2)),par.settings=noMargins())
  p3.2=barchart(value~cluster,data=neur,groups=condition,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','royalblue'),
                ylab='Avg Expr',main='Neurons', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','royalblue'),pch=15,cex=2)),par.settings=noMargins())
  p4.2=barchart(value~cluster,data=immu,groups=condition,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','orange3'),
                ylab='Avg Expr',main='Immune cells', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','orange3'),pch=15,cex=2)),par.settings=noMargins())
  p5.2=barchart(value~cluster,data=schw,groups=condition,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','olivedrab'),
                ylab='Avg Expr',main='Schwann cells', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','olivedrab'),pch=15,cex=2)),par.settings=noMargins())
  p6.2=barchart(value~cluster,data=astr,groups=condition,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','turquoise4'),
                ylab='Avg Expr',main='Astrocytes/Progenitors', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','turquoise4'),pch=15,cex=2)),par.settings=noMargins())
  p7.2=barchart(value~cluster,data=opc,groups=condition,
                scales=list(x=list(rot=45,cex=0.8),family='sans'),
                ylim=c(datamin,datamax),col=c('grey','springgreen4'),
                ylab='Avg Expr',main='Committed Oligo.', 
                key=list(text=list(as.character(unique(olig$condition)),
                                   family='sans'),space='right',
                         points=list(col=c('grey','springgreen4'),pch=15,cex=2)),par.settings=noMargins())
  hlay=rbind(c(NA,NA,1,1,1,1,1,NA,NA),
             c(NA,2,2,2,2,2,2,2,NA),
             c(3,3,3,3,3,3,3,3,3),
             c(NA,NA,4,4,4,4,4,NA,NA),
             c(NA,NA,5,5,5,5,5,NA,NA),
             c(NA,NA,6,6,6,6,6,NA,NA),
             c(NA,NA,7,7,7,7,7,NA,NA))
  lattice.options(
    layout.heights=list(bottom.padding=list(x=0), top.padding=list(x=0)),
    layout.widths=list(left.padding=list(x=0), right.padding=list(x=0))
  )
  grid.newpage()
  gl=list(as.grob(p1.2),as.grob(p2.2),as.grob(p3.2),as.grob(p4.2),as.grob(p5.2),as.grob(p6.2),as.grob(p7.2))
  grid.arrange(grobs=gl,layout_matrix=hlay,top=textGrob(paste0('Avg ',gene,' expression\n'),gp=gpar(fontface="bold",fontsize=24)))
}
