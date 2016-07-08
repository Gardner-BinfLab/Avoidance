require("fields")


retsign=function(a) {
  if(a[5] == 0) {
    a[1:4]=0 
    return(a)}
  else {
    a[1:4]=1
    return(a)}
  }

#SS="(((((((..((((...........)))).(((((.......))))).........................(((((.......))))))))))))."
#SS_chars=substring(SS,seq(1,nchar(SS),1),seq(1,nchar(SS),1))
tRNA_binding_heatmap=as.matrix(read.csv("~/Dropbox/tmp/tRNA.profile.thermus.combined.csv",sep="\t",header=F))
SRP_binding_heatmap=as.matrix(read.csv("~/Dropbox/tmp/thermus.SRP.binding.csv",sep="\t",header=F))
RnaseP_binding_heatmap=as.matrix(read.csv("~/Dropbox/tmp/thermus.RnaseP.binding.csv",sep="\t",header=F))
tmRNA_binding_heatmap=as.matrix(read.csv("~/Dropbox/tmp/thermus.tmRNA.binding.csv",sep="\t",header=F))


#tRNA_binding_heatmap[is.na(tRNA_binding_heatmap)]=1 #remove NAs
tRNA_order=hclust(dist(t(tRNA_binding_heatmap)))$order

tRNA_binding_heatmap=tRNA_binding_heatmap[,tRNA_order] #order after clustering

tRNA_binding_heatmap=tRNA_binding_heatmap < 0.001


tRNA_binding_heatmap[is.na(tRNA_binding_heatmap)]=-1



pdf("Figure_3.Supp.3.pdf",height = 12,width = 18)
par(mar=c(2,5,5,2),fig=c(0,1,0.5,1)) #c(bottom, left, top, right)
image(apply(tRNA_binding_heatmap,2,retsign),col=c("grey","white","red3"),axes=T,xaxt ="n",yaxt="n")


mtext("A",side=3,adj=-0.02,outer = F,font=2,cex=1.8,line=1)
axis(1,seq(0,1,1/95)[c(1,25,50,75,90)],labels=c(1,25,50,75,90),lwd.ticks = 0.5,cex.axis=1.5)
mtext("tRNAs",side=3,outer = F,cex=1.8,line=0)




par(mar=c(2,5,1,2),fig=c(0,1,0.35,0.45),new=T) #c(bottom, left, top, right)
image(as.matrix(tmRNA_binding_heatmap[,2] < 0.001),col=c("white","red3"),axes=T,xaxt ="n",yaxt="n" )
axis(1,seq(0,1,1/348)[c(1,100,200,300)],labels=c(1,100,200,300),lwd.ticks = 0.5,cex.axis=1.5)


mtext("B",side=3,adj=-0.02,outer = F,font=2,cex=1.8,line=1)
mtext("tmRNA",side=3,outer = F,cex=1.8,line=0)


par(mar=c(2,5,1,2),fig=c(0,1,0.2,0.3),new=T) #c(bottom, left, top, right)
image(as.matrix(RnaseP_binding_heatmap[,2] < 0.001),col=c("white","red3"),axes=T,xaxt ="n",yaxt="n" )
axis(1,seq(0,1,1/369)[c(1,100,200,300)],labels=c(1,100,200,300),lwd.ticks = 0.5,cex.axis=1.5)

mtext("C",side=3,adj=-0.02,outer = F,font=2,cex=1.8,line=1)
mtext("RNase P RNA",side=3,outer = F,cex=1.8,line=0)




par(mar=c(2,5,1,2),fig=c(0,1,0.05,0.15),new=T) #c(bottom, left, top, right)
image(as.matrix(SRP_binding_heatmap[,2] < 0.001),col=c("white","red3"),axes=T,xaxt ="n",yaxt="n")
axis(1,seq(0,1,1/95)[c(1,25,50,75,90)],labels=c(1,25,50,75,90),lwd.ticks = 0.5,cex.axis=1.5)

mtext("D",side=3,adj=-0.02,outer = F,font=2,cex=1.8,line=1)
mtext("SRP RNA",side=3,outer = F,cex=1.8,line=0)





dev.off()


save.image("Figure_3.Supp.3.RData")




