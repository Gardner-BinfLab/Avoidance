require("fields")


window_X=read.csv("/home/suu13/projects/antisense/kudla/Kudla_window21/window_21_distribution.extended.txt",header=F,sep="\t") #sort etmeyi unutma
names(window_X)=c("Type","StartPos","StopPos","CorrelationCoef","P-Value")


p_val_select=function(x) {if(x[1] >= 0) return(x[4]) else return(0)}

Correlation_DF=data.frame(Binding=window_X$CorrelationCoef[window_X$Type=="Binding"],StartPos=c(1:length(window_X$CorrelationCoef[window_X$Type=="Binding"])),Folding=window_X$CorrelationCoef[window_X$Type=="Folding"],Binding_logP=-log10(window_X$'P-Value'[window_X$Type=="Binding"]))


kudla21=read.table("/home/suu13/projects/antisense/kudla/Kudla_window21/ncRNAvsKudlaGFP.0.21_AM946981.2_200Shuffled_Native/GFP.0.21.binding.folding.txt",header = T)
kudla100=read.table("/home/suu13/projects/antisense/kudla/Kudla_window21/ncRNAvsKudlaGFP.100.121_AM946981.2_200Shuffled_Native/GFP.100.121.binding.folding.txt",header = T)
kudla300=read.table("/home/suu13/projects/antisense/kudla/Kudla_window21/ncRNAvsKudlaGFP.300.321_AM946981.2_200Shuffled_Native/GFP.300.321.binding.folding.txt",header = T)


plot_avoidance_regions = function() {
  #plot.new()
  image(t(matrix(apply(as.matrix(Correlation_DF),1,p_val_select),nrow=1,ncol=length(Correlation_DF$Binding_logP))),axes=T,col=two.colors(21,start="white",end="firebrick4",middle = "firebrick1"),xaxt="n",yaxt="n",bty="o")
  #mtext("a",side=3,adj=-0.1,outer = F,font=2,line=3)
  
  axis(1,at=c(seq(0,1,by=100*(1/700)),(21*(1/700))),labels=NA,line=0.2,cex.axis=0.8,las=3)
  text(c(seq(0,1,by=100*(1/700)),(21*(1/700))),-1.8,labels=c(1,100,"",300,rep("",4),21),xpd=NA,srt=45,cex=1.8)
  #text(c(seq(0,1,by=100*(1/700))),-2.5,labels=c(expression("R"[s]*"=0.65"),expression("R"[s]*"=0.18"),expression("R"[s]*"=0.03"),expression("R"[s]*"=0.1"),rep("",6)),xpd=NA,cex=1.8)
  rect(-0.2,-0.3,-0.001,0.3,xpd=NA)
  text(-0.1,0,labels=c("5'UTR"),xpd=NA,cex=1.8)
  rect(1.001,-0.3,1.2,0.3,xpd=NA)
  text(1.1,0,labels=c("3'UTR"),xpd=NA,cex=1.8)
  mtext("mRNA",side=3,line=3,outer=F,cex=1.8,font=2)
  mtext("ORF(CDS)",side=3,line=0.2,outer=F,cex=1.8)
  #text(-0.15,-2.5,labels=c("Correlation coef = "),xpd=NA,cex=1.8)
  

  
  
  
  
}



plot_generic_scatterplot_GFP=function(x,y,Colors,xlim) {
  
  pch=21
  colbinding="black"
  cex=1.8
  bg="grey"
  bgbottom="firebrick1"
  col="darkred"
  
  
  #correlation=cor.test(log10(y),x,method="spearman")
  lmodel=lm(log10(y)~x)
  
  plot(x ,log10(y),xlim=xlim,
       ann = F,pch=pch,cex=cex,bg=as.character(Colors),axes=F,ylim=c(2,4))
  axis(2,las=3)
  axis(1)
  
  abline(lmodel,col=col)
  
  
}

pdf("Figure_3.Supp.1.pdf",width = 18,height = 10)
par(oma=c(3,1,2,1),mar=c(5,1,5,1),las=1) #c(bottom, left, top, right)
par(lwd=1,lty=1,fig=c(0.15,0.8,0.65,1))#c(bottom, left, top, right)
plot_avoidance_regions()

par(lwd=1,lty=1,fig=c(0.15,0.97,0.65,1),new=T,cex.axis=1.5)#c(bottom, left, top, right)
mtext(expression('-log'[10]*italic(P)),side = 4,las=0,line=1.8,cex=1.8)
mtext("A",side=3,adj=-0.1,outer = F,font=2,cex=1.8,line=4)

image.plot(t(matrix(0:20,nrow=1,ncol=21)),
           col=two.colors(21,start="white",end="firebrick4",middle = "firebrick1"),axes=F,  legend.only = T, horizontal=F,
           axis.args = list(cex.axis = 1.5,lwd=1))


par(fig=c(0.04,0.34,0.05,0.6),mar=c(2,1,1,1),new=T)
plot_generic_scatterplot_GFP(kudla21$W_Binding_En,kudla21$GFP,"grey",c(-400,-100))
mtext(expression("log"[10]*"(fluorescence)"),side=2,cex=1.8,outer = F,las=3,line=3)

mtext("B",side=3,adj=-0.1,outer = F,font=2,cex=1.8,line=3)
mtext("1-21 region",side=3,outer = F,cex=1.8,line=1)
legend("bottomright",expression("R"[s]*"=0.65"),bty="n",cex=1.8)



par(fig=c(0.37,0.67,0.05,0.6),mar=c(2,1,1,1),new=T)
plot_generic_scatterplot_GFP(kudla100$W_Binding_En,kudla100$GFP,"grey",c(-400,-100))
mtext("100-121 region",side=3,outer = F,cex=1.8,line=1)
legend("bottomright",expression("R"[s]*"=0.13"),bty="n",cex=1.8)



par(fig=c(0.70,1,0.05,0.6),mar=c(2,1,1,1),new=T)
plot_generic_scatterplot_GFP(kudla300$W_Binding_En,kudla300$GFP,"grey",c(-400,-100))
mtext("300-321 region",side=3,outer = F,cex=1.8,line=1)
legend("bottomright",expression("R"[s]*"=-0.09"),bty="n",cex=1.8)



mtext("Binding Energy (kcal/mol)",side=1,cex=1.8,outer = T,las=1,line=1)


dev.off()

save.image(file="Figure_3.Supp.1.RData")






