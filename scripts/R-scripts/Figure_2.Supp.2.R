require(plotrix)


#antisense figure 4, GFP pilot study
GFP_df=read.table("~/projects/antisense/Avoidance_Git/data/Table3_GFP.csv",header = T)



kudla_raw=read.table("/home/suu13/projects/antisense/kudla/Kudla_window21/ncRNAvsKudlaGFP.0.21_AM946981.2_200Shuffled_Native/GFP.0.21.binding.folding.Energies.Folding37.csv",header=T)


plot_generic_scatterplot_GFP=function(x,y,Colors,xlim,ylim) {
  
  pch=21
  colbinding="black"
  cex=1.8
  bg="grey"
  bgbottom="firebrick1"
  col="darkred"
  
  
  #correlation=cor.test(log10(y),x,method="spearman")
  lmodel=lm(log10(y)~x)
    
  plot(x ,log10(y),xlim=xlim,ylim=ylim,
       ann = F,pch=pch,cex=cex,bg=as.character(Colors),axes=F)
  axis(2)
  axis(1)
  
  abline(lmodel,col=col)
  
  
}




pdf("Figure_2.Supp.2.pdf",width = 18,height = 10)
par(mfrow=c(2,3),oma=c(4,4,2,4),cex.axis=1.5,cex=1) #c(bottom, left, top, right) 

mcex=1.8

par(mar=c(2.5,2,2.5,2)) #c(bottom, left, top, right)




plot_generic_scatterplot_GFP(GFP_df$CAI,GFP_df$GFP,GFP_df$Color,c(0.40,0.85),c(4.2,4.7))
mtext("A",side=3,adj=-0.04,font=2,cex=1.8)
text(c(0.42,0.48),4.7,labels=c("hi","low"),xpd=T,cex=1.8)
points(c(0.42,0.42,0.42),c(4.65,4.60,4.55),pch=21,bg=c("firebrick4","dodgerblue4","darkorange4"),cex=1.8)
points(c(0.48,0.48,0.48),c(4.65,4.60,4.55),pch=21,bg=c("firebrick1","dodgerblue1","darkorange1"),cex=1.8)
points(c(0.48),c(4.50),pch=21,bg=c("black"),cex=1.8)
text(c(0.5,0.5,0.5,0.5),c(4.65,4.60,4.55,4.50),labels=c("Avoid","Fold","Codon","Optimal"),xpd=T,cex=1.8,pos=4)
legend("bottomright",expression("R"[s]*"=0.29"),bty="n",cex=1.8)


plot_generic_scatterplot_GFP(GFP_df$Folding_En,GFP_df$GFP,GFP_df$Color,c(-15,1),c(4.2,4.7))
mtext("B",side=3,adj=-0.04,font=2,cex=1.8)
legend("bottomright",expression("R"[s]*"=0.34"),bty="n",cex=1.8)



plot_generic_scatterplot_GFP(GFP_df$Binding_En,GFP_df$GFP,GFP_df$Color,c(-330,-120),c(4.2,4.7))
mtext("C",side=3,adj=-0.04,font=2,cex=1.8)
legend("bottomright",expression("R"[s]*"=0.56"),bty="n",cex=1.8)
mtext("Binding Energy (kcal/mol)",side = 1,outer=F,line=3,cex=mcex)







plot_generic_scatterplot_GFP(kudla_raw$CAI,kudla_raw$GFP,"grey",c(0.40,0.85),c(2,4))
mtext("D",side=3,adj=-0.04,font=2,cex=1.8)
mtext("CAI",side = 1,outer=F,line=3,cex=mcex)
legend("bottomright",expression("R"[s]*"=0.02"),bty="n",cex=1.8)



plot_generic_scatterplot_GFP(kudla_raw$Folding_En_0.37,kudla_raw$GFP,"grey",c(-15,1),c(2,4))
mtext("E",side=3,adj=-0.04,font=2,cex=1.8)
mtext("Folding Energy (kcal/mol)",side = 1,outer=F,line=3,cex=mcex)
legend("bottomright",expression("R"[s]*"=0.61"),bty="n",cex=1.8)

plot_generic_scatterplot_GFP(kudla_raw$W_Binding_En,kudla_raw$GFP,"grey",c(-330,-120),c(2,4))
mtext("F",side=3,adj=-0.04,font=2,cex=1.8)
mtext("Binding Energy (kcal/mol)",side = 1,outer=F,line=3,cex=mcex)
legend("bottomright",expression("R"[s]*"=0.65"),bty="n",cex=1.8)



mtext(expression("log"[10]*"(fluorescence)"),side=2,las=0,line=1,cex=mcex,outer=T)
dev.off()

























