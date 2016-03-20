

require("fields")


thermus_protein_protein=read.csv("TableS5_Thermus_SSU_profile.csv",sep="\t")



plot_avoidance_regions = function() {
  
  par(oma=c(1,1,2,1),las=1) #c(bottom, left, top, right)
  
  plot.new()
  mtext("A",side=3,adj=-0.05,outer = F,font=2,cex=1.8,line=4)
  mtext("B",side=1,adj=-0.05,outer = F,font=2,cex=1.8,line=5)
  par(lwd=1,lty=1,fig=c(0,1,0.7,0.95),mar=c(0,1,1,1),new=T)#c(bottom, left, top, right)
  image(t(matrix((thermus_protein_protein$Avoidance.P.Vals..0.001. < 0.001),nrow=1,ncol=nrow(thermus_protein_protein))),axes=T,col=c("white","red3"),xaxt ="n",yaxt="n" )
  text(0,1.5,"RNA Avoidance Regions",new=T,xpd=NA,pos=4,cex=1.5,font=2)
  
  par(lwd=1,lty=1,fig=c(0,1,0.30,0.55),mar=c(1,1,0,1),new=T)#c(bottom, left, top, right)
  image(t(matrix((thermus_protein_protein$Inaccessible.regions),nrow=1,ncol=nrow(thermus_protein_protein))),axes=T,col=c("white","royalblue4"),xaxt ="n",yaxt="n" )
  text(0,1.5,"Inaccessible Regions",new=T,xpd=NA,pos=4,cex=1.5,font=2)
  axis(1,at=c(seq(0,1,by=100*(1/1519))),labels=NA,line=0.2,cex.axis=0.8,las=3)
  text(c(seq(0,1,by=100*(1/1519))),-3,labels=c(1,seq(100,1500,by=100)),xpd=NA,srt=45,cex=1.8)
}


pdf("Figure_3.pdf",width = 18,height = 3)
plot_avoidance_regions()

dev.off()
