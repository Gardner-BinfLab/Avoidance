require(ggplot2)
require(RColorBrewer)

Colors=brewer.pal(6,"Dark2")


plot_density_negative_controls=function(EMBL,Colors){
  par(lwd=1,cex.main = 1,lty=1,las=1)

  
  Native_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.RNAup.csv"),header = F)$V1
  Shuffled_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".Negative.Controls.Shuffled.txt"),header = F)$V1
  Controls_Phylum=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".phylum.alternative.csv"),header = F)$V1
  Controls_Downstream=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.100-121.RNAup.csv"),header = F)$V1
  Controls_Reverse_Comp=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.revcomp.RNAup.csv"),header = F)$V1
  Controls_Intergenic=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".intergenic.csv"),header = F)$V1
  
  
  
  xlim=c(-18,1)
  ylim=c(0,0.25)
  bw=0.5
  lty=5
  plot(density(Native_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[1],lwd=lwd,main="",axes=F)
  
  lines(density(Shuffled_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[2],lwd=lwd,lty=lty)

  lines(density(Controls_Phylum,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[3],lwd=lwd,lty=lty)
  lines(density(Controls_Downstream,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[4],lwd=lwd,lty=lty)
  lines(density(Controls_Reverse_Comp,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[5],lwd=lwd,lty=lty)
  lines(density(Controls_Intergenic,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[6],lwd=lwd,lty=lty)
  

  axis(1)
  axis(2,las=3)
  

  
  #mtext("Binding Energy (kcal/mol)",1,cex=2,line=2.4,las=1)
  #mtext("Density",2,line=3,cex=2,las=0)
  
  mtext("Binding Energy (kcal/mol)",1,cex=1.8,line=4,las=1)

  
}



plot_density_negative_controls_intermolecular=function(EMBL,Colors){
  par(lwd=1,cex.main = 1,lty=1,las=1)
  
  
  Native_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.RNAup.intermolecular.csv"),header = F)$V1
  Shuffled_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".Negative.Controls.Shuffled.intermolecular.txt"),header = F)$V1
  Controls_Phylum=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".phylum.alternative.intermolecular.csv"),header = F)$V1
  Controls_Downstream=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.100-121.RNAup.intermolecular.csv"),header = F)$V1
  Controls_Reverse_Comp=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.revcomp.RNAup.intermolecular.csv"),header = F)$V1
  Controls_Intergenic=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".intergenic.intermolecular.csv"),header = F)$V1
  
  
  
  
  
  
  xlim=c(-27,1)
  ylim=c(0,0.16)
  bw=0.5
  lty=5
  plot(density(Native_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[1],lwd=lwd,main="",axes=F)
  
  lines(density(Shuffled_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[2],lwd=lwd,lty=lty)
  
  lines(density(Controls_Phylum,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[3],lwd=lwd,lty=lty)
  lines(density(Controls_Downstream,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[4],lwd=lwd,lty=lty)
  lines(density(Controls_Reverse_Comp,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[5],lwd=lwd,lty=lty)
  lines(density(Controls_Intergenic,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[6],lwd=lwd,lty=lty)
  
  
  axis(1)
  axis(2,las=3)
  
  
  
  #mtext("Binding Energy (kcal/mol)",1,cex=2,line=2.4,las=1)
  #mtext("Density",2,line=3,cex=2,las=0)
  
  mtext("Hybridization Energy (kcal/mol)",1,cex=1.8,line=4,las=1)

}


plot_density_negative_controls_intramolecular=function(EMBL,Colors){
  par(lwd=1,cex.main = 1,lty=1,las=1)
  
  
  Native_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.RNAfold.csv"),header = T,sep="\t")$Folding_En
  Shuffled_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.200Shuffled.RNAfold.csv"),header = T,sep="\t")$Folding_En
  Controls_Phylum=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".phylum.alternative.RNAfold.csv"),header = T,sep="\t")$Folding_En
  Controls_Downstream=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.100-121.RNAfold.csv"),header = T,sep="\t")$Folding_En
  Controls_Reverse_Comp=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.revcomp.RNAfold.csv"),header = T,sep="\t")$Folding_En
  Controls_Intergenic=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".intergenic.RNAfold.csv"),header = T,sep="\t")$Folding_En
  
  
  

  
  xlim=c(-10,1)
  ylim=c(0,0.65)
  bw=0.5
  lty=5
  plot(density(Native_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[1],lwd=lwd,main="",axes=F)
  
  lines(density(Shuffled_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[2],lwd=lwd,lty=lty)
  
  lines(density(Controls_Phylum,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[3],lwd=lwd,lty=lty)
  lines(density(Controls_Downstream,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[4],lwd=lwd,lty=lty)
  lines(density(Controls_Reverse_Comp,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[5],lwd=lwd,lty=lty)
  lines(density(Controls_Intergenic,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[6],lwd=lwd,lty=lty)
  
  
  axis(1)
  axis(2,las=3)

  mtext("Folding Energy (kcal/mol)",1,cex=1.8,line=4,las=1)
  mtext("Density",2,line=3,cex=1.8,las=0)
  

}




pdf("Figure_1.Supp.1.pdf",width = 18,height = 6)
par(mar=c(2,2,2,2),oma=c(4,5,2,2),cex.axis=1.5) #c(bottom, left, top, right) 
lwd=3
lty=2

layout(matrix(c(1:3),ncol=3,byrow = T))

plot_density_negative_controls_intramolecular("CP002633.1",Colors)
mtext("A",side=3,adj=-0.1,outer = F,font=2,line=0.5,cex=1.8,las=0)

legend("topleft",legend=c("Native",
                          expression("Shuffled ("*italic(P)*" = 0.9)"),
                          expression("Different phylum ("*italic(P)*" = 5.28"^-4*")"),
                          expression("Downstream ("*italic(P)*" = 7.74"^-5*")"),
                          expression("Rev. complement ("*italic(P)*" = 0.8)"),
                          expression("Intergenic ("*italic(P)*" = 0.95)")),
       col=Colors,cex=2,bty="n",lty=1,lwd=lwd)

plot_density_negative_controls_intermolecular("CP002633.1",Colors)
mtext("B",side=3,adj=-0.1,outer = F,font=2,line=0.5,cex=1.8,las=0)

legend("topleft",legend=c("Native",
                          expression("Shuffled ("*italic(P)*" = 3.55"^-15*")"),
                          expression("Different phylum ("*italic(P)*" = 0 )"),
                          expression("Downstream ("*italic(P)*" = 5.39"^-143*")"),
                          expression("Rev. complement ("*italic(P)*" = 3.54"^-36*")"),
                          expression("Intergenic ("*italic(P)*" = 7.49"^-32*")")),
       col=Colors,cex=2,bty="n",lty=1,lwd=lwd)

plot_density_negative_controls("CP002633.1",Colors)
mtext("C",side=3,adj=-0.1,outer = F,font=2,line=0.5,cex=1.8,las=0)

legend("topleft",legend=c("Native",
                          expression("Shuffled ("*italic(P)*" = 7.69"^-52*")"),
                          expression("Different phylum ("*italic(P)*" = 0 )"),
                          expression("Downstream ("*italic(P)*" = 2.66"^-124*")"),
                          expression("Rev. complement ("*italic(P)*" = 6.51"^-57*")"),
                          expression("Intergenic ("*italic(P)*" = 6.16"^-93*")")),
       col=Colors,cex=2,bty="n",lty=1,lwd=lwd)


dev.off()

