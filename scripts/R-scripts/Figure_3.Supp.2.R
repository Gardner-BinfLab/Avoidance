Avoidance_Table=read.csv("/home/suu13/projects/antisense/predictions/AvoidanceTableComplete.csv",sep="\t",header=T)
Avoidance_Table_2=read.csv("/home/suu13/projects/antisense/results/new-results.0.21/tablo.03.06.15.csv",sep="\t",header=T)



phylum_percent_model1=function(phylum,Avoidance_Table){
  tablo_av=table(Avoidance_Table[Avoidance_Table$Phylum == phylum,]$wxGreater.p.value.AVOIDANCEif.0.05. <= 0.05)
  return(round(100*as.numeric(tablo_av["TRUE"])/(as.numeric(tablo_av["FALSE"])+as.numeric(tablo_av["TRUE"]))))
}

phylum_n=function(phylum,Avoidance_Table){
  return(length(which(Avoidance_Table$Phylum == phylum)))
}

plot_barplot_bacteria=function(boxplot_of_pvalues) {
  
  boxplot(boxplot_of_pvalues,names=F,ylab="",axes=F,medcol="white",col=c("red","chartreuse4","yellow4","blue","purple","darkorange3","azure4"),ylim=c(-1,70))
  axis(1,at=1:7,labels=c("Act","Bac","Chl","Cya","Fir","Pro","Spi"),cex.axis=1.8,las=2)
  axis(2,las=3)
  abline(h=(-log10(0.05)),lty=4,lwd=1,col="black")
  
}

UTR5=sapply(as.vector(levels(Avoidance_Table$Phylum)),phylum_percent_model1,Avoidance_Table=Avoidance_Table)
UTR5=data.frame(Phylum=names(UTR5[!is.na(UTR5)]),Percent_Model1=as.numeric(UTR5[!is.na(UTR5)]),Count=as.numeric(sapply(names(UTR5[!is.na(UTR5)]),phylum_n,Avoidance_Table=Avoidance_Table)))
UTR5_20members=UTR5[UTR5$Count >= 20,]


CDS5=sapply(as.vector(levels(Avoidance_Table$Phylum)),phylum_percent_model1,Avoidance_Table=Avoidance_Table_2)
CDS5=data.frame(Phylum=names(CDS5[!is.na(CDS5)]),Percent_Model1=as.numeric(CDS5[!is.na(CDS5)]),Count=as.numeric(sapply(names(CDS5[!is.na(CDS5)]),phylum_n,Avoidance_Table=Avoidance_Table_2)))
CDS5_20members=CDS5[CDS5$Count >= 20,]





boxplot_of_pvalues_5UTR=list()

for(i in as.character(UTR5_20members$Phylum)){
  boxplot_of_pvalues_5UTR[[i]]=-log10(Avoidance_Table[Avoidance_Table$Phylum ==i,]$wxGreater.p.value.AVOIDANCEif.0.05.)
}

boxplot_of_pvalues_CDS=list()

for(i in as.character(CDS5_20members$Phylum)){
  boxplot_of_pvalues_CDS[[i]]=-log10(Avoidance_Table_2[Avoidance_Table_2$Phylum ==i,]$wxGreater.p.value.AVOIDANCEif.0.05.)
}

pdf("Figure_3.Supp.2.pdf",height = 9,width = 18)
layout(matrix(1:2,nrow=1,ncol=2))
par(oma=c(1,4,1,2),mar=c(5,3,2,1)) #c(bottom, left, top, right)

plot_barplot_bacteria(boxplot_of_pvalues_CDS)
mtext(expression('-log'[10]*italic(P)),cex=1.8,side=2,las=0,line=3)
mtext("A",side=3,adj=-0.1,outer = F,font=2,line=-0.5,cex=1.8,las=0)
mtext("5' CDS regions",side=3,outer = F,line=-0.5,cex=1.8,las=0)

plot_barplot_bacteria(boxplot_of_pvalues_5UTR)
mtext("B",side=3,adj=-0.1,outer = F,font=2,line=-0.5,cex=1.8,las=0)
mtext("5' UTR regions",side=3,outer = F,line=-0.5,cex=1.8,las=0)


dev.off()

save.image(file="Figure_3.Supp.2.RData")

