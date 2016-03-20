#barplot part
library(RColorBrewer)



phylum_n=function(phylum,Avoidance_Table){
  return(length(which(Avoidance_Table$Phylum == phylum)))
}

phylum_percent_model1=function(phylum,Avoidance_Table){
  tablo_av=table(factor(Avoidance_Table[Avoidance_Table$Phylum == phylum,]$wxGreater.p.value.AVOIDANCEif.0.05. <= 0.05,levels=c("TRUE","FALSE")))
  return(round(100*as.numeric(tablo_av["TRUE"])/(as.numeric(tablo_av["FALSE"])+as.numeric(tablo_av["TRUE"]))))
}

Avoidance_Table_2=read.csv("/home/suu13/projects/antisense/results/new-results.0.21/tablo.03.06.15.csv",sep="\t",header=T)
Avoidance_Table_Archaea=read.csv("/home/suu13/projects/antisense/new_results_for_0_21_archaea/archaea.table.txt",sep="\t",header=T)

PP_model1=sapply(as.vector(levels(Avoidance_Table_2$Phylum)),phylum_percent_model1,Avoidance_Table=Avoidance_Table_2)

PP_model_Archaea=sapply(as.vector(levels(Avoidance_Table_Archaea$Phylum)),phylum_percent_model1,Avoidance_Table=Avoidance_Table_Archaea)

Phylum_Percent_Table_2=data.frame(Phylum=names(PP_model1[!is.na(PP_model1)]),Percent_Model1=as.numeric(PP_model1[!is.na(PP_model1)]),Count=as.numeric(sapply(names(PP_model1[!is.na(PP_model1)]),phylum_n,Avoidance_Table=Avoidance_Table_2)))


Phylum_Percent_Table_Archaea=data.frame(Phylum=names(PP_model_Archaea[!is.na(PP_model_Archaea)]),Percent_Model1=as.numeric(PP_model_Archaea[!is.na(PP_model_Archaea)]),Count=as.numeric(sapply(names(PP_model_Archaea[!is.na(PP_model_Archaea)]),phylum_n,Avoidance_Table=Avoidance_Table_Archaea)))




plot_barplot=function(df,df_archaea) {
  
  par(oma=c(2,0,1,1),mar=c(0,15,1,2)) #c(bottom, left, top, right)
  
  par(fig=c(0,1,0.2,1))
  x=barplot(df$Percent_Model1,xlim=c(0,100),ylab="",axes = F,horiz = T,col="brown1") #x axisin yeri
  axis(2,at=x,labels=paste0(as.character(df$Phylum)," (n:",df$Count,")"),cex.axis=1,las=1)
  mtext("Bacteria",2,outer=F,line=13,cex=1.5)
  
  par(mar=c(0,15,0,2),fig=c(0,1,0.1,0.22),new=T)
  x=barplot(df_archaea$Percent_Model1,xlim=c(0,100),ylab="",axes = F,horiz = T,col="deepskyblue1") #x axisin yeri
  axis(2,at=x,labels=paste0(as.character(df_archaea$Phylum)," (n:",df_archaea$Count,")"),cex.axis=1,las=1)
  mtext("Archaea",2,outer=F,line=13,cex=1.5)
  
  axis(1,at=seq(0,100,by=20),labels=c("0","20","40","60","80","100"))
  mtext("Percentage",cex=1,side=1,las=1,line=2.3)
  

  
  
}

pdf("Figure_1.Supp.2.pdf",height = 10,width = 8)

plot_barplot(Phylum_Percent_Table_2[order(Phylum_Percent_Table_2$Percent_Model1,decreasing = T),],
             Phylum_Percent_Table_Archaea[order(Phylum_Percent_Table_Archaea$Percent_Model1,decreasing = T),])
dev.off()

