#!/usr/bin/Rscript
require("RColorBrewer")

#barplot part

phylum_n=function(phylum,Avoidance_Table){
  return(length(which(Avoidance_Table$Phylum == phylum)))
}

RfamGC.H=read.csv("/home/suu13/projects/antisense/core_fasta/CP002829.1/CP002829.1.Rfam.intrinsic.csv",header = F)$V1
CDSGC.H=read.csv("/home/suu13/projects/antisense/core_fasta/CP002829.1/CP002829.1.HMM_Core.0-21.intrinsic.csv",header = F)$V1

RfamGC.L=read.csv("/home/suu13/projects/antisense/core_fasta/CP002329.1/CP002329.1.Rfam.intrinsic.csv",header = F)$V1
CDSGC.L=read.csv("/home/suu13/projects/antisense/core_fasta/CP002329.1/CP002329.1.HMM_Core.0-21.intrinsic.csv",header = F)$V1



phylum_percent_model1=function(phylum,Avoidance_Table){
  tablo_av=table(Avoidance_Table[Avoidance_Table$Phylum == phylum,]$wxGreater.p.value.AVOIDANCEif.0.05. <= 0.05)
  return(round(100*as.numeric(tablo_av["TRUE"])/(as.numeric(tablo_av["FALSE"])+as.numeric(tablo_av["TRUE"]))))
}

Avoidance_Table_2=read.csv("/home/suu13/projects/antisense/results/new-results.0.21/tablo.03.06.15.csv",sep="\t",header=T)
Avoidance_Table_Archaea=read.csv("/home/suu13/projects/antisense/new_results_for_0_21_archaea/archaea.table.txt",sep="\t",header=T)


PP_model1=sapply(as.vector(levels(Avoidance_Table_2$Phylum)),phylum_percent_model1,Avoidance_Table=Avoidance_Table_2)

#PP_model1[!is.na(PP_model1)] #remove na values
Phylum_Percent_Table_2=data.frame(Phylum=names(PP_model1[!is.na(PP_model1)]),Percent_Model1=as.numeric(PP_model1[!is.na(PP_model1)]),Count=as.numeric(sapply(names(PP_model1[!is.na(PP_model1)]),phylum_n,Avoidance_Table=Avoidance_Table_2)))
Phylum_Percent_Table_2_20members=Phylum_Percent_Table_2[Phylum_Percent_Table_2$Count >= 20,]
Avoidance_Table_2_filtered=Avoidance_Table_2[Avoidance_Table_2$Phylum %in% as.character(Phylum_Percent_Table_2_20members$Phylum),]


boxplot_of_pvalues=list()

for(i in as.character(Phylum_Percent_Table_2_20members$Phylum)){
  boxplot_of_pvalues[[i]]=-log10(Avoidance_Table_2[Avoidance_Table_2$Phylum ==i,]$wxGreater.p.value.AVOIDANCEif.0.05.)
  #boxplot_of_pvalues[[i]]=Avoidance_Table_2[Avoidance_Table_2$Phylum ==i,]$wxGreater.p.value.AVOIDANCEif.0.05.
}



boxplot_of_pvalues[["Archaea"]]=-log10(read.csv("/home/suu13/projects/antisense/new_results_for_0_21_archaea/Archaea.pvalues.csv",header = F)$V1)

##############


Actinobacteria_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Actinobacteria.Native.txt")
Actinobacteria_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Actinobacteria.Shuffled.10Cols.txt")

Bacteroidetes_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Bacteroidetes.Native.txt")
Bacteroidetes_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Bacteroidetes.Shuffled.10Cols.txt")

Chlamydiae_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Chlamydiae.Native.txt")
Chlamydiae_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Chlamydiae.Shuffled.10Cols.txt")

Cyanobacteria_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Cyanobacteria.Native.txt")
Cyanobacteria_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Cyanobacteria.Shuffled.10Cols.txt")

Firmicutes_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Firmicutes.Native.txt")
Firmicutes_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Firmicutes.Shuffled.10Cols.txt")

Proteobacteria_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Proteobacteria.Native.txt")
Proteobacteria_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Proteobacteria.Shuffled.10Cols.txt")

Spirochaetes_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Spirochaetes.Native.txt")
Spirochaetes_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/Spirochaetes.Shuffled.10Cols.txt")

Archaea_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_archaea/Archaea.Native.txt")
Archaea_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_archaea/Archaea.Shuffled.10cols.txt")

#Negative.controls_Native=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/ncRNAvsmRNA.0-21.U00096.2_200Shuffled_Native_Core/U00096.2.Shuffled.1.txt")
#Negative.controls_Shuffled=read.table("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/ncRNAvsmRNA.0-21.U00096.2_200Shuffled_Native_Core/U00096.2.Shuffled.2.txt")



read_background_list=function(EMBL){
    Negative.controls_Native=read.table(paste0("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/ncRNAvsmRNA.0-21.",EMBL,"_200Shuffled_Native_Core/Negative.Controls.Shuffled5-45.txt"))
    Negative.controls_Shuffled=read.table(paste0("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/ncRNAvsmRNA.0-21.",EMBL,"_200Shuffled_Native_Core/Negative.Controls.Shuffled50-95.txt"))
    return(list(Native=Negative.controls_Native,Shuffled=Negative.controls_Shuffled,Color=rgb(0,0,0,alpha=0.3),lty=2,lwd=1))
  
}




phyla_list=list(Actinobacteria=list(Native=Actinobacteria_Native,Shuffled=Actinobacteria_Shuffled,Color="red",lty=1,lwd=4),
                Bacteroidetes=list(Native=Bacteroidetes_Native,Shuffled=Bacteroidetes_Shuffled,Color="chartreuse4",lty=1,lwd=4),
                Chlamydiae=list(Native=Chlamydiae_Native,Shuffled=Chlamydiae_Shuffled,Color="yellow4",lty=1,lwd=4),
                Cyanobacteria=list(Native=Cyanobacteria_Native,Shuffled=Cyanobacteria_Shuffled,Color="blue",lty=1,lwd=4),
                Firmicutes=list(Native=Firmicutes_Native,Shuffled=Firmicutes_Shuffled,Color="purple",lty=1,lwd=4),
                Proteobacteria=list(Native=Proteobacteria_Native,Shuffled=Proteobacteria_Shuffled,Color="darkorange3",lty=1,lwd=4),
                Spirochaetes=list(Native=Spirochaetes_Native,Shuffled=Spirochaetes_Shuffled,Color="skyblue4",lty=1,lwd=4),
                Archaea=list(Native=Archaea_Native,Shuffled=Archaea_Shuffled,Color="black",lty=1,lwd=4)
                )         


read_native_shuffled=function(EMBL){
  
  
  pattern_file=paste0("/home/suu13/projects/antisense/core_selected_submitted_finished_0_20_window/*",EMBL,"*/*/*sorted._R_.txt")
  Data_Matrix=read_files(Sys.glob(pattern_file))
  Native_Energies=as.vector(Data_Matrix[,1,])
  Shuffled_Energies=as.vector(Data_Matrix[,2:201,])
  return(list(Native=Native_Energies,Shuffled=Shuffled_Energies))
}


#read all calculation into a 3D matrix
read_files=function(R_Data_Files){
  i=1
  for (Single_File in R_Data_Files)
  {
    
    if (exists("DM")) {
      tmp_table_read_from_file=read.csv(Single_File,sep="\t",header=FALSE) #read.table yerine read.csv 
      
      DM[,,i]=as.matrix(tmp_table_read_from_file[,3:ncol(tmp_table_read_from_file)])
      
      
    } else {
      tmp_table_read_from_file=read.csv(Single_File,sep="\t",header=FALSE) #read.table yerine read.csv 
      DM=array(NA,c(nrow(tmp_table_read_from_file),201,length(R_Data_Files))) #define 3D matrix
      rownames(DM)=tmp_table_read_from_file[,1]
      DM[,,i]=as.matrix(tmp_table_read_from_file[,3:ncol(tmp_table_read_from_file)])
    }
    
    i=i+1;
    
  }
  
  return(DM)
}




heartbeat_paramaters=function(Phylum_List){
  
  NativeEnergyVal=Phylum_List$Native$V1
  ShuffledEnergyVal=Phylum_List$Shuffled$V1
  Color=Phylum_List$Color
  shuffled_counts=10 #default 200

  MinVal=-39 #def -30
  MaxVal=5 #def 5
  BreaksHist=seq(MinVal,MaxVal,length.out=150) #default 150
  
  NatHist=hist(+1*NativeEnergyVal,breaks=BreaksHist,plot=F)
  ShufHist=hist(+1*ShuffledEnergyVal,breaks=BreaksHist,plot=F)
  
  MaxPosition=max(NatHist$counts-(ShufHist$counts/shuffled_counts))
  MinPosition=min(NatHist$counts-(ShufHist$counts/shuffled_counts))
  par(new=T,fig=c(0.25,0.75,0,1),xpd=F)
  plot(NatHist$mids,(NatHist$density-ShufHist$density),ylab="",xlab="",xlim=c(-17,2),ylim=c(-0.010,0.015),axes=F,type="n")
  #lines(NatHist$mids,(NatHist$density-ShufHist$density),lwd=4,col=Color)
  lines(smooth.spline(NatHist$mids,(NatHist$density-ShufHist$density),spar = 0.45),lwd=Phylum_List$lwd,lty=Phylum_List$lty,col=Color)
  return(list(NatHist=NatHist,ShufHist=ShufHist,MaxPosition=MaxPosition,MinPosition=MinPosition))

}

plot_heartbeat=function(Lists_of_Phyla){
  
  par(mar=c(3,6,2,2),oma =c(2,4,1,0) ) #c(bottom, left, top, right)
  plot.new()
  lapply(lapply(as.character(Avoidance_Table_2_filtered$EMBL)[1:100],read_background_list), heartbeat_paramaters)
  
  lapply(Lists_of_Phyla,heartbeat_paramaters)
  axis(1)
  axis(2,line=3)
  mtext("Binding Energy (kcal/mol)",1,cex=1.8,line=3)
  mtext("Native Density - Shuffled Density",2,line=6,cex=1.8)
  mtext("B",side=3,adj=-0.2,outer = F,font=2,line=0.5,cex=1.8,las=0)
  
  abline(h=0,lwd=4,lty=2)
  par(fig=c(0.25,1,0,1),new=T,xpd=T)
  legend(-18,0.0175,legend=c(expression("Actinobacteria (n:163)"~italic(P)*" = 9.8x10"^-69),
                               expression("Bacteroidetes (n:60)"~italic(P)*" = 8.7x10"^-148),
                               expression("Chlamydiae (n:38)"~italic(P)*" = 1.4x10"^-193),
                               expression("Cyanobacteria (n:40)"~italic(P)*" = 3.8x10"^-11),
                               expression("Firmicutes (n:378)"~italic(P)*" = 0"),
                               expression("Proteobacteria (n:756)"~italic(P)*" = 0"),
                               expression("Spirochaetes (n:38)"~italic(P)*" = 1.6x10"^-98),
                               expression("Archaea (n:118)"~italic(P)*" = 4.2x10"^-177),
                               expression("Shuffled vs. Shuffled (n:100)")
                  )  
         ,pt.bg=c("red","chartreuse4","yellow4","blue","purple","darkorange3","skyblue4","black",rgb(0,0,0,alpha=0.3)                             )
         ,text.col=c("red","chartreuse4","yellow4","blue","purple","darkorange3","skyblue4","black",rgb(0,0,0,alpha=0.3)                             )
         ,pch=c(22),cex=1.5,
         y.intersp=1.5,
         bty="n")
  
  arrows(-13,-0.01,-17.3,-0.01,lwd=4,col="black")
  text(-15,-0.009,"More stable interactions",cex=1.5)
  
  
  arrows(-19.5,0.001,-19.5,0.012,lwd=4)
  text(-19.2,0.006,"Native interactions",cex=1.5,srt=90)
  
  arrows(-19.5,-0.001,-19.5,-0.012,lwd=4)
  text(-19.2,-0.006,"Shuffled interactions",cex=1.5,srt=90)
  
  par(mar=c(3,1,2,2), fig = c(0,0.23,0.50,1),new=T,xpd=F) #c(bottom, left, top, right)
  plot.new()
  plot_density_negative_controls("CP002633.1")
  par(mar=c(3,1,3,2),fig = c(0,0.23,0,0.50),new=T,xpd=F) #c(bottom, left, top, right)
  plot_barplot()
  
  par(mar=c(4,2,2,1), fig = c(0.75,1,0.42,1),new=T,xpd=F) #c(bottom, left, top, right)
  plot.new()
  plot_intrinsic_density()
  par(mar=c(3,2,3,1),fig = c(0.75,1,0,0.41),new=T,xpd=F) #c(bottom, left, top, right)
  plot_log_log_avoidance(Avoidance_Table_2_filtered,Avoidance_Table_Archaea)
  
  
  
}

plot_barplot=function() {
  
  boxplot(boxplot_of_pvalues,names=F,ylab="",axes=F,medcol="white",col=c("red","chartreuse4","yellow4","blue","purple","darkorange3","azure4","black"),ylim=c(-1,45))
  axis(1,at=1:8,labels=c("Act","Bac","Chl","Cya","Fir","Pro","Spi","Arc"),cex.axis=1.8,las=2)
  axis(2,las=3)
  mtext(expression('-log'[10]*italic(P)),cex=1.8,side=2,las=0,line=3)
  abline(h=(-log10(0.05)),lty=4,lwd=1,col="black")
  mtext("C",side=3,adj=-0.2,outer = F,font=2,line=-0.5,cex=1.8,las=0)
  
  
}

colfunction=function(phylum){
  col=c("red","chartreuse4","yellow4","blue","purple","darkorange3","azure4")
  names(col)=c("Actinobacteria","Bacteroidetes","Chlamydiae","Cyanobacteria","Firmicutes","Proteobacteria","Spirochaetes")
  
  return(as.character(col[phylum]))
  
}




plot_log_log_avoidance=function(df,df_archaea){
  
  ydata=c(-log10(df$Intrinsic),-log10(df_archaea$Intrinsic))
  xdata=c(-log10(df$wxGreater.p.value.AVOIDANCEif.0.05.),-log10(df_archaea$wxGreater.p.value.AVOIDANCEif.0.05.))
  coldata=c(colfunction(as.character(df$Phylum)),rep("black",nrow(df_archaea)))
  
  plot(xdata,ydata,xlab="",ylab="",pch=20,cex=0.5,col=coldata,yaxt="n")
  axis(2,las=3)
  abline(h=(-log10(0.05)),lty=4,lwd=2,col="black")
  abline(v=(-log10(0.05)),lty=4,lwd=2,col="black")
  
  mtext(expression('-log'[10]*italic(P)*'(Avoidance)'),1,cex=1.8,line=3,las=1)
  mtext(expression('-log'[10]*italic(P)*'(Intrinsic)'),2,line=3,cex=1.8,las=0)
  mtext("E",side=3,adj=-0.15,outer = F,font=2,line=1,cex=1.8,las=0)
  
  
}

plot_density_negative_controls=function(EMBL){
  
  Colors=brewer.pal(6,"Dark2")

  
  Native_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.0-21.RNAup.csv"),header = F)$V1
  Shuffled_Energies=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".Negative.Controls.Shuffled.txt"),header = F)$V1
  Controls_Phylum=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".phylum.alternative.csv"),header = F)$V1
  Controls_Downstream=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.100-121.RNAup.csv"),header = F)$V1
  Controls_Reverse_Comp=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".HMM_Core.revcomp.RNAup.csv"),header = F)$V1
  Controls_Intergenic=read.csv(paste0("/home/suu13/projects/antisense/two_window_test/",EMBL,"/",EMBL,".intergenic.csv"),header = F)$V1
  
  
  
  
  
  par(lwd=1,lty=1,las=1)
  lwd=4
  lty=2
  
  wxtest=wilcox.test(Native_Energies,Shuffled_Energies,alternative = "greater")
  
  xlim=c(-15,1)
  ylim=c(0,0.30)
  bw=0.25
  plot(density(Native_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[1],lwd=lwd,main="",axes=F)
  
  lines(density(Shuffled_Energies,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[2],lwd=lwd,lty=lty)
  
  lines(density(Controls_Phylum,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[3],lwd=lwd,lty=lty)
  lines(density(Controls_Downstream,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[4],lwd=lwd,lty=lty)
  lines(density(Controls_Reverse_Comp,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[5],lwd=lwd,lty=lty)
  lines(density(Controls_Intergenic,bw=bw),xlab="",ylab="",xlim=xlim,ylim=ylim,col=Colors[6],lwd=lwd,lty=lty)
  
  
  axis(1)
  axis(2,las=3)
  mtext("Density",side=2,line=3,cex=1.8,las=0)
  
  legend(-15,0.35,legend=c("Native","Shuffled","Different phylum","Downstream","Rev. complement","Intergenic"),
         col=Colors,cex=1.5,bty="n",lty=1,lwd=lwd,xpd=NA)
  
  
  mtext("Binding Energy (kcal/mol)",1,cex=1.8,line=3,las=1)
  mtext("A",side=3,adj=-0.2,outer = F,font=2,line=0.5,cex=1.8,las=0)
  

  
}


plot_intrinsic_density=function(){
  lwd=4
  
  plot(density(CDSGC.H,bw=0.05),xlim=c(0,1),ylim=c(0,10),main="",ylab="",xlab="",lwd=lwd,col="brown3",yaxt="n")
  axis(2,las=3)
  lines(density(RfamGC.H,bw=0.05),lty=2,lwd=lwd,col="brown3")
  lines(density(CDSGC.L,bw=0.05),lwd=lwd,col="deepskyblue3")
  lines(density(RfamGC.L,bw=0.05),lty=2,lwd=lwd,col="deepskyblue3")
  mtext("Density",2,line=3,cex=1.8,las=0)
  mtext("G+C Content",1,cex=1.8,line=3,las=1)
  legend("topleft",legend=c("ncRNA","mRNA"),lty=c(2,1),bty="n",cex=1.5,lwd=lwd)
  legend("topright",legend=c("High","Low"),fill = c("brown3","deepskyblue3"),bty="n",cex=1.5)
  mtext("D",side=3,adj=-0.15,outer = F,font=2,line=0.5,cex=1.8,las=0)
  
  
}
  

pdf("Figure_1.pdf",width = 18,height = 12)
par(cex.axis=1.5)
plot_heartbeat(phyla_list)
dev.off()

