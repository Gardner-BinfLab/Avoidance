

z_score=function(x,mean,sd){
  return ((x-mean)/sd)
}



proteomics_7_combined=read.table("/home/suu13/projects/antisense/protein_abundance_7/proteomics_7.csv",sep="\t",header = T)


extremes_7_protein=proteomics_7_combined[log10(proteomics_7_combined$Protein/proteomics_7_combined$mRNA) > 3.3 | log10(proteomics_7_combined$Protein/proteomics_7_combined$mRNA) < 1.495,]

extremes_7_protein=data.frame(extremes_7_protein,Group=log10(extremes_7_protein$Protein/extremes_7_protein$mRNA) <2,
                              Z_Score_Avoidance=z_score(extremes_7_protein$W_Binding_En.1.21,mean(extremes_7_protein$W_Binding_En.1.21),sd(extremes_7_protein$W_Binding_En.1.21)),
                              Z_Score_Folding=z_score(extremes_7_protein$Folding_En.1.37,mean(extremes_7_protein$Folding_En.1.37),sd(extremes_7_protein$Folding_En.1.37)),
                              Z_Score_CAI=z_score(extremes_7_protein$CAI,mean(extremes_7_protein$CAI),sd(extremes_7_protein$CAI)))

######################################################################

expPR<-read.table("/home/suu13/projects/antisense/protein_abundance_7/Laurent-E-coli-389-prot-mRNA-data.csv",sep="\t",header = T)
syy<-sort(log10( expPR$Protein/expPR$mRNA ))
breaks <- seq(from=floor(min(syy)),to=ceiling(max(syy)),length=15)
n <- 10
histL<-hist( syy[1:n]                             ,breaks=breaks, plot=F)
histM<-hist( syy[n:(length(syy)-n)]               ,breaks=breaks, plot=F)
histH<-hist( syy[(length(syy)-n+1):(length(syy))] ,breaks=breaks, plot=F)
histMAT<-matrix(c(histL$counts,histM$counts,histH$counts),nrow = length(breaks)-1, ncol = 3, dimnames=list(floor((10*histL$mids))/10,c("low exp.","","high exp.")))
histMAT=t(histMAT)

pdf("Figure_2.Supp.4.pdf",height = 12,width = 8)
cex.axis=1.5
par(mfrow=c(2,1),las=2,cex=1,mar=c(6,1,2,1),oma=c(2,5,0,0))#c(bottom, left, top, right)

greeny<-rgb(0.2,0.8,0.5,0.5)
barplot(histMAT,col=c(greeny,"azure","purple"),xlab="",ylab="Freq",main="") 
mtext("A",side=3,font=2,cex=1.8,las=1,adj = -0.15)

mtext(expression("log"[10]*"([Protein]/[mRNA])"),side = 1,F,line=3,cex=cex.axis,las=0)
mtext("Frequency",side = 2,outer=F,line=3,cex=cex.axis,las=0)
legend("topleft",c("low expression (n=10)","high expression (n=10)"),fill=c(greeny,"purple"),ncol=1,bty="n",cex=cex.axis)

boxplot(extremes_7_protein$Z_Score_Avoidance[extremes_7_protein$Group == TRUE],
        extremes_7_protein$Z_Score_CAI[extremes_7_protein$Group == TRUE],
        extremes_7_protein$Z_Score_Folding[extremes_7_protein$Group == FALSE],
        rnorm(n = 10000,mean = 0,sd = 1),
        extremes_7_protein$Z_Score_Folding[extremes_7_protein$Group == TRUE],
        extremes_7_protein$Z_Score_CAI[extremes_7_protein$Group == FALSE],
        extremes_7_protein$Z_Score_Avoidance[extremes_7_protein$Group == FALSE],
        names=c("Avoidance","CAI","Sec. Str.","Null","Sec. Str.","CAI","Avoidance"),
        xlab = "",ylab = "",axes = F,ylim=c(-2,2),col=c(greeny,greeny,"purple","white",greeny,"purple","purple"))
        mtext("B",side=3,font=2,cex=1.8,las=1,adj = -0.15)

lines(c(0.5,10),c(0,0))
axis(1,at=1:7,labels = c("Avoidance","Codon","Sec. Str.","Null","Sec. Str.","Codon","Avoidance"),cex.axis=cex.axis)
axis(2)
text(c(1,7),c(1,
            2),
     "*",cex=2.5)
mtext("Z score",side = 2,outer=F,line=3,cex=cex.axis,las=0)
legend("topleft",c("low expression (n=10)","high expression (n=10)"),bty="n",cex=cex.axis,fill=c(greeny,"purple"))
dev.off()



#determines order of box and whiskers:
sum(extremes_7_protein$Z_Score_Avoidance[extremes_7_protein$Group == TRUE])
sum(extremes_7_protein$Z_Score_CAI[extremes_7_protein$Group == TRUE])
sum(extremes_7_protein$Z_Score_Folding[extremes_7_protein$Group == FALSE])
sum(extremes_7_protein$Z_Score_Folding[extremes_7_protein$Group == TRUE])
sum(extremes_7_protein$Z_Score_CAI[extremes_7_protein$Group == FALSE])
sum(extremes_7_protein$Z_Score_Avoidance[extremes_7_protein$Group == FALSE])

save.image("Figure_2.Supp.4.RData")