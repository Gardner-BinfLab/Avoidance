CDS.core.exp=read.csv("/home/suu13/projects/antisense/core_gene_expression/combined.CDS.expressions.txt",header = F)
Rfam.core.exp=read.csv("/home/suu13/projects/antisense/core_gene_expression/combined.rfam.expressions.txt",header = F)


pdf("Figure_4.pdf",height = 9,width = 9)
par(oma=c(1,4,1,2),mar=c(5,3,2,1),xpd=NA) #c(bottom, left, top, right)
boxplot(log10(Rfam.core.exp$V1[Rfam.core.exp$V1 !=0]),log10(CDS.core.exp$V1[CDS.core.exp$V1 !=0]),
        names=c("Core ncRNA genes","Core protein coding genes"),col=c("firebrick","cornflowerblue"))

mtext(expression('log'[10]*"(Mean read depth)"),side=2,outer = F,line=3,cex=1.8)
mtext("Gene type",side=1,outer = F,line=3,cex=1.8)
#mtext("Mean read depths of core genes \n for Archaea (n=5) and Bacteria (n=37)",side=3,outer = F,line=1,cex=1.8)

dev.off()

save.image(file="Figure_4.RData")
