require("np")
require("car")
require("RColorBrewer")


explainedR2=read.csv(file="/home/suu13/projects/antisense/Avoidance_Git/files/Supplementary_file_5.csv",header = T)

explainedR2=explainedR2[with(explainedR2,order(Type,decreasing = T)),]

exporder=c("GFP reporter (n = 52(13))","GFP reporter (n = 154)","sfGFP-mCherry (n = 14234)"
           ,"Microarray-LC-MS (n = 389)","Microarray-APEX(MS) (n = 3301)"
           ,"Microarray-LC-MS (n = 5479)","Microarray-LC-MS (n = 1148)") 

Colors=rev(brewer.pal(7,"Accent"))


pdf("Figure_2.Supp.3.pdf",width = 14,height = 10)

cex.mtext=1.5
par(oma=c(5,5,5,8),cex.axis=1.5)  #c(bottom, left, top, right) 
layout(matrix(c(1,1,2,2,3,4,5,6),ncol = 4,byrow = T))


par(mar=c(5,1,3,1))
m1=explainedR2[explainedR2$Model == "Protein abundance ~ Avoidance + Folding Energy + CAI + mRNA abundance",]
barplot(rev(m1[match(exporder,m1$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ Avoidance + Sec. St. + Codon + mRNA ab.",cex=cex.mtext)
grid(col="black")
mtext(expression("R"^"2"),side = 1,outer=F,line=3,cex=cex.mtext)


m2=explainedR2[explainedR2$Model == "Protein abundance ~ Folding Energy + CAI + mRNA abundance",]
barplot(rev(m2[match(exporder,m2$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ Sec. St. + Codon + mRNA ab.",cex=cex.mtext)
legend(0.5,6.5,bty = "n",legend = exporder,fill = brewer.pal(7,"Accent") ,cex=1.8,xpd = NA)
grid(col="black")
mtext(expression("R"^"2"),side = 1,outer=F,line=3,cex=cex.mtext)


m3=explainedR2[explainedR2$Model == "Protein abundance ~ Avoidance",]
barplot(rev(m3[match(exporder,m3$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ Avoid.",cex=cex.mtext)
grid(col="black")

m4=explainedR2[explainedR2$Model == "Protein abundance ~ Folding Energy",]
barplot(rev(m4[match(exporder,m4$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ Sec. St.",cex=cex.mtext)
grid(col="black")

m5=explainedR2[explainedR2$Model == "Protein abundance ~ CAI",]
barplot(rev(m5[match(exporder,m5$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ Codon",cex=cex.mtext)
grid(col="black")

m6=explainedR2[explainedR2$Model == "Protein abundance ~ mRNA abundance",]
barplot(rev(m6[match(exporder,m6$Type),]$R2),beside = T,horiz = T,xlim = c(0,0.9),col = Colors)
mtext("[Protein] ~ mRNA ab.",cex=cex.mtext)
grid(col="black")


mtext(expression("R"^"2"),side = 1,outer=T,line=0,cex=cex.mtext)


dev.off()

save.image("Figure_2.Supp.3.RData")