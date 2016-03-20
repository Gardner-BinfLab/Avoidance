require("fields")


require("plotrix")


#antisense figure 4, GFP pilot study
GFP_df=read.table("Table3_GFP.csv",header = T)
kudla_raw=read.table("GFP.supp.csv",header=T)





#heatmap settings
cor.matrix=matrix(c(0.56,0.34,0.29,0.73, #Our GFPs : *,*,*,NA #
                             0.65,0.62,0.02,0.53, #Kudla 2009 : *,*,-,NA
                             0.11,0.26,-0.2,0.27, #Goodman 2013 : *,*,-,*       
                             0.11,-0.05,0.45,0.69, #Laurent 2010 389 : *,-,*,*
                             0.11,-0.04,0.535,0.68, #Lu 2007 : *,-,*,*
                             0.161,0.11,0.23,0.74, #Kwon 2014 : *,*,*,*
                             0.17,0.07,0.17,0.69 #Laurent 2010 1148 : *,-,-,*
),byrow=T,ncol=4)

colnames(cor.matrix)=c("Avoid.","Sec. St.","Codon","mRNA ab.")
rownames(cor.matrix)=c("GFP reporter\n(n = 52(13))\n",
                       "GFP reporter\n(n = 154)\n",
                       "sfGFP-mCherry\n(n = 14234)\n",
                       "Microarray-MS\n(n = 389)\n",
                       "Microarray-AP(MS)\n(n = 3301)\n",
                       "Microarray-MS\n(n = 5479)\n",
                       "Microarray-MS\n(n = 1148)\n"   
                       
                       
)


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
       ann = F,pch=pch,cex=cex,bg=as.character(Colors),axes=F,ylim=c(4.2,4.7))
  axis(2,las=3)
  axis(1)
  
  abline(lmodel,col=col)
  
  
}










pdf("Figure_2.pdf",width = 18,height = 8)

par(oma=c(0,5,1,2),mar=c(1,1,1,1),cex.axis=1.5) # c(bottom, left, top, right) 
plot.new()
mtext("A",side=3,adj=-0.04,font=2,cex=1.8,line=0)

par(fig = c(0.05,0.9,0.70,0.90),new=T) #c(bottom, left, top, right) 

image(cor.matrix[,ncol(cor.matrix):1],col=two.colors(n=14,"cadetblue","firebrick","coral1"),axes=F)

x_labels_place=seq(0,1,by=(1/6))
y_labels_place=seq(0,1,1/3)

axis(2,y_labels_place,labels=colnames(cor.matrix[,ncol(cor.matrix):1]),las=2,cex.axis=1.5)
par(xpd=NA)

axis(3,x_labels_place,labels=FALSE,las=3,cex.axis=1.8)
text(x_labels_place,1.7,labels=rownames(cor.matrix),cex=1.1)
par(xpd=NA)

star.cex=2.5
strain.cex=1.1
strain.srt=0

strain.y.place=1.4
strain.x.correction=0

text(c(x_labels_place[1:7]),y_labels_place[4],"*",cex=star.cex) #Avoid
#text(c(x_labels_place[2:7]),y_labels_place[6],"*",cex=star.cex) #Avoid.Norm
text(c(x_labels_place[c(1,2,3,6,7)]),y_labels_place[3],"*",cex=star.cex) #Secondary
#text(c(x_labels_place[c(2,3,6)]),y_labels_place[4],"*",cex=star.cex) #Secondary.Norm
text(c(x_labels_place[c(1,4,5,6,7)]),y_labels_place[2],"*",cex=star.cex) #Codon
#text(c(x_labels_place[c(4,5,6)]),y_labels_place[2],"*",cex=star.cex) #Codon.Norm
text(c(x_labels_place[1:7]),y_labels_place[1],"*",cex=star.cex) #mRNA

text(x_labels_place[7],strain.y.place,"P. aeruginosa",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[6],strain.y.place,"P. aeruginosa",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[5],strain.y.place,"E. coli",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[4],strain.y.place,"E. coli",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[3],strain.y.place,"E. coli",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[2],strain.y.place,"E. coli",cex=strain.cex,srt=strain.srt,font=4)
text(x_labels_place[1],strain.y.place,"E. coli",cex=strain.cex,srt=strain.srt,font=4)






text(1,-0.3,expression("*"*italic(P)*" < 0.05"),cex=1.5)

mtext("Correlation Coefficient",side=4,line=5,las=0,cex=1.5,adj=0.4)


par(mar=c(1,1,1,1),fig = c(0.3,1,0.65,0.95),xpd=F)

image.plot(cor.matrix,col=two.colors(n=12,"cadetblue","firebrick","coral1"),axes=F,legend.only = T,horizontal = F, axis.args = list(cex.axis = 1.5,las=1))


cex.mtext=1.8

par(mar=c(1,1,2,2),fig=c(0,0.3,0.15,0.70),new=T) #c(bottom, left, top, right)
plot_generic_scatterplot_GFP(GFP_df$CAI,GFP_df$GFP,GFP_df$Color,c(0.6,0.85))
mtext("B",side=3,adj=0,font=2,cex=1.8)
mtext("CAI",side = 1,outer=F,line=3,cex=cex.mtext)
mtext(expression("log"[10]*"(fluorescence)"),side=2,las=0,line=3,cex=cex.mtext,outer=F)
legend("bottomright",expression("R"[s]*"=0.29"),bty="n",cex=1.8)

par(fig=c(0.3,0.6,0.15,0.70),new=T)
plot_generic_scatterplot_GFP(GFP_df$Folding_En,GFP_df$GFP,GFP_df$Color,c(-15,1))
mtext("C",side=3,adj=0,font=2,cex=1.8)
mtext("Folding Energy (kcal/mol)",side = 1,outer=F,line=3,cex=cex.mtext)
legend("bottomright",expression("R"[s]*"=0.34"),bty="n",cex=1.8)


par(fig=c(0.6,0.9,0.15,0.70),new=T)
plot_generic_scatterplot_GFP(GFP_df$Binding_En,GFP_df$GFP,GFP_df$Color,c(-350,-100))
mtext("D",side=3,adj=0,font=2,cex=1.8)
mtext("Binding Energy (kcal/mol)",side = 1,outer=F,line=3,cex=cex.mtext)
legend("bottomright",expression("R"[s]*"=0.56"),bty="n",cex=1.8)


par(fig=c(0.80,0.98,0.2,0.6),new=T,xpd=NA)
plot.new()
text(c(0.7,0.9),0.5,labels=c("hi","low"),cex=1.8)
points(c(0.7,0.7,0.7),c(0.4,0.3,0.2),pch=21,bg=c("firebrick4","dodgerblue4","darkorange4"),cex=1.8)
points(c(0.9,0.9,0.9),c(0.4,0.3,0.2),pch=21,bg=c("firebrick1","dodgerblue1","darkorange1"),cex=1.8)
text(c(1,1,1,1),c(0.4,0.3,0.2,0.1),labels=c("Avoid","Fold","Codon","Optimal"),xpd=NA,cex=1.8,pos=4)

points(c(0.9),c(0.1),pch=21,bg=c("black"),cex=1.8)

par(mar=c(1,1,2,1),fig=c(0,1,0,0.1),new=T) #c(bottom, left, top, right)
plot.new()


mtext("E",side=3,adj=-0.04,font=2,cex=1.8)



dev.off()




