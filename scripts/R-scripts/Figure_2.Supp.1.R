
require("sfsmisc")

require("ggplot2")
require("cowplot")
require("grid")

#antisense figure boxplot, GFP pilot study
GFP_pilot_expression=read.csv("/home/suu13/projects/antisense/GFP_design/GFP_selected/results/gfp.expression.4_hpi.txt",sep="\t")
GFP_pilot_values=read.csv("/home/suu13/projects/antisense/GFP_design/GFP_selected/results/selected.combined.csv",sep=" ")
GFP_pilot_RNALfold=read.csv("/home/suu13/projects/antisense/GFP_design/GFP_selected/results/GFP_contructs.13.designed.RNALfold.csv",sep="\t")

GFP_pilot_combined=merge(GFP_pilot_expression,GFP_pilot_values,"ID")
GFP_pilot_combined=cbind(GFP_pilot_combined,GFP=apply(GFP_pilot_combined[,3:6],1,mean)) #take the mean of replicates
GFP_pilot_combined=merge(GFP_pilot_combined,GFP_pilot_RNALfold,"ID")


GFP_boxplot_df=data.frame(hi.Avoid=as.vector(as.matrix(GFP_pilot_expression[grep("hiAvoid",GFP_pilot_expression$ID),][,2:5])),
                          low.Avoid=as.vector(as.matrix(GFP_pilot_expression[grep("loAvoid",GFP_pilot_expression$ID),][,2:5])),
                          low.Fold=as.vector(as.matrix(GFP_pilot_expression[grep("hiFolding",GFP_pilot_expression$ID),][,2:5])),
                          hi.Fold=as.vector(as.matrix(GFP_pilot_expression[grep("loFolding",GFP_pilot_expression$ID),][,2:5])),
                          hi.Codon=as.vector(as.matrix(GFP_pilot_expression[grep("hiCAI",GFP_pilot_expression$ID),][,2:5])),
                          low.Codon=as.vector(as.matrix(GFP_pilot_expression[grep("loCAI",GFP_pilot_expression$ID),][,2:5])),
                          Optimal=c(as.vector(as.matrix(GFP_pilot_expression[grep("optimi",GFP_pilot_expression$ID),][,2:5])),rep(NA,4)))





#kruskal wallis test
kruskal.test(GFP_boxplot_df)

#Colors=c("firebrick4","firebrick1","dodgerblue1","dodgerblue4","darkorange4","darkorange1","black")

Colors=c("firebrick4","darkorange4","dodgerblue4","firebrick1","darkorange1","dodgerblue1","black")

GFP_ggboxplot_df=data.frame(Type=rep(names(GFP_boxplot_df),each=8),
                            Fluorescence=as.vector(as.matrix(GFP_boxplot_df)),
                            Colors=rep(Colors,each=8))



plot.boxplot=ggplot(GFP_ggboxplot_df,aes(x=reorder(Type,Fluorescence,FUN=median),y=Fluorescence)) + geom_boxplot(aes(fill=Type)) + scale_fill_manual(name="",values=Colors) + xlab("") + theme(axis.title=element_text(size=25)) + theme_bw(base_size = 25) + theme(legend.text=element_text(size=16),legend.key.size=unit(3,"lines"),legend.title=element_text(size=16),axis.text.x=element_blank())
save_plot("Figure_2.Supp.1.pdf",plot.boxplot,base_height = 10 ,base_aspect_ratio = 1.3)


save.image(file="Figure_2.Supp.1.RData")


