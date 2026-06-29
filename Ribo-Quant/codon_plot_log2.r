#!/usr/bin/Rscript
library(ggplot2)
options(warn=-1)
infilename=commandArgs()[5]
#print(infilename)
data=read.table(infilename)
foldname=commandArgs()[6]
colnames(data)=c("Codons","Samples","Class","Mean","Mean_nor","Sd_nor","Mean_nor_nor_wt","Sd_nor_nor_wt")

Quantile=ifelse(data$Mean>quantile(data$Mean,0.75),"4th quantile",ifelse(data$Mean>quantile(data$Mean,0.5),"3rd quantile",ifelse(data$Mean>quantile(data$Mean,0.25),"2nd quantile","1st quantile")))

data=cbind(data,Quantile)
Number_of_samples=length(unique(data$Samples))
Color=rainbow(Number_of_samples, s=.6, v=.9)[sample(1:Number_of_samples,Number_of_samples)]
data=data[data[,1]!="TAG" & data[,1]!="TAA" & data[,1]!="TGA",]

codons <- c("AAA", "AAG", "AAC", "AAT", "ACA", "ACG", "ACC", "ACT", "AGA", "AGG", "AGC", "AGT", "ATA", "ATG", "ATC", "ATT", "CAA", "CAG", "CAC", "CAT", "CCA", "CCG", "CCC", "CCT", "CGA", "CGG", "CGC", "CGT", "CTA", "CTG", "CTC", "CTT", "GAA", "GAG", "GAC", "GAT", "GCA", "GCG", "GCC", "GCT", "GGA", "GGG", "GGC", "GGT", "GTA", "GTG", "GTC", "GTT", "TAC", "TAT", "TCA", "TCG", "TCC", "TCT", "TGG", "TGC", "TGT", "TTA", "TTG", "TTC", "TTT")
#pdf("2.pdf")
pdf(paste(foldname,"/Codon_occupancy_log.pdf",sep=""),width=8,height=10)
ggplot(data,aes(x=factor(Codons,levels=rev(codons)),y=Mean_nor,color=factor(Samples),width=0, height=1.4))+geom_errorbar(aes(ymax=Mean_nor+Sd_nor,ymin=Mean_nor-Sd_nor),position = "dodge",colour="gray70", width=0.4, lwd=.5)+geom_point(aes(size=factor(Quantile)))+scale_size_discrete (range=c(2,4))+coord_flip()+scale_colour_manual(values= Color)+theme_bw()+ylab("Codon occupancy")+labs(size="Frequency Fraction",color="Sample")+xlab("Codons")
ggplot(data,aes(x=factor(Codons,levels=rev(codons)),y=log2(Mean_nor),color=factor(Samples),width=0, height=1.4))+geom_errorbar(aes(ymax=log2(Mean_nor+Sd_nor),ymin=log2(Mean_nor-Sd_nor)),position = "dodge",colour="gray70", width=0.4, lwd=.5)+geom_point(aes(size=factor(Quantile)))+scale_size_discrete (range=c(2,4))+coord_flip()+scale_colour_manual(values= Color)+theme_bw()+ylab("Codon occupancy")+labs(size="Frequency Fraction",color="Sample")+xlab("Codons")


data=data[data[,3]!="wt",]
#pdf(paste(foldname,"/Codon_occupancy_nor_wt.pdf",sep=""),width=8,height=10)
ggplot(data,aes(x=factor(Codons,levels=rev(codons)),y=Mean_nor_nor_wt,color=factor(Samples),width=0, height=1.4))+geom_errorbar(aes(ymax=Mean_nor_nor_wt+Sd_nor_nor_wt,ymin=Mean_nor_nor_wt-Sd_nor_nor_wt),position = "dodge",colour="gray70", width=0.4, lwd=.5)+geom_point(aes(size=factor(Quantile)))+scale_size_discrete (range=c(2,4))+coord_flip()+scale_colour_manual(values= Color)+theme_bw()+ylab("Codon occupancy")+labs(size="Frequency Fraction",color="Sample")+xlab("Codons")+theme(panel.grid.major.y=element_line(size=0.2,linetype = "dashed",color="gray"),text=element_text(family="Helvetica", size=12))+ggtitle("Codon Occupancy normalized to wt")
ggplot(data,aes(x=factor(Codons,levels=rev(codons)),y=log2(Mean_nor_nor_wt),color=factor(Samples),width=0, height=1.4))+geom_errorbar(aes(ymax=log2(Mean_nor_nor_wt+Sd_nor_nor_wt),ymin=log2(Mean_nor_nor_wt-Sd_nor_nor_wt)),position = "dodge",colour="gray70", width=0.4, lwd=.5)+geom_point(aes(size=factor(Quantile)))+scale_size_discrete (range=c(2,4))+coord_flip()+scale_colour_manual(values= Color)+theme_bw()+ylab("Codon occupancy")+labs(size="Frequency Fraction",color="Sample")+xlab("Codons")+theme(panel.grid.major.y=element_line(size=0.2,linetype = "dashed",color="gray"),text=element_text(family="Helvetica", size=12))+ggtitle("Codon Occupancy normalized to wt")
garbage <- dev.off()
options(warn=0)
