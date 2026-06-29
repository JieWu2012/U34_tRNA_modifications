# Generate dicodon plot

library(ggplot2)
library(reshape2)




a <- read.table("/Users/jwu/Downloads/dicodon_paper_ncs2elp6_wt.txt")

b <- aggregate(a$V3,list(a$V1,a$V2,a$V7),mean)

b_count=aggregate(a$V5,list(a$V1,a$V2),mean)

c <- dcast(b,b$Group.1+b$Group.2~b$Group.3,value.var = "x")

colnames(b_count) = c("Codon_P", "Codon_A", "Count")

colnames(c) = c("Codon_P", "Codon_A", "NcS2elp6", "wt")

d <- merge(c, b_count, by=c('Codon_P','Codon_A'))

head(d)

d <- d[d$Codon_P!="TGA" & d$Codon_P!="TAG" & d$Codon_P!="TAA" &  d$Codon_A!="TGA" & d$Codon_A!="TAG" & d$Codon_A!="TAA",]

d$MtoW <- apply(d,1,function(x){return(log2(as.numeric(x[3])/as.numeric(x[4])))})

d$Quantile <- apply(d,1,function(x){ifelse(as.numeric(x[5])>quantile(d$Count,0.75),"4th quantile",ifelse(as.numeric(x[5])>quantile(d$Count,0.5),"3rd quantile",ifelse(as.numeric(x[5])>quantile(d$Count,0.25),"2nd quantile","1st quantile")))})

color2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))

d[,8] <- apply(d,1,function(x){if(as.numeric(x[6])>2){return(as.numeric(2))}else{return(as.numeric(x[6]))}})

d$Codon_P <- gsub("T", "U", d$Codon_P)

d$Codon_A <- gsub("T", "U", d$Codon_A)

codons <- c("AAA", "AAG", "AAC", "AAU", "ACA", "ACG", "ACC", "ACU", "AGA", "AGG", "AGC", "AGU", "AUA", "AUG", "AUC", "AUU", "CAA", "CAG", "CAC", "CAU", "CCA", "CCG", "CCC", "CCU", "CGA", "CGG", "CGC", "CGU", "CUA", "CUG", "CUC", "CUU", "GAA", "GAG", "GAC", "GAU", "GCA", "GCG", "GCC", "GCU", "GGA", "GGG", "GGC", "GGU", "GUA", "GUG", "GUC", "GUU", "UAC", "UAU", "UCA", "UCG", "UCC", "UCU", "UGG", "UGC", "UGU", "UUA", "UUG", "UUC", "UUU")

pdf("dicodon_paper_ncs2elp6_wt.pdf",width = 12,height = 10)

ggplot(d)+geom_point(aes(x=factor(Codon_P,levels = codons),y=factor(Codon_A,levels=codons),color=as.numeric(V8),size=factor(Quantile)))+scale_color_gradientn(colours =rev(color2(20)),name=expression("("*italic("ncs2"*Delta*"elp6"*Delta)*"/wt)"))+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+xlab("Codons in the P sites")+ylab("Codons in the A sites")+scale_size_discrete (range=c(2,4))

dev.off()