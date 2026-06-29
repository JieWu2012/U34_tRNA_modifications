# R script for differential expression analysis of monsome and disome data.  

# For questions, please contact jie.wu@unibe.ch 

# Load required R package. 
library("ggplot2")
library("DESeq2")


# Set the working directory.
setwd("/Users/jwu/Modification_Project/")

# Read the monosome and disome data with raw read counts.  
mo=read.table("all_gene_monosome.txt", header = T, row.names = 1)
di=read.table("all_gene_disome.txt", header = T, row.names = 1)

# Rename the columns with monosome/disome and replicates. 
colnames(mo)=c("Ncs2elp6_1_mo","Ncs2elp6_2_mo","Ncs2elp6_3_mo","wt_1_mo","wt_2_mo","wt_3_mo","Ncs2elp6hel2_1_mo","Ncs2elp6hel2_2_mo","Hel2_1_mo","Hel2_2_mo")
colnames(di)=c("wt_1_di","wt_2_di","Hel2_1_di","Hel2_2_di","Ncs2elp6hel2_1_di","Ncs2elp6hel2_2_di","Ncs2elp6_1_di","Ncs2elp6_2_di")

# Merge monosome and disome counts and put the gene names in the rowname of data frame.  
rp=merge(mo,di,by=0,all=T)
rownames(rp)=rp[,1]
rp=rp[,-1]

# Extract raw counts in ncs2Delp6D and wt. 
countData=rp[,c(1:6,11:12,17:18)]

# Build sample information table "colData", similar to RPF vs mRNA in differential expression of translation efficiency:https://support.bioconductor.org/p/61509/.  .
colData <- data.frame (row.names = colnames(countData), condition = c("mutant","mutant","mutant","wt","wt","wt","wt","wt","mutant","mutant"),sampleType=c("mono","mono","mono","mono","mono","mono","di","di","di","di"))

# Construct DESeq2 object with interaction terms. 
dds <- DESeqDataSetFromMatrix (countData = countData,colData = colData, design = ~ sampleType+condition+sampleType:condition)

# Filter genes with sum of counts <= 5 in all samples. 
dds <- dds[ rowSums(counts(dds)) > 5, ]

# Set the reference level for both sampleType and condition. 
dds$sampleType <- relevel (dds$sampleType, ref = "mono")
dds$condition <- relevel (dds$condition, ref = "wt")

# Run DESeq 
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)

# Lists the coefficients. 
resultsNames(dds)

# Build a result table with threshold alpha = 0.05. 
res <- results(dds,name="sampleTypedi.conditionmutant",alpha = 0.05)

# Shrink the results
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)

# Defined genes with "Depleted" and "Enriched" disomes in mutant compared with wt. 
res1=data.frame(res2)
res1$Class=apply(res1,1,function(x){if(!is.na(x[5]) & as.numeric(x[5])<0.05){if(as.numeric(x[2]) < 0){return("Depleted")}else if(as.numeric(x[2])>0 ){return("Enriched")}}else{return("Nosig")}})

# Genes with |log2FC| > 4 are set to 4 or -4. Genes with -log10(p_adj) > 20 are set to 20 for better visualization.  
res1$log2FoldChange_1=apply(res1,1,function(x){if(!is.na(x[2]) & as.numeric(x[2])>4){return(4)}else if(as.numeric(x[2]) < -4){return(-4)}else{return(as.numeric(x[2]))}})
res1$padj_1=apply(res1,1,function(x){if(!is.na(x[5]) & -log10(as.numeric(x[5]))>20){return(20)}else{return(-log10(as.numeric(x[5])))}})

# Generate volcano plot for ncs2Delp6D vs wt. 
pdf("all_gene_di_vs_mo_ncs2elp6_vs_wt.pdf",width = 5.5,height = 4)
ggplot(res1,aes(x=log2FoldChange_1,y=padj_1,color=Class))+geom_abline(linetype="dashed",slope = 0,intercept = 0,color="gray")+geom_point(size=0.3)+theme_bw()+theme(panel.grid = element_blank(),text = element_text(size=12),axis.text=element_text(size=12),legend.text = element_text(size=12))+scale_color_manual(values = c("#4d9221","#c51b7d","gray"))+xlab("Log2 (fold change)")+ylab("-log10 (P value)")
dev.off()

# Write result table. 
res1 <- res1[order(res1$padj), ]
res1=data.frame(gene_name=rownames(res1),res1)
write.table(res1,"all_gene_di_vs_mo_ncs2elp6_vs_wt.txt",row.names=F,sep="\t",col.names=T,quote=F)


# Similar analysis is performed comparing ncs2Delp6Dhel2D vs ncs2Delp6D. 

countData=rp[,c(1,2,3,7,8,15:18)]
colData <- data.frame (row.names = colnames(countData), condition = c("wt","wt","wt","mutant","mutant","mutant","mutant","wt","wt"),sampleType=c("mono","mono","mono","mono","mono","di","di","di","di"))
dds <- DESeqDataSetFromMatrix (countData = countData,colData = colData, design = ~ sampleType+condition+sampleType:condition)
dds <- dds[ rowSums(counts(dds)) > 5, ]
dds$sampleType <- relevel (dds$sampleType, ref = "mono")
dds$condition <- relevel (dds$condition, ref = "wt")
dds <- DESeq(dds, test="LRT", reduced =  ~ sampleType+condition)
resultsNames(dds)
res <- results(dds,name="sampleTypedi.conditionmutant",alpha = 0.05)
res2 <- lfcShrink(dds, coef=4,type = "apeglm",res=res)
res1=data.frame(res2)
res1$Class=apply(res1,1,function(x){if(!is.na(x[5]) & as.numeric(x[5])<0.05){if(as.numeric(x[2]) < 0){return("Depleted")}else if(as.numeric(x[2])>0 ){return("Enriched")}}else{return("Nosig")}})
ggplot(res1,aes(x=log2FoldChange,y=-log10(padj),color=Class))+geom_abline(linetype="dashed",slope = 0,intercept = 0,color="gray")+geom_point(size=0.5)+theme_bw()+theme(panel.grid = element_blank(),legend.position = "non",text = element_text(size=12),axis.text=element_text(size=12),legend.text = element_text(size=12))+scale_color_manual(values = c("#4d9221","#c51b7d","gray"))+xlab("Mean of normalized counts")+ylab("log2(KO/wt)")
res1$log2FoldChange_1=apply(res1,1,function(x){if(!is.na(x[2]) & as.numeric(x[2])>4){return(4)}else if(as.numeric(x[2]) < -4){return(-4)}else{return(as.numeric(x[2]))}})
res1$padj_1=apply(res1,1,function(x){if(!is.na(x[5]) & -log10(as.numeric(x[5]))>10){return(10)}else{return(-log10(as.numeric(x[5])))}})
res1 <- res1[order(res1$padj), ]
res1=data.frame(gene_name=rownames(res1),res1)
write.table(res1,"all_gene_di_vs_mo_ncs2elp6hel2_vs_ncs2elp6.txt",row.names=F,sep="\t",col.names=T,quote=F)



# To make Figure 5C. 
n2e6_vs_wt=read.table("all_gene_di_vs_mo_ncs2elp6_vs_wt.txt",row.names = 1,header = T)
n2e6h2_vs_n2e6=read.table("all_gene_di_vs_mo_ncs2elp6hel2_vs_ncs2elp6.txt",row.names = 1,header = T)

# Merge two differential expression results.
all=merge(n2e6h2_vs_n2e6,n2e6_vs_wt,by=0,all=T)

# Remove empty values. 
all=na.omit(all)

# Generate volcano plot in Figure 5C. 
pdf("all_gene_di_vs_mo_ncs2elp6hel2_vs_ncs2elp6_coloredby_ncs2elp6_vs_wt.pdf",width = 5.5,height = 4)
ggplot(all,aes(x=log2FoldChange_1.x,y=padj_1.x,color=Class.y))+geom_abline(linetype="dashed",slope = 0,intercept = 0,color="gray")+geom_point(size=0.3)+theme_bw()+theme(panel.grid = element_blank(),text = element_text(size=12),axis.text=element_text(size=12),legend.text = element_text(size=12))+scale_color_manual(values = c("#4d9221","#c51b7d","gray"))+geom_point(data=all[all$Class.y!="Nosig",],aes(x=log2FoldChange_1.x,y=padj_1.x,color=Class.y),size=0.3)+xlab("Log2 (fold change)")+ylab("-log10 (P value)")
dev.off()


# RQC-target gene list. 
write.table(all[all$Class.x=="Depleted" & all$Class.y=="Enriched",],"RQC-target.txt",row.names=F,sep="\t",col.names=T,quote=F)

# RRT gene list. 
write.table(all[all$Class.x=="Enriched" & all$Class.y=="Depleted",],"RRT.txt",row.names=F,sep="\t",col.names=T,quote=F)