library(data.table)
library(fgsea)
library(ggplot2)

#Load gmt pathway
msigdb.hallmark<-gmtPathways("path/to/your/gmt/file")


##Read your Differential genes analysis file
Mydata=read.csv("path/to/your/csv/file")


#replacing zero p-value with 10^-6
Mydata_pval=Mydata$p_val_adj
Mydata_pval=ifelse(Mydata_pval!=0,Mydata_pval, 10^-6)

#make your rank file
Mydata_rank=ifelse(Mydata$avg_log2FC >=0, -log10(Mydata_pval),log10(Mydata_pval))
names(Mydata_rank)=Mydata$feature
Mydata_rank
Mydata_rank=sort(Mydata_rank)
Mydata_rank

#Run FGSEA FUNCTION
fgsea.Mydata<-fgsea(pathways = msigdb.hallmark, stats = Mydata_rank)


#Plot
plotEnrichment(msigdb.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],Mydata_rank) +labs(title="TNFA_SIGNALING_VIA_NFKB")