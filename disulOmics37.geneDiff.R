#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")


library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

expFile="disulGeneExp.txt"     
riskFile="risk.all.txt"        
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\37.geneDiff")     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

risk=read.table(riskFile, sep="\t", header=T, check.names=F, row.names=1)
sameSample=intersect(row.names(data),row.names(risk))
rt1=cbind(data[sameSample,,drop=F],risk[sameSample,"risk",drop=F])

rt1=melt(rt1,id.vars=c("risk"))
colnames(rt1)=c("risk","Gene","Expression")

group=levels(factor(rt1$risk))
rt1$risk=factor(rt1$risk, levels=c("low","high"))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill="risk",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Risk",
				  width=0.8,
				  palette = c("#0088FF", "#FF5555") )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=risk),
	method="wilcox.test",
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
	
pdf(file="genediff.pdf", width=7, height=5)
print(boxplot)
dev.off()
