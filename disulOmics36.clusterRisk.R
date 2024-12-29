#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#???ð?
library(limma)
library(ggpubr)
crgCluFile="DRGcluster.txt"        
geneCluFile="geneCluster.txt"     
scoreFile="risk.all.txt"           
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\36.clusterRisk")     

crgClu=read.table(crgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score$riskScore[score$riskScore>quantile(score$riskScore,0.99)]=quantile(score$riskScore,0.99)

twoCluster=cbind(crgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

data$DRGcluster=factor(data$DRGcluster, levels=levels(factor(data$DRGcluster)))
group=levels(factor(data$DRGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$CRGcluster)))]
	
boxplot=ggboxplot(data, x="DRGcluster", y="riskScore", color="DRGcluster",
			      xlab="DRGcluster",
			      ylab="Risk score",
			      legend.title="DRGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="DRGcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()

data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
boxplot=ggboxplot(data, x="geneCluster", y="riskScore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="Risk score",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
