#install.packages("survival")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")


#???ð?
library(limma)
library(survival)
library(ConsensusClusterPlus)

expFile="interGeneExp.txt"   
cliFile="time.txt"           
workDir="C:\\Users\\lexb\\Desktop\\disulOmics\\30.geneCluster"    
setwd(workDir)       

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
data2=data
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli,data)

sigGenes=c()
for(i in colnames(rt)[3:ncol(rt)]){
	cox=coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary=summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.05){ sigGenes=c(sigGenes,i) }
}
Tab=data2[,sigGenes]
outTab=cbind(id=row.names(outTab), outTab)
write.table(outTab, file="uniSigExp.txt", sep="\t", row.names=F, quote=F)
#???Tab=rt[,c("futime","fustat",sigGenes)]
outTab=cbind(id=row.names(outTab), outTab)
write.table(outTab, file="uniSigExpTime.txt", sep="\t", row.names=F, quote=F)

#???     #??????kdata2[,sigGenes])
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")


#??????clusterNum=3        =results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("geneCluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$geneCluster))
cluster$geneCluster=letter[match(cluster$geneCluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="geneCluster.txt", sep="\t", quote=F, col.names=F)


######