#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("VennDiagram")


#???ð?
library(limma) 
library(VennDiagram)

expFile="merge.txt"         
cluFile="DRGcluster.txt"    
logFCfilter=0.585          
adj.P.Val.Filter=0.05        
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\27.clusterDiff")      

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(colnames(data), row.names(cluster))
data=data[,sameSample]
cluster=cluster[sameSample,]

geneList=list()
Type=as.vector(cluster)
design=model.matrix(~0+factor(Type))
colnames(design)=levels(factor(Type))
comp=combn(levels(factor(Type)), 2)

allDiffGenes=c()
for(i in 1:ncol(comp)){
	fit=lmFit(data, design)
	contrast=paste0(comp[2,i], "-", comp[1,i])
	#print(contrast)
	cont.matrix=makeContrasts(contrast, levels=design)
	fit2=contrasts.fit(fit, cont.matrix)
	fit2=eBayes(fit2)
	

	allDiff=topTable(fit2,adjust='fdr',number=200000)
	allDiffOut=rbind(id=colnames(allDiff),allDiff)
	write.table(allDiffOut, file=paste0(contrast, ".all.txt"), sep="\t", quote=F, col.names=F)
	
	#?????????ԵĲ???????
	diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
	diffSigOut=rbind(id=colnames(diffSig),diffSig)
	write.table(diffSigOut, file=paste0(contrast, ".diff.txt"), sep="\t", quote=F, col.names=F)
	geneList[[contrast]]=row.names(diffSig)
}

#????vennͼ
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)) )
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#???潻??????
interGenes=Reduce(intersect,geneList)
write.table(file="interGene.txt",interGenes,sep="\t",quote=F,col.names=F,row.names=F)

#???潻???????ı???��
interGeneExp=data[interGenes,]
interGeneExp=rbind(id=colnames(interGeneExp), interGeneExp)
write.table(interGeneExp, file="interGeneExp.txt", sep="\t", quote=F, col.names=F)


###