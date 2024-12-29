#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)          
expFile="merge.txt"      
geneFile="gene.txt"      #?????б??ļ?
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\18.cuproExp")     #???ù???Ŀ¼

#??ȡ?????ļ??????????ݽ??д???
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#??ȡ?????б??ļ?,??ȡͭ?????????ı???��
gene=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
geneExp=data[sameGene,]

#???=rbind(ID=colnames(geneExp),geneExp)
write.table(out, file="cuprdisulExp.txt", sep="\t", quote=F, col.names=F)


###