#install.packages("pheatmap")


library(pheatmap)        
expFile="uniSigExp.txt"            
geneCluFile="geneCluster.txt"      
crgCluFile="DRGcluster.txt"       
cliFile="clinical.txt"            
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\32.geneHeatmap")    

read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
crgClu=read.table(crgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#?Ï²Sample=intersect(row.names(exp), row.names(crgClu))
exp=exp[sameSample,,drop=F]
expData=cbind(exp, geneCluster=geneClu[sameSample,], CRGcluster=crgClu[sameSample,])
Project=gsub("(.*?)\\_.*", "\\1", rownames(expData))
rownames(expData)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expData))
expData=cbind(expData, Project)

#?Ï²read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expData), row.names(cli))
expData=expData[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expData, cli)

#??È=data[order(data$geneCluster),]
Type=data[,((ncol(data)-2-ncol(cli)):ncol(data))]
data=t(data[,1:(ncol(expData)-3)])

#???ol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CRGcol=bioCol[1:length(levels(factor(Type$CRGcluster)))]
names(CRGcol)=levels(factor(Type$CRGcluster))
ann_colors[["CRGcluster"]]=CRGcol
GENEcol=bioCol[1:length(levels(factor(Type$geneCluster)))]
names(GENEcol)=levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]]=GENEcol

#??Í"heatmap.pdf", height=6, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


##