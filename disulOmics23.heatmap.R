#install.packages("pheatmap")


library(pheatmap)       #???ð?
expFile="disulGeneExp.txt"       #?????????ļ?
clusterFile="DRGcluster.txt"     #???͵Ľ????ļ?
cliFile="clinical.txt"           #?ٴ??????ļ?
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\23.heatmap")     #???ù???Ŀ¼

#??ȡ?????ļ?
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ??????ͷ???????
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)
Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
rownames(expCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expCluster))
expCluster=cbind(expCluster, Project)

#?ϲ??ٴ?????
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#??ȡ??ͼ????
data=data[order(data$CRGcluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

#??????ͼע?͵???ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$CRGcluster)))]
names(prgCluCol)=levels(factor(Type$CRGcluster))
ann_colors[["DRGcluster"]]=prgCluCol

#??ͼ???ӻ?
pdf("heatmap.pdf", width=7.5, height=5)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()