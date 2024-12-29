#install.packages("igraph")
#install.packages("psych")
#install.packages("reshape2")
#install.packages("RColorBrewer")

library(igraph)
library(psych)
library(reshape2)
library(RColorBrewer)

GeneExpfile="disulGeneExp.txt"     
Genefile="gene.txt"               
Coxfile="uniCox.txt"              
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\20.network")    

gene.group <- read.table(Genefile,header=T,sep="\t")
gene.exp <- read.table(GeneExpfile,header=T,sep="\t",row.names=1)
gene.cox <- read.table(Coxfile,header=T,sep="\t")

colnames(gene.group) <- c('id','group')
genelist <- intersect(gene.group$id, gene.cox$id)
genelist <- intersect(genelist, rownames(gene.exp))
gene.group <- gene.group[match(genelist,gene.group$id),]
gene.group <- gene.group[order(gene.group$group),]
gene.exp <- gene.exp[match(gene.group$id,rownames(gene.exp)),]
gene.cox <- gene.cox[match(gene.group$id,gene.cox$id),]

gene.cor <- corr.test(t(gene.exp))
gene.cor.cor <- gene.cor$r
gene.cor.pvalue <- gene.cor$p
gene.cor.cor[upper.tri(gene.cor.cor)] = NA
gene.cor.pvalue[upper.tri(gene.cor.pvalue)] = NA
gene.cor.cor.melt <- melt(gene.cor.cor)   #gene1 \t gene2 \t cor
gene.cor.pvalue.melt <- melt(gene.cor.pvalue)
gene.melt <- data.frame(from = gene.cor.cor.melt$Var2,to=gene.cor.cor.melt$Var1,cor=gene.cor.cor.melt$value,pvalue=gene.cor.pvalue.melt$value)
gene.melt <- gene.melt[gene.melt$from!=gene.melt$to&!is.na(gene.melt$pvalue),,drop=F]
gene.edge <- gene.melt[gene.melt$pvalue<0.0001,,drop=F]
gene.edge$color <- ifelse(gene.edge$cor>0,'pink','#6495ED')
gene.edge$weight <- abs(gene.edge$cor)*6

gene.node <- gene.group
group.color <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(gene.node$group)))
gene.node$color <- group.color[as.numeric(as.factor(gene.node$group))]
gene.node$shape <- "circle"
gene.node$frame <- ifelse(gene.cox$HR>1,'purple',"green")
gene.node$pvalue <- gene.cox$pvalue
#?ڵ??Ĵ?С
pvalue.breaks <- c(0,0.0001,0.001,0.01,0.05,1)
pvalue.size <- c(16,14,12,10,8)
cutpvalue <- cut(gene.node$pvalue,breaks=pvalue.breaks)
gene.node$size <- pvalue.size[as.numeric(cutpvalue)]
nodefile <- "network.node.txt"
edgefile <- "network.edge.txt"
write.table(gene.node, nodefile, sep="\t", col.names=T, row.names=F, quote=F)
write.table(gene.edge, edgefile, sep="\t", col.names=T, row.names=F, quote=F)


#????????ͼ
node = read.table(nodefile, header=T, sep="\t", comment.char="")
edge = read.table(edgefile, header=T, sep="\t", comment.char="")

g = graph.data.frame(edge,directed = FALSE)
node = node[match(names(components(g)$membership),node$id),]

if(!is.na(match('color',colnames(node)))) V(g)$color = node$color
if(!is.na(match('size',colnames(node)))) V(g)$size = node$size
if(!is.na(match('shape',colnames(node)))) V(g)$shape = node$shape
if(!is.na(match('frame',colnames(node)))) V(g)$frame = node$frame

pdf(file="network.pdf", width=10, height=8)
par(mar=c(0,0,0,0))
layout(matrix(c(1,1,4,2,3,4),nc=2),height=c(4,4,2),width=c(8,3))

coord = layout_in_circle(g)
degree.x = acos(coord[,1])
degree.y = asin(coord[,2])
degree.alpha = c()
for(i in 1:length(degree.x)){
	if(degree.y[i]<0) degree.alpha=c(degree.alpha,2*pi-degree.x[i]) else degree.alpha=c(degree.alpha,degree.x[i])
}
degree.cut.group = (0:8)/4*pi
degree.cut.group[1] = -0.0001
degree.cut = cut(degree.alpha,degree.cut.group)
degree.degree = c(-pi/4,-pi/4,-pi/2,-pi/2,pi/2,pi/2,pi/2,pi/4)
degree = degree.degree[as.numeric(degree.cut)]

values <- lapply(node$id,function(x)c(1,1))
V(g)$pie.color = lapply(1:nrow(node),function(x)c(node$color[x],node$frame[x]))
V(g)$frame = NA 

plot(g,layout=layout_in_circle,vertex.shape="pie",vertex.pie=values,
	vertex.label.cex=V(g)$lable.cex,edge.width = E(g)$weight,edge.arrow.size=0,
	vertex.label.color=V(g)$color,vertex.frame.color=V(g)$frame,edge.color=E(g)$color,
	vertex.label.cex=2.5,vertex.label.font=2.5,vertex.size=V(g)$size,edge.curved=0.4,
	vertex.color=V(g)$color,vertex.label.dist=1.35,vertex.label.degree=degree)
# label.degree : zero means to the right; and pi means to the left; up is -pi/2 and down is pi/2;  The default value is -pi/4
# label.dist If it is 0 then the label is centered on the vertex; If it is 1 then the label is displayed beside the vertex.

par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
groupinfo = unique(data.frame(group=node$group,color=node$color))
legend("left",legend=groupinfo$group,col=groupinfo$color,pch=16,bty="n",cex=3)
#???ƻ??????յ?ͼ??(??Щ?????Ǹ߷??յĻ???,??Щ?????ǵͷ??յĻ???)
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",ylab="",axes=F)
legend("left",legend=c('Risk factors','Favorable factors'),col=c('purple','green'),pch=16,bty="n",cex=2.5)
#????Ԥ??pvalue??ͼ??
par(mar=c(0,0,0,0))
plot(1,type="n",xlab="",axes=F,ylab="")
legend("top",legend=c('Postive correlation with P<0.0001','Negative correlation with P<0.0001'),lty=1,lwd=4,col=c('pink','#6495ED'),bty="n",cex=2.2)
legend('bottom',legend=c(0.0001,0.001,0.01,0.05,1),pch=16,pt.cex=c(1.6,1.4,1.2,1,0.8)*6,bty="n",ncol=5,cex=2.2,col="black",title="Cox test, pvalue")
dev.off()

