#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

riskFile="risk.all.txt"      
RNAssFile="StemnessScores.tsv"     
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\50.Stemness")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

RNAss=read.table(RNAssFile, header=T, sep="\t", check.names=F, row.names=1)
RNAss=t(RNAss[1,,drop=F])
rownames(RNAss)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(RNAss))
RNAss=avereps(RNAss)

sameSample=intersect(row.names(risk), row.names(RNAss))
risk=risk[sameSample,"riskScore",drop=F]
RNAss=RNAss[sameSample,,drop=F]
data=cbind(RNAss, risk)

xlab="riskScore"
ylab="RNAss"
outFile="RNAss.cor.pdf"
x=as.numeric(data[,xlab])
x[x>quantile(x,0.99)]=quantile(x,0.99)
y=as.numeric(data[,ylab])
df1=as.data.frame(cbind(x,y))
p1=ggplot(df1, aes(x, y)) + 
		  xlab("Risk score") + ylab(ylab)+ ylim(0,0.7)+
		  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))

pdf(file=outFile, width=5.2, height=5)
print(p2)
dev.off()
