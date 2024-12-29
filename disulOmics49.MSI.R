#install.packages("ggplot2")
#install.packages("ggpubr")


library(plyr)
library(ggplot2)
library(ggpubr)

riskFile="risk.all.txt"    
cliFile="MSI.txt"          
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\49.MSI")     

risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), row.names(cli))
rt=cbind(risk[sameSample,,drop=F], cli[sameSample,,drop=F])
rt$MSI=factor(rt$MSI, levels=c("MSS", "MSI-L", "MSI-H"))
rt$risk=factor(rt$risk, levels=c("low", "high"))

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,"MSI"]))]

rt1=rt[,c("MSI", "risk")]
df=as.data.frame(table(rt1))

df=ddply(df, .(risk), transform, percent = Freq/sum(Freq) * 100)

df=ddply(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")

p=ggplot(df, aes(x = factor(risk), y = percent, fill = MSI)) +
	   geom_bar(position = position_stack(), stat = "identity", width = .7) +
	   scale_fill_manual(values=bioCol)+
	   xlab("riskScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title="MSI"))+
	   geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
	   #coord_flip()+
	   theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()

rt2=rt[,c("MSI", "riskScore")]
type=levels(factor(rt2[,"MSI"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

boxplot=ggboxplot(rt2, x="MSI", y="riskScore", fill="MSI",
		          xlab="",
		          ylab="Risk score",
		          legend.title="MSI",
		          palette=bioCol
		          )+ 
	stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()
