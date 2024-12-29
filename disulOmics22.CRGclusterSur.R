#install.packages("survival")
#install.packages("survminer")


#???ð?
library(survival)
library(survminer)

clusterFile="DRGcluster.txt"     #???͵Ľ????ļ?
cliFile="time.txt"               #?????????ļ?
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\22.DRGclusterSur")      #???ù???Ŀ¼

#??ȡ?????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#???ݺϲ?
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

length=length(levels(factor(rt$DRGcluster)))
diff=survdiff(Surv(futime, fustat) ~ DRGcluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ CRGcluster, data = rt)
#print(surv_median(fit))

#????????????
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="DRGcluster",
		           legend.labs=levels(factor(rt[,"DRGcluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 2,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)

#????ͼ??
pdf(file="survival.pdf", width=6.5, height=5.25, onefile=FALSE)
print(surPlot)
dev.off()