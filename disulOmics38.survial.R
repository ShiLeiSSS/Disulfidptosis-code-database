#install.packages("survival")
#install.packages("survminer")


library(survival)
library(survminer)
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\38.survival")       

bioSurvival=function(inputFile=null, outFile=null){
	#??ȡ?????ļ?
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	rt$risk=factor(rt$risk, levels=c("low", "high"))

	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Risk",
		           legend.labs=c("Low risk", "High risk"),
		           xlab="Time(years)",
		           break.time.by = 2,
		           palette=c("#0088FF", "#FF5555"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)
	pdf(file=outFile, width=5, height=4.5, onefile=FALSE)
	print(surPlot)
	dev.off()
}

bioSurvival(inputFile="risk.train.txt", outFile="surv.train.pdf")
bioSurvival(inputFile="risk.test.txt", outFile="surv.test.pdf")
bioSurvival(inputFile="risk.all.txt", outFile="surv.all.pdf")

