#install.packages("survival")
#install.packages("survminer")
#install.packages("timeROC")


library(survival)
library(survminer)
library(timeROC)
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\39.ROC")      

bioROC=function(inputFile=null, rocFile=null){
	#??ȡ?????ļ?
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#??ȡROC???ߵĲ???
	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)
	
	pdf(file=rocFile, width=5, height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)

	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

bioROC(inputFile="risk.train.txt", rocFile="ROC.train.pdf")
bioROC(inputFile="risk.test.txt", rocFile="ROC.test.pdf")
bioROC(inputFile="risk.all.txt", rocFile="ROC.all.pdf")
