#install.packages("survival")
#install.packages("caret")
#install.packages("glmnet")
#install.packages("survminer")
#install.packages("timeROC")


library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)

setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\34.model")      
rt=read.table("uniSigExpTime.txt", header=T, sep="\t", check.names=F, row.names=1)     
rt$futime[rt$futime<=1]=1
rt$futime=rt$futime/365

n=1000     
for(i in 1:n){
	
	inTrain<-createDataPartition(y=rt[,2], p=0.5, list=F)
	train<-rt[inTrain,]
	test<-rt[-inTrain,]
	trainOut=cbind(id=row.names(train),train)
	testOut=cbind(id=row.names(test),test)
	
	x=as.matrix(train[,c(3:ncol(train))])
	y=data.matrix(Surv(train$futime,train$fustat))
	fit <- glmnet(x, y, family = "cox", maxit = 1000)
	cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
	coef <- coef(fit, s = cvfit$lambda.min)
	index <- which(coef != 0)
	actCoef <- coef[index]
	lassoGene=row.names(coef)[index]
	lassoSigExp=train[,c("futime", "fustat", lassoGene)]
	lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
	geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
	if(nrow(geneCoef)<2){next}

	multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
	multiCox=step(multiCox,direction = "both")
	multiCoxSum=summary(multiCox)
	
	outMultiTab=data.frame()
	outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	outMultiTab=outMultiTab[,1:2]
	
	riskScore=predict(multiCox,type="risk",newdata=train)      
	coxGene=rownames(multiCoxSum$coefficients)
	coxGene=gsub("`","",coxGene)
	outCol=c("futime","fustat",coxGene)
	medianTrainRisk=median(riskScore)
	risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
	trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
		
	riskScoreTest=predict(multiCox,type="risk",newdata=test)   
	riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
	testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
	
	diff=survdiff(Surv(futime, fustat) ~risk,data = train)
	pValue=1-pchisq(diff$chisq, df=1)
	diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
	pValueTest=1-pchisq(diffTest$chisq, df=1)

	predictTime=3    
	roc=timeROC(T=train$futime, delta=train$fustat,
	            marker=riskScore, cause=1,
	            times=c(predictTime), ROC=TRUE)
	rocTest=timeROC(T=test$futime, delta=test$fustat,
	            marker=riskScoreTest, cause=1,
	            times=c(predictTime), ROC=TRUE)	
	
	if((pValue<0.01) & (roc$AUC[2]>0.65) & (pValueTest<0.05) & (rocTest$AUC[2]>0.6)){
		
		write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
	    
	    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
		pdf("lasso.lambda.pdf")
		plot(fit, xvar = "lambda", label = TRUE)
		dev.off()
		pdf("lasso.cvfit.pdf")
		plot(cvfit)
		abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
		dev.off()
	   
		write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)
		write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)
		
		allRiskOut=rbind(trainRiskOut, testRiskOut)
		write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
		break
	}
}
