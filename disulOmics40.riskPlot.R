#install.packages("pheatmap")


library(pheatmap)       
setwd("C:\\Users\\lexb\\Desktop\\disulOmics\\40.riskPlot")    

bioRiskPlot=function(inputFile=null, project=null){
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    
	rt=rt[order(rt$riskScore),]      

	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	line[line>10]=10
	pdf(file=paste0(project, ".riskScore.pdf"), width=7, height=4)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)",
		 ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk","Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=paste0(project, ".survStat.pdf"), width=7, height=4)
	plot(rt$futime, pch=19,ylim=c(0,18),
		 xlab="Patients (increasing risk socre)",
		 ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead","Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
	
	ann_colors=list()
	bioCol=c("#0088FF", "#FF5555")
	names(bioCol)=c("low", "high")
	ann_colors[["Risk"]]=bioCol

	rt1=rt[c(3:(ncol(rt)-2))]
	rt1=t(rt1)
	annotation=data.frame(Risk=rt[,ncol(rt)])
	rownames(annotation)=rownames(rt)
	pdf(file=paste0(project, ".heatmap.pdf"), width=7, height=4)
	pheatmap(rt1, 
		     annotation=annotation,
		     annotation_colors = ann_colors, 
		     cluster_cols = FALSE,
		     cluster_rows = FALSE,
		     show_colnames = F,
		     scale="row",
		     color = colorRampPalette(c(rep("blue",3.5), "white", rep("red",3.5)))(50),
		     fontsize_col=7,
		     fontsize=7,
		     fontsize_row=8)
	dev.off()
}

bioRiskPlot(inputFile="risk.train.txt", project="train")

bioRiskPlot(inputFile="risk.test.txt", project="test")

bioRiskPlot(inputFile="risk.all.txt", project="all")
