#################
library(rms)
setwd("D:\\Two\\Disulfidptosis\\NEW\\0.700\\41")
data.nomogram <- read.csv("data.nomogram.csv", header = T)

data.nomogram$gender <- factor(data.nomogram$gender, levels = c(0,1), labels = c("Female", "Male"))
data.nomogram$T <- factor(data.nomogram$T, levels = c(1,2,3,4), labels = c("T1", "T2","T3", "T4"))
data.nomogram$N <- factor(data.nomogram$N, levels = c(0,1,2,3), labels = c("N0", "N1", "N2","N3"))
data.nomogram$M <- factor(data.nomogram$M, levels = c(0,1), labels = c("M0","M1"))
data.nomogram$pathologic_stage <- factor(data.nomogram$pathologic_stage, levels = c(1,2,3,4), labels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
data.nomogram$Risk <- factor(data.nomogram$risk, levels = c(0,1), labels = c("Low", "High"))

ddist <- datadist(data.nomogram)
options(datadist = "ddist")
f <- cph(Surv(OS.time, OS) ~ gender+T+N+M+Risk+riskScore, data = data.nomogram,  x = T, y = T, surv = T)
surv <- Survival(f)

pdf(file = "nomogram.All.pdf", width = 18, height = 9) 
nom <- nomogram(f, fun = list(function(x) surv(365, x), function(x) surv(1095,x), function(x) surv(1825, x)), funlabel = c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"))
plot(nom, xfrac = .3, lplabel="Linear Predictor", col.grid = gray(c(0.8, 0.95)), cex.var = 1.6,cex.axis = 1.05,lwd=5, label.every = 1)
dev.off()

####
#绘制calibration curve进行验证
####
f1 <- cph(formula = Surv(OS.time, OS) ~  gender+T+N+M+Risk+riskScore, data = data.nomogram, x = T, y = T, surv = T, time.inc = 365) 
cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=50, B=1000) 
pdf("calibration_1y.All.pdf",width = 8,height = 8)
plot(cal1,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

f3 <- cph(formula = Surv(OS.time, OS) ~  gender+T+N+M+Risk+riskScore, data = data.nomogram, x = T, y = T, surv = T, time.inc = 1095) 
cal2 <- calibrate(f3, cmethod="KM", method="boot", u=730, m=50, B=1000) 
pdf("calibration_3y.All.pdf",width = 8,height = 8)
plot(cal2,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal2[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

f5 <- cph(formula = Surv(OS.time, OS) ~  gender+T+N+M+Risk+riskScore, data = data.nomogram, x = T, y = T, surv = T, time.inc = 1825) 
cal3 <- calibrate(f5, cmethod="KM", method="boot", u=1095, m=50, B=1000) 
pdf("calibration_5y.All.pdf",width = 8,height = 8)
plot(cal3,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"),#error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小
lines(cal3[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("#2166AC")) #连线的颜色
mtext("")
box(lwd = 1) #边框粗细
abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
) 
dev.off()

pdf("calibration_compare.All.pdf",width = 8,height = 8)
plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2, lwd = 2, lty = 0, errbar.col = c("#669933"), 
     xlim = c(0,1),ylim= c(0,1), col = c("#669933"),add = T) 
lines(cal2[,c('mean.predicted',"KM")], 
      type = 'b', lwd = 1, col = c("#669933")) 

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year", "3-year", "5-year"), #图例文字
       col =c("#2166AC", "669933", "#B2182B"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()

rm(list = ls())
################

