
library(GD) 
library(glue)
library(purrr)
library(data.table)


for (i in 1:100){

  print(i)
  fileName<-paste("fishNet_NemerowCSV",i,sep = '')
  filepath<-paste("E:\\syc\\soil\\fishNet_NemerowCSV5000\\",fileName,'.csv',sep = '')
  filepath
  Data<-read.csv(filepath,encoding='UTF-8')
  factors<-Data[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]#取列
  name<-names(factors)
  print(name)#列名，不同因子名
  Yname<-names(Data[18])
  data<-Data[,c(18,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]

  
  #最优离散
  discmethod <- c("equal", "natural", "quantile", "geometric", "sd")#离散方法
  discitv <- c(4:10)#划分的组别数

  data.Nemerow <- data
  data.continuous <- data.Nemerow[,c(1,11:16)]
  data.continuous
  #寻找最优离散方法。Nemerow：Y因子名，data.continuous：连续性因子
  odc1<-optidisc(Nemerow~., data = data.continuous,discmethod, discitv)

  #划分组别
  data.continuous <- do.call(cbind, lapply(1:6, function(x)
    data.frame(cut(data.continuous[, -1][, x], unique(odc1[[x]]$itv), include.lowest = TRUE))))
  data.continuous
  data.Nemerow[,11:16] <- data.continuous
  data.Nemerow

  #
  factors_detector<-gd(Nemerow~.,data=data.Nemerow)
  factors_detector
  factors_detectorCSV=paste("E:\\syc\\soil\\geoDetector\\Nemerow_csv5000New\\",Yname,"\\factors_detectorTest",i,'.csv',sep = '')
  factors_detectorCSV
  #write.csv(factors_detector[1],factors_detectorCSV)
  dev.off()


  risk_detector<-gdrisk(Nemerow~.,data=data.Nemerow)
  risk_detector
  risk_detectorCSV=paste("E:\\syc\\soil\\geoDetector\\Nemerow_csv5000New\\",Yname,"\\risk_detectorTest",i,'.csv',sep = '')
  risk_detectorCSV
  #write.csv(risk_detector[1],risk_detectorCSV)
  dev.off()

  riskmean_detector<-riskmean(Nemerow~.,data=data.Nemerow)
  riskmean_detector
  riskmean_detectorCSV=paste("E:\\syc\\soil\\geoDetector\\Nemerow_csv5000New\\",Yname,"\\riskmean_detectorTest",i,'.csv',sep = '')
  riskmean_detectorCSV
  #write.csv(riskmean_detector[1],riskmean_detectorCSV)
  dev.off()

  eco_detector<-gdeco(Nemerow~.,data=data.Nemerow)
  eco_detector
  eco_detectorCSV=paste("E:\\syc\\soil\\geoDetector\\Nemerow_csv5000New\\",Yname,"\\eco_detectorTest",i,'.csv',sep = '')
  eco_detectorCSV
  #write.csv(eco_detector[1],eco_detectorCSV)
  dev.off()

  interact_detector<-gdinteract(Nemerow~.,data=data.Nemerow)
  interact_detector
  interact_detectorCSV=paste("E:\\syc\\soil\\geoDetector\\Nemerow_csv5000New\\",Yname,"\\interact_detectorTest",i,'.csv',sep = '')
  interact_detectorCSV
  #write.csv(interact_detector[1],interact_detectorCSV)
  dev.off()

  
  print('done')

}
