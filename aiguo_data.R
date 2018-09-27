library(readxl)
aiguo <- read_xlsx('F:/±äµãÎÄÏ×/aiguo.xlsx',sheet = 1)
aiguo <- aiguo[,-c(1:3,17)]
opar<-par(no.readonly=T)
par(mfrow=c(3,4))
for(i in 1:12){
  a <- unlist(aiguo[i])
  plot(a,type='l')
}
par(opar)


aiguo <- as.data.frame(aiguo)
times <- aiguo[,13]
aiguo <- aiguo[-13]
rownames(aiguo) <- times
library(xts)
my_xts <- as.xts(aiguo)
my_xts
par(mfrow=c(3,4))
plot(my_xts[,1])
plot(my_xts[,2])
plot(my_xts[,3])
plot(my_xts[,4])
plot(my_xts[,5])
plot(my_xts[,6])
plot(my_xts[,7])
plot(my_xts[,8])
plot(my_xts[,9])
plot(my_xts[,10])
plot(my_xts[,11])
plot(my_xts[,12])
par(opar)
