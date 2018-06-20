rm(list=ls())
library(readxl)
library(data.table)
library(reshape2)
library(glmnet)
options(scipen=3)
##生成时间序列数据
means <- c(-0.8,1.2,0,-1,0.6,1.5,-1)
len <- c(100,60,260,100,160,120,200)
vars <- c(0.05,0.5)
##y用于存储生成的序列数据，每一列为一段时间序列
y <- matrix(nrow = sum(len),ncol = length(vars))
for (j in 1:length(vars)) {
  t <- NULL
  for (i in 1:length(means)) {
    t <- c(t,rnorm(len[i],mean = means[i],sd = vars[j]))
  }
  y[,j] <- t
}
par(mfrow=c(1,2))
plot(y[,1],type = 'l')
plot(y[,2],type = 'l')
#生成下三角矩阵，非零元全为1
lowtri <- function(k){ 
  a <- matrix(1,k,k)
  a[upper.tri(a)] <- 0
  return(a)
}
x <- lowtri(sum(len))[,2:sum(len)]
biandian_real <- cumsum(len)[-length(len)]+1
biandian_real
par(mfrow=c(1,1))
#第一条时间序列
fit2 <- cv.glmnet(x,y[,1])
plot(fit2)
coefs1 <-as.matrix(coef(fit1,s = 'lambda.min'))
rownames(coefs1) <- paste('b',1:1000,sep = '')
##系数值作图
plot(coef(fit2,s = 'lambda.min'),type = 'l')
as.matrix(round(coefs1[which(coefs1!=0),],3))
#筛选系数值大于阈值的系数的位置
biandian1 <- which(abs(coef(fit1,s = 'lambda.min'))>0.3)
#序列开头和结尾，不当做变点
biandian_estimate_1 <- biandian1[!(biandian1<=20 | biandian1>=980)]
biandian_estimate_1
#第二条时间序列
fit2 <- cv.glmnet(x,y[,2])
coefs2 <-as.matrix(coef(fit2,s = 'lambda.min'))
rownames(coefs2) <- paste('b',1:1000,sep = '')
as.matrix(round(coefs2[which(coefs2!=0),],3))
##系数值作图
plot(coef(fit2,s = 'lambda.min'),type = 'l')
#筛选系数值大于阈值的系数的位置
biandian2 <- which(abs(coef(fit2,s = 'lambda.min'))>0.3)
#序列开头和结尾，不当做变点
biandian_estimate_2 <- biandian2[!(biandian2<=20| biandian2>=980)]
##如果两个变点距离太近，就选择其中估计值较大的点当做变点
biandian_estimate_2 <- biandian_estimate_2[-(which(diff(biandian_estimate_2)<10)+1)]
biandian_estimate_2



##行业电力使用数据
fh <- read_xls("F:/FUHE/2008.xls",sheet = 1,na = "NULL")
dt <- as.data.table(fh)
setnames(dt,c("户名"),c("huming"))
##选取户名为上海坦达钢结构厂(乙电源)的数据
dt <- dt[huming=="上海坦达钢结构厂(乙电源)"]
df <- data.frame(dt)
df[,8:103] <- apply(df[,8:103],2,as.numeric)
shuju <- df[,8:103]
fuhe_y <- as.vector(t(shuju))
n <- length(fuhe_y)
##完整数据太长，只选取其中一半做变点探测
fuhe_y<- fuhe_y[1:(n/2)]
plot(fuhe_y,type='l')
x <- lowtri(n/2)[,2:(n/2)]
fuhe_fit <- cv.glmnet(x,fuhe_y)
coefs<-as.matrix(coef(fuhe_fit,s = 'lambda.min'))
rownames(coefs) <- paste('b',1:(n/2),sep = '')
as.matrix(round(coefs[which(coefs!=0),],3))
##系数值作图
plot(coef(fuhe_fit,s = 'lambda.min'),type = 'l')
#筛选系数值大于阈值的系数的位置
biandian_fh <- which(abs(coef(fuhe_fit,s = 'lambda.min'))>20)
#序列开头和结尾，不当做变点
biandian_estimate_fh <- biandian_fh[!(biandian_fh<=10| biandian_fh>=1430)]
##如果两个变点距离太近，就选择其中估计值较大的点当做变点
biandian_estimate_fh <- biandian_estimate_fh[-(which(diff(biandian_estimate_fh)<12)+1)]
biandian_estimate_fh
plot(fuhe_y,type='l')
abline(v  = biandian_estimate_fh,lwd=1,col="blue")

z <- rnorm(100,0,1)
z <- c(z,rnorm(100,2,1))
plot(z)
a <- diag(1,nrow=200)
for(i in 3:199){
  a[i,] <- c(rep(0,i-3),c(-1,-1,1,1),rep(0,199-i))
}
x <- solve(a)
fit1 <- cv.glmnet(x,z)
coefs1 <-as.matrix(coef(fit1,s = 'lambda.min'))
rownames(coefs1) <- paste('b',1:1000,sep = '')
##系数值作图
plot(coef(fit1,s = 'lambda.min'),type = 'l')



