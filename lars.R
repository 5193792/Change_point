my_lars <- function(x,y){
  m <- ncol(x)
  n <- length(y)
  delta <- diag(cov(x))
  x <- scale(x)
  y <- scale(y,T,F)
  sigma_estimate <- sum((y-x%*%solve(t(x)%*%x)%*%t(x)%*%y)^2)/(n-m-1)
  miu_hat <- 0                  #当前估计值
  beta_hat <- matrix(0, nrow = m, ncol = m)      #系数估计值
  beta_hat[,m] <- solve(t(x)%*%x)%*%t(x)%*%y
  c <- t(x)%*%(y-miu_hat)       #相关系数向量
  c_max <- max(abs(c))          #相关系数最大值
  a_index <- which.max(abs(c))
  s <- sign(c[a_index])
  k <- 1
  x_current <- s*x[,a_index]
  cp <- NULL
  while (k<m) {
    g_current <- t(x_current)%*%x_current
    A_current <- 1/drop(sqrt(rep(1,k)%*%solve(g_current)%*%rep(1,k)))
    w_current <- A_current*solve(g_current)%*%rep(1,k)
    u_current <- x_current%*%w_current
    a <- drop(t(x)%*%u_current)
    y_1 <- (c(c_max)-c)/(A_current-a)
    y_2 <- (c(c_max)+c)/(A_current+a)
    y_1[y_1<0] <- 1/0
    y_1[a_index] <- 1/0
    y_2[y_2<0] <- 1/0
    y_2[a_index] <- 1/0
    eta <- min(y_1,y_2)
    miu_hat <- miu_hat+eta*u_current
    beta_hat[a_index,k] <- solve(t(x[,a_index])%*%x[,a_index])%*%t(x[,a_index])%*%miu_hat
    j <- which(c(y_1,y_2)==eta)
    j <- ifelse(j>m,j-m,j)
    a_index <- c(a_index,j)
    c_max <- c_max-eta*A_current
    c <- t(x)%*%(y-miu_hat)
    s=sign(c[a_index])
    cp[k] <- sum((y-miu_hat)^2)/sigma_estimate-n+2*(k+1)
    k <- k+1
    x_current <- cbind(x_current,s[k]*x[,j])
  }
  cp[k] <- sum((y-miu_hat)^2)/sigma_estimate-n+2*(k+1)
  beta_hat <- apply(beta_hat,2,function(x){x/sqrt(delta)})
  sum_abs_beta <- colSums(abs(beta_hat))
  gamr <- sum_abs_beta/sum_abs_beta[10]
  return(list(beta_hat=beta_hat,gamr=gamr,cp=cp))
}


#例子
library(lars)
library(ggplot2)
library(reshape2)
data(diabetes)
attach(diabetes)
w <- cbind(diabetes$x, diabetes$y, diabetes$x2)
x <- as.matrix(w[, 1:10])
y <- as.matrix(w[, 11])#响应变量
x2 <- as.matrix(w[, 12:21])
laa <- lars(x2, y,type = 'lasso') 
plot(laa)
cva <- cv.lars(x2, y, K = 10, plot.it = TRUE,type = 'lar')
best <- cva$index[which.min(cva$cv)]           #交叉验证选择最优的解
coef <- coef.lars(laa, mode = 'step', s = best) 
#lars解
coef

#自己写的函数my_lars得到的解
coef1 <- my_lars(x2,y)[[1]]
#系数轨迹图
gamr <- my_lars(x2,y)[[2]]
a <- as.data.frame(cbind(coef1,y= 1:ncol(x2)))
a1 <- melt(a,id.vars = 'y')
ggplot(a1,aes(rep(gamr,each=10),value,group=y))+geom_line(aes(col=y))
which.min(my_lars(x2,y)[[3]])