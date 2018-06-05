lasso <- function(x,y,t){
  k <- ncol(x)
  require(quadprog)
  D <- 2*t(x)%*%x
  d <- 2*t(y)%*%x
  b_ls <- solve(t(x)%*%x)%*%t(x)%*%y
  i <- 1
  delta_i <- drop(sign(b_ls))
  A <- cbind(-delta_i)
  b_hat <- solve.QP(Dmat = D,dvec = d, Amat = A,bvec = -t)$solution
  #判断两个向量是否相等
  iden <- function(x,y){
    if(all(x==y)){return(TRUE)}else{return(FALSE)}
  }
  #判断一个向量是否与矩阵的某行（列）相等
  vec_in_matrix <- function(x,m,marg){
    a <- apply(m, marg, iden,y=x)
    if(sum(a)==0){
      return(FALSE)
    }else{
      return(TRUE)
    }
  }
  while(sum(abs(round(b_hat,digits = 3)))>t+0.01){
   i <- i+1
   delta_i <- -sign(b_hat)
   while(vec_in_matrix(delta_i,A,2)==TRUE & i<2^k){
         delta_i <- sign(rnorm(k))
   }
   A <- cbind(A,delta_i)
   b_hat <- solve.QP(Dmat = D,dvec = d, Amat = A,bvec = rep(-t,ncol(A)))$solution
   round(b_hat,digits = 2)
  }
  return(round(b_hat,digits = 3))
}

#lars包中自带的diabetes数据集，对比lasso()的解和lars()的解
library(lars)
library(ggplot2)
library(reshape2)
data(diabetes)
attach(diabetes)
w <- cbind(diabetes$x, diabetes$y, diabetes$x2)
x <- as.matrix(w[, 1:10])
y <- as.matrix(w[, 11])#响应变量
x2 <- as.matrix(w[, 12:75])
laa <- lars(x2, y) #lars函数默认方法为lasso
cva <- cv.lars(x2, y, K = 10, plot.it = TRUE)
best <- cva$index[which.min(cva$cv)]
coef <- coef.lars(laa, mode = "fraction", s = best)
#lars解
coef
#lasso解
lasso(x2,y,t=best*sum(abs(solve(t(x2)%*%x2)%*%t(x2)%*%y)))

#系数轨迹图
b<- seq(0,1,length.out =100)
coef1 <- matrix(nrow =100 ,ncol = ncol(x2))
for (i in 1:100) {
  coef1[i,] <- lasso(x2,y,t=b[i]*sum(abs(solve(t(x2)%*%x2)%*%t(x2)%*%y)))
}
a <- as.data.frame(cbind(t(coef1),y= 1:ncol(x2)))
a1 <- melt(a,id.vars = 'y')
ggplot(a1,aes(variable,value,group=y))+geom_line(aes(col=y))
plot(laa)
