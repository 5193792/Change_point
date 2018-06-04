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
  #判断一个向量是否与矩阵的某行（列）相等
  vec_in_matrix <- function(x,m,marg){
    a <- apply(m, marg, identical,y=x)
    if(sum(a)==0){
      return(FALSE)
    }else{
      return(TRUE)
    }
  }
  while(sum(abs(b_hat))>t){
   i <- i+1
   delta_i <- -sign(b_hat)
   while(vec_in_matrix(delta_i,A,2)==TRUE & i<2^k){
         delta_i <- sign(rnorm(k))
   }
   A <- cbind(A,delta_i)
   b_hat <- solve.QP(Dmat = D,dvec = d, Amat = A,bvec = rep(-t,ncol(A)))$solution
  }
  return(b_hat)
}

data(diabetes)
attach(diabetes)
w <- cbind(diabetes$x, diabetes$y, diabetes$x2)
x <- as.matrix(w[, 1:10])
y <- as.matrix(w[, 11])#响应变量
x2 <- as.matrix(w[, 12:75])
laa <- lars(x2, y) #lars函数默认方法为lasso
plot(laa)
cva <- cv.lars(x2, y, K = 10, plot.it = TRUE)
best <- cva$index[which.min(cva$cv)]
coef <- coef.lars(laa, mode = "fraction", s = best)
sum(coef[coef != 0])
