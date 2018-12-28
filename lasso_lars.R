lasso_lars = function(x,y){
  n=length(y)
  m=ncol(x)
  delta=diag(cov(x))
  x=scale(x)
  y=scale(y,T,F)
  beta_ls=solve(t(x)%*%x)%*%t(x)%*%y   
  sigma_estimate=sum((y-x%*%beta_ls)^2)/(n-m-1)
  miu_hat =0                  #当前估计值
  beta_hat=diag(0,m)      #系数估计值
  c=t(x)%*%(y-miu_hat)       #相关系数向量
  c_max=max(abs(c))          #相关系数最大值
  a_index=which.max(abs(c))  #最先进入模型的变量
  s =sign(c[a_index]) 
  k = 1
  r=1
  x_current = s*x[,a_index]
  while(k<m+1){
    g_current = t(x_current)%*%x_current
    A_current = 1/drop(sqrt(rep(1,k)%*%solve(g_current)%*%rep(1,k)))
    w_current = A_current*solve(g_current)%*%rep(1,k)
    u_current = x_current%*%w_current
    a = drop(t(x)%*%u_current)
    ##计算下一个变量进入模型所需步长
    if(k<m){
      y_1 = (c(c_max)-c)/(A_current-a)
      y_2 = (c(c_max)+c)/(A_current+a)
      y_1[y_1<0] = 1/0
      y_1[a_index] = 1/0
      y_2[y_2<0] = 1/0
      y_2[a_index] = 1/0
      eta_hat = min(y_1,y_2)
      j = which(c(y_1,y_2)==eta_hat)
      j = ifelse(j>m,j-m,j)
    }else{
      eta_hat= max(abs(c))/A_current
    }
    ##计算eta_tilde，用于判断变量是否经过零点
    if(r>1){
      d=s*w_current
      gamma=-beta_hat[a_index,r-1]/d
      gamma[gamma<=0]=1/0
      eta_tilde=min(gamma)
      j_tilde=which.min(gamma)
    }else{
      eta_tilde=1/0
    }
    if(eta_tilde>=eta_hat){
      miu_hat=miu_hat+eta_hat*u_current
      beta_hat[a_index,r]=solve(t(x[,a_index])%*%x[,a_index])%*%t(x[,a_index])%*%miu_hat
      a_index=c(a_index,j)
      c_max=c_max-eta_hat*A_current
      c=t(x)%*%(y-miu_hat)
      s=sign(c[a_index])
      k=k+1
      x_current=cbind(x_current,s[k]*x[,j])
    }else{
      miu_hat=miu_hat+eta_tilde*u_current
      beta_hat[a_index,r]=solve(t(x[,a_index])%*%x[,a_index])%*%t(x[,a_index])%*%miu_hat
      beta_hat=cbind(beta_hat,0,0)
      a_index=a_index[-j_tilde]
      c_max=c_max-eta_tilde*A_current
      c=t(x)%*%(y-miu_hat)
      s=sign(c[a_index])
      k=k-1
      x_current=x_current[,-j_tilde]
    }
    r=r+1
  }
  beta_hat <- apply(beta_hat,2,function(x){x/sqrt(delta)})
  sum_abs_beta <- colSums(abs(beta_hat))
  gamr <- sum_abs_beta/max(sum_abs_beta)
  return(list(beta_hat=beta_hat,gamr=gamr))
}








#例子
library(lars)
library(ggplot2)
library(reshape2)
data(diabetes)
attach(diabetes)
w <- cbind(diabetes$x, diabetes$y, diabetes$x2)
y <- as.matrix(w[, 11])#响应变量
x <- as.matrix(w[, 12:21])
laa <- lars(x, y,type = 'lasso') 
plot(laa)

#自己写的函数my_lars得到的解
coef1 <- lasso(x,y)[[1]]
#系数轨迹图
gamr <- lasso(x,y)[[2]]
a <- as.data.frame(cbind(coef1,y= 1:ncol(x)))
a1 <- melt(a,id.vars = 'y')
ggplot(a1,aes(rep(gamr,each=10),value,group=y))+geom_line(aes(col=y))
