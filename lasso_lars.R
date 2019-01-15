lasso_lars = function(x,y){
  n=length(y)
  m=ncol(x)
  delta=diag(cov(x))
  x=scale(x)                  #标准化
  y=scale(y,T,F)              #中心化
  miu_hat =0                  #当前估计值
  beta_hat=diag(0,m)          #系数估计值
  c=t(x)%*%(y-miu_hat)        #相关系数向量
  c_max=max(abs(c))           #相关系数最大值
  a_index=which.max(abs(c))   #最先进入模型的变量
  s =sign(c[a_index])         #变量与残差向量相关系数的符号
  sigma_estimate=mean((y-x%*%solve(t(x)%*%x)%*%t(x)%*%y)^2)
  k = 1
  r=1
  x_current = s*x[,a_index]
  Cp=NULL                     #统计量Cp用于挑选最优参数
  BIC=NULL                    #BIC准则
  while(k<m+1){
    ##solution path求解路径
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
    ##计算eta_tilde，用于判断变量的系数是否经过零点
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
      error=y-miu_hat
      sse=sum(error^2)
      Cp[r]=sse/sigma_estimate-n+2*k
      BIC[r]=n*log(sse/n)+k*log(n)/2
      c=t(x)%*%error
      s=sign(c[a_index])
      k=k+1
      x_current=cbind(x_current,s[k]*x[,j])
    }else{
      miu_hat=miu_hat+eta_tilde*u_current
      beta_hat[a_index,r]=solve(t(x[,a_index])%*%x[,a_index])%*%t(x[,a_index])%*%miu_hat
      beta_hat=cbind(beta_hat,0,0)
      a_index=a_index[-j_tilde]
      c_max=c_max-eta_tilde*A_current
      error=y-miu_hat
      sse=sum(error^2)
      Cp[r]=sse/sigma_estimate-n+2*k
      BIC[r]=n*log(sse/n)+k*log(n)/2
      c=t(x)%*%error
      s=sign(c[a_index])
      k=k-1
      x_current=x_current[,-j_tilde]
    }
    r=r+1
  }
  beta_hat=apply(beta_hat,2,function(x){x/sqrt(delta)})
  sum_abs_beta=colSums(abs(beta_hat))
  lambda=sum_abs_beta/max(sum_abs_beta)
  return(list(beta_hat=beta_hat,lambda=lambda,Cp=Cp))
}

##交叉验证，默认10折交叉验证
cross_validation=function(x,y,k=10){
  n=nrow(x)
  r=sample(n,n)
  size=floor(n/k)
  lambda=seq(0,1,length=100)[-1]
  sse=matrix(nrow=k,ncol=99)
  for(i in 1:k){
    index=r[((i-1)*size+1):(i*size)]
    test_set_x=x[index,]
    test_set_y=y[index]
    train_set_x=x[-index,]
    train_set_y=y[-index]
    l=lasso_lars(train_set_x,train_set_y)
    train_beta=l$beta_hat
    train_lambda=l$lambda
    t=1
    for(j in 1:99){
      lambda_temp=lambda[j]
      while(lambda_temp>train_lambda[t]){
        t=t+1
      }
      if(t==1){
        beta_temp = train_beta[,t]*lambda_temp/train_lambda[t]
        sse[i,j]=mean((test_set_y - test_set_x%*%beta_temp)^2)
      }else{
        beta_temp = train_beta[,t-1]+
          (train_beta[,t]-train_beta[,t-1])*(lambda_temp-train_lambda[t-1])/(train_lambda[t]-train_lambda[t-1])
        sse[i,j]=mean((test_set_y - test_set_x%*%beta_temp)^2)
      }
    }
  }
  mse=colMeans(sse)  ##计算预测军方误差
  return(mse)
}

#例子
library(tidyverse)

##生成数据
t_j=c(20.1,35.1,58.1,80.2,130.2,150.1,162.1,176.1)
h_j=c(4,-2,-4,2.5,4.3,-3.1,2.1,-4.2)
n=200
t=1:n
f_real=NULL
for(i in 1:n){
  f_real[i]=h_j%*%(1+sign(t[i]-t_j))/2
}
f_t=tibble(t,f_t=f_real)
f_t_1=tibble(t,f_t=f_real+rnorm(n,sd=0.2))  
f_t_2=tibble(t,f_t=f_real+rnorm(n,sd=0.5))
f_t_3=tibble(t,f_t=f_real+rnorm(n,sd=0.7))
##时序图
ggplot(f_t,aes(t,f_t))+geom_line()+theme(axis.title=element_blank())
ggplot(f_t_1,aes(t,f_t))+geom_line()+theme(axis.title=element_blank())
ggplot(f_t_2,aes(t,f_t))+geom_line()+theme(axis.title=element_blank())
ggplot(f_t_3,aes(t,f_t))+geom_line()+theme(axis.title=element_blank())

##变点探测
###生成下三角矩阵X
X=diag(n)
X[lower.tri(X)]=1
###方差为0.2时
l_1=lasso_lars(X[,-1],f_t_1$f_t)
##选择Cp最小时的估计值
beta_h_1=l_1$beta_hat[,which.min(l_1$Cp)]
change_ponts_1=which(abs(beta_h_1)>.4)
ggplot(f_t_1,aes(t,f_t))+geom_line()+
  theme(axis.title=element_blank())+
  geom_vline(xintercept = change_ponts_1)

###方差为0.5时
l_2=lasso_lars(X[,-1],f_t_2$f_t)
##选择Cp最小时的估计值
beta_h_2=l_2$beta_hat[,which.min(l_2$Cp)]
change_ponts_2=which(abs(beta_h_2)>0.4)
ggplot(f_t_2,aes(t,f_t))+geom_line()+
  theme(axis.title=element_blank())+
  geom_vline(xintercept = change_ponts_2)

###方差为0.7时
l_3=lasso_lars(X[,-1],f_t_3$f_t)
##选择Cp最小时的估计值
beta_h_3=l_3$beta_hat[,which.min(l_3$Cp)]
change_ponts_3=which(abs(beta_h_3)>0.4)
ggplot(f_t_3,aes(t,f_t))+geom_line()+
  theme(axis.title=element_blank())+
  geom_vline(xintercept = change_ponts_3)
