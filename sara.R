sara=function(x,h,lamda){
  n=length(x)
  d=NULL
  local_max=NULL      #D(X)的局部最大值点
  for(i in (h+1):(n-h)){
    d[i-h]=(sum(x[(i-h):i])-sum(x[i:(h+i)]))/h
  }                   #计算D(x)
  d=abs(d)
  for(i in 1:length(d)){
    if(i<=h){
      if(all(d[i]>=d[1:i+h])){
        local_max=c(local_max,i+h)
      }
    }else if(i>=length(d)-h){
      if(all(d[i]>=d[(i-h):length(d)])){
        local_max=c(local_max,i+h)
      }
    }else{
      if(all(d[i]>=d[(i-h):(i+h)])){
        local_max=c(local_max,i+h)
      }
    }
  }
  return(local_max[d[local_max-h]>lamda])
}
y=c(rnorm(50),rnorm(50,1.2))
sara(y,10,1.5)
