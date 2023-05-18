iband<-function(x,alpha=.05,q = c(0.1, 0.25, 0.5, 0.75, 0.9),  method='BH', SW=FALSE, plotit=FALSE,SEED=TRUE,nboot=500,grp=c(1:4),
                xlab='X'){
  #
  # 2 by 2 design.
  #
  # For variables x1, x2, x3 and x4,
  # This function compares the quantiles of the distributions
  # d1=x1-x2 and d2=x3-x4
  #
  # SW=TRUE: switch rows and columns
  #
  if(SEED)set.seed(2)
  if(is.matrix(x) || is.data.frame(x))x<-listm(x)
  if(length(x)!=4)stop('Should be exactly 4 groups')
  for(j in 1:length(x))x[[j]]=elimna(x[[j]])
  if(SW)x=x[c(1,3,2,4)]
  n<-c(length(x[[1]]),length(x[[2]]),length(x[[3]]),length(x[[4]]))
  nq=length(q)
  output=matrix(NA,nrow=length(q),ncol=8)
  output[,1]=q
  for(j in 1:nq)output[j,2]=hd(outer(x[[1]],x[[2]],FUN='-'),q=q[j])
  for(j in 1:nq)output[j,3]=hd(outer(x[[3]],x[[4]],FUN='-'),q=q[j])
  output[,4]=output[,2]-output[,3]
  e=lapply(q,iband.sub,x=x,nboot=nboot)
  for(j in 1:nq)output[j,5]=e[[j]]$ci[1]
  for(j in 1:nq)output[j,6]=e[[j]]$ci[2]
  for(j in 1:nq)output[j,7]=e[[j]]$p.value
  output[,8]=p.adjust(output[,7],method=method)
  #output=as.numeric(format(round(output,3),nsmall = 3))
  dimnames(output)=list(NULL,c('Quant','Est.Lev 1','Est.Lev 2','Dif','ci.low','ci.up','p-value','p.adj'))
  if(plotit){
    g2plot(outer(x[[1]],x[[2]],FUN='-'),outer(x[[3]],x[[4]],FUN='-'),xlab=xlab)
  }
  output
}

iband.sub<-function(q,x,nboot=500,alpha=.05,SEED=FALSE){
  #
  #
  #
  if(SEED)set.seed(2)
  if(is.matrix(x))x<-listm(x)
  if(length(x)!=4)stop('There must be 4 groups')
  for(j in 1:length(x))x[[j]]=elimna(x[[j]])
  v1=NA
  v2=NA
  B=list()
  for(i in 1:nboot){
    for(j in 1:4)B[[j]]=sample(x[[j]],replace=TRUE)
    v1[i]=hd(outer(B[[1]],B[[2]],FUN='-'),q=q)
    v2[i]=hd(outer(B[[3]],B[[4]],FUN='-'),q=q)
  }
  p=mean(v1<v2)+.5*mean(v1==v2)
  pv=2*min(p,1-p)
  ilow<-round((alpha/2) * nboot)
  ihi<-nboot - ilow
  ilow<-ilow+1
  vs=sort(v1-v2)
  ci=vs[ilow]
  ci[2]=vs[ihi]
  list(ci=ci,p.value=pv)
}
