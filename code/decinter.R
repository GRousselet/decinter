decinter<-function(x,alpha=.05,q=c(1:9)/10,nboot=2000,SEED=TRUE,method='BH'){
  #
  # By default, use  all deciles when dealing with interactions in a 2-by-2 design.
  # The quantiles used can be altered via the argument q
  #
  if(SEED)set.seed(2)
  if(is.matrix(x))x=listm(x)
  x=elimna(x)
  bv1=matrix(NA,nrow=9,ncol=nboot)
  bv2=matrix(NA,nrow=9,ncol=nboot)
  bv3=matrix(NA,nrow=9,ncol=nboot)
  bv4=matrix(NA,nrow=9,ncol=nboot)
  data<-matrix(sample(x[[1]],size=length(x[[1]])*nboot,replace=TRUE),nrow=nboot)
  bv1=apply(data,1,hdmq,q=q)
  data<-matrix(sample(x[[2]],size=length(x[[2]])*nboot,replace=TRUE),nrow=nboot)
  bv2=apply(data,1,hdmq,q=q)
  data<-matrix(sample(x[[3]],size=length(x[[3]])*nboot,replace=TRUE),nrow=nboot)
  bv3=apply(data,1,hdmq,q=q)
  data<-matrix(sample(x[[4]],size=length(x[[4]])*nboot,replace=TRUE),nrow=nboot)
  bv4=apply(data,1,hdmq,q=q)
  be=bv1-bv2-bv3+bv4
  pv=NA
  nq=length(q)
  ilow<-round((alpha/2) * nboot)
  ihi<-nboot - ilow
  ilow<-ilow+1
  vs=sort(be)
  cilow=NA
  ciup=NA
  for(i in 1:nq){
    pv[i]=mean(be[i,]<0)
    pv[i]=2*min(pv[i],1-pv[i])
    bes=sort(be[i,])
    cilow[i]=bes[ilow]
    ciup[i]=bes[ihi]
  }
  output=matrix(NA,nrow=nq,ncol=8)
  dimnames(output)=list(NULL,c('Quant','Est.Lev 1','Est.Lev 2','Dif','ci.low','ci.up','p-value','p.adj'))
  output[,1]=q
  e=lapply(x,hdmq,q=q)
  est=e[[1]]-e[[2]]-e[[3]]+e[[4]]
  output[,2]=e[[1]]-e[[2]]
  output[,3]=e[[3]]-e[[4]]
  output[,4]=est
  output[,5]=cilow
  output[,6]=ciup
  output[,7]=pv
  output[,8]=p.adjust(pv,method=method)
  output
}

hdmq<-function(x,q=.5,tr=FALSE){
  #
  #
  # Estimate one or more quantiles.
  e=NA
  nq=length(q)
  if(!tr)for(i in 1:nq)e[i]=hd(x,q[i])
  if(tr)for(i in 1:nq)e[i]=thd(x,q[i])
  e
}

