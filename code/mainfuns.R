gskew <- function(g){
  #
  # skew and kurtosis of a g-and-h distribution when h=0
  #
  #
  v1=sqrt(3*exp(2*g^2)+exp(3*g^2)-4)
  v2=3*exp(2*g^2)+2*exp(3*g^2)+exp(4*g^2)-3  #Headrick has -6 not	-3, but	based on n=1000000, -3 works
  list(skew=v1,kurtosis=v2)
}

skew <- function(x){
  #
  # Compute skew and kurtosis
  #
  x=elimna(x)
  m1<-mean(x)
  m2<-var(x)
  m3<-sum((x-m1)^3)/length(x)
  m4<-sum((x-m1)^4)/length(x)
  sk<-m3/m2^1.5
  ku<-m4/m2^2
  list(skew=sk,kurtosis=ku)
}

keeporder <- function(x){
  x <- as.character(x)
  x <- factor(x, levels=unique(x))
  x
}

# save beta weights to compute HD
hd.w.calc <- function(nseq, qseq){
  nn <- length(nseq)
  hd.w <- vector("list", nn)
  nq <- length(qseq)
  for(N in 1:nn){
    n <- nseq[N]
    wq <- matrix(NA, nrow = nq, ncol = n)
    for(Q in 1:nq){
      q <- qseq[Q]
      m1 <- (n+1)*q
      m2 <- (n+1)*(1-q)
      vec <- seq(1,n)
      wq[Q,] <- pbeta(vec/n,m1,m2) - pbeta((vec-1)/n,m1,m2)
    }
    hd.w[[N]] <- wq
  }
  hd.w
}

hd <- function(x,q=.5,na.rm=TRUE,STAND=NULL,tr=FALSE){
  #
  #  Compute the Harrell-Davis estimate of the qth quantile
  #
  #  The vector x contains the data,
  #  and the desired quantile is q
  #  The default value for q is .5.
  #
  if(tr)e=thd(x,q=q)
  else{
    if(na.rm)x=elimna(x)
    n<-length(x)
    m1<-(n+1)*q
    m2<-(n+1)*(1-q)
    vec<-seq(along=x)
    w<-pbeta(vec/n,m1,m2)-pbeta((vec-1)/n,m1,m2)  # W sub i values
    y<-sort(x)
    e<-sum(w*y)
  }
  e
}

cnorm<-function(n,epsilon=.1,k=10){
  #
  # generate n observations from a contaminated normal
  # distribution
  # probability 1-epsilon from a standard normal
  # probability epsilon from normal with mean 0 and standard deviation k
  #
  if(epsilon>1)stop("epsilon must be less than or equal to 1")
  if(epsilon<0)stop("epsilon must be greater than or equal to 0")
  if(k<=0)stop("k must be greater than 0")
  val<-rnorm(n)
  uval<-runif(n)
  flag<-(uval<=1-epsilon)
  val[!flag]<-k*val[!flag]
  val
}

clnorm<-function(n,epsilon=.1,k=10){
  #
  # generate n observations from a contaminated lognormal
  # distribution
  #
  #  Using default values, median is approximately 1.14 and 20% trimmed mean is 1.33
  if(epsilon>1)stop('epsilon must be less than or equal to 1')
  if(epsilon<0)stop('epsilon must be greater than or equal to 0')
  if(k<=0)stop('k must be greater than 0')
  val<-rlnorm(n)
  uval<-runif(n)
  flag<-(uval<=1-epsilon)
  val[!flag]<-k*val[!flag]
  val
}

elimna <- function(m){
  #
  # remove any rows of data having missing values
  #
  DONE=FALSE
  if(is.list(m) && is.matrix(m)){
    z=pool.a.list(m)
    m=matrix(z,ncol=ncol(m))
    DONE=TRUE
  }
  if(!DONE){
    if(is.list(m) && is.matrix(m[[1]])){
      for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
      for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    e<-m[ikeep[ikeep>=1],]
  }
  e
}

listm <- function(x){
  #
  # Store the data in a matrix or data frame in a new
  # R variable having list mode.
  # Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
  #
  if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
  y<-list()
  for(j in 1:ncol(x))y[[j]]<-x[,j]
  y
}

deciles <- function(x,HD=TRUE,type=7){
  #
  #  Estimate the deciles for the data in vector x
  #  HD=TRUE: use the Harrell-Davis estimate of the qth quantile
  #   HD=FALSE:use R function quantile
  #
  x=elimna(x)
  if(HD){
    xs<-sort(x)
    n<-length(x)
    vecx<-seq(along=x)
    xq<-0
    for (i in 1:9){
      q<-i/10
      m1<-(n+1)*q
      m2<-(n+1)*(1-q)
      wx<-pbeta(vecx/n,m1,m2)-pbeta((vecx-1)/n,m1,m2)  # W sub i values
      xq[i]<-sum(wx*xs)
    }}
  if(!HD){
    pts=seq(.1,.9,.1)
    xq=quantile(x,probs=pts,type=type)
  }
  xq
}

ghdist<-function(n,g=0,h=0){
  #
  # generate n observations from a g-and-h dist.
  #
  x<-rnorm(n)
  if (g>0){
    ghdist<-(exp(g*x)-1)*exp(h*x^2/2)/g
  }
  if(g==0)ghdist<-x*exp(h*x^2/2)
  ghdist
}

bi2KMSv2<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),
                   x=NA,y=NA,nullval=0,alpha=.05){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success using method KMS.
  #
  #  Unlike the function bi2KMS, a p-value is returned
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  # Uses Kulinskaya et al. method American Statistician, 2010, 64, 350-
  #
  #  null value is the hypothesized value for p1-p2
  #
  alph<-c(1:99)/100
  for(i in 1:99){
    irem<-i
    chkit<-bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=x,alpha=alph[i])
    if(chkit$ci[1]>nullval || chkit$ci[2]<nullval)break
  }
  p.value<-irem/100
  if(p.value<=.1){
    iup<-(irem+1)/100
    alph<-seq(.001,iup,.001)
    for(i in 1:length(alph)){
      p.value<-alph[i]
      chkit<-bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=x,alpha=alph[i])
      if(chkit$ci[1]>nullval || chkit$ci[2]<nullval)break
    }}
  est=bi2KMS(r1=r1,n1=n1,r2=r2,n2=n2,x=x,y=y,alpha=alpha)
  list(ci=est$ci,p1=est$p1,p2=est$p2,est.dif=est$p1-est$p2,p.value=p.value)
}

bi2KMS<-function(r1=sum(elimna(x)),n1=length(elimna(x)),r2=sum(elimna(y)),n2=length(elimna(y)),
                 x=NULL,y=NULL,alpha=.05){
  #
  # Test the hypothesis that two independent binomials have equal
  # probability of success
  #
  # r1=number of successes in group 1
  # n1=number of observations in group 1
  #
  # Use Kulinskaya et al. method American Statistician, 2010, 64, 350-
  #
  #  This function was updated 1/10/19; results might differ very slightly compared to the
  #  original version.
  #
  N=n1+n2
  u=.5
  Dhat=(r1+.5)/(n1+1)-(r2+.5)/(n2+1)
  psihat=((r1+.5)/(n1+1)+(r2+.5)/(n2+1))/2
  nuhat=(1-2*psihat)*(.5-n2/N)
  what=sqrt(2*u*psihat*(1-psihat)+nuhat^2)
  se=qnorm(1-alpha/2)*sqrt(u/(2*n1*n2/N))
  val1=max(c(-1,(u*Dhat+nuhat)/what))
  ci=what*sin(asin(val1)-se)/u-nuhat/u
  val2=min(c(1,(u*Dhat+nuhat)/what))
  ci[2]=what*sin(asin(val2)+se)/u-nuhat/u
  list(ci=ci,p1=r1/n1,p2=r2/n2)
}

rbbinom <- function(n,nbin,r,s){
  #
  # Generate n values from a beta-binomial,
  # r and s are the parameters of the beta distribution.
  # nbin is for the binomial distribution,
  #  Example: nbin=10 means the sample space=c(0:10)
  #
  x<-NA
  for(i in 1:n){
    pval<-rbeta(1,r,s)
    x[i]<-rbinom(1,nbin,pval)
  }
  x
}

wmwpb<-function(x,y=NULL,est=median,alpha=.05,nboot=2000,SEED=TRUE,pr=TRUE,
                na.rm=TRUE,...){
  #
  #   Compute a bootstrap confidence interval for a
  #   measure of location associated with
  #   the distribution of x-y,
  #   est indicates which measure of location will be used
  #   x and y are possibly dependent
  #
  #   loc2dif.ci  computes a non-bootstrap confidence interval
  #
  if(is.null(y[1])){
    if(!is.matrix(x) & !is.data.frame(x))stop('With y missing, x should be a matrix')
    y=x[,2]
    x=x[,1]
  }
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  data1<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
  data2<-matrix(sample(length(y),size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvec<-NA
  for(i in 1:nboot)bvec[i]<-wmwloc(x[data1[i,]],y[data2[i,]],est=est,na.rm=na.rm,...)
  bvec<-sort(bvec)
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  estdiff=wmwloc(x,y,est=est,na.rm=na.rm,...)
  list(estimate=estdiff,ci=c(bvec[low],bvec[up]),p.value=sig.level)
}

wmwloc<-function(x,y,na.rm=TRUE,est=median,...){
  #
  # Estimate the median of the distribution of x-y
  #
  if(na.rm){
    x<-x[!is.na(x)]
    y<-y[!is.na(y)]
  }
  m<-outer(x,y,FUN="-")
  est=est(m,na.rm=TRUE,...)
  est
}

# ------------------------------------------------------
# ANOVA code -------------------------------------------
# ------------------------------------------------------
t2way<-function(J,K,x,tr=.2,grp=c(1:p),p=J*K,MAT=FALSE,
                lev.col=c(1:2),var.col=3,pr=TRUE,IV1=NULL,IV2=NULL){
  #  Perform a J by K  (two-way) ANOVA on trimmed means where
  #  all groups are independent.
  #
  #  The R variable x is assumed to contain the raw
  #  data stored in list mode, or a matrix with columns
  #  corresponding to groups. If stored in list mode, x[[1]] contains the data
  #  for the first level of both factors: level 1,1,.
  #  x[[2]] is assumed to contain the data for level 1 of the
  #  first factor and level 2 of the second factor: level 1,2
  #
  #  The default amount of trimming is tr=.2
  #
  #  It is assumed that x has length JK, the total number of
  #  groups being tested.
  #
  #  MAT=T, assumes x are stored in matrix with 3 columns
  #  with two of the columns indicated by the argument
  #  lev.col
  #  specifying the columns of x containing the values of the
  #  levels of the two factors.
  #  The outcome variable is in column
  #  var.col
  #  which defaults to column 3
  #  That is, this function calls selby2 for you.
  #
  #  IV1 and IV2: if specified, taken to be the independent variable
  #      That is, the group id values
  #      and x is assumed to be a vector containing all of the data
  #  EXAMPLE: t2way(x=data,IV1=iv1,IV2=iv2)
  #  would do a two-way ANOVA based on group id's in iv1 and iv2 and
  #  dependent variable data
  #
  if(is.data.frame(x))data=as.matrix(x)
  if(tr==.5){
    print("For medians, use med2way if there are no ties")
    print("With ties, use linear contrasts in conjunction with medpb")
    stop("")
  }
  if(MAT){
    if(!is.matrix(x))stop("With MAT=T, data must be a matrix")
    if(length(lev.col)!=2)stop("Argument lev.col should have 3 values")
    temp=selby2(x,lev.col,var.col)
    lev1=length(unique(temp$grpn[,1]))
    lev2=length(unique(temp$grpn[,2]))
    gv=apply(temp$grpn,2,rank)
    gvad=10*gv[,1]+gv[,2]
    grp=rank(gvad)
    if(pr){
      print(paste("Factor 1 has", lev1, "levels"))
      print(paste("Factor 2 has", lev2, "levels"))
    }
    if(J!=lev1)warning("J is being reset to the number of levels found")
    if(K!=lev2)warning("K is being reset to the number of levels found")
    J=lev1
    K=lev2
    x=temp$x
  }
  if(!is.null(IV1[1])){
    if(is.null(IV2[1]))stop("IV2 is NULL")
    if(pr)print("Assuming data is a vector containing all of the data; the dependent variable")
    xi=elimna(cbind(x,IV1,IV2))
    J=length(unique(xi[,2]))
    K=length(unique(xi[,3]))
    x=fac2list(xi[,1],xi[,2:3])
  }
  if(is.matrix(x))x=listm(x)
  if(!is.list(x))stop("Data are not stored in list mode")
  if(p!=length(x)){
    print("The total number of groups, based on the specified levels, is")
    print(p)
    print("The number of groups is")
    print(length(x))
    print("Warning: These two values are not equal")
  }
  tmeans<-0
  h<-0
  v<-0
  for (i in 1:p){
    x[[grp[i]]]=elimna(x[[grp[i]]])
    tmeans[i]<-mean(x[[grp[i]]],tr)
    h[i]<-length(x[[grp[i]]])-2*floor(tr*length(x[[grp[i]]]))
    #    h is the effective sample size
    if(winvar(x[[grp[i]]],tr)==0)print(paste('The Winsorized variance is zero for group',i))
    v[i]<-(length(x[[grp[i]]])-1)*winvar(x[[grp[i]]],tr)/(h[i]*(h[i]-1))
    #    v contains the squared standard errors
  }
  v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
  ij<-matrix(c(rep(1,J)),1,J)
  ik<-matrix(c(rep(1,K)),1,K)
  jm1<-J-1
  cj<-diag(1,jm1,J)
  for (i in 1:jm1)cj[i,i+1]<-0-1
  km1<-K-1
  ck<-diag(1,km1,K)
  for (i in 1:km1)ck[i,i+1]<-0-1
  #  Do test for factor A
  cmat<-kron(cj,ik)  # Contrast matrix for factor A
  alval<-c(1:999)/1000
  for(i in 1:999){
    irem<-i
    Qa<-johan(cmat,tmeans,v,h,alval[i])
    if(i==1)dfA=Qa$df
    if(Qa$teststat>Qa$crit)break
  }
  A.p.value=irem/1000
  # Do test for factor B
  cmat<-kron(ij,ck)  # Contrast matrix for factor B
  for(i in 1:999){
    irem<-i
    Qb<-johan(cmat,tmeans,v,h,alval[i])
    if(i==1)dfB=Qb$df
    if(Qb$teststat>Qb$crit)break
  }
  B.p.value=irem/1000
  # Do test for factor A by B interaction
  cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
  for(i in 1:999){
    irem<-i
    Qab<-johan(cmat,tmeans,v,h,alval[i])
    if(i==1)dfAB=Qab$df
    if(Qab$teststat>Qab$crit)break
  }
  AB.p.value=irem/1000
  tmeans=matrix(tmeans,J,K,byrow=TRUE)
  list(Qa=Qa$teststat,A.p.value=A.p.value, df.A=dfA,
       Qb=Qb$teststat,B.p.value=B.p.value,df.B=dfB,
       Qab=Qab$teststat,AB.p.value=AB.p.value,df.AB=dfAB,means=tmeans)
}

fac2list<-function(x,g,pr=TRUE){
  #
  # data are stored in x
  # information about the level of the value in x is stored in g,
  # which can be a matrix with up to 4 columns
  #
  # sort the data in x into groups based on values in g.
  # store results in list mode.
  #
  #  Example: fac2list(m[,2],m[,4]) would sort the values
  #  in column 2 of m according to the values in column 4 of m
  #
  g=as.data.frame(g)
  L=ncol(g)
  g=listm(g)
  for(j in 1:L)g[[j]]=as.factor(g[[j]])
  g=matl(g)
  Lp1=L+1
  if(L>4)stop("Can have at most 4 factors")
  if(L==1){
    res=selby(cbind(x,g),2,1)
    group.id=res$grpn
    res=res$x
  }
  if(L>1){
    res=selby2(cbind(x,g),c(2:Lp1),1)
    group.id=res$grpn
    res=res$x
  }
  if(pr)
  {print("Group Levels:")
    print(group.id)
  }
  res
}

listm<-function(x){
  #
  # Store the data in a matrix or data frame in a new
  # R variable having list mode.
  # Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
  #
  if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
  y<-list()
  for(j in 1:ncol(x))y[[j]]<-x[,j]
  y
}

selby<-function(m,grpc,coln){
  #
  #
  #  A commmon situation is to have data stored in an n by p matrix where
  #  one or more of the columns are  group identification numbers.
  #  This function groups  all values in column coln according to the
  #  group numbers in column grpc and stores the  results in list mode.
  #
  #  More than one column of data can sorted
  #
  # grpc indicates the column of the matrix containing group id number
  #
  if(is.null(dim(m)))stop("Data must be stored in a matrix or data frame")
  if(is.na(grpc[1]))stop("The argument grpc is not specified")
  if(is.na(coln[1]))stop("The argument coln is not specified")
  if(length(grpc)!=1)stop("The argument grpc must have length 1")
  x<-vector("list")
  grpn<-sort(unique(m[,grpc]))
  it<-0
  for (ig in 1:length(grpn)){
    for (ic in 1:length(coln)){
      it<-it+1
      flag<-(m[,grpc]==grpn[ig])
      x[[it]]<-m[flag,coln[ic]]
    }}
  list(x=x,grpn=grpn)
}

selby2<-function(m,grpc,coln=NA){
  # Create categories according to the grpc[1] and grpc[2] columns
  # of the matrix m. The function puts the values in column coln into
  # a vector having list mode.
  #
  if(is.na(coln))stop("The argument coln is not specified")
  if(length(grpc)>4)stop("The argument grpc must have length less than or equal to 4")
  x<-vector("list")
  ic<-0
  if(length(grpc)==2){
    cat1<-selby(m,grpc[1],coln)$grpn
    cat2<-selby(m,grpc[2],coln)$grpn
    for (i1 in 1:length(cat1)){
      for (i2 in 1:length(cat2)){
        temp<-NA
        it<-0
        for (i in 1:nrow(m)){
          if(sum(m[i,c(grpc[1],grpc[2])]==c(cat1[i1],cat2[i2]))==2){
            it<-it+1
            temp[it]<-m[i,coln]
          }
        }
        if(!is.na(temp[1])){
          ic<-ic+1
          x[[ic]]<-temp
          if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2]),1,2)
          if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2]))
        }
      }}
  }
  if(length(grpc)==3){
    cat1<-selby(m,grpc[1],coln)$grpn
    cat2<-selby(m,grpc[2],coln)$grpn
    cat3<-selby(m,grpc[3],coln)$grpn
    x<-vector("list")
    ic<-0
    for (i1 in 1:length(cat1)){
      for (i2 in 1:length(cat2)){
        for (i3 in 1:length(cat3)){
          temp<-NA
          it<-0
          for (i in 1:nrow(m)){
            if(sum(m[i,c(grpc[1],grpc[2],grpc[3])]==c(cat1[i1],cat2[i2],cat3[i3]))==3){
              it<-it+1
              temp[it]<-m[i,coln]
            }}
          if(!is.na(temp[1])){
            ic<-ic+1
            x[[ic]]<-temp
            if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3]),1,3)
            if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3]))
          }}}}
  }
  if(length(grpc)==4){
    cat1<-selby(m,grpc[1],coln)$grpn
    cat2<-selby(m,grpc[2],coln)$grpn
    cat3<-selby(m,grpc[3],coln)$grpn
    cat4<-selby(m,grpc[4],coln)$grpn
    x<-vector("list")
    ic<-0
    for (i1 in 1:length(cat1)){
      for (i2 in 1:length(cat2)){
        for (i3 in 1:length(cat3)){
          for (i4 in 1:length(cat4)){
            temp<-NA
            it<-0
            for (i in 1:nrow(m)){
              if(sum(m[i,c(grpc[1],grpc[2],grpc[3],grpc[4])]==c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))==4){
                it<-it+1
                temp[it]<-m[i,coln]
              }}
            if(!is.na(temp[1])){
              ic<-ic+1
              x[[ic]]<-temp
              if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]),1,4)
              if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))
            }}}}}
  }
  list(x=x,grpn=grpn)
}

matl<-function(x){
  #
  # take data in list mode and store it in a matrix
  #
  J=length(x)
  nval=NA
  for(j in 1:J)nval[j]=length(x[[j]])
  temp<-matrix(NA,ncol=J,nrow=max(nval))
  for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
  temp
}

winvar<-function(x,tr=.2,na.rm=FALSE,STAND=NULL){
  #
  #  Compute the gamma Winsorized variance for the data in the vector x.
  #  tr is the amount of Winsorization which defaults to .2.
  #
  remx=x
  x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  wv<-var(y)
  if(!na.rm)if(sum(is.na(remx)>0))wv=NA
  wv
}

kron<-function(m1,m2){
  #  compute the Kronecker product of the two matrices m1 and m2.
  #
  m1<-as.matrix(m1) # Vectors of length p are converted to a p by 1 matrix
  m2<-as.matrix(m2)
  kron<-vector(mode="numeric",length=0)
  for(i in 1:nrow(m1)){
    m3<-m1[i,1]*m2
    for(j in 2:ncol(m1))m3<-cbind(m3,m1[i,j]*m2)
    if(i==1)kron<-m3
    if(i>=2)kron<-rbind(kron,m3)
  }
  kron
}

johan<-function(cmat,vmean,vsqse,h,alpha=.05){
  #
  #  This function is used by other functions that come with this book,
  #  and it can be used to test hypotheses not covered in the text.
  #
  #  The function performs Johansen's test of C mu = 0 for p independent groups,
  #  where C is a k by p matrix of rank k and mu is a p by 1 matrix of
  #  of unknown trimmed means.
  #  The argument cmat contains the matrix C.
  #  vmean is a vector of length p containing the p trimmed means
  #  vsqe is a diagonal matrix containing the squared standard errors of the
  #  the trimmed means in vmean.
  #  h is a vector containing the effective sample sizes
  #
  yvec<-matrix(vmean,length(vmean),1)
  if(!is.matrix(vsqse))vsqse<-diag(vsqse)
  test<-cmat%*%vsqse%*%t(cmat)
  invc<-solve(test)
  test<-t(yvec)%*%t(cmat)%*%invc%*%cmat%*%yvec
  R<-vsqse%*%t(cmat)%*%invc%*%cmat
  A<-sum(diag((diag(R))^2/diag(h-1)))
  df<-nrow(cmat)
  crit<-qchisq(1-alpha,df)
  crit<-crit+(crit/(2*df))*A*(1+3*crit/(df+2))
  list(teststat=test[1],crit=crit[1],df=df)
}
# ---------------------------------------------

# ---------------------------------------------------------------------------
# SIMULATION CODE
# ---------------------------------------------------------------------------

decinter_wrap <- function(x, n, hd.w, nboot){
  res <- decinter_cpp(x, n, hd.w, nboot)
  diff <- res[["diff"]]
  mainA <- res[["mainA"]]
  mainB <- res[["mainB"]]
  bootdiff <- res[["bootdiff"]]
  bootmainA <- res[["bootmainA"]]
  bootmainB <- res[["bootmainB"]]
  # Compute p values using percentile bootstrap:
  temp <- colSums(bootdiff<0) / nboot + colSums(bootdiff==0) / (2*nboot)
  bootdiff.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainA<0) / nboot + colSums(bootmainA==0) / (2*nboot)
  bootmainA.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainB<0) / nboot + colSums(bootmainB==0) / (2*nboot)
  bootmainB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  list(AB = bootdiff.pv, A = bootmainA.pv, B = bootmainB.pv)
}

sim.counter <- function(S, nsim, inc){
  if(S == 1){
    # print(paste(nsim,"iterations:",S))
    cat(nsim,"iterations:",S)
    beep(2)
  }
  if(S %% inc == 0){
    # print(paste("iteration",S,"/",nsim))
    cat(" /",S)
    beep(2)
  }
}

# decinter: false positive simulation
# compare HD to QT7
# compare FDR to Hochberg
# Includes ANOVAs on means and 20% trimmed means for comparison
decinter_fp_sim <- function(rdist, ...){
  
  print("decinter: false positive simulation...")
  
  # declare variables  
  A.qt <- array(0, dim = c(nsim, nn, 9)) # Main effect of A
  B.qt <- array(0, dim = c(nsim, nn, 9)) # Main effect of B
  AB.qt <- array(0, dim = c(nsim, nn, 9)) # interaction
  A.hd <- array(0, dim = c(nsim, nn, 9)) # Main effect of A
  B.hd <- array(0, dim = c(nsim, nn, 9)) # Main effect of B
  AB.hd <- array(0, dim = c(nsim, nn, 9)) # interaction
  ANOVA.m <- array(0, dim = c(nsim, nn, 3))
  ANOVA.tm <- array(0, dim = c(nsim, nn, 3))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(rdist(nmax * 4), nrow = nmax)
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      # x <- matrix(rnorm(nx * 4), nrow = nx)
      pv <- decinter_hdqt7_wrap(x, nx, hd.w[[N]], nboot, qseq)
      # bootstrap p values
      A.qt[S,N,] <- pv$qt.A
      B.qt[S,N,] <- pv$qt.B
      AB.qt[S,N,] <- pv$qt.AB
      A.hd[S,N,] <- pv$hd.A
      B.hd[S,N,] <- pv$hd.B
      AB.hd[S,N,] <- pv$hd.AB

      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      res <- t2way(2,2,xLIST,tr=0)
      ANOVA.m[S,N,1] <- res$A.p.value
      ANOVA.m[S,N,2] <- res$B.p.value
      ANOVA.m[S,N,3] <- res$AB.p.value
      res <- t2way(2,2,xLIST,tr=0.2)
      ANOVA.tm[S,N,1] <- res$A.p.value
      ANOVA.tm[S,N,2] <- res$B.p.value
      ANOVA.tm[S,N,3] <- res$AB.p.value
    }
  }
  
  beep(8)
  
  list(A.qt = A.qt, B.qt = B.qt, AB.qt = AB.qt, 
       A.hd = A.hd, B.hd = B.hd, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# called by decinter_fp_sim and decinter_tp_sim
decinter_hdqt7_wrap <- function(x, n, w, nboot, probs){
  res <- decinter_hdqt7_cpp(x, n, w, nboot, probs)
  # QT7 ------------------------------------------
  AB <- res[["ABqt"]]
  A <- res[["Aqt"]]
  B <- res[["Bqt"]]
  bootAB <- res[["bootABqt"]]
  bootA <- res[["bootAqt"]]
  bootB <- res[["bootBqt"]]
  # Compute percentile bootstrap p values:
  temp <- colSums(bootAB<0) / nboot + colSums(bootAB==0) / (2*nboot)
  qt.bootAB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootA<0) / nboot + colSums(bootA==0) / (2*nboot)
  qt.bootA.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootB<0) / nboot + colSums(bootB==0) / (2*nboot)
  qt.bootB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  
  # HD ------------------------------------------
  diff <- res[["diffhd"]]
  mainA <- res[["mainAhd"]]
  mainB <- res[["mainBhd"]]
  bootAB <- res[["bootABhd"]]
  bootA <- res[["bootAhd"]]
  bootB <- res[["bootBhd"]]
  # Compute percentile bootstrap p values:
  temp <- colSums(bootAB<0) / nboot + colSums(bootAB==0) / (2*nboot)
  hd.bootAB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootA<0) / nboot + colSums(bootA==0) / (2*nboot)
  hd.bootA.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootB<0) / nboot + colSums(bootB==0) / (2*nboot)
  hd.bootB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  
  # Return results ------------------------------------------
  list(qt.AB = qt.bootAB.pv, qt.A = qt.bootA.pv, qt.B = qt.bootB.pv,
       hd.AB = hd.bootAB.pv, hd.A = hd.bootA.pv, hd.B = hd.bootB.pv)
}

# use results from decinter_fp_sim() to compute FWER
# also used to compute power in sim_tp.Rmd
decinter_fwer <- function(simres, ...){
  
  print("decinter: compute FWER...")
  
  # Main effect of A
  A.qt.raw <- array(0, dim = c(nsim, nn))
  A.qt.fdr <- array(0, dim = c(nsim, nn))
  A.qt.hoch <- array(0, dim = c(nsim, nn))
  
  # Main effect of B
  B.qt.raw <- array(0, dim = c(nsim, nn))
  B.qt.fdr <- array(0, dim = c(nsim, nn))
  B.qt.hoch <- array(0, dim = c(nsim, nn))
  
  # interaction
  AB.qt.raw <- array(0, dim = c(nsim, nn))
  AB.qt.fdr <- array(0, dim = c(nsim, nn))
  AB.qt.hoch <- array(0, dim = c(nsim, nn))
  
  # Main effect of A
  A.hd.raw <- array(0, dim = c(nsim, nn))
  A.hd.fdr <- array(0, dim = c(nsim, nn))
  A.hd.hoch <- array(0, dim = c(nsim, nn))
  
  # Main effect of B
  B.hd.raw <- array(0, dim = c(nsim, nn))
  B.hd.fdr <- array(0, dim = c(nsim, nn))
  B.hd.hoch <- array(0, dim = c(nsim, nn))
  
  # interaction
  AB.hd.raw <- array(0, dim = c(nsim, nn))
  AB.hd.fdr <- array(0, dim = c(nsim, nn))
  AB.hd.hoch <- array(0, dim = c(nsim, nn))
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    for(N in 1:nn){
      p.val <- simres$A.qt[S,N,]
      A.qt.raw[S,N] <- sum(p.val < alpha.val) >= 1
      A.qt.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      A.qt.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      p.val <- simres$B.qt[S,N,]
      B.qt.raw[S,N] <- sum(p.val < alpha.val) >= 1
      B.qt.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      B.qt.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      p.val <- simres$AB.qt[S,N,]
      AB.qt.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.qt.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      AB.qt.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      p.val <- simres$A.hd[S,N,]
      A.hd.raw[S,N] <- sum(p.val < alpha.val) >= 1
      A.hd.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      A.hd.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      p.val <- simres$B.hd[S,N,]
      B.hd.raw[S,N] <- sum(p.val < alpha.val) >= 1
      B.hd.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      B.hd.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      p.val <- simres$AB.hd[S,N,]
      AB.hd.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.hd.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      AB.hd.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
    }
  }
  list(AB.qt.raw = AB.qt.raw, 
       AB.qt.fdr = AB.qt.fdr,
       AB.qt.hoch = AB.qt.hoch,
       AB.hd.raw = AB.hd.raw, 
       AB.hd.fdr = AB.hd.fdr,
       AB.hd.hoch = AB.hd.hoch,
       A.qt.raw = A.qt.raw,
       A.qt.fdr = A.qt.fdr,
       A.qt.hoch = A.qt.hoch,
       A.hd.raw = A.hd.raw,
       A.hd.fdr = A.hd.fdr,
       A.hd.hoch = A.hd.hoch,
       B.qt.raw = B.qt.raw,
       B.qt.fdr = B.qt.fdr,
       B.qt.hoch = B.qt.hoch,
       B.hd.raw = B.hd.raw,
       B.hd.fdr = B.hd.fdr,
       B.hd.hoch = B.hd.hoch)
}

# decinter: false and true positive simulation for Poisson distributions
# compare HD to QT7
# compare FDR to Hochberg
# Includes ANOVAs on means and 20% trimmed means for comparison
decinter_pois_sim <- function(lambda, qseq, ...){
  
  print("decinter: poisson simulation...")
  
  # declare variables  
  A.qt <- array(0, dim = c(nsim, nn, 9)) # Main effect of A
  B.qt <- array(0, dim = c(nsim, nn, 9)) # Main effect of B
  AB.qt <- array(0, dim = c(nsim, nn, 9)) # interaction
  A.hd <- array(0, dim = c(nsim, nn, 9)) # Main effect of A
  B.hd <- array(0, dim = c(nsim, nn, 9)) # Main effect of B
  AB.hd <- array(0, dim = c(nsim, nn, 9)) # interaction
  ANOVA.m <- array(0, dim = c(nsim, nn, 3))
  ANOVA.tm <- array(0, dim = c(nsim, nn, 3))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(0, nrow = nmax, ncol = 4)
    for(C in 1:4){
      x.all[,C] <- rpois(nmax, lambda[C])
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      # x <- matrix(rnorm(nx * 4), nrow = nx)
      pv <- decinter_hdqt7_wrap(x, nx, hd.w[[N]], nboot, qseq)
      # bootstrap p values
      A.qt[S,N,] <- pv$qt.A
      B.qt[S,N,] <- pv$qt.B
      AB.qt[S,N,] <- pv$qt.AB
      A.hd[S,N,] <- pv$hd.A
      B.hd[S,N,] <- pv$hd.B
      AB.hd[S,N,] <- pv$hd.AB
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      res <- t2way(2,2,xLIST,tr=0)
      ANOVA.m[S,N,1] <- res$A.p.value
      ANOVA.m[S,N,2] <- res$B.p.value
      ANOVA.m[S,N,3] <- res$AB.p.value
      res <- t2way(2,2,xLIST,tr=0.2)
      ANOVA.tm[S,N,1] <- res$A.p.value
      ANOVA.tm[S,N,2] <- res$B.p.value
      ANOVA.tm[S,N,3] <- res$AB.p.value
    }
  }
  
  beep(8)
  
  list(A.qt = A.qt, B.qt = B.qt, AB.qt = AB.qt, 
       A.hd = A.hd, B.hd = B.hd, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# Beta-binomial simulations
# Used for false positive and power simulations
decinter_bb_sim <- function(rseq, s, nbin, qseq, nseq, ...){
  
  print("decinter: beta-binomial simulation...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of A
  A.qt <- array(0, dim = c(nsim, nn, 9))
  A.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of B
  B.qt <- array(0, dim = c(nsim, nn, 9))
  B.hd <- array(0, dim = c(nsim, nn, 9))
  
  # ANOVA
  ANOVA.m <- array(0, dim = c(nsim, nn, 3))
  ANOVA.tm <- array(0, dim = c(nsim, nn, 3))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(0, nrow = nmax, ncol = 4)
    for(C in 1:4){ # add constants
      x.all[,C] <- rbbinom(nmax, nbin = nbin, r = rseq[C], s = s)
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      pv <- decinter_hdqt7_wrap(x, nx, hd.w[[N]], nboot, qseq)
      # QT7 -------------------
      AB.qt[S,N,] <- pv$qt.AB
      A.qt[S,N,] <- pv$qt.A
      B.qt[S,N,] <- pv$qt.B
      # HD -------------------
      AB.hd[S,N,] <- pv$hd.AB
      A.hd[S,N,] <- pv$hd.A
      B.hd[S,N,] <- pv$hd.B
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      res <- t2way(2,2,xLIST,tr=0)
      ANOVA.m[S,N,1] <- res$A.p.value
      ANOVA.m[S,N,2] <- res$B.p.value
      ANOVA.m[S,N,3] <- res$AB.p.value
      res <- t2way(2,2,xLIST,tr=0.2)
      ANOVA.tm[S,N,1] <- res$A.p.value
      ANOVA.tm[S,N,2] <- res$B.p.value
      ANOVA.tm[S,N,3] <- res$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.qt = AB.qt, A.qt = A.qt, B.qt = B.qt, 
       AB.hd = AB.hd, A.hd = A.hd, B.hd = B.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}


# Power simulation for differences in means
# Includes ANOVAs on means and 20% trimmed means for comparison
# Used in sim_tp.Rmd
sim_mean_diff <- function(rdist, gpmeans, ...){
  
  print("True positive simulation: mean differences...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of A
  A.qt <- array(0, dim = c(nsim, nn, 9))
  A.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of B
  B.qt <- array(0, dim = c(nsim, nn, 9))
  B.hd <- array(0, dim = c(nsim, nn, 9))
  
  ANOVA.m <- array(0, dim = c(nsim, nn, 3))
  ANOVA.tm <- array(0, dim = c(nsim, nn, 3))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(rdist(nmax * 4), nrow = nmax)
    for(C in 1:4){ # add constants
      x.all[,C] <- x.all[,C] + gpmeans[C]  
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      pv <- decinter_hdqt7_wrap(x, nx, hd.w[[N]], nboot, qseq)
      # HD
      AB.hd[S,N,] <- pv$hd.AB
      A.hd[S,N,] <- pv$hd.A
      B.hd[S,N,] <- pv$hd.B
      # QT
      AB.qt[S,N,] <- pv$qt.AB
      A.qt[S,N,] <- pv$qt.A
      B.qt[S,N,] <- pv$qt.B
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      res <- t2way(2,2,xLIST,tr=0)
      ANOVA.m[S,N,1] <- res$A.p.value
      ANOVA.m[S,N,2] <- res$B.p.value
      ANOVA.m[S,N,3] <- res$AB.p.value
      res <- t2way(2,2,xLIST,tr=0.2)
      ANOVA.tm[S,N,1] <- res$A.p.value
      ANOVA.tm[S,N,2] <- res$B.p.value
      ANOVA.tm[S,N,3] <- res$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.hd = AB.hd, A.hd = A.hd, B.hd = B.hd, 
       AB.qt = AB.qt, A.qt = A.qt, B.qt = B.qt,
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# Power simulation for differences in g parameter of g-and-h distributions
# Includes ANOVAs on means and 20% trimmed means for comparison
sim_gh_diff <- function(gvals, hvals, ...){
  
  print("True positive simulation: g differences...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of A
  A.qt <- array(0, dim = c(nsim, nn, 9))
  A.hd <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of B
  B.qt <- array(0, dim = c(nsim, nn, 9))
  B.hd <- array(0, dim = c(nsim, nn, 9))
  
  ANOVA.m <- array(0, dim = c(nsim, nn, 3))
  ANOVA.tm <- array(0, dim = c(nsim, nn, 3))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(NA, nrow = nmax, ncol = 4)
    for(C in 1:4){ # add constants
      x.all[,C] <- ghdist(nmax, gvals[C], hvals[C])
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      # x <- matrix(rnorm(nx * 4), nrow = nx)
      pv <- decinter_hdqt7_wrap(x, nx, hd.w[[N]], nboot, qseq)
      # HD
      AB.hd[S,N,] <- pv$hd.AB
      A.hd[S,N,] <- pv$hd.A
      B.hd[S,N,] <- pv$hd.B
      # QT
      AB.qt[S,N,] <- pv$qt.AB
      A.qt[S,N,] <- pv$qt.A
      B.qt[S,N,] <- pv$qt.B
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      res <- t2way(2,2,xLIST,tr=0)
      ANOVA.m[S,N,1] <- res$A.p.value
      ANOVA.m[S,N,2] <- res$B.p.value
      ANOVA.m[S,N,3] <- res$AB.p.value
      res <- t2way(2,2,xLIST,tr=0.2)
      ANOVA.tm[S,N,1] <- res$A.p.value
      ANOVA.tm[S,N,2] <- res$B.p.value
      ANOVA.tm[S,N,3] <- res$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.hd = AB.hd, A.hd = A.hd, B.hd = B.hd, 
       AB.qt = AB.qt, A.qt = A.qt, B.qt = B.qt,
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}


# ----------------------------------------------
decinter_b1b9_wrap <- function(x, n, hd.w, nboot){
  res <- decinter_b1b9_cpp(x, n, hd.w, nboot)
  diff <- res[["diff"]]
  mainA <- res[["mainA"]]
  mainB <- res[["mainB"]]
  bootdiff_b1 <- res[["bootdiff_b1"]]
  bootmainA_b1 <- res[["bootmainA_b1"]]
  bootmainB_b1 <- res[["bootmainB_b1"]]
  bootdiff_b9 <- res[["bootdiff_b9"]]
  bootmainA_b9 <- res[["bootmainA_b9"]]
  bootmainB_b9 <- res[["bootmainB_b9"]]
  # Compute p values:
  # percentile bootstrap: boot1
  temp <- colSums(bootdiff_b1<0) / nboot + colSums(bootdiff_b1==0) / (2*nboot)
  bootdiffb1.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainA_b1<0) / nboot + colSums(bootmainA_b1==0) / (2*nboot)
  bootmainAb1.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainB_b1<0) / nboot + colSums(bootmainB_b1==0) / (2*nboot)
  bootmainBb1.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # percentile bootstrap: boot9
  temp <- colSums(bootdiff_b9<0) / nboot + colSums(bootdiff_b9==0) / (2*nboot)
  bootdiffb9.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainA_b9<0) / nboot + colSums(bootmainA_b9==0) / (2*nboot)
  bootmainAb9.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  temp <- colSums(bootmainB_b9<0) / nboot + colSums(bootmainB_b9==0) / (2*nboot)
  bootmainBb9.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # export
  list(b1.int = bootdiffb1.pv, b1.A = bootmainAb1.pv, b1.B = bootmainBb1.pv,
       b9.int = bootdiffb9.pv, b9.A = bootmainAb9.pv, b9.B = bootmainBb9.pv)
}

# Same as decinter_fp_sim but we compare two types
# of bootstrap sampling.
decinter_b1b9_sim <- function(rdist, ...){
  
  print("False positive simulation:")
  
  # interaction
  AB.b9 <- array(0, dim = c(nsim, nn, 9))
  AB.b1 <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of A
  A.b9 <- array(0, dim = c(nsim, nn, 9))
  A.b1 <- array(0, dim = c(nsim, nn, 9))
  
  # Main effect of B
  B.b9 <- array(0, dim = c(nsim, nn, 9))
  B.b1 <- array(0, dim = c(nsim, nn, 9))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(rdist(nmax * 4), nrow = nmax)
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      # x <- matrix(rnorm(nx * 4), nrow = nx)
      pv <- decinter_b1b9_wrap(x, nx, hd.w[[N]], nboot)
      # one bootstrap sample
      AB.b1[S,N,] <- pv$b1.int
      A.b1[S,N,] <- pv$b1.A
      B.b1[S,N,] <- pv$b1.B
      # nine bootstrap samples
      AB.b9[S,N,] <- pv$b9.int
      A.b9[S,N,] <- pv$b9.A
      B.b9[S,N,] <- pv$b9.B
    }
  }
  
  beep(8)
  
  list(AB.b1 = AB.b1, A.b1 = A.b1, B.b1 = B.b1, 
       AB.b9 = AB.b9, A.b9 = A.b9, B.b9 = B.b9)
}

# -----------------------------------------
decinter_b1b9_fwer <- function(simres, ...){
  
  print("compute fwer:")
  
  # interaction
  AB.b1.raw <- array(0, dim = c(nsim, nn))
  AB.b1.fdr <- array(0, dim = c(nsim, nn))
  AB.b9.raw <- array(0, dim = c(nsim, nn))
  AB.b9.fdr <- array(0, dim = c(nsim, nn))
  
  # Main effect of A
  A.b1.raw <- array(0, dim = c(nsim, nn))
  A.b1.fdr <- array(0, dim = c(nsim, nn))
  A.b9.raw <- array(0, dim = c(nsim, nn))
  A.b9.fdr <- array(0, dim = c(nsim, nn))
  
  # Main effect of B
  B.b1.raw <- array(0, dim = c(nsim, nn))
  B.b1.fdr <- array(0, dim = c(nsim, nn))
  B.b9.raw <- array(0, dim = c(nsim, nn))
  B.b9.fdr <- array(0, dim = c(nsim, nn))
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    for(N in 1:nn){
      # percentile bootstrap: boot 1
      p.val <- simres$AB.b1[S,N,]
      AB.b1.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.b1.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      p.val <- simres$A.b1[S,N,]
      A.b1.raw[S,N] <- sum(p.val < alpha.val) >= 1
      A.b1.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      p.val <- simres$B.b1[S,N,]
      B.b1.raw[S,N] <- sum(p.val < alpha.val) >= 1
      B.b1.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      
      # percentile bootstrap: boot 9
      p.val <- simres$AB.b9[S,N,]
      AB.b9.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.b9.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      p.val <- simres$A.b9[S,N,]
      A.b9.raw[S,N] <- sum(p.val < alpha.val) >= 1
      A.b9.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      p.val <- simres$B.b9[S,N,]
      B.b9.raw[S,N] <- sum(p.val < alpha.val) >= 1
      B.b9.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
    }
  }
  list(AB.b9.raw = AB.b9.raw, 
       AB.b9.fdr = AB.b9.fdr,
       AB.b1.raw = AB.b1.raw, 
       AB.b1.fdr = AB.b1.fdr,
       A.b9.raw = A.b9.raw,
       A.b9.fdr = A.b9.fdr,
       A.b1.raw = A.b1.raw,
       A.b1.fdr = A.b1.fdr,
       B.b9.raw = B.b9.raw,
       B.b9.fdr = B.b9.fdr,
       B.b1.raw = B.b1.raw,
       B.b1.fdr = B.b1.fdr)
}

# ----------------------------------------------
# APDINTER: COMPARE BOOT1 TO BOOT9 -- use HD and FDR only
apdinter_b1b9_wrap <- function(x, n, hd.w, nboot){
  res <- apdinter_b1b9_cpp(x, n, hd.w, nboot)
  bootdiff_b1 <- res[["bootdiff_b1"]]
  bootdiff_b9 <- res[["bootdiff_b9"]]
  # Compute p values:
  # percentile bootstrap: boot1
  temp <- colSums(bootdiff_b1<0) / nboot + colSums(bootdiff_b1==0) / (2*nboot)
  bootdiffb1.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # percentile bootstrap: boot9
  temp <- colSums(bootdiff_b9<0) / nboot + colSums(bootdiff_b9==0) / (2*nboot)
  bootdiffb9.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # export
  list(b1.int = bootdiffb1.pv, b9.int = bootdiffb9.pv)
}

apdinter_b1b9_sim <- function(rdist, n, hd.w, ...){
  
  print("False positive simulation:")
  
  # interaction
  AB.b9 <- array(0, dim = c(nsim, 9))
  AB.b1 <- array(0, dim = c(nsim, 9))
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x <- matrix(rdist(n * 4), nrow = n)
    pv <- apdinter_b1b9_wrap(x, n, hd.w, nboot)
    # percentile bootstrap
    AB.b1[S,] <- pv$b1.int
    AB.b9[S,] <- pv$b9.int
  }
  
  beep(8)
  
  list(AB.b1 = AB.b1, AB.b9 = AB.b9)
}

apdinter_b1b9_fwer <- function(simres, ...){
  
  print("compute fwer:")
  
  # interaction
  AB.b1.raw <- array(0, dim = c(nsim))
  AB.b1.fdr <- array(0, dim = c(nsim))
  AB.b9.raw <- array(0, dim = c(nsim))
  AB.b9.fdr <- array(0, dim = c(nsim))
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # percentile bootstrap: boot 1
    p.val <- simres$AB.b1[S,]
    AB.b1.raw[S] <- sum(p.val < alpha.val) >= 1
    AB.b1.fdr[S] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
    
    # percentile bootstrap: boot 9
    p.val <- simres$AB.b9[S,]
    AB.b9.raw[S] <- sum(p.val < alpha.val) >= 1
    AB.b9.fdr[S] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
    
  }
  list(AB.b9.raw = AB.b9.raw, 
       AB.b9.fdr = AB.b9.fdr,
       AB.b1.raw = AB.b1.raw, 
       AB.b1.fdr = AB.b1.fdr)
}

# ----------------------------------------------
# APDINTER SIMULATION WITH BOOT1 METHOD ONLY
apdinter_wrap <- function(x, n, w, nboot, qseq){
  res <- apdinter_hdqt7_cpp(x, n, w, nboot, qseq)
  # QT7 -----------------------
  bootdiff <- res[["bootABqt"]]
  # Compute percentile bootstrap p values:
  temp <- colSums(bootdiff<0) / nboot + colSums(bootdiff==0) / (2*nboot)
  qt.AB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # HD -----------------------
  bootdiff <- res[["bootABhd"]]
  # Compute percentile bootstrap p values:
  temp <- colSums(bootdiff<0) / nboot + colSums(bootdiff==0) / (2*nboot)
  hd.AB.pv <- 2*(apply(rbind(temp,1-temp), 2, min))
  # Return results ------------
  list(qt.AB = qt.AB.pv, hd.AB = hd.AB.pv)
}

apdinter_sim <- function(rdist, ...){
  
  print("apdinter: false positive simulation...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))

  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq*nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(rdist(nmax * 4), nrow = nmax)
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      w <- hd.w[[N]]
      pv <- apdinter_wrap(x, nx, w, nboot, qseq)
      AB.qt[S,N,] <- pv$qt.AB
      AB.hd[S,N,] <- pv$hd.AB
    }
  }
  
  beep(8)
  
  list(AB.qt = AB.qt, AB.hd = AB.hd)
}

# used for false positive and power simulations
apdinter_pois_sim <- function(lambda, qseq, ...){
  
  print("apdinter: Poisson simulation...")
 
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  ANOVA.m <- array(0, dim = c(nsim, nn))
  ANOVA.tm <- array(0, dim = c(nsim, nn))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq*nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(0, nrow = nmax, ncol = 4)
    for(C in 1:4){
      x.all[,C] <- rpois(nmax, lambda[C])
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      w <- hd.w[[N]]
      pv <- apdinter_wrap(x, nx, w, nboot, qseq)
      AB.qt[S,N,] <- pv$qt.AB
      AB.hd[S,N,] <- pv$hd.AB
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      ANOVA.m[S,N] <- t2way(2,2,xLIST,tr=0)$AB.p.value
      ANOVA.tm[S,N] <- t2way(2,2,xLIST,tr=0.2)$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.qt = AB.qt, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# used for false positive and power simulations
apdinter_bb_sim <- function(rseq, s, nbin, qseq, nseq, ...){
  
  print("apdinter: beta-binomial simulation...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  ANOVA.m <- array(0, dim = c(nsim, nn))
  ANOVA.tm <- array(0, dim = c(nsim, nn))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq*nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(0, nrow = nmax, ncol = 4)
    for(C in 1:4){ 
      x.all[,C] <- rbbinom(nmax, nbin = nbin, r = rseq[C], s = s)
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      w <- hd.w[[N]]
      pv <- apdinter_wrap(x, nx, w, nboot, qseq)
      AB.qt[S,N,] <- pv$qt.AB
      AB.hd[S,N,] <- pv$hd.AB
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      ANOVA.m[S,N] <- t2way(2,2,xLIST,tr=0)$AB.p.value
      ANOVA.tm[S,N] <- t2way(2,2,xLIST,tr=0.2)$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.qt = AB.qt, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# Compute FWER using FDR and Hochberg corrections
apdinter_fwer <- function(simres, ...){
  
  print("apdinter: compute fwer...")
  
  # interaction
  AB.qt.raw <- array(0, dim = c(nsim, nn))
  AB.qt.fdr <- array(0, dim = c(nsim, nn))
  AB.qt.hoch <- array(0, dim = c(nsim, nn))
  
  AB.hd.raw <- array(0, dim = c(nsim, nn))
  AB.hd.fdr <- array(0, dim = c(nsim, nn))
  AB.hd.hoch <- array(0, dim = c(nsim, nn))
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    for(N in 1:nn){
      # QT7
      p.val <- simres$AB.qt[S,N,]
      AB.qt.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.qt.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      AB.qt.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
      # HD
      p.val <- simres$AB.hd[S,N,]
      AB.hd.raw[S,N] <- sum(p.val < alpha.val) >= 1
      AB.hd.fdr[S,N] <- sum(p.adjust(p.val, method = "fdr") < alpha.val) >= 1
      AB.hd.hoch[S,N] <- sum(p.adjust(p.val, method = "hochberg") < alpha.val) >= 1
    }
  }
  list(AB.qt.raw = AB.qt.raw, 
       AB.qt.fdr = AB.qt.fdr,
       AB.qt.hoch = AB.qt.hoch,
       AB.hd.raw = AB.hd.raw, 
       AB.hd.fdr = AB.hd.fdr,
       AB.hd.hoch = AB.hd.hoch)
}

# Power simulation for differences in means
# HD + QT7
# Includes ANOVAs on means and 20% trimmed means for comparison
sim_apd_mean_diff <- function(rdist, gpmeans, ...){
  
  print("apdinter power simulation: mean differences...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  ANOVA.m <- array(0, dim = c(nsim, nn))
  ANOVA.tm <- array(0, dim = c(nsim, nn))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq*nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(rdist(nmax * 4), nrow = nmax)
    for(C in 1:4){ # add constants
      x.all[,C] <- x.all[,C] + gpmeans[C]  
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      # x <- matrix(rnorm(nx * 4), nrow = nx)
      pv <- apdinter_wrap(x, nx, hd.w[[N]], nboot, qseq)
      AB.qt[S,N,] <- pv$qt.AB
      AB.hd[S,N,] <- pv$hd.AB
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      ANOVA.m[S,N] <- t2way(2,2,xLIST,tr=0)$AB.p.value
      ANOVA.tm[S,N] <- t2way(2,2,xLIST,tr=0.2)$AB.p.value
    }
  }
  
  beep(8)

  list(AB.qt = AB.qt, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}

# Power simulation for differences in g parameter of g-and-h distributions
# Includes ANOVAs on means and 20% trimmed means for comparison
sim_apd_gh_diff <- function(gvals, hvals, ...){
  
  print("apdinter power simulation: g differences...")
  
  # interaction
  AB.qt <- array(0, dim = c(nsim, nn, 9))
  AB.hd <- array(0, dim = c(nsim, nn, 9))
  ANOVA.m <- array(0, dim = c(nsim, nn))
  ANOVA.tm <- array(0, dim = c(nsim, nn))
  
  # pre-compute the beta weights
  hd.w <- hd.w.calc(nseq*nseq, qseq)
  
  for(S in 1:nsim){
    
    sim.counter(S, nsim, inc = inc.step)
    
    # generate max size data
    x.all <- matrix(NA, nrow = nmax, ncol = 4)
    for(C in 1:4){ # add constants
      x.all[,C] <- ghdist(nmax, gvals[C], hvals[C])
    }
    
    for(N in 1:nn){
      # subsample
      nx <- nseq[N]
      x <- x.all[1:nx,]
      pv <- apdinter_wrap(x, nx, hd.w[[N]], nboot, qseq)
      AB.qt[S,N,] <- pv$qt.AB
      AB.hd[S,N,] <- pv$hd.AB
      
      xLIST <- list()
      for(C in 1:4){
        xLIST[[C]] <- x[,C] 
      }
      ANOVA.m[S,N] <- t2way(2,2,xLIST,tr=0)$AB.p.value
      ANOVA.tm[S,N] <- t2way(2,2,xLIST,tr=0.2)$AB.p.value
    }
  }
  
  beep(8)
  
  list(AB.qt = AB.qt, AB.hd = AB.hd, 
       ANOVA.m = ANOVA.m, ANOVA.tm = ANOVA.tm)
}
# ------------------------------------------------------


