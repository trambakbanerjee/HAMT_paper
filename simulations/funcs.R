library(irlba)

falt.nest=function(z,x,sig,weight,jacknife=T){
  #weighted NEST density function evaluated at z
  #In NEST the weight is just rep(1,length(x))
  #x are the observations (a vector)
  #sig are the corresponding standard deviations (a vector)
  #z is a 2-tuple,the first is the observation, the second is sd
  
  xds=density(x,from = min(x)-IQR(x)/2,to=max(x)+IQR(x)/2,n=length(x)/5)   #this line and the next are used to find bandwidth
  sigds=density(sig,from = min(sig)-IQR(sig)/2,to=max(sig)+IQR(sig)/2,n=length(x)/5)
  h=xds$bw
  hsig=sigds$bw
  if(jacknife==T){
    if(is.element(z[1],x)){
      z.ind=which(x==z[1])
      x=x[-z.ind[1]]
      sig=sig[-z.ind[1]]
      weight=weight[-z.ind[1]]/sum(weight[-z.ind[1]])
    }
  }
  obs=rep(z[1],length(x))
  obssig=rep(z[2],length(x))
  k=dnorm((obs-x)/h,0,sig)
  
  nestw=dnorm(obssig-sig,0,hsig)
  w=weight*nestw
  
  y=sum(w*k)/(sum(w)*h)
  return(y)
  
}


nest.func=function(x,sig,jacknife=T){
  #NEST procedure
  #Output is the density estimated at (x,sig) x and sig are both vectors 
  fm=function(a,b){
    res=falt.nest(c(a,b),x,sig,rep(1/length(x),length(x)),jacknife=jacknife)
    return(res)
  }
  denom=mapply(fm, x,sig)
  result=denom
  return(result)
}

nest.func.fast=function(x,sig,weight,n){
  
  xds=density(x,from = min(x)-IQR(x)/2,to=max(x)+IQR(x)/2,n=length(x)/5)   #this line and the next are used to find bandwidth
  sigds=density(sig,from = min(sig)-IQR(sig)/2,to=max(sig)+IQR(sig)/2,n=length(x)/5)
  h=xds$bw
  hsig=sigds$bw
  z=cbind(x,sig)
  y=matrix(0,n)
  for(i in 1:n){
    if(is.element(z[i,1],x)){
      z.ind=which(x==z[i,1])
      xx=x[-z.ind[1]]
      sigg=sig[-z.ind[1]]
      weightt=weight[-z.ind[1]]/sum(weight[-z.ind[1]])
    }
    obs=rep(z[i,1],length(xx))
    obssig=rep(z[i,2],length(xx))
    k=dnorm((obs-xx)/h,0,sigg)
    
    nestw=dnorm(obssig-sigg,0,hsig)
    w=weightt*nestw
    
    y[i]=sum(w*k)/(sum(w)*h)
  }
  return(y)
}

# Compute a truncated SVD approximation U*V' to rectangular matrix X,
# such that X is an n x m matrix, U is an n x k matrix, and V is an m
# x k matrix, where k <= min(n,m). The rank, k, is determined by the
# number of singular values surpassing the tolerance, "tol".
#
# Note that this function will only work if min(dim(X)) is 2 or
# greater.
#
# If irlba fails, the return value is NULL.
#
tsvd <- function (X, tol) {
  r <- min(dim(X))
  k <- 2
  
  # Iteratively increase the number of singular vectors in the SVD
  # until it is accurate enough (based on the tolerance setting,
  # "tol"), or until we hit an upper limit.
  while (TRUE) {
    # out <- tryCatch(irlba(X,k,tol = 1,svtol = 0.01*tol),
    #                 error   = function (e) NULL,
    #                 warning = function (e) NULL)
    out <- tryCatch(irlba(X,k,tol = 1,svtol = 0.01*tol),
                    error   = function (e) NULL)
    if (is.null(out))
      return(NULL)
    else if (k == r)
      break
    else if (min(out$d) < tol)
      break
    else
      k <- min(k+5,r)#min(2*k,r)
  }
  
  # Get the truncated SVD.
  i <- which(out$d > tol)
  if (length(i) < 2)
    i <- 1:2
  d <- out$d[i]
  U <- out$u[,i]
  V <- out$v[,i]
  U <- scale.cols(U,sqrt(d))
  V <- scale.cols(V,sqrt(d))
  return(list(U = U,V = V))
}

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b){
  t(t(A) * b)
}

#----------------- Functions for multiple testing --------------

bh.func<-function(pv, q)
{ 
  # the input 
  # pv: the p-values
  # q: the FDR level
  # the output 
  # nr: the number of hypothesis to be rejected
  # th: the p-value threshold
  # re: the index of rejected hypotheses
  # ac: the index of accepted hypotheses
  # de: the decision rule
  
  m=length(pv)
  st.pv<-sort(pv)   
  #print(length(st.pv))
  #print(m)
  pvi<-st.pv/1:m
  hps<-rep(0, m)
  if (max(pvi<=(q/m))==0)
  {
    k<-0
    pk<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k<-max(which(pvi<=(q/m)))
    pk<-st.pv[k]
    reject<-which(pv<=pk)
    accept<-which(pv>pk)
    hps[reject]<-1
  }
  y<-list(nr=k, th=pk, re=reject, ac=accept, de=hps)
  return (y)
}

sc.func<-function(lfdr, q)
{
  
  ## USAGE
  # mt.sc(lfdr, q)
  
  ## ARGUMENTS
  # lfdr: local false discovery rate sequence
  # q: the FDR level
  
  ## VALUES
  # nr: the number of rejected hypotheses
  # th: the threshold
  # re: the rejected hypotheses
  # ac: the accepted hypotheses
  # de: the decision rule
  
  m=length(lfdr)
  st.lfdr<-sort(lfdr)
  hps<-rep(0, m)
  #if (min(lfdr)>q)
  if (lfdr[which.min(lfdr)]>q)
  {
    k<-0
    threshold<-1
    reject<-NULL
    accept<-1:m
  }
  else
  {
    k=1
    while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
      k=k+1
    }
    k<-k-1
    threshold<-st.lfdr[k]
    reject<-which(lfdr<=threshold)
    accept<-which(lfdr>threshold)
    hps[reject]<-1
  }
  y<-list(nr=k, th=threshold, re=reject, ac=accept, de=hps)
  return (y)
}

clfdr.func=function(x,sig,mu0,w,gd,f_nest,type){
  
  n = length(sig)
  if(type=='npmle'){
    
    indm=max(which(gd<=mu0))
    gd.null=gd[1:indm]
    w.null=w[1:indm]
    # fnullest=function(x,s){
    #   #density at x for the null function
    #   res=dnorm(gd.null-x,0,s)
    #   return(res)
    # }
    # A=mapply(fnullest, x,sig)
    # fnull.dens=t(A)%*%w.null
    fnull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd.null-x[i],0,sig[i])*w.null)))
    ffull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd-x[i],0,sig[i])*w)))
    clfdr=pmin(fnull.dens/ffull.dens,1)
    clfdr=pmax(clfdr,0)
  }
  else if(type=='proposed'){
    
    indm=max(which(gd<=mu0))
    gd.null=gd[1:indm]
    w.null=w[1:indm,]
    # fnullest=function(x,s){
    #   #density at x for the null function
    #   res=dnorm(gd.null-x,0,s)
    #   return(res)
    # }
    # AA=mapply(fnullest, x,sig)
    # fnull.dens=sapply(1:n,function(i)t(AA[,i])%*%w.null[,i])
    fnull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd.null-x[i],0,sig[i])*w.null[,i])))
    ffull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd-x[i],0,sig[i])*w[,i])))
    clfdr=pmin(fnull.dens/ffull.dens,1)
    clfdr=pmax(clfdr,0)
    
  }
  return(clfdr)
}

clfdr.func.2=function(x,sig,mu0,mu1,w,gd,f_nest,type){
  
  n = length(sig)
  if(type=='npmle'){
    
    indm1=min(which(gd>=mu0))
    indm2=max(which(gd<=mu1))
    gd.null=gd[indm1:indm2]
    w.null=w[indm1:indm2]
    # fnullest=function(x,s){
    #   #density at x for the null function
    #   res=dnorm(gd.null-x,0,s)
    #   return(res)
    # }
    # A=mapply(fnullest, x,sig)
    # fnull.dens=t(A)%*%w.null
    fnull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd.null-x[i],0,sig[i])*w.null)))
    ffull.dens = unlist(sapply(1:n,function(i) sum(dnorm(gd-x[i],0,sig[i])*w)))
    clfdr=pmin(fnull.dens/ffull.dens,1)
    clfdr=pmax(clfdr,0)
  }
  else if(type=='proposed'){
    
    indm1=min(which(gd>=mu0))
    indm2=max(which(gd<=mu1))
    gd.null=gd[indm1:indm2]
    w.null=w[indm1:indm2,]
    # fnullest=function(x,s){
    #   #density at x for the null function
    #   res=dnorm(gd.null-x,0,s)
    #   return(res)
    # }
    # AA=mapply(fnullest, x,sig)
    # fnull.dens=sapply(1:n,function(i)t(AA[,i])%*%w.null[,i])
    fnull.dens = sapply(1:n,function(i) sum(dnorm(gd.null-x[i],0,sig[i])*w.null[,i]))
    ffull.dens = sapply(1:n,function(i) sum(dnorm(gd-x[i],0,sig[i])*w[,i]))
    clfdr=pmin(fnull.dens/ffull.dens,1)
    clfdr=pmax(clfdr,0)
    
  }
  return(clfdr)
}


#------------------------ Other functions ----------------
library("Rmosek")
library(Matrix)
ghat_knownsigma = function(f_pilot,A,B,n,m,K){
  
  f_pilot = matrix(f_pilot,n,1)
  
  #sparse matrix construction: A_sparse
  A_sparse = sparseMatrix(i=rep(1:n,each=m),j=1:(n*m),x=c(t(A)))
  cc = as.matrix(-2*t(A_sparse)%*%f_pilot)
  
  #sparse matrix construction: One_sparse
  One_sparse = sparseMatrix(i=rep(1:n,each=m),j=1:(n*m),x=rep(1,n*m))
  
  #sparse matrix construction: B_sparse
  BB = c(t(B))
  rr = c(sapply(1:n, function(i) rep((1+(i-1)*K):(i*K),m)))
  B_sparse = sparseMatrix(i=rep(1:(n*m),each = K),j=rep(1:(m*K),n),x=BB[rr])
  
  #------------------Mosek declarations-------------------------------
  #-------------------------------------------------------------------
  
  prob<-list(sense="min")
  
  #linear constraints in [u;w;t]
  prob$c = c(cc,rep(0,m*K),1)
  AA = rbind(cbind(One_sparse,
                   sparseMatrix(i=integer(0),j=integer(0),dims = list(n,m*K+1))),
             cbind(sparseMatrix(i=1:(n*m),j=1:(n*m),x = rep(1,n*m)),
                   -B_sparse,sparseMatrix(i=integer(0),j=integer(0),dims=list(n*m,1))))
  prob$A<- AA
  prob$bc<- rbind(blc=c(rep(1,n),rep(0,n*m)),
                  buc = c(rep(1,n),rep(0,n*m)))
  prob$bx<- rbind(blx=c(rep(0,n*m),rep(-Inf,m*K+1)),
                  bux = rep(Inf,n*m+m*K+1))
  
  #Affine conic constraint
  prob$F<- rbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,n*m+m*K+1)),
                 cbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,n*m+m*K)),1),
                 cbind(A_sparse,sparseMatrix(i=integer(0),j=integer(0),dims=list(n,m*K)),
                       sparseMatrix(i=integer(0),j=integer(0),dims=list(n,1))))
  # prob$F<- rbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,n*m+m*K+1)),
  #                cbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,n*m+m*K)),1),
  #                cbind(A_sparse,sparseMatrix(i=integer(0),j=integer(0),dims=list(n,m*K)),0))
  
  
  prob$g<-c(0.5,rep(0,n+1))
  prob$cones<- matrix(list("RQUAD",n+2,NULL))
  rownames(prob$cones)<- c("type","dim","conepar")
  prob$iparam <- list(PRESOLVE_USE = "PRESOLVE_MODE_OFF")
  
  out = mosek(prob)#mosek(prob,list(verbose=1))
  estimate = out$sol$itr$xx
  theta_mosek = matrix(estimate[1:(n*m)],ncol=m,byrow = T) 
  theta_mosek_clean = matrix(0,n,m)
  for(i in 1:n){
    theta_mosek_clean[i,] = theta_mosek[i,]/sum(theta_mosek[i,])
  }
  return(list('theta_hat'=theta_mosek_clean,'msg'=out$sol$itr$solsta))
  
}

hhat = function(f_s2_pilot,A,n,q){
  
  f_s2_pilot = matrix(f_s2_pilot,n,1)
  
  #sparse matrix construction: A_sparse
  A_sparse = as(A,"sparseMatrix")
  cc = -2*t(A)%*%f_s2_pilot
  
  #------------------Mosek declarations-------------------------------
  #-------------------------------------------------------------------
  
  prob<-list(sense="min")
  
  #linear constraints in [v;t]
  prob$c = c(cc,1)
  prob$A<- as(cbind(matrix(rep(1,q),1,q),0),"sparseMatrix")
  prob$bc<- rbind(blc=1,buc = 1)
  prob$bx<- rbind(blx=c(rep(0,q),-Inf),
                  bux = rep(Inf,q+1))
  
  #Affine conic constraint
  prob$F<- rbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,q+1)),
                 cbind(sparseMatrix(i=integer(0),j=integer(0),dims=list(1,q)),1),
                 cbind(A_sparse,sparseMatrix(i=integer(0),j=integer(0),dims=list(n,1))))
  
  prob$g<-c(0.5,rep(0,n+1))
  prob$cones<- matrix(list("RQUAD",n+2,NULL))
  rownames(prob$cones)<- c("type","dim","conepar")
  
  out = mosek(prob,list(verbose=1))
  estimate = out$sol$itr$xx
  v_mosek = matrix(estimate[1:q],nrow=q,ncol=1) 
  return(list('v_hat'=v_mosek,'msg'=out$sol$itr$solsta))
  
}

BDE <- function(y, sig,Tg = 300, df = 5, c0 = 0.1){
  # Bayesian Deconvolution Estimator: Efron (B'ka, 2016)
  require(splines)
  eps <- 1e-04
  if(length(Tg) == 1) Tg <- seq(min(y)-eps, max(y)+eps, length = Tg)
  X <- ns(Tg, df = df)
  a0 <- rep(0, ncol(X))
  A <- dnorm(outer(y,Tg,"-"), sd = sig)
  qmle <- function(a, X, A, c0){
    g <- exp(X %*% a)
    g <- g/sum(g)
    f <- A %*% g
    -sum(log(f)) + c0 * sum(a^2)^.5
  }
  ahat <- nlm(qmle, a0, X=X, A=A, c0 = c0)$estimate
  g <- exp(X %*% ahat)
  #g <- g/integrate(approxfun(T,g),min(T),max(T))$value
  g <- c(g/sum(g * diff(Tg)[1]))
  z <- list(x = Tg,y = g, sigma = sig)
  class(z) <- "GLmix"
  z
}

#----------------------------- Functions from Gu and Shen EJ ----------------
pwfwBH <- function(t, q = 0.1){ 
  pval<-1-pnorm(t)  
  Ng<-length(pval)
  gind=1:Ng
  rejpw=gind[pval<q]
  temp=data.frame(pval,gind)
  temp=temp[order(pval),]
  rejfw=sort(temp$gind[temp$pval<q/(Ng+1-seq(1:Ng))])
  rejBH=sort(temp$gind[temp$pval<seq(1:Ng)/Ng*q])	  #BH (1995paper)
  grank=temp$gind
  list(rejpw = rejpw, rejfw = rejfw, rejBH=rejBH,grank=grank)	
}
BH2 <- function(t, omega, q = 0.1){
  # modified BH procedure [based on Benjamini et al 2006paper]
  # omega is the non-null proportion 
  pval<-1-pnorm(t)  
  Ng<-length(pval)
  gind=1:Ng
  temp=data.frame(pval,gind)
  temp=temp[order(pval),]
  qstar = q/(1-omega)
  rejBH=sort(temp$gind[temp$pval<seq(1:Ng)/Ng*qstar])
  list(rejBH = rejBH)
}
Lfdr <- function(x, sigma, Acut = c(-1,1), g = NULL, method = "OR")
{
  if(method %in% c("OR","SM")){ # Oracle Lfdr Statistics (see Thm 1 of Sun-McLain), 
    T <- x
    if(method == "SM") { 
      g <- dNonNull(x, sigma)
      omega <- g$omega
      g$y[c(1, length(g$y))] <- 0 # Kill the tails
      g <- approxfun(g$x, g$y, rule = 2)
    }
    else {
      omega <- g$omega
      g <- g$g
    }
    for(i in 1:length(x)){ # NB (1-omegabar) cancels.
      h <- function(u) dnorm(x[i] - u, sd = sigma[i])*g(u)
      num <- ((1 - omega) * dnorm(x[i],sd = sigma[i]) + 
                omega * integrate(h,Acut[1],Acut[2], stop.on.error = FALSE)$value)
      denom <- (1 - omega) * dnorm(x[i],sd = sigma[i]) + 
        omega * integrate(h,-Inf, Inf, stop.on.error = FALSE)$value
      T[i] <- num/denom
    }
  }
  else if(method == "DE"){ # Simple Deconvolution Lfdr Statistics
    bw <- bw.dmise(x, sigma, error = "laplacian", kernel = "normal")
    f <- DeconPdf(x, sigma, bw = bw)
    dv <- diff(f$x)[1]
    A <- (dnorm(outer(x,f$x,"-"), sd = sigma) * dv * outer(rep(1,length(x)),f$y))
    T <- apply(A[, ((f$x > Acut[1]) & (f$x < Acut[2]))], 1, sum)/apply(A, 1, sum)
  }
  else if(method == "KW"){ # Kiefer-Wolfowitz Lfdr Statistics 
    f <- GLmix(x, sigma = sigma, v = 500, rtol = 1e-10)
    dv <- diff(f$x)[1]
    A <- (dnorm(outer(x,f$x,"-"), sd = sigma) * dv * outer(rep(1,length(x)),f$y))
    T <- apply(A[, ((f$x > Acut[1]) & (f$x < Acut[2]))], 1, sum)/apply(A, 1, sum)
  }
  else if(method == "KW.zero"){ # Kiefer-Wolfowitz Lfdr Statistics with a mass at zero
    eps = epsest.func(x/sigma,0,1)
    f <- GLmix.disc(x, omega =eps, sigma = sigma, v = 500, rtol = 1e-10)
    dv <- diff(f$x)[1]
    A <- eps * dnorm(outer(x,f$x,"-"), sd = sigma) * dv * outer(rep(1,length(x)),f$y)
    num <- ((1 - eps) * dnorm(x,sd = sigma) + 
              apply(A[,((f$x > Acut[1]) & (f$x < Acut[2]))], 1, sum))
    denom <- (1 - eps) * dnorm(x,sd = sigma) + 
      apply(A,1,sum)
    T <- num/denom
  }
  T
}
mFDR <- function(T, q = 0.1){ 
  # Adaptive step-up procedure Sun and Cai eq (10)
  # T is the vector of lfdr test statistics, q the chosen level
  cummean <- function(x) cumsum(x)/(1:length(x))
  S <- sort(T)
  k <- sum(cummean(S) < q)
  cval <- S[k]
  FDR <- mean((S < cval)*S)/mean(S < cval)
  list(FDR = FDR, cval = cval)
}
# Realized FDR and FNR for simulations
rFDR <- function(T, sample, Acut = c(-1,1), qFDR = 0.1){
  K <- 1:length(T)
  X <- K[T < mFDR(T, q = qFDR)$cval]
  O <- K[((sample$mu < Acut[1]) | (sample$mu > Acut[2]))]
  length(setdiff(X,O))/length(X)
}
rFNR <- function(T, sample, Acut = c(-1,1), qFDR = 0.1){
  K <- 1:length(T)
  X <- K[T < mFDR(T, q = qFDR)$cval]
  O <- K[((sample$mu < Acut[1]) | (sample$mu > Acut[2]))]
  length(setdiff(O,X))/(length(sample$x) - length(X))
}
rcover <- function(T, sample, Acut = c(-1,1), qFDR = 0.1){
  # computes coverage probability
  K <- 1:length(T)
  X <- K[T < mFDR(T, q = qFDR)$cval]
  O <- K[((sample$mu < Acut[1]) | (sample$mu > Acut[2]))]
  sum(X %in% O)/length(O)
}

epsest.func <- function(x,u,sigma)
{
  # x is a vector
  # u is the mean
  # sigma is the standard deviation
  
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.1)
  
  epsest=NULL
  
  for (j in 1:length(tt)) { 
    
    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi
    
    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}
GLmix.disc <- 
  function (x, omega, v = 300, sigma = 1, hist = FALSE, histm = 300, weights = NULL, 
            ...) 
  {
    n <- length(x)
    eps <- 1e-04
    if (length(v) == 1) 
      v <- seq(min(x) - eps, max(x) + eps, length = v)
    m <- length(v)
    weights <- weights/sum(weights)
    w <- weights
    if (hist) {
      histbin <- function(x, m = histm, eps = 1e-06, weights = weights) {
        u <- seq(min(x) - eps, max(x) + eps, length = m)
        xu <- findInterval(x, u)
        txu <- tabulate(xu)
        midu <- (u[-1] + u[-m])/2
        wnz <- (txu > 0)
        if (length(weights)) {
          if (length(weights) == length(x)) 
            w <- as.numeric(tapply(weights, xu, sum))
          else stop("length(weights) not equal length(x)")
        }
        else w <- txu[wnz]/sum(txu[wnz])
        list(midu = midu[wnz], w = w)
      }
      if (length(sigma) == 1) {
        h <- histbin(x, m = histm, weights = weights)
        x <- h$midu
        w <- h$w
      }
      else {
        sus <- sort(unique(sigma))
        us <- match(sigma, sus)
        nus <- table(us)
        if (min(nus) < 100) 
          stop("too few obs in some sigma bin")
        h <- as.list(1:length(nus))
        for (i in 1:length(sus)) {
          if (length(weights)) 
            h[[i]] <- histbin(x[us == i], m = histm, weights = weights[us == 
                                                                         i])
          else h[[i]] <- histbin(x[us == i], m = histm, 
                                 weights = weights)
        }
        x <- unlist(lapply(h, function(f) f$midu))
        w <- unlist(lapply(h, function(f) f$w))
        w <- w/sum(w)
        sigma <- rep(sus, unlist(lapply(h, function(f) length(f$midu))))
      }
    }
    if (!length(w)) 
      w <- rep(1, n)/n
    d <- diff(v)
    v <- (v[-1] + v[-m])/2
    A <- (1-omega) * dnorm(x, sd = sigma) + omega * dnorm(outer(x, v, "-"), sd = sigma)
    f <- KWDual(A, d, w, ...)
    y <- f$f
    g <- f$g
    logLik <- n * sum(w * log(g))
    dy <- as.vector((A %*% (y * d * v))/g)
    z <- list(x = v, y = y, g = g, logLik = logLik, sigma = sigma, 
              dx = x, dy = dy, status = f$status)
    class(z) <- c("GLmix", "density")
    return(z)
  }
