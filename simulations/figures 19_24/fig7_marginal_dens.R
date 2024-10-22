
library(Rmosek)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(Rfast)
library(ashr)

n<-10000
u = c(1.2,1.4,1.8)
mu0= 1
q=0.1
r = 42


#-----------------------------------------------------------------------------a


prior_hamt = prior_npmleb_1 = prior_npmleb_4 = array(0,c(50,n,length(u)))
prior_npmleb_2 = prior_npmleb_3 = array(0,c(100,n,length(u)))

fdp.npmleb = etp.npmleb = matrix(0,length(u),4)
fdp.dd = etp.dd = matrix(0,length(u),1)

sig_n = matrix(0,n,length(u))
gd_1 = gd_4 = matrix(0,50,length(u))
gd_2 = gd_3 = matrix(0,100,length(u))

for(nn in 1:length(u)){
  
  source('../funcs.R')
  set.seed(r)
  p = runif(n,0,1)
  sig = sapply(1:n, function(i) runif(1,0.5,1)*(p[i]<=0.9)+runif(1,1,u[nn])*(p[i]>0.9))
  mu = sapply(1:n, function(i) 0*(sig[i]<=1)+
                (2/sig[i])*(sig[i]>1))
  x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
  theta=1*(mu>mu0)
 
  f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
  
  gd=seq(from=min(x),to=max(x),length.out=50)
  
  hamt_w = hamt_opt(x,sig,f_hat,gd,length_out = 50,bases=10)
  
  clfdr.dd = clfdr.func(x,sig,mu0,hamt_w$probs,gd,type='proposed')
  decisions.dd = sc.func(clfdr.dd,q)
  fdp.dd[nn] = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
  etp.dd[nn]=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
  
  npmleb = npmleb_opt(x,sig,gd,length_out = 50,bases=10)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,1] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,1]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_hamt[,,nn] = hamt_w$probs
  prior_npmleb_1[,,nn] = npmleb$probs
  
  sig_n[,nn] = sig
  gd_1[,nn] = gd
  
  gd=seq(from=min(x),to=max(x),length.out=100)
  gd_2[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 100,bases=30)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,2] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,2]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_2[,,nn] = npmleb$probs
  
  gd=seq(from=min(x),to=max(x),length.out=100)
  gd_3[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 100,bases=10)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,3] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,3]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_3[,,nn] = npmleb$probs
  
  gd=seq(from=min(x),to=max(x),length.out=50)
  gd_4[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 50,bases=30)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,4] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,4]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_4[,,nn] = npmleb$probs
  
}

#--------------------------------------------------------------------

save.image('fig7_marginal_dens.RData')
