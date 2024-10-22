
library(REBayes)
library(Rmosek)
library(splines)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)

exp_1<-function(n,mu0,u,q,length_out,bases,reps){
  
  fdr.dd<-matrix(0,reps,length(u))
  fdr.npmle1<-fdr.npmle2<-matrix(0,reps,length(u))
  fdr.deconv<-matrix(0,reps,length(u))
  fdr.bh<-matrix(0,reps,length(u))
  fdr.or<-matrix(0,reps,length(u))

  etp.dd<-matrix(0,reps,length(u))
  etp.npmle1<-etp.npmle2<-matrix(0,reps,length(u))
  etp.deconv<-matrix(0,reps,length(u))
  etp.bh<-matrix(0,reps,length(u))
  etp.or<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(6)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("REBayes","MASS","isotone","adaptMT","splines"))%dopar%{
      source('../funcs.R')
      set.seed(r)
      p = runif(n,0,1)
      sig = sapply(1:n, function(i) 0.5*(p[i]<=0.33)+2*(p[i]>0.67)+
                     1*(p[i]>0.33 & p[i]<=0.67))
      set.seed(r+1)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+
                    u[nn]*sig[i]*(p[i]>0.9))
      set.seed(r+2)
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      
      f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
      ################## Oracle Clfdr ###############################
      clfdr.or = rep(0,n)
      for(i in 1:n){
        clfdr.or[i]=(0.9*dnorm(x[i],0,sig[i])+(u[nn]*sig[i]<=2)*0.1*dnorm(x[i],(u[nn]*sig[i]),sig[i]))/(0.9*dnorm(x[i], 0, sig[i])+
                                                    0.1*dnorm(x[i],(u[nn]*sig[i]),sig[i]))
      }
      decisions.or = sc.func(clfdr.or,q)
      fdp.or.r = sum((1-theta)*decisions.or$de)/max(sum(decisions.or$de),1)
      etp.or.r=sum((theta)*decisions.or$de)/max(sum(theta), 1)
      
      #### NPMLE based Clfdr ###################################
      tkw = Lfdr((x-mu0)/sig, sigma = 1, Acut = c(-Inf,0), method = "KW")   # KW
      tkw2 = Lfdr((x-mu0)/sig,sigma = 1,  Acut = c(-Inf,0), method = "KW.zero")   #KW-Hybrid
      decisions.npmle1 = sc.func(tkw,q)
      fdp.npmle1.r = sum((1-theta)*decisions.npmle1$de)/max(sum(decisions.npmle1$de),1)
      etp.npmle1.r=sum((theta)*decisions.npmle1$de)/max(sum(theta), 1)
      decisions.npmle2 = sc.func(tkw2,q)
      fdp.npmle2.r = sum((1-theta)*decisions.npmle2$de)/max(sum(decisions.npmle2$de),1)
      etp.npmle2.r=sum((theta)*decisions.npmle2$de)/max(sum(theta), 1)
      
      #### Efron's Deconv based Clfdr ###################################
      deconv.g<- tryCatch({BDE(x,sig)},error=function(e) return(NaN))
      if(!is.list(deconv.g)){
        
        fdp.deconv.r = NA
        etp.deconv.r = NA
        
      } else {
        w.deconv=deconv.g$y/sum(deconv.g$y)
        gd.deconv=deconv.g$x
        clfdr.deconv = clfdr.func(x,sig,mu0,w.deconv,gd.deconv,
                                  type='npmle')
        decisions.deconv = sc.func(clfdr.deconv,q)
        fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
        etp.deconv.r=sum((theta)*decisions.deconv$de)/max(sum(theta), 1)}
      ################# Proposed ####################
      m = length_out
      K = bases
      gd=seq(from=min(x),to=max(x),length.out=m)
      hamt = hamt_opt(x,sig,f_hat,gd,length_out = m,bases=K)
      w_mosek = hamt$probs
      
      clfdr.dd = clfdr.func(x,sig,mu0,w_mosek,gd,
                            type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      #################################### BH ###################################
      z=(x-mu0)/sig
      rejBH2 = BH2(t= z, omega = sum(theta)/n)$rejBH  # advanced BH based on Benjamini et al (2006paper)
      bh.res=list('de'=rep(0,n))
      bh.res$de[rejBH2]=1
      fdp.bh.r = sum((1-theta)*bh.res$de)/max(sum(bh.res$de),1)
      etp.bh.r=sum((theta)*bh.res$de)/max(sum(theta), 1)
      
      return(list("fdr.dd_r"=fdp.dd.r,"fdr.npmle1_r"=fdp.npmle1.r,"fdr.npmle2_r"=fdp.npmle2.r,
                  "fdr.deconv_r"=fdp.deconv.r,
                  "fdr.bh_r"=fdp.bh.r,"fdr.or_r"=fdp.or.r,
                  "etp.dd_r"=etp.dd.r,"etp.npmle1_r"=etp.npmle1.r,"etp.npmle2_r"=etp.npmle2.r,
                  "etp.deconv_r"=etp.deconv.r,
                  "etp.bh_r"=etp.bh.r,"etp.or_r"=etp.or.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.bh_r)
    fdr.npmle1[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle1_r)
    fdr.npmle2[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle2_r)
    fdr.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.deconv_r)
    fdr.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.dd_r)
    fdr.or[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.or_r)
    etp.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.bh_r)
    etp.npmle1[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle1_r)
    etp.npmle2[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle2_r)
    etp.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.deconv_r)
    etp.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.dd_r)
    etp.or[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.or_r)
    print(paste(u[nn],' -- ',
                mean(sapply(1:reps,function(i) result[[i]]$fdr.dd_r)),sep=''))
  }
  result<-list("fdr.dd"=fdr.dd,"fdr.npmle1"=fdr.npmle1,"fdr.npmle2"=fdr.npmle2,"fdr.deconv"=fdr.deconv,
               "fdr.bh"=fdr.bh,"fdr.or"=fdr.or,
               "etp.dd"=etp.dd,"etp.npmle1"=etp.npmle1,"etp.npmle2"=etp.npmle2,"etp.deconv"=etp.deconv,
               "etp.bh"=etp.bh,"etp.or"=etp.or)
}

#--------------------------------------------------------------------

n<-10000
u = c(2.2,2.3,2.4,2.5,2.6,2.7)
length_out = 50
bases = 10
mu0= 2
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,u,q,length_out,bases,reps)

save.image('fig5_set2.RData')
