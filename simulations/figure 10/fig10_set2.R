
library(REBayes)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(rmutil)
library(Rmosek)

exp_1<-function(n,mu0,mu1,u,q,length_out,bases,reps){
  
  fdr.dd<-matrix(0,reps,length(u))
  fdr.npmle<-matrix(0,reps,length(u))
  fdr.bh<-matrix(0,reps,length(u))
  fdr.or<-matrix(0,reps,length(u))
  fdr.deconv<-matrix(0,reps,length(u))
  
  etp.dd<-matrix(0,reps,length(u))
  etp.npmle<-matrix(0,reps,length(u))
  etp.bh<-matrix(0,reps,length(u))
  etp.or<-matrix(0,reps,length(u))
  etp.deconv<-matrix(0,reps,length(u))
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("REBayes","MASS","Rmosek","rmutil"))%dopar%{
      source('../funcs.R')
      set.seed(r)
      sig = runif(n,0.5,u[nn])
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+rnorm(1,3,sig[i])*(p[i]>0.9 & p[i]<=0.95)+
                    rnorm(1,-3,sig[i])*(p[i]>=0.95))
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu1)+1*(mu<mu0)
      
      f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
      ################## Oracle Clfdr ###############################
      clfdr.or = rep(0,n)
      for(i in 1:n){
        ff<-function(u,x,s,mm){dnorm(x,u,s)*dnorm(u,mm,s)
        } 
        denom = 0.9*dnorm(x[i],0,sig[i])+0.05*integrate(ff,-Inf,Inf,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,-Inf,Inf,x=x[i],s=sig[i],mm=-3)$value
        numer = 0.9*dnorm(x[i],0,sig[i])+0.05*integrate(ff,-2,0,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,-2,0,x=x[i],s=sig[i],mm=-3)$value+
          0.05*integrate(ff,0,2,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,0,2,x=x[i],s=sig[i],mm=-3)$value
        clfdr.or[i]=numer/denom
      }
      decisions.or = sc.func(clfdr.or,q)
      fdp.or.r = sum((1-theta)*decisions.or$de)/max(sum(decisions.or$de),1)
      etp.or.r=sum((theta)*decisions.or$de)/max(sum(theta), 1)
      
      #### NPMLE based Clfdr ###################################
      npmle.g<- GLmix(x,sigma=sig)
      w.npmle=npmle.g$y
      gd.npmle=npmle.g$x
      clfdr.npmle = clfdr.func.2(x,sig,mu0,mu1,w.npmle,gd.npmle,
                                 type='npmle')
      decisions.npmle = sc.func(clfdr.npmle,q)
      fdp.npmle.r = sum((1-theta)*decisions.npmle$de)/max(sum(decisions.npmle$de),1)
      etp.npmle.r=sum((theta)*decisions.npmle$de)/max(sum(theta), 1)
      
      #### Efron's Deconv based Clfdr ###################################
      deconv.g<- BDE(x,sig)
      w.deconv=deconv.g$y
      gd.deconv=deconv.g$x
      clfdr.deconv = clfdr.func.2(x,sig,mu0,mu1,w.deconv,gd.deconv,
                                  type='npmle')
      decisions.deconv = sc.func(clfdr.deconv,q)
      fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
      etp.deconv.r=sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
      
      ################# Proposed ####################
      m = length_out
      K = bases
      gd=seq(from=min(x),to=max(x),length.out=m)
      hamt = hamt_opt(x,sig,f_hat,gd,length_out = m,bases = K)
      w_mosek = hamt$probs
      
      clfdr.dd = clfdr.func.2(x,sig,mu0,mu1,w_mosek,gd,
                              type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      #################################### BH ###################################
      fdp.bh.r = 0
      etp.bh.r=0
      
      return(list("fdr.dd_r"=fdp.dd.r,"fdr.npmle_r"=fdp.npmle.r,"fdr.deconv_r"=fdp.deconv.r,
                  "fdr.bh_r"=fdp.bh.r,"fdr.or_r"=fdp.or.r,
                  "etp.dd_r"=etp.dd.r,"etp.npmle_r"=etp.npmle.r,"etp.deconv_r"=etp.deconv.r,
                  "etp.bh_r"=etp.bh.r,"etp.or_r"=etp.or.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.bh_r)
    fdr.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle_r)
    fdr.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.deconv_r)
    fdr.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.dd_r)
    fdr.or[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.or_r)
    etp.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.bh_r)
    etp.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle_r)
    etp.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.deconv_r)
    etp.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.dd_r)
    etp.or[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.or_r)
    print(paste(u[nn],' -- ',
                mean(sapply(1:reps,function(i) result[[i]]$fdr.dd_r)),sep=''))
  }
  result<-list("fdr.dd"=fdr.dd,"fdr.npmle"=fdr.npmle,"fdr.deconv"=fdr.deconv,
               "fdr.bh"=fdr.bh,"fdr.or"=fdr.or,
               "etp.dd"=etp.dd,"etp.npmle"=etp.npmle,"etp.deconv"=etp.deconv,
               "etp.bh"=etp.bh,"etp.or"=etp.or)
}

#--------------------------------------------------------------------

n<-10000
u = c(1,1.2,1.4,1.6,1.8,2)
length_out = 50
bases = 10
mu0 = -2
mu1 = 2
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,mu1,u,q,length_out, bases,reps)

save.image('fig10_set2.RData')
