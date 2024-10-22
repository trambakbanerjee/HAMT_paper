
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(ashr)
library(Rmosek)
library(REBayes)

exp_1<-function(n,mu0,mu1,u,q,reps){
  
  fdr.ash<-matrix(0,reps,length(u))
  fdr.or<-matrix(0,reps,length(u))
  etp.ash<-matrix(0,reps,length(u))
  fdr.ash.1<-matrix(0,reps,length(u))
  etp.ash.1<-matrix(0,reps,length(u))
  fdr.npmleb<-matrix(0,reps,length(u))
  etp.npmleb<-matrix(0,reps,length(u))
  fdr.dd<-matrix(0,reps,length(u))
  etp.dd<-matrix(0,reps,length(u))
  etp.or<-matrix(0,reps,length(u))
  fdr.npmle<-matrix(0,reps,length(u))
  fdr.deconv<-matrix(0,reps,length(u))
  etp.npmle<-matrix(0,reps,length(u))
  etp.deconv<-matrix(0,reps,length(u))
  fdr.adapt<-matrix(0,reps,length(u))
  etp.adapt<-matrix(0,reps,length(u))
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("ashr","Rmosek","REBayes"))%dopar%{
      source('../funcs.R')
      set.seed(r)
      sig = sample(c(0.25,0.75,1.5),n,replace = TRUE)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+(-u[nn]*sig[i])*(p[i]>0.9 & p[i]<=0.95)+
                    (u[nn]*sig[i])*(p[i]>0.95))
      x=sapply(1:n, function(i) mu[i]+rlogis(1,0,sig[i]))
      theta=1*(mu>mu1)+1*(mu<mu0)
      sig_x = sig*pi/sqrt(3)
      
      ################## Oracle Clfdr ###############################
      clfdr.or = rep(0,n)
      for(i in 1:n){
        denom = 0.9*dlogis(x[i],0,sig[i])+
          0.05*dlogis(x[i],-u[nn]*sig[i],sig[i])+
          0.05*dlogis(x[i],u[nn]*sig[i],sig[i])
        numer = 0.9*dlogis(x[i],0,sig[i])+
          (u[nn]*sig[i]<=mu1 & -u[nn]*sig[i]>=mu0)*(0.05*dlogis(x[i],u[nn]*sig[i],sig[i])+
                                                      0.05*dlogis(x[i],-u[nn]*sig[i],sig[i]))
        clfdr.or[i]=numer/denom
      }
      decisions.or = sc.func(clfdr.or,q)
      fdp.or.r = sum((1-theta)*decisions.or$de)/max(sum(decisions.or$de),1)
      etp.or.r=sum((theta)*decisions.or$de)/max(sum(theta), 1)
      
      ######################### AdaPTGMM ################################
      fdp.adapt.r = 0
      etp.adapt.r=0
      
      #### NPMLE based Clfdr ###################################
      fdp.npmle.r = 0
      etp.npmle.r=0
      
      #### Efron's Deconv based Clfdr ###################################
      fdp.deconv.r = 0
      etp.deconv.r= 0
      ########################## alpha = 0 ash ##############################
      fdp.ash.r = 0
      etp.ash.r= 0
      
      ####################################### alpha = 1 ##################
      fdp.ash.1.r = 0
      etp.ash.1.r= 0
      ################# NPMLE with basis expansion ####################
      m = 50
      K = 10
      gd=seq(from=min(x),to=max(x),length.out=m)
      npmleb = npmleb_opt(x,sig,gd,length_out = m,bases = K,likelihood = 'logistic')
      g_npmle = npmleb$probs
      
      clfdr.npmleb = clfdr.Logisfunc.2(x,sig,mu0,mu1,g_npmle,gd,
                                   type='proposed')
      decisions.npmleb = sc.func(clfdr.npmleb,q)
      fdp.npmleb.r = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
      etp.npmleb.r=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
      
      ################# Proposed ####################
      f_hat = nest.func.fast(x,sig_x,rep(1/length(x),length(x)),n)
      hamt = hamt_opt(x,sig,f_hat,gd,length_out = m,bases = K,likelihood = 'logistic')
      w_mosek = hamt$probs
      
      clfdr.dd = clfdr.Logisfunc.2(x,sig,mu0,mu1,w_mosek,gd,
                               type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      return(list("fdr.or_r"=fdp.or.r,
                  "fdr.npmle_r"=fdp.npmle.r,
                  "fdr.deconv_r"=fdp.deconv.r,
                  "fdr.ash_r"=fdp.ash.r,
                  "fdr.dd_r"=fdp.dd.r,
                  "fdr.ash.1_r"=fdp.ash.1.r,
                  "fdr.npmleb_r"=fdp.npmleb.r,
                  "etp.ash_r"=etp.ash.r,
                  "etp.dd_r"=etp.dd.r,
                  "etp.ash.1_r"=etp.ash.1.r,
                  "etp.npmleb_r"=etp.npmleb.r,
                  "etp.npmle_r"=etp.npmle.r,
                  "etp.deconv_r"=etp.deconv.r,
                  "etp.or_r"=etp.or.r,
                  "fdr.adapt_r"=fdp.adapt.r,
                  "etp.adapt_r"=etp.adapt.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.or[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.or_r)
    fdr.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle_r)
    fdr.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.deconv_r)
    fdr.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash_r)
    etp.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash_r)
    fdr.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash.1_r)
    etp.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash.1_r)
    fdr.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmleb_r)
    etp.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmleb_r)
    fdr.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.dd_r)
    etp.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.dd_r)
    etp.or[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.or_r)
    etp.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle_r)
    etp.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.deconv_r)
    fdr.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.adapt_r)
    etp.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.adapt_r)
    print(paste(u[nn],' -- ',
                mean(sapply(1:reps,function(i) result[[i]]$fdr.dd_r)),sep=''))
  }
  result<-list("fdr.or"=fdr.or,
               "fdr.npmle"=fdr.npmle,
               "fdr.deconv"=fdr.deconv,
               "fdr.ash"=fdr.ash,
               "fdr.ash.1"=fdr.ash.1,
               "fdr.npmleb"=fdr.npmleb,
               "fdr.dd"=fdr.dd,
               "etp.ash"=etp.ash,
               "etp.ash.1"=etp.ash.1,
               "etp.npmleb"=etp.npmleb,
               "etp.dd"=etp.dd,
               "etp.npmle"=etp.npmle,
               "etp.deconv"=etp.deconv,
               "etp.or"=etp.or,
               "fdr.adapt"=fdr.adapt,
               "etp.adapt"=etp.adapt)
}

#--------------------------------------------------------------------

n<-10000
u = c(5.5,5.6,5.7,5.8,5.9,6)
mu0 = -2
mu1 = 2
q=0.1
reps=500

out_exp1<-exp_1(n,mu0,mu1,u,q,reps)

save.image('fig13_set2.RData')
