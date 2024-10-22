
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(ashr)
library(Rmosek)
library(REBayes)

exp_1<-function(n,mu0,u,q,reps){
  
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
  fdr.npmle1<-fdr.npmle2<-matrix(0,reps,length(u))
  fdr.deconv<-matrix(0,reps,length(u))
  etp.npmle1<-etp.npmle2<-matrix(0,reps,length(u))
  etp.deconv<-matrix(0,reps,length(u))
  fdr.adapt<-matrix(0,reps,length(u))
  etp.adapt<-matrix(0,reps,length(u))
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("ashr","Rmosek","REBayes"))%dopar%{
      source('../funcs.R')
      mm = rep(u[nn],n)
      set.seed(r)
      p = runif(n,0,1)
      sig = sqrt(mm)*sapply(1:n, function(i) runif(1,0.5,1)*(p[i]<=0.9)+runif(1,1,1.5)*(p[i]>0.9))
      mu = sapply(1:n, function(i) 0*(sig[i]/sqrt(mm[i])<=1)+
                    (2*sqrt(mm[i])/sig[i])*(sig[i]/sqrt(mm[i])>1))
      x=sapply(1:n, function(i) mu[i]+sig[i]*rt(mm[i],df=10,ncp = 0)*sqrt(8/10))
      theta=1*(mu>mu0)
      xbar = colMeans(x)
      S =sqrt(sapply(1:n, function(i) var(x[,i])/mm[i]))
      f_hat = nest.func.fast(xbar,S,rep(1/length(xbar),length(xbar)),n)
      idd = which(is.nan(f_hat))
      f_hat[idd] = 0
      
      ################## Oracle Clfdr ###############################
      fdp.or.r = 0
      etp.or.r = 1
      ######################### AdaPTGMM ################################
      fdp.adapt.r = 0
      etp.adapt.r=0
      
      #### NPMLE based Clfdr ###################################
      tkw = Lfdr((xbar-mu0)/S, sigma = 1, Acut = c(-Inf,0), method = "KW")   # KW
      tkw2 = Lfdr((xbar-mu0)/S,sigma = 1,  Acut = c(-Inf,0), method = "KW.zero")   #KW-Hybrid
      decisions.npmle1 = sc.func(tkw,q)
      fdp.npmle1.r = sum((1-theta)*decisions.npmle1$de)/max(sum(decisions.npmle1$de),1)
      etp.npmle1.r=sum((theta)*decisions.npmle1$de)/max(sum(theta), 1)
      decisions.npmle2 = sc.func(tkw2,q)
      fdp.npmle2.r = sum((1-theta)*decisions.npmle2$de)/max(sum(decisions.npmle2$de),1)
      etp.npmle2.r=sum((theta)*decisions.npmle2$de)/max(sum(theta), 1)
      
      #### Efron's Deconv based Clfdr ###################################
      deconv.g<- tryCatch({BDE(xbar,S)},error=function(e) return(NaN))
      if(!is.list(deconv.g)){
        
        fdp.deconv.r = NA
        etp.deconv.r = NA
      } else{
        w.deconv=deconv.g$y/sum(deconv.g$y)
        gd.deconv=deconv.g$x
        clfdr.deconv = clfdr.func(xbar,S,mu0,w.deconv,gd.deconv,
                                  type='npmle')
        decisions.deconv = sc.func(clfdr.deconv,q)
        fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
        etp.deconv.r=sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
      }
      
      ################################ ash alpha = 0 #####################
      ash_out = ash(xbar,S,mixcompdist = "+uniform",mode="estimate",
                    outputlevel=3)
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=mu0))/3000)
      decisions.ash = sc.func(clfdr.ash,q)
      fdp.ash.r = sum((1-theta)*decisions.ash$de)/max(sum(decisions.ash$de),1)
      etp.ash.r=sum((theta)*decisions.ash$de)/max(sum(theta), 1)
      rm(post_sample_ash)
      
      ####################################### ash with alpha = 1 ##################
      ash_out = ash(xbar,S,mixcompdist = "+uniform",alpha=1,outputlevel=3,mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash.1 = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=(mu0/S[i])))/3000)
      decisions.ash.1 = sc.func(clfdr.ash.1,q)
      fdp.ash.1.r = sum((1-theta)*decisions.ash.1$de)/max(sum(decisions.ash.1$de),1)
      etp.ash.1.r=sum((theta)*decisions.ash.1$de)/max(sum(theta), 1)
      rm(post_sample_ash)
      
      ################# NPMLE with basis expansion ####################
      m = 50
      K = 10
      gd=seq(from=min(xbar),to=max(xbar),length.out=m)
      npmleb = npmleb_opt(xbar,S,gd,length_out = m,bases = K)
      g_npmle = npmleb$probs
      clfdr.npmleb = clfdr.func(xbar,S,mu0,g_npmle,gd,
                                type='proposed')
      decisions.npmleb = sc.func(clfdr.npmleb,q)
      fdp.npmleb.r = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
      etp.npmleb.r=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
      
      ################# Proposed ####################
      hamt = hamt_opt(xbar,S,f_hat, gd,length_out = m,bases = K)
      w_mosek = hamt$probs
      
      clfdr.dd = clfdr.func(xbar,S,mu0,w_mosek,gd,
                            type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      
      return(list("fdr.or_r"=fdp.or.r,
                  "fdr.npmle1_r"=fdp.npmle1.r,
                  "fdr.npmle2_r"=fdp.npmle2.r,
                  "fdr.deconv_r"=fdp.deconv.r,
                  "fdr.ash_r"=fdp.ash.r,
                  "fdr.dd_r"=fdp.dd.r,
                  "fdr.ash.1_r"=fdp.ash.1.r,
                  "fdr.npmleb_r"=fdp.npmleb.r,
                  "etp.ash_r"=etp.ash.r,
                  "etp.dd_r"=etp.dd.r,
                  "etp.ash.1_r"=etp.ash.1.r,
                  "etp.npmleb_r"=etp.npmleb.r,
                  "etp.npmle1_r"=etp.npmle1.r,
                  "etp.npmle2_r"=etp.npmle2.r,
                  "etp.deconv_r"=etp.deconv.r,
                  "etp.or_r"=etp.or.r,
                  "fdr.adapt_r"=fdp.adapt.r,
                  "etp.adapt_r"=etp.adapt.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.or[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.or_r)
    fdr.npmle1[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle1_r)
    fdr.npmle2[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle2_r)
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
    etp.npmle1[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle1_r)
    etp.npmle2[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle2_r)
    etp.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.deconv_r)
    fdr.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.adapt_r)
    etp.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.adapt_r)
    print(paste(u[nn],' -- ',
                mean(sapply(1:reps,function(i) result[[i]]$fdr.dd_r)),sep=''))
  }
  result<-list("fdr.or"=fdr.or,
               "fdr.npmle1"=fdr.npmle1,
               "fdr.npmle2"=fdr.npmle2,
               "fdr.deconv"=fdr.deconv,
               "fdr.ash"=fdr.ash,
               "fdr.ash.1"=fdr.ash.1,
               "fdr.npmleb"=fdr.npmleb,
               "fdr.dd"=fdr.dd,
               "etp.ash"=etp.ash,
               "etp.ash.1"=etp.ash.1,
               "etp.npmleb"=etp.npmleb,
               "etp.dd"=etp.dd,
               "etp.npmle1"=etp.npmle1,
               "etp.npmle2"=etp.npmle2,
               "etp.deconv"=etp.deconv,
               "etp.or"=etp.or,
               "fdr.adapt"=fdr.adapt,
               "etp.adapt"=etp.adapt)
}

#--------------------------------------------------------------------

n<-10000
u = c(20,40,60,80,100)
mu0= 1
q=0.1
reps=500

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig16_set2.RData')

