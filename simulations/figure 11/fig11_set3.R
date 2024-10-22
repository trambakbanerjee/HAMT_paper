
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(ashr)
library(Rmosek)
library(REBayes)
library(AdaPTGMM)


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
    
    result<-foreach(r = 1:reps,.packages=c("ashr","Rmosek","REBayes","AdaPTGMM"))%dopar%{
      source('../funcs.R')
      set.seed(r)
      sig = runif(n,0.5,u[nn])
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.8)+
                    (p[i]>0.8 & p[i]<=0.9)*rnorm(1,sig[i],1)+
                    (p[i]>0.9)*rnorm(1,sig[i],2))
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu<mu0)+1*(mu>mu1)
      
      ################## Oracle Clfdr ###############################
      clfdr.or = rep(0,n)
      for(i in 1:n){
        ff<-function(u,x,s_x,s_u,mm){dnorm(x,u,s_x)*dnorm(u,mm,s_u)
        } 
        denom = 0.8*dnorm(x[i],0,sig[i])+
          0.1*integrate(ff,-Inf,Inf,x=x[i],s_x=sig[i],s_u=1,mm=sig[i])$value+
          0.1*integrate(ff,-Inf,Inf,x=x[i],s_x=sig[i],s_u=2,mm=sig[i])$value
        numer = 0.8*dnorm(x[i],0,sig[i])+
          0.1*integrate(ff,-1,0,x=x[i],s_x=sig[i],s_u=1,mm=sig[i])$value+
          0.1*integrate(ff,-1,0,x=x[i],s_x=sig[i],s_u=2,mm=sig[i])$value+
          0.1*integrate(ff,0,1,x=x[i],s_x=sig[i],s_u=1,mm=sig[i])$value+
          0.1*integrate(ff,0,1,x=x[i],s_x=sig[i],s_u=2,mm=sig[i])$value
        clfdr.or[i]=numer/denom
      }
      decisions.or = sc.func(clfdr.or,q)
      fdp.or.r = sum((1-theta)*decisions.or$de)/max(sum(decisions.or$de),1)
      etp.or.r=sum((theta)*decisions.or$de)/max(sum(theta), 1)
      
      ######################### AdaPTGMM ################################
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = sig),z=x,se=sig,
                                     testing = "interval",
                                     rendpoint = mu1,
                                     lendpoint = mu0,
                                     alphas=q,beta_formulas = formulas,nfits=10)},
                          error=function(e) return(NaN))
      
      if(!is.list(adapt.res)){
        
        fdp.adapt.r = NA
        etp.adapt.r = NA
      } else{
        
        adapt.de=rep(0,n)
        adapt.de[adapt.res$rejs[[1]]]=1
        fdp.adapt.r = sum((1-theta)*adapt.de)/max(sum(adapt.de),1)
        etp.adapt.r=sum((theta)*adapt.de)/max(sum(theta), 1)
      }
      
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
      w.deconv=deconv.g$y/sum(deconv.g$y)
      gd.deconv=deconv.g$x
      clfdr.deconv = clfdr.func.2(x,sig,mu0,mu1,w.deconv,gd.deconv,
                                  type='npmle')
      decisions.deconv = sc.func(clfdr.deconv,q)
      fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
      etp.deconv.r=sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
      
      ########################## alpha = 0 ash ##############################
      ash_out = ash(x,sig,mixcompdist = "normal",
                    outputlevel=3,mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]>=mu0 &
                                                  post_sample_ash[,i]<=mu1))/3000)
      decisions.ash = sc.func(clfdr.ash,q)
      fdp.ash.r = sum((1-theta)*decisions.ash$de)/max(sum(decisions.ash$de),1)
      etp.ash.r=sum((theta)*decisions.ash$de)/max(sum(theta), 1)
      rm(post_sample_ash)
      
      ####################################### alpha = 1 ##################
      ash_out = ash(x,sig,mixcompdist = "normal",alpha=1,
                    outputlevel=3,mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash.1 = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]>=(mu0/sig[i]) &
                                                    post_sample_ash[,i]<=(mu1/sig[i])))/3000)
      decisions.ash.1 = sc.func(clfdr.ash.1,q)
      fdp.ash.1.r = sum((1-theta)*decisions.ash.1$de)/max(sum(decisions.ash.1$de),1)
      etp.ash.1.r=sum((theta)*decisions.ash.1$de)/max(sum(theta), 1)
      rm(post_sample_ash)
      ################# NPMLE with basis expansion ####################
      m = 50
      K = 10
      gd=seq(from=min(x),to=max(x),length.out=m)
      npmleb = npmleb_opt(x,sig,gd,length_out = m,bases = K)
      g_npmle = npmleb$probs
      clfdr.npmleb = clfdr.func.2(x,sig,mu0,mu1,g_npmle,gd,
                                  type='proposed')
      decisions.npmleb = sc.func(clfdr.npmleb,q)
      fdp.npmleb.r = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
      etp.npmleb.r=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
      
      ################# Proposed ####################
      f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
      m = 50
      K = 10
      gd=seq(from=min(x),to=max(x),length.out=m)
      hamt = hamt_opt(x,sig,f_hat,gd,length_out = m,bases = K)
      w_mosek = hamt$probs
      
      clfdr.dd = clfdr.func.2(x,sig,mu0,mu1,w_mosek,gd,
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
u = c(1,1.1,1.2,1.3,1.4,1.5)
mu0 = -1
mu1 = 1
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,mu1,u,q,reps)

save.image('fig11_set3.RData')
