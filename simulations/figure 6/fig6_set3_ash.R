

library(foreach)
library(doParallel)
library(ashr)
library(Rmosek)

exp_1<-function(n,mu0,u,q,reps){
  
  fdr.ash<-matrix(0,reps,length(u))
  etp.ash<-matrix(0,reps,length(u))
  fdr.ash.1<-matrix(0,reps,length(u))
  etp.ash.1<-matrix(0,reps,length(u))
  fdr.npmleb<-matrix(0,reps,length(u))
  etp.npmleb<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("ashr","Rmosek"))%dopar%{
      source("../funcs.R")
      set.seed(r)
      p<-runif(n)
      sig = runif(n,0.5,u[nn])
      mu = sapply(1:n,function(i)(p[i]<=0.9)*rnorm(1,-sig[i],0.5)+
                    (p[i]>0.9)*2*(sig[i])^2)
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      
      ash_out = ash(x,sig,mixcompdist = "halfuniform",outputlevel=3,
                    mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=mu0))/3000)
      decisions.ash = sc.func(clfdr.ash,q)
      fdp.ash.r = sum((1-theta)*decisions.ash$de)/max(sum(decisions.ash$de),1)
      etp.ash.r=sum((theta)*decisions.ash$de)/max(sum(theta), 1)
      
      ####################################### alpha = 1 ##################
      ash_out = ash(x,sig,mixcompdist = "halfuniform",alpha=1,outputlevel=3,
                    mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash.1 = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=(mu0/sig[i])))/3000)
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
      clfdr.npmleb = clfdr.func(x,sig,mu0,g_npmle,gd,
                                type='proposed')
      decisions.npmleb = sc.func(clfdr.npmleb,q)
      fdp.npmleb.r = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
      etp.npmleb.r=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
      
      return(list("fdr.ash_r"=fdp.ash.r,
                  "fdr.ash.1_r"=fdp.ash.1.r,
                  "fdr.npmleb_r"=fdp.npmleb.r,
                  "etp.ash_r"=etp.ash.r,
                  "etp.ash.1_r"=etp.ash.1.r,
                  "etp.npmleb_r"=etp.npmleb.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash_r)
    etp.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash_r)
    fdr.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash.1_r)
    etp.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash.1_r)
    fdr.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmleb_r)
    etp.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmleb_r)
    print(nn)
  }
  result<-list("fdr.ash"=fdr.ash,
               "fdr.ash.1"=fdr.ash.1,
               "fdr.npmleb"=fdr.npmleb,
               "etp.ash"=etp.ash,
               "etp.ash.1"=etp.ash.1,
               "etp.npmleb"=etp.npmleb)
}
#--------------------------------------------------------------------

n<-10000
u = c(1,1.1,1.2,1.3,1.4,1.5)
mu0= 1
q=0.1
reps=100

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig6_set3_ash.RData')

