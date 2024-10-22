
library(foreach)
library(doParallel)
library(AdaPTGMM)
library(ExtDist)

exp_1<-function(n,mu0,u,q,reps){
  
  fdr.adapt<-matrix(0,reps,length(u))
  etp.adapt<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(12)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("ExtDist","AdaPTGMM"))%dopar%{
      source('funcs.R')
      set.seed(r)
      sig = runif(n,0.3,u[nn])
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+rnorm(1,3,1)*(p[i]>0.9))
      x=sapply(1:n, function(i)rLaplace(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      sig_x = sqrt(2)*sig
      
      ######################### AdaPTGMM ################################
      pv = 1-pLaplace(x,mu0,sig)
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = sig_x),pvals=pv,alphas=q,
                                     beta_formulas = formulas,nfits=5)},
                          error=function(e) return(NaN))
      
      if(!is.list(adapt.res)){
        
        fdp.adapt.r = NA
        etp.adapt.r = NA
      } else{
        
        adapt.de=rep(0,n)
        adapt.de[adapt.res$rejs[[1]]]=1
        fdp.adapt.r = sum((1-theta)*adapt.de)/max(sum(adapt.de),1)
        etp.adapt.r = sum((theta)*adapt.de)/max(sum(theta), 1)
      }
      return(list("fdr.adapt_r"=fdp.adapt.r,
                  "etp.adapt_r"=etp.adapt.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.adapt_r)
    etp.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.adapt_r)
    print(nn)
  }
  result<-list("fdr.adapt"=fdr.adapt,
               "etp.adapt"=etp.adapt)
}

#--------------------------------------------------------------------

n<-10000
u = c(0.5,0.75,1,1.25,1.5,1.75,2)
mu0 = 2
q=0.1
reps=100

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig14_set3_adapt.RData')

