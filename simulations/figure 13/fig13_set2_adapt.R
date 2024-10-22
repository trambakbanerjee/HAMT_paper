
library(foreach)
library(doParallel)
library(AdaPTGMM)
library(ExtDist)

exp_1<-function(n,mu0,mu1,u,q,reps){
  
  fdr.adapt<-matrix(0,reps,length(u))
  etp.adapt<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(12)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("AdaPTGMM"))%dopar%{
      source('funcs.R')
      set.seed(r)
      sig = sample(c(0.25,0.75,1.5),n,replace = TRUE)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+(-u[nn]*sig[i])*(p[i]>0.9 & p[i]<=0.95)+
                    (u[nn]*sig[i])*(p[i]>0.95))
      x=sapply(1:n, function(i) mu[i]+rlogis(1,0,sig[i]))
      theta=1*(mu>mu1)+1*(mu<mu0)
      sig_x = sig*pi/sqrt(3)
      
      ######################### AdaPTGMM ################################
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = sig),z=x,se=sig_x,alphas=q,
                                     testing = "interval",
                                     rendpoint = mu1,
                                     lendpoint = mu0,
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
u = c(5.5,5.6,5.7,5.8,5.9,6)
mu0 = -2
mu1 = 2
q=0.1
reps=100

out_exp1<-exp_1(n,mu0,mu1,u,q,reps)

save.image('fig13_set2_adapt.RData')

