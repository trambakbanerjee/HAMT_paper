
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
      mm = rep(u[nn],n)
      p<-runif(n)
      sig = sqrt(mm)*runif(n,0.5,1.5)
      mu = sapply(1:n,function(i)(p[i]<=0.9)*rnorm(1,-sig[i]/sqrt(mm[i]),0.5)+
                    (p[i]>0.9)*2*((sig[i])^2/mm[i]))
      x=sapply(1:n, function(i) mu[i]+(1/sqrt(3))*runif(mm[i],-3*sig[i],3*sig[i]))
      theta=1*(mu>mu0)
      xbar = colMeans(x)
      S =sqrt(sapply(1:n, function(i) var(x[,i])/mm[i]))
      
      ######################### AdaPTGMM ################################
      pv = 1-pnorm(xbar,mu0,S)
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = S),pvals=pv,alphas=q,
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
u = c(30,40,50,60,70,80,90,100)
mu0= 1
q=0.1
reps=100

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig15_set1_adapt.RData')

