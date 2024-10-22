
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
    
    result<-foreach(r = 1:reps,.packages=c("ExtDist","AdaPTGMM"))%dopar%{
      source('funcs.R')
      set.seed(r)
      mm = rep(u[nn],n)
      sig = sqrt(mm)*sample(c(0.25,0.75,1.5),n,replace = TRUE)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+(-4.5*sig[i]/sqrt(mm[i]))*(p[i]>0.9 & p[i]<=0.95)+
                    (4.5*sig[i]/sqrt(mm[i]))*(p[i]>0.95))
      x=sapply(1:n, function(i) mu[i]+(sqrt(3)/pi)*rlogis(mm[i],0,sig[i]))
      
      theta=1*(mu>mu1)+1*(mu<mu0)
      xbar = colMeans(x)
      S =sqrt(sapply(1:n, function(i) var(x[,i])/mm[i]))
      
      ######################### AdaPTGMM ################################
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = S),z=xbar,se=S,
                                     testing = "interval",
                                     rendpoint = mu1,
                                     lendpoint = mu0,
                                     alphas=q,beta_formulas = formulas,nfits=5)},
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
u = c(20,30,40,50,60,70,80,90,100)
mu0 = -5
mu1 = 5
q=0.1
reps=100

out_exp1<-exp_1(n,mu0,mu1,u,q,reps)

save.image('fig18_set4_adapt.RData')

