

library(foreach)
library(doParallel)
library(CAMT)
library(splines)
library(AdaPTGMM)

exp_1<-function(n,mu0,u,q,reps){
  
  fdr.adapt<-matrix(0,reps,length(u))
  fdr.camt<-matrix(0,reps,length(u))
  
  etp.adapt<-matrix(0,reps,length(u))
  etp.camt<-matrix(0,reps,length(u))
  
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(8)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("CAMT","AdaPTGMM","splines"))%dopar%{
      set.seed(r)
      sig = runif(n,0.25,u[nn])
      mu = 3*sig
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      
      #################################### p-values ###################################
      z=(x-mu0)/sig
      pv=1-pnorm(z, 0, 1)
      
      ###################### CAMT ########################################
      out_camt = tryCatch({camt.fdr(pv,pi0.var=sig,f1.var=sig)},error=function(e) return(NaN))
      if(!is.list(out_camt)){
        
        fdp.camt.r = NA
        etp.camt.r = NA
      } else{
        camt.de = rep(0,n)
        camt.de[out_camt$fdr<=q]=1
        fdp.camt.r = sum((1-theta)*camt.de)/max(sum(camt.de),1)
        etp.camt.r=sum((theta)*camt.de)/max(sum(theta), 1)
      }
      
      ######################### AdaPTGMM ################################
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = sig),pvals=pv,alphas=q,
                                     beta_formulas = formulas,nfits=10)},
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
      return(list("fdr.camt_r"=fdp.camt.r,"fdr.adapt_r"=fdp.adapt.r,
                  "etp.camt_r"=etp.camt.r,"etp.adapt_r"=etp.adapt.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.camt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.camt_r)
    fdr.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.adapt_r)
    etp.camt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.camt_r)
    etp.adapt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.adapt_r)
    print(nn)
  }
  result<-list("fdr.camt"=fdr.camt,"fdr.adapt"=fdr.adapt,
               "etp.camt"=etp.camt,"etp.adapt"=etp.adapt)
}

#--------------------------------------------------------------------

n<-10000
u = c(1.5,1.6,1.7,1.8,1.9,2)
mu0 = 4
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig8_set5_adapt.RData')

