
library(foreach)
library(doParallel)
library(CAMT)
library(IHW)

exp_1<-function(n,mu0,u,q,reps){
  
  fdr.camt<-matrix(0,reps,length(u))
  etp.camt<-matrix(0,reps,length(u))
  fdr.ihw<-matrix(0,reps,length(u))
  etp.ihw<-matrix(0,reps,length(u))
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(6)
    registerDoParallel(cl)
    result<-foreach(r = 1:reps,.packages=c("CAMT","IHW"))%dopar%{
      set.seed(r)
      p = runif(n,0,1)
      sig = sapply(1:n, function(i) 0.5*(p[i]<=0.33)+2*(p[i]>0.67)+
                     1*(p[i]>0.33 & p[i]<=0.67))
      set.seed(r+1)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+
                    u[nn]*sig[i]*(p[i]>0.9))
      set.seed(r+2)
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      
      #################################### p-values ###################################
      z=(x-mu0)/sig
      pv=1-pnorm(z, 0, 1)
      
      ###################### CAMT ########################################
      out_camt = tryCatch({camt.fdr(pv,f1.var=sig)},error=function(e) return(NaN))
      if(!is.list(out_camt)){
        
        fdp.camt.r = NA
        etp.camt.r = NA
      } else{
        camt.de = rep(0,n)
        camt.de[out_camt$fdr<=q]=1
        fdp.camt.r = sum((1-theta)*camt.de)/max(sum(camt.de),1)
        etp.camt.r=sum((theta)*camt.de)/max(sum(theta), 1)
      }
      
      ####################### IHW #####################################
      out_ihw = tryCatch({ihw(pv,sig,q)},error=function(e) return(NaN))
      if(!isS4(out_ihw)){
        
        fdp.ihw.r = NA
        etp.ihw.r = NA
      } else{
        ihw.de = rep(0,n)
        ihw.de[adj_pvalues(out_ihw)<=q]=1
        fdp.ihw.r = sum((1-theta)*ihw.de)/max(sum(ihw.de),1)
        etp.ihw.r=sum((theta)*ihw.de)/max(sum(theta), 1)
      }
      
      return(list("fdr.camt_r"=fdp.camt.r,
                  "fdr.ihw_r"=fdp.ihw.r,
                  "etp.ihw_r"=etp.ihw.r,
                  "etp.camt_r"=etp.camt.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.camt[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.camt_r)
    etp.camt[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.camt_r)
    fdr.ihw[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ihw_r)
    etp.ihw[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ihw_r)
    print(nn)
  }
  result<-list("fdr.camt"=fdr.camt,
               "fdr.ihw"=fdr.ihw,
               "etp.camt"=etp.camt,
               "etp.ihw"=etp.ihw)
}

#--------------------------------------------------------------------

n<-10000
u = c(2.2,2.3,2.4,2.5,2.6,2.7)
mu0 = 2
q=0.1
reps = 100

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig5_set2_camt_ihw.RData')

