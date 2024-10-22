
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
   
    registerDoSEQ()
    result<-foreach(r = 1:reps,.packages=c("CAMT","IHW"))%dopar%{
     
      dof = 10
      mm = rep(u[nn],n)
      set.seed(r)
      p = runif(n,0,1)
      sig = sqrt(mm)*runif(n,0.5,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+
                    (p[i]>0.9)*rnorm(1,3,1))
      x=sapply(1:n, function(i)mu[i]+sig[i]*rt(mm[i],dof)*sqrt((dof-2)/dof))#sapply(1:n, function(i)runif(mm,mu[i]-3*sig[i],mu[i]+3*sig[i]))
      theta=1*(mu>mu0)
      xbar = colMeans(x)
      S =sqrt(sapply(1:n, function(i) var(x[,i])/mm[i]))
      #################################### p-values ###################################
      pv = 1-pnorm(xbar,mu0,S)
      ###################### CAMT ########################################
      out_camt = tryCatch({camt.fdr(pv,f1.var=S)},error=function(e) return(NaN))
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
      out_ihw = tryCatch({ihw(pv,S,q)},error=function(e) return(NaN))
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
u = c(20,40,60,80,100)
mu0= 2
q=0.1
reps = 100

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig17_set3_camt_ihw.RData')

