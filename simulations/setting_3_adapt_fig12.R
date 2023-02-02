

library(foreach)
library(doParallel)
library(AdaPTGMM)
library(splines)

exp_1<-function(n,mu0,mu1,u,q,reps){
  
  fdr.adapt<-matrix(0,reps,length(u))
  etp.adapt<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(8)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("CAMT","AdaPTGMM","splines"))%dopar%{
      mm = 100
      set.seed(r)
      sig = sqrt(mm)*sample(c(0.5,1,3),n,replace = TRUE)
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+(-u[nn]*sig[i]/sqrt(mm))*(p[i]>0.9 & p[i]<=0.95)+
                    (u[nn]*sig[i]/sqrt(mm))*(p[i]>0.95))
      x=sapply(1:n, function(i)rnorm(mm,mu[i],sig[i]))
      theta=1*(mu>mu1)+1*(mu<mu0)
      
      xbar = colMeans(x)
      S =sqrt(sapply(1:n, function(i) var(x[,i])/mm))
      
      #################################### p-values ###################################
      z=xbar/S
      pv=1-pnorm(abs(z)+5, 0, 1)+pnorm(-abs(z)+5,0,1)
      
      ######################### AdaptMT ################################
      # dist <- beta_family()
      # formulas <- paste0("ns(x, df = ", 6:10, ")")
      # adapt.res<-tryCatch({adapt_glm(x = data.frame(x = sig), pvals = pv, pi_formulas = formulas,alphas = q,
      #                                mu_formulas = formulas, dist = dist, nfits = 10)},
      #                     error=function(e) return(NaN))
      # 
      # if(!is.list(adapt.res)){
      #   
      #   fdp.adapt.r = NA
      #   etp.adapt.r = NA
      # } else{
      #   
      #   adapt.de=rep(0,n)
      #   adapt.de[adapt.res$rejs]=1
      #   fdp.adapt.r = sum((1-theta)*adapt.de)/max(sum(adapt.de),1)
      #   etp.adapt.r=sum((theta)*adapt.de)/max(sum(theta), 1)
      # }
      ######################### AdaPTGMM ################################
      formulas <- paste("ns(x, df = ",c(2,4,6)," )")
      adapt.res<-tryCatch({adapt_gmm(x=data.frame(x = S),z=xbar,se=S,
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
u = c(2,2.1,2.2,2.3,2.4,2.5)
mu0 = -5
mu1 = 5
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,mu1,u,q,reps)

save.image(paste(getwd(),'/setting_3_adapt_fig12.RData',sep=''))

