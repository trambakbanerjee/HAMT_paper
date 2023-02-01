
library(REBayes)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(rmutil)
library(Rmosek)



exp_1<-function(n,mu0,mu1,u,q,length_out,bases,reps){
  
  fdr.dd<-matrix(0,reps,length(u))
  fdr.npmle<-matrix(0,reps,length(u))
  fdr.bh<-matrix(0,reps,length(u))
  fdr.or<-matrix(0,reps,length(u))
  fdr.deconv<-matrix(0,reps,length(u))
  
  etp.dd<-matrix(0,reps,length(u))
  etp.npmle<-matrix(0,reps,length(u))
  etp.bh<-matrix(0,reps,length(u))
  etp.or<-matrix(0,reps,length(u))
  etp.deconv<-matrix(0,reps,length(u))
  
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("REBayes","MASS","Rmosek","rmutil"))%dopar%{
      source('funcs.R')
      set.seed(r)
      sig = runif(n,0.5,u[nn])
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.9)+rnorm(1,3,sig[i])*(p[i]>0.9 & p[i]<=0.95)+
                    rnorm(1,-3,sig[i])*(p[i]>=0.95))
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu1)+1*(mu<mu0)
      
      f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
      ################## Oracle Clfdr ###############################
      clfdr.or = rep(0,n)
      for(i in 1:n){
        ff<-function(u,x,s,mm){dnorm(x,u,s)*dnorm(u,mm,s)
        } 
        denom = 0.9*dnorm(x[i],0,sig[i])+0.05*integrate(ff,-Inf,Inf,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,-Inf,Inf,x=x[i],s=sig[i],mm=-3)$value
        numer = 0.9*dnorm(x[i],0,sig[i])+0.05*integrate(ff,-2,0,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,-2,0,x=x[i],s=sig[i],mm=-3)$value+
          0.05*integrate(ff,0,2,x=x[i],s=sig[i],mm=3)$value+
          0.05*integrate(ff,0,2,x=x[i],s=sig[i],mm=-3)$value
        clfdr.or[i]=numer/denom
      }
      decisions.or = sc.func(clfdr.or,q)
      fdp.or.r = sum((1-theta)*decisions.or$de)/max(sum(decisions.or$de),1)
      etp.or.r=sum((theta)*decisions.or$de)/max(sum(theta), 1)
      
      #### NPMLE based Clfdr ###################################
      npmle.g<- GLmix(x,sigma=sig)
      w.npmle=npmle.g$y
      gd.npmle=npmle.g$x
      clfdr.npmle = clfdr.func.2(x,sig,mu0,mu1,w.npmle,gd.npmle,
                                 f_hat,type='npmle')
      decisions.npmle = sc.func(clfdr.npmle,q)
      fdp.npmle.r = sum((1-theta)*decisions.npmle$de)/max(sum(decisions.npmle$de),1)
      etp.npmle.r=sum((theta)*decisions.npmle$de)/max(sum(theta), 1)
      
      #### Efron's Deconv based Clfdr ###################################
      deconv.g<- BDE(x,sig)
      w.deconv=deconv.g$y
      gd.deconv=deconv.g$x
      clfdr.deconv = clfdr.func.2(x,sig,mu0,mu1,w.deconv,gd.deconv,
                                  f_hat,type='npmle')
      decisions.deconv = sc.func(clfdr.deconv,q)
      fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
      etp.deconv.r=sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
      
      ################# Proposed ####################
      m = length_out
      K = bases
      gd=seq(from=min(x),to=max(x),length.out=m)
      B=matrix(0,n,K)
      for(i in 1:n){
        for(k in 1:K){
          B[i,k]=0.5*cos(k*sig[i])+0.5;
        }
      }
      A=matrix(0,n,m*K);
      for(i in 1:n){
        A[i,]=unlist(sapply(1:m,function(j) rep(dnorm(x[i]-gd[j],0,sig[i]),K)*B[i,]))
      }
      
      Q=matrix(0,n,m*K)
      for(i in 1:n){
        Q[i,] = rep(B[i,],m)
      }
      A_new = A
      H = crossprod(A_new)#t(A)%*%A
      H = H+10^{-5}*diag(m*K)
      #Sparse matrix construction
      Q_sparse = sparseMatrix(i=c(rep(1:n,each=m*K)),j=c(rep(1:(m*K),n)),x=c(t(Q)))
      BB = c(t(B))
      rr = c(unlist(sapply(1:n, function(i) rep((1+(i-1)*K):(i*K),m))))
      L_sparse = sparseMatrix(i=rep(1:(n*m),each = K),j=rep(1:(m*K),n),x=BB[rr])
      
      # Specify the non-quadratic part of the problem.
      prob <- list(sense="min")
      prob$c <- -t(A_new)%*%f_hat
      prob$A <- rbind(Q_sparse,L_sparse)
      prob$bc <- rbind(blc=c(rep(1,n),rep(0,n*m)), 
                       buc=c(rep(1,n),rep(Inf,n*m)))
      prob$bx <- rbind(blx=rep(-Inf,m*K), 
                       bux=rep(Inf,m*K))
      
      # Specify the quadratic objective matrix in triplet form.
      prob$qobj$i <- c(unlist(sapply(1:(m*K),function(i) i:(m*K))))
      prob$qobj$j <- rep(1:(m*K),(m*K):1)
      prob$qobj$v <- c(H[lower.tri(H,diag = T)])
      #prob$iparam <- list(PRESOLVE_USE = "PRESOLVE_MODE_OFF")
      
      # Solve the problem
      out_mosek_quad <- mosek(prob)
      if(out_mosek_quad$sol$itr$prosta=="PRIMAL_INFEASIBLE"){
        prob$bc <- rbind(blc=c(rep(0.9,n),rep(0,n*m)), 
                         buc=c(rep(1,n),rep(Inf,n*m)))
        out_mosek_quad <- mosek(prob)
        f = out_mosek_quad$sol$itr$xx
        f<- matrix(f,ncol=m,byrow=F)
        w_mosek_quad = t(B%*%f)
        for(i in 1:n){
          w_mosek_quad[,i] = w_mosek_quad[,i]/sum(w_mosek_quad[,i])
        }
      } else {
        f = out_mosek_quad$sol$itr$xx
        f<- matrix(f,ncol=m,byrow=F)
        w_mosek_quad = t(B%*%f)
      }
      clfdr.dd = clfdr.func.2(x,sig,mu0,mu1,w_mosek_quad,gd,
                              f_hat,type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      #################################### BH ###################################
      #z=(x-mu0)/sig
      #pv=1-pnorm(z, 0, 1)
      #bh.res<-bh.func(pv, q)
      fdp.bh.r = 0#sum((1-theta)*bh.res$de)/max(sum(bh.res$de),1)
      etp.bh.r=0#sum((theta)*bh.res$de)/max(sum(theta), 1)
      
      
      return(list("fdr.dd_r"=fdp.dd.r,"fdr.npmle_r"=fdp.npmle.r,"fdr.deconv_r"=fdp.deconv.r,
                  "fdr.bh_r"=fdp.bh.r,"fdr.or_r"=fdp.or.r,
                  "etp.dd_r"=etp.dd.r,"etp.npmle_r"=etp.npmle.r,"etp.deconv_r"=etp.deconv.r,
                  "etp.bh_r"=etp.bh.r,"etp.or_r"=etp.or.r))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.bh_r)
    fdr.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmle_r)
    fdr.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.deconv_r)
    fdr.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.dd_r)
    fdr.or[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.or_r)
    etp.bh[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.bh_r)
    etp.npmle[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmle_r)
    etp.deconv[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.deconv_r)
    etp.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.dd_r)
    etp.or[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.or_r)
    print(paste(u[nn],' -- ',
                mean(sapply(1:reps,function(i) result[[i]]$fdr.dd_r)),sep=''))
  }
  result<-list("fdr.dd"=fdr.dd,"fdr.npmle"=fdr.npmle,"fdr.deconv"=fdr.deconv,
               "fdr.bh"=fdr.bh,"fdr.or"=fdr.or,
               "etp.dd"=etp.dd,"etp.npmle"=etp.npmle,"etp.deconv"=etp.deconv,
               "etp.bh"=etp.bh,"etp.or"=etp.or)
}

#--------------------------------------------------------------------

n<-10000
u = c(1,1.2,1.4,1.6,1.8,2)
length_out = 50
bases = 10
mu0 = -2
mu1 = 2
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,mu1,u,q,length_out, bases,reps)

save.image(paste(getwd(),'/setting_2_fig11.RData',sep=''))

plotdata1<- as.data.frame(c(colMeans(out_exp1$fdr.npmle),
                            colMeans(out_exp1$fdr.deconv),
                            colMeans(out_exp1$fdr.dd,na.rm=TRUE),
                            colMeans(out_exp1$fdr.or)))
names(plotdata1)<-"fdr"
plotdata1$type<-as.factor(rep(c('NPMLE','DECONV','DD','OR'),each=length(u)))
plotdata1$n<- rep(u,4)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=n,y=fdr,color=type),size=1)+
  geom_point(data=plotdata1,aes(x=n,y=fdr,color=type,shape=type),size=4,fill=NA)+
  geom_hline(yintercept = 0.1,linetype='dotted',size=1, col = 'black')+
  scale_y_continuous(limits = c(0,0.2))+
  ylab('FDP')+theme_bw()+xlab('u')+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotdata2<- as.data.frame(c(colMeans(out_exp1$etp.npmle),
                            colMeans(out_exp1$etp.deconv),
                            colMeans(out_exp1$etp.dd,na.rm=TRUE),
                            colMeans(out_exp1$etp.or)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('NPMLE','DECONV','DD','OR'),each=length(u)))
plotdata2$n<- rep(u,4)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=n,y=etp,color=type),size=1)+
  geom_point(data=plotdata2,aes(x=n,y=etp,color=type,shape=type),size=4,fill=NA)+
  ylab('ETP')+theme_bw()+xlab('u')+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggarrange(g1,g2,ncol=2,nrow=1,common.legend = TRUE)


