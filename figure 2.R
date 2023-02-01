
library(REBayes)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(rmutil)
library(Rmosek)


n<-10000
u = c(2)
length_out = 50
bases = 10
mu0= 4
q=0.1
reps=1

cl <- makeCluster(6)
registerDoParallel(cl)

result<-foreach(r = 1:reps,.packages=c("REBayes","MASS","Rmosek","rmutil"))%dopar%{
  source('funcs.R')
  r=13
  set.seed(r)
  sig = runif(n,0.25,u[1])
  mu = 3*sig
  x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
  theta=1*(mu>mu0)
  
  f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
  
  #### NPMLE based Clfdr ###################################
  npmle.g<- GLmix(x,sigma=sig)
  w.npmle=npmle.g$y
  gd.npmle=npmle.g$x
  clfdr.npmle = clfdr.func(x,sig,mu0,w.npmle,gd.npmle,
                           f_hat,type='npmle')
  decisions.npmle = sc.func(clfdr.npmle,q)
  fdp.npmle.r = sum((1-theta)*decisions.npmle$de)/max(sum(decisions.npmle$de),1)
  etp.npmle.r=sum((theta)*decisions.npmle$de)/max(sum(theta), 1)
  
  #### Efron's Deconv based Clfdr ###################################
  deconv.g<- BDE(x,sig)
  w.deconv=deconv.g$y/sum(deconv.g$y)
  gd.deconv=deconv.g$x
  clfdr.deconv = clfdr.func(x,sig,mu0,w.deconv,gd.deconv,
                            f_hat,type='npmle')
  decisions.deconv = sc.func(clfdr.deconv,q)
  fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
  etp.deconv.r = sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
  
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
  H = crossprod(A_new)
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
  clfdr.dd = clfdr.func(x,sig,mu0,w_mosek_quad,gd,
                        f_hat,type='proposed')
  decisions.dd = sc.func(clfdr.dd,q)
  fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
  etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
  
  #################################### BH ###################################
  z=(x-mu0)/sig
  pv=1-pnorm(z, 0, 1)
  bh.res<-bh.func(pv, q)
  fdp.bh.r = sum((1-theta)*bh.res$de)/max(sum(bh.res$de),1)
  etp.bh.r=sum((theta)*bh.res$de)/max(sum(theta), 1)
  
  return(list("fdr.dd_r"=fdp.dd.r,"fdr.npmle_r"=fdp.npmle.r,"fdr.deconv_r"=fdp.deconv.r,
              "fdr.bh_r"=fdp.bh.r,
              "etp.dd_r"=etp.dd.r,"etp.npmle_r"=etp.npmle.r,"etp.deconv_r"=etp.deconv.r,
              "etp.bh_r"=etp.bh.r,"x"=x,"sig"=sig,'decisions.npmle'=decisions.npmle,
              'decisions.dd'=decisions.dd,'decisions.deconv'=decisions.deconv,'theta'=theta))
}
stopCluster(cl)
registerDoSEQ()

#--------------------------------------------------------------------

r=1
color.code=rep('gray',n)
for(i in 1:n){
  if(result[[r]]$decisions.npmle$de[i]==1){
    color.code[i]='#FF0000'
  }
}
plotdata1<- data.frame('x'=result[[r]]$x,'sig'=result[[r]]$sig,
                       'dec_true'=as.character(result[[r]]$decisions.npmle$de),
                       'cols'=color.code)
g4<-ggplot()+geom_point(data=plotdata1,aes(x=x,y=sig,color=dec_true,shape=dec_true))+
  scale_colour_manual(breaks = plotdata1$dec_true, 
                      values = plotdata1$cols)+
  geom_hline(aes(yintercept=(mu0/3)),linetype = "dashed")+
  theme_bw()+ggtitle('NPMLE')+
  ylab(TeX("$\\sigma$"))+
  theme(legend.position='none',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

color.code=rep('gray',n)
for(i in 1:n){
  if(result[[r]]$decisions.deconv$de[i]==1){
    color.code[i]='#FF0000'
  }
}
plotdata2<- data.frame('x'=result[[r]]$x,'sig'=result[[r]]$sig,
                       'dec_deconv'=as.character(result[[r]]$decisions.deconv$de),
                       'cols'=color.code)
g5<-ggplot()+geom_point(data=plotdata2,aes(x=x,y=sig,color=dec_deconv,
                                           shape=dec_deconv))+
  scale_colour_manual(breaks = plotdata2$dec_deconv, 
                      values = plotdata2$cols)+
  geom_hline(aes(yintercept=(mu0/3)),linetype = "dashed")+
  theme_bw()+ggtitle('DECONV')+ylab(TeX("$\\sigma$"))+
  theme(legend.position='none',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

color.code=rep('gray',n)
for(i in 1:n){
  if(result[[r]]$decisions.dd$de[i]==1){
    color.code[i]='#FF0000'
  }
}
plotdata3<- data.frame('x'=result[[r]]$x,'sig'=result[[r]]$sig,
                       'dec_dd'=as.character(result[[r]]$decisions.dd$de),
                       'cols'=color.code)
g6<-ggplot()+geom_point(data=plotdata3,aes(x=x,y=sig,color=dec_dd,
                                           shape=dec_dd))+
  scale_colour_manual(breaks = plotdata3$dec_dd, 
                      values = plotdata3$cols)+
  geom_hline(aes(yintercept=(mu0/3)),linetype = "dashed")+
  theme_bw()+ggtitle('HAMT')+ylab(TeX("$\\sigma$"))+
  theme(legend.position='none',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(g4,g5,g6,ncol=3,nrow=1,legend="none")


