
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
                           type='npmle')
  decisions.npmle = sc.func(clfdr.npmle,q)
  fdp.npmle.r = sum((1-theta)*decisions.npmle$de)/max(sum(decisions.npmle$de),1)
  etp.npmle.r=sum((theta)*decisions.npmle$de)/max(sum(theta), 1)
  
  #### Efron's Deconv based Clfdr ###################################
  deconv.g<- BDE(x,sig)
  w.deconv=deconv.g$y/sum(deconv.g$y)
  gd.deconv=deconv.g$x
  clfdr.deconv = clfdr.func(x,sig,mu0,w.deconv,gd.deconv,
                            type='npmle')
  decisions.deconv = sc.func(clfdr.deconv,q)
  fdp.deconv.r = sum((1-theta)*decisions.deconv$de)/max(sum(decisions.deconv$de),1)
  etp.deconv.r = sum((theta)*decisions.deconv$de)/max(sum(theta), 1)
  
  ################# Proposed ####################
  gd=seq(from=min(x),to=max(x),length.out=length_out)
  hamt = hamt_opt(x,sig,f_hat,gd,length_out,bases)

  clfdr.dd = clfdr.func(x,sig,mu0,hamt$probs,gd,
                        type='proposed')
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


