
library(REBayes)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(splines)


source('funcs.R')

n = 10000
reps=1

marg_deconv_1 = marg_npmle_1 =  marg_dd_1 =  marg_orc_1 = matrix(0,n,reps)
marg_deconv_2 = marg_npmle_2 =  marg_dd_2 =  marg_orc_2 = matrix(0,n,reps)
marg_deconv_3 = marg_npmle_3 =  marg_dd_3 =  marg_orc_3 = matrix(0,n,reps)
sig_val = matrix(0,reps,3)

for(r in 1:reps){
  
  set.seed(r)
  sig = runif(n,0.5,3)
  mu =  3*sig
  set.seed(r+2)
  x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
  
  # NPMLE 
  npmle = GLmix(x,sigma=sig)
  
  #### Efron's Deconv based Clfdr ###################################
  deconv.g<- BDE(x,sig)
  w.deconv=deconv.g$y/sum(deconv.g$y)
  gd.deconv=deconv.g$x

  
  f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
  ################# Proposed ####################
  gd=seq(from=min(x),to=max(x),length.out=50)
  hamt = hamt_opt(x,sig,f_hat,gd)
  w_mosek = hamt$probs

  #### Marginal density plot at sig = 1
  a= seq(-5,15,length.out=1000)
  id.sig = which.min(abs(sig-1))
  sig_val[r,1] = sig[id.sig]
  marg_deconv_1[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd.deconv,0,sig[id.sig])*w.deconv))
  marg_npmle_1[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-npmle$x,0,sig[id.sig])*npmle$y))
  marg_dd_1[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd,0,sig[id.sig])*w_mosek[,id.sig]))
  marg_orc_1[,r] = sapply(1:length(a),function(i) dnorm(a[i]-3*sig[id.sig],0,sig[id.sig]))
  
  #### Marginal density plot at sig = 1.5
  id.sig = which.min(abs(sig-1.5))
  sig_val[r,2] = sig[id.sig]
  marg_deconv_2[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd.deconv,0,sig[id.sig])*w.deconv))
  marg_npmle_2[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-npmle$x,0,sig[id.sig])*npmle$y))
  marg_dd_2[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd,0,sig[id.sig])*w_mosek[,id.sig]))
  marg_orc_2[,r] = sapply(1:length(a),function(i) dnorm(a[i]-3*sig[id.sig],0,sig[id.sig]))
  
  #### Marginal density plot at sig = 2
  a= seq(-5,15,length.out=1000)
  id.sig = which.min(abs(sig-2))
  sig_val[r,3] = sig[id.sig]
  marg_deconv_3[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd.deconv,0,sig[id.sig])*w.deconv))
  marg_npmle_3[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-npmle$x,0,sig[id.sig])*npmle$y))
  marg_dd_3[,r] = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd,0,sig[id.sig])*w_mosek[,id.sig]))
  marg_orc_3[,r] = sapply(1:length(a),function(i) dnorm(a[i]-3*sig[id.sig],0,sig[id.sig]))
  
  print(r)
}

library(ggplot2)
library(ggpubr)
library(latex2exp)
a= seq(-5,15,length.out=1000)
plotdata1 = data.frame('f'=c(rowMeans(marg_orc_1),rowMeans(marg_dd_1),rowMeans(marg_npmle_1),
                             rowMeans(marg_deconv_1)))
plotdata1$x = c(a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE','DECONV'),each=n))
g1<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma=1)$"))+geom_vline(xintercept = 3*1,color='black',linetype='dashed')+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="white"),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotdata1 = data.frame('f'=c(rowMeans(marg_orc_2),rowMeans(marg_dd_2),rowMeans(marg_npmle_2),
                             rowMeans(marg_deconv_2)))
plotdata1$x = c(a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE','DECONV'),each=n))
g2<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma=1.5)$"))+geom_vline(xintercept = 3*1.5,color='black',linetype='dashed')+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="white"),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

a= seq(-5,15,length.out=1000)
plotdata1 = data.frame('f'=c(rowMeans(marg_orc_3),rowMeans(marg_dd_3),rowMeans(marg_npmle_3),
                       rowMeans(marg_deconv_3)))
plotdata1$x = c(a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE','DECONV'),each=n))
g3<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma=2)$"))+geom_vline(xintercept = 3*2,color='black',linetype='dashed')+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="white"),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(g1,g2,g3,nrow=1,ncol=3,common.legend = T)


