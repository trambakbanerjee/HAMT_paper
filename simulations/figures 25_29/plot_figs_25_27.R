
library(Rmosek)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(Rfast)
library(ashr)

n<-10000
u = c(1.5,1.7,1.9)
mu0= 4
q=0.1
r = 42


#-----------------------------------------------------------------------------a


prior_hamt = prior_npmleb_1 = prior_npmleb_4 = array(0,c(50,n,length(u)))
prior_npmleb_2 = prior_npmleb_3 = array(0,c(100,n,length(u)))

fdp.npmleb = etp.npmleb = matrix(0,length(u),4)
fdp.dd = etp.dd = matrix(0,length(u),1)

sig_n = matrix(0,n,length(u))
gd_1 = gd_4 = matrix(0,50,length(u))
gd_2 = gd_3 = matrix(0,100,length(u))

for(nn in 1:length(u)){
  
  source('../funcs.R')
  set.seed(r)
  sig = runif(n,0.25,u[nn])
  mu = 3*sig
  x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
  theta=1*(mu>mu0)
 
  f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
  
  gd=seq(from=min(x),to=max(x),length.out=50)
  
  hamt_w = hamt_opt(x,sig,f_hat,gd,length_out = 50,bases=10)
  
  clfdr.dd = clfdr.func(x,sig,mu0,hamt_w$probs,gd,type='proposed')
  decisions.dd = sc.func(clfdr.dd,q)
  fdp.dd[nn] = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
  etp.dd[nn]=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
  
  npmleb = npmleb_opt(x,sig,gd,length_out = 50,bases=10)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,1] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,1]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_hamt[,,nn] = hamt_w$probs
  prior_npmleb_1[,,nn] = npmleb$probs
  
  sig_n[,nn] = sig
  gd_1[,nn] = gd
  
  gd=seq(from=min(x),to=max(x),length.out=100)
  gd_2[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 100,bases=30)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,2] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,2]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_2[,,nn] = npmleb$probs
  
  gd=seq(from=min(x),to=max(x),length.out=100)
  gd_3[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 100,bases=10)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,3] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,3]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_3[,,nn] = npmleb$probs
  
  gd=seq(from=min(x),to=max(x),length.out=50)
  gd_4[,nn] = gd
  npmleb = npmleb_opt(x,sig,gd,length_out = 50,bases=30)
  
  clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb,q)
  fdp.npmleb[nn,4] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn,4]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  
  prior_npmleb_4[,,nn] = npmleb$probs
  
}

#--------------------------------------------------------------------

save.image('plot_figs_25_27.RData')

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)

#-------------Figure 25 ----------------------
nn = 1

#### Marginal density plot and prior probs at sig = 0.7
a= seq(-3,12,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.5))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g1<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 0.5)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 1
id.sig = which.min(abs(sig_n[,nn]-1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g2<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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


#### Marginal density plot at sig = u[nn]-0.1
id.sig = which.min(abs(sig_n[,nn]-u[nn]+0.1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g3<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

ggarrange(g1,g2,g3,nrow=1,ncol=3,common.legend = TRUE)

#--------------
#-------------Figure 26 ----------------------
nn = 2

#### Marginal density plot and prior probs at sig = 0.7
a= seq(-3,12,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.5))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g1<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 0.5)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 1
id.sig = which.min(abs(sig_n[,nn]-1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g2<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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


#### Marginal density plot at sig = u[nn]-0.1
id.sig = which.min(abs(sig_n[,nn]-u[nn]+0.1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g3<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

ggarrange(g1,g2,g3,nrow=1,ncol=3,common.legend = TRUE)

#--------------
#-------------Figure 27 ----------------------
nn = 3

#### Marginal density plot and prior probs at sig = 0.7
a= seq(-3,12,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.5))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g1<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 0.5)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 1
id.sig = which.min(abs(sig_n[,nn]-1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g2<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = 1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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


#### Marginal density plot at sig = u[nn]-0.1
id.sig = which.min(abs(sig_n[,nn]-u[nn]+0.1))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-3*sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

plotdata1 = data.frame('f'=c(marg_orc,marg_dd,
                             marg_npmleb_1,marg_npmleb_2,
                             marg_npmleb_3,marg_npmleb_4))
plotdata1$x = c(a,a,a,a,a,a)
plotdata1$type = as.factor(rep(c('Oracle','HAMT','NPMLE B(S=50, K=10)',
                                 'NPMLE B(S=100, K=30)','NPMLE B(S=100, K=10)',
                                 'NPMLE B(S=50, K=30)'),each=length(a)))
g3<- ggplot(data=plotdata1)+geom_line(aes(x=x,y=f,color=type,linetype=type),size=1)+
  geom_point(aes(x=x,y=f,shape=type,color=type),size=1)+
  theme_bw()+xlab(expression(x))+
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 3*sig_n[id.sig,nn],color='black',linetype='dashed')+
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

ggarrange(g1,g2,g3,nrow=1,ncol=3,common.legend = TRUE)

#--------------