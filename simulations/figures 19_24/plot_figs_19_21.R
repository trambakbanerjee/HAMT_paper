
#--------------------------------------------------------------------

load('fig7_marginal_dens.RData')

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)

#--------Figure 19 ---------------------------
nn = 1

#### Marginal density plot and prior probs at sig = 0.75
a= seq(-8,10,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.7))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.7)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 0.95
id.sig = which.min(abs(sig_n[,nn]-0.95))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.95)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-2/sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 2/sig_n[id.sig,nn],color='black',linetype='dashed')+
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
#--------Figure 20 ---------------------------
nn = 2

#### Marginal density plot and prior probs at sig = 0.75
a= seq(-8,10,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.7))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.7)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 0.95
id.sig = which.min(abs(sig_n[,nn]-0.95))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.95)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-2/sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 2/sig_n[id.sig,nn],color='black',linetype='dashed')+
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
#--------Figure 21 ---------------------------
nn = 3

#### Marginal density plot and prior probs at sig = 0.75
a= seq(-8,10,length.out=1000)
id.sig = which.min(abs(sig_n[,nn]-0.7))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.7)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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

#### Marginal density plot and prior probs at sig = 0.95
id.sig = which.min(abs(sig_n[,nn]-0.95))
marg_npmleb_1 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_npmleb_1[,id.sig,nn]))
marg_npmleb_2 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_2[,nn],0,sig_n[id.sig,nn])*prior_npmleb_2[,id.sig,nn]))
marg_npmleb_3 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_3[,nn],0,sig_n[id.sig,nn])*prior_npmleb_3[,id.sig,nn]))
marg_npmleb_4 = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_4[,nn],0,sig_n[id.sig,nn])*prior_npmleb_4[,id.sig,nn]))
marg_dd = sapply(1:length(a),function(i) sum(dnorm(a[i]-gd_1[,nn],0,sig_n[id.sig,nn])*prior_hamt[,id.sig,nn]))
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-0,0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = 0.95)$"))+geom_vline(xintercept = 0,color='black',linetype='dashed')+
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
marg_orc = sapply(1:length(a),function(i) dnorm(a[i]-2/sig_n[id.sig,nn],0,sig_n[id.sig,nn]))

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
  ylab(TeX("$f(x|\\sigma = \\bar{\\sigma}-0.1)$"))+geom_vline(xintercept = 2/sig_n[id.sig,nn],color='black',linetype='dashed')+
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