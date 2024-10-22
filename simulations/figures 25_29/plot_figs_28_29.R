
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
length_out = 50
bases = 10
mu0= 4
q=0.1
r = 42


#-----------------------------------------------------------------------------a

marg_npmleb_1 =  marg_dd_1 =  marg_orc_1 = array(0,c(1000,length(u)))
marg_npmleb_2 =  marg_dd_2 =  marg_orc_2 = array(0,c(1000,length(u)))
prior_hamt = prior_npmleb = array(0,c(length_out,n,length(u)))

fdp.dd = fdp.npmleb = rep(0,length(u))
etp.dd = etp.npmleb = rep(0,length(u))

clfdr.dd = clfdr.npmleb = de.dd = de.npmleb = array(0,c(n,length(u)))

sig_n = matrix(0,n,length(u))
gd_n = matrix(0,length_out,length(u))

for(nn in 1:length(u)){
  
  source('../funcs.R')
  set.seed(r)
  sig = runif(n,0.25,u[nn])
  mu = 3*sig
  x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
  theta=1*(mu>mu0)

  f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
  
  gd=seq(from=min(x),to=max(x),length.out=length_out)
  
  hamt_w = hamt_opt(x,sig,f_hat,gd,length_out = length_out,bases=bases)
  
  clfdr.dd[,nn] = clfdr.func(x,sig,mu0,hamt_w$probs,gd,type='proposed')
  decisions.dd = sc.func(clfdr.dd[,nn],q)
  fdp.dd[nn] = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
  etp.dd[nn]=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
  de.dd[,nn] = decisions.dd$de
  
  npmleb = npmleb_opt(x,sig,gd,length_out = length_out,bases=bases)
  
  clfdr.npmleb[,nn] = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                            type='proposed')
  decisions.npmleb = sc.func(clfdr.npmleb[,nn],q)
  fdp.npmleb[nn] = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
  etp.npmleb[nn]=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
  de.npmleb[,nn] = decisions.npmleb$de
  
  prior_hamt[,,nn] = hamt_w$probs
  prior_npmleb[,,nn] = npmleb$probs
  
 
  sig_n[,nn] = sig
  gd_n[,nn] = gd
  
}

#--------------------------------------------------------------------

save.image('plot_figs_28_29.RData')

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)

#-------------Figure 28 ----------------------
nn = 2
sig_sort = order(sig_n[,nn])
w_ordered = prior_hamt[,sig_sort,nn]
rownames(w_ordered)=gd_n[,nn]
data_melt <- melt(w_ordered)
i1 = which(sig_n[sig_sort,nn]>4/3)
i2 = which(sig_n[sig_sort,nn]<=4/3)
pointxy.1 <- data.frame(
  Var1 = 3*sig_n[sig_sort[i1],nn],
  Var2 = i1)
pointxy.2 <- data.frame(
  Var1 = 3*sig_n[sig_sort[i2],nn],
  Var2 = i2)
g1<-ggplot(data_melt, aes(Var1, Var2)) +
  geom_raster(aes(fill = value), interpolate = TRUE)+
  geom_line(data = pointxy.1,size=0.5)+
  geom_line(data = pointxy.2,size=0.5,col='green')+
  scale_fill_gradientn(colours=c("white","blue","red"))+
  guides(fill = guide_colourbar(title = "Probability",barwidth = 1,
                                barheight = 12))+
  scale_x_continuous(n.breaks = 10)+
  xlab(TeX("$T$"))+ylab(TeX("$\\sigma_i$"))+
  theme_bw()+ggtitle('HAMT')+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

w_ordered = prior_npmleb[,sig_sort,nn]
rownames(w_ordered)=gd_n[,nn]
data_melt <- melt(w_ordered)
g2<-ggplot(data_melt, aes(Var1, Var2)) +
  geom_raster(aes(fill = value), interpolate = TRUE)+
  geom_line(data = pointxy.1,size=0.5)+
  geom_line(data = pointxy.2,size=0.5,col='green')+
  scale_fill_gradientn(colours=c("white","blue","red"))+
  guides(fill = guide_colourbar(title = "Probability",barwidth = 1,
                                barheight = 12))+
  scale_x_continuous(n.breaks = 10)+
  xlab(TeX("$T$"))+ylab(TeX("$\\sigma_i$"))+
  theme_bw()+ggtitle('NPMLE B(S=50, K=10)')+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(g1,g2,ncol=2)

#--------------
#-------------Figure 29 ----------------------
nn = 3
sig_sort = order(sig_n[,nn])
w_ordered = prior_hamt[,sig_sort,nn]
rownames(w_ordered)=gd_n[,nn]
data_melt <- melt(w_ordered)
i1 = which(sig_n[sig_sort,nn]>4/3)
i2 = which(sig_n[sig_sort,nn]<=4/3)
pointxy.1 <- data.frame(
  Var1 = 3*sig_n[sig_sort[i1],nn],
  Var2 = i1)
pointxy.2 <- data.frame(
  Var1 = 3*sig_n[sig_sort[i2],nn],
  Var2 = i2)
g1<-ggplot(data_melt, aes(Var1, Var2)) +
  geom_raster(aes(fill = value), interpolate = TRUE)+
  geom_line(data = pointxy.1,size=0.5)+
  geom_line(data = pointxy.2,size=0.5,col='green')+
  scale_fill_gradientn(colours=c("white","blue","red"))+
  guides(fill = guide_colourbar(title = "Probability",barwidth = 1,
                                barheight = 12))+
  scale_x_continuous(n.breaks = 10)+
  xlab(TeX("$T$"))+ylab(TeX("$\\sigma_i$"))+
  theme_bw()+ggtitle('HAMT')+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

w_ordered = prior_npmleb[,sig_sort,nn]
rownames(w_ordered)=gd_n[,nn]
data_melt <- melt(w_ordered)
g2<-ggplot(data_melt, aes(Var1, Var2)) +
  geom_raster(aes(fill = value), interpolate = TRUE)+
  geom_line(data = pointxy.1,size=0.5)+
  geom_line(data = pointxy.2,size=0.5,col='green')+
  scale_fill_gradientn(colours=c("white","blue","red"))+
  guides(fill = guide_colourbar(title = "Probability",barwidth = 1,
                                barheight = 12))+
  scale_x_continuous(n.breaks = 10)+
  xlab(TeX("$T$"))+ylab(TeX("$\\sigma_i$"))+
  theme_bw()+ggtitle('NPMLE B(S=50, K=10)')+
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggarrange(g1,g2,ncol=2)

#--------------