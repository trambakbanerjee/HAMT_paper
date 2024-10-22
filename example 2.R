
library(REBayes)
library(foreach)
library(doParallel)
require(ggplot2)
require(gridExtra)
library(ggpubr)
library(rmutil)
library(Rmosek)
library(ashr)


n<-10000        # number of hypotheses
u = 1.5         # parameter on the standard deviations
length_out = 50 # Grid size S
bases = 10      # no. of bases K
mu0= 1          # upper limit of the indifference region for one sided composite null
q=0.1           # FDR level

# --- A setting where $\mu$ is not independent of $\sigma$ ------ 
# \sigma ~ 0.9 Unif(0.5,1) + 0.1 Unif(1,u)
# \mu |\sigma_i = 0 (if sig[i]<=1), or 2/sig[i] (if sig[i]>1)
# X_i |\mu_i,\sigma_i ~ N(\mu_i,\sigma_i^2)
# \mathcal A = (-infinity, mu0]

source('funcs.R')
set.seed(42)
p = runif(n,0,1)
sig = sapply(1:n, function(i) runif(1,0.5,1)*(p[i]<=0.9)+runif(1,1,u)*(p[i]>0.9))
mu = sapply(1:n, function(i) 0*(sig[i]<=1)+
              (2/sig[i])*(sig[i]>1))
x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
theta=1*(mu>mu0)

################# HAMT ####################
f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
gd=seq(from=min(x),to=max(x),length.out=length_out)
hamt = hamt_opt(x,sig,f_hat,gd,length_out,bases)

clfdr.dd = clfdr.func(x,sig,mu0,hamt$probs,gd,
                      type='proposed')
decisions.dd = sc.func(clfdr.dd,q)
fdp.dd.r = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
etp.dd.r=sum((theta)*decisions.dd$de)/max(sum(theta), 1)

################### ASH (alpha = 0) #####################################

ash_out = ash(x,sig,mixcompdist = "normal",mode="estimate",
              outputlevel=3)

post_sample_ash=get_post_sample(ash_out,3000)
clfdr.ash = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=mu0))/3000)
decisions.ash = sc.func(clfdr.ash,q)
fdp.ash.r = sum((1-theta)*decisions.ash$de)/max(sum(decisions.ash$de),1)
etp.ash.r=sum((theta)*decisions.ash$de)/max(sum(theta), 1)
rm(post_sample_ash)

####################################### ASH (alpha = 1) ##################

ash_out = ash(x,sig,mixcompdist = "normal",alpha=1,outputlevel=3,mode="estimate")

post_sample_ash=get_post_sample(ash_out,3000)
clfdr.ash.1 = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=(mu0/sig[i])))/3000)
decisions.ash.1 = sc.func(clfdr.ash.1,q)
fdp.ash.1.r = sum((1-theta)*decisions.ash.1$de)/max(sum(decisions.ash.1$de),1)
etp.ash.1.r=sum((theta)*decisions.ash.1$de)/max(sum(theta), 1)
rm(post_sample_ash)

################# NPMLE with basis expansion ####################
m = length_out
K = bases
gd=seq(from=min(x),to=max(x),length.out=m)
npmle = npmleb_opt(x,sig,gd,length_out = m,bases=K)
g_npmle = npmle$probs

clfdr.npmleb = clfdr.func(x,sig,mu0,g_npmle,gd,
                          type='proposed')
decisions.npmleb = sc.func(clfdr.npmleb,q)
fdp.npmleb.r = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
etp.npmleb.r=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)


