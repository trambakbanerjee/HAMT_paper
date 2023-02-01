

####### simple example ###########################
q=0.1
pii = 0.1
del = 1.5
a = 0.5
b = 4

t_seq =seq(0,10,length.out=1000)
val = matrix(0,length(t_seq))
for(i in 1:length(t_seq)){
  
    integrand <- function(sig,t,d) {1-pnorm(t-sig^{d-1})}
    val[i] = ((1-pii)*(1-pnorm(t_seq[i]))/((1-pii)*(1-pnorm(t_seq[i]))+
                                    (pii/(b-a))*integrate(integrand,lower=a,
                                            upper=b,t=t_seq[i],d=del)$value))-q
}
tmin<-t_seq[min(which(val<=0))]
integrand <- function(sig,pii,d) {1-pnorm(tmin-sig^{d-1})}
power_standardized = (1/(b-a))*integrate(integrand,lower=a,
                                       upper=b,pii=0.1,d=del)$value

lam_seq = seq(0.001,0.999,length.out=1000)
val = matrix(0,length(lam_seq))
for(i in 1:length(lam_seq)){
  
  int.1<-function(sig,pii,lam,d){
    
    1-pnorm((0.5*sig^{2*(d-1)}-log((lam*pii)/((1-lam)*(1-pii))))/sig^{d-1})
  }
  int.2<-function(sig,pii,lam,d){
    
    1-pnorm((0.5*sig^{2*(d-1)}-log((lam*pii)/((1-lam)*(1-pii))))/sig^{d-1}-sig^{d-1})
  }
  val.1 = integrate(int.1, lower = a, upper =b,pii=0.1,lam=lam_seq[i],d=del)$value/(b-a)
  val.2 = integrate(int.2, lower = a, upper =b,pii=0.1,lam=lam_seq[i],d=del)$value/(b-a)
  val[i] = ((1-pii)*val.1/((1-pii)*val.1+pii*val.2))-q
  
}
lam_star<-lam_seq[which.max(val[val<=0])+1]
int.2<-function(sig,pii,lam,d){
  
  1-pnorm((0.5*sig^{2*(d-1)}-log((lam*pii)/((1-lam)*(1-pii))))/sig^{d-1}-sig^{d-1})
}
power_orc = integrate(int.2, lower = a, upper =b,pii=0.1,lam=lam_star,d=del)$value/(b-a)

n = 10000
set.seed(1)
sig = runif(n,a,b)
p = runif(n,0,1)
mu = sapply(1:n, function(i) (p[i]<=0.9)*0+(p[i]>0.9)*(sig[i])^{1.5})
x = sapply(1:n,function(i) rnorm(1,mu[i],sig[i]))

z = x/sig
temp = (lam_star/(1-lam_star))*(0.1/0.9)
thresh2 = (1/sig^{del-1})*(-log(temp)+0.5*sig^{2*del-2})

library(ggplot2)
library(ggpubr)
library(latex2exp)
plotdata1 = data.frame('sig'=sig,'z'=z,'thresh2'=thresh2,'p'=p,'tmin'=rep(tmin,n))
ss = sig[which.min(abs(tmin-thresh2))]
x1_z = c(seq(0.5,ss,length.out=50))
y1_z = rep(tmin,50)
x2_z = c(seq(0.5,ss,length.out=50))
y2_z = thresh2[sapply(1:50,function(i) which.min(abs(sig-x1_z[i])))]
x1_or = c(seq(ss,4,length.out=50))
y1_or = thresh2[sapply(1:50,function(i) which.min(abs(sig-x1_or[i])))]
x2_or = c(seq(ss,4,length.out=50))
y2_or = rep(tmin,50)
x1_none = seq(0.5,4,length.out=400)
y1_none = rep(2.5,400)
x2_none = x1_none
temp_none = sapply(1:400,function(i) which.min(abs(sig-x1_none[i])))#c(seq(ss,4,length.out=200))
y2_none = pmin(tmin,thresh2[temp_none])
ggplot()+geom_line(data=plotdata1,aes(x=sig,y=thresh2),col='blue',size=1.5)+
  geom_line(data=plotdata1,aes(x=sig,y=tmin),col='red',size=1.5)+xlab(TeX("$\\sigma$"))+
  ylab(TeX('$Z$'))+ylim(2.5,6)+
  geom_segment(aes(x = x1_z, y = y1_z, xend = x2_z, yend = y2_z),col='red')+
  geom_segment(aes(x = x1_or, y = y1_or, xend = x2_or, yend = y2_or),col='blue')+
  geom_segment(aes(x = x1_none, y = y1_none, xend = x2_none, yend = y2_none),col='gray')+
  geom_point(data=plotdata1[plotdata1$p>0.9 & plotdata1$z>=pmin(thresh2,tmin),],aes(x=sig,y=z),col='black')+
  theme_bw()+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
######################################################################################
######################################## Another example #################################
q=0.1
pii = 0.1
del = 1.5
a = 0.5
b = 4
a0 = 0.9*b+a*(1-0.9) # 90th percentile of U(a,b)
t_seq =seq(0,10,length.out=1000)
val = matrix(0,length(t_seq))
for(i in 1:length(t_seq)){
  
  integrand <- function(sig,t,d) {1-pnorm(t-sig^{d-1})}
  val[i] = ((1-pii)*(1-pnorm(t_seq[i]))/((1-pii)*(1-pnorm(t_seq[i]))+
                                           (pii/(b-a))*integrate(integrand,lower=a0,
                                                               upper=b,t=t_seq[i],d=del)$value))-q
}
tmin<-t_seq[min(which(val<=0))]
integrand <- function(sig,pii,d) {1-pnorm(tmin-sig^{d-1})}
power_standardized = (1/(b-a))*integrate(integrand,lower=a0,
                                       upper=b,pii=0.1,d=del)$value


#######################################################################################
n = 10000
set.seed(1)
sig = runif(n,a,b)
mu = sapply(1:n, function(i) (sig[i]<=a0)*0+(sig[i]>a0)*(sig[i])^{1.5})
x = sapply(1:n,function(i) rnorm(1,mu[i],sig[i]))

z = x/sig
thresh2 = rep(a0,n)
plotdata2 = data.frame('sig'=sig,'z'=z,'thresh2'=thresh2,'tmin'=rep(tmin,n))
ss = a0
x1_z = c(seq(3,ss-0.01,length.out=50))
y1_z = rep(tmin,50)
x2_z = c(seq(3,ss-0.01,length.out=50))
y2_z = rep(max(z)+1,50)
x1_or = c(seq(ss,4,length.out=50))#rep(a0,50)
y1_or = rep(min(z[sig>a0]),50)#seq(min(z[sig>a0]),tmin,length.out=50)
x2_or = c(seq(ss,4,length.out=50))#rep(4,50)
y2_or = rep(tmin,50)#y1_or
x1_none = c(seq(3,ss-0.005,length.out=200))
y1_none = rep(min(z[sig>a0]),200)
x2_none = x1_none
y2_none = rep(tmin,200)
ggplot()+geom_line(data=plotdata2[plotdata2$sig>=3,],aes(x=sig,y=tmin),col='red',size=1.5)+
  geom_segment(aes(x = a0,y=min(z[sig>a0]),xend = a0, yend = max(z)+1),col='blue',size=1.5)+
  geom_segment(aes(x = x1_z, y = y1_z, xend = x2_z, yend = y2_z),col='red')+
  geom_segment(aes(x = x1_or, y = y1_or, xend = x2_or, yend = y2_or),col='blue')+
  geom_segment(aes(x = x1_none, y = y1_none, xend = x2_none, yend = y2_none),col='gray')+
  geom_point(data=plotdata2[plotdata2$sig>a0,],
             aes(x=sig,y=z),col='black')+xlab(TeX("$\\sigma$"))+
  ylab(TeX('$Z$'))+ylim(min(z[sig>a0]),max(z)+1)+
  theme_bw()+
  theme(legend.position='top',legend.title=element_blank(),
        legend.background = element_rect(fill="white",
                                         size=1, linetype="solid", colour ="black"),
        legend.text=element_text(size=12),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


