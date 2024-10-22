

library(foreach)
library(doParallel)
library(ashr)
library(Rmosek)

exp_1<-function(n,mu0,u,q,reps){
  
  fdr.ash<-matrix(0,reps,length(u))
  etp.ash<-matrix(0,reps,length(u))
  fdr.ash.1<-matrix(0,reps,length(u))
  etp.ash.1<-matrix(0,reps,length(u))
  fdr.npmleb<-matrix(0,reps,length(u))
  etp.npmleb<-matrix(0,reps,length(u))
  fdr.dd<-matrix(0,reps,length(u))
  etp.dd<-matrix(0,reps,length(u))
 
  for(nn in 1:length(u)){
    
    cl <- makeCluster(20)
    registerDoParallel(cl)
    
    result<-foreach(r = 1:reps,.packages=c("ashr","Rmosek"))%dopar%{
      source("../funcs.R")
      set.seed(r)
      sig = runif(n,0.25,u[nn])
      p = runif(n,0,1)
      mu = sapply(1:n, function(i) 0*(p[i]<=0.8)+
                    (p[i]>0.8)*rnorm(1,0,4))
      x=sapply(1:n, function(i)rnorm(1,mu[i],sig[i]))
      theta=1*(mu>mu0)
      
      ash_out = ash(x,sig,mixcompdist = "normal",outputlevel=3,
                    mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=mu0))/3000)
      decisions.ash = sc.func(clfdr.ash,q)
      fdp.ash.r = sum((1-theta)*decisions.ash$de)/max(sum(decisions.ash$de),1)
      etp.ash.r=sum((theta)*decisions.ash$de)/max(sum(theta), 1)
      
      ####################################### alpha = 1 ##################
      ash_out = ash(x,sig,mixcompdist = "normal",alpha=1,outputlevel=3,
                    mode="estimate")
      post_sample_ash=get_post_sample(ash_out,3000)
      clfdr.ash.1 = sapply(1:n,function(i) sum(1*(post_sample_ash[,i]<=(mu0/sig[i])))/3000)
      decisions.ash.1 = sc.func(clfdr.ash.1,q)
      fdp.ash.1.r = sum((1-theta)*decisions.ash.1$de)/max(sum(decisions.ash.1$de),1)
      etp.ash.1.r=sum((theta)*decisions.ash.1$de)/max(sum(theta), 1)
      rm(post_sample_ash)
      
      f_hat = nest.func.fast(x,sig,rep(1/length(x),length(x)),n)
      gd=seq(from=min(x),to=max(x),length.out=50)
      hamt_w = hamt_opt(x,sig,f_hat,gd,length_out = 50,bases=10)
      clfdr.dd = clfdr.func(x,sig,mu0,hamt_w$probs,gd,type='proposed')
      decisions.dd = sc.func(clfdr.dd,q)
      fdp.dd = sum((1-theta)*decisions.dd$de)/max(sum(decisions.dd$de),1)
      etp.dd=sum((theta)*decisions.dd$de)/max(sum(theta), 1)
      
      ################# NPMLE with basis expansion ####################
      npmleb = npmleb_opt(x,sig,gd,length_out = 50,bases=10)
      
      clfdr.npmleb = clfdr.func(x,sig,mu0,npmleb$probs,gd,
                                     type='proposed')
      decisions.npmleb = sc.func(clfdr.npmleb,q)
      fdp.npmleb = sum((1-theta)*decisions.npmleb$de)/max(sum(decisions.npmleb$de),1)
      etp.npmleb=sum((theta)*decisions.npmleb$de)/max(sum(theta), 1)
    
      return(list("fdr.ash_r"=fdp.ash.r,
                  "fdr.ash.1_r"=fdp.ash.1.r,
                  "fdr.npmleb_r"=fdp.npmleb,
                  "fdr.dd"=fdp.dd,
                  "etp.ash_r"=etp.ash.r,
                  "etp.ash.1_r"=etp.ash.1.r,
                  "etp.npmleb_r"=etp.npmleb,
                  "etp.dd"=etp.dd))
    }
    stopCluster(cl)
    registerDoSEQ()
    
    fdr.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash_r)
    etp.ash[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash_r)
    fdr.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.ash.1_r)
    etp.ash.1[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.ash.1_r)
    fdr.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.npmleb_r)
    etp.npmleb[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.npmleb_r)
    fdr.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$fdr.dd)
    etp.dd[,nn]<-sapply(1:reps,function(i) result[[i]]$etp.dd)
    print(nn)
  }
  result<-list("fdr.ash"=fdr.ash,
               "fdr.ash.1"=fdr.ash.1,
               "fdr.npmleb"=fdr.npmleb,
               "fdr.dd"=fdr.dd,
               "etp.ash"=etp.ash,
               "etp.ash.1"=etp.ash.1,
               "etp.npmleb"=etp.npmleb,
               "etp.dd"=etp.dd)
}
#--------------------------------------------------------------------

n<-10000
u = c(1,1.1,1.2,1.3,1.4,1.5)
mu0 = 1
q=0.1
reps=200

out_exp1<-exp_1(n,mu0,u,q,reps)

save.image('fig31.RData')

library(ggplot2)
library(ggpubr)
library()
plotdata1<- as.data.frame(c(colMeans(out_exp1$fdr.dd),
                            colMeans(out_exp1$fdr.ash),
                            colMeans(out_exp1$fdr.ash.1),
                            colMeans(out_exp1$fdr.npmleb)))
names(plotdata1)<-"fdr"
plotdata1$type<-as.factor(rep(c('HAMT','ASH','ASH 1','NPMLE B'),each=length(u)))
plotdata1$n<- rep(u,4)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=n,y=fdr,color=type),size=1.5)+
  geom_point(data=plotdata1,aes(x=n,y=fdr,color=type,shape=type),size=4)+
  scale_shape_manual(values=1:nlevels(plotdata1$type)) +
  geom_hline(yintercept = 0.1,linetype='dotted', size=2,col = 'black')+
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

plotdata2<- as.data.frame(c(colMeans(out_exp1$etp.dd),
                            colMeans(out_exp1$etp.ash),
                            colMeans(out_exp1$etp.ash.1),
                            colMeans(out_exp1$etp.npmleb)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('HAMT','ASH','ASH 1','NPMLE B'),each=length(u)))
plotdata2$n<- rep(u,4)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=n,y=etp,color=type),size=1.5)+
  geom_point(data=plotdata2,aes(x=n,y=etp,color=type,shape=type),size=4)+
  scale_shape_manual(values=1:nlevels(plotdata1$type))+
  ylab('PTP')+theme_bw()+xlab('u')+
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