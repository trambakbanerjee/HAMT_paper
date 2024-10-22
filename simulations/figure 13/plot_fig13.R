
require(ggplot2)
library(latex2exp)
require(gridExtra)
library(ggpubr)

# Setting 2: Logistic likelihood (two-sided)

load("fig13_set2.RData")
out_exp1_others = out_exp1
load("fig13_set2_adapt.RData")

plotdata1<- as.data.frame(c(colMeans(out_exp1$fdr.adapt,na.rm=TRUE),
                            colMeans(out_exp1_others$fdr.dd,na.rm=TRUE),
                            colMeans(out_exp1_others$fdr.or),
                            colMeans(out_exp1_others$fdr.npmleb)))
names(plotdata1)<-"fdr"
plotdata1$type<-as.factor(rep(c('AdaPTGMM','HAMT',
                                'OR','NPMLE B'),each=length(u)))
plotdata1$n<- rep(u,4)
g1<-ggplot()+geom_line(data=plotdata1,aes(x=n,y=fdr,color=type),size=1.5)+
  geom_point(data=plotdata1,aes(x=n,y=fdr,color=type,shape=type),size=4)+
  scale_shape_manual(values=1:nlevels(plotdata1$type)) +
  geom_hline(yintercept = 0.1,linetype='dotted', size=2,col = 'black')+
  scale_y_continuous(limits = c(0,0.2))+
  ylab('FDP')+theme_bw()+xlab(TeX("$\\bar{\\sigma}$"))+
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

plotdata2<- as.data.frame(c(colMeans(out_exp1$etp.adapt,na.rm=TRUE),
                            colMeans(out_exp1_others$etp.dd,na.rm=TRUE),
                            colMeans(out_exp1_others$etp.or),
                            colMeans(out_exp1_others$etp.npmleb)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('AdaPTGMM','HAMT',
                                'OR','NPMLE B'),each=length(u)))
plotdata2$n<- rep(u,4)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=n,y=etp,color=type),size=1.5)+
  geom_point(data=plotdata2,aes(x=n,y=etp,color=type,shape=type),size=4)+
  scale_shape_manual(values=1:nlevels(plotdata1$type))+
  ylab('PTP')+theme_bw()+xlab(TeX("$\\bar{\\sigma}$"))+
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
