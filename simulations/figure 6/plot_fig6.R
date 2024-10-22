
require(ggplot2)
library(latex2exp)
require(gridExtra)
library(ggpubr)



# Figure 6 Setting 3

load("fig6_set3.RData")
out_exp1_others = out_exp1
load("fig6_set3_ash.RData")
out_exp1_ashes = out_exp1
load("fig6_set3_adapt.RData")
out_exp1_adapt = out_exp1
load("fig6_set3_camt_ihw.RData")

plotdata1<- as.data.frame(c(colMeans(out_exp1$fdr.camt,na.rm=TRUE),
                            colMeans(out_exp1$fdr.ihw,na.rm=TRUE),
                            colMeans(out_exp1_adapt$fdr.adapt,na.rm=TRUE),
                            colMeans(out_exp1_others$fdr.bh),
                            colMeans(out_exp1_others$fdr.npmle1),
                            colMeans(out_exp1_others$fdr.npmle2),
                            colMeans(out_exp1_others$fdr.deconv,na.rm=TRUE),
                            colMeans(out_exp1_others$fdr.dd,na.rm=TRUE),
                            colMeans(out_exp1_others$fdr.or),
                            colMeans(out_exp1_ashes$fdr.ash),
                            colMeans(out_exp1_ashes$fdr.ash.1),
                            colMeans(out_exp1_ashes$fdr.npmleb)))
names(plotdata1)<-"fdr"
plotdata1$type<-as.factor(rep(c('CAMT','IHW','AdaPTGMM','BH','GS 1','GS 2','DECONV','HAMT',
                                'OR','ASH','ASH 1','NPMLE B'),each=length(u)))
plotdata1$n<- rep(u,12)
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

plotdata2<- as.data.frame(c(colMeans(out_exp1$etp.camt,na.rm=TRUE),
                            colMeans(out_exp1$etp.ihw,na.rm=TRUE),
                            colMeans(out_exp1_adapt$etp.adapt,na.rm=TRUE),
                            colMeans(out_exp1_others$etp.bh),
                            colMeans(out_exp1_others$etp.npmle1),
                            colMeans(out_exp1_others$etp.npmle2),
                            colMeans(out_exp1_others$etp.deconv,na.rm=TRUE),
                            colMeans(out_exp1_others$etp.dd,na.rm=TRUE),
                            colMeans(out_exp1_others$etp.or),
                            colMeans(out_exp1_ashes$etp.ash),
                            colMeans(out_exp1_ashes$etp.ash.1),
                            colMeans(out_exp1_ashes$etp.npmleb)))
names(plotdata2)<-"etp"
plotdata2$type<-as.factor(rep(c('CAMT','IHW','AdaPTGMM','BH','GS 1','GS 2','DECONV','HAMT',
                                'OR','ASH','ASH 1','NPMLE B'),each=length(u)))
plotdata2$n<- rep(u,12)
g2<-ggplot()+geom_line(data=plotdata2,aes(x=n,y=etp,color=type),size=1.5)+
  geom_point(data=plotdata2,aes(x=n,y=etp,color=type,shape=type),size=4)+
  scale_shape_manual(values=1:nlevels(plotdata2$type))+
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

