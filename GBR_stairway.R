stfa=read.table("~/Dropbox/GBR_pooledStairway/GBR.final.summary",sep="\t",header=T)
library(ggplot2)
head(stfa)
plot(stfa$mutation_per_site)
ciribbon75=aes(ymax= stfa$Ne_87.5.,ymin= stfa$Ne_12.5.)
ciribbon95=aes(ymax= stfa$Ne_97.5.,ymin= stfa$Ne_2.5.)

stfa$year2= stfa$year*2

#ggplot(data=stfa,aes(x=year,y=Ne_median))+geom_line()+geom_ribbon(ciribbon95,alpha=0.15)+geom_ribbon(ciribbon75,alpha=0.2)+scale_x_log10(limits=c(10,1e+6),labels=c("10","100","1k","10k","100k","1M"),breaks=c(10,1e+2,1e+3,1e+4,1e+5,1e+6))+xlab("years before present")+ylab("effective population size")+scale_y_log10(breaks=c(1e+4,3e+4,1e+5,3e+5),labels=c(1e+4,3e+4,1e+5,3e+5))+theme_bw()

st=read.table("~/Dropbox/GBR_pooledStairway/KEP.final.summary",sep="\t",header=T)
st$year2=st$year*2
head(st)
stfa$population="GBR"
st$population="Keppel"
sts=rbind(stfa,st)
ciribbon75=aes(ymax=sts$Ne_87.5.,ymin=sts$Ne_12.5.)
ciribbon95=aes(ymax=sts$Ne_97.5.,ymin=sts$Ne_2.5.)

str(sts)
ggplot(data=sts,aes(x=year,y=Ne_median,group=population,fill= population))+geom_line(aes(colour= population))+geom_ribbon(ciribbon95,alpha=0.15)+geom_ribbon(ciribbon75,alpha=0.2)+scale_x_log10(limits=c(10,1.3e+6),labels=c("100","1k","10k","100k","1M"),breaks=c(1e+2,1e+3,1e+4,1e+5,1e+6))+xlab("years before present")+ylab("effective population size")+scale_y_log10(breaks=c(1e+4,3e+4,1e+5,3e+5),labels=c(1e+4,3e+4,1e+5,3e+5))+theme_bw()
