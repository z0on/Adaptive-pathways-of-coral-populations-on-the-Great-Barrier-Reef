# see "detailed_walkthrough...": "analyzing one-pop growth bootstraps" about how to create input files from dadi output

#setwd("~/Documents/allRAD_august2015/dadi7/noBS_fdr05/")
library(ggplot2)

#---------------------

plotNuGG=function(pop,mutRate,genTime,Tmax,likesCut=0.1,CI=0.95,col="black") {
#pop="W";Tmax=500000;mutRate=0.018;genTime=5;likesCut=0.1;CI=0.95;col="hotpink";SD=FALSE
	# mutRate - number of detectable mutations per DIPLOID individual per generation
	require(ggplot2)
	para=read.table("2grBootParams.txt",sep=" ")
	th=read.table("2grBootTheta.txt",sep=" ")
	th$V3=as.numeric(as.character(th$V3))
	li=read.table("2grBootLikes.txt",sep=" ")
	for(i in 3:6) { para[,i]=as.numeric(as.character(para[,i])) }
	para=subset(para,V2==pop)
	row.names(para)=para$V1
	li=subset(li,V2==pop)
	row.names(li)=li$V1
	th=subset(th,V2==pop)
	row.names(th)=th$V1
	liCut=quantile(li[,3],prob=likesCut)
	check=apply(para[,3:6],1,mean)
	check2=row.names(th)[which(!is.na(th$V3))]
	check3=row.names(th)[which(th$V3>1)]
	check4=row.names(para)[which(!is.na(check) & check!=Inf)]
	check5=row.names(li)[which(li$V3>liCut)]
	check6=intersect(check4,check5)
	check7=intersect(check2,check3)
	Check=intersect(check6,check7)
	para=para[Check,]
	li=li[Check,]
	th=th[Check,]
	for(i in 3:6) { para[,i]=as.numeric(as.character(para[,i])) }
	Ne=th[,3]/(2*mutRate)
#	hist(log(para$V5))
#	hist(Ne)
	para[,3:6]=para[,3:6]*Ne
	para[,5:6]=para[,5:6]*2*genTime
	b=1
	nus=c()
#	hist(log(para$V6,10))
	for (b in 1:nrow(para)){
		nu1=para[b,3]
		nuF=para[b,4]
		Tg=para[b,5]
		Tf=para[b,6]
		nu=c()
		for (t in exp(seq(0,log(Tmax),length.out=50))){
			if (t>(Tg+Tf)) { 
				nu=append(nu,Ne[b]) 
				next
			} 
			if (t<(Tg+Tf) & t>Tf ) { 
				t0=Tg-(t-Tf)
				nu=append(nu,Ne[b]*exp(log(nu1/Ne[b]) * t0/Tg))
				next
			} else {  
				t0=Tf-t
				nu=append(nu,nu1*exp(log(nuF/nu1) * t0/Tf))
			} 
		}
		nus=data.frame(cbind(nus,nu))
	}
	time=1e-3*exp(seq(0,log(Tmax),length.out=50))
	mn=c();low=c();hi=c();sd=c();se=c()
	mn=apply(nus,1,mean,na.rm=TRUE)
#	plot(mn)
	sdd=apply(nus,1,sd,na.rm=TRUE)
#	se=sd/sqrt(ncol(nus))
	lop=(1-CI)/2
	hip=1-(1-CI)/2
	means=data.frame(cbind("mean"=mn,"time"=time))
#	plot(mean~time,means,log="x")
	means$lo=apply(nus,1,quantile,prob=lop)
	means$hi=apply(nus,1,quantile,prob=hip)
	ciribbon=aes(ymax=hi,ymin=lo)
	ggp=ggplot(data=means,aes(x=time,y=mean))+geom_line(color=col)+geom_ribbon(ciribbon,alpha=0.3,fill=col)+ylab("effective population size")
	return(list("plot"=ggp,"means"=means))
#	ggp+scale_x_log10()
}
#---------------------

tmax=7e+5;mutRate=0.018;genTime=5;lcut=0.1;CI=0.95 # mutRate= number of detectable mutations er diploid individual per generation = 2x dadi's mut rate

#plotNuGG("K",mutRate,genTime,tmax,tstep,lcut,CI=CI)+theme_bw()+ggtitle("Keppels")+xlab("time, KY")
library(gridExtra)
quartz()

mm=c()
pop=c()
times=c()

ppp=plotNuGG("W",mutRate,genTime,tmax,lcut,CI=CI)
Wi= ppp$plot+theme_bw()+ggtitle("Wilkie")+xlab("years before present")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+ylim(1000,56000)
ppp$means$pop="W"
mm=append(mm, ppp$means$mean)
pop=append(pop,rep("W",nrow(ppp$means)))
times=append(times,ppp$means$time)

ppp=plotNuGG("S",mutRate,genTime,tmax,lcut,CI=CI)
Su= ppp$plot+theme_bw()+ggtitle("Sudbury")+xlab("years before present")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+ylim(1000,56000)
ppp$means$pop="S"
mm=append(mm, ppp$means$mean)
pop=append(pop,rep("S",nrow(ppp$means)))
times=append(times,ppp$means$time)

ppp=plotNuGG("O",mutRate,genTime,tmax,lcut,CI=CI)
Or= ppp$plot+theme_bw()+ggtitle("Orpheus")+xlab("years before present")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+ylim(1000,56000)
ppp$means$pop="O"
mm=append(mm, ppp$means$mean)
pop=append(pop,rep("O",nrow(ppp$means)))
times=append(times,ppp$means$time)

ppp=plotNuGG("M",mutRate,genTime,tmax,lcut,CI=CI)
Ma= ppp$plot+theme_bw()+ggtitle("Magnetic")+xlab("years before present")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+ylim(1000,56000)
ppp$means$pop="M"
mm=append(mm, ppp$means$mean)
pop=append(pop,rep("M",nrow(ppp$means)))
times=append(times,ppp$means$time)

ppp=plotNuGG("K",mutRate,genTime,tmax,lcut,CI=CI)
Ke= ppp$plot+theme_bw()+ggtitle("Keppel")+xlab("years before present")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+ylim(1000,56000)
ppp$means$pop="K"
mm=append(mm, ppp$means$mean)
pop=append(pop,rep("K",nrow(ppp$means)))
times=append(times,ppp$means$time)

allmeans=data.frame(cbind("psize"=mm,"population"=pop,times))
allmeans$psize=as.numeric(as.character(allmeans$psize))
allmeans$times=as.numeric(as.character(allmeans$times))
allmeans$population=factor(allmeans$population,levels=c("W","S","O","M","K"))
ggplot(allmeans,aes(x=times,y=psize,colour=population))+geom_step(stat="identity")+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k","500k"),breaks=c(0.1,1,10,100,500))+xlab("years before present")+ylab("effective population size")+theme_bw()

slo=read.csv("seaLevel_O18.csv")
ciribbon.sl=aes(ymax=slo$sl.upper,ymin=slo$sl.lower)
ciribbon.o=aes(ymax=-slo$o18.upper,ymin=-slo$o18.lower)
quartz()
sl=ggplot(data=slo,aes(x=ka,y=sl))+geom_line(colour="blue")+geom_ribbon(ciribbon.sl,alpha=0.6,fill="skyblue")+ylab("sea level, m")+theme_bw()+scale_x_log10(limits=c(0.05,tmax/1000),labels=c("100","1k","10k","100k"),breaks=c(0.1,1,10,100))+xlab("years before present")

grid.arrange(Wi,Su,Or,Ma,Ke,sl,nrow=2)


