# extracting mutation data from slim output
# ls e*out | perl -pe 's/(\S+)/grep \"#OUT\" $1 \| perl -pe \"s\/#OUT\:\\s\/\/\" >$1\.muts/' | bash

par(mfrow=c(3,3))

infiles=c(
"e0s0.5a.out.muts",
"e1s0.5a.out.muts",
"e2s0.5a.out.muts",
"e0a.out.muts",
"e1a.out.muts",
"e2a.out.muts",
"e0s2a.out.muts",
"e1s2a.out.muts",
"e2s2a.out.muts",
"e0s0.5b.out.muts",
"e1s0.5b.out.muts",
"e2s0.5b.out.muts",
"e0b.out.muts",
"e1b.out.muts",
"e2b.out.muts",
"e0s2b.out.muts",
"e1s2b.out.muts",
"e2s2b.out.muts"
)
colors=c("firebrick","khaki3","green3","coral","skyblue")
alpha=0.15
#colors=c(rgb(0.7,0,0,alpha=alpha),rgb(0,0.7,0,alpha=alpha),rgb(0,1,0,alpha=alpha),rgb(1,0,0,alpha=alpha),rgb(.3,0.3,1,alpha=alpha))
infi=infiles[1]
nmuts=c();ceff=c()
for (infi in infiles) {
	mm=read.table(infi)
	mm=mm[,-c(1,2,5,6,8,10)]
	names(mm)=c("pop","id","beta","pop.ori","n")
	cols=rep(colors[1],nrow(mm))
	p=1
	nm=c();ce=c()
#	str(mm)
	for (p in 1:length(levels(mm$pop))) {
		Pop=levels(mm$pop)[p]
		s=subset(mm,pop==Pop)
		cols[which(mm$pop==levels(mm$pop)[p])]=colors[p]
		nm=append(nm,sum(mm$pop==Pop))
		ce=append(ce,round(sum(abs(s$beta)*s$n),0))
		print(paste(infi,p,Pop,sum(mm$pop==Pop),round(sum(abs(s$beta)*s$n),0)))
	}
	nmuts=data.frame(cbind(nmuts,nm))
	ceff=data.frame(cbind(ceff,ce))	
	par(mai=c(0.3,0.3,0.2,0.2))
	p=plot(n~beta,mm,ylim=c(1,11000),xlim=c(-0.8,0.8),log="y",cex=0.2,pch=19,col=cols,mgp=c(2.1,1,0),xlab=NA,ylab=NA)
}
nmuts
ceff

# plotting legend
quartz()
plot(n~beta,mm,ylim=c(1,11000),xlim=c(-0.8,0.8),log="y",cex=0.3,pch=19,col="white",mgp=c(2.1,1,0),xlab=NA,ylab=NA)
alpha=1
#colors=c(rgb(0.7,0,0,alpha=alpha),rgb(0,0.7,0,alpha=alpha),rgb(0,1,0,alpha=alpha),rgb(1,0,0,alpha=alpha),rgb(.3,0.3,1,alpha=alpha))
legend("topleft",legend=c("W","S","O","M","K"),col=colors, pch=19)

nmuts=nmuts[-5,]
snm=stack(nmuts)
snm$e=c(rep("0",4),rep("1",4),rep("2",4))
snm$s=c(rep("0.5",12),rep("1",12),rep("2",12))
snm$env=c("hot","cool","cool","hot")
snm=subset(snm,env=="cool")
#snm$h2=factor(snm$h2,levels=c("1","0.5","0.25"))
library(ggplot2)
pd=position_dodge(0)
quartz()
ggplot(snm,aes(e,values))+geom_jitter(width=0.1)+theme_bw()+ylab("number of mutations")+theme(legend.position="top")
#ggplot(snm,aes(h2,values,colour=s))+geom_boxplot(position=pd,width=1.5)+theme_bw()+ylab("number of mutations")
ggplot(snm,aes(e,values,colour=s))+geom_boxplot(width=1)+theme_bw()+ylab("number of mutations")
+theme(legend.position="top")

ceff=ceff[-5,]
snm=stack(ceff)
snm$e=c(rep("0",4),rep("1",4),rep("2",4))
snm$s=c(rep("0.5",12),rep("1",12),rep("2",12))
snm$env=c("hot","cool","cool","hot")
snm=subset(snm,env=="cool")
#snm$h2=factor(snm$h2,levels=c("1","0.5","0.25"))
library(ggplot2)
pd=position_dodge(0)
quartz()
ggplot(snm,aes(e,values))+geom_jitter(width=0.1)+theme_bw()+ylab("number of mutations")+theme(legend.position="top")
#ggplot(snm,aes(h2,values,colour=s))+geom_boxplot(position=pd,width=1.5)+theme_bw()+ylab("number of mutations")
ggplot(snm,aes(e,values,colour=s))+geom_boxplot(width=0.7)+theme_bw()+ylab("cumulative effect")
+theme(legend.position="top")

