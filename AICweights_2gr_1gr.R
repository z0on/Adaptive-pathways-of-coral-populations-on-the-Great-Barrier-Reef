# AIC weights of two 2d models: one growth, or two growths (in all cases, growth then decline).

#number of parameters:
npar=c(2,4)

nosplit=read.table("~/Documents/allRAD_august2015/dadi7/noBS_fdr05/1grBootlikes.txt")
nosplit$V3=as.numeric(as.character(nosplit$V3))
row.names(nosplit)=paste(nosplit$V1,nosplit$V2,sep=".")
nstheta=read.table("~/Documents/allRAD_august2015/dadi7/noBS_fdr05/1grThetas.txt")
row.names(nstheta)=paste(nstheta$V1,nstheta$V2,sep=".")
c1=row.names(nstheta)[nstheta$V3>1]
c2=row.names(nosplit)[!is.na(nosplit$V3)]
cc1=intersect(c1,c2)
names(nosplit)=c("dset","pop","like")

split=read.table("~/Documents/allRAD_august2015/dadi7/noBS_fdr05/2grBootlikes.txt")
split$V3=as.numeric(as.character(split$V3))
row.names(split)=paste(split$V1,split$V2,sep=".")
st=read.table("~/Documents/allRAD_august2015/dadi7/noBS_fdr05/2grBootTheta.txt")
row.names(st)=paste(st$V1,st$V2,sep=".")
c3=row.names(st)[st$V3>1]
c4=row.names(split)[!is.na(split$V3)]
cc2=intersect(c3,c4)
names(split)=c("dset","pop","like")

goods=intersect(cc1,cc2)

split=split[goods,]
nosplit=nosplit[goods,]

b.nosplit=nosplit[grep("boot",nosplit[,1]),]
b.split = split[grep("boot", split[,1]),]

#------------

npar=c(2,4)
likeCut=0.1
par(mfrow=c(2,3))
pops=c("W","S","O","M","K")
p1=5
for (p1 in 1:(length(pops))){
		pp1=pops[p1]
		s1=subset(b.nosplit,pop==pp1 )
		lc1=quantile(s1$like,prob=likeCut)
		s1=s1[s1$like>lc1,]
		s2=subset(b.split,pop==pp1)
		lc2=quantile(s2$like,prob=likeCut)
		s2=s2[s2$like>lc2,]
		goods.p=intersect(row.names(s1),row.names(s2))
		s1=s1[goods.p,]
		s2=s2[goods.p,]
#		lims=c(min(s2$like,s1$like),max(s2$like,s1$like))	
#		plot(s2$like~s1$like,ylim=lims,xlim=lims)
#		abline(0,1)
		aic1=2*npar[1]-2*s1$like
		aic2=2*npar[2]-2*s2$like
		delta.aic=aic1-aic2
		mac=median(delta.aic)
		sdac=sd(delta.aic)
#		delta.aic=delta.aic[delta.aic>(mac-10*mac) & delta.aic<(mac+10*mac)]
		print(paste(pops[p1],round(mac,1),signif(exp(mac/2),1),round(sum(delta.aic>0)/length(delta.aic),2)))
		hist(delta.aic,breaks=10,main=paste(pp1),mgp=c(2.3,1,0))
		mtext(round(100*sum(delta.aic>0)/length(delta.aic),0),cex=0.8)
}
# [1] "W 2.5 4 0.83"
# [1] "S 0.4 1 0.59"
# [1] "O 2.2 3 0.8"
# [1] "M -2.6 0.3 0.07"
# [1] "K 37.7 2e+08 1"
