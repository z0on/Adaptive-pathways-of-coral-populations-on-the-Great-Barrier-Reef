# AIC weights of two 2d models: one with split (s2m), another without (te2d).
# test for significance of population split

# see "detailed_walkthrough...": "combining all likelihoods for AIC-bootstrap" about how to create input files from dadi output

#---------------------

#number of parameters:
npar=c(2,5)

nosplit=read.table("te2d_allLikes.txt")
row.names(nosplit)=paste(nosplit[,1],nosplit[,2],nosplit[,3],sep=".")
names(nosplit)=c("dset","pop1","pop2","like")
nosplit$like=as.numeric(as.character(nosplit$like))
nosplit=nosplit[!(is.na(nosplit$like)),]

split=read.table("s2m_allLikes.txt")
row.names(split)=paste(split[,1], split[,2],split[,3],sep=".")
names(split)=c("dset","pop1","pop2","like")
split$like=as.numeric(as.character(split$like))
split=split[!(is.na(split$like)),]

goods=intersect(row.names(split),row.names(nosplit))
split=split[goods,]
nosplit=nosplit[goods,]

b.nosplit=nosplit[grep("boot",nosplit[,1]),]
b.split = split[grep("boot", split[,1]),]

#------------

npar=c(2,5)
likeCut=0.15
par(mfrow=c(4,3))
# plot(c(1:10)~1,col="white",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
# plot(c(1:10)~1,col="white",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
pops=c("W","S","O","M","K")
p1=1;p2=2
for (p1 in 1:(length(pops)-1)){
	for (p2 in (p1+1):(length(pops))){
		pp1=pops[p1];pp2=pops[p2]
		s1=subset(b.nosplit,pop1==pp1 & pop2==pp2)
		lc1=quantile(s1$like,prob=likeCut)
		s1=s1[s1$like>lc1,]
		s2=subset(b.split,pop1==pp1 & pop2==pp2)
		lc2=quantile(s2$like,prob=likeCut)
		s2=s2[s2$like>lc2,]
		goods.p=intersect(row.names(s1),row.names(s2))
		s1=s1[goods.p,]
		s2=s2[goods.p,]	
		aic1=2*npar[1]-2*s1$like
		aic2=2*npar[2]-2*s2$like
		delta.aic=aic1-aic2
		mac=median(delta.aic)
		sdac=sd(delta.aic)
		delta.aic=delta.aic[delta.aic>(mac-2*sdac) & delta.aic<(mac+2*sdac)]
		length(delta.aic)
		print(paste(pops[p1],pops[p2],round(mac,0),signif(exp(mac/2),1),round(sum(delta.aic>0)/length(delta.aic),2)))
		hist(delta.aic,breaks=10,main=paste(pp1,pp2),mgp=c(2.3,1,0))
		mtext(100*round(sum(delta.aic>0)/length(delta.aic),2),cex=0.8)
	}
}
# [1] "W S 12 500 0.98"
# [1] "W O 16 4000 1"
# [1] "W M 23 80000 1"
# [1] "W K 372 7e+80 1"
# [1] "S O 3 5 0.62"
# [1] "S M -1 0.6 0.43"
# [1] "S K 293 3e+63 1"
# [1] "O M 9 70 0.81"
# [1] "O K 384 2e+83 1"
# [1] "M K 205 4e+44 1"

