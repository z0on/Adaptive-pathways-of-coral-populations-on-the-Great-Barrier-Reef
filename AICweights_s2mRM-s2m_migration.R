# AIC weights of two 2d models: s2m, and s2m with recent (T=0.01) migration change 

#number of parameters:
npar=c(5,7)

oneM=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2m_boots/s2m_allLikes2.txt")
row.names(oneM)=paste(oneM[,1],oneM[,2],oneM[,3],sep=".")
names(oneM)=c("dset","pop1","pop2","like")
oneM$like=as.numeric(as.character(oneM$like))
oneM=oneM[!(is.na(oneM$like)),]

twoM=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/s2mRM_allLikes.txt")
row.names(twoM)=paste(twoM[,1], twoM[,2],twoM[,3],sep=".")
names(twoM)=c("dset","pop1","pop2","like")
twoM$like=as.numeric(as.character(twoM$like))
twoM=twoM[!(is.na(twoM$like)),]

goods=intersect(row.names(twoM),row.names(oneM))
twoM=twoM[goods,]
oneM=oneM[goods,]
nrow(twoM)
b.oneM=oneM[grep("boot",oneM[,1]),]
b.twoM = twoM[grep("boot", twoM[,1]),]

#------------

npar=c(5,7)
likeCut=0.25
par(mfrow=c(3,4))
pops=c("W","S","O","M","K")
p1=1;p2=2
for (p1 in 1:(length(pops)-1)){
	for (p2 in (p1+1):(length(pops))){
		pp1=pops[p1];pp2=pops[p2]
		s1=subset(b.oneM,pop1==pp1 & pop2==pp2)
		lc1=quantile(s1$like,prob=likeCut)
		s1=s1[s1$like>lc1,]
		s2=subset(b.twoM,pop1==pp1 & pop2==pp2)
		lc2=quantile(s2$like,prob=likeCut)
		s2=s2[s2$like>lc2,]
		goods.p=intersect(row.names(s1),row.names(s2))
		s1=s1[goods.p,]
		s2=s2[goods.p,]	
		aic1=2*npar[1]-2*s1$like
		aic2=2*npar[2]-2*s2$like
		delta.aic=aic1-aic2
		delta.aic=delta.aic[abs(delta.aic)<30]
		mac=median(delta.aic)
		sdac=sd(delta.aic)
		delta.aic=delta.aic[delta.aic>(mac-2*sdac) & delta.aic<(mac+2*sdac)]
		print(paste(pops[p1],pops[p2],round(mac,0),signif(exp(mac/2),1),round(sum(delta.aic>0)/length(delta.aic),2)))
		length(delta.aic)
		hist(delta.aic,breaks=10,main=paste(pp1,pp2),mgp=c(2.3,1,0))
	}
}

# [1] "W S -2 0.4 0.22"
# [1] "W O 4 6 0.85"
# [1] "W M 0 1 0.6"
# [1] "W K -3 0.3 0.02"
# [1] "S O -1 0.6 0.44"
# [1] "S M -1 0.7 0.41"
# [1] "S K -3 0.2 0"
# [1] "O M -2 0.4 0.18"
# [1] "O K 0 0.9 0.46"
# [1] "M K -2 0.4 0.22"


# boxplot of recent change in migration rates: WO
wo=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/WO_s2mRM_params.txt")
row.names(wo)=wo[,1]
wot=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/WO_s2mRM_theta.txt")
row.names(wot)=wot[,1]
wol=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/WO_s2mRM_logLike.txt")
row.names(wol)=wol[,1]
lc=quantile(wol$V4,prob=likeCut)
goods=row.names(wol)[wol$V4>lc]
wo=wo[goods,]
wot=wot[goods,]
mutRate=0.64
Ne=round(wot[,4]/(2*mutRate),0)
for (i in 5:8) { wo[,i]=wo[,i]/(2*Ne) }
names(wo)[5:8]=c("N","S","N:recent","S:recent")
wos=stack(wo[,5:8])
wos$ind=factor(wos$ind,levels=c("N","S","N:recent","S:recent"))
boxplot(values~ind,wos,main="Wilkie-Orpheus",xlab="direction : time period",ylab="immigration rate",mgp=c(2.1,1,0))

# boxplot of recent change in migration rates: WM
wo=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/MK_s2mRM_params.txt")
row.names(wo)=wo[,1]
wot=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/MK_s2mRM_theta.txt")
row.names(wot)=wot[,1]
wol=read.table("~/Documents/allRAD_august2015/digitifera/dadi5/s2mRM/MK_s2mRM_logLike.txt")
row.names(wol)=wol[,1]
lc=quantile(wol$V4,prob=0.25)
goods=row.names(wol)[wol$V4>lc]
wo=wo[goods,]
wot=wot[goods,]
mutRate=0.018
Ne=round(wot[,4]/(2*mutRate),0)
for (i in 5:8) { wo[,i]=wo[,i]/(2*Ne) }
names(wo)[5:8]=c("N","S","Nr","Sr")
wos=stack(wo[,5:8])
wos$ind=factor(wos$ind,levels=c("N","S","Nr","Sr"))
wos$logm=log(wos$values,10)
boxplot(logm~ind,wos,main="Magnetic-Keppel",xlab="direction : time period",ylab="log10(immigration rate)",mgp=c(2.1,1,0))
