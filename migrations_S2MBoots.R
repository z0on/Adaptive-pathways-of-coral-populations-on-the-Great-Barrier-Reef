# see "detailed_walkthrough...": "collecting results: s2m" section about how to create input files from dadi output


setwd("/Users/c-monstr/Documents/allRAD_august2015/dadi7/noBS_fdr05")
pops=c("W","S","O","M","K")

likesCut=0.15;mutRate=0.018;gentime=5
popp1=c();popp2=c();deltas=c();weights=c();reps=c();M12=c();M21=c();TT=c();S2Mpara=c()
p1=1;p2=3
for(p1 in 1:(length(pops)-1)){
	for (p2 in (p1+1):length(pops)) {
		pop1=pops[p1];pop2=pops[p2]
		para=read.table(paste(pop1,pop2,"_s2m_params.txt",sep=""),sep=" ")
		para=para[,-c(9:10)]
		for(i in 4:8) { para[,i]=as.numeric(as.character(para[,i])) }
		check=apply(para[,4:8],1,mean)
		para=para[!is.na(check) & check!=Inf,]
		para=unique(para)
		goods=row.names(para)=para$V1
		para=para[,-1]
		th=read.table(paste(pop1,pop2,"_s2m_theta.txt",sep=""),sep=" ")
		row.names(th)=th$V1
		th=th[intersect(goods,th$V1),]
		li=read.table(paste(pop1,pop2,"_s2m_logLike.txt",sep=""),sep=" ")
		li[,4]=as.numeric(as.character(li[,4]))
		li=li[!(is.na(li[,4])),]
		row.names(li)=li$V1
		li=li[intersect(goods,li$V1),]
		liCut=quantile(li[,4],prob=likesCut)
		selects=row.names(li[li[,4]>=liCut[1],])
		sel2=row.names(th)[th$V4>1]
		selects=intersect(selects,sel2)
		para=para[selects,]
		names(para)=c("pop1","pop2","nu1","nu2","T","m12","m21")
		th=th[selects,]
		Ne=round(th[,4]/(4*mutRate),0)
		para$th=round(th[,4],0)
		para$nu1=round(para$nu1*Ne,0)
		para$nu2=round(para$nu2*Ne,0)
		para$logLike=li[selects,4]
		para$Nref=Ne
		para$M12=para$m12
		para$M21=para$m21
		para$m12=signif(para$m12/(2*Ne),2)
		para$m21=signif(para$m21/(2*Ne),2)
		para$TT=round(gentime*para$T*2*Ne,0)
		M12=rbind(M12,quantile(para$m12,prob=c(0.5,0.025,0.975)))
		M21=rbind(M21,quantile(para$m21,prob=c(0.5,0.025,0.975)))
		TT=rbind(TT,round(quantile(para$TT,prob=c(0.5,0.025,0.975)),0))
		popp1=append(popp1,pop1)
		popp2=append(popp2,pop2)
		S2Mpara=data.frame(rbind(S2Mpara,para))
	}
}

S2Mpara$nu1.0=S2Mpara$nu1/S2Mpara$Nref
S2Mpara$nu2.0=S2Mpara$nu2/S2Mpara$Nref
save(S2Mpara,file="s2m_allParameters_noBS.RData")
nus=data.frame(cbind(
	"ne"=c(S2Mpara$nu1,S2Mpara$nu2),
	"pop"=c(as.character(S2Mpara$pop1),as.character(S2Mpara$pop2))
))
nus$pop=factor(nus$pop,levels=c("W","S","O","M","K"))
nus$ne=as.numeric(as.character(nus$ne))
boxplot(log(ne,10)~pop,nus,	xlab="population",ylab="effective population size",
	yaxt="n"
)
axis(2,labels=c(1000,5000,25000,100000),at=log(c(1000,5000,25000,100000),10))

nus0=data.frame(cbind(
	"ne"=c(S2Mpara$nu1.0,S2Mpara$nu2.0),
	"pop"=c(as.character(S2Mpara$pop1),as.character(S2Mpara$pop2))
))
nus0$pop=factor(nus0$pop,levels=c("W","S","O","M","K"))
nus0$ne=as.numeric(as.character(nus0$ne))
summary(lm(ne~0+pop,nus0))
# popW  1.71168    0.03926  43.599   <2e-16 ***
# popS  1.81823    0.04855  37.452   <2e-16 ***
# popO  1.54246    0.04364  35.347   <2e-16 ***
# popM  2.01317    0.04709  42.750   <2e-16 ***
# popK  0.33642    0.03938   8.544   <2e-16 ***

# projections
# W 32  7108
# S 34  7197
# O 52  8172
# M 38  7567
# K 36  7000

res=data.frame(cbind(M12,M21,TT))
res$pop1=popp1
res$pop2=popp2
names(res)=c("m12.median","m12.lo","m12.hi","m21.median","m21.lo","m21.hi","T.median","T.lo","T.hi","pop1","pop2")

plot(density(res$T.median),yaxt="n",bty="n",ylab="",xlab="years before present",main="s2m: pop split time")

ibd0=read.csv("rad_Fst_dist_adig.csv")
ibd0=ibd0[-c(2,6,10,11,12),]
res$dist=ibd0$dist
row.names(res)=paste(res$pop1,res$pop2,sep=".")
save(res,file="migrationBoot_dadi7_S2M_noBS_fdr05.RData")

#----------------------------
ll=load("migrationBoot_dadi7_S2M_noBS_fdr05.RData")

# migration vs distance plot
pops=c("W","S","O","M","K")
res
plot(log(m21.median,10)~dist,res,xlim=c(10,1500),ylim=c(-4,-1.7),xlab="distance, km",mgp=c(2.3,1,0),ylab="immigration rate",yaxt="n")
axis(2,at=c(-4,-3,-2),labels=c(1e-4,1e-3,1e-2))
points(log(m12.median,10)~dist,res,pch=2,log="y")
legend("topright",pch=c(1,2),legend=c("S-ward","N-ward"),cex=0.9)

# producing migration matrix
migrat=matrix(0,nrow=length(pops),ncol=length(pops))
colnames(migrat)=pops
row.names(migrat)=pops
for(p1 in 1:(length(pops)-1)){
	for (p2 in (p1+1):length(pops)) {
 		pop1=pops[p1];pop2=pops[p2]
		migrat[pop1,pop2]=res[res$pop1==pop1 & res$pop2==pop2,"m12.median"]
		migrat[pop2,pop1]=res[res$pop1==pop1 & res$pop2==pop2,"m21.median"]
	}
}
#colnames(migrat)=c("W","S","O","M","K")
#row.names(migrat)=c("W","S","O","M","K")
write.csv(migrat,file="migrationMatrix_dadi7_s2m_noBS_fdr05b.csv")

library(pheatmap)
colf=colorRampPalette(c("white","lightyellow","orange","red","firebrick"))
cols=colf(100)
pheatmap(log(migrat+2e-4,10),cluster_rows=F,cluster_cols=F, color=cols,main="immigration rate",border_col="white")

library(ggplot2)

# migration rates barplot with whiskers
 
 migration=c(res[,4],res[,1])
 hi=c(res[,5],res[,2])
 lo=c(res[,6],res[,3])
 dist=rep(res$dist,2)
 r2=data.frame(cbind(migration,hi,lo,dist))
 r2$populations=rep(row.names(res),2)
 r2$direction=c(rep("S",nrow(res)),rep("N",nrow(res)))
 r2=r2[order(r2$dist),]
 r2$populations =factor(r2$populations,levels=unique(r2$populations))
 r2$direction=factor(r2$direction,levels=c("S","N"))
 #r2=r2[-grep("C",r2$populations),]
 pd = position_dodge(0.7)
 ggplot(r2,aes(populations,migration,fill=direction))+ylab("immigration rate")+geom_bar(stat="identity",position=pd)+ geom_errorbar(aes(ymin = lo, ymax = hi),position=pd,lwd = 0.3, width = 0.5, colour = "grey50")+theme_bw()+scale_fill_manual(values=c("grey20","grey80"))

