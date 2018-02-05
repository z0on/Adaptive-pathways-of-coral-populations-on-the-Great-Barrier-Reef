# run this line first:
# ls *slim*out | perl -pe 's/(\S+)/grep \"#Total\" $1 \| perl -pe \"s\/#TotalGenVar\\s\/\/\" >$1\.G/' | bash

setwd("~/Dropbox/Documents/milleporaRAD_paper_2015/milleporaRAD_v2/SLIM_simulation_GBR/equilibrium_v5.2/m1e-6/") 
infiles=c(
"e0_p0.5.slim.a2.out",
"e0_p1.slim.a2.out",
"e0_p2.slim.a2.out",
"e1_p0.5.slim.a2.out",
"e1_p1.slim.a2.out",
"e1_p2.slim.a2.out",
"e2_p0.5.slim.a2.out",
"e2_p1.slim.a2.out",
"e2_p2.slim.a2.out"
)
# infiles=c(
# "e0_p0.5.slim.b2.out",
# "e0_p1.slim.b2.out",
# "e0_p2.slim.b2.out",
# "e1_p0.5.slim.b2.out",
# "e1_p1.slim.b2.out",
# "e1_p2.slim.b2.out",
# "e2_p0.5.slim.b2.out",
# "e2_p1.slim.b2.out",
# "e2_p2.slim.b2.out"
# )
# infiles=c(
# "e0_p0.5.slim.c2.out",
# "e0_p1.slim.c2.out",
# "e0_p2.slim.c2.out",
# "e1_p0.5.slim.c2.out",
# "e1_p1.slim.c2.out",
# "e1_p2.slim.c2.out",
# "e2_p0.5.slim.c2.out",
# "e2_p1.slim.c2.out",
# "e2_p2.slim.c2.out"
# )


infiles=infiles[c(1,4,7,2,5,8,3,6,9)]
quartz()
par(mfrow=c(3,3),mai=c(0.4,0.4,0.4,0.4))

 # infiles=c(
 # "e1_p1.slim.a2.out",
 # "50k.slim.a2.out",
 # "c100.slim.a2.out",
 # "c100_50k.slim.a2.out"
 # )
 # par(mfrow=c(2,2),mai=c(0.4,0.4,0.4,0.4))

for (i in infiles) {
	slimfile=i
	slim=read.table(slimfile,sep="\t",skip=18)
	tail(slim)
	names(slim)=c("generation","population","fitness","phenotype","nmuts","G")
#	slim$G=2*slim$G
	slim= slim[slim $generation>19500 & slim $generation<20200,]
	slim$generation=	slim$generation-20000
	head(slim)
	g2=paste(slimfile,".G",sep="")
	g2s=read.table(g2,sep="\t",skip=18)
	names(g2s)=c("generation","nmuts","G")
#	g2s$G=2*g2s$G
	g2s = g2s[g2s$generation>19500 & g2s$generation<20200,]
	g2s$generation= g2s$generation-20000
	plot(G~generation,type="l",g2s,col="white",lwd=1,lty=2,ylim=c(0,0.4),main=sub(".slim.a2.out","",i),mgp=c(2.1,1,0),xlim=c(-500,200))
	colors=c("firebrick","khaki3","green3","coral","skyblue")
	if (i==infiles[1]) { legend("topleft",lty=c(1,1,1,1,1,2),lwd=c(1,1,1,1,1,2),col=c(colors,"grey80"),c("W","S","O","M","K","All"),bty="n",cex=0.8) }
	for (p in c(0:4)){ lines(G~generation,type="l",subset(slim,population==p),col=colors[p+1]) }
	lines(G~generation,type="l",g2s,col="grey80",lwd=2,lty=2)
	abline(v=0,lty=3,col="grey50")
}


