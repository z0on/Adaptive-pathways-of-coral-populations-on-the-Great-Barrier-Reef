# This scipt contains a function to produce fitness and phenotype plots
# based on SLiM model output, and its example uses to produce Fig. 3 panels.

#------------------------
slimplots=function(slimfile,burnin,start=-10,title=" ",cf=1,max.generation=400,yscale=NULL) {
	require(gridExtra)
	require(ggplot2)
	slim=read.table(slimfile,sep="\t",skip=77)
#head(slim)
	names(slim)=c("generation","T","W","S","O","M","K","pW","pS","pO","pM","pK")
	slim$T=slim$T++29.93
	slim=slim[slim$generation>(burnin+start),]
	slim$generation=slim$generation-burnin
	if(is.null(max.generation)){ max.generation=max(slim$generation)}
	slim=slim[slim$generation<=max.generation,]
	fits=stack(slim[,3:7])
	names(fits)=c("fitness","population")
	fits$fitness=fits$fitness*cf
	fits$generation=slim$generation
	fits$population=factor(fits$population,levels=c("W","S","O","M","K"))
	denoise=c();p="K"
	for (p in levels(fits$population)) {
		pf=subset(fits,population==p & generation>100)
		for (i in 2:nrow(pf)){
			if (pf[i,"fitness"]/pf[i-1,"fitness"]>100 & pf[i-1,"fitness"]<0.001) {
				pf[i,"fitness"]=pf[i-1,"fitness"]
			}
		}
		denoise=data.frame(rbind(denoise,pf))
	}
	fits=rbind(subset(fits,generation<101),denoise)
		
	names(slim)=c("generation","T","pW","pS","pO","pM","pK","W","S","O","M","K")
	pheno=stack(slim[,8:12])
	names(pheno)=c("phenotype","population")
	pheno$generation=slim$generation
#	head(pheno)
	pheno$phenotype=pheno$phenotype+29.93
	pheno$population=factor(pheno$population,levels=c("W","S","O","M","K"))
	
	W=slim$T+1.6
	O=slim$T
	K=slim$T-1.8
	t=data.frame(cbind("generation"=slim$generation,W,O,K))
	st=stack(t[,2:4])
	names(st)=c("phenotype","population")
	st$generation=t$generation
	st$population=factor(st$population,levels=c("W","O","K"))
		colors=c("firebrick","khaki3","green3","coral","skyblue")
	fit=ggplot(fits,aes(generation,fitness,colour=population))+
		geom_line()+
		scale_colour_manual(values=colors)+
		ylab("relative fitness")+
		ggtitle(title) 
	if (!is.null(yscale)) { fit=fit+scale_y_continuous(limits=yscale) }
	phen=ggplot()+
		geom_line(data=pheno,aes(generation,phenotype,colour=population))+
		geom_line(data=st,aes(generation,phenotype,colour=population),size=0.5,alpha=0.25)+
		ylim(27,50)+
		ggtitle(" ")+
		ylab(expression("thermal tolerance,"^o*C))+
		scale_colour_manual(values=colors)

	grid.arrange(fit,phen)
}

#-------------
# plasticity =0.5, with noise

slimplots("h1p0.5.out2",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("h0.5p0.5.out1",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.5 ; p = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("h0.25p0.5.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.25 ; p = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("h1e-5p0.5.out1",burnin=2000,start=-100,title=expression(paste(h^2,"= 1e-5 ; p = 0.5")),max.generation=300,yscale=c(0,1))

#-------------
# plasticity =1, with noise

slimplots("h1p1.out2",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 1")),max.generation=300,yscale=c(0,1))

slimplots("h0.5p1.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.5 ; p = 1")),max.generation=300,yscale=c(0,1))

quartz()
slimplots("h0.25p1.out1",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.25 ; p = 1")),max.generation=300,yscale=c(0,1))

# no southward migration
slimplots("h0.25p1noS.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.25 ; p = 1, no S migr.")),max.generation=300,yscale=c(0,1))

# no migration at all
slimplots("h0.25p1noM.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.25 ; p = 2, no migr.")),max.generation=300,yscale=c(0,1))

#-------------------
# plasticity =2, with noise

slimplots("h1p2.out2",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 2")),max.generation=300,yscale=c(0,1))

slimplots("h0.5p2.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.5 ; p = 2")),max.generation=300,yscale=c(0,1))

slimplots("h0.25p2.out2",burnin=2000,start=-100,title=expression(paste(h^2,"= 0.25 ; p = 2")),max.generation=300,yscale=c(0,1))

slimplots("h1e-5p2.out1",burnin=2000,start=-100,title=expression(paste(h^2,"= 1e-5 ; p = 1")),max.generation=300,yscale=c(0,1))


#-----------
# noise amplitude
setwd("")

	slim0=read.table("h0.25p1s.out1",sep="\t",skip=77)
	names(slim0)=c("generation","T","W","S","O","M","K","pW","pS","pO","pM","pK")
	slim0=slim0[slim0$generation>1950 & slim0$generation<2151,]
	slim0$generation=	slim0$generation-2000
	
	slim=read.table("h0.25p1.out3",sep="\t",skip=77)
	names(slim)=c("generation","T","W","S","O","M","K","pW","pS","pO","pM","pK")
	slim=slim[slim$generation>1950 & slim$generation<2151,]
	slim$generation=	slim$generation-2000

plot(-0.1*(slim$T-slim0$T)~slim$generation,type="l",col="grey40",ylim=c(-0.18,0.15),mgp=c(2.1,1,0),xlab="generations since warming started",ylab="fitness deviation",yaxt="n")
tlabs=c(1,0.5,0,-0.5,-1)
flabs=c(-0.1,-0.05,0, 0.05,0.1)
axis(2,labels=flabs,at=flabs)
axis(4,labels=tlabs,at=flabs)

lo=loess(W~generation,slim,span=0.25)
slim$Wlo=predict(lo,newdata=slim)
#lines(slim$W~slim$generation,type="l")
#lines(slim$Wlo~slim$generation,type="l",col="red")
lines((slim$W-slim$Wlo)~slim$generation,col="coral")
lo=loess(K~generation,slim,span=0.25)
slim$Klo=predict(lo,newdata=slim)
# plot(slim$K~slim$generation,type="l")
# lines(slim$Klo~slim$generation,type="l",col="red")
lines((slim$K-slim$Klo)~slim$generation,type="l",col="cyan2")
abline(v=0,lty=3,col="grey50")
abline(h=0,lty=3,col="grey50")

lo=loess(O~generation,slim,span=0.25)
slim$Olo=predict(lo,newdata=slim)
# plot(slim$O~slim$generation,type="l")
# lines(slim$Olo~slim$generation,type="l",col="red")
lines((slim$O-slim$Olo)~slim0$generation,col="green2")

x=seq(-6,6,0.01)
max05=dnorm(0,0,0.5)
max1=dnorm(0,0,1)
max2=dnorm(0,0,2)
plot(dnorm(x,0,2)/max2~x,type="l")
lines(dnorm(x,0,1)/max1~x,type="l",col="red")

1-(max1-dnorm(1.17,0,1))/max1


