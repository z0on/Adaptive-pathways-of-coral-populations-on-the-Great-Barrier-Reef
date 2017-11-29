#------------------------
slimplots=function(slimfile,burnin,start=-10,title=" ",max.generation=400,yscale=NULL,envfile="~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/reslim/sim_environment.txt",popnames=c("W","S","O","M","K")) {
#	slimfile="h1p0.5.a.out";burnin=2000;cf=1;start=-100;title=expression(paste(h^2,"=1 ; p = 0.5"));max.generation=300;yscale=c(0,1);popnames=c("W","S","O","M","K"); envfile="~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/sim_environment.txt"
	require(gridExtra)
	require(ggplot2)
	envi=data.frame(t(read.table(envfile,header=F)))
	names(envi)=c("generation",popnames)
	envi[,2:6]=envi[,2:6]+30
	slim=read.table(slimfile,sep="\t",skip=18)
	names(slim)=c("generation","population","fitness","phenotype")
	for (i in unique(slim$population)) {
		slim$population[slim$population==i]=popnames[i+1]
	}
	slim$population=factor(slim$population, levels=popnames)
	slim=slim[slim$generation>(burnin+start),]
	slim$generation=slim$generation-burnin
	if(is.null(max.generation)){ max.generation=max(slim$generation)}
	slim=slim[slim$generation<=max.generation,]
	slim$phenotype=slim$phenotype+30
	head(slim)
	denoise=c();p="K";eg=c()
	for (p in levels(slim$population)) {
		pf=subset(slim,population==p & fitness<0.011)
		extinct=min(pf$generation)
		eg=append(eg,extinct)
		pf=subset(slim,population==p)
		pf$phenotype[pf$generation>extinct]=NA
		# for (i in 2:nrow(pf)){
			# if (pf[i,"fitness"]/pf[i-1,"fitness"]>100 & pf[i-1,"fitness"]<0.001) {
				# pf[i,"fitness"]=pf[i-1,"fitness"]
			# }
		# }
		denoise=data.frame(rbind(denoise,pf))
	}
	mm=max(denoise$generation[!(is.na(denoise$phenotype))])
	slim=denoise[denoise$generation<=mm,]
	# for (p in levels(slim$population)) {
		# pf=subset(slim,population==p & generation>100)
		# for (i in 2:nrow(pf)){
			# if (pf[i,"fitness"]/pf[i-1,"fitness"]>100 & pf[i-1,"fitness"]<0.001) {
				# pf[i,"fitness"]=pf[i-1,"fitness"]
			# }
		# }
		# denoise=data.frame(rbind(denoise,pf))
	# }
#	slim=rbind(subset(slim,generation<101),denoise)
	st=stack(data.frame(envi[,c(2,4,6)]))
	names(st)=c("phenotype","population")
	st$generation=envi$generation-burnin
	st=st[st$generation %in% slim$generation,]
	st$population=factor(st$population,levels=c("W","O","K"))
		colors=c("firebrick","khaki3","green3","coral","skyblue")
	fit=ggplot(slim,aes(generation,fitness,colour=population))+
		geom_line(size=0.3)+
		scale_colour_manual(values=colors)+
		xlim(start,max.generation)+
		ylab("relative fitness")+
		ggtitle(title) 
	if (!is.null(yscale)) { fit=fit+scale_y_continuous(limits=yscale) }
	phen=ggplot()+
		geom_line(data=slim,aes(generation,phenotype,colour=population),size=0.3)+
		geom_line(data=st,aes(generation,phenotype,colour=population),size=0.3,alpha=0.5)+
		ylim(27,39)+
		ggtitle(" ")+
		xlim(start,max.generation)+
		ylab(expression("thermal tolerance,"^o*C))+
		scale_colour_manual(values=colors)
#		theme(legend.position="bottom")
	grid.arrange(fit,phen)
}
#------------------------
setwd("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/reslim") 

slimplots("e0a.out",burnin=2000,start=-100,title=expression(paste("e=0; s=1")),max.generation=250,yscale=c(0,1))
slimplots("e1a.out",burnin=2000,start=-100,title=expression(paste("e=1; s=1")),max.generation=250,yscale=c(0,1))
slimplots("e2a.out",burnin=2000,start=-100,title=expression(paste("e=2; s=1")),max.generation=250,yscale=c(0,1))

slimplots("e0s0.5a.out",burnin=2000,start=-100,title=expression(paste("e=0; s=0.5")),max.generation=250,yscale=c(0,1))
slimplots("e1s0.5a.out",burnin=2000,start=-100,title=expression(paste("e=1; s=0.5")),max.generation=250,yscale=c(0,1))
slimplots("e2s0.5a.out",burnin=2000,start=-100,title=expression(paste("e=2; s=0.5")),max.generation=250,yscale=c(0,1))

slimplots("e0s2a.out",burnin=2000,start=-100,title=expression(paste("e=0; s=2")),max.generation=250,yscale=c(0,1))
slimplots("e1s2a.out",burnin=2000,start=-100,title=expression(paste("e=1; s=2")),max.generation=250,yscale=c(0,1))
slimplots("e2s2a.out",burnin=2000,start=-100,title=expression(paste("e=2; s=2")),max.generation=250,yscale=c(0,1))

slimplots("e1s1_noS_a.out",burnin=2000,start=-100,title=expression(paste("e=1; s=1, no S-ward")),max.generation=250,yscale=c(0,1))
slimplots("e1s1_noMig_a.out",burnin=2000,start=-100,title=expression(paste("e=1; s=1, no migration")),max.generation=250,yscale=c(0,1))

# with mutations
setwd("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/reslim/EwithMut") 
slimplots("e1s1.a.out",burnin=2000,start=-100,title=expression(paste("e=0; s=1")),max.generation=250,yscale=c(0,0.75))
slimplots("e1s1_mut6.a.out",burnin=2000,start=-100,title=expression(paste("e=0; s=1; mu=1e-6")),max.generation=250,yscale=c(0,.75))
slimplots("e1s1_mut5.c.out",burnin=2000,start=-100,title=expression(paste("e=0; s=1; mu=1e-5")),max.generation=250,yscale=c(0,.75))
slimplots("e1s1_mut4.b.out",burnin=2000,start=-100,title=expression(paste("e=0; s=1; mu=1e-4")),max.generation=250,yscale=c(0,.75))



#-----------
# noise amplitude

setwd("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/reslim/")

T0=data.frame(t(read.table("sim_environment_noNoise.txt",sep="\t")))
names(T0)=c("generation","W","S","O","M","K")
T0=T0[T0$generation>1950 & T0$generation<2061,]
T0$generation=	T0$generation-2000
head(T0)

T1=data.frame(t(read.table("sim_environment_c.txt",sep="\t")))
names(T1)=c("generation","W","S","O","M","K")
T1=T1[T1$generation>1950 & T1$generation<2061,]
T1$generation=	T1$generation-2000
head(T1)

plot(-0.1*(T1$S-T0$S)~T1$generation,type="l",col="grey40",ylim=c(-0.18,0.15),mgp=c(2.1,1,0),xlab="generations since warming started",ylab="fitness deviation",yaxt="n")
tlabs=c(1,0.5,0,-0.5,-1)
flabs=c(-0.1,-0.05,0, 0.05,0.1)
axis(2,labels=flabs,at=flabs)
axis(4,labels=tlabs,at=flabs)


slimfile="e1b.out"
slim=read.table(slimfile,sep="\t",skip=18)
names(slim)=c("generation","population","fitness","phenotype")
slim= slim[slim $generation>1950 & slim $generation<2061,]
slim$generation=	slim$generation-2000
head(slim)

pop=0
sp=subset(slim,population==pop)
lo=loess(fitness~generation,sp,span=0.25)
sp$Wlo=predict(lo,newdata=sp)
lines((sp$fitness-sp$Wlo)~sp$generation,col="coral")

pop=4
sp=subset(slim,population==pop)
lo=loess(fitness~generation,sp,span=0.25)
sp$Klo=predict(lo,newdata=sp)
lines((sp$fitness-sp$Klo)~sp$generation,col="cyan2")

abline(v=0,lty=3,col="grey50")
abline(h=0,lty=3,col="grey50")

#---------------






