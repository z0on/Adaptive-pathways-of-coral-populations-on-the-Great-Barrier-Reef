#load('~/Documents/milleporaRAD_paper_2015/water_quality_WSOMK_2015/WSOMK_environment.RData')


# library(ggplot2)
# cols=c("cyan3","green2","orange","red","purple")
# ggplot(st,aes(generation,T,colour=pop))+geom_line()+theme_bw()+scale_colour_manual(values=cols)


#------------------------
slimplots=function(slimfile,burnin,start=-10,title=" ",max.generation=400,yscale=NULL,envfile="~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/sim_environment.txt",popnames=c("W","S","O","M","K")) {
#	slimfile="~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.25s1.out";burnin=2000;cf=1;start=-100;title=expression(paste(h^2,"=1 ; p = 0.5"));max.generation=300;yscale=c(0,1);popnames=c("W","S","O","M","K"); envfile="~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/sim_environment.txt"
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
	denoise=c();p="K"
	for (p in levels(slim$population)) {
		pf=subset(slim,population==p & generation>100)
		for (i in 2:nrow(pf)){
			if (pf[i,"fitness"]/pf[i-1,"fitness"]>100 & pf[i-1,"fitness"]<0.001) {
				pf[i,"fitness"]=pf[i-1,"fitness"]
			}
		}
		denoise=data.frame(rbind(denoise,pf))
	}
	slim=rbind(subset(slim,generation<101),denoise)
	st=stack(data.frame(envi[,c(2,4,6)]))
	names(st)=c("phenotype","population")
	st$generation=envi$generation-burnin
	st=st[st$generation %in% slim$generation,]
	st$population=factor(st$population,levels=c("W","O","K"))
		colors=c("firebrick","khaki3","green3","coral","skyblue")
	fit=ggplot(slim,aes(generation,fitness,colour=population))+
		geom_line()+
		scale_colour_manual(values=colors)+
		ylab("relative fitness")+
		ggtitle(title) 
	if (!is.null(yscale)) { fit=fit+scale_y_continuous(limits=yscale) }
	phen=ggplot()+
		geom_line(data=slim,aes(generation,phenotype,colour=population))+
		geom_line(data=st,aes(generation,phenotype,colour=population),size=0.5,alpha=0.25)+
		ylim(27,50)+
		ggtitle(" ")+
		ylab(expression("thermal tolerance,"^o*C))+
		scale_colour_manual(values=colors)
	grid.arrange(fit,phen)
}
#------------------------


#-------------
# plasticity =0.5, with noise

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h1s0.5.out",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; s = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.5s0.5.out",burnin=2000,start=-100,title=expression(paste(h^2,"=0.5 ; s = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.25s0.5.out",burnin=2000,start=-100,title=expression(paste(h^2,"=0.25 ; s = 0.5")),max.generation=300,yscale=c(0,1))


#-------------
# plasticity =1, with noise

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h1s1.out",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; s = 1")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.5s1.out",burnin=2000,start=-100,title=expression(paste(h^2,"=0.5 ; s = 1")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.25s1.out",burnin=2000,start=-100,title=expression(paste(h^2,"=0.25 ; s = 1")),max.generation=300,yscale=c(0,1))

#-------------------
# plasticity =2, with noise

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h1s2.out",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.5s2.out",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 0.5")),max.generation=300,yscale=c(0,1))

slimplots("~/Documents/SLiM_adaptation_metapop_variableSelection_Dec2016/SLiM_fileReading/slimQTL_v2_h0.25s2.out",burnin=2000,start=-100,title=expression(paste(h^2,"=1 ; p = 0.5")),max.generation=300,yscale=c(0,1))

#-------------------



