pops=c("W","S","O","M","K") # population names
e.local=c(1.57,0,0,1.6,-1.8) # local environmental settting
e.sd=0.25 # SD of random environmentsl fluctuations (common to all populations)
e.increment=0.05 # increment in environmental setting per generation past burn-in period

sin.p=10 # period of sinusoidal climate fluctuations, in generations
sin.amp=0 # amplitude of sinusoidal climate fluctuations

burnin=2000 # length of burnin period
gmax=2500 # maximum number of generations

envs=c()
for (g in 1:gmax){
	if (g<=burnin){
		envs=cbind(envs,round(e.local+rnorm(1,0,e.sd)+sin(g*2*pi/sin.p)*0.5*sin.amp,3))
	}
	else {
		envs=cbind(envs,round(e.local+rnorm(1,0,e.sd)+sin(g*2*pi/sin.p)*0.5*sin.amp+(g-burnin)*e.increment,3));
	}
}
colnames(envs)=seq(gmax)
E=data.frame(t(envs))
names(E)=pops
te=stack(E)
names(te)=c("environment","population")
te$generation=seq(gmax)
library(ggplot2)

ggplot(te,aes(generation,environment,colour=population))+geom_line()

write.table(envs,row.names=F,quote=F,sep="\t",file="sim_environment.txt")
