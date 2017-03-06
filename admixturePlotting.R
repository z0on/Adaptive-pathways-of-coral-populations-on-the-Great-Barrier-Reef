#install.packages("RColorBrewer")

# execute the chunk between ------- lines to load the function
#-----------------
plotAdmixture=function(data,npops,colors=NA,ord=NULL,angle=0,space=0,...) {
#data=z;npops=4;colors=cols;angle=0;ord=NULL
	tbl=data
	require(RColorBrewer)
	if (is.na(colors[1])){ colors=brewer.pal(n=npops,name="Set3") }
	p=levels(data$pop)[1]
	tsort=c()
	if (is.null(ord)) {
		for(p in levels(tbl$pop)){
			s=subset(tbl,pop==p)
			dd=cor(t(s[,1:npops]))
			dd[dd<0]=0
			dd=1-dd
			if (length(dd)>2){
				cl=hclust(as.dist(dd),method="single")
				sorted=cl$labels[cl$order]
				s=s[sorted,]
			}
			tsort=data.frame(rbind(tsort,s))
		}
	} else {
		tsort=tbl[ord,]
	}
	midpts=barplot(t(as.matrix(tsort[1:npops])), col=colors,xlab="population", space=space,ylab="Ancestry", border=NA, xaxt="n",mgp=c(2.3,1,0))
	pops=levels(tbl$pop)
	np=0;lp=0;lpoints=c()
	abline(v=0)
	abline(h=1,lwd=1)
	abline(h=0,lwd=1)
	for( p in pops) {
		np0=table(tbl$pop==p)[2]
		lp0=np0/2
		lp=np+lp0
		np=np+np0
		lpoints=append(lpoints,lp)
		abline(v=np)
	}
	text(as.numeric(lpoints)+2.5,par("usr")[3]-0.05,labels=pops,xpd=TRUE, srt=angle, pos=2)
	row.names(tsort)=tsort$ind
#	return(tsort)
	return(row.names(tsort))
}

#---------------------
# assembling the input table
dir="my/dir/with/Qfiles_and_inds2pops_table/"
inName="bestA.3.Q" # name of the input file to plot, output of ADMIXTURE run
npops=as.numeric(strsplit(inName,".",fixed=T)[[1]][length(strsplit(inName,".",fixed=T)[[1]])-1])
inds2pops="inds2pops" # two-column tab-delimited file listing sample names and their population assignemnts, in the order corresponding to the ADMIXTURE output.

tbl=read.table(paste(dir,inName,sep=""),header=F)
i2p=read.table(paste(dir,inds2pops,sep=""),header=F)
names(i2p)=c("ind","pop")
tbl=cbind(tbl,i2p)
tbl$pop=factor(tbl$pop,levels=c("W","S","O","M","K"))
row.names(tbl)=tbl$ind

head(tbl,20) # this is how the resulting dataset must look

# bright palette:
cols=c("gold", "purple1", "darkorange", "green3", "royalblue", "hotpink", "red3")  # or use your own colors
cols=c("green3", "purple1","gold", "darkorange",  "royalblue", "hotpink", "red3")  # or use your own colors
cols=c( "gold", "darkorange","green3",  "purple1", "royalblue", "hotpink", "red3")  # or use your own colors

# subdued palette (except gold):
cols=c("gold","skyblue", "palegreen2","plum2","grey20","slateblue","coral","sandybrown","grey50")  # or use your own colors
cols=c("skyblue", "plum","gold","palegreen2","khaki2","slateblue","coral","sandybrown")  # or use your own colors
cols=c("gold","palegreen2","skyblue","plum2", "khaki2","slateblue","coral","sandybrown")
cols=c("sandybrown","khaki1","pink")
quartz()
ords=plotAdmixture(data=tbl,npops=npops,colors=cols,angle=90) # ords now contain the names of individuals in order in which they are plotted, useful to plot something else about them and line up with admixture plot

