# function plots xkcd style (i.e., cartoonish with a hand drawn look) trees
# written by Liam J. Revell 2012, 2013

xkcdTree<-function(tree,file=NULL,gsPath=NULL,fsize=2,lwd=7,color="blue",dim=c(8.5,11),jitter=0.001,waver=c(0.1,0.1),tilt=0,right=TRUE){
	margin<-0.1
	message(paste(c("**** NOTE: use in Windows requires the installation of Ghostscript to\n", 
					"**** embed 'xkcd' font in PDF plot")))
	# if gsPath not specified
	if(is.null(gsPath)) gsPath="C:/Program Files/gs/gs9.10/bin/gswin32.exe"
	if(.Platform$OS.type=="windows") Sys.setenv(R_GSCMD=gsPath)
	# try to load extrafont
	if(is.na(match("xkcd", extrafont::fonts()))){
		message(paste(c("**** cannot find 'xkcd' font, searching for 'xkcd.ttf'\n",
						"**** (this may take a moment)\n",
						"**** if this fails, it is probably because 'xkcd' font is not installed\n"))) 
		 extrafont::font_import(prompt=FALSE,pattern="xkcd.ttf")
	}
	if(is.na(match("xkcd", extrafont::fonts()))) stop(paste(c("\n",
								"**** cannot find 'xkcd' font\n",
								"**** locate & install 'xkcd.ttf'")))
	 extrafont::loadfonts(quiet=TRUE)
	# if filename not specified	
	if(is.null(file)) file<-"xkcdTree.pdf"
	# check filename ending
	temp<-strsplit(file,"")[[1]]
	if(paste(temp[(length(temp)-3):length(temp)],collapse="")!=".pdf") file<-paste(file,".pdf",collapse="",sep="")
	# check tree
	if(class(tree)!="phylo") stop("tree should be object of class 'phylo.'")
	if(is.null(tree$edge.length)) tree<-compute.brlen(tree)
	# swap out "_" character for spaces (assumes _ is a place holder)
	tree$tip.label<-gsub("_"," ",tree$tip.label)
	tree$tip.label<-paste("",tree$tip.label)
	# reorder
	cw<-reorder(tree)
	pw<-reorder(tree,"pruningwise")
	# get coordinates for nodes & edges
	n<-length(cw$tip); m<-cw$Nnode
	Y<-matrix(NA,m+n,1)
	Y[cw$edge[cw$edge[,2]<=length(cw$tip),2]]<-1:n-1
	nodes<-unique(pw$edge[,1])
	for(i in 1:m){
		desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
		Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
	}
	Y<-Y/max(Y)
	H<-nodeHeights(tree); if(!is.null(tree$root.edge)) H<-H+tree$root.edge
	H<-H/max(H)
	# open plot
	pdf(file,family="xkcd",width=dim[1],height=dim[2])
	par(mar=rep(margin,4))
	plot.new()
	if(max(strwidth(cw$tip.label,units="figure",cex=fsize,family="xkcd"))<1.0){
		dx<-abs(max(strwidth(cw$tip.label,units="figure",cex=fsize,family="xkcd"))*cos(deg2rad(tilt)))
		dy<-abs(max(strwidth(cw$tip.label,units="figure",cex=fsize,family="xkcd"))*sin(deg2rad(tilt)))
		if(!right){ 
			dx<-dx*dim[1]/dim[2]
			dy<-dy*dim[2]/dim[1]
		}
		H<-H*(1-dx)
		Y<-Y*(1-dy)
		if((tilt>0&&!right)||(tilt<0&&right)) Y<-Y+dy
	} else message("Font size too large to properly rescale tree to window.")
	plot.window(xlim=c(0,1),ylim=c(0,1))
	for(i in 1:m){
		x<-H[which(cw$edge[,1]==nodes[i]),1]; y<-Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]]
		if(nodes[i]==(n+1)) rootPosn<-mean(y)
		k<-max(ceiling(abs(y[2]-y[1])/((max(Y)-min(Y))*waver[2]))-1,0)
		x<-c(x[1],rep(x[1],k),x[2])
		y<-y[1]+0:(k+1)/(k+1)*(y[2]-y[1])
 		if(right) xkcdLine(x,y,color=color,lwd=lwd,scale=1/jitter)
		else xkcdLine(y,x,color=color,lwd=lwd,scale=1/jitter)
	}
	if(!is.null(tree$root.edge)){
		x<-c(0,min(H)); y<-rep(rootPosn,2)
		k<-max(ceiling(abs(x[2]-x[1])/((max(H)-min(H))*waver[1]))-1,0)
		x<-x[1]+0:(k+1)/(k+1)*(x[2]-x[1])
		y<-c(y[1],rep(y[1],k),y[2])
		if(right) xkcdLine(x,y,color=color,lwd=lwd,scale=1/jitter)
		else xkcdLine(y,x,color=color,lwd=lwd,scale=1/jitter)
	}
	for(i in 1:nrow(cw$edge)){
		x<-H[i,]; y<-c(Y[cw$edge[i,2]],Y[cw$edge[i,2]])
		k<-max(ceiling(abs(x[2]-x[1])/((max(H)-min(H))*waver[1]))-1,0)
		x<-x[1]+0:(k+1)/(k+1)*(x[2]-x[1])
		y<-c(y[1],rep(y[1],k),y[2])
		if(right) xkcdLine(x,y,color=color,lwd=lwd,scale=1/jitter)
		else xkcdLine(y,x,color=color,lwd=lwd,scale=1/jitter)
	}
	for(i in 1:n){
		if(right) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],adj=c(0,0.5),cex=fsize,srt=tilt)
		else text(Y[i],H[which(cw$edge[,2]==i),2],cw$tip.label[i],adj=c(0,0.5),cex=fsize,srt=tilt+90)
	}
	dev.off()
	message(paste(c("**** NOTE: an 'embed_fonts' error most likely means that Ghostscript\n",
					"**** is not installed or that the path is incorrect")))
	extrafont::embedFonts(file)
	extrafont::embed_fonts(file)
	# reset margin
	par(mar=c(5,4,4,2)+0.1)
}

# function plots xkcd style line
# modified from function by user295691 on http://stackoverflow.com/questions/12675147

xkcdLine<-function(x,y,color,lwd,scale=1000){
	len<-length(x);
	rg<-par("usr");
	yjitter<-(rg[4]-rg[3])/scale;
	xjitter<-(rg[2]-rg[1])/scale;
	x_mod<-x+rnorm(len)*xjitter;
	y_mod<-y+rnorm(len)*yjitter;
	lines(x_mod,y_mod,col=color,lwd=lwd);
}

# function converts deg to rad
# written by Liam J. Revell 2012

deg2rad<-function(x) x*pi/180

