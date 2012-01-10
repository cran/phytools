# function
# written by Liam J. Revell 2012

sampleFrom<-function(xbar=0,xvar=1,n=1,randn=NULL,type="norm"){
	if(length(xvar)==1&&length(xbar)!=length(xvar)) xvar<-rep(xvar,length(xbar))
	if(!is.null(randn))
		for(i in 1:length(xbar)) n[i]<-floor(runif(n=1,min=randn[1],max=(randn[2]+1)))
	x<-vector()
	for(i in 1:length(xbar)){
		y<-rnorm(n=n[i],mean=xbar[i],sd=sqrt(xvar[i]))
   		names(y)<-rep(names(xbar)[i],length(y))
   		x<-c(x,y)
	}
	return(x)
}
		