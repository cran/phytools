# function returns the MRCA for a set of taxa (in tips)
# written by Liam Revell 2012

findMRCA<-function(tree,tips=NULL){
	if(is.null(tips)) return(mrca(tree))
	else {
		H<-nodeHeights(tree)
		X<-matrix(NA,length(tips),length(tips),dimnames=list(tips,tips))
		for(i in 1:length(tips)) for(j in i:length(tips)) X[i,j]<-X[j,i]<-fastMRCA(tree,tips[i],tips[j])
		n<-length(tips)
		nodes<-height<-vector(); k<-1
		for(i in 1:(n-1)) for(j in (i+1):n){
			nodes[k]<-X[tips[i],tips[j]]
			height[k]<-H[match(nodes[k],tree$edge[,1]),1]
			k<-k+1
		}
		z<-match(min(height),height)
		return(nodes[z])
	}
}

# fast pairwise MRCA function
# written by Liam Revell 2012

fastMRCA<-function(tree,sp1,sp2){
	x<-match(sp1,tree$tip.label)
	y<-match(sp2,tree$tip.label)
	a<-Ancestors(tree,x)
	b<-Ancestors(tree,y)
	z<-a%in%b
	return(a[min(which(z))])
}
