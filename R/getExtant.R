# get all the extant/extinct tip names
# written by Liam J. Revell 2012

getExtant<-function(tree,tol=1e-8){
	H<-nodeHeights(tree)
	tl<-max(H)
	x<-which(H[,2]>=(tl-tol))
	y<-tree$edge[x,2]
	y<-y[y<=length(tree$tip)]
	z<-tree$tip.label[y]
	return(z)
}

getExtinct<-function(tree,tol=1e-8) setdiff(tree$tip.label,getExtant(tree,tol))
