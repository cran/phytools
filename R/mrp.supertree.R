# function for Matrix Representation Parsimony supertree estimation in R
# uses the parsimony ratchet, pratchet() from the "phangorn" package for parsimony tree estimation
# written by Liam J. Revell 2011

mrp.supertree<-function(phy,weights=NULL){
	
	# some minor error checking
	if(!class(phy)=="multiPhylo") stop("phy must be object of class 'multiPhylo.'")

	X<-list() # list of bipartitions
	characters<-0 # number of characters

	for(i in 1:length(phy)){
		temp<-prop.part(phy[[i]]) # find all bipartitions
		# create matrix representation of phy[[i]] in X[[i]]
		X[[i]]<-matrix(0,nrow=length(phy[[i]]$tip),ncol=length(temp)-1)
		for(j in 1:ncol(X[[i]])) X[[i]][c(temp[[j+1]]),j]<-1
		rownames(X[[i]])<-attr(temp,"labels") # label rows
		if(i==1) species<-phy[[i]]$tip.label
		else species<-union(species,phy[[i]]$tip.label) # accumulate labels
		characters<-characters+ncol(X[[i]]) # count characters
	}

	XX<-matrix(data="?",nrow=length(species),ncol=characters,dimnames=list(species))
	
	j<-1
	for(i in 1:length(X)){
		# copy each of X into supermatrix XX
		XX[rownames(X[[i]]),c(j:((j-1)+ncol(X[[i]])))]<-X[[i]][1:nrow(X[[i]]),1:ncol(X[[i]])]
		j<-j+ncol(X[[i]])
	}

	# compute contrast matrix for phangorn
	contrast<-matrix(data=c(1,0,0,1,1,1),3,2,dimnames=list(c("0","1","?"),c("0","1")),byrow=TRUE)

	XX<-phyDat(XX,type="USER",contrast=contrast) # convert XX to phyDat object

	supertree<-pratchet(XX,trace=0,all=TRUE) # infer supertree(s) via pratchet()

	if(class(supertree)=="phylo")
		message(paste("The MRP supertree, optimized via pratchet(),\nhas a parsimony score of ",attr(supertree,"pscore")," (minimum ",characters,")",sep=""))
	else if(class(supertree)=="multiPhylo")
		message(paste("pratchet() found ",length(supertree)," supertrees\nwith a parsimony score of ",attr(supertree[[1]],"pscore")," (minimum ",characters,")",sep=""))

	return(supertree)

}
