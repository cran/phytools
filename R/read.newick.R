# function reads Newick style tree with branch lengths into memory as an {ape} "phylo" object
# written by Liam J. Revell 2011

read.newick<-function(file="",text){

	# check to see if reading from file
	if(file!=""){
		text<-scan(file,sep="\n",what="character") # read from file
	}

	tree<-unlist(strsplit(text, NULL)) # split string 'text'
	Ntip<-1; Nnode<-0
	for(i in 1:length(tree)){ 
		if(tree[i]==",") Ntip<-Ntip+1 # count tips
		if(tree[i]=="(") Nnode<-Nnode+1 # count nodes
	}
	tip.label<-vector(mode="character") # create tip label vector
	edge<-matrix(data=0,Nnode+Ntip-1,2) # create edge matrix
	edge.length<-rep(0,Nnode+Ntip-1) # create edge.length vector
	ntip<-vector(mode="numeric")
	
	currnode<-Ntip+1; nodecount<-currnode # start with root node
	i<-1; j<-1; k<-1
	while(tree[i]!=";"){
		if(tree[i]=="("){
			edge[j,1]<-currnode; i<-i+1
			# is the next element a label?
			if(is.na(match(tree[i],c("(",")",",",":",";")))){
				l<-1; temp<-vector()
				while(is.na(match(tree[i],c(",",":",")")))){
					temp[l]<-tree[i]; l<-l+1; i<-i+1
				}
				tip.label[k]<-paste(temp,collapse=""); edge[j,2]<-k; k<-k+1 # create tip label
				ntip[j]<-1 # count of the number of descendants (one in this case, for a tip)
				# is there a branch length?
				if(tree[i]==":"){
					i<-i+1; l<-1; temp<-vector()
					while(is.na(match(tree[i],c(",",")")))){
						temp[l]<-tree[i]; l<-l+1; i<-i+1
					}
					branch<-as.numeric(paste(temp,collapse="")); edge.length[j]<-branch # create branch length
				}	
			} else if(tree[i]=="("){
				nodecount<-nodecount+1 # creating a new internal node
				currnode<-nodecount; edge[j,2]<-currnode; # move to new internal node
			}
			j<-j+1;
		} else if(tree[i]==")"){
			i<-i+1
			# is there a branch length?
			if(tree[i]==":"){
				i<-i+1; l<-1; temp<-vector()
				while(is.na(match(tree[i],c(",",")")))){
					temp[l]<-tree[i]; l<-l+1; i<-i+1
				}
				branch<-as.numeric(paste(temp,collapse="")); edge.length[match(currnode,edge[,2])]<-branch # assign branch length
			}
			ntip[match(currnode,edge[,2])]<-sum(ntip[which(edge[,1]==currnode)])
			currnode<-edge[match(currnode,edge[,2]),1] # move down the tree
		} else if(tree[i]==","){
			edge[j,1]<-currnode; i<-i+1
			# is the next element a label?
			if(is.na(match(tree[i],c("(",")",",",":",";")))){
				l<-1; temp<-vector()
				while(is.na(match(tree[i],c(",",":",")")))){
					temp[l]<-tree[i]; l<-l+1; i<-i+1
				} 
				tip.label[k]<-paste(temp,collapse=""); edge[j,2]<-k; k<-k+1 # create tip label
				ntip[j]<-1 # count of the number of descendants (one in this case, for a tip)
				# is there a branch length?				
				if(tree[i]==":"){
					i<-i+1; l<-1; temp<-vector()
					while(is.na(match(tree[i],c(",",")")))){
						temp[l]<-tree[i]; l<-l+1; i<-i+1
					}
					branch<-as.numeric(paste(temp,collapse="")); edge.length[j]<-branch # create branch length
				}
			} else if(tree[i]=="("){
				nodecount<-nodecount+1 # creating a new internal node
				currnode<-nodecount; edge[j,2]<-currnode; # move to internal node
			}
			j<-j+1
		}
	}
	
	# assemble into "phylo" object
	if(sum(edge.length)>1e-8)		
		phy<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,Ndesc=ntip) # with branch lengths
	else
		phy<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,Ndesc=ntip) # without branch lengths
	
	class(phy)<-"phylo" # assign class

	return(phy) # return object

}
