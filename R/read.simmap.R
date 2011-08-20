# Function reads SIMMAPv1.0 tree into memory as modified {ape} "phylo" or "multiPhylo" object.
# Tree format is:
# text<-"((A:{0,0.3:1,0.4:0,0.3},B:{0,1.0}):{0,1.0},(C:{0,0.25:1,0.75},D:{0,1.0}):{0,1.0});"
# Can also read modified "nexus" style files (produced by SIMMAP).
# States need not be binary and can be one or several characters long (e.g., "aqua" vs. "terr"):
# text<-"((A:{aqua,0.3:terr,0.4:aqua,0.3},B:{aqua,1.0}):{aqua,1.0},(C:{aqua,0.25:terr,0.75},D:{aqua,1.0}):{aqua,1.0});".
# Times spent on each edge in each state are stored in the list of named vectors, phy$maps.
# Total time spent in each state on each edge is stored in matrix phy$mapped.edge.
# Written by Liam J. Revell 2011

read.simmap<-function(file="",text,format="phylip",rev.order=TRUE){

	# check to see if reading from file and "nexus"
	if(file!=""&&format=="nexus"){
		# read modified nexus block
		# this code is adapted from read.nexus() in the ape package (Paradis et al. 2004)
		X<-scan(file=file,what="",sep="\n",quiet=TRUE)
		left<-grep("\\[",X)
		right<-grep("\\]",X)
		if(length(left)){
			w<-left==right
			if(any(w)){
				s<-left[w]
				X[s]<-gsub("\\[[^]]*\\]","",X[s])
			}
			w<-!w
			if(any(w)){
			s<-left[w]
			X[s]<-gsub("\\[.*","",X[s])
			sb<-right[w]
			X[sb]<-gsub(".*\\]","",X[sb])
			if(any(s<sb-1))
				X<-X[-unlist(mapply(":",(s+1),(sb-1)))]
			}
		}
		endblock<-grep("END;|ENDBLOCK;",X,ignore.case=TRUE)
		semico<-grep(";",X)
		i1<-grep("begin smptrees;",X,ignore.case=TRUE)
		i2<-grep("translate",X,ignore.case=TRUE)
		translation<-if(length(i2)==1 && i2>i1) TRUE else FALSE
		if(translation){
			end<-semico[semico>i2][1]
			x<-X[(i2+1):end]
			x<-unlist(strsplit(x,"[,;\t]"))
			x<-unlist(strsplit(x," ")) # this is an addition
			x<-x[nzchar(x)]
			trans<-matrix(x,ncol=2,byrow=TRUE)
			trans[,2]<-gsub("['\"]","",trans[,2])
			n<-dim(trans)[1]
		}
		start<-if(translation) semico[semico>i2][1]+1 else semico[semico>i1][1]
		end<-endblock[endblock>i1][1]-1
		tree<-X[start:end]
		rm(X)
		tree<-gsub("^.*= *","",tree)
		tree<-tree[tree!=""]
		semico<-grep(";",tree)
		Ntree<-length(semico)
		if(Ntree==1&&length(tree)>1){
			STRING<-paste(tree,collapse="")
		} else {
			if(any(diff(semico)!=1)){
				STRING<-character(Ntree)
				s<-c(1,semico[-Ntree]+1)
				j<-mapply(":",s,semico)
				if(is.list(j)){
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[[i]]],collapse="")
				} else {
					for(i in 1:Ntree) STRING[i]<-paste(tree[j[,i]],collapse="")
				}
			} else STRING<-tree
		}
		rm(tree); text<-STRING
		if(translation==TRUE){
			rownames(trans)<-trans[,1]; trans<-trans[,2]
		}
	}
	# end adapted "ape" code here

	# check to see if reading from file and "phylip"
	if(file!=""&&format=="phylip"){
		text<-scan(file,sep="\n",what="character")
		Ntree<-length(text)
	}

	# if reading from text
	if(file=="") Ntree<-length(text)

	# create "multiPhylo" object
	phy<-list(); class(phy)<-"multiPhylo"

	# loop over Ntree
	for(h in 1:Ntree){

		tree<-unlist(strsplit(text[h], NULL)) # split string 'text'
		Ntip<-1; Nnode<-0; i<-1
		while(tree[i]!=";"){
			while(tree[i]=="{") while(tree[i]!="}") i<-i+1 # skip maps
			if(tree[i]==",") Ntip<-Ntip+1 # count tips
			if(tree[i]=="(") Nnode<-Nnode+1 # count nodes
			i<-i+1
		}
		tip.label<-vector(mode="character") # create tip label vector
		edge<-matrix(data=0,Nnode+Ntip-1,2) # create edge matrix
		edge.length<-rep(0,Nnode+Ntip-1) # create edge.length vector
		maps<-list()
	
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
					# is there a branch length?
					if(tree[i]==":"){
						i<-i+1; l<-1
						if(tree[i]=="{"){
							maps[[j]]<-vector()
							i<-i+1
							while(tree[i]!="}"){
								temp<-vector(); m<-1
								while(tree[i]!=","){
									temp[m]<-tree[i]
									i<-i+1; m<-m+1
								}
								state<-paste(temp,collapse="")
								i<-i+1; temp<-vector(); m<-1
								while(tree[i]!=":"&&tree[i]!="}"){
									temp[m]<-tree[i]
									i<-i+1; m<-m+1
								}
								length<-as.numeric(paste(temp,collapse=""))
								maps[[j]][l]<-length; names(maps[[j]])[l]<-as.character(state)
								l<-l+1
								if(tree[i]==":") i<-i+1
							}
						}
						edge.length[j]<-sum(maps[[j]]) # create branch length
						i<-i+1
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
					i<-i+1; l<-1
					if(tree[i]=="{"){
						maps[[match(currnode,edge[,2])]]<-vector()
						i<-i+1
						while(tree[i]!="}"){
							temp<-vector(); m<-1
							while(tree[i]!=","){
								temp[m]<-tree[i]
								i<-i+1; m<-m+1
							}
							state<-paste(temp,collapse="")
							i<-i+1; temp<-vector(); m<-1
							while(tree[i]!=":"&&tree[i]!="}"){
								temp[m]<-tree[i]
								i<-i+1; m<-m+1
							}
							length<-as.numeric(paste(temp,collapse=""))
							maps[[match(currnode,edge[,2])]][l]<-length; names(maps[[match(currnode,edge[,2])]])[l]<-as.character(state)
							l<-l+1
							if(tree[i]==":") i<-i+1
						}
					}
					edge.length[match(currnode,edge[,2])]<-sum(maps[[match(currnode,edge[,2])]]) # create branch length
					i<-i+1
				}	
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
					# is there a branch length?
					if(tree[i]==":"){
						i<-i+1; l<-1
						if(tree[i]=="{"){
							maps[[j]]<-vector()
							i<-i+1
							while(tree[i]!="}"){
								temp<-vector(); m<-1
								while(tree[i]!=","){
									temp[m]<-tree[i]
									i<-i+1; m<-m+1
								}
								state<-paste(temp,collapse="")
								i<-i+1; temp<-vector(); m<-1
								while(tree[i]!=":"&&tree[i]!="}"){
									temp[m]<-tree[i]
									i<-i+1; m<-m+1
								}
								length<-as.numeric(paste(temp,collapse=""))
								maps[[j]][l]<-length; names(maps[[j]])[l]<-as.character(state)
								l<-l+1
								if(tree[i]==":") i<-i+1
							}
						}
						edge.length[j]<-sum(maps[[j]]) # create branch length
						i<-i+1
					}
				} else if(tree[i]=="("){
					nodecount<-nodecount+1 # creating a new internal node
					currnode<-nodecount; edge[j,2]<-currnode; # move to internal node
				}
				j<-j+1
			}
		}

		# construct the matrix mapped.edge (for backward compatibility)
		allstates<-vector()
		for(i in 1:nrow(edge)) allstates<-c(allstates,names(maps[[i]]))
		allstates<-unique(allstates)
		mapped.edge<-matrix(data=0,length(edge.length),length(allstates),dimnames=list(edge=apply(edge,1,function(x) paste(x,collapse=",")),state=allstates))
		for(i in 1:length(maps)) for(j in 1:length(maps[[i]])) mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]

		# assemble into modified "phylo" object
		temp<-list(edge=edge,Nnode=as.integer(Nnode),tip.label=tip.label,edge.length=edge.length,maps=maps,mapped.edge=mapped.edge); class(temp)<-"phylo"
		phy[[h]]<-temp

		# untranslate
		if(format=="nexus"&&translation==TRUE)
			phy[[h]]$tip.label<-trans[phy[[h]]$tip.label]

		# if rev.order==TRUE
		if(rev.order){
			phy[[h]]$maps<-lapply(phy[[h]]$maps,function(x) x<-x[length(x):1])
			attr(phy[[h]],"map.order")<-"right-to-left"
		} else
			attr(phy[[h]],"map.order")<-"left-to-right"

	}

	# if only one tree
	if(length(phy)==1)
		phy<-phy[[1]]

	return(phy) # return modified "phylo" or "multiPhylo" object

}
