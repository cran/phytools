# function conducts phylogenetic ANOVA & posthoc tests
# some code from phy.anova() in "geiger"
# written by Liam Revell 2011, 2015, 2017, 2021

phylANOVA<-function(tree,x,y,nsim=1000,posthoc=TRUE,p.adj="holm"){
	if(is.null(names(x))){
		cat("Warning: no labels for x. Assuming order of tree$tip.label.\n\n")
		names(x)<-tree$tip.label
	}
	x<-x[tree$tip.label]
	if(is.null(names(y))){
		cat("Warning: no labels for y. Assuming order of tree$tip.label.\n\n")
		names(y)<-tree$tip.label
	}	
	y<-y[tree$tip.label]
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	sig2<-mean(pic(y,multi2di(tree,random=FALSE))^2) # compute BM rate for y
	x<-as.factor(x) # change x to factor
	m<-length(levels(x))
	aov<-anova(lm(y~x))
	F.obs<-aov[1,4] # F on empirical data
	if(posthoc) T.obs<-tTests(x,y) # empirical ts
	sims<-fastBM(tree,sig2=sig2,nsim=(nsim-1)) # simulate
	F.null<-vector()
	if(posthoc) T.null<-array(NA,dim=c(m,m,nsim),dimnames=list(levels(x),levels(x),NULL))
	F.null[1]<-F.obs
	if(posthoc) T.null[,,1]<-T.obs
	for(i in 2:nsim){
		F.null[i]<-anova(lm(sims[,i-1]~x))[1,4]
		if(posthoc) T.null[,,i]<-tTests(x,sims[,i-1])
	}
	P.F<-sum(F.null>=F.obs)/nsim # p-value for F-test
	if(posthoc){
		P.T<-matrix(NA,m,m,dimnames=list(levels(x),levels(x)))
		# uncorrected p-values from simulation
		for(i in 1:m) for(j in i:m){
			P.T[i,j]<-sum(abs(T.null[i,j,])>=abs(T.obs[i,j]))/nsim
			P.T[j,i]<-P.T[i,j]
		}
		# control for multiple tests (if p.adj!="none")
		P.T[lower.tri(P.T)]<-p.adjust(P.T[lower.tri(P.T)],method=p.adj)
		for(i in 1:m) for(j in i:m) P.T[i,j]<-P.T[j,i]
		obj<-list(F=F.obs,Pf=P.F,T=T.obs,method=p.adj,Pt=P.T,
			"Sum Sq"=aov$"Sum Sq","Mean Sq"=aov$"Mean Sq")
	} else obj<-list(F=F.obs,Pf=P.F,"Sum Sq"=aov$"Sum Sq",
		"Mean Sq"=aov$"Mean Sq")
	class(obj)<-"phylANOVA"
	obj
}

# computes pairwise t-statistics with a pooled SD
# some code from pairwise.t.test() in "stats"
# written by Liam Revell 2011

tTests<-function(x,y){
	if(!is.factor(x)) x<-as.factor(x)
	ybar<-tapply(y,x,mean,na.rm=TRUE)
	s<-tapply(y,x,sd,na.rm=TRUE)
	n<-tapply(!is.na(y),x,sum)
	m<-length(levels(x))
	degf<-n-1
	total.degf<-sum(degf)
	pooled.sd<-sqrt(sum(s^2*degf)/total.degf)
	compare.levels<-function(i,j){
		dif<-ybar[i]-ybar[j]
		se.dif<-pooled.sd*sqrt(1/n[i]+1/n[j])
		t.val<-dif/se.dif
		return(t.val)
	}
	T<-matrix(NA,m,m,dimnames=list(levels(x),levels(x)))
	for(i in 1:m) for(j in 1:m) T[i,j]<-compare.levels(levels(x)[i],levels(x)[j])
	return(T)
}

## S3 print method

print.phylANOVA<-function(x,digits=6,...){
	cat("ANOVA table: Phylogenetic ANOVA\n\n")
	cat("Response: y\n")
	object<-data.frame(round(x$"Sum Sq",digits),
		round(x$"Mean Sq",digits),
		c(round(x$F,digits),""),
		c(round(x$Pf,digits),""))
	colnames(object)<-c("Sum Sq","Mean Sq","F value",
		"Pr(>F)")
	rownames(object)<-c("x","Residual")
	print(object)
	cat("\nP-value based on simulation.\n")
	cat("---------\n")
	cat("\n")
	if(!is.null(x$T)){
		cat(paste("Pairwise posthoc test using method = \"",
			x$method,"\"\n\n",sep=""))
		cat("Pairwise t-values:\n")
		print(round(x$T,digits))
		cat("\nPairwise corrected P-values:\n")
		print(round(x$Pt,digits))
		cat("---------\n")
		cat("\n")
	}
}
