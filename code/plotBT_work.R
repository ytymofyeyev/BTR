
#### plot 2-treatment case
plotBT.2D <- function(btObj,nGen=NULL){
	if (is.null(nGen)) nGen<-length(btObj)
	if(length(btObj)<nGen) stop("Need to define more genereration of BT obeject to be plotted\n")
	w.norm <- attr(btObj,"w")/sum(attr(btObj,"w"))
	plot(rbind(rep(0,2),ceiling(w.norm*nGen)),type="n",xlab="TRT A",ylab="TRT B")
	for (g in 1:nGen){
		for (nodeI in 1:dim(btObj[[g]]$nodeMat)[2]){
			node <- btObj[[g]]$nodeMat[,nodeI]
			for ( i in 1:2){   # connect node to next generation nodes
				cnodeI <- btObj[[g]]$mapToNextGenerNode[2*(nodeI-1)+i] # index of connected node
				if (is.na(cnodeI)) next
				cnode <- btObj[[g]]$nextGenerNodeMat[,cnodeI]
				x <- c( node[1], cnode[1])
				y <- c( node[2], cnode[2])
				points(x,y,col="blue",type="l")
			}
		}
	} 
}

# plot schedule on 2D plane
plotAlloc <- function(s, jigleSd=0.05,col=1){
	len <- length(s)
	x <- cumsum(s==1)
	y <- cumsum(s==2)
	points(x+rnorm(len,0,jigleSd),y+rnorm(len,0,jigleSd),pch="*",cex=0.75,col=col)
}


#### plot 3-treatment case 
require(scatterplot3d)
btObj <- constructBT(c(14,21,25),optimType="minNA",secondaryProperty="symmetry" )


nGen<-11
w <- attr(btObj,"w")
w.e <- w*nGen/sum(w)
r3d <- scatterplot3d(
	x = c(0,w.e[1]), y = c(0,w.e[3]),z = c(0,w.e[3]),
	col.axis="grey",type="l", angle=60, color="red",
	xlim=c(0,w.e[1]+1),ylim=c(0,w.e[2]+1),zlim=c(0,w.e[3]+1),
	xlab="Treatment A", ylab="Treatment B", zlab="Treatment C",
	cex.axis=1.5, cex.lab=1.5
)

for (g in 1:nGen){
	for (nodeI in 1:dim(btObj[[g]]$nodeMat)[2]){
		node <- btObj[[g]]$nodeMat[,nodeI]
		for ( i in 1:3){   # connect node to next generation nodes
			cnodeI <- btObj[[g]]$mapToNextGenerNode[3*(nodeI-1)+i] # index of connected node
			if (is.na(cnodeI)) next
			cnode <- btObj[[g]]$nextGenerNodeMat[,cnodeI]
			x <- c( node[1], cnode[1])
			y <- c( node[2], cnode[2])
			z <- c( node[3], cnode[3])
			if (g<nGen)
				r3d$points3d(x,y,z,col="blue",type="l")
			else 
				r3d$points3d(node[1],node[2],node[3],col="red",type="p",pch=16,cex=2)
		}
	}
} 
