# experiment with definition of tunnel geometry


require(Rglpk)

# compute geometry of BT, i.e all feasible nodes in each generations
getBTgeometry <- function(w,numExtraLayers=1){
   if (! all(is.numeric(w)) && all(round(w, 0) == w)){
       stop("Randomization weights must an integer vector")
   }
   K <- length(w)
   allSubsets <- sapply(1:K,function(x)combn(1:K,x,simplify=F))
   cube <- sapply(unlist(allSubsets,F),function(x){ 
       a<-rep(0,K) 
       a[x] <- 1
       return(a)
   })
   cube <- cbind(rep(0,K),cube)
   r <- c(0,sort(unique(unlist(lapply(w,function(x){ (1:(x-1))/x})))))
   cubeConers <- sapply(r,function(x){ floor(x*w) })
   BTnodes <- unique(matrix(apply(cubeConers,2,function(x){ x+cube }),nrow=K),MARGIN=2)
   if (numExtraLayers>0){
	   # put additional layer recoded as coordinates of cubes closest to the origin coners
	   extraLayerUp <- unique(do.call(cbind,
	   		lapply(1:1,function(dm){ ## TODO change 1:1 to 1:K
				   shiftedInADim <- BTnodes
				   shiftedInADim[row(shiftedInADim)==dm] <- shiftedInADim[row(shiftedInADim)==dm] + 1
				   shiftedInADim
			   })
	   ),MARGIN=2)
		extraLayerDown <- unique(do.call(cbind,
			lapply(1:1,function(dm){ ## TODO change 1:1 to 1:K
					shiftedInADim <- BTnodes
					shiftedInADim[row(shiftedInADim)==dm] <- shiftedInADim[row(shiftedInADim)==dm] - 1
					shiftedInADim
			})
		),MARGIN=2)
		extraLayerDown <- extraLayerDown[,which( apply(extraLayerDown,2, function(x) all(x>=0)))]

	   BTnodes <- unique(cbind(BTnodes,extraLayerUp,extraLayerDown),MARGIN=2)
   }
   gen <- apply(BTnodes,2,sum)
   list(gen=sort(gen),bt=BTnodes[,order(gen)])
}


constructBT <-function(w, m=NULL, optimType="minED", distType="Euclidean", secondaryProperty="none"){
   if (! all(is.numeric(w)) && all(round(w, 0) == w)){
       stop("Randomization weights must an integer vector")
   }
   if (is.null(m)) {m<-sum(w)+1}
   if (!is.element(optimType, c("minED", "minFarthest", "minNA","maxED"))) {
       stop("Character specification of optimization type may only be minED, minminFarthest or minNA")
   }
   if (!is.element(secondaryProperty, c("symmetry", "poolToOneHalf", "none"))) {
       stop("Character specification for additional properties of randomization probabilaties type may only be symmetry, poolToOneHalf or none")
   }
   if (!is.element(distType, c("Euclidean", "norm1"))) {
       stop("Character specification of optimization type may only be Euclidean or norm1 ")
   }
   distance <- switch (distType, 
       Euclidean = function(x,r){
          proj <- crossprod(x,r)/crossprod(r) * r
          D <- sqrt(crossprod(x - proj))
       },  
       norm1 = function(x,r){
          perfPoint <- r*sum(x)/sum(r) 
          D <- sum(abs(x - perfPoint))
       }  
   )
   tg <- getBTgeometry(w)
 
   getNextGeneration <-function(curGener, nodeMat, nodePr, w, nodeNumInPaths){
   # w - randomization ratio
   # curGener -(subject) for nextGenerNodeMat and nextGenerNodePr are going to be computed
   # for the next stage (subject)
   # each column of nodeMat represent a node defined by number of subject at 
   # each treatment arm ( j-th row corresponds to the j-th treantment arm)                             
   
                                # compute possible nodeMat for next stage
   # for each column of the nodeMat generatate K new nodes (b/c the next subject 
   # is assigned to one of the K possible treatments
 
   K <-length(w)       # number of treatment
   numNodes <- dim(nodeMat)[2]  # number of nodes at the current stage
   # matrix of all possible nodes for the next generation with the same structure as nodeMat 
   candidateNodeMat <- matrix(apply(nodeMat,2, function(x){
      matrix(rep(x,K),nrow=K)+diag(K)
   }),nrow=K)
   notFeasibleNode <- is.na(match(
      apply(candidateNodeMat,2,function(x)paste(x,collapse='.')),
      apply(tg$bt[,tg$gen == (curGener%%sum(w))+1,drop=F]+w*(curGener%/%sum(w)),2,function(x)paste(x,collapse='.'))
   ))
#############################################################################################
   # compute transition probabilities to the candidate nodes at the next stage
   # by setting system of linear equations A * (trans.prob) = b
   
   # transition probabilities from each current node should add to 1
   A.sumTo1 <- kronecker( diag(numNodes) ,matrix(rep(1,K),nrow=1))  
   b.sumTo1 <- rep(1,numNodes)
   # set impossible transition b/c of the tunnel restriction
   ind.p0 <- which(notFeasibleNode)
   if (length(ind.p0)>0){
      A.p0 <- matrix(0,ncol= numNodes * K, nrow =length(ind.p0) )
      A.p0[cbind(1:length(ind.p0),ind.p0) ] <- 1
      b.p0 <- rep(0,length(ind.p0))
   } else {
      A.p0 <- NULL
      b.p0 <- NULL
   }   
   # transition prob. should be such that unconditionally requirement as for w is met
   A.m <- matrix(sapply(nodePr,function(x,K){x*diag(K)},K=K),nrow=K)
   b.m <- w/sum(w)
   # find trans. probabilities    
   A <- rbind( A.sumTo1, A.m, A.p0 )
   b <- matrix(c(b.sumTo1, b.m, b.p0))
   
   # use linear programming to solve under-defined system  
   nv <- dim(A)[2] 
  
   distToRay <- sapply(as.data.frame(candidateNodeMat), distance, r = w)
   distToRay[notFeasibleNode]<- -1
   farthestToRay <- distToRay==max(distToRay)
 
   if (secondaryProperty=='symmetry'){
       Apd <- cbind( rep(1,nv-1), -diag(nv-1))
       if (nv>2)
       for(u in 2:(nv-1)){
          Apd <- rbind(Apd, cbind(  matrix(0,ncol=u-1,nrow=nv-u),rep(1,nv-u),-diag(nv-u))  )
       }
       nvApd <- nv*(nv-1)/2
       A.abs.1 <- cbind(   Apd, -diag(nvApd) )
       A.abs.2 <- cbind(  -Apd, -diag(nvApd) ) 
       b.abs.1 <-  rep(0,nvApd)
       b.abs.2 <-  rep(0,nvApd)
       A.sec <- rbind(A.abs.1, A.abs.2)
       b.sec <- c(b.abs.1, b.abs.2)
       dir.sec <- rep("<=", 2*nvApd) 
       numAuxVar <- nvApd
   } else if (secondaryProperty=='poolToOneHalf') {
       A.abs.1 <- cbind(  diag(nv), -diag(nv) )
       A.abs.2 <- cbind( -diag(nv), -diag(nv) ) 
       b.abs.1 <- rep( 0.5, length.out=nv)
       b.abs.2 <- rep(-0.5, length.out=nv)
       A.sec <- rbind(A.abs.1, A.abs.2)
       b.sec <- c(b.abs.1, b.abs.2)
       dir.sec <- rep("<=", 2*nv) 
       numAuxVar <- nv
   } else {
       A.sec <- NULL
       b.sec <- NULL
       dir.sec <- NULL 
       numAuxVar <- 0
   }
   # define constraint matrix for linear optimization program
   A.lp <- rbind(
       cbind(A, matrix(0,nrow=dim(A)[1],ncol=numAuxVar)),  # tunel constraints, probabilities sums to 1 and uncondinional prob. requirements
       cbind(diag(nv), matrix(0,nrow=nv,ncol=numAuxVar)),  # to restrict probab. be less then 1 
       A.sec                                               # to incorporate secondary properties, e.g. symmetry
   ) 
   # corresponding right hand side of A.lp
   rhs.lp <- c( b,
      rep(1,nv),
      b.sec
   )
   # corresponding directions
   dir.lp <- c( rep("==",length(b)),  
      rep("<=",nv),         
      dir.sec       
   )
   # define objective for lp  
   objective<-switch(optimType,
       minED = c( 10000*rep(nodePr,each=K)*distToRay, rep(1,numAuxVar )),
       minFarthest = c( 10000*rep(nodePr,each=K) * farthestToRay, rep(1,numAuxVar )),
       minNA = c( rep(0,nv), rep(1,numAuxVar )),
       maxED = c( -10000*rep(nodePr,each=K)*distToRay, rep(1,numAuxVar ))

   )
   lpRes <- Rglpk_solve_LP( obj=objective, mat=A.lp, dir=dir.lp, rhs=rhs.lp,
       types = NULL, max = FALSE, bounds = NULL, verbose = FALSE) 
   if ( lpRes$status == 0){
       trProb <- lpRes$solution[1:nv]
   } else {
       stop(paste("Solution was not found by linear programming solver. \n",
           "Attepted to solve problem for", nv, "transition probabilities, plus", numAuxVar,"auxlilary variables."))
   }

   # compute probability of getting to a particular node at next stage (subject)
   # and number of different paths that lead to the node at next stage
   aux.nextGenerNodePr <- rep(nodePr,each=K)*c(trProb)
   aux.nextGenerNodeNumInPaths <- rep(nodeNumInPaths,each=K)*c(trProb>0.00001)
   nextGenerNodeMat <- matrix(unique(candidateNodeMat,MARGIN=2),nrow=K)
   nextGenerNodePr <- numeric(dim(nextGenerNodeMat)[2])
   nextGenerNodeNumInPaths <- numeric(dim(nextGenerNodeMat)[2]) # count number of paths that lead to node
   for (j in 1:dim(nextGenerNodeMat)[2]){
       un.j <- apply(candidateNodeMat,2,function(x,st){
           all(x == st)
       }, st = nextGenerNodeMat[,j])  # get indexes that map to the same actual node
       nextGenerNodePr[j] <- sum(aux.nextGenerNodePr[un.j])
       nextGenerNodeNumInPaths[j] <- sum( aux.nextGenerNodeNumInPaths[un.j] )
   }
   nextGenerNodeMat <- matrix(nextGenerNodeMat[ ,nextGenerNodePr > 0.00001],nrow=K)
   nextGenerNodePr <- nextGenerNodePr[nextGenerNodePr > 0.00001]
   nextGenerNodeNumInPaths <- nextGenerNodeNumInPaths[nextGenerNodeNumInPaths>0]
   return(list(curGener=curGener, nodeMat=nodeMat, 
      nodePr=nodePr, nodeNumInPaths = nodeNumInPaths,
      curNodeRandProb= round(matrix(trProb,nrow=K),6),
      nextGenerNodeMat=nextGenerNodeMat, nextGenerNodePr=nextGenerNodePr, 
      nextGenerNodeNumInPaths = nextGenerNodeNumInPaths,
#     nextGener_ED=  sum( nextGenerNodePr *sapply(as.data.frame(nodeMat),distance,r=w)),
      nextGenerCandidateNodeMat =  candidateNodeMat,
#      nextGenerCandidateNode_distToRay = distToRay,
      # add to information about mapping of candidate nodes to actual nodes for next generation  
      # to speed up actual allocation generation
      mapToNextGenerNode = match( 
         apply(t(candidateNodeMat),1,function(y)paste(y,collapse='')),
         apply(t(nextGenerNodeMat),1,function(y)paste(y,collapse=''))
      )
   ))
   } # getNextGeneration

   # build tunnel 
   tunnel <- list(   # null (initial) gener
     list( generation=0, nodeMat=matrix(rep(0,length(w)),nrow=length(w)),nodePr = 1,nodeNumInPaths=1)
   )
   for( i in 1:m){
     temp <- getNextGeneration(curGener=i-1, tunnel[[i]]$nodeMat, tunnel[[i]]$nodePr, w, tunnel[[i]]$nodeNumInPaths )
     tunnel[[i]]$randProb <- temp$curNodeRandProb
     tunnel[[i]]$nextGenerNodeMat <- temp$nextGenerNodeMat     
	 tunnel[[i]]$nextGenerCandidateNodeMat <- temp$nextGenerCandidateNodeMat  # todo     
     tunnel[[i]]$mapToNextGenerNode <- temp$mapToNextGenerNode   # used for generation of alloc schedule instanses             
     tunnel[[i+1]] <- list(generation=i, nodeMat=temp$nextGenerNodeMat, 
        nodePr = temp$nextGenerNodePr, nodeNumInPaths=temp$nextGenerNodeNumInPaths)
   }
   tunnel<-tunnel[1:m]
   class(tunnel)<-'btRand'
   attr(tunnel,'w')<-w
   attr(tunnel,'distance')<- distance
   return(tunnel)
}

##########################################################################################################
# print btRand object as table 
print.btRand <- function(obj, printAllNodes=T, checkSymmInBlk=F, ...){
    w <- attr(obj,'w')
    K<- length(w)
    if (printAllNodes) BT<- getBTgeometry(attr(obj,'w'))
    tbl<-c()
    for (g in 1:length(obj) )
        tbl<-with(obj[[g]],{ 
            distToRay <- sapply(as.data.frame(nodeMat), attr(obj,'distance'), r=w)
            for ( inode in 1:dim(nodeMat)[2]){
                tblrow <- c(generation, inode,nodeMat[,inode], distToRay[inode], nodePr[inode],
                   nodeNumInPaths[inode],randProb[,inode])
                tbl <- rbind(tbl,tblrow)
            }
            if (printAllNodes){
                allNodes <- BT$bt[,BT$gen == g-1 %% sum(w) ]
                idroppedNodes <- which(!duplicated( cbind(nodeMat,allNodes ),MARGIN=2)[-c(1:dim(nodeMat)[2])])
                if (length(idroppedNodes)>0){
                    distToRay.dropped <- sapply(as.data.frame(allNodes[,idroppedNodes]),attr(obj,'distance'),r=w)
                    for ( j in 1:length(idroppedNodes)){
                        tblrow <- c(generation, inode+j, allNodes[,idroppedNodes[j]],distToRay.dropped[j],0,NA,rep(NA,K))
                        tbl <- rbind(tbl,tblrow)
                    }
                }
            }
            tbl
        }) # with
    colnames(tbl) <- c("gen","node",paste("nodeCoord",1:K,sep='.'),"distToAR","residProb","numPaths",paste("randProb",1:K,sep='.'))
    rownames(tbl) <- NULL
    if (checkSymmInBlk){
        coord <- tbl[,paste("nodeCoord",1:K,sep=".")]
        blkSize_coord <- t(w-t(coord))
        imatched <- match(apply(tbl[,paste("nodeCoord",1:K,sep=".")],1,function(x)paste(x,collapse='.')),  
           apply( blkSize_coord,1, function(x)paste(x,collapse=".")))
        tbl<- cbind( tbl, symmInBlk= round(tbl[,"residProb"][imatched],5)==round(tbl[,"residProb"],5))
    }
#  ED <- tapply( tbl[,"residProb"]*tbl[,"distToAR"],tbl[,"gen"],sum)
    cat("Target allocation ratio",paste(w,collapse=":"),'\n')
    print(tbl,...)
}

# function to generate a rialization of random treatment allocation sequence
generateAllocSeq <- function(allocRule,m=NULL){
    w <- attr(allocRule,'w')
    if (is.null(m)) 
      m <- length(allocRule)-1
    K <- dim(allocRule[[1]]$randProb)[1]
    trtSeq <-numeric(m)
    curNodeInd <- 1
    for (ii in 1:m){
        i <- (ii-1)%%sum(w)+1
        trtSeq[ii] <- findInterval(  runif(1), c(0,cumsum(allocRule[[i]]$randProb[,curNodeInd])) )
          #      curNodeInd <- allocRule[[i]]$mapToNextGenerNode[ (curNodeInd-1)*K + trtSeq[i] ] bug fixed Mar 19 2013
          curNodeInd <- allocRule[[i]]$mapToNextGenerNode[ (curNodeInd-1)*K + trtSeq[ii] ]
          
    }
    return(trtSeq)
}

 # Example 
if(0){
 set.seed(777)
 #allocRule <- constructBT(w=c(4,2,1,1,1),optimType='minNA' )
 allocRule <- constructBT(w=c(4,2,1,1,1),optimType='minED')  
 numSchedules <- 10000
 schedules <- array(0,c(numSchedules,length(allocRule)-1))
 for (j in 1:dim(schedules)[1])
    schedules[j,] <- generateAllocSeq(allocRule,length(allocRule)-1) 

 # look at prelances of rand sequenses
 apply(schedules[,],1,function(x)paste(x,collapse='')) |> 
   table() |>
   sort(decreasing = TRUE) |> data.frame()
 
 # check that ARP holds
 apply(schedules,2,table)/numSchedules 
}

 