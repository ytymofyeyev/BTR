require(Rglpk)
setwd("C:/Projects/BTR")
source("constructBT_v1.R")
source("plotBT_work.R")

C1<-2
C2<-3
b <- 4.5
seqLenthCap <- 15
bBT <- (C2-1)/C1+1
# number of layers is L (WBT might also include L+1 layer)
L <- floor(b-bBT)
m <- 1:seqLenthCap
Ytop_m <- ceiling(m*C2/C1) 
Ybot_m <- floor( (m-1)*C2/C1)
pos <- layer <- bound <- numeric(0)
for( j in 0:L){
	switchPosTop <- m+Ytop_m + j 
	switchPosBot <- m-1+Ybot_m - j
	if (j==L) {# check if node withing WBT
	  checkTop <- Ytop_m+L+1 - (m-1)*C2/C1 <= b
	  checkBot <- m*C2/C1-(Ybot_m-L-1) <= b
	  switchPosTop <- switchPosTop[checkTop]
	  switchPosBot <- switchPosBot[checkBot]
	}
	switchPos    <- c( switchPosTop,switchPosBot) # position of switch
	boundInd <- c(rep("t",length(switchPosTop)),rep("b",length(switchPosBot)) )
    badPairCond <- ceiling((m+1)*C2/C1) - ceiling(m*C2/C1)  >= 2
	
	cutInd <- switchPos < seqLenthCap & switchPos > 1 
	switchPos <- switchPos[cutInd]
	boundInd <- boundInd[cutInd]
	
	# save results from current layer
	pos <-   c(pos, switchPos)
	layer <- c(layer, rep(j+1,length(switchPos)))
	bound <- c(bound, boundInd)
}

dat0 <- data.frame(layer,pos,bound)
dat0 <- dat0[order(dat0$layer,dat0$pos,dat0$bound),]
dat0$isConseq <- dat0$pos+1 == c(dat0$pos[2:length(dat0$pos)],NA) 
dat0$isTB <- dat0$bound == "t" & as.character(dat0$bound[2:length(dat0$bound)],NA)=="b" 
dat0$badPair <- (dat0$isConseq & dat0$isTB) # 1-st element of bad pair

switchRule <- unique(subset(dat0,select=c("layer","pos","badPair")))
switchRule

# apply switches
applySwitches <- function(sch, sr, delta=0.5){
	len <- dim(sr)[1]
	switchDelta <- rep(delta,len)
	switchDelta[ c(which(sr$badPair),which(sr$badPair)+1)] <- 0 # reset bad pair to 0 
	pickElementFromBadPair <- which(sr$badPair)+ ifelse( runif( sum(sr$badPair) )<=0.5,0,1) 
	switchDelta[pickElementFromBadPair] <- delta*ifelse(delta>=0.5,1,2)
	sr$swOn <- runif(len) < switchDelta
	for (l in 1:max(sr$layer)){
		layerSR <- subset(sr, layer==l)
		for (i in 1:length(layerSR$pos) ){
			toSwitch <- c( layerSR$pos[i], layerSR$pos[i]+1 ) 
			if( layerSR$swOn[i]){
				sch[toSwitch] <- sch[rev(toSwitch)]
			}
		}
	}
	return(sch)
}

# generate BT
allocRule <- constructBT(c(C1,C2), m=seqLenthCap)
numSchedules <- 10000
sch <- array(0,c(numSchedules,length(allocRule)))
for (j in 1:dim(sch)[1])
	sch[j,] <- generateAllocSeq(allocRule) 

plotBT.2D(allocRule)
grid()
abline(b=C2/C1,a=0)
#apply( sch,1,plotAlloc,col=1)
sch.WBT.1 <- t(apply( sch,1,applySwitches,sr=switchRule,delta=0.8))
#apply( sch.WBT.1, 1, plotAlloc,col=1)
sch.WBT.2 <- t(apply( sch,1,applySwitches,sr=switchRule,delta=0.8))
#apply( sch.WBT.2, 1, plotAlloc,col=2)
#legend("topleft",c("delta=0.9","delta=0.3"),col=1:2,pch='*')


apply(sch.WBT.1, 2,table)/numSchedules  #marginalProb 
apply(sch.WBT.2, 2,table)/numSchedules  #marginalProb 
 

getResidentPr <- function(S){
	nsim <- dim(S)[1]
	len  <- dim(S)[2]
	A <- S==1
	cA <- apply(A,1,cumsum)
	rProb <- lapply(apply(cA,1,table),function(x){
		x/nsim
	})
#	len <- length(x)	
#	sgn <- sign( (1:len) - mean(1:len))
#	names(x)<-sgn*ceiling(abs((1:len) - mean(1:len)))
}

rPr.BT <- getResidentPr(sch)
rPr.WBT.1 <- getResidentPr(sch.WBT.1)
rPr.WBT.2 <- getResidentPr(sch.WBT.2)



########################
# apply switches
if(0){
# merge unique and duplicated enries wrt 'layer' and 'pos'
uInd <- which(!duplicated(subset(dat0,select=c("layer","pos"))))
switch <- merge(
		dat0[ uInd, c("layer","pos","badPair")],
		dat0[ setdiff(1:dim(dat0)[1] ,uInd),c("layer","pos","badPair")], by=c("layer","pos"),all=T,
		suffix=c("",".2rec")
)
# over-write "badPair" indicator if same 'pos' appeared twice
switch$badPair[ switch$badPair.2rec==TRUE] <- TRUE
switch
}





#switchPosTop <- m+Ytop_m 
#switchPosBot <- m-1+Ybot_m
#sort(union(switch1Top,switch1Bot))

 # Example 
if(0){
 set.seed(777)
 allocRule <- constructBT(w=c(4,2,1,1,1),optimType='minNA' ) 
 numSchedules <- 10000
 schedules <- array(0,c(numSchedules,length(allocRule)))
 for (j in 1:dim(schedules)[1])
    schedules[j,] <- generateAllocSeq(allocRule) 

 # to report number of unique schedules and number of occurances
 summary(as.data.frame(apply(schedules[,],1,function(x)paste(x,collapse=''))),maxsum=20) 
 ( marginalProb <- apply(schedules,2,table)/numSchedules )
}


