require(Rglpk)
setwd( "/home/yt/Projects/BTR")
source("./code/constructBT_v1.R")

##############################################################
numSchedules <- 5
C1 <- 2
C2 <- 3
b  <- 4.5
delta <- c(1)
ver <- "WBT_2_3_delta_1"
sLen <- 15
simResDir <- "WBTsimRes"

#############################################################
defineSwitchRules <- function(C1,C2,b,seqLenthCap){
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
	dat0$isTB <- dat0$bound == "t" & c(as.character(dat0$bound[2:length(dat0$bound)]),NA)=="b" 
	dat0$badPair <- (dat0$isConseq & dat0$isTB) # 1-st element of bad pair
	switchRule <- unique(subset(dat0,select=c("layer","pos","badPair")))
	return(switchRule)
	# (layer, pos, badPair), where badPair is the indicator of the 1-st element of 
	# the 'bad Pair' of swithces (i.e. only one switch from the pair should be executed 
	# to avoid leaving bound of the WBT)
}

# apply switch rules to BT schedule
applySwitches <- function(
	sch,      # input BT schedule
	sr,   	  # dataset that defines switch rules
	delta = 0.5 # switch probability
){
	if (length(delta)>1 ) {
		if (length(delta)!=length(unique(sr$layer))) {
			stop( "Mismatch between length of delta and number of layers")
		}
	} else {
		delta <- rep(delta,length(unique(sr$layer)))			
	}
	len <- dim(sr)[1]
	switchDelta <- delta[sr$layer]
	allBadPairInd <- c(which(sr$badPair),which(sr$badPair)+1)
	# double delta in bad pairs if delta is less than 1/2 and set detla to 1 otherwise
	switchDelta[allBadPairInd] <- ifelse(switchDelta[allBadPairInd]>=0.5,1,2*switchDelta[allBadPairInd])
	# execute at most one switch from the pair
	pickElementFromBadPair <- which(sr$badPair)+ ifelse( runif( sum(sr$badPair) )<=0.5,0,1)
	switchDelta[pickElementFromBadPair] <- 0

	#len <- dim(sr)[1]
	#switchDelta <- rep(delta,len)
	#switchDelta[ c(which(sr$badPair),which(sr$badPair)+1)] <- 0 # reset bad pair to 0 
	#pickElementFromBadPair <- which(sr$badPair)+ ifelse( runif( sum(sr$badPair) )<=0.5,0,1) 
	#switchDelta[pickElementFromBadPair] <- delta*ifelse(delta>=0.5,1,2)
	
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

#############################################################

# pick scenario to be run by a processer (worker)
job.id <- as.integer(Sys.getenv("SGE_TASK_ID"))
if(is.na(job.id)) job.id <- 1

switchRule <-defineSwitchRules(C1,C2,b=b,sLen)
allocRule <- constructBT(c(C1,C2))

set.seed( 777 + job.id )
# generate BT
S <- array(0,c(numSchedules, sLen ))
for (j in 1:dim(S)[1]) S[j,] <- generateAllocSeq(allocRule,sLen) 
# generate WBT and save results
S.wbt <- t(apply( S,1,applySwitches,sr=switchRule,delta=delta))

write.table( S.wbt,file=file.path(simResDir,sprintf("%s_job_%d.out",ver,job.id)),
	row.names=FALSE,col.names=FALSE, eol="\r\n")
if (job.id==1) 
	save(C1,C2,b,delta,sLen,ver,simResDir,file = file.path(simResDir,sprintf("%s_param.Rdata",ver)))
