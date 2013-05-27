## lading package
library(parallel)

##get file names from args
args=commandArgs(TRUE)
count_file=args[1]
gtf_file=args[2]
path_outR=args[3]
path_outS=args[4]
threadNum=as.integer(args[5])
if (threadNum > 19) {
	threadNum=19
}
DoSat=args[6]
##Do sat values 0=RPKM ; 1=RPKM+SATURATION ; 2=SATURATION

## open files
co=read.table(count_file,header=T,sep="\t")
gtf=read.table(gtf_file,header=F,sep="\t")

##rpkm
leT=cbind(gtf[,2])
rownames(leT)=gtf[,1]
coA=cbind(co[,3:dim(co)[2]])
rownames(coA)=co[,1]
leO=match(rownames(coA),rownames(leT))
letO=cbind(leT[leO,]/1000)
coK=cbind(coA/letO)
Mr=colSums(coA)/1000000
coKM=NULL
for (i in 1:dim(coK)[2]) {
    coKM=cbind(coKM,coK[,i]/Mr[i])
}
rownames(coKM)=rownames(coK)
colnames(coKM)=colnames(coK)
write.table(round(coKM,3),paste(path_outR,"matrixRPKM.tsv",sep="/"),sep="\t",quote=F,col.names=T,row.names=T)


##saturation
if (DoSat == "1") {
    satModel=function(te,s,qs,qf,q1,q2,q3,q4) {
	readNum=round(length(te)*s)
	res=NULL
	resP=NULL
	for (k in 1:100) { ## generate 100 replicats by sub coverage
	repV=table(sample(te,s=readNum,replace=F))
	repR=rep(0,length(qs))
	repR[as.numeric(names(repV))]=repV
	repRK=repR/qs
	repRKM=repRK/(readNum/1000000)
	repPRE=((abs(repRKM-qf)/(qf+0.000001))*100)
	res=cbind(res,repRKM)
	resP=cbind(resP,repPRE)
	}
	res1=apply(res[q1,],2,median)
	res2=apply(res[q2,],2,median)
	res3=apply(res[q3,],2,median)
	res4=apply(res[q4,],2,median)
	res1p=colMeans(resP[q1,])
	res2p=colMeans(resP[q2,])
	res3p=colMeans(resP[q3,])
	res4p=colMeans(resP[q4,])
	return(cbind(res1,res2,res3,res4,res1p,res2p,res3p,res4p))
    }
    satM=seq(0.05,0.95,by=0.05)
    satMinv=seq(0.95,0.05,by=-0.05)
    steps=floor(19/threadNum)
    subS=list()
    cpS=1
    for (i in 1:steps) {
        subS[[i]]=vector()
        for (j in 1:threadNum) {
             subS[[i]]=c(subS[[i]],cpS)
             cpS=cpS+1
        }
    }
    if (cpS <= 19) {
        subS[[steps+1]]=cpS:19
	steps=steps+1
    }
    for (i in 1:dim(coKM)[2]) { ## iterate over samples
        ## generate saturation empty results vectors
        satRQ1=NULL
        satRQ2=NULL
        satRQ3=NULL
        satRQ4=NULL
	satRQ1p=NULL
        satRQ2p=NULL
        satRQ3p=NULL
        satRQ4p=NULL
	#pos=1:dim(co2)[1]
        ## generate quantile contigs/transcripts set
	coKM2=as.numeric(as.vector(coKM[coKM[,i] > 0,i]))
	coKMZ=as.numeric(as.vector(coKM[coKM[,i] == 0,i]))
	co2=as.numeric(as.vector(co[coKM[,i] > 0,i+2]))
	letO2=as.numeric(as.vector(letO[coKM[,i] > 0,1]))
	pos=1:length(co2)
	fpkmZP=round(round(length(coKMZ)/(length(coKMZ)+length(coKM2)),2)*100)
        limRPKM=quantile(coKM2,p=c(0.25,0.5,0.75))
        ql1=pos[coKM2 <= limRPKM[1]]
        ql2=pos[coKM2 > limRPKM[1] & coKM2 <= limRPKM[2]]
        ql3=pos[coKM2 > limRPKM[2] & coKM2 <= limRPKM[3]]
        ql4=pos[coKM2 > limRPKM[3]]
        ## generate presampling vector ... Warnings: it needs at least 10G of ram
        testV=rep(as.vector(1:length(co2)),co2)
        ## resampling
	resAll=NULL
	for (k in 1:steps) {
		for (j in subS[[k]]) { ## iterate over sub-coverage
			cmd=expression(satModel(testV,satM[j],letO2,coKM2,ql1,ql2,ql3,ql4))
			mcparallel(cmd)
		}
		tmpres=c(mccollect(),resAll)
		resAll=tmpres
	}
	for (j in 1:length(satM)) {
		satRQ1=cbind(satRQ1,resAll[[j]][,1])
		satRQ2=cbind(satRQ2,resAll[[j]][,2])
		satRQ3=cbind(satRQ3,resAll[[j]][,3])
		satRQ4=cbind(satRQ4,resAll[[j]][,4])
		satRQ1p=cbind(satRQ1p,resAll[[j]][,5])
		satRQ2p=cbind(satRQ2p,resAll[[j]][,6])
		satRQ3p=cbind(satRQ3p,resAll[[j]][,7])
		satRQ4p=cbind(satRQ4p,resAll[[j]][,8])
	}
	colnames(satRQ1)=as.character(round(satMinv*100))
	colnames(satRQ2)=as.character(round(satMinv*100))
	colnames(satRQ3)=as.character(round(satMinv*100))
	colnames(satRQ4)=as.character(round(satMinv*100))
	colnames(satRQ1p)=as.character(round(satMinv*100))
	colnames(satRQ2p)=as.character(round(satMinv*100))
	colnames(satRQ3p)=as.character(round(satMinv*100))
	colnames(satRQ4p)=as.character(round(satMinv*100))
        jpeg(paste(path_outS,paste(colnames(coKM)[i],"saturation.jpeg",sep="_"),sep="/"),1000,1000)
	par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1,1,3,1))
        boxplot(satRQ4,main="Q4 saturation",xlab="Resampling precentage",ylab="median RPKM")
	boxplot(satRQ2,main="Q2 saturation",xlab="Resampling precentage",ylab="median RPKM")
        boxplot(satRQ3,main="Q3 saturation",xlab="Resampling precentage",ylab="median RPKM")
	boxplot(satRQ1,main="Q1 saturation",xlab="Resampling precentage",ylab="median RPKM")
	mtext(paste("Saturation estimate of the mean FPKM - Excluding",as.character(fpkmZP),"perc of genes with fpkm = 0",sep=" "),NORTH<-3, line=1, adj=0.5, cex=1.2, outer=TRUE, col="black")
        dev.off()
	jpeg(paste(path_outS,paste(colnames(coKM)[i],"PRE_saturation.jpeg",sep="_"),sep="/"),1000,1000)
        par(mfcol=c(2,2), mar=c(4,4,0.5,0.5), oma=c(1,1,4,1))
	yl1=yl2=yl3=yl4=1
	yl1=max(satRQ1p)+5
	yl2=max(satRQ2p)+5
	yl3=max(satRQ3p)+5
	yl4=max(satRQ4p)+5
	print(max(satRQ1p))
	print(max(satRQ2p))
	print(max(satRQ3p))
	print(max(satRQ4p))
        boxplot(satRQ4p,main="Q4 saturation",xlab="Resampling precentage",ylab="mean PRE",ylim=c(0,yl4))
	boxplot(satRQ2p,main="Q2 saturation",xlab="Resampling precentage",ylab="mean PRE",ylim=c(0,yl2))
        boxplot(satRQ3p,main="Q3 saturation",xlab="Resampling precentage",ylab="mean PRE",ylim=c(0,yl3))
	boxplot(satRQ1p,main="Q1 saturation",xlab="Resampling precentage",ylab="mean PRE",ylim=c(0,yl1))
	mtext(paste("Saturation estimate of the Percent Relative Error - Excluding",as.character(fpkmZP),"perc of genes with fpkm = 0",sep=" "),NORTH<-3, line=1, adj=0.5, cex=1.2, outer=TRUE, col="black")
        dev.off()
    }
}

