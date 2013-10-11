## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

args=commandArgs(TRUE)
fileDir=args[1]
outputFile=args[2]
type=args[3]

name=NULL
align=NULL
duplicate=NULL
paternFile=c(".flagstat",".insert_size_metrics",".all.coverage.sample_summary",".CCDS.coverage.sample_summary")
## get flagstat metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paternFile[1],recursive=T))
sampleNum=length(listFile)
nameSample=strsplit(basename(listFile),".",fixed=T)
for(i in 1:length(listFile)) {
	name=c(name,nameSample[[i]][1])
	stats=scan(file=listFile[i],what="character",sep="\n")
	align=c(align,strsplit(stats[grep("mapped (",stats,fixed=T)]," + ",fixed=T)[[1]][1])
	duplicate=c(duplicate,strsplit(stats[grep("duplicates",stats,fixed=T)]," + ",fixed=T)[[1]][1])
}
pairOrient=rep("NA",sampleNum)
medianInsS=rep("NA",sampleNum)
meanInsS=rep("NA",sampleNum)
averageDev=rep("NA",sampleNum)
standD=rep("NA",sampleNum)
wgMeanCov=rep("NA",sampleNum)
wgbase10=rep("NA",sampleNum)
wgbase25=rep("NA",sampleNum)
wgbase50=rep("NA",sampleNum)
wgbase75=rep("NA",sampleNum)
wgbase100=rep("NA",sampleNum)
wgbase500=rep("NA",sampleNum)
ccdsMeanCov=rep("NA",sampleNum)
ccdsbase10=rep("NA",sampleNum)
ccdsbase25=rep("NA",sampleNum)
ccdsbase50=rep("NA",sampleNum)
ccdsbase75=rep("NA",sampleNum)
ccdsbase100=rep("NA",sampleNum)
ccdsbase500=rep("NA",sampleNum)
## get insert size metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paternFile[2],recursive=T))
nameSample=strsplit(basename(listFile),".",fixed=T)
if (sampleNum == length(listFile)) {
	sampleNum=length(listFile)
	for(i in 1:length(listFile)) {
		nameLoc=match(nameSample[[i]][1],name)
		stats=scan(file=listFile[i],what="character",sep="\n")
		statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
		statsValue=strsplit(statsLigne,"\t",fixed=T)
		pairOrient[nameLoc]=statsValue[[1]][8]
		medianInsS[nameLoc]=statsValue[[1]][1]
		meanInsS[nameLoc]=statsValue[[1]][5]
		averageDev[nameLoc]=statsValue[[1]][2]
		standD[nameLoc]=statsValue[[1]][6]
	}
} else {
	sampleNumTmp=length(listFile)
	for(i in 1:length(listFile)) {
		if (nameSample[[i]][1] %in% name) {
			nameLoc=match(nameSample[[i]][1],name)
			stats=scan(file=listFile[i],what="character",sep="\n")
			statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			pairOrient[nameLoc]=statsValue[[1]][8]
			medianInsS[nameLoc]=statsValue[[1]][1]
			meanInsS[nameLoc]=statsValue[[1]][5]
			averageDev[nameLoc]=statsValue[[1]][2]
			standD[nameLoc]=statsValue[[1]][6]
		} else {
			sampleNum=sampleNum+1
			stats=scan(file=listFile[i],what="character",sep="\n")
			name=c(name,nameSample[[i]][1])
			align=c(align,"NA")
			duplicate=c(duplicate,"NA")
			statsLigne=stats[grep("WIDTH_OF_10_PERCENT",stats,fixed=T)+1]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			pairOrient=c(pairOrient,statsValue[[1]][8])
			medianInsS=c(medianInsS,statsValue[[1]][1])
			meanInsS=c(meanInsS,statsValue[[1]][5])
			averageDev=c(averageDev,statsValue[[1]][2])
			standD=c(standD,statsValue[[1]][6])
			wgMeanCov=c(wgMeanCov,"NA")
			wgbase10=c(wgbase10,"NA")
			wgbase25=c(wgbase25,"NA")
			wgbase50=c(wgbase50,"NA")
			wgbase75=c(wgbase75,"NA")
			wgbase100=c(wgbase100,"NA")
			wgbase500=c(wgbase500,"NA")
			ccdsMeanCov=c(ccdsMeanCov,"NA")
			ccdsbase10=c(ccdsbase100,"NA")
			ccdsbase25=c(ccdsbase25,"NA")
			ccdsbase50=c(ccdsbase50,"NA")
			ccdsbase75=c(ccdsbase75,"NA")
			ccdsbase100=c(ccdsbase100,"NA")
			ccdsbase500=c(ccdsbase500,"NA")
		}
	}
}

## get WG coverage metrics
if (type == "wholeGenome") {
	listFile=file.path(fileDir,list.files(fileDir,pattern=paternFile[3],recursive=T))
	nameSample=strsplit(basename(listFile),".",fixed=T)
	if (sampleNum == length(listFile)) {
		sampleNum=length(listFile)
		for(i in 1:length(listFile)) {
			nameLoc=match(nameSample[[i]][1],name)
			stats=scan(file=listFile[i],what="character",sep="\n")
			statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			wgMeanCov[nameLoc]=statsValue[[1]][3]
			wgbase10[nameLoc]=statsValue[[1]][7]
			wgbase25[nameLoc]=statsValue[[1]][8]
			wgbase50[nameLoc]=statsValue[[1]][9]
			wgbase75[nameLoc]=statsValue[[1]][10]
			wgbase100[nameLoc]=statsValue[[1]][11]
			wgbase500[nameLoc]=statsValue[[1]][12]
		}
	} else {
		sampleNumTmp=length(listFile)
		for(i in 1:length(listFile)) {
			if (nameSample[[i]][1] %in% name) {
				nameLoc=match(nameSample[[i]][1],name)
				stats=scan(file=listFile[i],what="character",sep="\n")
				statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
				statsValue=strsplit(statsLigne,"\t",fixed=T)
				wgMeanCov[nameLoc]=statsValue[[1]][3]
				wgbase10[nameLoc]=statsValue[[1]][7]
				wgbase25[nameLoc]=statsValue[[1]][8]
				wgbase50[nameLoc]=statsValue[[1]][9]
				wgbase75[nameLoc]=statsValue[[1]][10]
				wgbase100[nameLoc]=statsValue[[1]][11]
				wgbase500[nameLoc]=statsValue[[1]][12]
			} else {
				sampleNum=sampleNum+1
				stats=scan(file=listFile[i],what="character",sep="\n")
				name=c(name,nameSample[[i]][1])
				align=c(align,"NA")
				duplicate=c(duplicate,"NA")
				pairOrient=c(pairOrient,"NA")
				medianInsS=c(medianInsS,"NA")
				meanInsS=c(meanInsS,"NA")
				averageDev=c(averageDev,"NA")
				standD=c(standD,"NA")
				statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
				statsValue=strsplit(statsLigne,"\t",fixed=T)
				wgMeanCov=c(wgMeanCov,statsValue[[1]][3])
				wgbase10=c(wgbase10,statsValue[[1]][7])
				wgbase25=c(wgbase25,statsValue[[1]][8])
				wgbase50=c(wgbase50,statsValue[[1]][9])
				wgbase75=c(wgbase75,statsValue[[1]][10])
				wgbase100=c(wgbase100,statsValue[[1]][11])
				wgbase500=c(wgbase500,statsValue[[1]][12])
			}
		}
	}
}
## get CCDS coverage metrics
listFile=file.path(fileDir,list.files(fileDir,pattern=paternFile[4],recursive=T))
nameSample=strsplit(basename(listFile),".",fixed=T)
if (sampleNum == length(listFile)) {
	sampleNum=length(listFile)
	for(i in 1:length(listFile)) {
		nameLoc=match(nameSample[[i]][1],name)
		stats=scan(file=listFile[i],what="character",sep="\n")
		statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
		statsValue=strsplit(statsLigne,"\t",fixed=T)
		ccdsMeanCov[nameLoc]=statsValue[[1]][3]
		ccdsbase10[nameLoc]=statsValue[[1]][7]
		ccdsbase25[nameLoc]=statsValue[[1]][8]
		ccdsbase50[nameLoc]=statsValue[[1]][9]
		ccdsbase75[nameLoc]=statsValue[[1]][10]
		ccdsbase100[nameLoc]=statsValue[[1]][11]
		ccdsbase500[nameLoc]=statsValue[[1]][12]
	}
} else {
	sampleNumTmp=length(listFile)
	for(i in 1:length(listFile)) {
		if (nameSample[[i]][1] %in% name) {
			nameLoc=match(nameSample[[i]][1],name)
			stats=scan(file=listFile[i],what="character",sep="\n")
			statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			ccdsMeanCov[nameLoc]=statsValue[[1]][3]
			ccdsbase10[nameLoc]=statsValue[[1]][7]
			ccdsbase25[nameLoc]=statsValue[[1]][8]
			ccdsbase50[nameLoc]=statsValue[[1]][9]
			ccdsbase75[nameLoc]=statsValue[[1]][10]
			ccdsbase100[nameLoc]=statsValue[[1]][11]
			ccdswgbase500[nameLoc]=statsValue[[1]][12]
			
		} else {
			sampleNum=sampleNum+1
			stats=scan(file=listFile[i],what="character",sep="\n")
			name=c(name,nameSample[[i]][1])
			align=c(align,"NA")
			duplicate=c(duplicate,"NA")
			pairOrient=c(pairOrient,"NA")
			medianInsS=c(medianInsS,"NA")
			meanInsS=c(meanInsS,"NA")
			averageDev=c(averageDev,"NA")
			standD=c(standD,"NA")
			wgMeanCov=c(wgMeanCov,"NA")
			wgbase10=c(wgbase10,"NA")
			wgbase25=c(wgbase25,"NA")
			wgbase50=c(wgbase50,"NA")
			wgbase75=c(wgbase75,"NA")
			wgbase100=c(wgbase100,"NA")
			wgbase500=c(wgbase500,"NA")
			statsLigne=stats[grep(nameSample[[i]][1],stats,fixed=T)]
			statsValue=strsplit(statsLigne,"\t",fixed=T)
			ccdsMeanCov=c(ccdsMeanCov,statsValue[[1]][3])
			ccdsbase10=c(ccdsbase10,statsValue[[1]][7])
			ccdsbase25=c(ccdsbase25,statsValue[[1]][8])
			ccdsbase50=c(ccdsbase50,statsValue[[1]][9])
			ccdsbase75=c(ccdsbase75,statsValue[[1]][10])
			ccdsbase100=c(ccdsbase100,statsValue[[1]][11])
			ccdsbase500=c(ccdsbase500,statsValue[[1]][12])
		}
	}
}
finalTable=cbind(name,align,duplicate,pairOrient,medianInsS,meanInsS,averageDev,standD,wgMeanCov,wgbase10,wgbase25,wgbase50,wgbase75,wgbase100,wgbase500,ccdsMeanCov,ccdsbase10,ccdsbase25,ccdsbase50,ccdsbase75,ccdsbase100,ccdsbase500)
colnames(finalTable)=c("SampleName","Aligned","Duplicates","Pair Orientation","Median Insert Size","Mean Insert Size","Average Deviation","Standard Deviation","WG Mean Coverage","WG %_bases_above_10","WG %_bases_above_25","WG %_bases_above_50","WG %_bases_above_75","WG %_bases_above_100","WG %_bases_above_500","CCDS Mean Coverage","CCDS %_bases_above_10","CCDS %_bases_above_25","CCDS %_bases_above_50","CCDS %_bases_above_75","CCDS %_bases_above_100","CCDS %_bases_above_500")
write.table(finalTable,file=outputFile,sep="\t",row.names=F,col.names=T,quote=F)
