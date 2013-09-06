## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

args=commandArgs(TRUE)
paternFile=args[1]
fileDir=args[2]
outputFile=args[3]
type=args[4]

listFile=list.files(fileDir,pattern=paternFile,recursive=T)
nameSample=strsplit(listFile,"/")
name=NULL
info=NULL
rawR=NULL
trimR=NULL
for(i in 1:length(listFile)) {
	name=c(name,nameSample[[i]][1])
	stats<-read.table(file.path(fileDir,listFile[i]), header=F, sep=",")
	if (type == "paired") {
		ligne=grep("Fragment Surviving",stats$V1)
	} else if (type == "single") {
		ligne=grep("Single Surviving",stats$V1)
	} else {
		print("unrocognized library type use of the paired mode")
		ligne=grep("Fragment Surviving",stats$V1)
	}
	rawR=c(rawR,stats[1,2])
	trimR=c(trimR,stats[ligne,2])
	info=c(info,nameSample[[i]][2])
}
finalTable=cbind(name,info,rawR,trimR)
write.table(finalTable,file=outputFile,sep="\t",row.names=F,col.names=F,quote=F)
