## merge all sample readStat
## Mathieu bourgey 
## 2013/01/31

args=commandArgs(TRUE)
paternFile=args[1]
fileDir=args[2]
outputFile=args[3]

resTable=NULL
listFile=list.files(fileDir,pattern=paternFile,recursive=F)
for(i in 1:length(listFile)) {
	stats<-read.table(paste(fileDir,listfile[i],sep="/"), header=F, sep=",")
	resTable=rbind(resTable,stats[1,])
}

finalTable=cbind(resTable[,1:3],round((resTable[,3]/resTable[,2])*100,2),resTable[,4],round((resTable[,4]/resTable[,3])*100,2))
colnames(finaleTable)=c("Samples","Raw","Filtered","Filtered%","Aligned","Aligned%")

write.table(finaleTable,file=outputFile,sep="\t",row.names=F,col.names=T,quote=F)