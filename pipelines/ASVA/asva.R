suppressMessages(library(dada2, lib.loc="/home/egonzale/dada2_test/libraries_ASVA"))
paste("Dada2 version: ", packageVersion("dada2"))
suppressMessages(library(ShortRead))
paste("ShortRead version: ", packageVersion("ShortRead"))
suppressMessages(library(ggplot2))
paste("ggplot2 version: ", packageVersion("ggplot2"))
suppressMessages(library(phyloseq, lib.loc="/home/egonzale/dada2_test/libraries_ASVA"))
paste("Phyloseq version: ", packageVersion("phyloseq"))
suppressMessages(library(biomformat, lib.loc="/home/egonzale/dada2_test/libraries_ASVA"))
paste("biomformat version: ", packageVersion("biomformat"))

usage=function(errM) {
  cat("\nUsage : Rscript asva.R [option] <Value>\n")
  cat("       -r        : raw reads folder\n")
  cat("       -d        : design file\n")
  cat("       -o        : output directory\n")
  cat("       -tr        : database trainset\n")
  cat("       -tax        : taxonomy file\n")
  cat("       -h        : this help\n\n")
  stop(errM)
}



ARG = commandArgs(trailingOnly = T)
pipeDir=getwd()
## default arg values
rawReadsFolder=""
output_directory=""
designFile_path=""
trainset_path=""
taxonomy_path=""
ampliconLengthFile = paste(pipeDir, "/metrics/FlashLengths.tsv", sep="")


## get arg variables
for (i in 1:length(ARG)) {
  if (ARG[i] == "-r") {
    rawReadsFolder=ARG[i+1]
  } else if (ARG[i] == "-o") {
    output_directory=ARG[i+1]
  } else if (ARG[i] == "-d") {
    designFile_path=ARG[i+1]
  } else if (ARG[i] == "-tr") {
    trainset_path=ARG[i+1]
  } else if (ARG[i] == "-tax") {
    taxonomy_path=ARG[i+1]
  } else if (ARG[i] == "-h") {
    usage("")
  }
}
## check arg consitency
if (!(file.exists(designFile_path))) {
  usage("Error : Design file not found")
}
if (!(file.exists(trainset_path))) {
  usage("Error : trainset file not found")
}
if (!(file.exists(taxonomy_path))) {
  usage("Error : taxonomy file not found")
}
if (output_directory == "") {
  usage("Error : Output directory not specified")
}
if (!(file.exists(ampliconLengthFile))) {
  usage("Error : Parsed amplicon length file not found in metrics/FlashLengths.tsv")
}

#Parse ampliconLengthFile to get the minumunm and maximum lengths
ampLen = read.table(file =ampliconLengthFile, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
ampLen<- as.data.frame(ampLen)
minAmpLen=min(ampLen$`Minimum Amplicon Length`)
maxAmpLen=min(ampLen$`Maximum Amplicon Length`)

#Parse raw reads folder to get the read lengths (uses ShortRead package)

gzReadsFolder <- list.files(rawReadsFolder, pattern = "pair[12].fastq.gz$", full = TRUE, recursive = TRUE)
length(gzReadsFolder)



reads <- readFastq(gzReadsFolder)
readsLength = floor(quantile(reads@quality@quality@ranges@width,.01))

#parse sample names
sample.names <- unique(sapply(strsplit(basename(gzReadsFolder), ".pair[12].fastq"), `[`, 1))

#Create a job check folder (check files .ok will be there)
subDir = "job_checks"
dir.create(file.path(output_directory, subDir), showWarnings=FALSE)
job_checks = paste(output_directory, "/job_checks", sep="")

cat("\n---\nLoading input files:\n")
cat("Raw reads folder:", rawReadsFolder)
cat("\nDesign file:", designFile_path)
cat("\nTrainset file:", trainset_path)
cat("\nTaxonomy file:", taxonomy_path)
cat("\nOutput folder:", output_directory)
cat("\nMinimum amplicon length:",minAmpLen,"\nMaximum amplicon length:", maxAmpLen)
cat("\nRead lengths (1st quantile):",readsLength)
cat("\nNumber of samples:",length(sample.names))
cat("\nSample Names:\n", gzReadsFolder)
cat("\n---\n")


cat("\n###########################UNZIPPING READS#################################\n")



#Unzip reads into a new directory
subDir = sprintf("dada2_rawReads")
dir.create(file.path(output_directory, subDir), showWarnings=FALSE)
rawReadsFolder = paste(output_directory, "/dada2_rawReads", sep="")


checkfile=sprintf("%s/UnzipppingReads.ok",job_checks)
if(!file.exists(checkfile)){
  start_time <- Sys.time()
  cat("\nUnzipping reads\n")
  for (i in gzReadsFolder) {
    destfile<-basename(strsplit(i, ".gz")[[1]])
    if (!file.exists(destfile)) {
      cat(".")
      system(sprintf("gunzip < %s > %s/%s",i, rawReadsFolder, destfile))
    }
  }
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
  try(system(sprintf("touch %s/UnzipppingReads.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}

# Forward and reverse fastq filenames have the following format: SampleName.pair[12].fastq
fnFs <- sort(list.files(rawReadsFolder, pattern = ".pair1.fastq$", full = TRUE, recursive = TRUE))
fnRs <- sort(list.files(rawReadsFolder, pattern = ".pair2.fastq$", full = TRUE, recursive = TRUE))
length(fnFs)
length(fnRs)


cat("\n###########################QUALITY CONTROL#################################\n")

cat("\nExamine quality profiles graphs of forward and reverse reads")
subDir = sprintf("qualityProfilesRawReads")
dir.create(file.path(output_directory, subDir), showWarnings=FALSE)
savefolder = paste(output_directory, "/qualityProfilesRawReads", sep="")

checkfile=sprintf("%s/QCGraphs.ok",job_checks)
if(!file.exists(checkfile)){
  start_time <- Sys.time()
  cat("Forward reads\n")
  for (i in fnFs) {
    j<-sapply(strsplit(i, sprintf("%s/",rawReadsFolder)), `[`, -1)
    destfile=sprintf("%s/qualityProfile_%s.pdf",savefolder,j)
    if (!file.exists(destfile)) {
      cat(".")
      pdf(file=destfile, width=10, height=7)
      print(plotQualityProfile(i))
      dev.off()
    }
  }
  cat("\nReverse reads\n")
  for (i in fnRs) {
    j <- sapply(strsplit(i, sprintf("%s/",rawReadsFolder)), `[`, -1)
    destfile=sprintf("%s/qualityProfile_%s.pdf",savefolder,j)
    if (!file.exists(destfile)) {
      cat(".")
      pdf(file=destfile, width=10, height=7)
      print(plotQualityProfile(i))
      dev.off()
    }
  }
  cat("\n-----\ndone! Check pdf files in: ", savefolder,"\n", sep = "")
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  try(system(sprintf("touch %s/QCGraphs.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}


checkfile=sprintf("%s/QCStats.ok",job_checks)
if(!file.exists(checkfile)){

  cat("\nExamine quality statistics of forward and reverse reads at different quantiles")
  srqa=qa(rawReadsFolder,"fastq")
  df <- srqa[["perCycle"]]$quality
  sampleList=unique(df$lane)
  for (sample in sampleList) {
    cat("\n",sample)
    df_temp=df[df[,"lane"] == sample,]
    df_temp$lane <- NULL
    df_temp$Quality <- NULL
    means <- rowsum(df_temp$Score * df_temp$Count, df_temp$Cycle)/rowsum(df_temp$Count,df_temp$Cycle)
    means=data.frame(means)
    temp=data.frame("cycle"=rownames(means),"means"=means$means)
    write.table(temp ,file =sprintf("%s/%s_ReadsQC.txt",savefolder,sample), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  #We'll check the quality at different quantile
  cat("\nQuality of all nucleotide in forward reads at different quantiles")
  srqa=qa(rawReadsFolder,".pair1.fastq$")
  df <- srqa[["perCycle"]]$quality
  qualPerCycleFwd <- matrix(ncol=(55/5), nrow=(max(df$Cycle)+1)-min(df$Cycle))
  colnames(qualPerCycleFwd) <- colnames(qualPerCycleFwd, do.NULL = FALSE)
  for (i in min(df$Cycle):max(df$Cycle)) {
    temp=df[df$Cycle == i,]
    counter=0
    for (q in seq(from=0, to=50, by=5)) {
      if (q == 0) {q = 1}
      q=q/100
      counter=counter+1
      colnames(qualPerCycleFwd)[counter]<-paste("q", q, sep = "")
      qualPerCycleFwd[i,counter] = quantile(temp$Score,q)[[1]]
    }
  }  
  qualPerCycleFwd <- data.frame(qualPerCycleFwd)
  qualPerCycleFwd$Cycle <- rownames(qualPerCycleFwd)
  qualPerCycleFwd = qualPerCycleFwd[,c(which(colnames(qualPerCycleFwd)=="Cycle"),which(colnames(qualPerCycleFwd)!="Cycle"))]
  write.table(qualPerCycleFwd ,file =sprintf("%s/ForwardReads_QC.txt",savefolder), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n-----\ndone! Check file: ", savefolder,"/ForwardReads_QC.txt", sep = "")
  
  cat("\nQuality of all nucleotide in reverse reads at different quantiles")
  srqa=qa(rawReadsFolder,".pair2.fastq$")
  df <- srqa[["perCycle"]]$quality
  qualPerCycleRev <- matrix(ncol=(55/5), nrow=(max(df$Cycle)+1)-min(df$Cycle))
  colnames(qualPerCycleRev) <- colnames(qualPerCycleRev, do.NULL = FALSE)
  for (i in min(df$Cycle):max(df$Cycle)) {
    temp=df[df$Cycle == i,]
    counter=0
    for (q in seq(from=0, to=50, by=5)) {
      if (q == 0) {q = 1}
      q=q/100
      counter=counter+1
      colnames(qualPerCycleRev)[counter]<-paste("q", q, sep = "")
      qualPerCycleRev[i,counter] = quantile(temp$Score,q)[[1]]
    }
  }  
  qualPerCycleRev <- data.frame(qualPerCycleRev)
  qualPerCycleRev$Cycle <- rownames(qualPerCycleRev)
  qualPerCycleRev = qualPerCycleRev[,c(which(colnames(qualPerCycleRev)=="Cycle"),which(colnames(qualPerCycleRev)!="Cycle"))]
  write.table(qualPerCycleRev ,file =sprintf("%s/ReverseReads_QC.txt",savefolder), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n-----\ndone! Check file: ", savefolder,"/ReverseReads_QC.txt", sep = "")
  
  
  
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
  try(system(sprintf("touch %s/QCStats.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}




# Make directory and filenames for the filtered fastqs
start_time <- Sys.time()

filt_path <- file.path(output_directory, "Filtered_reads")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
length(filtFs)
length(filtRs)

checkfile=sprintf("%s/FilteringAndTrimming.ok",job_checks)
if(!file.exists(checkfile)){
    cat("\n###########################FILTERING AND TRIMMING#################################\n")
    cat("\nGetting the right-side average minimum quality position for forward and reverse reads\nRules: the minimum quality is set when 2 consecutive positions are below 40, 30 or 20\n")
    #Here the goal is to have a table where for each phred33 quality (from 10 to 40, increment of 5) we select at which position the forward and revertse reads should be trimmed at each lower quantile  of the data (from 1st to 50th, increment of 5)
    #for example:
    #Quantile    Quality    ForwardReads    ReverseReads
    #0.01           20           80             80
    #0.50           20          229            229
    #This table says:
    #  - To keep 99% (~0.01 quantile) of the data at a quality >=20, we should trim the forward and reverse reads at position 80
    #  - To keep 50% of the data at a quality >=20, we should trim the forward and reverse reads at position 229
    
    #dimension of the final table: 4 columns and the number of rows is set by the range 10 to 40 (incremewnt of 5) and by the number of quantile analyzed (11: 0.01, 0.05, 0.1 ... 0.5)
    qualTrimOptions <- matrix(ncol=4, nrow=(45-10)/5*11)
    colnames(qualTrimOptions) <- c("Quantile", "Quality", "ForwardReads", "ReverseReads")
    counter=0
    #Check if qualPerCycleFwd exists, if not upload it form previous run
    if (!exists("qualPerCycleRev")) {
      qualPerCycleFwdLoc = sprintf("%s/qualityProfilesRawReads/ForwardReads_QC.txt",output_directory)
      qualPerCycleFwd = read.table(file =qualPerCycleFwdLoc, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
      qualPerCycleFwd <- as.data.frame(qualPerCycleFwd)
      qualPerCycleRevLoc = sprintf("%s/qualityProfilesRawReads/ReverseReads_QC.txt",output_directory)
      qualPerCycleRev = read.table(file =qualPerCycleRevLoc, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
      qualPerCycleRev <- as.data.frame(qualPerCycleRev)
    }
    for (qua in seq(from=10, to=40, by=5)) {
      for (perc in 2:length(colnames(qualPerCycleFwd))) {
        counter = counter + 1
        fwdQC = qualPerCycleFwd[,c(1,perc)]
        quantile=gsub('q', '', colnames(fwdQC)[2])
        qualTrimOptions[counter,1]=quantile
        qualTrimOptions[counter,2]=qua
        goodQualTestFwd=dim(fwdQC)[1] + 1 - head(which(sapply((seq(from=dim(fwdQC)[[1]], to=1, by=-1)),function(i){fwdQC[,2][i]>qua && fwdQC[,2][i-1]>qua})),1)
        if (length(goodQualTestFwd)==0) {
          qualTrimOptions[counter,3]=1
        } else {
          qualTrimOptions[counter,3]=goodQualTestFwd
        }
        
        #cat("\n",quantile,"quantile.",qua,"quality. Forward reads should be trimmed from position:", goodQualTestFwd)
        revQC = qualPerCycleRev[,c(1,perc)]
        #Here the goal is to look at which position we have 2 consecutive records that are lower than the quality threshold. We start from the end of the reads so we don't stop in the beginning if the quality is lower
        goodQualTestRev=dim(revQC)[1] + 1 - head(which(sapply((seq(from=dim(revQC)[[1]], to=1, by=-1)),function(i){revQC[,2][i]>qua && revQC[,2][i-1]>qua})),1)
        if (length(goodQualTestRev)==0) {
          qualTrimOptions[counter,4]=1
        } else {
          qualTrimOptions[counter,4]=goodQualTestRev
        }
      }
    }
    qualTrimOptions = na.omit(qualTrimOptions)
    qualTrimOptions <- data.frame(qualTrimOptions)
    write.table(qualTrimOptions, file =sprintf("%s/qualityProfilesRawReads/qualTrimOptions.txt",output_directory), sep = "\t", quote = FALSE, row.names = FALSE)

    cat("\ndone!\n")
    ####################################################################################
    #We need now to take the best trim values for the forward and reverse reads. We'll do this in several steps:
    #    1. Based on the qualTrimOptions table just created, we will pick all possible forward (posFwd) and reverse (posRev) read lengths 
    #    2. We test all potential read assemblies (step 1) and keep only the ones with an overlap of at least 25bp and a length = to the longest amplicon (maxAmpLen)
    #    3. For all good read lengths (step 2), we will screen the ones where both have a quality at least >25.
    #        When multiple read length satisfy the 25 threshold, we will keep the ones with minimum quantile (largest proportion of the data) and maximum quality
    #    4. Now we have all possible pairs of reads that assemble to give at least the longest amplicon at a quality at least of 25.
    #        We then pick the lengths that will allow to have most data (lowest quantiles) with the highest quality and we will favours longer forward read lengths eventually in case of multiple hits
    #    5. #If possible, we will shorten the selected read length (step4) so a 25bp overlap will give the longest amplicon (maxAmpLen)
    
    #step 1
    posFwd = unique(qualTrimOptions$ForwardReads)
    posRev = unique(qualTrimOptions$ReverseReads)
    #step2
    dfTest <- matrix(ncol=2, nrow=length(posFwd)*length(posRev))
    colnames(dfTest) <- c("posFwd", "posRev")
    counter = 0
    for (i in posFwd) {
      for (j in posRev) {
        if ((as.numeric(as.character(i)) + as.numeric(as.character(j)) - 25 >= maxAmpLen) && (as.numeric(as.character(i)) > 200) && (as.numeric(as.character(j)) > 160)){
          counter = counter + 1
          dfTest[counter,1]=i
          dfTest[counter,2]=j
        }
      }
    }

    dfTest = na.omit(dfTest)
    dfTest <- data.frame(dfTest)
    write.table(dfTest, file =sprintf("%s/qualityProfilesRawReads/potentialTrimValues.txt",output_directory), sep = "\t", quote = FALSE, row.names = FALSE)

    #step3
    df <- matrix(ncol=6, nrow=nrow(dfTest))
    colnames(df) <- c("quantileFwd", "quantileRev","qualityFwd", "qualityRev","ForwardReads", "ReverseReads")
    counter = 0
    for (i in 1:nrow(dfTest)) {
      posFwd=as.numeric(as.character(dfTest[i,1]))
      posRev=as.numeric(as.character(dfTest[i,2]))
      indicesfwd <- which(qualTrimOptions$ForwardReads == posFwd)
      verifyFwd = qualTrimOptions[indicesfwd,]
      indicesrev <- which(qualTrimOptions$ReverseReads == posRev)
      verifyRev = qualTrimOptions[indicesrev,]
      if ((length(which(as.numeric(as.character(verifyFwd$Quality)) >= 20)) != 0) && (length(which(as.numeric(as.character(verifyRev$Quality)) >= 20)) != 0)) {
        counter = counter + 1
        verifyFwd$ReverseReads <- NULL
        #down to a single value: minimum quantile and maximum quality
        verifyFwd = verifyFwd[order(as.numeric(as.character(verifyFwd$Quantile)), -as.numeric(as.character(verifyFwd$Quality))),]
        rownames(verifyFwd) <- NULL
        verifyFwd = verifyFwd[min(which(as.numeric(as.character(verifyFwd$Quality)) >= 20)),]
        verifyRev$ForwardReads <- NULL
        #down to a single value: minimum quantile and maximum quality
        verifyRev = verifyRev[order(as.numeric(as.character(verifyRev$Quantile)), -as.numeric(as.character(verifyRev$Quality))),]
        #reset index
        rownames(verifyRev) <- NULL
        verifyRev = verifyRev[min(which(as.numeric(as.character(verifyRev$Quality)) >= 20)),]
        df[counter,1] <- as.numeric(as.character(verifyFwd$Quantile))
        df[counter,2] <- as.numeric(as.character(verifyRev$Quantile))
        df[counter,3] <- as.numeric(as.character(verifyFwd$Quality))
        df[counter,4] <- as.numeric(as.character(verifyRev$Quality))
        df[counter,5] <- as.numeric(as.character(verifyFwd$ForwardReads))
        df[counter,6] <- as.numeric(as.character(verifyRev$ReverseReads))
      }
    }
    df = na.omit(df)
    df <- data.frame(df)
    #step4
    df$quantSum = df$quantileFwd + df$quantileRev
    df$qualSum = df$qualityFwd + df$qualityRev
    df = df[order(as.numeric(as.character(df$quantSum)), -as.numeric(as.character(df$qualSum)), -as.numeric(as.character(df$ForwardReads))),]
    write.table(df, file =sprintf("%s/qualityProfilesRawReads/bestTrimValues.txt",output_directory), sep = "\t", quote = FALSE, row.names = FALSE)
    #we will pick top value by iterating (so we don't have small values of truncLenFwd or truncLenRev)
    check = 0
    for (i in 1:nrow(df)) {
      trimFwd=as.numeric(as.character(df[i,"ForwardReads"]))
      trimRev=as.numeric(as.character(df[i,"ReverseReads"]))
      if (check != 1) {
        truncLenFwd = trimFwd
        truncLenRev = trimRev
        quantFwd = as.numeric(as.character(df[i,"quantileFwd"]))
        quantRev = as.numeric(as.character(df[i,"quantileRev"]))
        check=1
      }
    }

    Summary_path <- file.path(output_directory, "Summary")
    if(!file_test("-d", Summary_path)) dir.create(Summary_path)
    
    cat("\n---\nFilter and trim parameters:\n",file=sprintf("%s/FilterAndTrimming.txt",Summary_path),sep="\n")
    cat("\nMinimum amplicon length:",minAmpLen,"\nMaximum amplicon length:", maxAmpLen,"\nForward reads trimming up to position:", truncLenFwd, "\nReverse reads trimming up to position:", truncLenRev,"\nQuantile of total reads kept (lower the better):", quantFwd,"(forward)", quantRev, "(reverse)\n",file=sprintf("%s/FilterAndTrimming.txt",Summary_path),append=TRUE)
    cat("\nStarting filtering process\n---\n")
    
    #trunclen parameter: If a length 1 vector is provided, the same parameter value is used for the forward and reverse reads. 
    #                    If a length 2 vector is provided, the first value is used for the forward reads, and the second for the reverse reads.
    #truncLen must be large enough to maintain at least 25 nucleotides of overlap between them.
    cat("\nfnFs:", length(fnFs), "\n")
    cat("filtFs:", length(filtFs), "\n")
    cat("fnRs:", length(fnRs), "\n")
    cat("filtRs:", length(filtRs), "\n")
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(truncLenFwd,truncLenRev), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,, compress=FALSE, multithread=6) 

    cat("\nCreating a summary table of filterAndTrim process")
    Summary_path <- file.path(output_directory, "Summary")
    if(!file_test("-d", Summary_path)) dir.create(Summary_path)
    track <- out
    # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
    colnames(track) <- c("input", "filtered")
    rownames(track) <- sample.names
    write.table(track ,file =sprintf("%s/Summary1.txt",Summary_path), sep = "\t", quote = FALSE)
    try(system(sprintf("touch %s/Summary1.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))

    try(system(sprintf("touch %s/FilteringAndTrimming.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
    end_time <- Sys.time()
    cat("\nLast process time:",end_time - start_time,"secs")
    cat("\n---\ndone!\n")
}




errorLearning_path <- file.path(output_directory, "Error_Rate_Learning")
if(!file_test("-d", errorLearning_path)) dir.create(errorLearning_path)

checkfile=sprintf("%s/ErrorLearning.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n###########################ERROR RATES LEARNING#################################\n")
  cat("\nError rates are learned by alternating between sample inference and error rate estimation until convergence.\n")
  errF <- learnErrors(filtFs, verbose = FALSE, MAX_CONSIST = 50, randomize = TRUE, multithread = 6)
  errR <- learnErrors(filtRs, verbose = FALSE, MAX_CONSIST = 50, randomize = TRUE, multithread = 6)
  cat("\ndone!\n")
  
  cat("Creating error rates graphs\n")
  
  
  pdf(file=sprintf("%s/ErrorLearning_ForwardReads.pdf",errorLearning_path), width=10, height=7)
  print(plotErrors(errF, nominalQ=TRUE))
  dev.off()
  
  pdf(file=sprintf("%s/ErrorLearning_ReverseReads.pdf",errorLearning_path), width=10, height=7)
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
  cat("\n-----\ndone! Check pdf files: ",errorLearning_path,"\n")
  #Saving Error Learning Matrices in case we have to rerun the process
  errorLearningMat_path <- file.path(output_directory, "ErrorLearningMatrices")
  if(!file_test("-d", errorLearningMat_path)) dir.create(errorLearningMat_path)
  #the number of matrices inside errR$err_in varies, so I'll use a loop
  errFFilestoOutput=list()
  for (i in 1:length(errF$err_in)) {
    print(i)
    errFFilestoOutput[[i]] = errF$err_in[[i]]
  }
  i=length(errF$err_in)+1
  errFFilestoOutput[[i]] = errF$err_out
  i=length(errF$err_in)+2
  errFFilestoOutput[[i]] = errF$trans

  #Same thing for errRFilestoOutput
  errRFilestoOutput=list()
  for (i in 1:length(errR$err_in)) {
    print(i)
    errRFilestoOutput[[i]] = errR$err_in[[i]]
  }
  i=length(errR$err_in)+1
  errRFilestoOutput[[i]] = errR$err_out
  i=length(errR$err_in)+2
  errRFilestoOutput[[i]] = errR$trans

  for (i in 1:length(errFFilestoOutput)) {
    output= sprintf("%s/errF_%s.txt",errorLearningMat_path,i)
    write.table(errFFilestoOutput[[i]] ,file = output, sep = "\t", quote = FALSE, row.names = TRUE)
  }
  for (i in 1:length(errRFilestoOutput)) {
    output= sprintf("%s/errR_%s.txt",errorLearningMat_path,i)
    write.table(errRFilestoOutput[[i]] ,file = output, sep = "\t", quote = FALSE, row.names = TRUE)
  }
  try(system(sprintf("touch %s/ErrorLearning.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}


checkfile=sprintf("%s/MergePairedReads.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n########################### DEREPLICATION #################################\n")
  
  start_time <- Sys.time()
  cat("\nStarting dereplication  process")
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  try(system(sprintf("touch %s/Dereplication.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
  
  
  
  cat("\n###########################Sample Inference#################################\n")
  start_time <- Sys.time()
  cat("Apply the core sequence-variant inference algorithm to the dereplicated data.\n")
  #pool = TRUE: the algorithm will pool together all samples prior to sample inference it may give better results for low sampling depths at the cost of increased computation time
  #    Note about pool: f you are dealing with datasets that approach or exceed available RAM, it is preferable to process samples one-by-one
  
  #Check if the errR and errF lists from ERROR RATES LEARNING exist (in case an analysis is re-running )
  if (!exists("errF")) {
      #We need to know how many files were output before:
      errF_path <- file.path(output_directory, "ErrorLearningMatrices")
      errF_files <- list.files(errF_path, pattern = "errF_*", full = TRUE, recursive = TRUE)
      #I'll read all files first
      MAT=list()
      for (i in 1:length(errF_files)) {
        inputfile= sprintf("%s/errF_%s.txt",errF_path,i)
        MAT[[i]] = as.matrix(read.table(inputfile, header=T, sep = "\t", row.names=1, com='', check.names=FALSE))
      }
      #now that the matrices are in a list, we'll create a new list for the matrices inside dada2's err_in
      errInMat = list() 
      for (i in 1:(length(errF_files)-2)){
        errInMat[[i]] = MAT[[i]]
      }
      #Reconstruct the err object
      errF = list(err_out = MAT[[length(errF_files)-1]], err_in = errInMat, trans = MAT[[length(errF_files)]])
  }
  
  if (!exists("errR")) {
      #We need to know how many files were output before:
      errR_path <- file.path(output_directory, "ErrorLearningMatrices")
      errR_files <- list.files(errR_path, pattern = "errR_*", full = TRUE, recursive = TRUE)
      #I'll read all files first
      MAT=list()
      for (i in 1:length(errR_files)) {
        inputfile= sprintf("%s/errR_%s.txt",errR_path,i)
        MAT[[i]] = as.matrix(read.table(inputfile, header=T, sep = "\t", row.names=1, com='', check.names=FALSE))
      }
      #now that the matrices are in a list, we'll create a new list for the matrices inside dada2's err_in
      errInMat = list() 
      for (i in 1:(length(errR_files)-2)){
        errInMat[[i]] = MAT[[i]]
      }
      #Reconstruct the err object
      errR = list(err_out = MAT[[length(errR_files)-1]], err_in = errInMat, trans = MAT[[length(errR_files)]])
  }
  
  
  
  dadaFs <- dada(derepFs, err=errF, multithread=6, pool=FALSE)
  dadaRs <- dada(derepRs, err=errR, multithread=6, pool=FALSE)

  try(system(sprintf("touch %s/SampleInference.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
  
  
  cat("\n###########################MERGE PAIRED READS#################################\n")
  start_time <- Sys.time()
  cat("\nAssembling paired-end reads into ASVs")
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=20, verbose=TRUE)
  # Inspect the merger data.frame from the first sample
  #print(head(mergers[[1]]))
  cat("\ndone!\n")
  
  cat("\nConstructing sequence table")
  
  seqtab <- makeSequenceTable(mergers)
  print(dim(seqtab)[-1])
  cat("Distribution of sequence lengths:\n")
  #print(table(nchar(getSequences(seqtab))))
  

  # Make directory and filenames for sequence inference
  sequenceInference_path <- file.path(output_directory, "sequenceInference")
  if(!file_test("-d", sequenceInference_path)) dir.create(sequenceInference_path)

  write.table(seqtab, file =sprintf("%s/seqtab.txt",sequenceInference_path), sep = "\t", quote = FALSE)

  cat("\nCreating sequence inferrence histogram\n")
  pdf(file=sprintf("%s/Sequence_inference_histogram.pdf",sequenceInference_path), width=10, height=7)
  print(hist(nchar(getSequences(seqtab)),  
       main=paste("Amplicon lengths overview:", dim(seqtab)[-1],"generated amplicons in total"),
       xlab="Amplicon length",
       border="black", col="gold",
       breaks=dim(seqtab)[-1]))
  dev.off()
  #Inspecting the dada-class object returned by dada
  cat("\nCreating a summary table of Sample Inference process")
  Summary_path <- file.path(output_directory, "Summary")
  if(!file_test("-d", Summary_path)) dir.create(Summary_path)
  getN <- function(x) sum(getUniques(x))
  track <- cbind(sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab))
  colnames(track) <- c("denoised", "merged", "rawOTUs")
  rownames(track) <- sample.names
  write.table(track, file = sprintf("%s/Summary2.txt",Summary_path), sep = "\t", quote = FALSE)
  try(system(sprintf("touch %s/Summary2.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))

  try(system(sprintf("touch %s/MergePairedReads.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
}
  
  #Sequences that are much longer or shorter than expected may be the result of non-specific priming, 
  #and may be worth removing (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)]). 
  #This is analogous to “cutting a band” in-silico to get amplicons of the targeted length.
  
  
  
  
  
checkfile=sprintf("%s/RemoveChimeras.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n###########################REMOVE CHIMERAS#################################\n")
  if (!exists("seqtab")) {
    sequenceInference_path <- file.path(output_directory, "sequenceInference")
    seqtab = read.table(file = sprintf("%s/seqtab.txt",sequenceInference_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
    seqtab <- as.matrix(seqtab)
  }
  start_time <- Sys.time()
  cat("\nRemoving chimeras")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=6, verbose=TRUE)
  cat("\nFraction of chimeras:", (1-sum(seqtab.nochim)/sum(seqtab))*100,"% of the total sequence reads.\n")
  #Most of your reads should remain after chimera removal (it is not uncommon for a majority of sequence 
  #variants to be removed though). If most of your reads were removed as chimeric, upstream processing 
  #may need to be revisited. In almost all cases this is caused by primer sequences with ambiguous 
  #nucleotides that were not removed prior to beginning the DADA2 pipeline.
  write.table(seqtab.nochim, file =sprintf("%s/seqtab.nochim.txt",sequenceInference_path), sep = "\t", quote = FALSE)

  cat("\nCreating a summary table of remove chimera process")
  Summary_path <- file.path(output_directory, "Summary")
  if(!file_test("-d", Summary_path)) dir.create(Summary_path)
  track <- rowSums(seqtab.nochim)
  write.table(track ,file =sprintf("%s/Summary3.txt",Summary_path), sep = "\t", col.names="Chimera", quote = FALSE)
  try(system(sprintf("touch %s/Summary3.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))

  try(system(sprintf("touch %s/RemoveChimeras.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}


checkfile=sprintf("%s/Summary.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n###########################SUMMARY TABLES #################################\n")
  #we will gather the three summary tables into one
  Summary_path <- file.path(output_directory, "Summary")
  if(!file_test("-d", Summary_path)) dir.create(Summary_path)
  df1 = read.table(file = sprintf("%s/Summary1.txt",Summary_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
  df2 = read.table(file = sprintf("%s/Summary2.txt",Summary_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
  df3 = read.table(file = sprintf("%s/Summary3.txt",Summary_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
  Statsondf = merge(df1,df2,by="row.names",all.x=TRUE)
  rownames(Statsondf) <- Statsondf$Row.names
  Statsondf$Row.names <- NULL
  Statsondf = merge(Statsondf,df3,by="row.names",all.x=TRUE)
  rownames(Statsondf) <- Statsondf$Row.names
  Statsondf$Row.names <- NULL
  write.table(Statsondf ,file =sprintf("%s/Summary.txt",Summary_path), sep = "\t", quote = FALSE)
  system(sprintf("rm -f %s/Summary[123].txt",Summary_path), intern = FALSE, ignore.stderr = FALSE)
  try(system(sprintf("touch %s/Summary.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}





checkfile=sprintf("%s/OTUTable.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n###########################OTU TABLES #################################\n")
  if (!exists("seqtab.nochim")) {
    sequenceInference_path <- file.path(output_directory, "sequenceInference")
    seqtab.nochim = read.table(file = sprintf("%s/seqtab.nochim.txt",sequenceInference_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
    seqtab.nochim <- as.matrix(seqtab.nochim)
  }

  cat("\nCreating OTU table in 2 different format (txt and biom)\n")
  dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
  if(!file_test("-d", dada2outputFiles_path)) dir.create(dada2outputFiles_path)
  write.table(t(seqtab.nochim), file =sprintf("%s/otu_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)
  biomFile <- make_biom(t(seqtab.nochim))
  write_biom(biomFile, sprintf("%s/otu_table.biom",dada2outputFiles_path))
  

  try(system(sprintf("touch %s/OTUTable.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
}


###################################TAXONOMY ASSIGNMENT##############################################

checkfile=sprintf("%s/TaxAssignment.ok",job_checks)
if(!file.exists(checkfile)){
  cat("\n###########################TAXONOMY ASSIGNMENT#################################\n")
  start_time <- Sys.time()
  
  
  if (!exists("seqtab.nochim")) {
    sequenceInference_path <- file.path(output_directory, "sequenceInference")
    seqtab.nochim = read.table(file = sprintf("%s/seqtab.nochim.txt",sequenceInference_path), header = T, sep = "\t", row.names=1, com='', check.names=FALSE)
    seqtab.nochim <- as.matrix(seqtab.nochim)
  }
  cat("\nReading taxonomy databases")
  taxa <- assignTaxonomy(seqtab.nochim, trainset_path)
  taxa <- addSpecies(taxa, taxonomy_path)
  cat("\ndone!\n")
  
  cat("\nWriting taxonomy file")
  dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
  if(!file_test("-d", dada2outputFiles_path)) dir.create(dada2outputFiles_path)
  write.table(taxa ,file =sprintf("%s/tax_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)
  try(system(sprintf("touch %s/TaxAssignment.ok",job_checks), intern = FALSE, ignore.stderr = FALSE))
  end_time <- Sys.time()
  cat("\nLast process time:",end_time - start_time,"secs")
  cat("\n---\ndone!\n")
}


cat("\n##################################### MICROBIOME ANALYST #################################\n")
start_time <- Sys.time()
MA_path <- file.path(output_directory, "MicrobiomeAnalyst")
if(!file_test("-d", MA_path)) dir.create(MA_path)

cat("\nLoading Sample data (design)\n")
presampleTable = read.table(file =designFile_path, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
presampleTable<- as.data.frame(presampleTable)


cat("\nLoading OTU table\n")
dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
myOtuTable_path = sprintf("%s/otu_table.txt",dada2outputFiles_path)
preOtuTable <- data.matrix(read.table(file =myOtuTable_path, header=T, row.names=1, com='', stringsAsFactors=FALSE, check.names=FALSE)) #check.names=FALSE is used to force R not to interpret numbers as numbers (as my Sample id are numbers)


cat("\nLoading taxonomy table\n")
mytaxtable_path = sprintf("%s/tax_table.txt",dada2outputFiles_path)
preTaxTable <- as.matrix(read.table(file =mytaxtable_path, header=T, row.names=1, com='', sep="\t"))


#create MA otu table
MA_otu = as.data.frame(preOtuTable)
#add an otu column
MA_otu[,"#NAME"] <- paste0("OTU_",seq.int(nrow(MA_otu),sep=""))
#move last column to first
MA_otu = MA_otu[,c(which(colnames(MA_otu)=="#NAME"),which(colnames(MA_otu)!="#NAME"))]
#output this table as OTU table
write.table(MA_otu ,file =sprintf("%s/OTU_table.txt",MA_path), sep = "\t", quote = FALSE, row.names = FALSE)

#create MA sample table
MA_sample = presampleTable
#add # name column which is the row.names
MA_sample[,"#NAME"] = row.names(presampleTable)
#move last column to first
MA_sample = MA_sample[,c(which(colnames(MA_sample)=="#NAME"),which(colnames(MA_sample)!="#NAME"))]
#output this table as OTU table
write.table(MA_sample ,file =sprintf("%s/Metadata_file.txt",MA_path), sep = "\t", quote = FALSE, row.names = FALSE)

#create MA taxonomy table
MA_tax = preTaxTable
#create a crib between 
crib=data.frame("OTU"=MA_otu[,"#NAME"], "sequences" = row.names(MA_otu))
#merge crib with tax table
MA_tax = merge(x=crib,y=MA_tax, by.x ="sequences", by.y ="row.names", all.x = TRUE, all.y = TRUE)
#remove sequences
MA_tax$sequences <- NULL
#here I will join genera and species when species is not NA (NS for New Species)
MA_tax$NS <- paste(MA_tax$Genus,"_",MA_tax$Species, sep = "")
MA_tax$NS <- gsub('.*_NA.*', 'NA', MA_tax$NS)
MA_tax$Species = MA_tax$NS
MA_tax$NS <- NULL

#rename OTU col
names(MA_tax)[names(MA_tax)=="OTU"] <- "#TAXONOMY"
#output this table as Taxonomy table
write.table(MA_tax ,file =sprintf("%s/Taxonomy_table .txt",MA_path), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\ndone!")




cat("\n##################################### Excel-ready #################################\n")
excel_path <- file.path(output_directory, "ExcelReadyOutput")
if(!file_test("-d", excel_path)) dir.create(excel_path)
#create the otu table
countTable = as.data.frame(preOtuTable)
#add an ASV_name column (handling sequence as unique identifier can be hard)
countTable[,"ASV_name"] <- paste0("ASV_",seq.int(nrow(countTable),sep=""))
#create the taxonomy table
taxTable = preTaxTable
excelTable=merge(taxTable,countTable, by="row.names", all.x = FALSE, all.y = TRUE)
names(excelTable)[names(excelTable)=="Row.names"] <- "ASV_sequences"
#here I will join genera and species when species is not NA (NS for New Species)
excelTable$NS <- paste(excelTable$Genus,"_",excelTable$Species, sep = "")
excelTable$NS <- gsub('.*_NA.*', 'NA', excelTable$NS)
excelTable$Species = excelTable$NS
excelTable$NS <- NULL
#move ASV_name column to first column
excelTable = excelTable[,c(which(colnames(excelTable)=="ASV_name"),which(colnames(excelTable)!="ASV_name"))]
write.table(excelTable ,file =sprintf("%s/ASV_table .txt",excel_path), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\ndone!")



cat("\n##################################### SOME GRAPHS#################################\n")

start_time <- Sys.time()
graphs_path <- file.path(output_directory, "graphics")
if(!file_test("-d", graphs_path)) dir.create(graphs_path)


cat("\nUsing Phyloseq package to create some graphs about the computed data")
sampleTable = sample_data(presampleTable)
OtuTable = otu_table(preOtuTable, taxa_are_rows=TRUE)

TaxTable = tax_table(preTaxTable)


cat("\nLoading a phyloseq object")
# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(OtuTable, TaxTable, sampleTable)
ps
cat("\ndone!\n")

# Ouput the design in the same folder as otu and taxonomy table
dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
if(!file_test("-d", dada2outputFiles_path)) dir.create(dada2outputFiles_path)
write.table(sample_data(ps) ,file =sprintf("%s/sample_data.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)


#temporary fix a bug in phyloseq when you only have one condition: see https://github.com/joey711/phyloseq/issues/541
sample_data(ps)[ , 2] <- sample_data(ps)[ ,1]


cat("\nAlpha diversity plots")
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
listTaxLevels=c("Phylum", "Class", "Order", "Family", "Genus")
destfile=sprintf("%s/barCharts_top20.pdf",graphs_path)
pdf(file=destfile, width=12, height=9)
for(l in listTaxLevels) {
  p = plot_bar(ps.top20, x=sprintf("%s",l), fill=sprintf("%s",l), facet_grid=~Condition, title=sprintf("%s level of the 20 most abundant ASVs", l))
  q= p + geom_bar(aes_string(color=sprintf("%s",l), fill=sprintf("%s",l)), stat="identity", position="stack")
  print(q)
}
dev.off()
cat("\ndone!\n")


cat("\nRichness plots")
destfile=sprintf("%s/richness.pdf",graphs_path)
pdf(file=destfile, width=12, height=12)
plot_richness(ps, x = "Condition", color = "Condition", measures=c("Observed", "Shannon", "Simpson", "Fisher", "Chao1", "InvSimpson"), nrow=3) + 
  geom_boxplot() + geom_point(size = 0.5, alpha = 0.5) +
  ggtitle("Richness plot") +
  theme_classic(base_size = 16)
dev.off()
end_time <- Sys.time()
cat("\nLast process time:",end_time - start_time,"secs")
cat("\n---\ndone!\n")




cat("\nOrdination plots")
destfile=sprintf("%s/ordination.pdf",graphs_path)
pdf(file=destfile, width=10, height=7)
print("PCoA")
ps.pcoa <- ordinate(ps, method="PCoA", distance="bray")
p=plot_ordination(ps, ps.pcoa, type="samples", color="Condition")
df=p$data
ggplot(df, aes(df[,1], df[,2], color = Condition)) + ggtitle("PCoA plot") + xlab(colnames(df[1])) + ylab(colnames(df[2])) +
  geom_point(size = 3, alpha = 0.8) +geom_text(aes(label=row.names(df)),hjust=0, vjust=0, size = 3, color = 'grey') +  theme_classic(base_size = 16)
print("CCA")
ps.cca <- ordinate(ps, method="CCA", distance="bray")
p=plot_ordination(ps, ps.cca, type="Sample", color="Condition")+
  ggtitle("CCA")
df=p$data
ggplot(df, aes(df[,1], df[,2], color = Condition)) + 
  ggtitle("CCA plot") + xlab(colnames(df[1])) + ylab(colnames(df[2])) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label=row.names(df)),hjust=0, vjust=0, size = 3, color = 'grey') +
  theme_classic(base_size = 16)
cat("\ndone!\n")


cat("\nASVA is exiting normally\n Well done dude!\n")
