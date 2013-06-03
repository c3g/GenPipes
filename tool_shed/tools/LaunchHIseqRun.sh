#! /bin/bash

## usage: LaunchHIseqRun.sh is.indexed(0: no 1:yes =>$2) indexStringValue($3)
## to do usage

#Index string dictionnary
Istring[0]="null"
Istring[1]="Y100,I6,y100"
Istring[2]="Y100,I6n1,y100"
Istring[3]="Y100,n2I3n1,y100"
Istring[4]="Y50,I6"
Istring[5]="Y50,I6n1"
Istring[6]="Y100,I1n5,y100"
Istring[7]="Y100,I1n6,y100"
Istring[8]="Y50,I6,y50"
Istring[9]="Y50,I6n1,y50"
Istring[10]="Y100,I4n2,y100"
Istring[11]="Y150,I6,y150"
Istring[12]="Y100,I6"
Istring[13]="Y50,n1I5,y50"


#export necessary env variables
export fol=$(pwd)
export PATH=/opt/runTools/OLB/bin/:/opt/runTools/CASAVA/bin/:/data/solexa/tools/bioinformatics-svn/tools/:$PATH
export RUN_ID=$( echo $fol | cut -d\/ -f 4)
RUN_NUM=$( echo $RUN_ID | cut -d_ -f 5 |  cut -dH -f 1 )
export RUN_LABEL="run$RUN_NUM"
export RUN_FOLDER=$(echo "$fol")
export INDEXED=$(echo "$1")
export INDEX_STRING=${Istring[$2]}

#Generate adequate SampleSheet
cd $RUN_FOLDER
#if [ $INDEXED == "0" ]
#then
#mv SampleSheet.csv nanuqSampleSheet.csv
#export INDEX_STRING=${Istring[0]}
#fi

#if [ $INDEXED == "1" ]
#then
#if [ $3 == "2" ]
#then
#cat SampleSheet.csv | sed 's/[ATGC][ATGC],,N/,,N/g' | sed 's/nanuq,[ATGC][ATGC]/nanuq,/g' > nanuqSampleSheet.csv
#elif [ $3 == "7" ]
#then
#cat SampleSheet.csv | sed 's/[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC],,N/,,N/g'  > nanuqSampleSheet.csv
#else
#cat SampleSheet.csv | sed 's/[ATGC],,N/,,N/g' > nanuqSampleSheet.csv
#fi
#fi

cat SampleSheet.csv > nanuqSampleSheet.csv

#remove the Unaligned folder otherwise casava will stop
#if [ -d $RUN_FOLDER/Unaligned ]
#then
#rm -rf $RUN_FOLDER/Unaligned
#fi

#run the processing
echo "process_hiseq_run.sh $RUN_LABEL $RUN_ID $RUN_FOLDER $INDEXED $INDEX_STRING" | qsub -N $RUN_LABEL -cwd -V -j y -p 612
