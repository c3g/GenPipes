#!/bin/bash

# This script was written by Louis Letourneau and Maxime Caron
# Processes a hiseq run


if [ $# != 5 ]; then
echo "usage: $0 RUN_LABEL RUN_ID RUN_FOLDER INDEXED INDEX_STRING"
exit
fi

set -o
set -e
RUN_LABEL=$1
RUN_ID=$( echo $2 | sed 's/\/$//' )
RUN_FOLDER=$( echo $3 |sed 's/\/$//' )
INDEXED=$4
INDEX_STRING=$5

echo "Values:"
echo $RUN_LABEL;
echo $RUN_ID;
echo $RUN_FOLDER;
echo $INDEXED;
echo $INDEX_STRING;

START=`date`
END=0

echo "$RUN_LABEL in $RUN_FOLDER"

send_email() {
  /data/solexa/tools/smtp-cli-2.8.pl --verbose --host=smtp.mcgill.ca --port=587 --creds=/data/analysis/bioinformatics/misc/smtp.bioinformatics.credentials.txt --from bioinformatics.genome@mail.mcgill.ca --to louis.letourneau@mail.mcgill.ca --to etienne.sdicu@mail.mcgill.ca --to maxime.caron2@mail.mcgill.ca --to mathieu.bourgey@mail.mcgill.ca --to etienne.sdicu@mail.mcgill.ca --subject "$RUN_LABEL done $? start: $START end: $END"
}

trap send_email INT TERM EXIT

if [ $INDEXED == 0 ]; then
  echo "No Indexes, *NOT* Generating indexes `date`"
else
  echo "Generating indexes `date`"
  #CONVERTED_MASK=`echo "$INDEX_STRING" | sed 's/Y\([0-9]\+\)/\1T/g' | sed 's/y\([0-9]\+\)/\1T/g' | sed 's/I\([0-9]\+\)/\1B/g' | sed 's/n\([0-9]\+\)/\1S/g' | sed 's/,//g'`
  OIFS=$IFS
  IFS=','
  CONVERTED_MASK=""
  for x in $INDEX_STRING
  do
      if [[ $x == [yY]* ]] ; then CONVERTED_MASK="${CONVERTED_MASK}`echo \"$x\" | sed 's/[yY]\([0-9]\+\)/\1T/g'`" ; fi
      if [[ $x == n*I*n* ]]
      then
        toEval=`echo "$x" | sed 's/n\([0-9]\+\)I\([0-9]\+\)n\([0-9]\+\)/\1 \+ \2 \+ \3/g'`
        CONVERTED_MASK="${CONVERTED_MASK}$(($toEval))B"
      else
        if [[ $x == n*I* ]] || [[ $x == I*n* ]]
        then
          toEval=`echo "$x" | sed 's/[In]\([0-9]\+\)[In]\([0-9]\+\)/\1 + \2/g'`
          CONVERTED_MASK="${CONVERTED_MASK}$(($toEval))B"
        else
          if [[ $x == I* ]]
          then
            CONVERTED_MASK="${CONVERTED_MASK}`echo \"$x\" | sed 's/[I]\([0-9]\+\)/\1B/g'`";
          fi
        fi
      fi
  done
  IFS=$OIFS
  echo "New mask: $CONVERTED_MASK"
  for i in 1 2 3 4 5 6 7 8
  do
    echo "java -Xmx5G -jar /data/solexa/tools/CountIlluminaBarcodes-0.2-jar-with-dependencies.jar BASECALLS_DIR=${RUN_FOLDER}/Data/Intensities/BaseCalls LANE=${i} READ_STRUCTURE=${CONVERTED_MASK} METRICS_FILE=${RUN_FOLDER}/${RUN_ID}_${i}.metrics MAX_MISMATCHES=1 NUM_PROCESSORS=\$NSLOTS BARCODE_FILE=/data/solexa/tools/barcodes.txt" | qsub -j y -N idxs_${i}_${RUN_LABEL} -pe multicore 20,10 -cwd -v PATH -p 512
  done
fi

echo "Generate Unaligned directory `date`"
if [ $INDEXED == 0 ]; then
  echo "configureBclToFastq.pl --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls --sample-sheet ${RUN_FOLDER}/nanuqSampleSheet.csv --fastq-cluster-count 0"
  configureBclToFastq.pl --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls --sample-sheet ${RUN_FOLDER}/nanuqSampleSheet.csv --fastq-cluster-count 0
else
  echo "configureBclToFastq.pl --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls --sample-sheet ${RUN_FOLDER}/nanuqSampleSheet.csv --fastq-cluster-count 0 --mismatches 1 --use-bases-mask ${INDEX_STRING}"
  configureBclToFastq.pl --input-dir ${RUN_FOLDER}/Data/Intensities/BaseCalls --sample-sheet ${RUN_FOLDER}/nanuqSampleSheet.csv --fastq-cluster-count 0 --mismatches 1 --use-bases-mask ${INDEX_STRING}
fi

echo "Generate fastq.gz `date`"
cd ${RUN_FOLDER}/Unaligned
echo "make -j \$NSLOTS" | qsub -j y -N fastq${RUN_LABEL} -pe multicore 20,15 -cwd -sync yes -v PATH -p 1024

echo "Run qc graphs `date`"
cd ${RUN_FOLDER};
IS_PAIRED=`ls -l Unaligned/Project_nanuq/*/*_R2_001.fastq.gz 2> /dev/null | wc -l`

JOB_IDS="";
if [ ${IS_PAIRED} == 0 ]; then
  # Single
  echo "Is single: ${IS_PAIRED}";
  for i in Unaligned/Project_nanuq/*/*_R1_001.fastq.gz
  do
    SNAME=`echo $i | sed 's/.*\/\([^/]\+\)_R1_001.fastq.gz/\1/g'`;
    DIR=`echo $i | sed 's/\(.*\)\/[^/]*_001.fastq.gz/\1/g'`;
    JOB_IDS=${JOB_IDS},`echo "mkdir -p $DIR/qc;time java -Xmx24g -Djava.awt.headless=true -jar /data/solexa/tools/mps-tools.jar -i ${DIR}/${SNAME}_R1_001.fastq.gz -Q33 -t FASTQ -o ${DIR}/qc/ -n 10 -outputIlmnHandler -regionName ${SNAME}" | qsub -cwd -v PATH -j y -N qc${RUN_LABEL} -pe multicore 10 -p 1023 | sed 's/^Your job \([0-9]\+\) .*/\1/g'`;
  done;
else
#Paired
  echo "Is paired: ${IS_PAIRED}";
  for i in Unaligned/Project_nanuq/*/*_R1_001.fastq.gz
  do
    SNAME=`echo $i | sed 's/.*\/\([^/]\+\)_R1_001.fastq.gz/\1/g'`
    DIR=`echo $i | sed 's/\(.*\)\/[^/]*_001.fastq.gz/\1/g'`
    JOB_IDS=${JOB_IDS},`echo "mkdir -p $DIR/qc;time java -Xmx24g -Djava.awt.headless=true -jar /data/solexa/tools/mps-tools.jar -i ${DIR}/${SNAME}_R1_001.fastq.gz -i2 ${DIR}/${SNAME}_R2_001.fastq.gz -Q33 -t FASTQ -o ${DIR}/qc/ -n 10 -outputIlmnHandler -regionName ${SNAME}" | qsub -cwd -v PATH -j y -N qc${RUN_LABEL} -pe multicore 10 -p 1023 | sed 's/^Your job \([0-9]\+\) .*/\1/g'`;
  done;
fi

ls ${RUN_FOLDER}/Unaligned/Project_nanuq/*/*.gz | sed 's/\(.*\/\([^\/]\+L00.\).*gz\)/\1,\2/g' > $RUN_FOLDER/blast_sampleSheet.txt
ls ${RUN_FOLDER}/Unaligned/Undetermined_indices/*/*.gz | sed 's/\(.*\/\([^\/]\+L00.\).*gz\)/\1,\2/g' >> $RUN_FOLDER/blast_sampleSheet.txt

echo "/data/solexa/tools/bioinformatics-svn/tools/testSample.sh $RUN_FOLDER/blast_sampleSheet.txt $RUN_FOLDER/Unaligned/Blast_sample ${IS_PAIRED}" | qsub -N $RUN_LABEL.Blast -cwd -V -j y -p 612

# TODO: Missing error checking. If a qc job fails, we'll never know...but it's not so bad we'll see in nanuq that stats are missing.
echo "echo \"Finished running\"" | qsub -cwd -j y -N qcWait${RUN_LABEL} -hold_jid $JOB_IDS,$RUN_LABEL.Blast -p 1023 -sync yes

echo "rsync `date`"
rsync -avP --include "**/*onfig*" --exclude "Unaligned/**.old" --include "Unaligned/**" --exclude "Thumbnail_Images/" --exclude "Data/Intensities/B*/*" --include "Data/Intensities/B*/" --exclude "Data/Intensities/*" $RUN_FOLDER /data/newrobot/hiSeqSequencer/hiSeqRuns-drop/

echo "setfacl `date`"
setfacl -R -m mask:rwx /data/newrobot/hiSeqSequencer/hiSeqRuns-drop/$RUN_ID

END=`date`
trap - INT TERM EXIT

send_email
echo "Done ${END}"
