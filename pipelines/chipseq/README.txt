# This is an example of experiment using the chipseq pipeline
### Set your working directory
#DIR=/path/to/my/working/dir/


#########################################################################################################################
### 1- Copy Perl wrapper in your working dir. from .../mugqic_pipeline/pipeline/chipseq/chipSeq.pl
###    Copy Perl mugqic_pipeline/lib folder to your working dir.
###    Copy chipSeq.mammouth.ini file to your working dir and edit it to enter proper values.


cp ~/repos/mugqic_pipeline/pipelines/chipseq/chipSeq.pl bin
cp ~/repos/mugqic_pipeline/pipelines/chipseq/chipSeq.mammouth.ini bin
mkdir -p bin/../../lib
cp ~/repos/mugqic_pipeline/lib/* bin/../../lib/

#########################################################################################################################
### 2- Prepare NANUQ sample sheet. Here are 2 examples of how to run the sampleSetup.pl script.

module load mugqic/tools/1.0 && perl sampleSetup.pl --projectId 8968  --nanuqAuthFile ~/nanuq_sample_setup.txt   
module load mugqic/tools/1.0 && perl sampleSetup.pl --nanuqAuthFile ~/nanuq_sample_setup.txt --sampleSheet ./project.nanuq.csv 

#########################################################################################################################
### 3- Prepare the design file, either one column by pair (background, CHIP) or a column by group of chip samples to be 
### compared to one background sample (background, CHIP1, CHIP2, CHIP3, etc) or a column with a Chip sample withoutbackground


#########################################################################################################################
### 4- Generate either all commands on one file:
# perl bin/chipSeq.pl -c bin/chipSeq.mammouth.ini -n project.nanuq.csv -d design.csv -w `pwd` -s 2 -e 11 > toRun.sh

## Or run each step separately:

#for s in `seq 1 11`;
#do
#  perl bin/chipSeq.pl -c bin/chipSeq.mammouth.ini -n project.nanuq.csv -d design.csv -w `pwd` -s $s -e $s > commands_step$s.sh
#  chmod u+rwx commands_step$s.sh
#done;

#########################################################################################################################
## 5- Then source these .sh files

## either all commands

#sh toRun.sh

## or step by step
sh commands_step1.sh
# check that all steps were completed with exit code 0 and continue
sh commands_step2.sh
# and so on

