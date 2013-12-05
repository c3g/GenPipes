### Set your working dir.
#DIR=/path/to/my/working/dir/

#########################################################################################################################
### 1- Copy Perl wrapper in your working dir. from .../mugqic_pipeline/pipeline/pacBioAssembly/pacBioAssembly.pl
###    Copy Perl mugqic_pipeline/lib folder to your working dir.
###    Copy pacBioAssembly.abacus.ini file to your working dir and edit it to enter proper values.      
###    Don't forget to update outdir path in the .ini file! Then generate commands for the pipeline
###    Copy the sampleSetupPacBio.pl from .../mugqic_tools/perl-tools/


#########################################################################################################################
### 2- Prepare PacBio sample sheet. Here are 3 examples of how to run the sampleSetupPacBio.pl script.

#~/build/mugqic_tools/perl-tools/sampleSetupPacBio.pl --runName Run037  --nanuqAuthFile ~/nanuq_sample_setup.txt > ./samples.txt 
#~/build/mugqic_tools/perl-tools/sampleSetupPacBio.pl --projectId 8968  --nanuqAuthFile ~/nanuq_sample_setup.txt > ./samples.txt  
#~/build/mugqic_tools/perl-tools/sampleSetupPacBio.pl --nanuqAuthFile ~/nanuq_sample_setup.txt --sampleSheet ./project.nanuq.8968.csv > ./samples.txt

## Then manually add the estimated genome size in the samples.txt file (create a last columns with a \t and enter an integer.) 


#########################################################################################################################
### 2- Either all commands on one file:
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 1 -e 8  -w $DIR/ -n ./samples.txt > commands_all.sh

## Or run each steps separately:

#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 1 -e 1 -w $DIR/ -n ./samples.txt > commands_step1.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 2 -e 2 -w $DIR/ -n ./samples.txt > commands_step2.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 3 -e 3 -w $DIR/ -n ./samples.txt > commands_step3.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 4 -e 4 -w $DIR/ -n ./samples.txt > commands_step4.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 5 -e 5 -w $DIR/ -n ./samples.txt > commands_step5.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 6 -e 6 -w $DIR/ -n ./samples.txt > commands_step6.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 7 -e 7 -w $DIR/ -n ./samples.txt > commands_step7.sh
#~/build/mugqic_pipeline/pipelines/pacBioAssembly/pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 8 -e 8 -w $DIR/ -n ./samples.txt > commands_step8.sh

#########################################################################################################################
## 3- Then source these .sh files. Pray the almighty that all jobs run with no error.

