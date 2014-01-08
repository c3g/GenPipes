MUGQIC - PacBio Assembly pipeline README.

#########################################################################################################################
### 1- Copy Perl wrapper in your working dir. from .../mugqic_pipeline/pipeline/pacBioAssembly/pacBioAssembly.pl
###    Copy Perl mugqic_pipeline/lib folder to your working dir.
###    Copy pacBioAssembly.abacus.ini file to your working dir and edit it to enter proper values.      
###    Don't forget to update outdir path in the .ini file! Then generate commands for the pipeline
###    Copy the sampleSetupPacBio.pl from .../mugqic_tools/perl-tools/


#########################################################################################################################
### 2- Prepare PacBio sample sheet. Here are 3 examples of how to run the sampleSetupPacBio.pl script.

#./sampleSetupPacBio.pl --runName Run037  --nanuqAuthFile /path/to/nanuq_sample_setup.txt > ./samples.txt 
#./sampleSetupPacBio.pl --projectId 8968  --nanuqAuthFile /path/to/nanuq_sample_setup.txt > ./samples.txt  
#./sampleSetupPacBio.pl --nanuqAuthFile ~/nanuq_sample_setup.txt --sampleSheet ./project.nanuq.8968.csv > ./samples.txt

## Then manually add the estimated genome size in the samples.txt file (create a last columns with a \t and enter an integer.)
## Here (in the ./samples.txt file), you can also edit the file to just keep the samples you want to process. For instance, in 
## project 8968, there might be 4 samples associated with the project and you could want to only process sample #1 or sample #1 and #2.

#########################################################################################################################
### 2- Either all commands on one file:
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 1 -e 8 -n ./samples.txt > commands_all.sh

## Or run each steps separately:

#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 1 -e 1  -n ./samples.txt > commands_step1.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 2 -e 2  -n ./samples.txt > commands_step2.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 3 -e 3  -n ./samples.txt > commands_step3.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 4 -e 4  -n ./samples.txt > commands_step4.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 5 -e 5  -n ./samples.txt > commands_step5.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 6 -e 6  -n ./samples.txt > commands_step6.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 7 -e 7  -n ./samples.txt > commands_step7.sh
#./pacBioAssembly.pl -c  ./pacBioAssembly.abacus.ini -s 8 -e 8  -n ./samples.txt > commands_step8.sh

#########################################################################################################################
## 3- Then source these .sh files. Pray the almighty that all jobs run with no error.

