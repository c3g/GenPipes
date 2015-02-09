#!/bin/bash

#TODO: Not sure if absolutely need to install PDF::API2 in mirdeep2 folder instead of any other $PERL5LIB location.
#      also not sure if will work with Vienna1.8,3 and newer versions of randfold,squid,PDF::A2...

## mirDeep2 is a real pain to install.
## Anyway, it depends on bowtie1 (module),  ViennaRNA (decided to make a separate module), and
## on squid and randfold. Those two last are installed within this module. It also needs a specific perl module, installed to
## $MUGIQC/software/perl5libs/

#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#
# 
# cd ~ ; rm -rf /mnt/lustre03/bourque/bourque_group/opt/software/mirdeep2


# rm -rf /gs/project/mugqic/analyste_dev/phase2/software/mirdeep2
export PERL5LIBPREFIX=$MUGQIC_INSTALL_HOME_DEV/software/perl5libs
export PERL5LIB=$PERL5LIBPREFIX/share/perl5/:$PERL5LIBPREFIX/lib64/perl5/

SOFTWARE=mirdeep2 
VERSION="2_0_0_7"  
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download
wget "https://www.mdc-berlin.de/43969303/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/mirdeep$VERSION.zip"
unzip mirdeep$VERSION.zip
rm -rf mirdeep$VERSION.zip
cd  mirdeep$VERSION

# Download two obscure dependencies: SQUID AND RANFOLD
wget http://selab.janelia.org/software/squid/squid.tar.gz
tar -xvf squid.tar.gz ; rm squid.tar.gz
ln -s squid* squid
wget http://bioinformatics.psb.ugent.be/supplementary_data/erbon/nov2003/downloads/randfold-2.0.tar.gz
tar -xvf randfold-2.0.tar.gz  ; rm randfold-2.0.tar.gz
ln -s randfold* randfold
cd squid 
./configure
make
cd ../randfold

# MANUAL FIX: replace "INCLUDE=-I." by "$MUGQIC_INSTALL_HOME_DEV/software/mirdeep2/mirdeep2-2_0_0_5/squid"
echo $MUGQIC_INSTALL_HOME_DEV/software/mirdeep2/mirdeep$VERSION/squid/ #
# /gs/project/mugqic/analyste_dev/phase2/software/mirdeep2/mirdeep2_0_0_7/squid/
# INCLUDE=-I. -I/gs/project/mugqic/analyste_dev/phase2/software/mirdeep2/mirdeep2_0_0_7/squid/ -L/gs/project/mugqic/analyste_dev/phase2/software/mirdeep2/mirdeep2_0_0_7/squid/
nano Makefile
# INCLUDE=-I. -I/mnt/lustre03/bourque/bourque_group/opt/software/mirdeep2/mirdeep2-2_0_0_5/squid/ -L/mnt/lustre03/bourque/bourque_group/opt/software/mirdeep2/mirdeep2-2_0_0_5/squid/
# sed -e s/INCLUDE=-I./$MUGQIC_INSTALL_HOME_DEV/g temp
# 2.c.ii.3)  edit Makefile e.g. emacs Makefile
#            change line with INCLUDE=-I. to           
#            INCLUDE=-I. -Iyour_path_to_squid-1.9g -Lyour_path_to_squid-1.9g
#       e.g. INCLUDE=-I. -I/home/Pattern/squid-1.9g/ -L/home/Pattern/squid-1.9g/ 
make
cd ..

# PDF-API2 perl module (cpan not working for some reason)


wget http://cpan.metacpan.org/authors/id/G/GA/GAAS/IO-String-1.08.tar.gz
tar -xvf IO-String-1.08.tar.gz


wget http://search.cpan.org/CPAN/authors/id/M/MH/MHOSKEN/Font-TTF-1.05.tar.gz
tar -xvf Font-TTF-1.05.tar.gz 


wget http://search.cpan.org/CPAN/authors/id/S/SS/SSIMMS/PDF-API2-2.023.tar.gz
tar -xvf PDF-API2-2.023.tar.gz

cd IO-String-1.08
perl Makefile.PL PREFIX=$PERL5LIBPREFIX
make
make test
make install
cd ..

cd Font-TTF-1.05
perl Makefile.PL PREFIX=$PERL5LIBPREFIX
make
make test
make install
cd ..

cd PDF-API2-2.023
perl Makefile.PL PREFIX=$PERL5LIBPREFIX
make
make test
make install
cd ../.. 
 
# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .

# Determine perl version


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  
}
module-whatis \"$SOFTWARE  \" ;  
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/mirdeep$VERSION ;  
prepend-path    PATH                \$root ; 
prepend-path    PATH                \$root/randfold ;

module          load                mugqic/bowtie/1.0.0
module          load                mugqic_dev/ViennaRNA/1.8.3

prepend-path    PERL5LIB            $PERL5LIB ; 

" > $VERSION


# (echo 'export PERL5LIB=PERL5LIB:your_path_to_miRDeep2/lib/perl5/site_perl/5.x/' >> ~/.bashrc)

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
