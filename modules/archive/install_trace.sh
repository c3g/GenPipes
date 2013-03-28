# TODO:
# Commit to bitbucket
# Figure out a way to mke stubs accessigle from bitbucket
# Reflect on whther or not to use MUGQIC_HOME.. how installtion script or within script guess with enforced structure
# make java jars executable with stubs? --> then need to solve java_args problem, flag...

#cd ~
#echo "hello world" > temp.txt
#mv temp.txt /sb/programs/analyste/
#cd /sb/programs/analyste/
#chmod g+rw temp.txt
#chgrp analyste temp.txt
#ls -alrth temp.txt
#
#chmod -R g+rw ./*
#chgrp -R analyste ./*
#chmod -R a+r ./*



#### Epilogue
umask 0002
## mugic module
HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;
if [[ $HOST == abacus* ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $HOST == lg-* ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $HOST == grid1.genome.mcgill.ca ]]; then
 export MUGQIC_INSTALL_HOME=/home/nfs/flefebvr
elif [[ $HOST == ip03 ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/scratch_mp2/bourque/bourque_group/opt
elif [[ $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $DNSDOMAIN == guillimin.clumeq.ca ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $DNSDOMAIN == m ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/scratch_mp2/bourque/bourque_group/opt
fi
module use $MUGQIC_INSTALL_HOME/modulefiles

# 


###################
################### MUGQIC TOOLS (hosted on svn for now)
###################
VERSION="0.1"
screen -S svn
ssh -l flefebvr  -L9443:esx-svn.genome.mcgill.ca:443 gallium.genome.mcgill.ca
# ctl+A +D, then enter when asking pwd
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mugqic_tools # where to install..
mkdir -p $INSTALL_PATH
cp -r bioinformatics/R-tools bioinformatics/perl-tools bioinformatics/java-tools bioinformatics/tools $INSTALL_PATH 

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - MUGQIC developped tools \"
}
module-whatis \"MUGQIC - MUGQIC developped tools \"
                       
set             root            \$::env(MUGQIC_INSTALL_HOME)/software/mugqic_tools
prepend-path    PATH            \$root/tools
prepend-path    PATH            \$root/perl-tools
setenv          R_TOOLS         \$root/R-tools 

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tools


# mugqic/tools



###################
################### R
###################

## Install R itself (libcairo must be installed?s)
VERSION="2.15.2"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/R/R-$VERSION # where to install..
mkdir -p $INSTALL_PATH
wget http://cran.r-project.org/src/base/R-${VERSION:0:1}/R-$VERSION.tar.gz
tar -xvf R-$VERSION.tar.gz
cd R-$VERSION
./configure --prefix=$INSTALL_PATH  # TEMP s--with-readline=yes --with-readline=no
make -j8
make install
# cleanup

## Install prefered add on packages
DEP_PATH="https://bitbucket.org/mugqic/rpackages/raw/6608eeca25d1ca70352400aed045d73fe1c94820/DEPENDENCIES/package_dependencies_list.txt" # path to list of prefered packages (if possible change that to bitbucket link in the future...)
wget $DEP_PATH
$INSTALL_PATH/bin/Rscript -e "source(\"http://bioconductor.org/biocLite.R\");deps=readLines(\"package_dependencies_list.txt\");biocLite(deps,lib=.Library);biocLite(deps,lib=.Library)" # runs twice, sometimes mult attempts necessary



## chmod after installtion
chmod -R g+rw $INSTALL_PATH/lib64/R/library

#bin/Rscript -e "library(BiocInstaller);deps=readLines(\"$DEP_PATH\");biocLite(deps,lib=.Library);biocLite(deps,lib=.Library)"  # This leverages biocLite argh readLines() not fetching bitbucket... why???
# http://dl.dropbox.com/u/2528754/gqRsuite/dependencies_list.txt
#bin/Rscript -e "deps=readLines(\"\");install.packages(deps,lib=.Library,repos=\"http://cran.us.r-project.org\")" 

# Module def file..
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds R to your environment \"
}
module-whatis \"MUGQIC - Adds R to your environment \"
                       
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/R/R-$VERSION
setenv          R_LIBS             \$root/lib64/R/library
prepend-path    MANPATH            \$root/share              
prepend-path    PATH               \$root/bin
prepend-path    LD_LIBRARY_PATH    \$root/lib64:/software/libraries/GotoBLAS_LAPACK/shared
#prepend-path   LD_LIBRARY_PATH    \$root/lib64:\$root/standalone:/software/libraries/GotoBLAS_LAPACK/shared
#prepend-path   CPATH              \$root/include

" > $VERSION

## THEN--> Move module definition manually, and edit .version
# set             root               \$::env(MUGQIC_INSTALL_HOME)/software/R/R-$VERSION


# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R






###################
################### tophat
###################
VERSION="2.0.7"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/tophat/tophat-$VERSION.Linux_x86_64 # where to install..
mkdir -p $INSTALL_PATH
wget http://tophat.cbcb.umd.edu/downloads/tophat-$VERSION.Linux_x86_64.tar.gz
tar -xvf tophat-$VERSION.Linux_x86_64.tar.gz
mv tophat-$VERSION.Linux_x86_64/* $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds tophat to your environment \"
}
module-whatis \"MUGQIC - Adds tophat to your environment \"
                       
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/tophat/tophat-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/tophat



###################
################### RNASeqC
###################
VERSION="1.1.7" # Tue 29 Jan 10:10:21 2013 
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/rnaseqc/RNA-SeQC_v$VERSION # where to install..
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v$VERSION.jar

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - RNAseq QC software. Depends on bwa \"
}
module-whatis \"MUGQIC -Java tool to generate RNA QC html report \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/rnaseqc/RNA-SeQC_v$VERSION/RNA-SeQC_v$VERSION.jar
setenv          RNASEQC_JAR        \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/rnaseqc
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/rnaseqc


###################
################### FASTQC
###################
VERSION="0.10.1"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/fastqc/fastqc_v"$VERSION" # where to install..
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v"$VERSION".zip
unzip fastqc_v"$VERSION".zip
rm fastqc_v"$VERSION".zip
chmod +x FastQC/fastqc

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - FASTQC \"
}
module-whatis \"MUGQIC -FASTQC \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/fastqc/fastqc_v"$VERSION"/FastQC
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastqc
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastqc/



###################
################### Cufflinks
###################
VERSION="2.0.2"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/cufflinks
mkdir -p $INSTALL_PATH
wget http://cufflinks.cbcb.umd.edu/downloads/cufflinks-$VERSION.Linux_x86_64.tar.gz
tar -xvf cufflinks-$VERSION.Linux_x86_64.tar.gz
mv cufflinks-$VERSION.Linux_x86_64 $INSTALL_PATH


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Cufflinks \"
}
module-whatis \"MUGQIC - Cufflinks \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/cufflinks/cufflinks-$VERSION.Linux_x86_64
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/cufflinks
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/cufflinks/




###################
################### sqlite3
###################
VERSION="3071502"
INSTALL_PATH="$MUGQIC_INSTALL_HOME/software/sqlite-shell-linux-x86-"$VERSION
mkdir -p $INSTALL_PATH
FILE="sqlite-shell-linux-x86-"$VERSION".zip" #http://www.sqlite.org/sqlite-shell-linux-x86-3071502.zip
wget http://www.sqlite.org/$FILE
unzip $FILE
mv sqlite3  $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - sqlite3 \"
}
module-whatis \"MUGQIC - sqlite3 shell \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/sqlite-shell-linux-x86-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/sqlite3
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/sqlite3/



###################
################### HMMER3
###################
VERSION="3.0"
INSTALL_PATH="$MUGQIC_INSTALL_HOME/software/hmmer/hmmer-"$VERSION
mkdir -p $INSTALL_PATH
FILE="hmmer-"$VERSION"-linux-intel-x86_64.tar.gz" 
wget "http://selab.janelia.org/software/hmmer3/$VERSION/$FILE" 
tar -xvf $FILE
cd "hmmer-$VERSION-linux-intel-x86_64"
./configure
make
make check
cd ..
mv "hmmer-$VERSION-linux-intel-x86_64" $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - hmmer3 \"
}
module-whatis \"MUGQIC - hmmer3 \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/hmmer/hmmer-$VERSION/hmmer-$VERSION-linux-intel-x86_64/binaries
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/hmmer/










###################
################### fastx
###################
# /sb/programs/analyste/software/fastx_toolkit-0.0.13.2
# http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit-0.0.13.2.tar.bz2
# http://hannonlab.cshl.edu/fastx_toolkit/libgtextutils-0.6.1.tar.bz2
# http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13.2_binaries_Linux_2.6_amd64.tar.bz2
VERSION="0.0.13.2"
LIBVERSION="0.6.1"
NAME=fastx_toolkit-$VERSION
LIBNAME=libgtextutils-$LIBVERSION
wget http://hannonlab.cshl.edu/fastx_toolkit/$NAME.tar.bz2
wget http://hannonlab.cshl.edu/fastx_toolkit/$LIBNAME.tar.bz2
tar -xvf $NAME.tar.bz2
tar -xvf $LIBNAME.tar.bz2

INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/fastx/$NAME
mkdir -p $INSTALL_PATH

cd $LIBNAME
./configure --prefix=$INSTALL_PATH
make
make install
cd ..

export PKG_CONFIG_PATH=$INSTALL_PATH/lib/pkgconfig
cd $NAME
./configure --prefix=$INSTALL_PATH
make
make install
cd ..

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - fastx \"
}
module-whatis \"MUGQIC - fastx \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/fastx/fastx_toolkit-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastx
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/fastx/




###################
################### BLAT
###################
VERSION="35x1"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/blat/blat-$VERSION
mkdir -p $INSTALL_PATH
wget "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat"
chmod +x blat
mv blat $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BLAT \"
}
module-whatis \"MUGQIC - BLAT \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/blat/blat-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blat
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/blat/




###################
################### Bowtie2
###################
VERSION="2.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bowtie/bowtie-$VERSION
mkdir -p $INSTALL_PATH
# Download and extract
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$VERSION/bowtie2-$VERSION-source.zip/download
unzip bowtie2-$VERSION-source.zip
# Compile
cd bowtie2-$VERSION
make -j8
mv bowtie2* $INSTALL_PATH
cd ..


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Bowtie2 aligner \"
}
module-whatis \"MUGQIC - Bowtie2 aligner \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bowtie/bowtie-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bowtie
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bowtie/



###################
################### bedGraphToBigWig
###################
VERSION="v4"
NAME=bedGraphToBigWig # same could apply to all ucsc tools
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME/$NAME-$VERSION
mkdir -p $INSTALL_PATH
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/$NAME
chmod +x $NAME
mv $NAME $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $NAME \"
}
module-whatis \"MUGQIC - $NAME \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/$NAME/$NAME-$VERSION
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$NAME/





###################
################### BEDTools
###################
VERSION="2.17.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/bedtools/bedtools-$VERSION
mkdir -p $INSTALL_PATH
wget "http://bedtools.googlecode.com/files/BEDTools.v$VERSION.tar.gz"
tar -xvf BEDTools.v$VERSION.tar.gz
cd bedtools-$VERSION
make -j8
mv ./* $INSTALL_PATH



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BEDtools \"
}
module-whatis \"MUGQIC - BEDtools  \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/bedtools/bedtools-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/bedtools/






###################
################### Trimmomatic
###################
VERSION="0.25"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$VERSION.zip 
unzip Trimmomatic-$VERSION.zip
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/trimmomatic
mkdir -p $INSTALL_PATH
mv Trimmomatic-$VERSION $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - BEDtoolsTrimmomatic to trim fastqs \"
}
module-whatis \"Trimmomatic to trim fastqs  \"
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/trimmomatic/Trimmomatic-$VERSION
setenv          TRIMMOMATIC_JAR     \$root/trimmomatic-$VERSION.jar

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/trimmomatic/




###################
################### IDBA
###################
VERSION="1.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/idba/
mkdir -p $INSTALL_PATH


# IDBA : compile and test whether it segfaults
wget http://hku-idba.googlecode.com/files/idba-$VERSION.tar.gz
tar -xvf idba-$VERSION.tar.gz
cd idba-$VERSION



./configure --prefix=$INSTALL_PATH/idba-$VERSION
make -j8
make install
cp -rf bin/* $INSTALL_PATH/idba-$VERSION/bin


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IDBA assembler \"
}
module-whatis \"MUGQIC - IDBA assembler  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/idba/idba-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba/






###################
################### IDBA LONG
###################
VERSION="1.0.9"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/idba_long/
mkdir -p $INSTALL_PATH


# IDBA : compile and test whether it segfaults
wget http://hku-idba.googlecode.com/files/idba_ud-$VERSION.tar.gz
tar -xvf idba_ud-$VERSION.tar.gz
cd idba_ud-$VERSION

#wget http://hku-idba.googlecode.com/files/idba-$VERSION.tar.gz
#tar -xvf idba-$VERSION.tar.gz
#cd idba-$VERSION

# TODO MANUAL: Tweak for longer reads https://groups.google.com/forum/?fromgroups=#!topic/hku-idba/ShB95FPswN8
# https://groups.google.com/forum/?fromgroups=#!topic/hku-idba/NE2JXqNvTFY


./configure --prefix=$INSTALL_PATH/idba-$VERSION
make -j8
make install
cp -rf bin/* $INSTALL_PATH/idba-$VERSION/bin


# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - IDBA assembler \"
}
module-whatis \"MUGQIC - IDBA assembler  \"
            
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/idba_long/idba-$VERSION/bin
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba_long
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/idba_long/





###################
################### A5
###################
VERSION="20120518"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/a5/
mkdir -p $INSTALL_PATH
wget http://ngopt.googlecode.com/files/ngopt_a5pipeline_linux-x64_$VERSION.tar.gz
tar -xvf ngopt_a5pipeline_linux-x64_$VERSION.tar.gz
# sh ngopt_a5pipeline_linux-x64_$VERSION/test.a5.sh # only tested on login node


# TODO manual : modify the perl script to support long reads... -l + -r with empty file
	# PATCH: create dumy emptyfile	
	#my $idba_cmd = "idba -r $reads -o $WD/$outbase --mink ".IDBA_MIN_K." --maxk $maxrdlen";
	#my @dummy = `echo "" > empty.fa`;
	#my $idba_cmd = "idba -r empty.fa -l $reads -o $WD/$outbase --mink ".IDBA_MIN_K." --maxk $maxrdlen";
	

mv ngopt_a5pipeline_linux-x64_$VERSION $INSTALL_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - A5 microbial genomes assembly pipeline \"
}
module-whatis \"MUGQIC - A5  \"
 
conflict        mugqic/bwa mugqic/fish mugqic/idba mugqic/samtools mugqic/sga mugqic/sspace mugqic/tagdust mugqic/bowtie mugqic/SSPACE             
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/a5/ngopt_a5pipeline_linux-x64_$VERSION/bin
prepend-path    PATH               \$root
prepend-path    PATH               \$root/SSPACE
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/a5
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/a5/










###################
################### Python
###################
VERSION="2.7.3"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/python/Python-$VERSION
mkdir -p $INSTALL_PATH
wget "http://www.python.org/ftp/python/$VERSION/Python-$VERSION.tgz"
tar -xvf Python-$VERSION.tgz
cd Python-$VERSION
./configure --prefix=$INSTALL_PATH
make -j8
make install


# Module file (TODO: BLAS and LAPACK are not portable)
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - Adds Python $VERSION to your environment \"
}
module-whatis \"Adds Python $VERSION to your environment  \"

setenv  BLAS  			/software/libraries/GotoBLAS_LAPACK/shared/libblas.so
setenv  LAPACK		   /software/libraries/GotoBLAS_LAPACK/shared/liblapack.so

set             root               \$::env(MUGQIC_INSTALL_HOME)/software/python/Python-$VERSION
prepend-path    MANPATH            \$root/share/man              
prepend-path    PATH               \$root/bin
prepend-path    LD_LIBRARY_PATH    \$root/lib
prepend-path    CPATH              \$root/include

#setenv		PYTHONPATH	\$root/lib
#setenv		PYTHONHOME	\$root/lib
# setenv                INTEL_LICENSE_FILE      $root/licenses
# setenv                FC                      gfortran
# setenv                F77                     gfortran
# setenv                CC                      gcc
#prepend-path    PATH               \$root/bin:/software/tools/swig-2.0.4/bin:/software/tools/wx-2.8.12/bin
#prepend-path    CPATH              $root/include:/software/tools/wx-2.8.12/include
#prepend-path    LD_LIBRARY_PATH    /software/libraries/GotoBLAS_LAPACK/shared:$root/lib:/software/tools/wx-2.8.12/lib


" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/python/


## Install numpy
module load mugqic/python



## Install HTSeq
module load mugqic/python










###################
################### Ray
###################
#module load compat-openmpi-x86_64 # 4abacus
#module load compat-openmpi-psm-x86_64
VERSION="2.1.0"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/ray/Ray-v$VERSION
NAME=Ray-v$VERSION
wget http://sourceforge.net/projects/denovoassembler/files/$NAME.tar.bz2/download
tar -xvf $NAME.tar.bz2
cd $NAME
# ...... TODO: meh
make PREFIX=ray-build
make install
ls ray-build
mpiexec -n 1 ray-build/Ray -o test -p test_1.fastq test_2.fastq -k 31
#mkdir -p $INSTALL_PATH
Then, type
== Options ==
You can provide compilation options to the Makefile.
MPICXX                  The path to the C++ compiler wrapper (usually called mpicxx)
PREFIX                  Where to install stuff
MAXKMERLENGTH           maximum k-mer length, default is MAXKMERLENGTH=32
FORCE_PACKING           save memory by not aligning addresses, default is FORCE_PACKING=n
ASSERT                  run assertions too, default is ASSERT=n
For other options, read the Makefile header.




########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################



#!/bin/sh

## F Lefebvre
## Tue  1 Jan 11:50:34 2013
##  
## Download, compile and install external tools..
## Versions are tracked here in file versions.csv which can be used by reporting mecanism
## The guiding principle should be to avoid external tools as much as possible and stick to suitable
## equivalent tools from the R/Bioconductor ecosystem.
##
## This script is run upon package installation. It downloads the source, compiles, registers the 
## tool versions to a file as well as the path to binaries that should be appended to PATH
## upon package load (see onLoad function)
##
## 
##
## JAR are made executables by pre-pending a shell script (inst/scripts/stub.sh). Reference:
## https://coderwall.com/p/ssuaxa
##  
## 
UNAME=`uname`
rm -rf inst/tools/* # pre-cleanup, because devtools run this script twice!
cd inst/tools

# rm -rf ./*
### Initialize version track file and the paths lists
echo "Tool,Version,Note" > versions.csv
echo -n "" > paths.txt


#########  Trimmomatic (java jar classpath) ##################
VERSION="0.22"
# Download and extract
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$VERSION.zip 
unzip Trimmomatic-$VERSION.zip
# Create executable by prepending a shell script (java_classpath_stub.sh) to the binary..
cat ../scripts/java_classpath_stub.sh Trimmomatic-$VERSION/trimmomatic-$VERSION.jar > Trimmomatic-$VERSION/trimmomatic
chmod +x Trimmomatic-$VERSION/trimmomatic # usage: ./trimmomatic org.usadellab.trimmomatic.TrimmomaticPE 
# Add to paths list
echo  Trimmomatic-$VERSION >> paths.txt
# Register tool version
echo "Trimmomatic,$VERSION,http://www.usadellab.org/cms/index.php?page=trimmomatic" >> versions.csv


#%Module1.0

proc ModulesHelp { } {
        puts stderr "\tadd  Trimmomatic to trim fastqs"
}

module-whatis "Trimmomatic to trim fastqs"

set             root                $::env(MUGQIC_INSTALL_HOME)/software/Trimmomatic-0.25
setenv          TRIMMOMATIC_JAR     $root/trimmomatic-0.25.jar



######### samtools ##################
VERSION="0.1.18"
# Download and extract
wget http://sourceforge.net/projects/samtools/files/samtools/$VERSION/samtools-$VERSION.tar.bz2/download
tar xvjf samtools-$VERSION.tar.bz2
# Compile
cd samtools-$VERSION
make -j8
cd ..
# Add to paths list
echo samtools-$VERSION >> paths.txt
echo samtools-$VERSION/bcftools >> paths.txt
# Register tool version
echo "samtools,$VERSION,http://samtools.sourceforge.net" >> versions.csv


######### picard (java jars) ##################
VERSION="1.82"
# Download and extract
wget http://sourceforge.net/projects/picard/files/picard-tools/$VERSION/picard-tools-$VERSION.zip/download
unzip picard-tools-$VERSION.zip
# Create executable by prepending a shell script (jar_stub.sh) to the binary..
for f in `ls picard-tools-$VERSION/*.jar `; do
	filename=$(basename "$f")
	toolname="${filename%.*}"
	cat ../scripts/jar_stub.sh $f >  picard-tools-$VERSION/$toolname
	chmod +x picard-tools-$VERSION/$toolname
done
# Add to paths list
echo picard-tools-$VERSION >> paths.txt
# Register tool version
echo "picard,$VERSION,http://picard.sourceforge.net" >> versions.csv



######### bwa ##################
VERSION="0.6.2"
# Download and extract
wget http://sourceforge.net/projects/bio-bwa/files/bwa-$VERSION.tar.bz2/download
tar xvjf bwa-$VERSION.tar.bz2
# Compile
cd bwa-$VERSION
make -j8
cd ..
# Append to path
echo bwa-$VERSION >> paths.txt
# Register tool version
echo "bwa,$VERSION,http://bio-bwa.sourceforge.net" >> versions.csv



######### bowtie2 ##################
VERSION="2.1.0"
# Download and extract
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$VERSION/bowtie2-$VERSION-source.zip/download
unzip bowtie2-$VERSION-source.zip
# Compile
cd bowtie2-$VERSION
make -j8
cd ..
# append to path
echo bowtie2-$VERSION >> paths.txt
# Register tool version
echo "bowtie2,$VERSION,http://bowtie-bio.sourceforge.net/bowtie2" >> versions.csv





######### tophat2 ################## 
VERSION="2.0.6"
# Download and extract (tophat is a pain in the a** to compile)
if [ "$UNAME" == "Linux" ]; then
   os="Linux"
fi
if [ "$UNAME" == "Darwin" ]; then 
   os="OSX"
fi
wget http://tophat.cbcb.umd.edu/downloads/tophat-$VERSION."$os"_x86_64.tar.gz
tar -xvf tophat-$VERSION."$os"_x86_64.tar.gz
# paths
echo tophat-$VERSION."$os"_x86_64 >> paths.txt
# Register tool version
echo "tophat,$VERSION,http://tophat.cbcb.umd.edu" >> versions.csv



######### cufflinks and co. ##################
VERSION="2.0.2"
# Download and extract 
if [ "$UNAME" == "Linux" ]; then
   os="Linux"
fi
if [ "$UNAME" == "Darwin" ]; then 
   os="OSX"
fi
wget http://cufflinks.cbcb.umd.edu/downloads/cufflinks-$VERSION."$os"_x86_64.tar.gz
tar -xvf cufflinks-$VERSION."$os"_x86_64.tar.gz
# paths
echo cufflinks-$VERSION."$os"_x86_64 >> paths.txt
# Register tool version
echo "cufflinks,$VERSION,http://cufflinks.cbcb.umd.edu" >> versions.csv




######### tophat2-fusion (already in tophat2!)##################
# ...


######### STAR ##################
VERSION="2.2.0"
# Download and extract
wget http://rna-star.googlecode.com/files/STAR_"$VERSION"c.tgz
tar xvzf STAR_"$VERSION"c.tgz
# Compile
cd STAR_"$VERSION"c
make -j8
cd ..
# paths
echo STAR_"$VERSION"c >> paths.txt
# Register tool version
echo "STAR,$VERSION,http://code.google.com/p/rna-star/" >> versions.csv


######### RNASeQC (java jar) ##################
VERSION="1.1.7"
# Download and extract
mkdir RNA-SeQC_v$VERSION
cd RNA-SeQC_v$VERSION
wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v$VERSION.jar
cd ..
# Create executable by prepending a shell script (jar_stub.sh) to the binary..
cat ../scripts/jar_stub.sh RNA-SeQC_v$VERSION/RNA-SeQC_v$VERSION.jar > RNA-SeQC_v$VERSION/RNA-SeQC
chmod +x RNA-SeQC_v$VERSION/RNA-SeQC
# version abstracted link
echo RNA-SeQC_v$VERSION >> paths.txt
# Register tool version
echo "RNA-SeQC,$VERSION,https://confluence.broadinstitute.org/display/CGATools/RNA-SeQC" >> versions.csv




######### FastQC (java jar) ##################
VERSION="0.10.1"
# Download and extract
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v"$VERSION".zip
unzip fastqc_v"$VERSION".zip
# link to executable
chmod +x FastQC/fastqc
# paths
echo FastQC >> paths.txt
# Register tool version
echo "FastQC,$VERSION,http://www.bioinformatics.babraham.ac.uk/projects/fastqc" >> versions.csv





######### Final Cleanup
rm *.zip *.gz *.bz2 *.jar *.tgz
cd ..






######### RSeQC ##################
#### What else???? whatever other tools out there..., BLAST, BFAST, Bowtie1, MMummer ####
# BLAST
#BFAST
#LastZ
#Novoalign
#Blat
# Mummer







