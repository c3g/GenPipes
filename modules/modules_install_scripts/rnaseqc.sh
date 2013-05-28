
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


