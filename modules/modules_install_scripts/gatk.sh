

###################
################### samtools
###################
echo "You must download and install GATK manually";
echo "It can be found here:"
echo "http://www.broadinstitute.org/gatk/download"
VERSION="2.5-2-gf57256b"
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/GenomeAnalysisTK/GenomeAnalysisTK-${VERSION}
mkdir -p $INSTALL_PATH
echo "Install it here: ${INSTALL_PATH}"

echo "#%Module1.0

proc ModulesHelp { } {
        puts stderr "\tadd  GenomeAnalysisTK"
}

module-whatis "The Broads toolsuite to work with resequencing"

set             root         $::env(MUGQIC_INSTALL_HOME)/software/GenomeAnalysisTK/GenomeAnalysisTK-${VERSION}
setenv          GATK_JAR     $root/GenomeAnalysisTK.jar
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK/

echo "Module is installed here: $MUGQIC_INSTALL_HOME/modulefiles/mugqic/GenomeAnalysisTK"
