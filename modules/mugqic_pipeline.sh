###################
################### MUGQIC pipeline 
###################
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/pipeline/tmp/untar/
cd $MUGQIC_INSTALL_HOME/modulefiles/mugqic/pipeline/tmp
VERSION="1.0"
wget https://bitbucket.org/mugqic/mugqic_pipeline/get/${VERSION}.tar.gz
tar -xvf ${VERSION}.tar.gz -C untar/
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/mugqic_pipeline # where to install..
ARCHIVE_PATH=$MUGQIC_INSTALL_HOME/archive/mugqic_pipeline 
mkdir -p $INSTALL_PATH $ARCHIVE_PATH
cp -r untar/mugqic-mugqic_pipeline-*/*  $INSTALL_PATH
chmod -R 775 $INSTALL_PATH 
mv ${VERSION}.tar.gz $ARCHIVE_PATH

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - MUGQIC developped tools \"
}
module-whatis \"MUGQIC - MUGQIC developped tools \"
                       
set             root                   \$::env(MUGQIC_INSTALL_HOME)/software/mugqic_pipeline
setenv          MUGQIC_PIPELINE_HOME   \$root
prepend-path    PATH                   \$root/pipelines/chipseq
prepend-path    PATH                   \$root/pipelines/dnaseq
prepend-path    PATH                   \$root/pipelines/rnaseq
prepend-path    PERL5LIB               \$root/lib

" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"

" > .version

mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/pipeline
cd ..
rm -rf tmp
