#!/bin/bash

#
# MUGQIC pipeline
#

SOFTWARE=mugqic_pipeline
VERSION=1.4
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
ARCHIVE=$MUGQIC_INSTALL_HOME/archive/$SOFTWARE
mkdir -p $INSTALL_DOWNLOAD $ARCHIVE
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://bitbucket.org/mugqic/$SOFTWARE/get/$VERSION.tar.gz
tar zxvf $VERSION.tar.gz
mv mugqic-$SOFTWARE-* v$VERSION

# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -i v$VERSION $INSTALL_PATH
mv -i $VERSION.tar.gz $ARCHIVE

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                   \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/v$VERSION ;
setenv          MUGQIC_PIPELINE_HOME   \$root
prepend-path    PATH                   \$root/pipelines/chipseq
prepend-path    PATH                   \$root/pipelines/dnaseq
prepend-path    PATH                   \$root/pipelines/rnaseq
prepend-path    PATH                   \$root/pipelines/rnaseq_denovo
prepend-path    PERL5LIB               \$root/lib
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Well... here, module directory is named "pipeline" instead of "mugqic_pipeline" for aesthetical reasons

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
MODULE_DIR=$MUGQIC_INSTALL_HOME/modulefiles/mugqic/pipeline
mkdir -p $MODULE_DIR
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MODULE_DIR

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
