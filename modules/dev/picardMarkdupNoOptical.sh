#!/bin/sh

#
# Picard
#

SOFTWARE=picard
VERSION="1.90MarkDupNoOpt"
DOWNLOADVERSION=1.90
INSTALL_PATH=$MUGQIC_INSTALL_HOME_DEV/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
wget http://downloads.sourceforge.net/project/picard/picard-tools/$DOWNLOADVERSION/$SOFTWARE-tools-$DOWNLOADVERSION.zip
unzip $SOFTWARE-tools-$DOWNLOADVERSION.zip

# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
mv $SOFTWARE-tools-$DOWNLOADVERSION  $SOFTWARE-tools-$VERSION
mv -i $SOFTWARE-tools-$VERSION $INSTALL_PATH
mv -i $SOFTWARE-tools-$DOWNLOADVERSION.zip $MUGQIC_INSTALL_HOME_DEV/archive

cp /lb/project/mugqic/projects/bayfield_PRJBFX_389/MarkDuplicates-noOptical-v1.90.jar $INSTALL_PATH/$SOFTWARE-tools-$VERSION/MarkDuplicates.jar

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;
}
module-whatis \"$SOFTWARE  \" ;
                      
set             root                \$::env(MUGQIC_INSTALL_HOME_DEV)/software/$SOFTWARE/$SOFTWARE-tools-$VERSION ;
setenv          PICARD_HOME         \$root ;
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE
chmod -R ug+rwX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
