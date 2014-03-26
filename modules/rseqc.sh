#!/bin/bash

# Prerequisite: gcc; python2.7; numpy; R
# 
# Install RSeQC (Example):
# 
# tar zxf RSeQC-VERSION.tar.gz
# 
# cd RSeQC-VERSION
# 
# #type 'python setup.py install --help' to see options
# python setup.py install
# 
# #This is only an example. Change path according to your system configuration
# export PYTHONPATH=/home/user/lib/python2.7/site-packages:$PYTHONPATH
# 
# #This is only an example. Change path according to your system configuration
# export PATH=/home/user/bin:$PATH
# Finally, type: python -c ‘from qcmodule import SAM’. If no error message comes out, RSeQC modules have been installed successfully.
# 


#
# Software_name  ## TO BE MODIFIED WITH e.g. BLAST, HMMER, SAMtools, etc.
#
SOFTWARE=rseqc  
VERSION=2.3.8  ## TO BE MODIFIED WITH e.g. 2.2.28+, 3.0, 0.1.19, etc.
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
INSTALL_DOWNLOAD=$INSTALL_PATH/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
# Write here the specific commands to download, extract, build the software, typically similar to:
wget https://downloads.sourceforge.net/project/rseqc/RSeQC-$VERSION.tar.gz
tar xvf RSeQC-$VERSION.tar.gz
cd RSeQC-$VERSION
module load mugqic/python/2.7.5 # wouldn't work with 2.7.8s
python setup.py install 


# Add permissions and install software
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX .
chmod -R o+rX .
mv -i $SOFTWARE-$VERSION $INSTALL_PATH  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
mv -i $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) $MUGQIC_INSTALL_HOME/archive  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION.(zip|tar.gz|tar.bz2) IF DIFFERENT

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
}
module-whatis \"$SOFTWARE  \" ;  ## TO BE MODIFIED WITH DETAILED DESCRIPTION IF ANY
                      
set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  ## TO BE MODIFIED WITH SPECIFIC $SOFTWARE-$VERSION IF DIFFERENT
prepend-path    PATH                \$root/bin ;  ## TO BE ADDED IF NECESSARY
prepend-path    PATH                \$root/other_tools/bin ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
setenv          ${SOFTWARE}_JAR     \$root/$SOFTWARE-$VERSION.jar ;  ## TO BE ADDED AND MODIFIED IF NECESSARY
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Add permissions and install module
mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
chmod -R ug+rwX $VERSION .version
chmod -R o+rX $VERSION .version
mv $VERSION .version $MUGQIC_INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

# Clean up temporary installation files if any
rm -rf $INSTALL_DOWNLOAD
