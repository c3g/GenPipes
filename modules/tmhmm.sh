#!/bin/bash

## NOTE (FL): 
# - Tools from cbs.dtu.dk under some sort of danish license and need to be downloaded manually.
# - Assumes gnuplot, gs and are avail on the system
# - Empty output problem: http://sourceforge.net/p/trinotate/mailman/message/31228028/
# 

SOFTWARE=tmhmm
VERSION=2.0c
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$SOFTWARE
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

# Download, extract, build
wget https://www.dropbox.com/s/2lgbq3ci3zvvm2d/tmhmm-$VERSION.Linux.tar.gz -O tmhmm-$VERSION.Linux.tar.gz
tar xvf tmhmm-$VERSION.Linux.tar.gz
cd $SOFTWARE-$VERSION 

# ppm2gif, could not find this on the web. Link to ppm2tiff may work instead
ln -s `which ppm2tiff` bin/ppm2gif

# Stupid file tweaks again
sed -i s,"#\!/usr/local/bin/perl,#\!/usr/bin/env perl,g" bin/tmhmm
sed -i s,"#\!/usr/local/bin/perl -w,#\!/usr/bin/env perl,g" bin/tmhmmformat.pl

 # You need an executable of the program decodeanhmm that runs under
 # Unix. The program may already be in bin/decodeanhmm.
 # 
 # The scripts require perl 5.x
 # 
 # For plotting gnuplot is needed (making postscript plots).
 # 
 # When generating html output the postscript plots are converted to
 # gif, and for this you need the programs ghostscript (gs) and ppmtogif.
 # 
 # After unpacking the directory you should
 # 
 # 1. Insert the correct path for perl 5.x in the first line of the scripts
 #    bin/tmhmm and bin/tmhmmformat.pl (if not /usr/local/bin/perl).
 # 2. Make sure you have an executable version of decodeanhmm in the bin
 #    directory.
 # 3. Include the directory containing tmhmm in your path.
 # 4. Read the TMHMM2.0.guide.html.
 # 5. Run the program.
 # 
 # 

cd ..


# Add permissions and install software
chmod -R ug+rwX .
chmod -R o+rX .
mv -if tmhmm-$VERSION.Linux.tar.gz $MUGQIC_INSTALL_HOME/archive/ 



# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \" ; 
}
module-whatis \"$SOFTWARE  \" ; 

set             root                \$::env(MUGQIC_INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE-$VERSION ;  
prepend-path    PATH                \$root/bin ;
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


