#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#
# Perl
#

SOFTWARE=perl
VERSION=5.18.2

# 'MUGQIC_INSTALL_HOME_DEV' for development, 'MUGQIC_INSTALL_HOME' for production (don't write '$' before!)
INSTALL_HOME=MUGQIC_INSTALL_HOME

# Indirection call to use $INSTALL_HOME value as variable name
INSTALL_DIR=${!INSTALL_HOME}/software/$SOFTWARE

# Create install directory with permissions if necessary
if [[ ! -d $INSTALL_DIR ]]
then
  mkdir $INSTALL_DIR
  chmod ug+rwX,o+rX-w $INSTALL_DIR
fi

INSTALL_DOWNLOAD=$INSTALL_DIR/tmp
mkdir -p $INSTALL_DOWNLOAD
cd $INSTALL_DOWNLOAD

# Download, extract, build
ARCHIVE=$SOFTWARE-$VERSION.tar.gz
# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f ${!INSTALL_HOME}/archive/$ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in ${!INSTALL_HOME}/archive/: using it..."
  cp -a ${!INSTALL_HOME}/archive/$ARCHIVE .
else
  echo "Archive $ARCHIVE not in ${!INSTALL_HOME}/archive/: downloading it..."
  wget http://www.cpan.org/src/5.0/$ARCHIVE
fi
tar zxvf $ARCHIVE

SOFTWARE_DIR=$SOFTWARE-$VERSION
cd $SOFTWARE_DIR
./Configure -des -Dusethreads -Dprefix=$INSTALL_DIR/$SOFTWARE_DIR
make
make install

# Install Perl modules from CPAN

# MakeMaker's prompt function will always return the default without waiting for user input.
export PERL_MM_USE_DEFAULT=1
# Install module prerequisites automatically
export PERL_EXTUTILS_AUTOINSTALL="--defaultdeps"

for MODULE in \
Config::Simple \
File::Slurp \
File::Which \
Filesys::Df \
Parse::Range \
PDF::Create \
PerlIO::gzip \
Proc::ParallelLoop \
Statistics::Descriptive \
Text::CSV \
Text::CSV::Encoded \
XML::Simple \
; do
$INSTALL_DIR/$SOFTWARE_DIR/bin/perl -MCPAN -e"CPAN::Shell->force(qw(install $MODULE))"
# Test if module is properly installed
$INSTALL_DIR/$SOFTWARE_DIR/bin/perl -e "use $MODULE"
done

# Add permissions
cd $INSTALL_DOWNLOAD
chmod -R ug+rwX,o+rX-w . $INSTALL_DIR/$SOFTWARE_DIR
# Store archive if not already present or if different from the previous one
if [[ ! -f ${!INSTALL_HOME}/archive/$ARCHIVE || `diff ${!INSTALL_HOME}/archive/$ARCHIVE $ARCHIVE` ]]
then
  mv -i $ARCHIVE ${!INSTALL_HOME}/archive/
fi

# Module file
echo "#%Module1.0
proc ModulesHelp { } {
  puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"$SOFTWARE\"

set             root                \$::env($INSTALL_HOME)/software/$SOFTWARE/$SOFTWARE_DIR
prepend-path    PATH                \$root/bin
prepend-path    PERL5LIB            \$root/lib
prepend-path    PERL5LIB            \$root/lib/$VERSION
prepend-path    PERL5LIB            \$root/lib/site_perl/$VERSION
" > $VERSION

################################################################################
# Everything below this line should be generic and not modified

# Default module version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"" > .version

# Set module directory path by lowercasing $INSTALL_HOME and removing '_install_home' in it
MODULE_DIR=${!INSTALL_HOME}/modulefiles/`echo ${INSTALL_HOME,,} | sed 's/_install_home//'`/$SOFTWARE

# Create module directory with permissions if necessary
if [[ ! -d $MODULE_DIR ]]
then
  mkdir $MODULE_DIR
  chmod ug+rwX,o+rX-w $MODULE_DIR
fi

# Add permissions and install module
chmod ug+rwX,o+rX-w $VERSION .version
mv $VERSION .version $MODULE_DIR/

# Clean up temporary installation files if any
cd
rm -rf $INSTALL_DOWNLOAD
