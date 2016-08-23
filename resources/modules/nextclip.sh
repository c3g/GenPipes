#!/bin/bash
set -e
umask 0002
me=`basename $0`

## Default arg values
SOFTWARE="nextclip"
WHATIS="NextClip: an analysis and read preparation tool for Nextera long mate pair libraries"
VERSION="master"
MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE"
INSTALL_DIR="$MUGQIC_INSTALL_HOME_TMP/software/$SOFTWARE"
FORCE_INSTALL=false

## Parse arguments
usage()
{
cat << EOF
usage: $0 options

NextClip installation

OPTIONS:
   -v      Force a specifc commit instead of latest
   -f      Force re-install, i.e. overwrite module  and installation even if module already exists
   -m      Where to write the module file
   -i	   Installation directory. Actual install dir will be reassigned to <-i>/$SOFTWARE-version
   -h      Print this message

NOTES: 
   ...
EOF
}

while getopts “v:m:i:fh” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         v)
             VERSION=$OPTARG
             ;;
	     m)
	         MODULEFILE_DIR=$OPTARG
	         ;;
         i)
             INSTALL_DIR=$OPTARG
             ;;
         f)
             FORCE_INSTALL=true
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

## Tmp dir to work in
#SOFTWARE="nextclip" ; VERSION="master" ; MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_TMP/modulefiles/mugqic/$SOFTWARE" ;INSTALL_DIR="$MUGQIC_INSTALL_HOME_TMP/software/$SOFTWARE"; FORCE_INSTALL=false
TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR 
echo "Working in $TEMPDIR"

## Get commit
git clone https://github.com/richardmleggett/nextclip.git
cd nextclip
git checkout $VERSION
VERSION=`git log --pretty=format:'%h' -n 1` # rename to pretty commit hash number

# Compile
make all

## Paths, mkdirs
INSTALL_DIR=$INSTALL_DIR/$SOFTWARE-$VERSION
MODULEFILE="$MODULEFILE_DIR/$VERSION"
MODULEVERSIONFILE="$MODULEFILE_DIR/.version"

## Install if required by force or absence of module files
if  [ ! -f $MODULEFILE ] || $FORCE_INSTALL
then
	# Prelim. cleanup
	rm -rf $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
	mkdir -p $MODULEFILE_DIR $INSTALL_DIR
 
	# move
	cp -rf ./* $INSTALL_DIR/
	
	# Create module files
	cat > $MODULEFILE <<-EOF	
		#%Module1.0
		proc ModulesHelp { } {
		puts stderr "MUGQIC - $WHATIS"
		}
		module-whatis "MUGQIC - $WHATIS"
		
		set             root               $INSTALL_DIR          
		prepend-path    PATH               \$root/bin
		EOF
	cat > $MODULEVERSIONFILE <<-EOF
		#%Module1.0
		set ModulesVersion $VERSION
		EOF
fi

## Adjust permissions
chmod -R ug+rwX  $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
chmod -R o+rX    $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE

exit
