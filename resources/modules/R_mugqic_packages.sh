#!/bin/bash
set -e
umask 0002
me=`basename $0`

## Repo infos
REPO="rpackages" # in case the name of the repo changes
PCKGS="gqData gqUtils gqMicroarrays gqSeqUtils" # determines the order of installation (dependencies)
SOFTWARE="mugqic_R_packages" # in case we don't like the name

## Default arg values
#REF="master"
REF=1.0.6
INSTALL_PREFIX_ENV_VARNAME=""
MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev"
INSTALL_DIR="$MUGQIC_INSTALL_HOME_DEV/software"
R_MODULE="mugqic_dev/R_Bioconductor"

## Parse arguments
usage()
{
cat << EOF
usage: $0 options

This script installs a given mugqic in-house package from bitbucket.

OPTIONS:
   -v	The desired git reference to the repository. Defaults to master.
   -R	The R module name with which to build and prereq this package installation. Defaults to mugqic_dev/R_Bioconductor
   -p	Name of an environment variable which defines a prefix path to -m and -i. E.g. MUGQIC_INSTALL_HOME
   -m  	Path the module files directory. Defaults to \$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev. Actual module file  dir will be <-m>/$SOFTWARE/<-v>
   -i	library installation directory, defaults to \$MUGQIC_INSTALL_HOME_DEV/software. Actual install will be <-i>/$SOFTWARE/$SOFTWARE-<-v>/
   -h 	Print this message
   
EXAMPLE USAGE:

R_mugqic_packages.sh  -v master -R mugqic_dev/R_Bioconductor/3.1.2_3.0 -p MUGQIC_INSTALL_HOME_DEV -m modulefiles/mugqic_dev -i software

NOTES: 
...

EOF
}

while getopts “r:v:R:p:m:i:h” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
		 v)
		     REF=$OPTARG
		     ;;
	     R)
	         R_MODULE=$OPTARG
	         ;;
         p)
             INSTALL_PREFIX_ENV_VARNAME=$OPTARG
             ;;
	     m)
	         MODULEFILE_DIR=$OPTARG
	         ;;
         i)
             INSTALL_DIR=$OPTARG
             ;;
         ?)
             usage
             exit
             ;;
     esac
done



## Tmp dir to work in
TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR 
echo "Working in $TEMPDIR"


## Paths, mkdirs, prefix
export INSTALL_DIR="$INSTALL_DIR/$SOFTWARE/$SOFTWARE-$REF"
MODULEFILE_DIR="$MODULEFILE_DIR/$SOFTWARE"
TCLROOT=$INSTALL_DIR
if [[ $INSTALL_PREFIX_ENV_VARNAME != "" ]]
then
	echo "prefixing..."
	TCLROOT="\$::env($INSTALL_PREFIX_ENV_VARNAME)/$INSTALL_DIR"
	INSTALL_DIR=${!INSTALL_PREFIX_ENV_VARNAME}/$INSTALL_DIR
	MODULEFILE_DIR=${!INSTALL_PREFIX_ENV_VARNAME}/$MODULEFILE_DIR
fi
MODULEFILE="$MODULEFILE_DIR/$REF"
MODULEVERSIONFILE="$MODULEFILE_DIR/.version"

echo "The software install location is $INSTALL_DIR"
echo "The module file directory is $MODULEFILE_DIR"
echo "The module file is $MODULEFILE"
echo "The module version file is $MODULEVERSIONFILE"


## Dir creation
mkdir -p $MODULEFILE_DIR $INSTALL_DIR

# Download from bitbucket, install
wget "https://bitbucket.org/mugqic/$REPO/get/$REF.tar.gz" -O "$SOFTWARE-$REF.tar.gz" 
tar -zxf "$SOFTWARE-$REF.tar.gz"  
rm "$SOFTWARE-$REF.tar.gz"
cd mugqic*
module load $R_MODULE
R CMD INSTALL -l $INSTALL_DIR $PCKGS

# roxygenize the packages?


# Now write the module file
cat > $MODULEFILE <<-EOF
#%Module1.0
proc ModulesHelp { } {
        puts stderr "MUGQIC - Adds the R mugqic packages to your R_LIBS"
}
module-whatis "MUGQIC -  Adds the R mugqic packages to your R_LIB"
set root $TCLROOT
prepend-path R_LIBS \$root
EOF

# We'll have the last installed version to be default
cat > $MODULEVERSIONFILE <<-EOF
#%Module1.0
set ModulesVersion $REF
EOF


## Adjust permissions
chmod -R ug+rwX  $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
chmod -R o+rX    $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE



echo "done."

exit 0 ;



