#!/bin/bash
set -e
umask 0002
me=`basename $0`


## Annoying cluster-specific actions
# Guillimin phase 1 gcc is outdated
GCC_MODULE_CALL=""
if [[ `hostname` == lg-* ]] && [[ ! $(cat /proc/cpuinfo | grep -Ec 'E5-2650|E5-2670|E5-4620') -gt 0 ]]  ;then # Guillimin phase1
	GCC_MODULE_CALL="module load gcc/4.7.2" # otherwise Rarmadillo will not install
	$GCC_MODULE_CALL
	#echo $GCC_MODULE_CALL
fi


## Neutralize $R_LIBS
export R_LIBS=

## Default arg values
VERSION="latest" 
INSTALL_PREFIX_ENV_VARNAME=""
MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/R"
INSTALL_DIR="$MUGQIC_INSTALL_HOME_DEV/software/R"
FORCE_INSTALL=false

## Parse arguments
usage()
{
cat << EOF
usage: $0 options

This script installs R, installs the corresponding module and other dependencies.

OPTIONS:
   -v      Force a specifc R version different than latest, e.g. 3.0.0
   -f      Force re-install, i.e. overwrite module and R installation even if module already exists
   -p	   Name of an environnment variable which defines a prefix path to -m and -i. E.g. MUGQIC_INSTALL_HOME
   -m      Where to write the module file, defaults to modulesfiles/mugqic_dev/R/
   -i	   R installation directory, defaults to software/R. Actual install dir will be reassigned to <-i>/R-version
   -h      Print this message

EXAMPLE USAGE:

R.sh -f -v 3.0.2 -p MUGQIC_INSTALL_HOME -m modulefiles/mugqic/R -i software/R

NOTES: 
   - On Mammouth, byte compiling of the base package might fail if an R module is already loaded.
   The problem was not investiguated further except that it might be related to the gcc compiler there.
   - Frequently missing libraries: libcurl-devel, cairo...
   - ()...)
EOF
}

while getopts “v:p:m:i:fh” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         v)
             VERSION=$OPTARG
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
TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR 
echo "Working in $TEMPDIR"

## If latest is requested, determine version number. Unfort., only way seem to download R-latest tar.gz!
if [[  $VERSION == "latest" ]]
then
	wget --no-verbose http://cran.r-project.org/src/base/R-latest.tar.gz
	tar -xf R-latest.tar.gz
	VERSION=`cat ./*/VERSION`
	rm -r R*
fi

## Paths, mkdirs
INSTALL_DIR="$INSTALL_DIR/R-$VERSION"
MODULEFILE="$MODULEFILE_DIR/$VERSION"
MODULEVERSIONFILE="$MODULEFILE_DIR/.version"
TCLROOT=$INSTALL_DIR
if [[ $INSTALL_PREFIX_ENV_VARNAME != "" ]]
then
	echo "prefixing..."
	TCLROOT="\$::env($INSTALL_PREFIX_ENV_VARNAME)/$INSTALL_DIR"
	INSTALL_DIR=${!INSTALL_PREFIX_ENV_VARNAME}/$INSTALL_DIR
	MODULEFILE=${!INSTALL_PREFIX_ENV_VARNAME}/$MODULEFILE
	MODULEVERSIONFILE=${!INSTALL_PREFIX_ENV_VARNAME}/$MODULEVERSIONFILE
fi
mkdir -p $MODULEFILE_DIR $INSTALL_DIR
# NOTE: this is somewhat complicated because we want the ROOT dir MUGQIC_INSTALL_HOME to be resolved at module execution.
# TCLROOT is just a variable holding the TCL script value for the 'root' variable in the module file.

## Install if required by force or absence of module files
if  [ ! -f $MODULEFILE ] || $FORCE_INSTALL
then
	# Prelim. cleanup
	rm -rf $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
 
	# Download, compile, install
	wget --no-verbose http://cran.r-project.org/src/base/R-${VERSION:0:1}/R-$VERSION.tar.gz
	# wget --no-verbose http://probability.ca/cran/src/base/R-${VERSION:0:1}/R-$VERSION.tar.gz
	#wget https://dl.dropboxusercontent.com/u/2528754/R-$VERSION.tar.gz # TODO TEMP
	tar -xf R-$VERSION.tar.gz
	cd R-$VERSION
	
	# hack to force umask 0002: in a large group, want pkgs installed with write perm.
	sed -i 's/Sys.umask("022")/Sys\.umask("002")/g' src/library/tools/R/build.R 

	./configure --prefix=$INSTALL_DIR  # TEMP s--with-readline=yes --with-readline=no
	make -j8
	make install
	cd $TEMPDIR
	
	# Create mechoodule files
	cat > $MODULEFILE <<-EOF	
		#%Module1.0
		proc ModulesHelp { } {
		puts stderr "MUGQIC - Adds R to your environment"
		}
		module-whatis "MUGQIC - Adds R to your environment"	
		set             root               $TCLROOT
		setenv          R_LIBS             \$root/lib64/R/library
		#prepend-path   MANPATH            \$root/share              
		prepend-path    PATH               \$root/bin
		prepend-path    LD_LIBRARY_PATH    \$root/lib64:/software/libraries/GotoBLAS_LAPACK/shared
		#prepend-path   LD_LIBRARY_PATH    \$root/lib64:\$root/standalone:/software/libraries/GotoBLAS_LAPACK/shared
		#prepend-path   CPATH              \$root/include
		$GCC_MODULE_CALL
	EOF
	cat > $MODULEVERSIONFILE <<-EOF
		#%Module1.0
		set ModulesVersion $VERSION
		EOF
		
		## Define Rprofile.site: 
		# On headless nodes, can only use cairo backend OR have Xvb running
		# so we can have cairo-devel.x86_64 installed by admins, and then re-compile R. It should
		# then have cairo support. There is a weird behavior though. Accoring to the R-doc, cairo
		# is supposed to be set as the default backend -> getOption("bitmapType") on linux instead
		# of Xlib if cairo is available.
		# However it is not always the case. Not sure if this is a bug, might have something to do with pango present or not:
		# http://r.789695.n4.nabble.com/cairo-is-not-the-default-when-available-td4666691.html
		# http://r.789695.n4.nabble.com/Xll-options-td3725879.html
		#  so the only other way is to set options(bitmapType="cairo")
		# This can be set in R/lib64/etc/Rprofile.site, but then R should not be invked with R --vanilla
		# because vanilla will ignore Rprofile.site.
		cat > $INSTALL_DIR/lib64/R/etc/Rprofile.site <<-EOF
			if(capabilities()["cairo"]){ options(bitmapType="cairo") }
			Sys.umask("002")
		EOF
fi


## Finally, update/install library!
$INSTALL_DIR/bin/R  --no-save --no-restore  <<-'EOF'

	#' This script:
	#' 1) Installs or update all packages hard-coded below
	#' 2) Installs or updates all Bioconductor org.* packages 
	#' 3) Installs one non-CRAN/Biocuctor packages (Vennerable)
	#' 4) Installs "gqUtils","gqSeqUtils","gqData","gqMicroarrays" from the *master* branch of Rpackages i.e.
	#'	 https://bitbucket.org/mugqic/rpackages/get/master.zip
	#' 5) Grants +rX permission to all and +w to group on package library (.Library)

	## Install library path
	.libPaths(.Library) # useful because e.g. devtools::install() installs in .libPaths()[1], and the latter will be ~/R/... if user library exists...

	## biocLite
	source("http://bioconductor.org/biocLite.R")

	## Define the list of packages to standard packages to install.
	deps = c("affxparser","affy","affyio","affyPLM","akima","annotate","AnnotationDbi"
	,"AnnotationForge","ape","ash","base","BatchExperiments","BatchJobs","beanplot","Biobase","BiocGenerics"
	,"BiocInstaller","bioDist","biomaRt","Biostrings","biovizBase","bit"
	,"bitops","boot","brew","BSgenome","caTools","charm","charmData","class"
	,"cluster","clusterProfiler","codetools","colorspace","compiler","corpcor","crlmm","ctc"
	,"cummeRbund","datasets","DBI","DESeq","devtools","dichromat","digest"
	,"edgeR","ellipse","evaluate","fastcluster","ff","fields"
	,"foreach","foreign","gcrma","gdata","genefilter","GenomicFeatures"
	,"GenomicRanges","genoset","GEOquery","ggplot2","googleVis","goseq"
	,"gplots","graph","graphics","grDevices","grid","gtable","gtools"
	,"Gviz","hdrcde","Hmisc","hwriter","HTSFilter","igraph","IlluminaHumanMethylation450k.db"
	,"IlluminaHumanMethylation450kmanifest","impute","IRanges","iterators"
	,"KernSmooth","ks","labeling","lattice","latticeExtra","limma","locfit"
	,"lumi","LVSmiRNA","maps","markdown","MASS","Matrix","matrixStats","mclust"
	,"memoise","methods","methyAnalysis","methylumi","mgcv","minfi","mirbase.db","misc3d"
	,"multicore","multtest","munsell","mvtnorm","NBPSeq","nleqslv","nlme"
	,"nnet","nor1mix","Nozzle.R1","oligo","oligoClasses","outliers","parallel"
	,"pd.charm.hg18.example","pheatmap","plotrix","plyr","plyr","preprocessCore"
	,"proto","quantreg","R2HTML","RBGL","RColorBrewer","Rcpp","RcppEigen","RCurl"
	,"ReportingTools","reshape","reshape2","rgl","RJSONIO","R.methodsS3","roxygen2"
	,"rpart","Rsamtools","RSQLite","rtracklayer","scales","sendmailR","ShortRead","siggenes","snow"
	,"SNPchip","SortableHTMLTables","spam","SparseM","spatial","splines","SQN"
	,"statmod","stats","stats4","stringr","survival","sva","tcltk","testthat"
	,"tools","TxDb.Hsapiens.UCSC.hg19.knownGene","utils","Vennerable","vioplot","vsn"
	,"WriteXLS","XML","xtable","zlibbioc")

	## Programmatically add all the org packages (excluding MeSH mess which takes too long)
	contribUrl = contrib.url(biocinstallRepos(), type = 'source')
	availPkgs  = available.packages(contribUrl, type = 'source')	
	org.packages = rownames(availPkgs)[grepl("^org", rownames(availPkgs))]
	org.packages = org.packages[!grepl("^org.MeSH.",org.packages)]
	deps = c(deps,org.packages)

	## Install pkgs not already installed, with ask=FALSE biocLite() takes care of updating if necessary
	biocLite(ask=FALSE)
	deps = setdiff(deps,rownames(installed.packages())) # Define packages that need actual install
	biocLite(deps,lib=.Library,ask=FALSE)
	deps = setdiff(deps,rownames(installed.packages()))
	biocLite(deps,lib=.Library,ask=FALSE) # twice, just to make sure

	## Install Vennerable, since not yet in CRAN
	install.packages("Vennerable", repos="http://R-Forge.R-project.org",lib=.Library, type='source')

	## Install homebrew packages
#	require(roxygen2)
	require(devtools)
	deps = c("gqUtils","gqSeqUtils","gqData","gqMicroarrays")
	download.file("https://bitbucket.org/mugqic/rpackages/get/master.zip",destfile='.packages.zip',method='wget')
	unzip(".packages.zip",exdir='.packages')
	deps = file.path( list.files(".packages",full.names=TRUE), deps )
#	sapply(deps,roxygenize) # msg sent to R-help; roxygen2 not available R 3.0.0 or 3.0.1 !!!
	install_local(deps)
	unlink(c(".packages.zip",".packages"),recursive=TRUE)

	## chmod
	#system(paste("chmod -R ug+rwX",  .Library))
	#system(paste("chmod -R o+rX", .Library))
	EOF


echo "R packages installation done."

## Adjust permissions
chmod -R ug+rwX  $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
chmod -R o+rX    $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE

exit






# ## Test
# me="R.sh"
# export R_LIBS=
# VERSION="3.0.1"
# INSTALL_PREFIX_ENV_VARNAME="MUGQIC_INSTALL_HOME_DEV"
# MODULEFILE_DIR="modulefiles/mugqic_dev/R"
# INSTALL_DIR="software/R"
# #MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/R"
# #INSTALL_DIR="$MUGQIC_INSTALL_HOME_DEV/software/R"
# FORCE_INSTALL=true
# 

# ### Blurb to test graphics
# # module load mugqic/R/3.0.2
# # R --no-save <<'EOF'
#  capabilities()
#  getOption("bitmapType") # returns default bitmap device: abacus=cairo,guil=cairo,mp2=Xlib
#  jpeg("temp.jpeg")#,type='cairo')
#  plot(1)
#  dev.off()
# # EOF

