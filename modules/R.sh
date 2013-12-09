#!/bin/bash
set -e
umask 0002
me=`basename $0`

## Neutralize $R_LIBS
export R_LIBS=

## Default arg values
VERSION="latest"
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
   -m      Where to write the module file, defaults to $MUGQIC_INSTALL_HOME_DEV/modulesfiles/mugqic_dev/R/
   -i	   R installation directory, defaults to $MUGQIC_INSTALL_HOME_DEV/software/R. Actual install dir will be reassigned to <-i>/R-version
   -h      Print this message

NOTES: 
   On Mammouth, byte compiling of the base package might fail if an R module is already loaded.
   The problem was not investiguated further except that it might be related 
   to the gcc compiler there.
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
TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR 
echo "Working in $TEMPDIR"

## If latest is requested, determine version number. Unfort., only way seem to download R-latest tar.gz!
if [  $VERSION = "latest" ]
then
	wget --no-verbose http://cran.r-project.org/src/base/R-latest.tar.gz
	tar -xf R-latest.tar.gz
	VERSION=`cat ./*/VERSION`
	rm -r R*
fi

## Paths, mkdirs
INSTALL_DIR=$INSTALL_DIR/R-$VERSION
MODULEFILE="$MODULEFILE_DIR/$VERSION"
MODULEVERSIONFILE="$MODULEFILE_DIR/.version"
mkdir -p $MODULEFILE_DIR $INSTALL_DIR

## Install if required by force or absence of module files
if  [ ! -f $MODULEFILE ] || $FORCE_INSTALL
then
	# Prelim. cleanup
	rm -rf $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
 
	# Download, compile, install
	wget --no-verbose http://cran.r-project.org/src/base/R-${VERSION:0:1}/R-$VERSION.tar.gz
	tar -xf R-$VERSION.tar.gz
	cd R-$VERSION
	sed -i 's/Sys.umask("022")/Sys\.umask("002")/g' src/library/tools/R/build.R # hack to force umask 0002
	./configure --prefix=$INSTALL_DIR  # TEMP s--with-readline=yes --with-readline=no
	make -j8
	make install
	cd $TEMPDIR
	
	# Create module files
	cat > $MODULEFILE <<-EOF	
		#%Module1.0
		proc ModulesHelp { } {
		puts stderr "MUGQIC - Adds R to your environment"
		}
		module-whatis "MUGQIC - Adds R to your environment"
		
		set             root               $INSTALL_DIR
		setenv          R_LIBS             \$root/lib64/R/library
		#prepend-path   MANPATH            \$root/share              
		prepend-path    PATH               \$root/bin
		prepend-path    LD_LIBRARY_PATH    \$root/lib64:/software/libraries/GotoBLAS_LAPACK/shared
		#prepend-path   LD_LIBRARY_PATH    \$root/lib64:\$root/standalone:/software/libraries/GotoBLAS_LAPACK/shared
		#prepend-path   CPATH              \$root/include
		EOF
	cat > $MODULEVERSIONFILE <<-EOF
		#%Module1.0
		set ModulesVersion $VERSION
		EOF
fi

## Finally, update/install library!
$INSTALL_DIR/bin/R --vanilla  <<-'EOF'

	#' This script:
	#' 1) Installs or update all packages hard-coded below
	#' 2) Installs or updates all Bioconductor org.* packages 
	#' 3) Installs one non-CRAN/Biocuctor packages (Vennerable)
	#' 4) Installs "gqUtils","gqSeqUtils","gqData","gqMicroarrays" from the *master* branch of Rpackages i.e.
	#'	 https://bitbucket.org/mugqic/rpackages/get/master.zip
	#' 5) Grants +rX permission to all and +w to group on package library (.Library)

	## biocLite
	source("http://bioconductor.org/biocLite.R")

	## Define the list of packages to standard packages to install.
	deps = c("affxparser","affy","affyio","affyPLM","akima","annotate","AnnotationDbi"
	,"AnnotationForge","ape","ash","base","beanplot","Biobase","BiocGenerics"
	,"BiocInstaller","bioDist","biomaRt","Biostrings","biovizBase","bit"
	,"bitops","boot","brew","BSgenome","caTools","charm","charmData","class"
	,"cluster","codetools","colorspace","compiler","corpcor","crlmm","ctc"
	,"cummeRbund","datasets","DBI","DESeq","devtools","dichromat","digest"
	,"edgeR","ellipse","evaluate","fastcluster","ff","fields"
	,"foreach","foreign","gcrma","gdata","genefilter","GenomicFeatures"
	,"GenomicRanges","genoset","GEOquery","ggplot2","googleVis","goseq"
	,"gplots","graph","graphics","grDevices","grid","gtable","gtools"
	,"Gviz","hdrcde","Hmisc","hwriter","IlluminaHumanMethylation450k.db"
	,"IlluminaHumanMethylation450kmanifest","impute","IRanges","iterators"
	,"KernSmooth","ks","labeling","lattice","latticeExtra","limma","locfit"
	,"lumi","LVSmiRNA","maps","markdown","MASS","Matrix","matrixStats","mclust"
	,"memoise","methods","methyAnalysis","methylumi","mgcv","minfi","misc3d"
	,"multicore","multtest","munsell","mvtnorm","NBPSeq","nleqslv","nlme"
	,"nnet","nor1mix","Nozzle.R1","oligo","oligoClasses","outliers","parallel"
	,"pd.charm.hg18.example","pheatmap","plotrix","plyr","plyr","preprocessCore"
	,"proto","quantreg","R2HTML","RBGL","RColorBrewer","Rcpp","RcppEigen","RCurl"
	,"ReportingTools","reshape","reshape2","rgl","RJSONIO","R.methodsS3","roxygen2"
	,"rpart","Rsamtools","RSQLite","rtracklayer","scales","siggenes","snow"
	,"SNPchip","SortableHTMLTables","spam","SparseM","spatial","splines","SQN"
	,"statmod","stats","stats4","stringr","survival","sva","tcltk","testthat"
	,"tools","TxDb.Hsapiens.UCSC.hg19.knownGene","utils","Vennerable","vsn"
	,"WriteXLS","XML","xtable","zlibbioc")

	## Programmatically add all the org pacakges
	contribUrl = contrib.url(biocinstallRepos(), type = 'source')
	availPkgs  = available.packages(contribUrl, type = 'source')	
	org.packages = rownames(availPkgs)[grepl("^org", rownames(availPkgs))]
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

## Adjust permissions
chmod -R ug+rwX  $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
chmod -R o+rX    $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE

exit
