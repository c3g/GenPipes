#!/bin/bash
set -e
umask 0002
me=`basename $0`

# Constants
SOFTWARE="R_Bioconductor" # used for module and install dir names

## Neutralize $R_LIBS
export R_LIBS=

## Default arg values
R_VERSION="latest"
INSTALL_PREFIX_ENV_VARNAME=""
MODULEFILE_DIR="$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev"
INSTALL_DIR="$MUGQIC_INSTALL_HOME_DEV/software"
UPDATE_MODE=1

## Parse arguments
usage()
{
cat << EOF
usage: $0 options

This script installs R, installs the corresponding module and other dependencies.

OPTIONS:
   -v      Force a specifc R version different than latest, e.g. 3.1.1. Note that the Bioconductor version installed will always be the latest one for this R version.
   -p      Name of an environnment variable which defines a prefix path to -m and -i. E.g. MUGQIC_INSTALL_HOME
   -m      The root folder for the module files, defaults to \$MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev. 
   -h      Print this message

EXAMPLE USAGE:

R_Bioconductor.sh -v 3.1.1 -p MUGQIC_INSTALL_HOME     -m modulefiles/mugqic     -i software
R_Bioconductor.sh -v 3.1.1 -p MUGQIC_INSTALL_HOME_DEV -m modulefiles/mugqic_dev -i software

NOTES: 
   - On Mammouth, byte compiling of the base package might fail if an R module is already loaded.
   The problem was not investiguated further except that it might be related to the gcc compiler there.
   - Frequently missing libraries: libcurl-devel, cairo...
   - ()...)
EOF
}

while getopts "v:b:p:m:i:rh" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         v)
             R_VERSION=$OPTARG
             ;;
         b)
             BIOCVERSION=$OPTARG
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
         r)
             UPDATE_MODE=0
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

#module load mugqic/openssl/1.1.0 mugqic/autoconf/2.69

## Tmp dir to work in
TEMPDIR=`mktemp -d -t $me.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX` && cd $TEMPDIR
echo "Working in $TEMPDIR"

## If latest is requested, determine version number. Unfort., only way seem to download R-latest tar.gz!
if [[  $R_VERSION == "latest" ]]
then
        wget --no-verbose https://cran.r-project.org/src/base/R-latest.tar.gz
        tar -xf R-latest.tar.gz
        R_VERSION=`cat ./*/VERSION`
        rm -r R*
        echo "Latest R version appears to be $R_VERSION"
        BIOCVERSION=`wget -qO- https://bioconductor.org/packages/release/bioc/ | grep 'Bioconductor version: Release ' | grep -oE '[0-9]*\.[0-9]*'`
        echo "Latest Bioconductor version appears to be $BIOCVERSION"
        VERSION="$R_VERSION""_""$BIOCVERSION"
elif [[ $BIOCVERSION != "" ]]
then
        VERSION="$R_VERSION""_""$BIOCVERSION"
else
        VERSION=$R_VERSION
fi

## Paths, mkdirs
INSTALL_DIR="$INSTALL_DIR/$SOFTWARE/$SOFTWARE-$VERSION"
MODULEFILE_DIR="$MODULEFILE_DIR/$SOFTWARE"
TCLROOT=$INSTALL_DIR
if [[ $INSTALL_PREFIX_ENV_VARNAME != "" ]]
then
        echo "prefixing..."
        TCLROOT="\$::env($INSTALL_PREFIX_ENV_VARNAME)/$INSTALL_DIR"
        INSTALL_DIR=${!INSTALL_PREFIX_ENV_VARNAME}/$INSTALL_DIR
        MODULEFILE_DIR=${!INSTALL_PREFIX_ENV_VARNAME}/$MODULEFILE_DIR
fi
MODULEFILE="$MODULEFILE_DIR/$VERSION"
MODULEVERSIONFILE="$MODULEFILE_DIR/.version"
# NOTE: this is somewhat complicated because we want the ROOT dir MUGQIC_INSTALL_HOME to be resolved at module execution.
# TCLROOT is just a variable holding the TCL script value for the 'root' variable in the module file.

if [ `lsb_release -i | cut -f 2` == "Ubuntu" ]
then
  echo "Ubuntu" > /dev/null
  C3G_SYSTEM_LIBRARY=/cvmfs/soft.mugqic/apt/ubuntu1604/1.0
  LIB=lib
  INTERPRETER=$C3G_SYSTEM_LIBRARY/$LIB/x86_64-linux-gnu/ld-linux-x86-64.so.2
  LIBDIR=$C3G_SYSTEM_LIBRARY/$LIB/x86_64-linux-gnu:$C3G_SYSTEM_LIBRARY/usr/$LIB/x86_64-linux-gnu:$C3G_SYSTEM_LIBRARY/$LIB:$C3G_SYSTEM_LIBRARY/usr/$LIB:$INSTALL_DIR/$LIB/R/lib
elif [ `lsb_release -i | cut -f 2` == "CentOS" ]
then
  echo "CentOS" > /dev/null
  C3G_SYSTEM_LIBRARY=/cvmfs/soft.mugqic/yum/centos7/1.0
  LIB=lib64
  INTERPRETER=$C3G_SYSTEM_LIBRARY/$LIB/ld-linux-x86-64.so.2
  LIBDIR=$INSTALL_DIR/$LIB/R/lib:$INSTALL_DIR/$LIB/R/library:$C3G_SYSTEM_LIBRARY/usr/local/c3g/rpm/usr/lib64:$C3G_SYSTEM_LIBRARY/usr/local/c3g/compile/lib:$C3G_SYSTEM_LIBRARY/usr/local/lib64:$C3G_SYSTEM_LIBRARY/usr/lib64:$C3G_SYSTEM_LIBRARY/usr/lib:$C3G_SYSTEM_LIBRARY/usr/lib64/mysql:$C3G_SYSTEM_LIBRARY/$LIB/mysql
else
  echo "*** ERROR ***"
  echo "'"`lsb_release -i | cut -f 2`"' OS detected... should be either 'Ubuntu' neither 'CentOS'..."
  exit 1
fi

echo "The software install location is $INSTALL_DIR"
echo "The module file directory is $MODULEFILE_DIR"
echo "The module file is $MODULEFILE"
echo "The module version file is $MODULEVERSIONFILE"

## Dir creation
mkdir -p $MODULEFILE_DIR $INSTALL_DIR

## Install if required by force or absence of module files
if [[ ( ! ( -f "$MODULEFILE" && -f "$INSTALL_DIR/bin/R") ) || "$UPDATE_MODE" = 0   ]]
then
        echo "Installing base R..."

        # Prelim. cleanup
        rm -rf $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE

        # Download, compile, install
        wget --no-verbose https://cran.r-project.org/src/base/R-${R_VERSION:0:1}/R-$R_VERSION.tar.gz
        tar -xf R-$R_VERSION.tar.gz
        cd R-$R_VERSION

        # hack to force umask 0002: in a large group, want pkgs installed with write perm.
        sed -i 's/Sys.umask("022")/Sys\.umask("002")/g' src/library/tools/R/build.R

        ./configure --enable-R-shlib --prefix=$INSTALL_DIR  # --enable-R-shlib  is for Rpy
        make -j12
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
setenv          R_LIBS             \$root/$LIB/R/library           
prepend-path    PATH               \$root/bin
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
# This can be set in R/lib64(or lib)/etc/Rprofile.site, but then R should not be invked with R --vanilla
# because vanilla will ignore Rprofile.site.

cat > $INSTALL_DIR/$LIB/R/etc/Rprofile.site <<-'EOF'
if(capabilities()["cairo"]){ options(bitmapType="cairo") }
Sys.setenv(PAGER="/usr/bin/less")
Sys.umask("002")
EOF

fi

## Finally, update/install library!
$INSTALL_DIR/bin/R  --no-save --no-restore  <<-'EOF'

    #' This script:
    #' 1) Installs or update all packages hard-coded below
    #' 2) Installs or updates all Bioconductor org.* packages 
    #' 3) Installs one non-CRAN/Biocuctor packages (Vennerable)
    #' 4) Installs "gqUtils", "gqSeqUtils", "gqData", "gqMicroarrays" from the *master* branch of Rpackages i.e.
    #'       https://bitbucket.org/mugqic/rpackages/get/master.zip
    #' 5) Grants +rX permission to all and +w to group on package library (.Library)

    ## Install library path
    .libPaths(.Library) # useful because e.g. devtools::install() installs in .libPaths()[1], and the latter will be ~/R/... if user library exists...

    ## Install udunits2 as it is needed as a dependency for some of the below packages
    install.packages("udunits2", configure.args='--with-udunits2-include=/usr/include/udunits2/', repos='http://cran.rstudio.org', lib=.Library)

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos='http://cran.rstudio.org', lib=.Library)

    ## Define which installer to iuse : biocLite (for 'old' versions) vs. BiocManager (for newer versions)
    if (version$major < 3 && version$minor < 5.0) {
        source("https://bioconductor.org/biocLite.R")
        INSTALLER=biocLite
        REPO=biocinstallRepos
    } else {
        INSTALLER=BiocManager::install
        REPO=BiocManager::repositories
    }

    ## RcppArmadillo temporary patch: CentOS 6 is old and using an old gcc (4.4.7) which is incompatible with the latest armadillo libs. The easiest workaround seems to force install
    # of an archived version of RcppArmadillo. Note that update attempts below will fail with ERROR.
    # https://github.com/RcppCore/RcppArmadillo/issues/30
    # http://stackoverflow.com/questions/27296522/rcpparmadillo-failing-to-install-on-centos
    INSTALLER("RcppArmadillo", ask=FALSE) # despite failing, this will install Rcpparmadillo dependencies for us
    rcpp.armadillo.archive="RcppArmadillo_0.4.500.0.tar.gz"
    download.file(sprintf("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/%s", rcpp.armadillo.archive), destfile=rcpp.armadillo.archive)
    install.packages(rcpp.armadillo.archive, repos=NULL, type="source", lib=.Library)
    install.packages("slam", repos="http://cran.us.r-project.org/", lib=.Library)
    # Since R3.5 'udunits2' needs to be installed so that 'clusterProfiler' and its dependencies get installed properly
    # Also, 'udunits2' installation requires hard-coding the 'configure.args' parameter
    install.packages('udunits2', repos="http://cran.us.r-project.org/", configure.args='--with-udunits2-include=/usr/include/udunits2', lib=.Library)

    ## Define the list of packages to standard packages to install.
    deps = c("affxparser", "affy", "affyio", "affyPLM", "akima", "allgown", "annotate", "AnnotationDbi", "AnnotationForge", "annotatr", "ape", "ash", "ASCAT",
    "ballgown", "BatchExperiments", "BatchJobs", "batchtools", "beanplot", "Biobase", "BiocGenerics", "BiocInstaller", "bioDist", "biomaRt", "biomformat", "Biostrings" "biovizBase", "bit",
    "biovizBase", "bit", "bit64", "bitops", "boot", "brew", "BSgenome", "BSgenome.Hsapiens.1000genomes.hs37d5", "BSgenome.Hsapiens.UCSC.hg19", "bumphunter",
    "Cairo", "caTools", "charm", "charmData", "ChIPseeker", "circlize", "class", "cluster", "clusterStab", "clusterProfiler", "cn.mops", "codetools", "colorspace", "ConsensusClusterPlus",
    "ConsensusClusterPlus", "corpcor", "crlmm", "ctc", "cummeRbund",
    "dada2", "datasets", "data.table", "DBI", "deconstructSigs", "DESeq", "devtools", "dendextend", "dichromat", "diffHic", "digest", "dmrseq", "DNAcopy", "dplyr", "DSS",
    "edgeR", "ellipse", "evaluate",
    "farver", "fastcluster", "fdrtool", "ff", "fields", "FDb.InfiniumMethylation.hg19", "foreach", "foreign",
    "gcrma", "gdata", "genefilter", "genomation", "GenomicFeatures", "GenomicRanges", "GenomeInfoDb", "genoset", "GEOquery", "ggplot2", "ggvis", "googleVis", "goseq",
    "gplots", "graph", "gsalib", "gtable", "gtools", "Gviz",
    "hdrcde", "Hmisc", "hwriter", "HTqPCR", "HTSFilter", "hopach",
    "igraph", "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "impute", "InteractionSet", "IRanges", "iterators",
    "KernSmooth", "ks",
    "labeling", "lattice", "latticeExtra", "limma", "lme4", "locfit", "lumi", "LVSmiRNA",
    "magrittr", "maps", "markdown", "MASS", "MAST", "Matrix", "matrixStats", "mclust", "memoise", "methyAnalysis", "methylKit", "methylumi", "mgcv", "minfi",
    "mirbase.db", "misc3d", "monocle", "multtest", "munsell", "mvtnorm",
    "NBPSeq", "nleqslv", "nlme", "NMF", "nnet", "nondetects", "nor1mix", "Nozzle.R1",
    "oligo", "oligoClasses", "optparse", "outliers",
    "pd.charm.hg18.example", "pheatmap", "phyloseq", "plotrix", "plyr", "preprocessCore", "proto",
    "qqman", "quantreg",
    "R2HTML", "RBGL", "RColorBrewer", "Rcpp", "RcppEigen", "RCurl", "remotes", "rhdf5", "ReportingTools", "reshape", "reshape2", "rgl", "RJSONIO", "Rmisc", "R.methodsS3",
    "rmarkdown", "roxygen2", "rpart", "Rsamtools", "RSQLite", "rtracklayer", "Rtsne",
    "scales", "sendmailR", "sequenza", "Seurat", "shiny", "ShortRead", "siggenes", "slam", "snow", "SNPchip", "SomaticCancerAlterations", "SomaticSignatures",
    "SortableHTMLTables", "spam", "SparseM", "spatial", "spp", "SQN", "statmod", "stringi", "stringr", "survival", "sva",
    "testthat", "tidyr", "tidyverse", "tweenr", "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "vioplot", "vsn",
    "WriteXLS",
    "XML", "xtable",
    "zlibbioc")

    ## Programmatically add all the org packages (excluding MeSH mess which takes too long)
    contribUrl = contrib.url(REPO(), type = 'source')
    availPkgs  = available.packages(contribUrl, type = 'source')    
    org.packages = rownames(availPkgs)[grepl("^org", rownames(availPkgs))]
    org.packages = org.packages[!grepl("^org.MeSH.", org.packages)]
    deps = c(deps, org.packages)

    ## Install pkgs not already installed, with ask=FALSE INSTALLER() takes care of updating if necessary
    INSTALLER(ask=FALSE)
    deps = setdiff(deps, rownames(installed.packages())) # Define packages that need actual install
    INSTALLER(deps, lib=.Library, ask=FALSE)
    deps = setdiff(deps, rownames(installed.packages()))
    INSTALLER(deps, lib=.Library, ask=FALSE) # twice, just to make sure
    deps = setdiff(deps, rownames(installed.packages()))
    INSTALLER(deps, lib=.Library, ask=FALSE) # and why not a third time !

    ## Install Vennerable, since not yet in CRAN
    install.packages("Vennerable", repos="http://R-Forge.R-project.org", lib=.Library)
    ## Force Rmarkdown and knitr, not available fot R 3.2
    install.packages('knitr', repos='http://cran.rstudio.org', lib=.Library)
    install.packages('rmarkdown', repos='http://cran.rstudio.org', lib=.Library)

    ## Sleuth : needs devtools and remotes to be installed (both done above) 
    INSTALLER("pachterlab/sleuth", lib=.Library, ask=FALSE)

    require(devtools)
    ## PopSV
    devtools::install_github("jmonlong/PopSV")
    ## ASCAT
    devtools::install_github("Crick-CancerGenomics/ascat/ASCAT")
    ## ChIAnalysis (with its eric.utils dependency)
    #devtools::install_bitbucket("ericfournier2/sb_lab/eric.utils")
    #devtools::install_github("ArnaudDroitLab/ChIAnalysis")
    ## DoubltFinder
    devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
    ## binless
    devtools::install_github("3DGenomes/binless",subdir="binless")

    ## Print the list of the installed packages along with their version into a file
    library(data.table)
    write.table(data.table(installed.packages())[,c(1,3)], paste(c(Sys.getenv('R_HOME')),"/../../installed.packages.txt", sep=""), sep="\t", row.names=F, quote=F)
EOF

echo "Patching C3G executables..."
for i in `find $INSTALL_DIR/ -type f -executable -exec file {} \; | grep ELF | cut -d":" -f1`; do
  if readelf -l $i | grep go.build > /dev/null
  then
    echo "GO Done" > /dev/null
  elif [ ${i##*.} == "so" ] || [[ ${i##*/} =~ "so"*(\.[0-9]{1,2})*$ ]]
  then
    $MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --set-rpath $($MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --print-rpath $i):$LIBDIR $i
  else
    echo $i
    $MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --set-interpreter $INTERPRETER --set-rpath $($MUGQIC_INSTALL_HOME/software/patchelf/patchelf-0.9/bin/patchelf --print-rpath $i):$LIBDIR $i
  fi
done

echo "Building C3G wrappers for executables..."
sed -i "s,R_binary=\"\${R_HOME}/bin/exec\${R_ARCH}/R,R_binary=\"\${R_HOME}/bin/exec\${R_ARCH}.wrap/R," $INSTALL_DIR/bin/R
sed -i "s,R_binary=\"\${R_HOME}/bin/exec\${R_ARCH}/R,R_binary=\"\${R_HOME}/bin/exec\${R_ARCH}.wrap/R," $INSTALL_DIR/$LIB/R/bin/R
mkdir $INSTALL_DIR/$LIB/R/bin/exec.wrap
echo "$INTERPRETER --library-path $LIBDIR $INSTALL_DIR/$LIB/R/bin/exec/R \${args} \${@}" > $INSTALL_DIR/$LIB/R/bin/exec.wrap/R
chmod 775 $INSTALL_DIR/$LIB/R/bin/exec.wrap/R
for i in $INSTALL_DIR/bin/Rscript $INSTALL_DIR/$LIB/R/bin/Rscript; do
  mv $i $i.raw;
  echo "$INTERPRETER --library-path $LIBDIR $i.raw \${@}" > $i;
  chmod 775 $i;
done

echo "Adjusting permissions..."
## Adjust permissions
chmod -R ug+rwX  $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE
chmod -R o+rX    $INSTALL_DIR $MODULEFILE $MODULEVERSIONFILE

echo "R packages installation done."
exit 0 ;

# ### Blurb to test graphics
# # module load mugqic/R/3.0.2
# # R --no-save <<'EOF'
#  capabilities()
#  getOption("bitmapType") # returns default bitmap device: abacus=cairo, guil=cairo, mp2=Xlib
#  jpeg("temp.jpeg")#, type='cairo')
#  plot(1)
#  dev.off()
# # EOF


## Install homebrew packages to their own separate module
# # Possible directly from devtools?
# Define tag version variablw

# #     require(roxygen2)
#       require(devtools)
#       deps = c("gqUtils", "gqSeqUtils", "gqData", "gqMicroarrays")
#       download.file("https://bitbucket.org/mugqic/rpackages/get/master.zip", destfile='.packages.zip', method='wget')
#       unzip(".packages.zip", exdir='.packages')
#       deps = file.path( list.files(".packages", full.names=TRUE), deps )
# #     sapply(deps, roxygenize) # msg sent to R-help; roxygen2 not available R 3.0.0 or 3.0.1 !!!
#       install_local(deps)
#       unlink(c(".packages.zip", ".packages"), recursive=TRUE)

