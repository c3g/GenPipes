




# Usage: Rscript update_R_and_Bioconductor.r
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
# Fetch  list of standard packages from bitbucket
# download.file("https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_and_Bioconductor_packages.txt"
# 	,destfile='.packages',method='wget')
# deps=readLines(".packages")
# file.remove(".packages")
# writeLines(paste(sort(deps),collapse="\",\""))

## Programmatically add all the org pacakges
contribUrl = contrib.url(biocinstallRepos(), type = 'source')
availPkgs  = available.packages(contribUrl, type = 'source')	
org.packages = rownames(availPkgs)[grepl("^org", rownames(availPkgs))]
deps = c(deps,org.packages)

## Define packages that need actual install
deps = setdiff(deps,rownames(installed.packages()))

## Install pkgs not already installed, with ask=FALSE biocLite() takes care of updating if necessary
biocLite(deps,lib=.Library,ask=FALSE)
biocLite(deps,lib=.Library,ask=FALSE) # twice, just to make sure

## Install Vennerable, since not yet in CRAN
install.packages("Vennerable", repos="http://R-Forge.R-project.org",lib=.Library, type='source')

## Install homebrew packages
require(roxygen2)
require(devtools)
deps = c("gqUtils","gqSeqUtils","gqData","gqMicroarrays")
download.file("https://bitbucket.org/mugqic/rpackages/get/master.zip",destfile='.packages.zip',method='wget')
unzip(".packages.zip",exdir='.packages')
deps = file.path( list.files(".packages",full.names=TRUE), deps )
sapply(deps,roxygenize)
install_local(deps)
unlink(c(".packages.zip",".packages"),recursive=TRUE)

## chmod
system(paste("chmod -R a+rX",  .Library))
system(paste("chmod -R g+w", .Library))

