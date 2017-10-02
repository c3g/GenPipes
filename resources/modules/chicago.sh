#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

SOFTWARE=chicago
VERSION=1.1.5
SOFTWARE_DIR=${SOFTWARE}_$VERSION

INSTALL_HOME=$MUGQIC_INSTALL_HOME

MODULE_R=mugqic/R_Bioconductor/3.2.3_3.2

INSTALL_PATH=$INSTALL_HOME/software/$SOFTWARE/$SOFTWARE_DIR
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH

git clone git@bitbucket.org:chicagoTeam/chicago.git
mv $SOFTWARE/* .
rm -fr $SOFTWARE/

chmod 775 $INSTALL_PATH/chicagoTools/*.R
chmod 775 $INSTALL_PATH/chicagoTools/*.py

## module file
mkdir -p $INSTALL_HOME/modulefiles/mugqic/$SOFTWARE

module load $MODULE_R
R  --no-save --no-restore  <<-'EOF'
## will need to install R libraries:
  install.packages("httr", repos="http://cran.rstudio.com/")
  install.packages("devtools", repos="http://cran.rstudio.com/")
  install.packages("argparser", repos="http://cran.rstudio.com/")
  library(httr)
  set_config(config(ssl_verifypeer = 0L))
  library(devtools)
  install_github("trevorld/findpython", dependencies = TRUE)
  install.packages("rjson", repos="http://cran.rstudio.com/")
  install_github("trevorld/getopt", dependencies = TRUE)
  install_github("trevorld/argparse", dependencies = TRUE)
  install_bitbucket("aadler/delaporte")
  install_bitbucket("chicagoTeam/Chicago", subdir="Chicago", dependencies = TRUE)
EOF

#Module definition to use
echo "\
#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - $SOFTWARE \"
}
module-whatis \"Chicago R package for capture Hi-C analysis\"

set             root                $INSTALL_PATH
prepend-path    PATH                \$root
prepend-path    PATH                \$root/chicagoTools

" > $VERSION


# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mv .version $VERSION $INSTALL_HOME/modulefiles/mugqic/$SOFTWARE
