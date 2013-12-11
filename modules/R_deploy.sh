#!/bin/bash
set -e

## This script calls R.sh on known cluster. It assumes being run on abacus, and that the user can ssh
## directly without password.
wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R.sh -O R.sh

## Abacus
sh R.sh >& abacus.R.log $@

## Guillimin
ssh flefebvre1@guillimin.clumeq.ca "bash -l -s" -- >& guillimin.R.log < R.sh $@

## Mammouth MP2
ssh lefebvr3@bourque-mp2.rqchp.ca  "bash -l -s" -- >& mammouth.R.log  < R.sh $@

exit

# wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_deploy.sh -O R_deploy.sh && wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R.sh -O R.sh

# sh R_deploy.sh -f -v 3.0.0
# sh R_deploy.sh -f -v 3.0.2

# sh R.sh -f -v 3.0.0
# sh R.sh -f -v 3.0.2
# sh R.sh -f -m $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R -i $MUGQIC_INSTALL_HOME/software/R -v 3.0.0
# sh R.sh -f -m $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R -i $MUGQIC_INSTALL_HOME/software/R -v 3.0.2

# ls -alrth $MUGQIC_INSTALL_HOME_DEV/software/R/R-3.0.0/lib64/R/library/
# ls -alrth $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/R/3.0.0
# 
# ls -alrth $MUGQIC_INSTALL_HOME_DEV/software/R/R-3.0.2/lib64/R/library/
# ls -alrth $MUGQIC_INSTALL_HOME_DEV/modulefiles/mugqic_dev/R/3.0.2
# 
# ls -alrth  $MUGQIC_INSTALL_HOME/software/R/R-3.0.0/lib64/R/library/
# ls -alrth  $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R/3.0.0
# 
# ls -alrth  $MUGQIC_INSTALL_HOME/software/R/R-3.0.2/lib64/R/library/
# ls -alrth  $MUGQIC_INSTALL_HOME/modulefiles/mugqic/R/3.0.2



