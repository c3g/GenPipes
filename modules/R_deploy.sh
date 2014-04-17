#!/bin/bash
set -e

## This script calls R.sh on known cluster. It assumes being run on abacus, and that the user can ssh
## directly without password.
wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R.sh -O R.sh

## Abacus
sh R.sh >& abacus.R.log $@

## Mammouth MP2
ssh lefebvr3@bourque-mp2.rqchp.ca  "bash -l -s" -- >& mammouth.R.log  < R.sh $@


## Guillimin phase 1
#ssh flefebvre1@guillimin.clumeq.ca "bash -l -s" -- >& guillimin.R.log < R.sh $@

## Guillimin phase 2
ssh flefebvre1@guillimin-p2.hpc.mcgill.ca "bash -l -s" -- >& guillimin2.R.log < R.sh $@

exit

# wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_deploy.sh -O R_deploy.sh && wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R.sh -O R.sh

# sh -l R_deploy.sh -f -v 3.0.0
# sh -l R_deploy.sh -f -v 3.0.2


# sh -l R.sh -f -v 3.0.2 -p MUGQIC_INSTALL_HOME_DEV -i software/R -m modulefiles/mugqic_dev/R >& logdev
# sh -l R.sh -f -v 3.0.2 -p MUGQIC_INSTALL_HOME -i software/R -m modulefiles/mugqic/R >& logprod



