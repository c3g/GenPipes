#!/bin/bash
set -e


## This script calls R.sh on known cluster. It assumes being run on abacus, and that the user can ssh
## directly without password.
wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_Bioconductor.sh -O R_Bioconductor.sh
wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_mugqic_packages.sh -O R_mugqic_packages.sh

## Abacus
sh R_Bioconductor.sh >& abacus.R_Bioconductor.log $@
sh R_mugqic_packages.sh >& abacus.R_mugqic_packages.log $@

## Guillimin 
ssh flefebvre1@guillimin-p2.hpc.mcgill.ca "bash -l -s" -- >& guillimin.R_Bioconductor.log < R_Bioconductor.sh $@
ssh flefebvre1@guillimin-p2.hpc.mcgill.ca "bash -l -s" -- >& guillimin.R_mugqic_packages.log < R_mugqic_packages.sh $@

## Mammouth MP2
ssh lefebvr3@bourque-mp2.rqchp.ca  "bash -l -s" -- >& mammouth.R_Bioconductor.log  < R_Bioconductor.sh $@
ssh lefebvr3@bourque-mp2.rqchp.ca  "bash -l -s" -- >& mammouth.R_mugqic_packages.log  < R_mugqic_packages.sh $@



exit

# wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_deploy.sh -O R_deploy.sh && wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R.sh -O R.sh

# sh -l R_deploy.sh -f -v 3.0.0
# sh -l R_deploy.sh -f -v 3.0.2


# sh -l R.sh -f -v 3.0.2 -p MUGQIC_INSTALL_HOME_DEV -i software/R -m modulefiles/mugqic_dev/R >& logdev
# sh -l R.sh -f -v 3.1.0 -p MUGQIC_INSTALL_HOME -i software/R -m modulefiles/mugqic/R >& logprod

# sh -l R.sh -f -v 3.1.1 -i /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_dev/software/R -m /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_dev/modulefiles/mugqic_dev/R >& logdev &
#
# sh -l R.sh -f -v 3.1.1 -i /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod/software/R -m /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod/modulefiles/mugqic/R >& logprod &

