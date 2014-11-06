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

# Deploy prod
# wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_Bioconductor.sh -O R_Bioconductor.sh
# wget https://bitbucket.org/mugqic/mugqic_resources/raw/master/modules/R_mugqic_packages.sh -O R_mugqic_packages.sh
# sh R_Bioconductor.sh -r -p MUGQIC_INSTALL_HOME -m modulesfiles/mugqic -i software &> R_Bioconductor_prod.log
# sh R_mugqic_packages.sh -v 0.1 -R mugqic/R_Bioconductor/3.1.2_3.0 -p MUGQIC_INSTALL_HOME -m modulesfiles/mugqic -i software &> R_mugqic_packages_prod.log
#
#

