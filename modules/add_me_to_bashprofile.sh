#!/bin/bash
#
#
#
# This script runs operations common to all installtions, such as setting the mask and defining $MUGQIC_INSTALL_HOME. It should also be added to the user's bash profile

# Set umask
umask 0002

# Guess cluster and set $MUGQIC_INSTALL_HOME accordingly
HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;
if [[ $HOST == abacus* ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $HOST == lg-* ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $HOST == ip03 ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
elif [[ $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
 export MUGQIC_INSTALL_HOME=/sb/programs/analyste
elif [[ $DNSDOMAIN == guillimin.clumeq.ca ]]; then
 export MUGQIC_INSTALL_HOME=/software/areas/genomics
elif [[ $DNSDOMAIN == m ]]; then
 export MUGQIC_INSTALL_HOME=/mnt/lustre03/bourque/bourque_group/opt
fi

# Module use
module use $MUGQIC_INSTALL_HOME/modulefiles


