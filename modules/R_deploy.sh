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