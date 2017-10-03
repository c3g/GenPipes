#!/bin/bash
# Exit immediately on error
set -ex

if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

ROOT="${!INSTALL_HOME}/genomes/mirbase/"; mkdir -p $ROOT ; cd $ROOT
wget ftp://mirbase.org/pub/mirbase/CURRENT/README -O mirbase_README # to trace version number
wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O hairpin.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz -O mature.fa.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/miFam.dat.gz -O miFam.dat.gz
wget ftp://mirbase.org/pub/mirbase/CURRENT/organisms.txt.gz -O organisms.txt.gz
gunzip -f *.gz

