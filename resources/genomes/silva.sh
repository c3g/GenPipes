#!/bin/bash
# Exit immediately on error
set -ex

if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

ROOT="${!INSTALL_HOME}/genomes/silva_db/"; mkdir -p $ROOT ; cd $ROOT

#VERSION=119
VERSION=123 # The last version is 123
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_${VERSION}_release.zip
unzip Silva_${VERSION}_release.zip

rm -rf __MACOSX
mv SILVA${VERSION}_QIIME_release $VERSION
rm -f Silva_${VERSION}_provisional_release.zip

## Rename files for configuration pattern.
cd 123/taxonomy/16S_only/
cp 90/taxonomy_all_levels.txt 90_otu_taxonomy.txt
cp 94/taxonomy_all_levels.txt 94_otu_taxonomy.txt
cp 97/taxonomy_all_levels.txt 97_otu_taxonomy.txt
cp 99/taxonomy_all_levels.txt 99_otu_taxonomy.txt

cd $ROOT
cd 123/rep_set/rep_set_16S_only/
cp 90/90_otus_16S.fasta 90_otus.fasta
cp 94/94_otus_16S.fasta 94_otus.fasta
cp 97/97_otus_16S.fasta 97_otus.fasta
cp 99/99_otus_16S.fasta 99_otus.fasta

