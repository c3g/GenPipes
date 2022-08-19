#!/bin/bash
# Exit immediately on error
set -ex

if [[ ${1:-} == MUGQIC_INSTALL_HOME ]]
then
  INSTALL_HOME=MUGQIC_INSTALL_HOME
else
  INSTALL_HOME=MUGQIC_INSTALL_HOME_DEV
fi

ROOT="${!INSTALL_HOME}/genomes/silva_db/";
mkdir -p $ROOT; 
cd $ROOT

VERSION=132 # The last version is 132
ARCHIVE=Silva_${VERSION}_release.zip
ARCHIVE_URL=https://www.arb-silva.de/fileadmin/silva_databases/qiime/${ARCHIVE}

# If archive was previously downloaded, use the local one, otherwise get it from remote site
if [[ -f $ARCHIVE ]]
then
  echo "Archive $ARCHIVE already in $ROOT/: using it..."
else
  echo "Archive $ARCHIVE not in $ROOT/: downloading it..."
  wget --no-check-certificate ${ARCHIVE_URL} --output-document=$ROOT/$ARCHIVE
fi

unzip $ARCHIVE

rm -rf __MACOSX
mv SILVA_${VERSION}_QIIME_release $VERSION
rm -f Silva_${VERSION}_provisional_release.zip

## Rename files for configuration pattern.
for i in 16S 18S
do
  cd $ROOT/${VERSION}/taxonomy/${i}_only/
  for j in 90 94 97 99
  do
    cp ${j}/taxonomy_all_levels.txt ${j}_otu_taxonomy.txt
  done

  cd $ROOT/${VERSION}/rep_set/rep_set_${i}_only/
  for j in 90 94 97 99
  do
    file=`ls $j/*`
    cp $file ${j}_otus.fasta
  done

done 
