##
ROOT="$MUGQIC_INSTALL_HOME/genomes/unite_db/"; mkdir -p $ROOT ; cd $ROOT

## The last version is 1211 
wget https://github.com/downloads/qiime/its-reference-otus/its_12_11_otus.tar.gz
tar -zxvf its_12_11_otus.tar.gz
mv its_12_11_otus 1211
rm -f its_12_11_otus.tar.gz

## Rename files for configuration pattern.
cd 1211/taxonomy/
gunzip 97_otu_taxonomy.txt.gz
gunzip 99_otu_taxonomy.txt.gz

cd $ROOT
cd 1211/rep_set/
gunzip 97_otus.fasta.gz
gunzip 99_otus.fasta.gz

