
######### virus of rins for VirusFinder

cd $MUGQIC_INSTALL_HOME/genomes/blast_db/
wget http://khavarilab.stanford.edu/downloads/rins_core.tar.gz
tar -xvf rins_core.tar.gz
mv rins_core/indexes/virus.fa.gz ./
rm -r rins_core
rm rins_core.tar.gz
zcat virus.fa.gz > virus.fa
rm rins_core.tar.gz

makeblastdb -in virus.fa -dbtype nucl -out virus
rm virus.fa.gz 
