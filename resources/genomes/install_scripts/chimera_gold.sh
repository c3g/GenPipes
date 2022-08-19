##
ROOT="$MUGQIC_INSTALL_HOME/genomes/chimera_gold_db/"; mkdir -p $ROOT ; cd $ROOT

## The last version is the 20110519 from Broad Microbiome Utilities. 
wget http://drive5.com/uchime/gold.fa
mkdir 20110519
mv gold.fa 20110519/gold.fasta

