##
######### galGal3

module load mugqic/picard
module load mugqic/java/oracle-jdk1.6.0_38

## Fetch from iGenomes
SPECIE="Gallus_gallus"
GENOME="galGal3"
mkdir -p $MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}/fasta
mkdir -p $MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}/fasta/byChr
mkdir -p $MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}/annotations
cd $MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}/fasta/byChr
wget "http://hgdownload.soe.ucsc.edu/goldenPath/galGal3/bigZips/chromFa.tar.gz"
tar xzvf chromFa.tar.gz
rm chromFa.tar.gz
cat chr[0-9].fa chr[1-9][0-9].fa chr[WZM].fa chrE64.fa chr*_*.fa > ../${GENOME}.fasta

# Reference dictionary (not done for iGenomes mm10)
cd $MUGQIC_INSTALL_HOME/genomes/${SPECIE}/${GENOME}/fasta
java -Xmx1G -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=${GENOME}.fasta  OUTPUT=${GENOME}.fasta.dict GENOME_ASSEMBLY=${GENOME} SPECIES=${SPECIE}

