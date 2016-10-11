##
ROOT="$MUGQIC_INSTALL_HOME/genomes/chimera_unite_db/"; mkdir -p $ROOT ; cd $ROOT

## The 20150311 version
#mkdir 20150311
#mkdir tmp/; cd tmp/
#wget https://unite.ut.ee/sh_files/uchime_reference_dataset_11.03.2015.zip
#unzip uchime_reference_dataset_11.03.2015.zip
#mv uchime_sh_refs_dynamic_original_985_11.03.2015.fasta ../20150311/unite.fasta
#cd ..; rm -rf tmp/

## the lastest version
mkdir 20160101
mkdir tmp/; cd tmp/
wget https://unite.ut.ee/sh_files/uchime_reference_dataset_01.01.2016.zip
unzip uchime_reference_dataset_01.01.2016.zip
mv uchime_reference_dataset_01.01.2016/uchime_reference_dataset_01.01.2016.fasta ../20160101/unite.fasta
cd ..; rm -rf tmp/
