##
ROOT="$MUGQIC_INSTALL_HOME/genomes/silva_db/"; mkdir -p $ROOT ; cd $ROOT

## The last version is 111 
wget http://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_111_release.tgz
tar -zxvf Silva_111_release.tgz
mv Silva_111_post 111
rm -f Silva_111_release.tgz

## Rename files for configuration pattern.
cd 111/taxonomy/
mv 90_Silva_111_taxa_map.txt 90_otu_taxonomy.txt
mv 94_Silva_111_taxa_map.txt 94_otu_taxonomy.txt
mv 97_Silva_111_taxa_map.txt 97_otu_taxonomy.txt
mv 99_Silva_111_taxa_map.txt 99_otu_taxonomy.txt

cd $ROOT
cd 111/rep_set/
gunzip 90_Silva_111_rep_set.fasta.gz
mv 90_Silva_111_rep_set.fasta 90_otus.fasta
gunzip 94_Silva_111_rep_set.fasta.gz
mv 94_Silva_111_rep_set.fasta 94_otus.fasta
gunzip 97_Silva_111_rep_set.fasta.gz
mv 97_Silva_111_rep_set.fasta 97_otus.fasta
gunzip 99_Silva_111_rep_set.fasta.gz
mv 99_Silva_111_rep_set.fasta 99_otus.fasta


