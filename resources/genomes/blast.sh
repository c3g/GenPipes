#BLAST db ids and descriptions have been created in tab-delimited format using blastdbcmd:

#nr
#--
blastdbcmd -entry all -db nr -outfmt '%i %t' > nr.id.title

#uniprot_sprot_2013_11
#---------------------
#ids are not readily available in UniProtKBr/Swiss-Prot db but present in the title therefore:

blastdbcmd -entry all -db uniprot_sprot_2013_11 -outfmt '%i %t' > uniprot_sprot_2013_11.id.title.tmp
perl -pe 's/^No ID available (\S+) /\1\t/' uniprot_sprot_2013_11.id.title.tmp > uniprot_sprot_2013_11.id.title
rm uniprot_sprot_2013_11.id.title.tmp

# Trinotate UniProt resources (adjust file names for newer Trinotate version)
TRINOTATE_VERSION=2.0.1
cd $MUGQIC_INSTALL_HOME/genomes/blast_db
SPROT=uniprot_sprot.trinotate_v${TRINOTATE_VERSION%.*}.pep
UNIREF90=uniprot_uniref90.trinotate_v${TRINOTATE_VERSION%.*}.pep
for RESOURCE in \
  $SPROT \
  $UNIREF90 \
; do
  wget "ftp://ftp.broadinstitute.org/pub/Trinity/Trinotate_v${TRINOTATE_VERSION%.*}_RESOURCES/$RESOURCE.gz"
  gunzip $RESOURCE.gz
done

module load mugqic/blast/2.2.29+
makeblastdb -in $SPROT -dbtype prot
makeblastdb -in $UNIREF90 -dbtype prot
