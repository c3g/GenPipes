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
