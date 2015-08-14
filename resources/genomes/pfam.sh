##
ROOT="$MUGQIC_INSTALL_HOME/genomes/pfam_db"; mkdir -p $ROOT ; cd $ROOT
module load mugqic/hmmer/3.1b1 # Assuming any version will do...

## pfam-A
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm


## pfam-AB, as suggested by Trinotate authors for non-model species
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.gz
gunzip -c Pfam-A.hmm.gz Pfam-B.hmm.gz > Pfam-AB.hmm
hmmconvert -b Pfam-AB.hmm > Pfam-AB.hmm.bin
rm Pfam-A.hmm.gz Pfam-B.hmm.gz 
hmmpress Pfam-AB.hmm


# wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz ftp://ftp.sanger.ac.uk/pub/databases/Pfam/current_release/Pfam-B.hmm.gz
# gunzip -c Pfam-A.hmm.gz Pfam-B.hmm.gz > Pfam-AB.hmm
# hmmconvert -b Pfam-AB.hmm > pfam/Pfam-AB.hmm.bin
# rm Pfam-A.hmm.gz Pfam-B.hmm.gz Pfam-AB.hmm pfam/Pfam-AB.hmm.bin.* ## ???
# make prep_pfam ## ???

# The following HMMER algorithms/programs are supported by this server:
# phmmer - used to search one or more query protein sequences against a protein sequence database.
# hmmscan - search protein sequences against collections of profiles, e.g. Pfam. In HMMER2 this was called hmmpfam.
# hmmsearch - used to search one or more profiles against a protein sequence database.
# jackhmmer - iteratively search a query protein sequence, multiple sequence alignment or profile HMM against the target protein sequence database.
# Longer term, we will support nhmmer that will enable searches with nucleotide sequences. This software has been released as part of the HMMER software package (version 3.1).
# Other Programs
# The following is a brief description of the other programs in the HMMER suite. These are only available from downloaded distributions. However, they are used indirectly when performing the searches on the server.
# hmmalign performs a multiple sequence alignment of all the sequences (usually identified by running an hmmsearch) in the input, by aligning them individually to the profile HMM.
# hmmbuild builds a profile HMM for each multiple sequence alignment in the input multiple sequence alignment file, and saves it to a new file.
# hmmconvert utility converts an input profile file to different HMMER formats.
# hmmfetch retrieves one or more profile HMMs from a profile database (e.g. Pfam).
# hmmpress takes a profile database in standard HMMER3 format and constructs binary compressed data files for hmmscan.
# hmmstat utility prints out a tabular file of summary statistics for each profile.
# 
# 
# 
# 