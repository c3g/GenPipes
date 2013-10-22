

######### VirusFinder ##################
VERSION="1.0"
NAME="virus_finder" # same could apply to all ucsc tools
INSTALL_PATH=$MUGQIC_INSTALL_HOME/software/$NAME
mkdir -p $INSTALL_PATH
cd $INSTALL_PATH
# Download and extract 
wget http://bioinfo.mc.vanderbilt.edu/VirusFinder/VirusFinder-$VERSION.tgz
tar -xvzf VirusFinder-$VERSION.tgz
mv VirusFinder VirusFinder-$VERSION

module add mugqic/blat/35x1
module add mugqic/blast/2.2.27+
module add mugqic/samtools/0.1.18
module add mugqic/bowtie/2.1.0
module add mugqic/bwa/0.7.4
module add mugqic/trinity/2013-02-25
module add mugqic/svdetect/0.8b
#module add mugqic/blat/35x1 mugqic/blast/2.2.27+ mugqic/samtools/0.1.18 mugqic/bowtie/2.1.0 mugqic/bwa/0.7.4 mugqic/trinity/2013-02-25 mugqic/svdetect/0.8b

echo "###########################################################################" > VirusFinder-1.0/template-config.txt
echo "#" >> VirusFinder-1.0/template-config.txt
echo "# Configuration file for VirusFinder" >> VirusFinder-1.0/template-config.txt
echo "#" >> VirusFinder-1.0/template-config.txt
echo "###########################################################################" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo "## Input NGS data can be either of the following two options:" >> VirusFinder-1.0/template-config.txt
echo "##  (a) an alignment file (in BAM format)," >> VirusFinder-1.0/template-config.txt
echo "##  (b) two paired-end FASTQ files. " >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "alignment_file = /scratch/kingw/VirusFinder/simulation/simulation.bam" >> VirusFinder-1.0/template-config.txt
echo "#fastq1        = /scratch/kingw/VirusFinder/simulation/seq_1.fastq.gz" >> VirusFinder-1.0/template-config.txt
echo "#fastq2        = /scratch/kingw/VirusFinder/simulation/seq_2.fastq.gz" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "mailto         = \$MAIL" >> VirusFinder-1.0/template-config.txt
echo "thread_no      = 8" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo "## The full paths to the following third-party tools are required by VirusFinder:" >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "samtools_bin    = samtools" >> VirusFinder-1.0/template-config.txt
echo "blat_bin        = blat" >> VirusFinder-1.0/template-config.txt
echo "blastn_bin      = blastn" >> VirusFinder-1.0/template-config.txt
echo "bowtie_bin      = bowtie2" >> VirusFinder-1.0/template-config.txt
echo "bwa_bin         = bwa" >> VirusFinder-1.0/template-config.txt
echo "trinity_script  = $MUGQIC_INSTALL_HOME/software/trinity/trinityrnaseq_r2013-02-25/Trinity.pl" >> VirusFinder-1.0/template-config.txt
echo "SVDetect_dir    = $MUGQIC_INSTALL_HOME/software/sv_detect/SVDetect_r0.8b/" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo "## Reference files indexed for Bowtie2, BWA, and BLAST" >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "virus_database     = $MUGQIC_INSTALL_HOME/genomes/blast_db/virus.fa" >> VirusFinder-1.0/template-config.txt
echo "bowtie_index_human = $MUGQIC_INSTALL_HOME/genomes/Homo_sapiens/hg19/fasta/bowtie2/" >> VirusFinder-1.0/template-config.txt
echo "blastn_index_human = $MUGQIC_INSTALL_HOME/genomes/Homo_sapiens/hg19/fasta/" >> VirusFinder-1.0/template-config.txt
echo "blastn_index_virus = $MUGQIC_INSTALL_HOME/genomes/blast_db/virus" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo "## Modifiable parameters (users are suggested not to modify them)" >> VirusFinder-1.0/template-config.txt
echo "################################" >> VirusFinder-1.0/template-config.txt
echo ""  >> VirusFinder-1.0/template-config.txt
echo "min_contig_length  = 300" >> VirusFinder-1.0/template-config.txt
echo "blastn_evalue_thrd = 0.05" >> VirusFinder-1.0/template-config.txt
echo "similarity_thrd    = 0.8"  >> VirusFinder-1.0/template-config.txt
echo "chop_read_length   = 25" >> VirusFinder-1.0/template-config.txt
echo "minIdentity        = 80" >> VirusFinder-1.0/template-config.txt






# Module file
echo "#%Module1.0
proc ModulesHelp { } {
       puts stderr \"\tMUGQIC - VirusFinder \"
}
module-whatis \"MUGQIC - VirusFinder \"
                      
set             root               \$::env(MUGQIC_INSTALL_HOME)/software/virus_finder/VirusFinder-$VERSION/bin/
prepend-path    PATH               \$root
" > $VERSION

# version file
echo "#%Module1.0
set ModulesVersion \"$VERSION\"
" > .version

mkdir -p $MUGQIC_INSTALL_HOME/modulefiles/mugqic/svdetect
mv .version $VERSION $MUGQIC_INSTALL_HOME/modulefiles/mugqic/svdetect/

