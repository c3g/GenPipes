# Install genome annotations (to be sure we'll have the same sequence names)
# SNPeff
module load mugqic/java/oracle-jdk1.7.0_15 mugqic/snpEff/3.3 && java -Xmx4g -jar $SNPEFF_HOME/snpEff.jar download -c $SNPEFF_HOME/snpEff.config -v athalianaTair10

# Reference genomes were generated using a bash script

wget ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas -O tmp.fasta

# change the sequence names for chloroplast and mitochondria to fit the SNPeff annotation

cat tmp.fasta | sed -e 's/[>]mitochondria/\>M/g' | sed -e 's/[>]chloroplast/\>C/g' > tmp
mv tmp tmp.fasta

module add mugqic/tools/0.1  mugqic/bowtie/2.1.0 mugqic/samtools/0.1.19 && install_genomes.sh genomes_list.tsv

# Gene annotations
mkdir -p annotations
wget ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff -O annotations/TAIR10_GFF3_genes.gff 
cat annotations/TAIR10_GFF3_genes.gff  | sed -e 's/^Chr//g' > tmp
mv tmp annotations/TAIR10_GFF3_genes.gff 

# Convert TAIR10_GFF3_genes.gff to .gtf
module load mugqic/cufflinks && gffread -E TAIR10_GFF3_genes.gff -T -o- > TAIR10_GFF3_genes.gtf

# genesize

~/src/mugqic_tools/tools/gtf2geneSize.awk TAIR10_GFF3_genes.gtf geneSize.tsv

awk ' BEGIN {FS="\""} {print $4 "\t" $6} ' TAIR10_GFF3_genes.gtf | sort -u > TAIR10_GFF3_genesID2genename.tsv

# Split by chromosome

mkdir fasta/byChro
load mugqic/exonerate/2.2.0 && fastaexplode Arabidopsis_Thaliana_TAIR10.fasta -d fasta/byChro/
