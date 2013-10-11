
# Instructions used to obtain genome reference

wget ftp://ftp.sanger.ac.uk/pub/1000genomes/tk2/main_project_reference/human_g1k_v37.fasta.gz 
gunzip human_g1k_v37.fasta.gz
mkdir fasta
mv human_g1k_v37.fasta fasta/hg1k_v37.fasta
ln -s hg1k_v37.fasta fasta/hg1k_v37.fasta.fa

# Run bowtie-build

/software/areas/genomics/aligners/bowtie-0.12.8/bowtie-build /lb/project/mugqic/epigenome/genome_files/hg1k_v37/fasta/hg1k_v37.fasta /lb/project/mugqic/epigenome/genome_files/hg1k_v37/fasta/hg1k_v37.fasta

# Download rnaseqc rRNA fasta files

mkdir rnaseqc
wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/examples/rRNA/rRNA.tar.gz
tar -zxvf rRNA.tar.gz -C rnaseqc/
rm rnaseqc/human_all_rRNA.fasta.*
/software/areas/genomics/aligners/bwa-0.6.1/bwa index rnaseqc/human_all_rRNA.fasta

# Index reference with samtools

/software/areas/genomics/tools/samtools-0.1.18/samtools faidx fasta/hg1k_v37.fasta

# Create reference dictionary

java -jar /software/areas/genomics/tools/picard-tools-1.70/CreateSequenceDictionary.jar REFERENCE=fasta/hg1k_v37.fasta OUTPUT=fasta/hg1k_v37.fasta.dict

# Generate tab delimited chromosome size file

sed 1d fasta/hg1k_v37.fasta.dict | awk '{print $2"\t"$3}' | sed -e 's/SN://g' | sed -e 's/LN://g' > hg1k_v37.chromsize.txt
grep -v "GL" /sb/programs/analyste/genomes/hg1k_v37/hg1k_v37.chromsize.txt > /sb/programs/analyste/genomes/hg1k_v37/hg1k_v37.AutosomeSize.txt

# Create BWA indexes

/software/areas/genomics/aligners/bwa-0.6.1/bwa index -a bwtsw fasta/bwa/hg1k_v37.fasta
ln -s fasta/bwa/hg1k_v37.fasta fasta/bwa/

# Download Ensembl v68 annotations

wget ftp://ftp.ensembl.org/pub/release-68/gtf/homo_sapiens/Homo_sapiens.GRCh37.68.gtf.gz
gunzip Homo_sapiens.GRCh37.68.gtf.gz
mv Homo_sapiens.GRCh37.68.gtf annotations/
ln -s /lb/project/mugqic/epigenome/genome_files/hg1k_v37/annotations/Homo_sapiens.GRCh37.68.gtf annotations/transcripts.gtf

# Download lambda reference sequence

mkdir lambda
wget http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi --post-dat 'db=nuccore&id=NC_001416.1&rettype=fasta&retmode=text' -O lambda/lambda.fa

# Append lambda sequence to hg1k fasta, create converted C->T and G->A genome sequences

mkdir fasta/Bisulfite_Genome
cp fasta/hg1k_v37.fasta fasta/Bisulfite_Genome/
echo ">lambda" >> fasta/Bisulfite_Genome/hg1k_v37.fasta
grep -vP '^>' lambda/lambda.fa >> fasta/Bisulfite_Genome/hg1k_v37.fasta
echo -e "\n" >> fasta/Bisulfite_Genome/hg1k_v37.fasta

cat fasta/Bisulfite_Genome/hg1k_v37.fasta | perl -pe 'if(!(/^>/)){ tr/cC/tT/ }' > fasta/Bisulfite_Genome/hg1k_v37.ct.fasta
cat fasta/Bisulfite_Genome/hg1k_v37.fasta | perl -pe 'if(!(/^>/)){ tr/gG/aA/ }' > fasta/Bisulfite_Genome/hg1k_v37.ga.fasta
cp fasta/Bisulfite_Genome/hg1k_v37.fasta fasta/Bisulfite_Genome/hg1k_v37.fasta.1
gzip fasta/Bisulfite_Genome/hg1k_v37.fasta fasta/Bisulfite_Genome/hg1k_v37.ct.fasta fasta/Bisulfite_Genome/hg1k_v37.ga.fasta
mv fasta/Bisulfite_Genome/hg1k_v37.fasta.1 fasta/Bisulfite_Genome/hg1k_v37.fasta

# Index with bwa

/software/areas/genomics/aligners/bwa-0.6.1/bwa index -a bwtsw fasta/Bisulfite_Genome/hg1k_v37.ct.fasta.gz
/software/areas/genomics/aligners/bwa-0.6.1/bwa index -a bwtsw fasta/Bisulfite_Genome/hg1k_v37.ga.fasta.gz

# Create exon bed files per transcript category

mkdir bed_files
rm bed_files/*
awk '$3=="exon"' annotations/transcripts.gtf > forbed.gtf
perl splitAnnotation.pl transcript_types.csv forbed.gtf 

# Sort and merge bed files

for files in bed_files/*.bed
do
a=`echo $files | sed -e 's/.bed//g'`
sort -k1n,1n -k2n,2n $files > $a.sorted.bed
done

for files in bed_files/*.sorted.bed
do
a=`echo $files | sed -e 's/.sorted.bed//g'`
/software/areas/genomics/tools/BEDTools-Version-2.16.2/bin/mergeBed -nms -i $files | sort -k1n,1n -k2,2n > $a.merged.bed
done

# generate genesize file (need for saturation)

gtf2geneSize.awk Homo_sapiens.GRCh37.66.gtf Homo_sapiens.GRCh37.66_geneSize.tsv

# get dbSnp variant file
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz.tbi -O dbSnp-138.vcf.gz.tbi
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz -O dbSnp-138.vcf.gz

# get the dbnsfp annotations
wget http://dbnsfp.houstonbioinformatics.org/dbNSFPzip/dbNSFP2.0.zip
unzip dbNSFP2.0.zip
(head -n 1 dbNSFP2.0_variant.chr1 ; cat dbNSFP2.0_variant.chr* | grep -v "^#" ) > dbNSFP2.0.txt

sed 's/.*SN:\([^	]\+\).*LN:\([0-9]\+\).*/\1	\2/g' hg1k_v37.fasta.dict > hg1k_v37.fasta.length
