# TODO:
# Improve rRNA definition
# Bisulfite stuff...
# Small RNA
# Variation
# mappability



module load mugqic/R
module load mugqic/bwa
module load mugqic/samtools
module load mugqic/picard


##
######### mm10

## Fetch from iGenomes
GENOME="mm10"
cd $MUGQIC_INSTALL_HOME/genomes
wget "ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz"
tar -xvf Mus_musculus_UCSC_mm10.tar.gz
ROOT="$MUGQIC_INSTALL_HOME/genomes/Mus_musculus/UCSC/mm10"

# Bowtie1, Bowtie2, BWA indices
# (already done ...  iGenomes)

# fasta index samtools faidx
# (already done ...  iGenomes)

# Reference dictionary (not done for iGenomes mm10)
cd $ROOT/Sequence/WholeGenomeFasta
java -jar $PICARD_HOME/CreateSequenceDictionary.jar REFERENCE=genome.fa  OUTPUT=genome.dict


## Link WholeGenomeFasta/* in aligner indices  dirs
for dir in $ROOT/Sequence/{BWA,Bowtie,Bowtie2}Index
do
 cd $dir
 ln -s "../WholeGenomeFasta/genome.fa.fai"
 ln -s "../WholeGenomeFasta/genome.dict"
done

## Chrom size file, simple tab-delimied 1st col chr , 2nd col chr length. (already in iGenomes Annotation/Genes/ChromInfo.txt)
# OR cut -f1,2 $ROOT/Sequence/WholeGenomeFasta/genome.fa.fai > $ROOT/Sequence/WholeGenomeFasta/genome.chromsize.txt
# OR  sed 1d fasta/hg19.fasta.dict | awk '{print $2"\t"$3}' | sed -e 's/SN://g' | sed -e 's/LN://g' > hg19.chromsize.txt # example Maxime

## rRNA sequence + BWA index them (either from the Broad OR from iGenomes) (ln -s to). Put into $ROOT/Sequence/AbundantSequences/rRNA.fa
# cat $ROOT/Sequence/AbundantSequences/hum5SrDNA.fa $ROOT/Sequence/AbundantSequences/humRibosomal.fa > $ROOT/Sequence/AbundantSequences/rRNA.fa
# OR :
cd $ROOT/Sequence/AbundantSequences
ln -s musRibosomal.fa rRNA.fa
bwa index rRNA.fa

## Gene models aka gtf.   (Already done for iGenomes, Annotation/Genes/genes.gtf)
# ...

## Do a bunch of stuff in R :Gene size file for the gene model + gene id to symbol mapping + gene annotation file 
cd $ROOT/Annotation/Genes
R --vanilla <<'EOF'
	require(org.Mm.eg.db) # ORGANISM SPECIFIC!!!!
	org.prefix = 'org.Mm.eg'
	require(stringr)
	require(rtracklayer)
	options(stringsAsFactors=FALSE)
	x = import.gff('genes.gtf')
	x$geneid = gsub('((gene_id |;)|\")','', str_extract(as.character(x$group),'gene_id.*?;') ) # MIGHT DIFFER FOR NON-UCSC GENOMES

	# Write the gene_id to symbol file
	g2s  = data.frame('ensembl'=unique(x$geneid),'geneSymbol'=unique(x$geneid))
	write.table(g2s,file='genes_geneid2Symbol.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

	# Prepare the gene lengths file
	gl = gl[gl$type=='exon',]
	gl =  as(gl,'GRanges') 		
	gl = split(gl,f=factor(x$geneid))	 # split for GRanges signature turns into GRangesList!
	gl = reduce(gl) # reduce eliminates overlaps
	gl = sum ( width(gl) )
	gl = as.data.frame(gl)
	write.table(gl,file='genes_lengths.txt',col.names=FALSE,row.names=TRUE,sep='\t',quote=FALSE)

	# Gene Annoation file (Using Bioc USE a different STRATEGY if not org package available!!!)
	# NOTE: this strategy is not perfect for UCSC iGenomes (some genes unmapped). But then nothing important depends on these annotation
	ann = data.frame("featureID"=unique(x$geneid),row.names=unique(x$geneid))
	ann[['Entrez Gene ID']] =  mget(ann$featureID,get( paste0(org.prefix,'SYMBOL2EG') ), ifnotfound=NA)
	# which(sapply(ann[,2],length)>=2) # ok checked manually
	ann[['Entrez Gene ID']] = sapply(ann[['Entrez Gene ID']] , function(z)z[1])
	ann[['Entrez Gene ID']][is.na(ann[['Entrez Gene ID']])] = "NA"
	for(enviro in c("SYMBOL","CHR","GENENAME","UNIGENE","ENZYME")) # MIGT be specific see e.g. ls("package:org.Mm.eg.db")
	{
		print(enviro)
		ann[[enviro]] =  mget( ann[['Entrez Gene ID']],  get( paste0(org.prefix,enviro) )  , ifnotfound=NA)
		print(max(sapply(ann[[enviro]],length)))
		ann[[enviro]] = sapply(ann[[enviro]] ,paste, collapse=', ') # this is kind of slow
	}
	write.table(ann,file='genes.txt',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

EOF
cd -
# OR Alternative is MAthieu's:
# module load mugqic/tools
# gtf2geneSize.awk <...>










