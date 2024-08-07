.. _docs_steps_dnaseq:

.. spelling::

    Haplotype
    Mpileup
    GC
    FastQ
    config
    param
    bwa
    picard
    sam   
    BVATools
    vcf
    gvcfs
    recalibrator
    recalibrate
    Omni
    Biomart
    dbNSFP
    dbSNP
    Polyphen
    PhyloP
    transversion
    gz
    nb
    varfilter
    snpAndIndelBCF
    covarying
    bcf
    Mappability
    translocations
    bgzipped
    vcfs
    indv
    missingness
    idepth
    readgroups
    imiss

**Step Details**

Following are the various steps that are part of GenPipes DNA Sequencing genomic analysis pipeline:

.. _Picard SAM to FastQ:

**Picard SAM to FastQ**

Convert SAM/BAM files from the input readset file into FASTQ format if FASTQ files are not already specified in the readset file. Do nothing otherwise.

.. _Sym Link FastQ:

**Sym Link FastQ**

Automated symbolic linking of FASTQ files.

.. _Step Trimmomatic:

**Trimmomatic**

If available, this step takes FastQ file from the readset as input. Otherwise, FastQ output files from the previous `Picard SAM to FastQ`_ step, where BAM files are converted to FastQ format, are utilized in this step.

Raw reads quality trimming and removing of Illumina adapters is performed using `Trimmomatic`_. If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'), it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets, only Adapter1 is used and left unchanged.

.. _Merge Trimmomatic Stats:

**Merge Trimmomatic Stats**

The trim statistics per readset file are merged in this step.

.. _GATK Indel Re-aligner:

**GATK Indel Re-aligner**

Insertion and deletion realignment is performed on regions where multiple base mismatches are preferred over indels by the aligner since it can appear to be less costly by the algorithm. Such regions will introduce false positive variant calls which may be filtered out by realigning those regions properly. Realignment is done using `GATK Software`_. The reference genome is divided by a number regions given by the nb_jobs parameter.

.. _Fix Mate by Coordinate:

**Fix Mate by Coordinate**

This step fixes the read mates. Once local regions are realigned, the read mate coordinates of the aligned reads need to be recalculated since the reads are realigned at positions that differ from their original alignment. Fixing the read mate positions is done using`BVATools`_.

.. _Picard Mark Duplicates:

**Picard Mark Duplicates**

This step marks duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score) will be marked as a duplicate in the BAM file. Marking duplicates is done using `Picard`_.

.. _Recalibration:

**Recalibration**

This step is used to recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome. Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so, provides not only more accurate quality scores but also more widely dispersed ones.

.. note::

   The unmapped reads are extracted from the merged/sorted bam files. After the base recalibration step, the unmapped reads are re-injected into the FASTQ. This ensures that the unmapped reads are available for further analysis, if needed. Without this process, the base recalibration may eventually remove the unmapped reads.

.. _Sym Link Final BAM:

**Sym Link Final BAM**

Creates symbolic links of final BAM for delivery of data to the clients.

.. _Metrics DNA Picard:

**Metrics DNA Picard**

This step generates metrics with picard, including `CollectMultipleMetrics <https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard->`_, `CollectOxoGMetrics <https://gatk.broadinstitute.org/hc/en-us/articles/360037428231-CollectOxoGMetrics-Picard->`_, and `CollectGcBiasMetrics <https://gatk.broadinstitute.org/hc/en-us/articles/360036481572-CollectGcBiasMetrics-Picard->`_.

.. _Metrics DNA Sample Quality Map:

**Metrics DNA Sample Quality Map**

Generates metrics with `Qualimap BAM QC <http://qualimap.conesalab.org/doc_html/analysis.html#bamqc>`_. BAM QC reports information for the evaluation of the quality of the provided alignment data (for e.g., a BAM file). In short, it provides the basic statistics of the alignment  such as the number of reads, coverage, GC-content, etc., in a summarized form along with a number of useful graphs are produced.

.. _Metrics DNA SAMBAM Flag Stats:

**Metrics DNA SAMBAM Flag Stats**

Outputs `flag statistics <https://lomereiter.github.io/sambamba/docs/sambamba-flagstat.html>`_ from the BAM file.

.. _Metrics DNA FastQC:

**Metrics DNA FastQC**

This step performs quality control with `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`.

.. _Picard Calculate HS Metrics:

**Picard Calculate HS Metrics**

Compute on target percent of hybridization based capture happens in this step.

.. _GATK Callable Loci:

**GATK Callable Loci**

Computes the callable region or the genome as a bed track.

.. _Extract Common SNP Frequencies:

**Extract Common SNP Frequencies**

Extracts allele frequencies of possible variants across the genome.

.. _BAF Plot:

**BAF Plot**

Plots DepthRatio and B allele frequency of previously extracted alleles.

.. _GATK Haplotype Caller:

**GATK Haplotype Caller**

GATK Haplotype Caller is used for SNPs and small indels.

.. _Merge and call individual GVCF:

**Merge and call individual GVCF**

Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.

.. _Combine GVCF:

**Combine GVCF**

Combine the per sample gvcfs of haplotype caller into one main file for all sample.

.. _Merge and call combined GVCF:

**Merge and call combined GVCF**

Merges the combined gvcfs and also generates a general vcf containing genotypes.

.. _Variant Recalibrator:

**Variant Recalibrator**

This step involves GATK variant recalibrator. The purpose of the variant recalibrator is to assign a well-calibrated probability to each variant call in a call set. You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call. The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship between SNP call annotations (QD, MQ, HaplotypeScore, and ReadPosRankSum, for example) and the probability that a SNP is a true genetic variant versus a sequencing or data processing artifact. This model is determined adaptively based on "true sites" provided as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP chip array. This adaptive error model can then be applied to both known and novel variation discovered in the call set of interest to evaluate the probability that each call is real. The score that gets added to the INFO field of each variant is called the VQSLOD. It is the log odds ratio of being a true variant versus being false under the trained Gaussian mixture model. Using the tranche file generated by the previous step the ApplyRecalibration walker looks at each variant's VQSLOD value and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level have their filter field annotated with its tranche level. This will result in a call set that simultaneously is filtered to the desired level but also has the information necessary to pull out more variants for a higher sensitivity but a slightly lower quality level.

.. _Haplotype caller decompose and normalize:

**Haplotype caller decompose and normalize**

Performs decompose/normalization for variant comparison at the haplotype level. Replay the variants from the VCF into the reference and determine whether variants match by whether the resulting Haplotype match.

.. _Haplotype caller flag mappability:

**Haplotype caller flag mappability**

Mappability annotation applied to haplotype caller vcf. An in-house database identifies regions in which reads are confidently mapped to the reference genome.

.. _Haplotype caller SNP ID annotation:

**Haplotype caller SNP ID annotation**

dbSNP annotation applied to haplotype caller vcf. The .vcf files are annotated for dbSNP using the software SnpSift (from the `SNPEff Suite`_).

.. _Haplotype caller SNP Effect:

**Haplotype caller SNP Effect**

Variant effect annotation applied to haplotype caller vcf. The .vcf files are annotated for variant effects using the `SNPEff Suite`_ software. SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

.. _Haplotype caller dbNSFP annotation:

**Haplotype caller dbNSFP annotation**

Additional SVN annotations applied to haplotype caller vcf. Provides extra information about SVN by using numerous published databases. Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information) and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy) and other function annotations).

.. _Haplotype caller Gemini annotation:

**Haplotype caller Gemini annotation**

Uses `Haplotype caller <https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller>`_ to call germline SNPs and indels via local re-assembly of Haplotype for exploring genetic variations using the `Gemini annotations`_. 

.. _Haplotype caller metrics VCF stats:

**Haplotype caller metrics VCF stats**

Metrics SNV applied to haplotype caller vcf. Multiple metrics associated to annotations and effect prediction are generated at this step: change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect, counts by genomic region, SNV quality, coverage, InDel lengths, base changes, transition-transversion rates, summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

.. _Run MultiQC:

**Run MultiQC**

Aggregate results from bioinformatics analyses across many samples into a single report. `MultiQC`_ searches a given directory for analysis logs and compiles a HTML report. It's a general use tool, perfect for summarizing the output from numerous bioinformatics tools.

.. _Raw MPileup:

**Raw MPileup**

Full pileup (optional). A raw Mpileup file is created using samtools mpileup and compressed in gz format. One packaged Mpileup file is created per sample/chromosome.

.. _Compress Raw MPileup:

**Compress Raw MPileup**

Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

.. _SNP and indel BCF:

**SNP and indel BCF**

Mpileup and Variant calling. Variants (SNPs and INDELs) are called using `SAMTools software package`_ Mpileup. bcftools view is used to produce binary bcf files.

.. _Merge Filter BCF:

**Merge Filter BCF**

bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step. The output of bcftools is fed to varfilter, which does an additional filtering of the variants and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls for all samples in the experiment.

.. _MPileup decompose and normalize:

**MPileup decompose and normalize**

Uses multi-way pileup that produces genotype likelihoods for decomposing and normalizing complex variants in VCF files.

.. _MPileup Flag Mappability:

**MPileup Flag Mappability**

Mappability annotation applied to mpileup vcf. An in-house database identifies regions in which reads are confidently mapped to the reference genome.

.. _Mpileup SNP ID annotation:

**Mpileup SNP ID annotation**

dbSNP annotation applied to `MPileup`_ vcf. The .vcf files are annotated for dbSNP using the software SnpSift (from the `SNPEff Suite`_).

.. _MPileup SNP Effect:

**MPileup SNP Effect**

Variant effect annotation applied to mpileup vcf. The .vcf files are annotated for variant effects using the SnpEff software. SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

.. _MPileup dbNSFP annotation:

**MPileup dbNSFP annotation**

Additional SVN annotations applied to mpileup vcf. Provides extra information about SVN by using numerous published databases. Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information) and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy) and other function annotations).

.. _MPileup Gemini annotation:

**MPileup Gemini annotation**

Explores genetic variations in the VCF file generated after Mpileup using `Gemini annotations`_.

.. _MPileup Metrics VCF stats:

**MPileup Metrics VCF stats**

Metrics SNV applied to mpileup caller vcf. Multiple metrics associated to annotations and effect prediction are generated at this step: change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect, counts by genomic region, SNV quality, coverage, InDel lengths, base changes, transition-transversion rates, summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

.. _SAMBAM Mark Duplicates:

**SAMBAM Mark Duplicates**

In this step duplicates are marked. Aligned reads per sample are duplicates if they have the same 5' alignment positions (for both mates in case of paired-end reads). All but the best pair (based on alignment score) will be marked as duplicates in the BAM file. Marking duplicates is done using `Picard`_.

.. _Delly2 Call Filter:

**Delly2 Call Filter**

This step uses normal and tumor final BAMs as input and generates a binary variant call format (BCF) file as output. It utilizes `Delly2`_, an integrated structural variant (SV) prediction method that can discover genotype and visualize deletions, tandem duplications, inversions and translocations at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends, split-reads and read-depth to sensitively and accurately delineate genomic rearrangements throughout the genome. Structural variants can be annotated using `Delly-sansa`_ and visualized using `Delly-maze`_ or `Delly-suave`_.

.. _Delly2 SV Annotation:

**Delly2 SV Annotation**

SV Annotation step utilizes the BCF file generated in previous `Delly2 Call Filter`_ step and performs genome annotation at various levels. At the nucleotid level it tries to identify the physical location of the SV dna sequences. Next, at the protein level the annotation process tries to determine the possible functions of the SV genes. Lastly, at the process-level annotation, it tries to identify the pathways and processes in which different SV genes interact, assembling an efficient functional annotation. For more details on annotation see `Genome Annotations`_. 

.. _Manta SV Calls:

**Manta SV Calls**

`Manta`_ calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.  Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a single efficient workflow.

Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format

.. _Manta SV Annotation:

** Manta SV Annotation**

This step uses the VCF file generated in previous step and performs SV annotations to compares types and breakpoints for candidate SVs from different callsets and enables fast comparison of SVs to genomic features such as genes and repetitive regions, as well as to previously established SV datasets.

.. _Lumpy Paired SV:

**Lumpy Paired SV**

This step uses `Lumpy`_ for structural variant discovery in the given input file. The output is available in BAM format.

Comprehensive discovery of structural variation (SV) from whole genome sequencing data requires multiple detection signals including read-pair, split-read, read-depth and prior knowledge. Owing to technical challenges, extant SV discovery algorithms either use one signal in isolation, or at best use two sequentially. Lumpy is a novel SV discovery framework that naturally integrates multiple SV signals jointly across multiple samples. It yields improved sensitivity, especially when SV signal is reduced owing to either low coverage data or low intra-sample variant allele frequency.

.. _Lumpy SV Annotation:

**Lumpy SV Annotation**

This step performs LumPy SV Annotation for mapping and characterization of SVs.

.. _Wham SV Call:

**Wham SV Call**

`Wham`_ (Whole-genome Alignment Metrics) provides a single, integrated framework for both structural variant calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts to employ SVs in association testing.  This step returns a VCF file.

.. _Wham SV Annotation:

**Wham SV Annotation**

This step uses the VCF file generated in the previous step `Wham SV Call`_ and performs SV annotations. 

.. _CNVkit Batch:

**CNVkit Batch**

A copy number variation (`CNV`_) is when the number of copies of a particular gene varies from one individual to the next. Copy-number variation (CNV) is a large category of structural variation, which includes insertions, deletions and duplications. For copy number variation analysis, GenPipes DNA Sequencing pipeline (-t sv option) uses CNVkit that allows for CNV calling on single samples (e.g., tumor samples). 

CNVkit provides an advantageous way to run the entire pipeline using the batch option to run various stages in copy number calling pipeline such as:

* Create target/anti-target bed files
* Gather read depths for those regions
* Compile a copy number reference
* Correct biases in tumor samples while calculating copy ratios
* Mark copy number segments

.. _CNVkit SV Annotation:

**CNVkit SV Annotation**

This step performs CNVkit SV annotation.

.. _Run BreakSeq2:

**Run BreakSeq2**

In this step, `BreakSeq2`_ is used to combine DNA double-strand breaks (DSBs) labeling with next generation sequencing (NGS) to map chromosome breaks with improved sensitivity and resolution. It is an ultra fast and accurate nucleotide-resolution analysis of structural variants.

.. _Ensemble MetaSV:

**Ensemble MetaSV**

`MetaSV`_ uses highly effective ensemble approach for calling SVs. It is an integrated SV caller which leverages multiple orthogonal SV signals for high accuracy and resolution. MetaSV proceeds by merging SVs from multiple tools for all types of SVs. It also analyzes soft-clipped reads from alignment to detect insertions accurately since existing tools underestimate insertion SVs. Local assembly in combination with dynamic programming is used to improve breakpoint resolution. Paired-end and coverage information is used to predict SV genotypes.

.. _MetaSV Annotation:

**MetaSV Annotation**

This step uses output from previous step and performs SV annotations.

.. _Metrics:

**Metrics**

This step computes metrics and generates coverage tracks per sample. Multiple metrics are computed at this stage:

* Number of raw reads,
* Number of filtered reads,
* Number of aligned reads,
* Number of duplicate reads,
* Median, mean and standard deviation of insert sizes of reads after alignment,
* Percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
* Whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads).

A TDF (.tdf) coverage track is also generated at this step for easy visualization of coverage in the IGV browser.

.. _Metrics NGSCheckmate:

**Metrics NGSCheckmate**

`NGSCheckMate <https://github.com/parklab/NGSCheckMate>`_ is a software package for identifying next generation sequencing (NGS) data files from the same individual. It analyzes various types of NGS data files including (but not limited to) whole genome sequencing (WGS), whole exome sequencing (WES), RNA-seq, ChIP-seq, and targeted sequencing of various depths. Data types can be mixed (e.g. WES and RNA-seq, or RNA-seq and ChIP-seq). It takes BAM (reads aligned to the genome), VCF (variants) or FASTQ (unaligned reads) files as input. NGSCheckMate uses depth-dependent correlation models of allele fractions of known single-nucleotide polymorphisms (SNPs) to identify samples from the same individual. 

This step takes as input the file containing all vcfs in project output.

.. _Metrics Verify BAM ID:

**Metrics Verify BAM ID**

In this step, `VerifyBAMID <https://genome.sph.umich.edu/wiki/VerifyBamID>`_ software is used to verify whether the reads in particular file match previously known genotypes for an individual (or group of individuals), and checks whether the reads are contaminated as a mixture of two samples. 

.. _Metrics VCFTools Missing Individual:

**Metrics VCFTools Missing Individual**

This step uses `VCFtools <https://vcftools.github.io/index.html>`_ and --missing_indv option to generate a file reporting the `missingness` factor in the analysis on a per-individual basis. It takes bgzipped .vcf file as input and creates .imiss flat file indicating missingness.

.. _Metrics VCFTools Depth Individual:

**Metrics VCFTools Depth Individual**

This step uses `VCFtools <https://vcftools.github.io/index.html>`_ and --depth option to generate a file containing the mean depth per individual. It takes as input bgzipped .vcf file and creates a .idepth flat file.

.. _Metrics GATK Sample Fingerprint:

**Metrics GATK Sample Fingerprint**

`CrosscheckFingerprints (Picard) <https://gatk.broadinstitute.org/hc/en-us/articles/360040098592-CrosscheckFingerprints-Picard->`_ functionality in GATK toolkit is used to cross-check readgroups, libraries, samples, or files to determine if all data in the set of input files appears to come from the same individual. In this step, sample SAM/BAM or VCF file is taken as input and a fingerprint file is generated using CrosscheckFingerprints (Picard) in GATK. It checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) against a set of known genotypes in the supplied genotype file (in VCF format).

.. _Metrics GATK Cluster Fingerprint:

**Metrics GATK Cluster Fingerprint**

In this step, ClusterCrosscheckMetrics function from GATK is used as a follow-up step to running CrosscheckFingerprints created in the `Metrics GATK Sample Fingerprint`_ step earlier. There are cases where one would like to identify a few groups out of a collection of many possible groups (say to link a bam to it's correct sample in a multi-sample vcf. In this step, sample SAM/BAM or VCF file is taken as input and a fingerprint file is generated.

.. include:: steps_dnaseq_common.inc
.. include:: steps_cram_op.inc

.. Following are the html link references:

.. _Trimmomatic: http://www.usadellab.org/cms/index.php?page=trimmomatic
.. _GATK Software: https://gatk.broadinstitute.org/hc/en-us
.. _BVATools: https://bitbucket.org/mugqic/bvatools/src/master/
.. _SNPEff Suite: http://snpeff.sourceforge.net
.. _Delly2: https://github.com/dellytools/delly
.. _Delly-sansa: https://github.com/dellytools/sansa
.. _Delly-maze: https://github.com/dellytools/maze
.. _Delly-suave: https://github.com/dellytools/suave
.. _Genome Annotations: https://www.sciencedirect.com/topics/biochemistry-genetics-and-molecular-biology/genome-annotation
.. _Manta: https://github.com/Illumina/manta
.. _Lumpy: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84
.. _Wham: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004572
.. _CNV: https://www.genome.gov/genetics-glossary/Copy-Number-Variation
.. _BreakSeq2: http://bioinform.github.io/breakseq2/
.. _MetaSV: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4528635/
.. _Gemini annotations: https://gemini.readthedocs.io
