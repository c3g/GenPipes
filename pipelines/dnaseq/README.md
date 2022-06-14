[TOC]



Usage
-----
```
#!text

usage: dnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {wrapper, singularity} <IMAGE PATH>]
                 [--genpipes_file GENPIPES_FILE]
                 [-t {mugqic,mpileup,light,sv}] [-r READSETS] [-v]

Version: 4.3.0

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -t {mugqic,mpileup,light,sv}, --type {mugqic,mpileup,light,sv}
                        DNAseq analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![dnaseq mugqic workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mugqic.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mugqic.png)
```
mugqic:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_extract_unmapped
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- picard_mark_duplicates
8- recalibration
9- gatk_haplotype_caller
10- merge_and_call_individual_gvcf
11- combine_gvcf
12- merge_and_call_combined_gvcf
13- variant_recalibrator
14- haplotype_caller_decompose_and_normalize
15- haplotype_caller_flag_mappability
16- haplotype_caller_snp_id_annotation
17- haplotype_caller_snp_effect
18- haplotype_caller_dbnsfp_annotation
19- haplotype_caller_gemini_annotations
20- metrics_dna_picard_metrics
21- metrics_dna_sample_qualimap
22- metrics_dna_fastqc
23- picard_calculate_hs_metrics
24- metrics
25- gatk_callable_loci
26- extract_common_snp_freq
27- baf_plot
28- run_multiqc
29- cram_output
30- sym_link_fastq
31- sym_link_final_bam
32- metrics_ngscheckmate
33- metrics_verify_bam_id
34- metrics_vcftools_missing_indiv
35- metrics_vcftools_depth_indiv
36- metrics_gatk_sample_fingerprint
37- metrics_gatk_cluster_fingerprint
----
```
![dnaseq mpileup workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mpileup.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mpileup.png)
```
mpileup:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_extract_unmapped
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- picard_mark_duplicates
8- recalibration
9- rawmpileup
10- rawmpileup_cat
11- snp_and_indel_bcf
12- merge_filter_bcf
13- mpileup_decompose_and_normalize
14- mpileup_flag_mappability
15- mpileup_snp_id_annotation
16- mpileup_snp_effect
17- mpileup_dbnsfp_annotation
18- mpileup_gemini_annotations
19- mpileup_metrics_vcf_stats
20- cram_output
21- metrics_dna_picard_metrics
22- metrics_dna_sample_qualimap
23- metrics_dna_fastqc
24- picard_calculate_hs_metrics
25- gatk_callable_loci
26- extract_common_snp_freq
27- baf_plot
28- run_multiqc
29- sym_link_fastq
30- sym_link_final_bam
----
```
![dnaseq light workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_light.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_light.png)
```
light:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_extract_unmapped
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- picard_mark_duplicates
8- recalibration
9- sym_link_final_bam
10- metrics_dna_picard_metrics
11- metrics_dna_sample_qualimap
12- metrics_dna_sambamba_flagstat
13- metrics_dna_fastqc
14- picard_calculate_hs_metrics
15- gatk_callable_loci
16- extract_common_snp_freq
17- baf_plot
18- gatk_haplotype_caller
19- merge_and_call_individual_gvcf
20- combine_gvcf
21- merge_and_call_combined_gvcf
22- variant_recalibrator
23- haplotype_caller_decompose_and_normalize
24- haplotype_caller_flag_mappability
25- haplotype_caller_snp_id_annotation
26- haplotype_caller_snp_effect
27- haplotype_caller_dbnsfp_annotation
28- haplotype_caller_gemini_annotations
29- run_multiqc
30- cram_output
----
```
![dnaseq sv workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_sv.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_sv.png)
```
sv:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_extract_unmapped
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- picard_mark_duplicates
8- recalibration
9- gatk_haplotype_caller
10- merge_and_call_individual_gvcf
11- metrics_dna_picard_metrics
12- delly_call_filter
13- delly_sv_annotation
14- manta_sv_calls
15- manta_sv_annotation
16- lumpy_paired_sv
17- lumpy_sv_annotation
18- wham_call_sv
19- wham_sv_annotation
20- cnvkit_batch
21- cnvkit_sv_annotation
22- run_breakseq2
23- ensemble_metasv
24- metasv_sv_annotation

```

picard_sam_to_fastq
-------------------
Converts SAM/BAM files from the input readset file into FASTQ format.
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

skewer_trimming
---------------
Trimming using [skewer](https://sourceforge.net/projects/skewer/)

bwa_mem_sambamba_sort_sam
-------------------------
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html)
This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

sambamba_merge_sam_extract_unmapped
-----------------------------------
BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

gatk_indel_realigner
--------------------
Insertion and deletion realignment is performed on regions where multiple base mismatches
are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
Such regions will introduce false positive variant calls which may be filtered out by realigning
those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
The reference genome is divided by a number regions given by the `nb_jobs` parameter.

sambamba_merge_realigned
------------------------
BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

picard_mark_duplicates
----------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).

recalibration
-------------
Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
the quality scores in the QUAL field in each read in the output BAM are more accurate in that
the reported quality score is closer to its actual probability of mismatching the reference genome.
Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
and sequence context, and by doing so, provides not only more accurate quality scores but also
more widely dispersed ones.

gatk_haplotype_caller
---------------------
GATK haplotype caller for snps and small indels.

merge_and_call_individual_gvcf
------------------------------
Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.

combine_gvcf
------------
Combine the per sample gvcfs of haplotype caller into one main file for all sample.

merge_and_call_combined_gvcf
----------------------------
Merges the combined gvcfs and also generates a general vcf containing genotypes.

variant_recalibrator
--------------------
GATK VariantRecalibrator.
The purpose of the variant recalibrator is to assign a well-calibrated probability to each variant call in a call set.
You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
between SNP call annotations (QD, MQ, HaplotypeScore, and ReadPosRankSum, for example) and the probability that a SNP
is a true genetic variant versus a sequencing or data processing artifact. This model is determined adaptively based
on "true sites" provided as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP
chip array. This adaptive error model can then be applied to both known and novel variation discovered in the call set
of interest to evaluate the probability that each call is real. The score that gets added to the INFO field of each variant
is called the VQSLOD. It is the log odds ratio of being a true variant versus being false under the trained Gaussian mixture model.
Using the tranche file generated by the previous step the ApplyRecalibration walker looks at each variant's VQSLOD value
and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
have their filter field annotated with its tranche level. This will result in a call set that simultaneously is filtered
to the desired level but also has the information necessary to pull out more variants for a higher sensitivity but a
slightly lower quality level.

haplotype_caller_decompose_and_normalize
----------------------------------------
haplotype_caller_flag_mappability
---------------------------------
Mappability annotation applied to haplotype caller vcf.
An in-house database identifies regions in which reads are confidently mapped
to the reference genome.

haplotype_caller_snp_id_annotation
----------------------------------
dbSNP annotation applied to haplotype caller vcf.
The .vcf files are annotated for dbSNP using the software SnpSift (from the [SnpEff suite](http://snpeff.sourceforge.net/)).

haplotype_caller_snp_effect
---------------------------
Variant effect annotation applied to haplotype caller vcf.
The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

haplotype_caller_dbnsfp_annotation
----------------------------------
Additional SVN annotations applied to haplotype caller vcf.
Provides extra information about SVN by using numerous published databases.
Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
(SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
and other function annotations).

haplotype_caller_gemini_annotations
-----------------------------------
metrics_dna_picard_metrics
--------------------------
metrics_dna_sample_qualimap
---------------------------
metrics_dna_fastqc
------------------
picard_calculate_hs_metrics
---------------------------
Compute on target percent of hybridisation based capture.

metrics
-------
Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
for easy visualization of coverage in the IGV browser.

gatk_callable_loci
------------------
Computes the callable region or the genome as a bed track.

extract_common_snp_freq
-----------------------
Extracts allele frequencies of possible variants accross the genome.

baf_plot
--------
Plots DepthRatio and B allele frequency of previously extracted alleles.

run_multiqc
-----------
cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format.
Using this function will include the orginal final bam file into the  removable file list.

sym_link_fastq
--------------
:return:

sym_link_final_bam
------------------
metrics_ngscheckmate
--------------------
NGSCheckMate is a software package for identifying next generation sequencing (NGS) data files from the same individual.
It analyzes various types of NGS data files including (but not limited to) whole genome sequencing (WGS), whole exome
sequencing (WES), RNA-seq, ChIP-seq, and targeted sequencing of various depths. Data types can be mixed (e.g. WES and
RNA-seq, or RNA-seq and ChIP-seq). It takes BAM (reads aligned to the genome), VCF (variants) or FASTQ (unaligned reads)
files as input. NGSCheckMate uses depth-dependent correlation models of allele fractions of known single-nucleotide
polymorphisms (SNPs) to identify samples from the same individual.
input: file containing all vcfs in project
output:

metrics_verify_bam_id
---------------------
:param self:
:return:

metrics_vcftools_missing_indiv
------------------------------
vcftools: --missing_indv: Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
input: bgzipped vcf file
ouput: missingness flat file

metrics_vcftools_depth_indiv
----------------------------
vcftools: --depth: Generates a file containing the mean depth per individual. This file has the suffix ".idepth".
input: bgzipped vcf file
ouput: idepth flat file

metrics_gatk_sample_fingerprint
-------------------------------
		CheckFingerprint (Picard)
        Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) against a set of known genotypes in the supplied genotype file (in VCF format).
        input: sample SAM/BAM or VCF
        output: fingerprint file

metrics_gatk_cluster_fingerprint
--------------------------------
CheckFingerprint (Picard)
Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) against a set of known genotypes in the supplied genotype file (in VCF format).
input: sample SAM/BAM or VCF
output: fingerprint file

rawmpileup
----------
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

rawmpileup_cat
--------------
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

snp_and_indel_bcf
-----------------
Mpileup and Variant calling. Variants (SNPs and INDELs) are called using
[SAMtools](http://samtools.sourceforge.net/) mpileup.
bcftools view is used to produce binary bcf files.

merge_filter_bcf
----------------
bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
The output of bcftools is fed to varfilter, which does an additional filtering of the variants
and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
for all samples in the experiment.

mpileup_decompose_and_normalize
-------------------------------
mpileup_flag_mappability
------------------------
Mappability annotation applied to mpileup vcf.
An in-house database identifies regions in which reads are confidently mapped
to the reference genome.

mpileup_snp_id_annotation
-------------------------
dbSNP annotation applied to mpileyp vcf.
The .vcf files are annotated for dbSNP using the software SnpSift (from the [SnpEff suite](http://snpeff.sourceforge.net/)).

mpileup_snp_effect
------------------
Variant effect annotation applied to mpileup vcf.
The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

mpileup_dbnsfp_annotation
-------------------------
Additional SVN annotations applied to mpileup vcf.
Provides extra information about SVN by using numerous published databases.
Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
(SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
and other function annotations).

mpileup_gemini_annotations
--------------------------
mpileup_metrics_vcf_stats
-------------------------
Metrics SNV applied to mpileup caller vcf.
Multiple metrics associated to annotations and effect prediction are generated at this step:
change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect,
counts by genomic region, SNV quality, coverage, InDel lengths, base changes,  transition-transversion rates,
summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

metrics_dna_sambamba_flagstat
-----------------------------
delly_call_filter
-----------------
Delly2 is an integrated structural variant prediction method that can
discover, genotype and visualize deletions, tandem duplications, inversions and translocations
at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends
and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.
Structural variants can be visualized using Delly-maze and Delly-suave.
Input: normal and tumor final bams
Returns: bcf file

delly_sv_annotation
-------------------
manta_sv_calls
--------------
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
single efficient workflow.
Returns:Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format.

manta_sv_annotation
-------------------
lumpy_paired_sv
---------------
A probabilistic framework for structural variant discovery.
Lumpy traditional with paired ends and split reads on tumor normal pair.
Returns: bams.

lumpy_sv_annotation
-------------------
wham_call_sv
------------
Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
to employ SVs in association testing.
Returns: vcf.

wham_sv_annotation
------------------
cnvkit_batch
------------

cnvkit_sv_annotation
--------------------
run_breakseq2
-------------
BreakSeq2: Ultrafast and accurate nucleotide-resolution analysis of structural variants.

ensemble_metasv
---------------

metasv_sv_annotation
--------------------

