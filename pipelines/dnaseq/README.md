[TOC]



Usage
-----
```
#!text

usage: dnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                 [-t {mugqic,mpileup,light}] [-r READSETS] [-v]

Version: 3.1.5

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

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
  --container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}
                        run pipeline inside a container providing a container
                        image path or accessible docker/singularity hub path
  -t {mugqic,mpileup,light}, --type {mugqic,mpileup,light}
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
2- sym_link_fastq
3- trimmomatic
4- merge_trimmomatic_stats
5- skewer_trimming
6- bwa_mem_picard_sort_sam
7- sambamba_merge_sam_files
8- gatk_indel_realigner
9- sambamba_merge_realigned
10- fix_mate_by_coordinate
11- picard_mark_duplicates
12- recalibration
13- sym_link_final_bam
14- metrics_dna_picard_metrics
15- metrics_dna_sample_qualimap
16- metrics_dna_sambamba_flagstat
17- metrics_dna_fastqc
18- picard_calculate_hs_metrics
19- gatk_callable_loci
20- extract_common_snp_freq
21- baf_plot
22- gatk_haplotype_caller
23- merge_and_call_individual_gvcf
24- combine_gvcf
25- merge_and_call_combined_gvcf
26- variant_recalibrator
27- haplotype_caller_decompose_and_normalize
28- haplotype_caller_flag_mappability
29- haplotype_caller_snp_id_annotation
30- haplotype_caller_snp_effect
31- haplotype_caller_dbnsfp_annotation
32- haplotype_caller_gemini_annotations
33- haplotype_caller_metrics_vcf_stats
34- run_multiqc
35- cram_output
----
```
![dnaseq mpileup workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mpileup.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_mpileup.png)
```
mpileup:
1- picard_sam_to_fastq
2- sym_link_fastq
3- trimmomatic
4- merge_trimmomatic_stats
5- skewer_trimming
6- bwa_mem_picard_sort_sam
7- sambamba_merge_sam_files
8- gatk_indel_realigner
9- sambamba_merge_realigned
10- fix_mate_by_coordinate
11- picard_mark_duplicates
12- recalibration
13- sym_link_final_bam
14- metrics_dna_picard_metrics
15- metrics_dna_sample_qualimap
16- metrics_dna_sambamba_flagstat
17- metrics_dna_fastqc
18- picard_calculate_hs_metrics
19- gatk_callable_loci
20- extract_common_snp_freq
21- baf_plot
22- rawmpileup
23- rawmpileup_cat
24- snp_and_indel_bcf
25- merge_filter_bcf
26- mpileup_decompose_and_normalize
27- mpileup_flag_mappability
28- mpileup_snp_id_annotation
29- mpileup_snp_effect
30- mpileup_dbnsfp_annotation
31- mpileup_gemini_annotations
32- mpileup_metrics_vcf_stats
33- run_multiqc
34- cram_output
----
```
![dnaseq light workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_light.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_dnaseq_light.png)
```
light:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_picard_sort_sam
4- sambamba_merge_sam_files
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- sambamba_mark_duplicates
8- metrics_dna_picard_metrics
9- metrics_dna_sample_qualimap
10- metrics_dna_sambamba_flagstat
11- metrics_dna_fastqc
12- picard_calculate_hs_metrics
13- gatk_callable_loci
14- extract_common_snp_freq
15- baf_plot
16- gatk_haplotype_caller
17- merge_and_call_individual_gvcf
18- combine_gvcf
19- merge_and_call_combined_gvcf
20- variant_recalibrator
21- haplotype_caller_decompose_and_normalize
22- haplotype_caller_flag_mappability
23- haplotype_caller_snp_id_annotation
24- haplotype_caller_snp_effect
25- haplotype_caller_dbnsfp_annotation
26- haplotype_caller_gemini_annotations
27- run_multiqc
28- cram_output

```
picard_sam_to_fastq
-------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

sym_link_fastq
--------------

:return:

trimmomatic
-----------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
only Adapter1 is used and left unchanged.

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

merge_trimmomatic_stats
-----------------------
The trim statistics per readset are merged at this step.

skewer_trimming
---------------

bwa_mem_picard_sort_sam
-----------------------
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

sambamba_merge_sam_files
------------------------
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

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

fix_mate_by_coordinate
----------------------
Fix the read mates. Once local regions are realigned, the read mate coordinates of the aligned reads
need to be recalculated since the reads are realigned at positions that differ from their original alignment.
Fixing the read mate positions is done using [BVATools](https://bitbucket.org/mugqic/bvatools).

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

sym_link_final_bam
------------------
metrics_dna_picard_metrics
--------------------------
metrics_dna_sample_qualimap
---------------------------
metrics_dna_sambamba_flagstat
-----------------------------
metrics_dna_fastqc
------------------
picard_calculate_hs_metrics
---------------------------
Compute on target percent of hybridisation based capture.

gatk_callable_loci
------------------
Computes the callable region or the genome as a bed track.

extract_common_snp_freq
-----------------------
Extracts allele frequencies of possible variants accross the genome.

baf_plot
--------
Plots DepthRatio and B allele frequency of previously extracted alleles.

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
haplotype_caller_metrics_vcf_stats
----------------------------------
Metrics SNV applied to haplotype caller vcf.
Multiple metrics associated to annotations and effect prediction are generated at this step:
change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect,
counts by genomic region, SNV quality, coverage, InDel lengths, base changes,  transition-transversion rates,
summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.

run_multiqc
-----------
cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list 

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
[SAMtools](http://samtools.sourceforge.net/) mpileup. bcftools view is used to produce binary bcf files.

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

sambamba_mark_duplicates
------------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).


