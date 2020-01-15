[TOC]


Tumor Pair Pipeline
=================

The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline with the exception of the
indel realignment (IR) step. In the tumor pipeline the IR step utilizes both the normal and tumor bam to further reduce
false positives (FPs) in and around indels. The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step.
At this point, a paired caller is used to call SNVs and Indels from the pairs given as input. Additional, muliple cancer callers
are utilized using an ensemble approach and SNVs and Indels seen in at least 2 different callers are retained for further
investigation.

Example command:
python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv

-c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin

-r readset: derived from GQ lims or made yourself. See : https://bitbucket.org/mugqic/mugqic_pipelines#markdown-header-readset-file

-p pairs : format - patient_name,normal_sample_name,tumor_sample_name


Usage
-----
```
#!text

usage: tumor_pair.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                     [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                     [--no-json] [--report] [--clean]
                     [-l {debug,info,warning,error,critical}] [--sanity-check]
                     [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                     [-p PAIRS] [-t {mugqic,mpileup,light}] [-r READSETS] [-v]

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
  -p PAIRS, --pairs PAIRS
                        pairs file
  -t {mugqic,mpileup,light}, --type {mugqic,mpileup,light}
                        DNAseq analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- bwa_mem_picard_sort_sam
5- sambamba_merge_sam_files
6- gatk_indel_realigner
7- sambamba_merge_realigned
8- sambamba_mark_duplicates
9- recalibration
10- conpair_concordance_contamination
11- rawmpileup_panel
12- paired_varscan2_panel
13- merge_varscan2_panel
14- preprocess_vcf_panel
15- snp_effect_panel
16- gemini_annotations_panel
17- metrics
18- picard_calculate_hs_metrics
19- gatk_callable_loci
20- extract_common_snp_freq
21- baf_plot
22- rawmpileup
23- paired_varscan2
24- merge_varscan2
25- paired_mutect2
26- merge_mutect2
27- samtools_paired
28- merge_filter_paired_samtools
29- vardict_paired
30- merge_filter_paired_vardict
31- ensemble_somatic
32- gatk_variant_annotator_somatic
33- merge_gatk_variant_annotator_somatic
34- compute_cancer_effects_somatic
35- combine_tumor_pairs_somatic
36- all_pairs_compute_effects_somatic
37- gemini_annotations_somatic
38- ensemble_germline_loh
39- gatk_variant_annotator_germline
40- merge_gatk_variant_annotator_germline
41- compute_cancer_effects_germline
42- combine_tumor_pairs_germline
43- all_pairs_compute_effects_germline
44- gemini_annotations_germline
45- cram_output

```
picard_sam_to_fastq
-------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

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

Note: modified to use both normal and tumor bams to reduce FPs around indels


sambamba_merge_realigned
------------------------
BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

sambamba_mark_duplicates
------------------------
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

conpair_concordance_contamination
---------------------------------


rawmpileup_panel
----------------
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

paired_varscan2_panel
---------------------
Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing

merge_varscan2_panel
--------------------
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

preprocess_vcf_panel
--------------------
Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
vcf FORMAT modification for correct loading into gemini

snp_effect_panel
----------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).

gemini_annotations_panel
------------------------
Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html

metrics
-------
Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
for easy visualization of coverage in the IGV browser.

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

rawmpileup
----------
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

paired_varscan2
---------------
Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
SSC INFO field remove to prevent collison with Samtools output during ensemble

merge_varscan2
--------------
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

paired_mutect2
--------------
GATK MuTect2 caller for SNVs and Indels.

merge_mutect2
-------------
Merge SNVs and indels for mutect2
Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
Generate a somatic vcf containing only PASS variants

samtools_paired
---------------
Samtools caller for SNVs and Indels using verison 0.1.19.

merge_filter_paired_samtools
----------------------------
bcftools is used to merge the raw binary variants files created in the snpAndIndelBCF step.
The output of bcftools is fed to varfilter, which does an additional filtering of the variants
and transforms the output into the VCF (.vcf) format. One vcf file contain the SNP/INDEL calls
for all samples in the experiment.
Additional somatic filters are performed to reduce the number of FPs:
1. vcflibs vcfsamplediff tags each variant with <tag>={germline,somatic,loh} to specify the type
of variant given the genotype difference between the two samples.
2. bcftools filter is used to retain only variants with CLR>=15 and have STATUS=somatic from
vcfsamplediff
3. bcftools filter is used to retain only variants that have STATUS=germline or STATUS=loh from
vcfsamplediff

vardict_paired
--------------
vardict caller for SNVs and Indels.
Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is AUPAC nomenclature

merge_filter_paired_vardict
---------------------------
The fully merged vcf is filtered using following steps:
1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
2. Somatics identified in step 1 must have PASS filter

ensemble_somatic
----------------
Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls
Filter ensemble calls to retain only calls overlapping 2 or more callers

gatk_variant_annotator_somatic
------------------------------
Add vcf annotations to ensemble vcf: Standard and Somatic annotations

merge_gatk_variant_annotator_somatic
------------------------------------
Merge annotated somatic vcfs

compute_cancer_effects_somatic
------------------------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
Modified arguments to consider paired cancer data.

combine_tumor_pairs_somatic
---------------------------
Combine numerous ensemble vcfs into one vcf for gemini annotations

all_pairs_compute_effects_somatic
---------------------------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
Modified arguments to consider paired cancer data.
Applied to all tumor pairs.

gemini_annotations_somatic
--------------------------
Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html

ensemble_germline_loh
---------------------
Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls
Filter ensemble calls to retain only calls overlapping 2 or more callers

gatk_variant_annotator_germline
-------------------------------
Add vcf annotations to ensemble vcf: most importantly the AD field

merge_gatk_variant_annotator_germline
-------------------------------------
Merge annotated germline and LOH vcfs

compute_cancer_effects_germline
-------------------------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
Modified arguments to consider paired cancer data.

combine_tumor_pairs_germline
----------------------------
Combine numerous ensemble vcfs into one vcf for gemini annotations

all_pairs_compute_effects_germline
----------------------------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
Modified arguments to consider paired cancer data.
Applied to all tumor pairs.

gemini_annotations_germline
---------------------------
Load functionally annotated vcf file into a mysql lite annotation database : http://gemini.readthedocs.org/en/latest/index.html

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list 


