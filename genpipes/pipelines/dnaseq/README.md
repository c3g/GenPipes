<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [DNA-Seq Pipeline](#dna-seq-pipeline)
  - [The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage](#the-pipeline-is-designed-to-be-run-on-a-cluster-and-is-configured-using-a-configuration-file-the-pipeline-can-be-run-in-a-single-step-or-in-multiple-steps-the-pipeline-can-also-be-run-in-parallel-to-process-multiple-samples-simultaneously%0Ausage)
  - [gatk_sam_to_fastq](#gatk_sam_to_fastq)
  - [trim_fastp](#trim_fastp)
  - [bwa_mem2_samtools_sort](#bwa_mem2_samtools_sort)
  - [gatk_mark_duplicates](#gatk_mark_duplicates)
  - [set_interval_list](#set_interval_list)
  - [gatk_haplotype_caller](#gatk_haplotype_caller)
  - [merge_and_call_individual_gvcf](#merge_and_call_individual_gvcf)
  - [combine_gvcf](#combine_gvcf)
  - [merge_and_call_combined_gvcf](#merge_and_call_combined_gvcf)
  - [variant_recalibrator](#variant_recalibrator)
  - [haplotype_caller_decompose_and_normalize](#haplotype_caller_decompose_and_normalize)
  - [haplotype_caller_flag_mappability](#haplotype_caller_flag_mappability)
  - [haplotype_caller_snp_id_annotation](#haplotype_caller_snp_id_annotation)
  - [haplotype_caller_snp_effect](#haplotype_caller_snp_effect)
  - [haplotype_caller_dbnsfp_annotation](#haplotype_caller_dbnsfp_annotation)
  - [haplotype_caller_gemini_annotations](#haplotype_caller_gemini_annotations)
  - [metrics_dna_picard_metrics](#metrics_dna_picard_metrics)
  - [metrics_dna_sample_mosdepth](#metrics_dna_sample_mosdepth)
  - [metrics_picard_calculate_hs](#metrics_picard_calculate_hs)
  - [metrics_verify_bam_id](#metrics_verify_bam_id)
  - [run_multiqc](#run_multiqc)
  - [sym_link_fastq](#sym_link_fastq)
  - [sym_link_final_bam](#sym_link_final_bam)
  - [metrics_vcftools_missing_indiv](#metrics_vcftools_missing_indiv)
  - [metrics_vcftools_depth_indiv](#metrics_vcftools_depth_indiv)
  - [metrics_gatk_sample_fingerprint](#metrics_gatk_sample_fingerprint)
  - [metrics_gatk_cluster_fingerprint](#metrics_gatk_cluster_fingerprint)
  - [delly_call_filter](#delly_call_filter)
  - [delly_sv_annotation](#delly_sv_annotation)
  - [germline_manta](#germline_manta)
  - [manta_sv_annotation](#manta_sv_annotation)
  - [lumpy_paired_sv](#lumpy_paired_sv)
  - [lumpy_sv_annotation](#lumpy_sv_annotation)
  - [wham_call_sv](#wham_call_sv)
  - [wham_sv_annotation](#wham_sv_annotation)
  - [cnvkit_batch](#cnvkit_batch)
  - [cnvkit_sv_annotation](#cnvkit_sv_annotation)
  - [run_breakseq2](#run_breakseq2)
  - [ensemble_metasv](#ensemble_metasv)
  - [metasv_sv_annotation](#metasv_sv_annotation)
  - [samtools_merge_files](#samtools_merge_files)
  - [gatk_fixmate](#gatk_fixmate)
  - [germline_varscan2](#germline_varscan2)
  - [preprocess_vcf](#preprocess_vcf)
  - [snp_effect](#snp_effect)
  - [gemini_annotations](#gemini_annotations)
  - [cram_output](#cram_output)
  - [split_tumor_only](#split_tumor_only)
  - [filter_tumor_only](#filter_tumor_only)
  - [report_cpsr](#report_cpsr)
  - [report_pcgr](#report_pcgr)
  - [sequenza](#sequenza)
  - [rawmpileup](#rawmpileup)
  - [paired_varscan2](#paired_varscan2)
  - [merge_varscan2](#merge_varscan2)
  - [filter_germline](#filter_germline)
  - [filter_somatic](#filter_somatic)
  - [conpair_concordance_contamination](#conpair_concordance_contamination)
  - [sym_link_report](#sym_link_report)
  - [sym_link_fastq_pair](#sym_link_fastq_pair)
  - [sym_link_panel](#sym_link_panel)
  - [manta_sv_calls](#manta_sv_calls)
  - [strelka2_paired_somatic](#strelka2_paired_somatic)
  - [strelka2_paired_germline](#strelka2_paired_germline)
  - [strelka2_paired_snpeff](#strelka2_paired_snpeff)
  - [purple](#purple)
  - [paired_mutect2](#paired_mutect2)
  - [merge_mutect2](#merge_mutect2)
  - [vardict_paired](#vardict_paired)
  - [merge_filter_paired_vardict](#merge_filter_paired_vardict)
  - [ensemble_somatic](#ensemble_somatic)
  - [gatk_variant_annotator_somatic](#gatk_variant_annotator_somatic)
  - [merge_gatk_variant_annotator_somatic](#merge_gatk_variant_annotator_somatic)
  - [ensemble_germline_loh](#ensemble_germline_loh)
  - [gatk_variant_annotator_germline](#gatk_variant_annotator_germline)
  - [merge_gatk_variant_annotator_germline](#merge_gatk_variant_annotator_germline)
  - [report_djerba](#report_djerba)
  - [sym_link_ensemble](#sym_link_ensemble)
  - [gridss_paired_somatic](#gridss_paired_somatic)
  - [purple_sv](#purple_sv)
  - [linx_annotations_somatic](#linx_annotations_somatic)
  - [linx_annotations_germline](#linx_annotations_germline)
  - [linx_plot](#linx_plot)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

[TOC]

DNA-Seq Pipeline
================

A pipeline to process DNA sequencing data. The pipeline uses Trimmomatic for quality control, BWA for alignment to a reference genome, Picard for marking duplicates, GATK for indel realignment and variant calling, SnpEff for variant annotation, SnpSift for filtering variants and MultiQC for aggregate reports.

The pipeline contains protocols for processing both germline and somatic sequencing data; high-coverage data from targeted sequencing experiments; SNVs and structural variants.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.
Usage
-----

```
#!text
usage: genpipes dnaseq [-h] [--clean] -c CONFIG [CONFIG ...]
                       [--container {wrapper, singularity} <IMAGE PATH>] [-f]
                       [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                       [--genpipes_file GENPIPES_FILE]
                       [-j {pbs,batch,daemon,slurm}] [--json-pt]
                       [-l {debug,info,warning,error,critical}]
                       [-o OUTPUT_DIR] [--sanity-check] [-s STEPS]
                       [--wrap [WRAP]] -r READSETS_FILE [-d DESIGN_FILE] [-v]
                       [-p PAIRS] [--profyle]
                       [-t {germline_snv,germline_sv,germline_high_cov,somatic_tumor_only,somatic_fastpass,somatic_ensemble,somatic_sv}]

For more documentation, visit our website: https://genpipes.readthedocs.io

options:
  -h, --help            show this help message and exit
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -c, --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --genpipes_file, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -j, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  -l, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -o, --output-dir OUTPUT_DIR
                        output directory (default: current)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  -s, --steps STEPS     step range e.g. '1-5', '3,6,7', '2,4-8'
  --wrap [WRAP]         Path to the genpipes cvmfs wrapper script. Default is 
                        genpipes/ressources/container/bin/container_wrapper.sh
                        . This is a convenience option for using genpipes in a
                        container
  -r, --readsets READSETS_FILE
                        readset file
  -d, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit
  -p, --pairs PAIRS     pairs file
  --profyle             adjust deliverables to PROFYLE folder conventions
                        (Default: False)
  -t, --type {germline_snv,germline_sv,germline_high_cov,somatic_tumor_only,somatic_fastpass,somatic_ensemble,somatic_sv}
                        DNAseq analysis type

Steps:

Protocol germline_snv
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 set_interval_list
6 gatk_haplotype_caller
7 merge_and_call_individual_gvcf
8 combine_gvcf
9 merge_and_call_combined_gvcf
10 variant_recalibrator
11 haplotype_caller_decompose_and_normalize
12 haplotype_caller_flag_mappability
13 haplotype_caller_snp_id_annotation
14 haplotype_caller_snp_effect
15 haplotype_caller_dbnsfp_annotation
16 haplotype_caller_gemini_annotations
17 metrics_dna_picard_metrics
18 metrics_dna_sample_mosdepth
19 metrics_picard_calculate_hs
20 metrics_verify_bam_id
21 run_multiqc
22 sym_link_fastq
23 sym_link_final_bam
24 metrics_vcftools_missing_indiv
25 metrics_vcftools_depth_indiv
26 metrics_gatk_sample_fingerprint
27 metrics_gatk_cluster_fingerprint

Protocol germline_sv
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 sym_link_final_bam
6 set_interval_list
7 gatk_haplotype_caller
8 merge_and_call_individual_gvcf
9 metrics_dna_picard_metrics
10 metrics_dna_sample_mosdepth
11 metrics_picard_calculate_hs
12 run_multiqc
13 delly_call_filter
14 delly_sv_annotation
15 germline_manta
16 manta_sv_annotation
17 lumpy_paired_sv
18 lumpy_sv_annotation
19 wham_call_sv
20 wham_sv_annotation
21 cnvkit_batch
22 cnvkit_sv_annotation
23 run_breakseq2
24 ensemble_metasv
25 metasv_sv_annotation

Protocol germline_high_cov
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 samtools_merge_files
5 gatk_fixmate
6 metrics_dna_picard_metrics
7 metrics_dna_sample_mosdepth
8 metrics_picard_calculate_hs
9 metrics_verify_bam_id
10 germline_varscan2
11 preprocess_vcf
12 snp_effect
13 gemini_annotations
14 run_multiqc
15 cram_output

Protocol somatic_tumor_only
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 sym_link_final_bam
6 metrics_dna_picard_metrics
7 metrics_dna_sample_mosdepth
8 metrics_picard_calculate_hs
9 metrics_verify_bam_id
10 run_multiqc
11 set_interval_list
12 gatk_haplotype_caller
13 merge_and_call_individual_gvcf
14 combine_gvcf
15 merge_and_call_combined_gvcf
16 variant_recalibrator
17 haplotype_caller_decompose_and_normalize
18 cnvkit_batch
19 split_tumor_only
20 filter_tumor_only
21 report_cpsr
22 report_pcgr

Protocol somatic_fastpass
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 set_interval_list
6 sequenza
7 rawmpileup
8 paired_varscan2
9 merge_varscan2
10 preprocess_vcf
11 cnvkit_batch
12 filter_germline
13 report_cpsr
14 filter_somatic
15 report_pcgr
16 conpair_concordance_contamination
17 metrics_dna_picard_metrics
18 metrics_dna_sample_mosdepth
19 run_multiqc
20 sym_link_report
21 sym_link_fastq_pair
22 sym_link_panel
23 cram_output

Protocol somatic_ensemble
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 set_interval_list
6 conpair_concordance_contamination
7 metrics_dna_picard_metrics
8 metrics_dna_sample_mosdepth
9 sequenza
10 manta_sv_calls
11 strelka2_paired_somatic
12 strelka2_paired_germline
13 strelka2_paired_snpeff
14 purple
15 rawmpileup
16 paired_varscan2
17 merge_varscan2
18 paired_mutect2
19 merge_mutect2
20 vardict_paired
21 merge_filter_paired_vardict
22 ensemble_somatic
23 gatk_variant_annotator_somatic
24 merge_gatk_variant_annotator_somatic
25 ensemble_germline_loh
26 gatk_variant_annotator_germline
27 merge_gatk_variant_annotator_germline
28 cnvkit_batch
29 filter_germline
30 report_cpsr
31 filter_somatic
32 report_pcgr
33 report_djerba
34 run_multiqc
35 sym_link_fastq_pair
36 sym_link_final_bam
37 sym_link_report
38 sym_link_ensemble
39 cram_output

Protocol somatic_sv
1 gatk_sam_to_fastq
2 trim_fastp
3 bwa_mem2_samtools_sort
4 gatk_mark_duplicates
5 set_interval_list
6 manta_sv_calls
7 strelka2_paired_somatic
8 gridss_paired_somatic
9 purple_sv
10 linx_annotations_somatic
11 linx_annotations_germline
12 linx_plot
13 run_multiqc
14 cram_output
```

gatk_sam_to_fastq 
-----------------
 
Converts SAM/BAM files from the input readset file into FASTQ format,
if FASTQ files are not already specified in the readset file. 
Do nothing otherwise.

trim_fastp 
----------
 
[Fastp](https://github.com/OpenGene/fastp): A tool designed to provide fast all-in-one preprocessing for FastQ
files. This tool is developed in C++ with multithreading supported to afford high performance.

bwa_mem2_samtools_sort 
----------------------
 
The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
The alignment software used is [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) with algorithm: bwa mem2.
BWA output BAM files are then sorted by coordinate using [Samtools](https://www.htslib.org/doc/samtools.html)
This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

gatk_mark_duplicates 
--------------------
 
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [GATK](http://broadinstitute.github.io/picard/).

set_interval_list 
-----------------
 
Create an interval list with ScatterIntervalsByNs from GATK: [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360041416072-ScatterIntervalsByNs-Picard).
Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. 
By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering 
the variants with the resulting interval list will be the same as calling with one large region.

gatk_haplotype_caller 
---------------------
 
[GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) haplotype caller for snps and small indels.

merge_and_call_individual_gvcf 
------------------------------
 
Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.

combine_gvcf 
------------
 
Combine the per sample gvcfs of haplotype caller into one main file for all samples.

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
 
Variants with multiple alternate alleles will not be handled correctly by gemini (or by the tools used to annotate the variants).
To reduce the number of false negatives, the authors of gemini strongly recommend that gemini users split, left-align, and trim their variants.
For more info on preprocessing, see the gemini docs: https://gemini.readthedocs.io/en/latest/content/preprocessing.html
The tool used for decomposing and normalizing VCFs is vt: https://github.com/atks/vt

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
The .vcf files are annotated for variant effects using the [SnpEff](http://snpeff.sourceforge.net/) software.
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
 
Load functionally annotated vcf file into a mysql lite annotation database :
http://gemini.readthedocs.org/en/latest/index.html

metrics_dna_picard_metrics 
--------------------------
 
Generates metrics with picard, including:
    [CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard-)
    [CollectOxoGMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037428231-CollectOxoGMetrics-Picard-)
    [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036481572-CollectGcBiasMetrics-Picard-)
    [CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectWgsMetrics-Picard-)
    [CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectHsMetrics-Picard-)
    [CollectInsertSizeMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectInsertSizeMetrics-Picard-)
    [CollectSequencingArtifactMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectSequencingArtifactMetrics-Picard-)
    [CollectQualityYieldMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectQualityYieldMetrics-Picard-)
    [CollectQualityByCycle](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectQualityByCycle-Picard-)
    [CollectBaseDistributionByCycle](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectBaseDistributionByCycle-Picard-)

metrics_dna_sample_mosdepth 
---------------------------
 
Calculate depth stats for captured regions with [Mosdepth](https://github.com/brentp/mosdepth)

metrics_picard_calculate_hs 
---------------------------
 
Compute on target percent of hybridisation based capture with [Picard CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard)

metrics_verify_bam_id 
---------------------
 
[VerifyBamID](https://github.com/Griffan/VerifyBamID) is a software that verifies whether the reads in particular file match previously known
genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
as a mixture of two samples. VerifyBamID can detect sample contamination and swaps when external
genotypes are available. When external genotypes are not available, verifyBamID still robustly
detects sample swaps.

run_multiqc 
-----------
 
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. 
It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.
https://multiqc.info/

sym_link_fastq 
--------------
 
Create sym link of raw reads fastq files.

sym_link_final_bam 
------------------
 
Create sym link of final bam for delivery of data to clients.

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
Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) 
against a set of known genotypes in the supplied genotype file (in VCF format).
input: sample SAM/BAM or VCF
output: fingerprint file

metrics_gatk_cluster_fingerprint 
--------------------------------
 
CheckFingerprint (Picard). Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) 
against a set of known genotypes in the supplied genotype file (in VCF format).
input: sample SAM/BAM or VCF
output: fingerprint file

delly_call_filter 
-----------------
 
Delly2 is an integrated structural variant prediction method that can
discover, genotype and visualize deletions, tandem duplications, inversions and translocations
at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends
and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.
Structural variants can be visualized using Delly-maze and Delly-suave.
Input: normal and tumor final bams
Output: bcl file

delly_sv_annotation 
-------------------
 
Preprocess and annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

germline_manta 
--------------
 
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
single efficient workflow.
Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format.

manta_sv_annotation 
-------------------
 
Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

lumpy_paired_sv 
---------------
 
A probabilistic framework for structural variant discovery.
Lumpy traditional with paired ends and split reads on tumor normal pair.
Outputs: bams

lumpy_sv_annotation 
-------------------
 
Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

wham_call_sv 
------------
 
Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
to employ SVs in association testing.
Outputs: vcf

wham_sv_annotation 
------------------
 
Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

cnvkit_batch 
------------
 
[CNVkit](https://cnvkit.readthedocs.io/en/stable/index.html) is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data.

cnvkit_sv_annotation 
--------------------
 
Annotate VCF with SnpEff.
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
https://pcingola.github.io/SnpEff/se_introduction/

run_breakseq2 
-------------
 
[BreakSeq2](https://bioinform.github.io/breakseq2/): Ultrafast and accurate nucleotide-resolution analysis of structural variants.

ensemble_metasv 
---------------
 
[MetaSV](http://bioinform.github.io/metasv/) is an integrated SV caller which leverages multiple orthogonal SV signals for high accuracy and resolution.
MetaSV proceeds by merging SVs from multiple tools for all types of SVs.

metasv_sv_annotation 
--------------------
 
Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).

samtools_merge_files 
--------------------
 
BAM readset files are merged into one file per sample. Merge is done using [Samtools](https://www.htslib.org/).

This step takes as input files:

1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
2. Else, BAM files from the readset file

gatk_fixmate 
------------
 
Verify mate-pair information between mates and fix if needed.
This ensures that all mate-pair information is in sync between each read and its mate pair.
Fix is done using [Picard](http://broadinstitute.github.io/picard/).

germline_varscan2 
-----------------
 
[VarScan](https://dkoboldt.github.io/varscan/) caller for insertions and deletions.

preprocess_vcf 
--------------
 
Preprocess vcf for loading into an annotation database - Gemini : http://gemini.readthedocs.org/en/latest/index.html
Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
vcf FORMAT modification for correct loading into Gemini.

snp_effect 
----------
 
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
Arguments:
    input_file (str): The input vcf file to annotate for variant effects. Default is allSamples.hc.vqsr.vt.mil.snpId.vcf.gz.
    output (str): The output vcf file. Default is allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.
    job_name (str): The name of the job. Default is snp_effect.hc.

gemini_annotations 
------------------
 
Load functionally annotated vcf file into a mysql lite annotation database :
http://gemini.readthedocs.org/en/latest/index.html
Arguments:
    input_file (str): The input vcf file to load into a mysql lite annotation database. Default is allSamples.merged.flt.vt.mil.snpId.snpeff.dbnsfp.vcf.gz.
    output (str): The output database file. Default is allSamples.gemini.db.
    job_name (str): The name of the job. Default is gemini_annotations.

cram_output 
-----------
 
Generate long term storage version of the final alignment files in CRAM format.
Using this function will add the orginal final bam file to the removable file list.

split_tumor_only 
----------------
 
Splits the merged VCF produced in previous steps to generate a report on a per-patient basis.
The merged VCF is split using the bcftools +split function with the removal of homozygous reference calls.
Creates one VCF per patient to be used for downstream reporting.

filter_tumor_only 
-----------------
 
Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
the filter on those generated fields.

report_cpsr 
-----------
 
Creates a cpsr germline report (https://sigven.github.io/cpsr/)
input: annotated/filter vcf
output: html report and addtional flat files

report_pcgr 
-----------
 
Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
input: filtered somatic vcf
output: html report and addtional flat files

sequenza 
--------
 
Sequenza is a novel set of tools providing a fast Python script to genotype cancer samples,
and an R package to estimate cancer cellularity, ploidy, genome-wide copy number profile and infer
for mutated alleles.

rawmpileup 
----------
 
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

paired_varscan2 
---------------
 
Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
Varscan2 thresholds based on DREAM3 results generated by the author see: https://github.com/dkoboldt/varscan/releases
SSC INFO field remove to prevent collision with Samtools output during ensemble.

merge_varscan2 
--------------
 
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

filter_germline 
---------------
 
Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
the filter on those generated fields.

filter_somatic 
--------------
 
Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
the filter on those generated fields.

conpair_concordance_contamination 
---------------------------------
 
Conpair is a fast and robust method dedicated to human tumor-normal studies to perform concordance verification
(= samples coming from the same individual), as well as cross-individual contamination level estimation in
whole-genome and whole-exome sequencing experiments. Importantly, the method of estimating contamination in
the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.

sym_link_report 
---------------
 
Create a sym link of the MultiQC report for delivery to clients.

sym_link_fastq_pair 
-------------------
 
Create sym links and md5 sums for tumor and normal fastq files.

sym_link_panel 
--------------
 
Create sym links of panel variants for deliverables to the clients.

manta_sv_calls 
--------------
 
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
the analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
single efficient workflow.
Outputs: Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
in VCF 4.1 format.

strelka2_paired_somatic 
-----------------------
 
[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
cohorts and somatic variation in tumor/normal sample pairs
This implementation is optimized for somatic calling.

strelka2_paired_germline 
------------------------
 
[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
cohorts and somatic variation in tumor/normal sample pairs.
This implementation is optimized for germline calling in cancer pairs.

strelka2_paired_snpeff 
----------------------
 
[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
cohorts and somatic variation in tumor/normal sample pairs.
This implementation is optimized for germline calling in cancer pairs.

purple 
------
 
PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.

paired_mutect2 
--------------
 
GATK MuTect2 caller for SNVs and Indels.

merge_mutect2 
-------------
 
Merge SNVs and indels for mutect2.
Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
Generate a somatic vcf containing only PASS variants.

vardict_paired 
--------------
 
vardict caller for SNVs and Indels.
Note: variants are filtered to remove the instance where REF == ALT and REF is modified to 'N' when REF is
AUPAC nomenclature.

merge_filter_paired_vardict 
---------------------------
 
The fully merged vcf is filtered using following steps:
1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
2. Somatics identified in step 1 must have PASS filter

ensemble_somatic 
----------------
 
Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls.
Filter ensemble calls to retain only calls overlapping 2 or more callers.

gatk_variant_annotator_somatic 
------------------------------
 
Add vcf annotations to ensemble vcf: Standard and Somatic annotations.

merge_gatk_variant_annotator_somatic 
------------------------------------
 
Merge annotated somatic vcfs.

ensemble_germline_loh 
---------------------
 
Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls.
Filter ensemble calls to retain only calls overlapping 2 or more callers.

gatk_variant_annotator_germline 
-------------------------------
 
Add vcf annotations to ensemble vcf: most importantly the AD field.

merge_gatk_variant_annotator_germline 
-------------------------------------
 
Merge annotated germline and LOH vcfs.

report_djerba 
-------------
 
Produce Djerba report.

sym_link_ensemble 
-----------------
 
Create sym link of ensemble output for delivery of data to clients.

gridss_paired_somatic 
---------------------
 
Performs joint variant calling on tumor/normal samples using [GRIDSS](https://github.com/PapenfussLab/gridss),
followed by filtering with [GRIPSS](https://github.com/hartwigmedical/hmftools/tree/master/gripss).
GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce 
a high confidence set of somatic SV for a tumor sample. GRIPSS processes the GRIDSS output and produces a somatic vcf.

purple_sv 
---------
 
Runs PURPLE with the optional structural variant input VCFs.
PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.

linx_annotations_somatic 
------------------------
 
[Linx](https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md) is an annotation, interpretation and visualisation tool for structural variants.
The primary function of Linx is grouping together individual SV calls into distinct events 
and properly classify and annotating the event to understand both its mechanism and genomic impact.

linx_annotations_germline 
-------------------------
 
Runs [Linx](https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md) in germline mode.
Linx is an annotation, interpretation and visualisation tool for structural variants.
The primary function of Linx is grouping together individual SV calls into distinct events 
and properly classify and annotating the event to understand both its mechanism and genomic impact.

linx_plot 
---------
 
Generate Linx Plot of the tumor pair analysis.

