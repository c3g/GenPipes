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
                     [--container {wrapper, singularity} <IMAGE PATH>]
                     [--genpipes_file GENPIPES_FILE] [-p PAIRS] [--profyle]
                     [-t {fastpass,ensemble,sv}] [-r READSETS] [-v]

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
  -p PAIRS, --pairs PAIRS
                        pairs file
  --profyle             adjust deliverables to PROFYLE folder conventions
                        (Default: False)
  -t {fastpass,ensemble,sv}, --type {fastpass,ensemble,sv}
                        Tumor pair analysis type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![tumor_pair fastpass workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_fastpass.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_fastpass.png)
```
fastpass:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_files
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- sambamba_mark_duplicates
8- recalibration
9- manta_sv_calls
10- rawmpileup_panel
11- paired_varscan2_panel
12- merge_varscan2_panel
13- preprocess_vcf_panel
14- snp_effect_panel
15- gemini_annotations_panel
16- conpair_concordance_contamination
17- metrics_dna_picard_metrics
18- metrics_dna_sample_qualimap
19- metrics_dna_fastqc
20- sequenza
21- run_pair_multiqc
22- sym_link_report
23- sym_link_fastq_pair
24- sym_link_panel
----
```
![tumor_pair ensemble workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_ensemble.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_ensemble.png)
```
ensemble:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_files
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- sambamba_mark_duplicates
8- recalibration
9- conpair_concordance_contamination
10- metrics_dna_picard_metrics
11- metrics_dna_sample_qualimap
12- metrics_dna_fastqc
13- sequenza
14- manta_sv_calls
15- strelka2_paired_somatic
16- strelka2_paired_germline
17- strelka2_paired_germline_snpeff
18- purple
19- rawmpileup
20- paired_varscan2
21- merge_varscan2
22- paired_mutect2
23- merge_mutect2
24- vardict_paired
25- merge_filter_paired_vardict
26- ensemble_somatic
27- gatk_variant_annotator_somatic
28- merge_gatk_variant_annotator_somatic
29- ensemble_germline_loh
30- gatk_variant_annotator_germline
31- merge_gatk_variant_annotator_germline
32- cnvkit_batch
33- filter_ensemble_germline
34- filter_ensemble_somatic
35- report_cpsr
36- report_pcgr
37- run_pair_multiqc
38- sym_link_fastq_pair
39- sym_link_final_bam
40- sym_link_report
41- sym_link_ensemble
----
```
![tumor_pair sv workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_sv.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_tumor_pair_sv.png)
```
sv:
1- picard_sam_to_fastq
2- skewer_trimming
3- bwa_mem_sambamba_sort_sam
4- sambamba_merge_sam_files
5- gatk_indel_realigner
6- sambamba_merge_realigned
7- sambamba_mark_duplicates
8- recalibration
9- strelka2_paired_somatic
10- strelka2_paired_germline
11- metrics_dna_picard_metrics
12- sequenza
13- delly_call_filter
14- delly_sv_annotation
15- manta_sv_calls
16- manta_sv_annotation
17- lumpy_paired_sv
18- lumpy_sv_annotation
19- wham_call_sv
20- wham_sv_annotation
21- cnvkit_batch
22- cnvkit_sv_annotation
23- scones
24- svaba_assemble
25- svaba_sv_annotation
26- ensemble_metasv_somatic
27- ensemble_metasv_germline
28- metasv_sv_annotation
29- sym_link_sequenza
30- sym_link_metasv
31- sym_link_delly
32- sym_link_manta
33- sym_link_lumpy
34- sym_link_wham
35- sym_link_cnvkit

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

sambamba_merge_sam_files
------------------------
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

Note: modified to use both normal and tumor bams to reduce FPs around indels.

sambamba_merge_realigned
------------------------
BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

sambamba_mark_duplicates
------------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using Sambamba](http://lomereiter.github.io/sambamba/index.html).

recalibration
-------------
Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
the quality scores in the QUAL field in each read in the output BAM are more accurate in that
the reported quality score is closer to its actual probability of mismatching the reference genome.
Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
and sequence context, and by doing so, provides not only more accurate quality scores but also
more widely dispersed ones.

manta_sv_calls
--------------
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
single efficient workflow.
Returns:Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
in VCF 4.1 format.

rawmpileup_panel
----------------
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

paired_varscan2_panel
---------------------
Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing.

merge_varscan2_panel
--------------------
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

preprocess_vcf_panel
--------------------
Preprocess vcf for loading into a annotation database - gemini : http://gemini.readthedocs.org/en/latest/index.html
Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
vcf FORMAT modification for correct loading into gemini.

snp_effect_panel
----------------
Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
[SnpEff](https://pcingola.github.io/SnpEff/) annotates and predicts the effects of variants on genes (such as amino acid changes).

gemini_annotations_panel
------------------------
Load functionally annotated vcf file into a mysql lite annotation database [Gemini] (http://gemini.readthedocs.org/en/latest/index.html).

conpair_concordance_contamination
---------------------------------
Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification
(= samples coming from the same individual), as well as cross-individual contamination level estimation in
whole-genome and whole-exome sequencing experiments. Importantly, the method of estimates contamination in
the tumor samples not affected by copy number changes and is able to detect contamination levels as low as 0.1%.

metrics_dna_picard_metrics
--------------------------
Runs specific QC metrics on DNA data.
Functions: collect_multiple_metrics, CollectOxoGMetrics and collect_sequencing_artifacts_metrics
[Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html).

metrics_dna_sample_qualimap
---------------------------
QC alignment metrics generated by [Qualimap](http://qualimap.conesalab.org/).

metrics_dna_fastqc
------------------
QCing metrics generated on the read level using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

sequenza
--------
Sequenza is a novel set of tools providing a fast python script to genotype cancer samples,
and an R package to estimate cancer cellularity, ploidy, genome wide copy number profile and infer
for mutated alleles.

run_pair_multiqc
----------------
Aggregate results from bioinformatics analyses across many samples into a single report
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

sym_link_report
---------------
sym_link_fastq_pair
-------------------
sym_link_panel
--------------
Create sym links of panel variants for deliverables to the clients.

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

strelka2_paired_germline_snpeff
-------------------------------
[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
cohorts and somatic variation in tumor/normal sample pairs.
This implementation is optimized for germline calling in cancer pairs.

purple
------
PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.

rawmpileup
----------
Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
One packaged mpileup file is created per sample/chromosome.

paired_varscan2
---------------
Variant calling and somatic mutation/CNV detection for next-generation sequencing data. 
Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
Varscan2 thresholds based on DREAM3 results generated by author see: https://github.com/dkoboldt/varscan/releases
SSC INFO field remove to prevent collison with Samtools output during ensemble.

merge_varscan2
--------------
Merge mpileup files per sample/chromosome into one compressed gzip file per sample.

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
Note: variants are filtered to remove instantance where REF == ALT and REF modified to 'N' when REF is
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

cnvkit_batch
------------
CNVkit is a Python library and command-line software toolkit to infer and visualize copy number from
high-throughput DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and
custom target panels, and short-read sequencing platforms such as Illumina and Ion Torrent.

filter_ensemble_germline
------------------------
Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
the filter on those generated fields.

filter_ensemble_somatic
-----------------------
Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
the filter on those generated fields.

report_cpsr
-----------
Creates a cpsr gremline report (https://sigven.github.io/cpsr/)
input: filtered ensemble gremline vcf
output: html report and addtionalflat files

report_pcgr
-----------
Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
input: filtered ensemble gremline vcf
output: html report and addtionalflat files

sym_link_final_bam
------------------
Create sym link of final bam for delivery of data to clients.

sym_link_ensemble
-----------------
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
cnvkit_sv_annotation
--------------------
scones
------
This step aims to estimate somatic Copy Number Variation using BVAtools and SCoNEs. BVAtools generate the bined Depth ratio values from the
tumor and normal BAM files. SCoNEs is tool to deconvolution the logR signal of the tumor-normal coverage into a mixture of baysian sub-signal
for each copy number state. The result is a set of several deconvolution using  0-7 sub-signal. As each tumor sample is unique the choice of
the best final model (number of sub-signal) needs to be manually evaluated using the log ratio graphical representation.

svaba_assemble
--------------
SvABA - Structural variation and indel analysis by assembly.

svaba_sv_annotation
-------------------
ensemble_metasv_somatic
-----------------------
MetaSV: An accurate and integrative structural-variant caller for next generation sequencing.

ensemble_metasv_germline
------------------------
MetaSV: An accurate and integrative structural-variant caller for next generation sequencing.

metasv_sv_annotation
--------------------
sym_link_sequenza
-----------------
Sym link of sequenza outputs.

sym_link_metasv
---------------
sym_link_delly
--------------
sym_link_manta
--------------
sym_link_lumpy
--------------
sym_link_wham
-------------
sym_link_cnvkit
---------------

