[TOC]

Usage
-----
```
#!text

usage: rnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--json-pt] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {wrapper, singularity} <IMAGE PATH>]
                 [--genpipes_file GENPIPES_FILE] [-t {stringtie,cufflinks}]
                 [-d DESIGN] [-r READSETS] [-v]

Version: 4.6.0

For more documentation, visit our website: https://genpipes.readthedocs.io/en/latest/user_guide/user_guide.html

For source code, visit our bitbucket repository : https://bitbucket.org/mugqic/genpipes/

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
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
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
  -t {stringtie,cufflinks}, --type {stringtie,cufflinks}
                        RNAseq analysis type
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![rnaseq stringtie workflow diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_stringtie.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_stringtie.png)
```
stringtie:
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
<<<<<<< HEAD:pipelines/rnaseq/README.md
4- sortmerna
5- star
6- picard_merge_sam_files
7- picard_sort_sam
8- mark_duplicates
9- picard_rna_metrics
10- estimate_ribosomal_rna
11- rnaseqc2
12- wiggle
13- raw_counts
14- raw_counts_metrics
15- stringtie
16- stringtie_merge
17- stringtie_abund
18- ballgown
19- differential_expression
20- multiqc
21- cram_output
----
```
![rnaseq variants workflow diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_variants.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_variants.png)
=======
4- star
5- picard_merge_sam_files
6- picard_sort_sam
7- picard_mark_duplicates
8- picard_rna_metrics
9- estimate_ribosomal_rna
10- rnaseqc
11- wiggle
12- raw_counts
13- raw_counts_metrics
14- stringtie
15- stringtie_merge
16- stringtie_abund
17- ballgown
18- differential_expression
19- multiqc
20- cram_output
----
```
![rnaseq cufflinks workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_cufflinks.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq_cufflinks.png)
>>>>>>> 9d3195661 (squashing commits on package branch):genpipes/pipelines/rnaseq/README.md
```
cufflinks:
1- picard_sam_to_fastq
2- skewer_trimming
<<<<<<< HEAD:pipelines/rnaseq/README.md
3- sortmerna
4- star
5- picard_merge_sam_files
6- mark_duplicates
7- split_N_trim
8- sambamba_merge_splitNtrim_files
9- gatk_indel_realigner
10- sambamba_merge_realigned
11- recalibration
12- gatk_haplotype_caller
13- merge_hc_vcf
14- run_vcfanno
15- variant_filtration
16- decompose_and_normalize
17- compute_snp_effects
18- gemini_annotations
19- picard_rna_metrics
20- estimate_ribosomal_rna
21- rnaseqc2
22- gatk_callable_loci
23- wiggle
24- multiqc
25- cram_output
=======
3- star
4- picard_merge_sam_files
5- mark_duplicates
6- split_N_trim
7- sambamba_merge_splitNtrim_files
8- gatk_indel_realigner
9- sambamba_merge_realigned
10- recalibration
11- gatk_haplotype_caller
12- merge_hc_vcf
13- run_vcfanno
14- variant_filtration
15- decompose_and_normalize
16- compute_snp_effects
17- gemini_annotations
18- picard_rna_metrics
19- estimate_ribosomal_rna
20- rnaseqc2
21- gatk_callable_loci
22- wiggle
23- multiqc
24- cram_output
>>>>>>> 9d3195661 (squashing commits on package branch):genpipes/pipelines/rnaseq/README.md
----
```
![rnaseq cancer workflow diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_cancer.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_rnaseq_cancer.png)
```
cancer:
1- picard_sam_to_fastq
2- skewer_trimming
<<<<<<< HEAD:pipelines/rnaseq/README.md
3- sortmerna
4- star
5- picard_merge_sam_files
6- mark_duplicates
7- split_N_trim
8- sambamba_merge_splitNtrim_files
9- gatk_indel_realigner
10- sambamba_merge_realigned
11- recalibration
12- gatk_haplotype_caller
13- merge_hc_vcf
14- run_vcfanno
15- decompose_and_normalize
16- filter_gatk
17- report_cpsr
18- report_pcgr
19- run_star_fusion
20- run_arriba
21- run_annofuse
22- picard_rna_metrics
23- estimate_ribosomal_rna
24- rnaseqc2
25- rseqc
26- gatk_callable_loci
27- wiggle
28- multiqc
29- cram_output
=======
3- star
4- picard_merge_sam_files
5- mark_duplicates
6- split_N_trim
7- sambamba_merge_splitNtrim_files
8- gatk_indel_realigner
9- sambamba_merge_realigned
10- recalibration
11- gatk_haplotype_caller
12- merge_hc_vcf
13- run_vcfanno
14- decompose_and_normalize
15- filter_gatk
16- report_cpsr
17- report_pcgr
18- run_star_fusion
19- run_arriba
20- run_annofuse
21- picard_rna_metrics
22- estimate_ribosomal_rna
23- rnaseqc2
24- rseqc
25- gatk_callable_loci
26- wiggle
27- multiqc
28- cram_output
>>>>>>> 9d3195661 (squashing commits on package branch):genpipes/pipelines/rnaseq/README.md

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

sortmerna
---------
Calculation of ribosomal RNA per read based on known ribosomal sequences from archea, bacteria and eukaryotes.
Using [sortmeRNA] (https://github.com/sortmerna/sortmerna)

Taking trimmed fastqs and reporting on each read, either paired-end or single end.

star
----
The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

This step takes as input files:

1. Trimmed FASTQ files if available
2. Else, FASTQ files from the readset file if available
3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

picard_merge_sam_files
----------------------
BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

picard_sort_sam
---------------
The alignment file is reordered (QueryName) using [Picard](http://broadinstitute.github.io/picard/). The QueryName-sorted bam files will be used to determine raw read counts.

picard_mark_duplicates
----------------------
Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
(for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).

picard_rna_metrics
------------------
Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
metrics are collected using [Picard](http://broadinstitute.github.io/picard/).

estimate_ribosomal_rna
----------------------
Use bwa mem to align reads on the rRNA reference fasta and count the number of read mapped
The filtered reads are aligned to a reference fasta file of ribosomal sequence. The alignment is done per sequencing readset.
The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

This step takes as input files:

readset Bam files

<<<<<<< HEAD:pipelines/rnaseq/README.md
rnaseqc2
--------
=======
rnaseqc
-------
>>>>>>> 9d3195661 (squashing commits on package branch):genpipes/pipelines/rnaseq/README.md
Computes a series of quality control metrics using [RNA-SeQC](https://www.broadinstitute.org/cancer/cga/rna-seqc).

wiggle
------
Generate wiggle tracks suitable for multiple browsers.

raw_counts
----------
Count reads in features using [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html).

raw_counts_metrics
------------------
Create rawcount matrix, zip the wiggle tracks and create the saturation plots based on standardized read counts.

stringtie
---------
Assemble transcriptome using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: Still in testing.

stringtie_merge
---------------
Merge assemblies into a master teranscriptome reference using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: still in testing

stringtie_abund
---------------
Assemble transcriptome and compute RNA-seq expression using [stringtie](https://ccb.jhu.edu/software/stringtie/index.shtml).
Warning: Still in testing.

ballgown
--------
[Ballgown](https://bioconductor.org/packages/release/bioc/html/ballgown.html) is used to calculate differential transcript and gene expression levels and test them for significant differences.

Warning: still in testing

differential_expression
-----------------------
Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file.

multiqc
-------
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format.
Using this function will include the orginal final bam file into the  removable file list.

skewer_trimming
---------------
[Skewer](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182) is used mainly for
detection and trimming adapter sequences from raw fastq files. Other features of Skewer is listed
[here](https://github.com/relipmoc/skewer).

split_N_trim
------------
SplitNtrim. A [GATK](https://software.broadinstitute.org/gatk/) tool called SplitNCigarReads developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions.

sambamba_merge_splitNtrim_files
-------------------------------
BAM readset files are merged into one file per sample. Merge is done using [Sambamba] (http://lomereiter.github.io/sambamba/docs/sambamba-merge.html).

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

merge_hc_vcf
------------
Merges vcfs from haplotype caller to generate a sample level vcf

run_vcfanno
-----------
vcfanno is used to annotate VCF files with preferred INFO fields from anu number of VCF or BED files. For more
information [visit](https://github.com/brentp/vcfanno)

variant_filtration
------------------
GATK VariantFiltration.
VariantFiltration is a GATK tool for hard-filtering variant calls based on certain criteria. Records are hard-filtered
by changing the value in the FILTER field to something other than PASS.

cufflinks
---------
Compute RNA-Seq data expression using [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/).
Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format

cuffmerge
---------
Merge assemblies into a master transcriptome reference using [cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/).

cuffquant
---------
Compute expression profiles (abundances.cxb) using [cuffquant](http://cole-trapnell-lab.github.io/cufflinks/cuffquant/).
Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format

<<<<<<< HEAD:pipelines/rnaseq/README.md
gatk_callable_loci
------------------
Computes the callable region or the genome as a bed track.
=======
cuffdiff
--------
[Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/) is used to calculate differential transcript expression levels and test them for significant differences.

cuffnorm
--------
Global normalization of RNA-Seq expression levels using [Cuffnorm](http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/).
>>>>>>> 9d3195661 (squashing commits on package branch):genpipes/pipelines/rnaseq/README.md

fpkm_correlation_matrix
-----------------------
Compute the pearson corrleation matrix of gene and transcripts FPKM. FPKM data are those estimated by cuffnorm.

gq_seq_utils_exploratory_analysis_rnaseq
----------------------------------------
Exploratory analysis using the gqSeqUtils R package.

differential_expression_goseq
-----------------------------
Gene Ontology analysis for RNA-Seq using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
Generates GO annotations for differential gene expression analysis.

ihec_metrics
------------
Generate IHEC's standard metrics.


