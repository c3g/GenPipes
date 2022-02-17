[TOC]


Hi-C Pipeline
==============

Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques.
This pipeline analyzes both Hi-C experimental data (-t hic) and capture Hi-C data (-t capture).
The Hi-C pipeline, selected using the "-t hic" parameter, starts by trimming adaptors and low quality bases.
It then maps the reads to a reference genome using HiCUP.
HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
Samples from different lanes are merged and a tag directory is created by Homer, which is also used to produce the interaction
matrices and compartments. TopDom is used to predict topologically associating domains (TADs) and homer is used to identify
significant interactions.

The capture Hi-C pipeline, selected using the "-t capture" parameter, starts by trimming adaptors and low quality bases.
It then maps the reads to a reference genome using HiCUP.
HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates.
Samples from different lanes are merged and CHiCAGO is then used to filter capture-specific artifacts and call significant
interactions. This pipeline also identifies enrichement of regulatory features when provided with ChIP-Seq marks. It can also
return bed interctions with significant baits (bait_intersect step) or with captured interactions (capture_intersect step).

An example of the Hi-C report for an analysis on public data (GM12878 Rao. et al.) is available for illustration purpose only:
[Hi-C report](<url>).

[Here](<url>) is more information about Hi-C pipeline that you may find interesting.


Usage
-----
```
#!text

usage: hicseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                 [--no-json] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [--sanity-check]
                 [--container {wrapper, singularity} <IMAGE PATH>]
                 [--genpipes_file GENPIPES_FILE] -e
                 {DpnII,HindIII,NcoI,MboI,Arima} [-t {hic,capture}]
                 [-r READSETS] [-v]

Version: 4.1.2

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
  -e {DpnII,HindIII,NcoI,MboI,Arima}, --enzyme {DpnII,HindIII,NcoI,MboI,Arima}
                        Restriction Enzyme used to generate Hi-C library
                        (default DpnII)
  -t {hic,capture}, --type {hic,capture}
                        Hi-C experiment type (default hic)
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
```
![hicseq hic workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_hicseq_hic.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_hicseq_hic.png)
```
hic:
1- samtools_bam_sort
2- picard_sam_to_fastq
3- trimmomatic
4- merge_trimmomatic_stats
5- fastq_readName_Edit
6- hicup_align
7- samtools_merge_bams
8- homer_tag_directory
9- interaction_matrices_Chr
10- interaction_matrices_genome
11- identify_compartments
12- identify_TADs_TopDom
13- identify_TADs_RobusTAD
14- identify_peaks
15- create_hic_file
16- reproducibility_scores
17- quality_scores
18- cram_output
19- multiqc_report
----
```
![hicseq capture workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_hicseq_capture.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_hicseq_capture.png)
```
capture:
1- samtools_bam_sort
2- picard_sam_to_fastq
3- trimmomatic
4- merge_trimmomatic_stats
5- fastq_readName_Edit
6- hicup_align
7- samtools_merge_bams
8- create_rmap_file
9- create_baitmap_file
10- create_design_files
11- create_input_files
12- runChicago
13- runChicago_featureOverlap
14- bait_intersect
15- capture_intersect
16- create_hic_file
17- multiqc_report
18- cram_output

```
samtools_bam_sort
-----------------
Sorts bam by readname prior to picard_sam_to_fastq step in order to minimize memory consumption.
If bam file is small and the memory requirements are reasonable, this step can be skipped.

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

fastq_readName_Edit
-------------------
Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER

hicup_align
-----------
Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and
duplicated reads. It is ready for use as input for downstream analysis.

For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)

samtools_merge_bams
-------------------
BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/).

This step takes as input files the aligned bams/sams from the hicup_align step

homer_tag_directory
-------------------
The bam file produced by HiCUP is used to create a tag directory using HOMER for further analysis that includes interaction matrix generation,
compartments and identifying significant interactions.

For more detailed information about the HOMER process visit: [HOMER] (http://homer.ucsd.edu/homer/interactions/index.html)

interaction_matrices_Chr
------------------------
IntraChromosomal interaction matrices are produced by Homer at resolutions defined in the ini config file and plotted by HiCPlotter.
For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)
For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)

interaction_matrices_genome
---------------------------
Genomewide interaction matrices are produced by Homer at resolutions defined in the ini config file
For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)
For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)

identify_compartments
---------------------
Genomic compartments are identified using Homer at resolutions defined in the ini config file
For more detailed information about the HOMER compartments visit: [HOMER compartments] (http://homer.ucsd.edu/homer/interactions/HiCpca.html)

identify_TADs_TopDom
--------------------
Topological associating Domains (TADs) are identified using TopDom at resolutions defined in the ini config file.
For more detailed information about the TopDom visit: [TopDom] (https://www.ncbi.nlm.nih.gov/pubmed/26704975)

identify_TADs_RobusTAD
----------------------
Topological associating Domain (TAD) scores are calculated using RobusTAD for every bin in the genome.
RobusTAD is resolution-independant and will use the first resolution in "resolution_TADs"  under [identify_TADs] in the ini file.
For more detailed information about the RobusTAD visit: [RobusTAD] (https://github.com/rdali/RobusTAD)

identify_peaks
--------------
Significant intraChromosomal interactions (peaks) are identified using Homer.
For more detailed information about the Homer peaks visit: [Homer peaks] (http://homer.ucsd.edu/homer/interactions/HiCinteractions.html)

create_hic_file
---------------
A .hic file is created per sample in order to visualize in JuiceBox, WashU epigenome browser or as input for other tools.
For more detailed information about the JuiceBox visit: [JuiceBox] (http://www.aidenlab.org/software.html)

reproducibility_scores
----------------------
hic-rep is a R package for calculating the inter-chromosomal reproducibility score.
Pairwise reproducibility scores for each chromosome pair in each sample pair are calculated using
hic-rep at resolutions (bin size) defined in interaction_matrices_Chr step
and other parameters defined in reproducibility_scores step of ini config file. All the scores are finally merged
together and output a csv file with parameters used in analysis with chromosome number, reproducibility scores,
standard deviation and smoothing value used for the analysis (in order to compare samples, smoothing value
and the sequencing depth should be similar across samples).
Down-sampling of samples can be performed using the down_sampling parameter in the ini config file.
Correlation matrices and weight matrices can be saved using  corr=TRUE and weights=TRUE in ini config file
for more information visit: [https://bioconductor.org/packages/release/bioc/html/hicrep.html]

quality_scores
--------------
Quality score per chromosome for each sample is calculated using QUASAR-QC at all resolutions
and sequencing depths (coverages) and down_sampling value (coverage) defined in quality_scores step of ini config file
QUASAR-QC is a part of the hifive hic-seq analysis suite
for more information visit: [http://hifive.docs.taylorlab.org/en/latest/quasar_scoring.html]

cram_output
-----------
Generate long term storage version of the final alignment files in CRAM format
Using this function will include the orginal final bam file into the  removable file list

multiqc_report
--------------
A quality control report for all samples is generated.
For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)

create_rmap_file
----------------
rmap file for Chicago capture analysis is created using the hicup digestion file.

create_baitmap_file
-------------------
baitmap file for Chicago capture analysis is created using the created rmap file and the probe capture bed file.

create_design_files
-------------------
design files (nperbin file (.npb), nbaitsperbin file (.nbpb), proxOE file (.poe)) for Chicago capture analysis are created using the rmap file and the baitmap file.

create_input_files
------------------
input file (sample.chinput) for Chicago capture analysis is created using the rmap file, the baitmap file and the hicup aligned bam.

runChicago
----------
Chicago is run on capture data. Chicago will filter capture hic artifacts and identify significant interactions. It will output data as a bed file and will also output SeqMonk and WashU tracks.
For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago] https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html

runChicago_featureOverlap
-------------------------
Runs the feature enrichement of chicago significant interactions.
For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago] https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html

bait_intersect
--------------
provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the baits that have significant interactions.
Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.

capture_intersect
-----------------
provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the captured ends ("Other Ends") that have significant interactions.
Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.


