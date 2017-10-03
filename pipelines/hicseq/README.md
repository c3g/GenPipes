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
  -j {pbs,batch}, --job-scheduler {pbs,batch}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
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
  -e {DpnII,HindIII,NcoI,MboI}, --enzyme {DpnII,HindIII,NcoI,MboI}
                        Restriction Enzyme used to generate Hi-C library
  -t {hic,capture}, --type {hic,capture}
                        Hi-C experiment type
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------

----
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
12- identify_TADs
13- identify_peaks
14- create_hic_file
15- multiqc_report
----
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
16- multiqc_report


```
1- samtools_bam_sort
----------------------
Sorts bam by readname prior to picard_sam_to_fastq step in order to minimize memory consumption.
If bam file is small and the memory requirements are reasonable, this step can be skipped.


2- picard_sam_to_fastq
----------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

3- trimmomatic
--------------
Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
only Adapter1 is used and left unchanged.

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

4- merge_trimmomatic_stats
--------------------------
The trim statistics per readset are merged at this step.


5- fastq_readName_Edit
----------------------
Removes the added /1 and /2 by picard's sam_to_fastq transformation to avoid issues with downstream software like HOMER
      

6- hicup_align
--------------------------
Paired-end Hi-C reads are truncated, mapped and filtered using HiCUP. The resulting bam file is filtered for Hi-C artifacts and duplicated reads. It is ready for use as input for downstream analysis.

For more detailed information about the HICUP process visit: [HiCUP] (https://www.bioinformatics.babraham.ac.uk/projects/hicup/overview/)
   

7- samtools_merge_bams
-----------------------
BAM readset files are merged into one file per sample. Merge is done using [samtools](http://samtools.sourceforge.net/). This step takes as input files the aligned bams/sams from the hicup_align step.
    


8- homer_tag_directory
-----------------------
The bam file produced by HiCUP is used to create a tag directory using HOMER for further analysis that includes interaction matrix generation, compartments and identifying significant interactions.

For more detailed information about the HOMER process visit: [HOMER] (http://homer.ucsd.edu/homer/interactions/index.html)
   

9- interaction_matrices_Chr
---------------------------

IntraChromosomal interaction matrices are produced by Homer at resolutions defined in the ini config file and plotted by HiCPlotter.

For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)

For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)

10- interaction_matrices_genome
--------------------------------

Genomewide interaction matrices are produced by Homer at resolutions defined in the ini config file

For more detailed information about the HOMER matrices visit: [HOMER matrices] (http://homer.ucsd.edu/homer/interactions/HiCmatrices.html)

For more detailed information about HiCPlotter visit: [HiCPlotter] (https://github.com/kcakdemir/HiCPlotter)
        

11- identify_compartments
--------------------------

Genomic compartments are identified using Homer at resolutions defined in the ini config file

For more detailed information about the HOMER compartments visit: [HOMER compartments] (http://homer.ucsd.edu/homer/interactions/HiCpca.html)
        

12- identify_TADs
------------------

Topological associating Domains (TADs) are identified using TopDom at resolutions defined in the ini config file.

For more detailed information about the TopDom visit: [TopDom] (https://www.ncbi.nlm.nih.gov/pubmed/26704975)



13- identify_peaks
-------------------

Significant intraChromosomal interactions (peaks) are identified using Homer.

For more detailed information about the Homer peaks visit: [Homer peaks] (http://homer.ucsd.edu/homer/interactions/HiCinteractions.html)
       

14- create_hic_file
-------------------

A .hic file is created per sample in order to visualize in JuiceBox, WashU epigenome browser or as input for other tools.

For more detailed information about the JuiceBox visit: [JuiceBox] (http://www.aidenlab.org/software.html)


15- multiqc_report
-------------------
A quality control report for all samples is generated.

For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)


----------------------------------------- Capture specific steps:


8- create_rmap_file
--------------------

rmap file for Chicago capture analysis is created using the hicup digestion file.


9- create_baitmap_file
-----------------------

baitmap file for Chicago capture analysis is created using the created rmap file and the probe capture bed file.


10- create_design_files
------------------------

design files (nperbin file (.npb), nbaitsperbin file (.nbpb), proxOE file (.poe)) for Chicago capture analysis are created using the rmap file and the baitmap file.


11- create_input_files
-----------------------

input file (sample.chinput) for Chicago capture analysis is created using the rmap file, the baitmap file and the hicup aligned bam.


12- runChicago
---------------

Chicago is run on capture data. Chicago will filter capture hic artifacts and identify significant interactions. It will output data as a bed file and will also output SeqMonk and WashU tracks.
For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago] https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html


13- runChicago_featureOverlap
------------------------------

Runs the feature enrichement of chicago significant interactions.
For more detailed information about the Chicago, including how to interpret the plots, please visit: [Chicago] https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html

14- bait_intersect
-------------------

provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the baits that have significant interactions. 
Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.

15- capture_intersect
----------------------

provided with a bed file, for example a bed of GWAS snps or features of interest, this method returns the lines in the bed file that intersect with the captured ends ("Other Ends") that have significant interactions. 
Input bed must have 4 columns (<chr> <start> <end> <annotation>) and must be tab separated.
