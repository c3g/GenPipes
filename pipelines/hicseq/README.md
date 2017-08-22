[TOC]


Hi-C Pipeline
==============

Hi-C experiments allow researchers to understand chromosomal folding and structure using proximity ligation techniques. 
The pipeline starts by trimming adaptors and low quality bases. It then maps the reads to a reference genome using HiCUP. HiCUP first truncates chimeric reads, maps reads to the genome, then it filters Hi-C artifacts and removes read duplicates. Samples from different lanes are merged and a tag directory is created by Homer, which is also used to produce the interaction matrices and compartments. TopDom is used to predict topologically associating domains (TADs) and homer is used to identify significant interactions.

    An example of the Hi-C report for an analysis on public data (GM12878 Rao. et al.) is available for illustration purpose only:
    [Hi-C report](<url>).

    [Here](<url>) is more information about Hi-C pipeline that you may find interesting.


Usage
-----
```
#!text

usage: hicseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS] [-e Restriction_Enzyme {DpnII, MboI, HindIII, NcoI}]
                  [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report] [--clean]
                  [-l {debug,info,warning,error,critical}] [-d DESIGN]
                  [-r READSETS] [-v]

Version: 1.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -e Restriction_Enzyme --enzyme Restriction_Enzyme
                        restriction enzyme used in generating Hi-C library
                        enzymes currently supported: DpnII, MboI, HindIII, NcoI
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
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
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
14- multiqc_report

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
       

14- multiqc_report
-------------------
A quality control report for all samples is generated.

For more detailed information about the MultiQc visit: [MultiQc] (http://multiqc.info/)

