[TOC]

('AmpliconSeq',)

Usage
-----
```
#!text

Protocol default
0 trimmomatic16S
1 merge_trimmomatic_stats16S
2 flash_pass1
3 ampliconLengthParser
4 flash_pass2
5 merge_flash_stats
6 asva
7 multiqc
```

trimmomatic16S 
--------------
 
MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

This step takes as input files:
1. MiSeq paired-End FASTQ files from the readset file

merge_trimmomatic_stats16S 
--------------------------
 
The trim statistics per readset are merged at this step.

flash_pass1 
-----------
 
Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). Overlapping regions between paired-end reads are found and 
then merged into a continuous strand.

ampliconLengthParser 
--------------------
 
Looks at FLASH output statistics to set input amplicon lengths for dada2. Minimum lengths are set by ensuring that they represent 
at least 1% of the total number of amplicons.

flash_pass2 
-----------
 
Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). The second pass uses statistics obtained from the first pass
to adjust merging.

merge_flash_stats 
-----------------
 
Merges statistics from both flash passes.

asva 
----
 
Checks for design files required for PCA plots, sets up directories, links readset fastq files, and initiates 
[DADA2](https://benjjneb.github.io/dada2/). 

DADA2 is used to infer sequence variants of microbial communities.

multiqc 
-------
 
A quality control report for all samples is generated.
For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)


