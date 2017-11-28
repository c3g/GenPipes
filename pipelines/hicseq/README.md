usage: hicseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch,daemon}] [-f] [--report]
                 [--clean] [-l {debug,info,warning,error,critical}] -e
                 {DpnII,HindIII,NcoI,MboI} [-t {hic,capture}] [-r READSETS]
                 [-v]

Version: 3.0.0

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
  -j {pbs,batch,daemon}, --job-scheduler {pbs,batch,daemon}
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
