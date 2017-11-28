usage: rnaseq_denovo_assembly.py [-h] [--help] [-c CONFIG [CONFIG ...]]
                                 [-s STEPS] [-o OUTPUT_DIR]
                                 [-j {pbs,batch,daemon}] [-f] [--report]
                                 [--clean]
                                 [-l {debug,info,warning,error,critical}]
                                 [-d DESIGN] [-r READSETS] [-v]

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
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- insilico_read_normalization_readsets
5- insilico_read_normalization_all
6- trinity
7- exonerate_fastasplit
8- blastx_trinity_uniprot
9- blastx_trinity_uniprot_merge
10- transdecoder
11- hmmer
12- rnammer_transcriptome
13- blastp_transdecoder_uniprot
14- signalp
15- tmhmm
16- trinotate
17- align_and_estimate_abundance_prep_reference
18- align_and_estimate_abundance
19- gq_seq_utils_exploratory_analysis_rnaseq_denovo
20- differential_expression
21- filter_annotated_components
22- gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
23- differential_expression_filtered
