usage: ampliconseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                      [-o OUTPUT_DIR] [-j {pbs,batch,daemon}] [-f] [--report]
                      [--clean] [-l {debug,info,warning,error,critical}]
                      [-r READSETS] [-v]

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
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- trimmomatic
2- merge_trimmomatic_stats
3- flash
4- merge_flash_stats
5- catenate
6- uchime
7- merge_uchime_stats
8- otu_picking
9- otu_rep_picking
10- otu_assigning
11- otu_table
12- otu_alignment
13- filter_alignment
14- phylogeny
15- qiime_report
16- multiple_rarefaction
17- alpha_diversity
18- collate_alpha
19- sample_rarefaction_plot
20- qiime_report2
21- single_rarefaction
22- css_normalization
23- rarefaction_plot
24- summarize_taxa
25- plot_taxa
26- plot_heatmap
27- krona
28- plot_to_alpha
29- beta_diversity
30- pcoa
31- pcoa_plot
32- plot_to_beta
