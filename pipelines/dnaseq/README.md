usage: dnaseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
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
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- bwa_mem_picard_sort_sam
5- picard_merge_sam_files
6- gatk_indel_realigner
7- merge_realigned
8- fix_mate_by_coordinate
9- picard_mark_duplicates
10- recalibration
11- metrics
12- picard_calculate_hs_metrics
13- gatk_callable_loci
14- extract_common_snp_freq
15- baf_plot
16- gatk_haplotype_caller
17- merge_and_call_individual_gvcf
18- combine_gvcf
19- merge_and_call_combined_gvcf
20- variant_recalibrator
21- dna_sample_metrics
22- haplotype_caller_filter_nstretches
23- haplotype_caller_flag_mappability
24- haplotype_caller_snp_id_annotation
25- haplotype_caller_snp_effect
26- haplotype_caller_dbnsfp_annotation
27- haplotype_caller_metrics_vcf_stats
28- haplotype_caller_metrics_snv_graph_metrics
29- rawmpileup
30- rawmpileup_cat
31- snp_and_indel_bcf
32- merge_filter_bcf
33- mpileup_filter_nstretches
34- mpileup_flag_mappability
35- mpileup_snp_id_annotation
36- mpileup_snp_effect
37- mpileup_dbnsfp_annotation
38- mpileup_metrics_vcf_stats
39- mpileup_metrics_snv_graph_metrics
40- verify_bam_id
