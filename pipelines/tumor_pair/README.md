usage: tumor_pair.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                     [-o OUTPUT_DIR] [-j {pbs,batch,daemon}] [-f] [--report]
                     [--clean] [-l {debug,info,warning,error,critical}]
                     [-p PAIRS] [-r READSETS] [-v]

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
  -p PAIRS, --pairs PAIRS
                        pairs file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
1- picard_sam_to_fastq
2- trimmomatic
3- merge_trimmomatic_stats
4- bwa_mem_picard_sort_sam
5- sambamba_merge_sam_files
6- gatk_indel_realigner
7- sambamba_merge_realigned
8- sambamba_mark_duplicates
9- recalibration
10- conpair_concordance_contamination
11- rawmpileup_panel
12- paired_varscan2_panel
13- merge_varscan2_panel
14- preprocess_vcf_panel
15- snp_effect_panel
16- gemini_annotations_panel
17- metrics
18- picard_calculate_hs_metrics
19- gatk_callable_loci
20- extract_common_snp_freq
21- baf_plot
22- rawmpileup
23- paired_varscan2
24- merge_varscan2
25- paired_mutect2
26- merge_mutect2
27- samtools_paired
28- merge_filter_paired_samtools
29- vardict_paired
30- merge_filter_paired_vardict
31- ensemble_somatic
32- gatk_variant_annotator_somatic
33- merge_gatk_variant_annotator_somatic
34- compute_cancer_effects_somatic
35- combine_tumor_pairs_somatic
36- all_pairs_compute_effects_somatic
37- gemini_annotations_somatic
38- ensemble_germline_loh
39- gatk_variant_annotator_germline
40- merge_gatk_variant_annotator_germline
41- compute_cancer_effects_germline
42- combine_tumor_pairs_germline
43- all_pairs_compute_effects_germline
44- gemini_annotations_germline
