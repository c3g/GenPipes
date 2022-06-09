35 tags, 9446 commits

HEAD        Thu Jun 9 13:15:16 2022 -0400        0 commits

4.2.1        Thu Jun 9 17:48:47 2022 +0000        107 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      16 commits

       4e90741 GenPipes - DNASeq : fixing dependencies in recalibratino step for unmapped_reads
       46042e1 GenPipes - remove useless sections in ini
       40d6fe2 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7dff261 GenPipes - move 'sambamba_merge_sam_files' into common.py + renamed 'sambamba_merge_sam_files' to 'sambamba_merge_sam_extract_unmapped' in DNASeq
       34ea28a Cleaning bfx/sambamba from samtools dependencies
       adf5e7f GenPipes - DNASeq : revisited recal jobs to avoid useless restart of the step
       ed49153 GenPipes - DNASeq : Fixing dependencies for recal .bam.bai
       b7ef560 GenPipes - DNASeq : Recreating of the index of the recal bam file after merging of the unmapped reads
       373fb13 GenPipes - DNASeq : correcting the job which merges the unmapped reads to the recalibrated fastqs
       d30127d GenPipes - DNASeq - Fix job queueing for merge unmapped to recalibrated
       1fc5243 fix typo
       afffee2 GenPipes - DNASeq : corrected unmapped_reads job queueing
       379dfdf GenPipes - DNASeq : removing useless part of code...
       4b897f0 GenPipes - DNASeq : extracting unmapped reads from the merged bam, then re-instering them into the recalibrated bam
       460b6a3 Version bump to 4.2.0
       5cc640c Version bump to 4.1.4-beta

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       5690df4 GenPipes - C3G Software Stack : showing only the 3 top versions of the softwares and adding '...' when necessary to indicate more older versions are also available
       34b06e6 GenPipes - C3G Software Stack : reverse-sorting the software versions when combining the metadata JSONs
       416fe4a GenPipes - C3G Software Stack : improved sorting of JSON output
       864465e GenPipes - C3G Software Stack : sort JSON output
       d3f03f6 GenPipes - C3G Software Stack : parametrized the redirection of the output to a file or to stdout
       701298b GenPipes - C3G Software Stack : removed the sending of the JSON
       85b7009 GenPipes - C3G Software Stack : corrected typo
       833d3ca GenPipes - C3G Software Stack : corrected message when file is sent
       3303a59 GenPipes - C3G Software Stack : Added sending of the combined JSON to url + reformating the module helpers
       0b8b7d5 GenPipes - module helper : updating the pythonSearcher class

  ehenrion <edouard.henrion@mcgill.ca>      6 commits

       5f690ff Merged in soft_jsondb_gsoc2020_eh (pull request #181)
       ba84d8c dnaseq.base.ini edited online with Bitbucket
       943e103 Merged in unmapped_reads (pull request #355)
       ef817ea fixed typo
       ecdcb0a dnaseq.base.ini edited online with Bitbucket
       9f80a8c Merged in release_4.2.0 (pull request #354)

  Moonshroom <yatharthrai16@ducic.ac.in>      4 commits

       a51bba0 Some documentation edits, README updates, fixed bad variable naming
       8d223ca Migration to argparse
       c6bde0f Added -h/--help flag
       810321f 1. Pushing pre-final code 2. Added CLI args 3. Added CLI Documentation

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      20 commits

       cd2845a Merged in HotFix_dev (pull request #357)
       a5caedf General - General - Fixing UnboundLocalError: local variable 'job_name_prefix' referenced before assignment
       100567c Merged in HotFix_dev (pull request #356)
       a35809f General - General - Fixing typo
       edc01ae General - General - Fixing typo
       20801c4 General - General - Adding step_wrapper to slurm and batch mode (fixing regression)
       e786922 General - general - Standardizing inis + fixing methylseq.base.ini
       343ec56 Merged in HotFix_dev (pull request #350)
       401be79 ampliconseq - General - Upgrading to latest release of mugqic_tools
       a52c6d6 covseq - General - Upgrading to latest release of covseq_tools
       92e0a5f ampliconseq - flash - Debug previous flash commit
       d73c35d ampliconseq - ampliconLengthParser - Fixing head with pipefail
       c74fe1c ampliconseq - flash - Fixing None appearing in cmd
       720962b covseq - qualimap - Fixing qualimap parameter for threads
       92081f8 covseq - prepare_report - Fixing pipefail issue for ncovtools: we always want ncvotools error to be ignored and continue the job
       e40e7cc covseq - prepare_report - Fixing pipefail issue when grep doesn't find anything
       2c6b46c covseq - prepare_table - Adding threads param
       5387d29 covseq - prepare_table - Adding samtools as module required
       b8cfc5e methylseq - metrics - Changing ini path of file
       f20a986 General - general - Standardizing inis + fixing methylseq.base.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       cf3d727 Merged in release_4.2.1 (pull request #358)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       c7b3aa6 Merge remote-tracking branch 'origin/dev' into release_4.2.1
       b092002 Update readme and version for release
       5469a90 remove -nt completely from indel realigner
       ccf8529 abacus ini typo
       720fd57 typo in tabix option, and no threads in gatk_indel_realigner
       68183a8 always force overwrite on tabix
       ef425e7 tweak and fix base ini

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      39 commits

       79ef518 addressed comments from Paul
       2011bfb removed methylseqraw class
       a726b97 Fix issues in input and output dependencies
       66e55e3 resolved differences in dev branch files
       9c8dd41 Modified dragen.ini file
       c5acace fixed documentation typo
       5949d7b fixed issues in output file tracking added documentation
       fea7a92 fixed issue in dragen hybrid protocol in single-pass
       8fca6d1 added dragen_bedgraph step
       454abe6 added split_dragen_methylation_report
       313ce78 fixed issues in dragen_methylation_call step working pipeline
       e2a952e trying to fix dependencies
       0c48154 added missing job to the dragen_aling
       15dda51 fixed a bug
       d3e9874 Add changes in common and methylseq
       c853189 Implement jsonator disabling at job level
       4d64d73 Debug dragen protocol
       d1c7bf2 update dragen protocol
       8ed0ab8 Implement Modified dependency of dragen command to allo mv commands Added dragen function
       2ce9490 Implement Modified dependency of dragen command to allo mv commands Added dragen function
       dfa349d fixing issues with the done file tracking system
       8612011 Trying to fix dependency issue in dragen
       61945dc resolved merge conflicts
       f22f0f6 resolved merge conflicts
       6664a85 added dragen mark_duplication
       427d7a0 fixed a bug
       33caf90 update dragen protocol
       3cb1965 added missing job to the dragen_aling
       e0e21a7 Resolved incorrect tool link
       93a084b fixed a bug
       30ef76c fixed a bug
       fa534de Add changes in common and methylseq
       e8833f5 Removing debug outputdir log info
       8485ae7 Implement jsonator disabling at job level
       98cc972 Debug dragen protocol
       5ec4bfc update dragen protocol
       c32117b pipelines/methylseq/methylseq.dragen.ini
       5e71366 Implement Modified dependency of dragen command to allo mv commands Added dragen function
       179cabf Implement Modified dependency of dragen command to allo mv commands Modify dependency of bash commands Added bfx/dragen.py

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       d82bab0 added comments and fix hybrid mode single pass directional-complement
       865a65f fixed ihec metric by adding a job to create an empty metric file for estimated_library_size
       cd94f55 changed dragen structure and added a new protocol

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       7faffc5 Merged in dragen_update_version2 (pull request #346)

4.2.0        Wed Jun 1 20:01:30 2022 +0000        230 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      16 commits

       695798b Generating READMEs for 4.2.0
       9947d05 Version bump to 4.1.3
       d7f1942 Merge remote-tracking branch 'origin/dev' into release_4.1.3
       e12f9ed in prep for a release
       e135fce GenPipes - Tumor Pair : corrected location of the coverage_bed in strelka jobs + prepend strelka job with 'rm' command to ensure folder is clean
       b7ed614 GenPipes - Tumor Pair : some code cleaning before correcting strelka errors
       1d34130 GenPipes - DNASseq : updating base.ini
       20453e2 GenPipes - DNASEQ : fixing base.ini for SV protocol
       fe3ccfc GenPipes - Scheduler : Re-using the default python, from ini, to call job2json.py
       61218ac GenPipes - DNASeq SV : fix sambamba_merge_realigned
       fa7cc65 GenPipes - DNASeq SV : using Python2 for breakseq2
       b790c91 GenPipes - BFX : CNVkit now uses the muqgic/CNVkit module, instead of a pyhton module which would
       9db4af2 GenPipes - DNASeq SV : set use of Python2 for MantaSV in base.ini
       97513bd GenPipes - DNASeq SV : fixing resource dependencies
       5685c7b GenPipes - DNASeq SV : fixed delly call
       1ef5a9e Version bump to 4.1.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      4 commits

       da1d31b Resources - install scripts : update to mugqic_tools 2.10.4 + install_modules updates for improved installation in dev space
       5616bb7 GenPipes - Resources : ensuring the use of rpath when patching binaries
       b15769c GenPipes - Resources : adding some new installation scripts
       802c02d GenPipes - Resources : Updating some install scripts with newer version + adeind bioinfokit and lofreq install scripts

  Édouard Henrion <henrione@narval1.narval.calcul.quebec>      3 commits

       13b4a69 GenPipes - Tumor Pair : setting Python2 for Strelka2 in the config file
       69201a3 GenPipes - Tumor Pair : set the output_dirs global variable to better handle outputs, thus dependencies
       a9bf30c GenPipes - Resourcess : adding pycoqc install script + updating other scripts

  Édouard Henrion <henrione@narval2.narval.calcul.quebec>      4 commits

       70c9b89 GenPipes - Tumor Pair : fixed dependency issues by using relative path everywhere in the pipeline, no more absolute path. Also setting 'samples' for all the jobs.
       1156e33 GenPipes - Tumor Pair : updating mugqic_tools to latest version 2.10.2
       2a541de GenPipes - BFX : updated GATK wrappers to remove print_reads outputs before running
       576724b GenPipes - Tumor Pair : fixing Strelka & Manta SV step

  Édouard Henrion <henrione@narval3.narval.calcul.quebec>      1 commits

       bab3a85 GenPipes - Tumor Pair : Fixing conpair with specific call to Python2

  ehenrion <edouard.henrion@mcgill.ca>      16 commits

       927e40e Merged in release_4.2.0 (pull request #353)
       a546428 Merged in release_4.1.3 (pull request #352)
       f1e6716 Merged in release_4.1.3 (pull request #351)
       43e87b7 dnaseq.base.ini edited online with Bitbucket
       cdcebcb Merged in tumor_pair_hotfix_dev (pull request #339)
       0205057 chipseq.base.ini edited online with Bitbucket
       f52b3e2 Merged in tumorpair_hotfix_strelka (pull request #331)
       4db4e81 Merged in dnaseqsv_fix_eh (pull request #329)
       66eca37 Merged in fix_dnaseqsv_eh (pull request #320)
       e634cdd Merged in fix_delly_call_dnaseqsv_eh (pull request #317)
       2c20036 GenPipes - Scheduler : stop using python from ini, use default mugqic python3 instead to ensure job2json is working in al cases
       2ee5d49 GenPipes - COVSEQ : corrected covseq.graham.ini for dna_sample_qualimap
       0b5b896 GenPipes - COVSEQ : corrected covseq.cedar.ini for dna_sample_qualimap
       41b755e GenPipes - COVSEQ : corrected covseq.beluga.ini for dna_sample_qualimap
       2364a6c GenPipes - Readset : corrected "writer" call in checkDuplicateReadets
       a1be2b9 Merged in release_4.1.2 (pull request #314)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       97a0701 Merged in JoseHector-Galvez-Lopez/nanopore_covseqbaseini-edited-online-wit-1645559267177 (pull request #316)
       7c1b1f5 Final adjustment to the nanopore_covseq base ini to optimize job submission on Abacus.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      54 commits

       9c24099 General - general - Adding cpulimit to be used by default on Abacus for all jobs (discussed with PO)
       65971dc Merged in covseq_stretenp (pull request #341)
       a03cb25 Merged in HotFix_dev (pull request #342)
       10b13c0 covseq - prepare_report - Adding pandoc from cvmfs
       9f12d3e General - ini - Updating rnaseq and dnaseq ini
       058e9d0 rnaseq - star_index - Auto-Calculating star genomeSAindexNbases based on genome length
       8e4d4d5 rnaseq - bedtools - Not normalizing if input reads too low
       d1261c0 Merged in HotFix_dev (pull request #340)
       129c709 methylseq - bismark - Hardcoding the multicore value for now
       4c787bd General - ini - Standardizing ini base section
       5c02400 Merged in Paul-Stretenowich/covseqbaseini-edited-online-with-bitbuck-1649436403200 (pull request #336)
       e29820b covseq - ini - Fixing bam2fq number of threads hard coded
       75863a0 Merged in HotFix_dev (pull request #334)
       15069c5 chipseq - macs2_callpeak + homer_annotate - making genomre_size not required in ini
       2456dec Merged in HotFix_dev (pull request #330)
       3f46d5c chipseq - homer_annotate_peaks - Adding way to customize genome size in ini for cit
       47ab87a chipseq - macs2_callpeak - Adding way to customize genome size in ini for cit
       f987882 General - scheduler - Addressing PO's comment
       835b0f5 chipseq - homer_annotate_peaks - Debug genome_size
       2942015 chipseq - Misc - Adding genome size + log to .o file in batch mode
       a57654c chipseq - homer_annotate_peaks + homer_find_motifs_genome - Adding way of using genome fasta and not only uscs naming
       ee79ed5 chipseq - homer_make_ucsc_file_bigWig - Adding ini_section variable + IMPORTANT: adding way to skip exit code 141 due to (z)cat | head sending SIGPIPE despite not being an error
       4c9fbc2 chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       4fbf274 chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       7c972fe General - scheduler - Fixing batch mode job2json
       81ceffe General - scheduler - Fixing batch mode job2json
       484353c General - scheduler - Fixing batch mode
       a9d439c General - scheduler - Fixing batch mode
       8118756 chipseq - base.ini - Fixing cluster_cpu old way
       3d5c6a3 Merged in covseq_nanopore (pull request #325)
       5d8cadb nanopore_covseq - General - Code cleaning
       564fc03 nanopore_covseq - General - Code cleaning
       06c695b covseq - rename_consensus_header - Updating date from 2021 to 2022 by default
       20df473 nanopore_covseq - General - Adding step_wrapper to slurm and batch mode + Adding pipefail in jobs to catch errors in pipes
       5231d85 nanopore_covseq - General - Changing default resources
       9db1946 nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       cba39ac nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       f91a244 nanopore_covseq - jsonator - Adding NanoporeCoVSeq for json Nanopore
       e818ccf nanopore_covseq - General - Upgrading to python3 ini by default
       27453fb nanopore_covseq - prepare_report - Re-implementing cpulimit to be set at ini level and not automatically
       73b5b5b Merged in HotFix_dev (pull request #328)
       2d75e43 dnaseq - General - HotFix
       01c2b49 Merged in chipseq_urgent_fix (pull request #326)
       96a5f32 chipseq - General - Upgrading to python3
       17f337f chipseq - differential_binding - Adding pandoc to module load
       da53a9d Merged in chipseq_urgent_fix (pull request #322)
       97dfd74 Merged in covseq_nanopore (pull request #321)
       0215ccf chipseq - macs2_callpeak - Fixing typo error
       32d6175 nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       c4e2ad5 nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       a8017e3 nanopore_covseq - prepare_report - Adding cpulimit for Abacus
       20f047e Merged in Fixing-Issue-#144 (pull request #315)
       2b3c2b5 chipseq - macs2_callpeak - Fixing other_options to be used via ini
       38d0c9a chipseq - macs2_callpeak - Fixing other_options to be used via ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       868b49d Merged in tp_reports (pull request #348)
       3649818 Merged in tp_reports (pull request #347)
       ea6260f Merged in new_ini (pull request #319)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      20 commits

       98bc236 tabix ovewright existing index in tumor pair
       166b252 fix ram in depth coverage gatk
       433f7a3 fix format2pcgr config
       4230791 up mugqic tool version in tumor pair
       5a9f72c fix compute_effect with cancer
       e42f778 revert R for gq_seq_utils_exploratory_analysis_rnaseq
       4eb6fd6 python 3.10.4 for filter ensemble
       edc5569 remove sed from manta_sv
       de2d425 text file in not compress
       7fabc99 gatk_indel_realigner does not support multithread
       34d2e11 Cleanup dev files
       8e0e9fe typo
       cee0321 typo ini files
       d40a229 let cpulimit do its stuff on abacus
       4c6514f let cpulimit do its stuff on abacus
       0e25a16 fix rnaseq denovo crashes
       3a88478 disable reporting to mgcill.hpc with NO_MUGQIC_REPORT variable set
       2924b41 fix nodes set in pbs
       6abbcc4 more robust log_report.py
       d93e094 new ini file setup with common ini for clusters

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       9ea6e29 fixed a bug in input file list of haplotypecaller
       81cce4e Added new parameters to change the FDR and P-value of the differential binding analysis (chip-seq)
       c7ff997 fixed a bug in input file list of haplotypecaller
       9002783 fixed a bug in input file list of haplotypecaller
       75c58e4 fixed a bug in input file list of haplotypecaller
       4d7ba2c fixed a bug in input file list of haplotypecaller
       c3b27d8 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       193ef5a Added variant calling using GATK4 to chipseq pipeline

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar1.cedar.computecanada.ca>      6 commits

       c5d8d3e fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       7605e7c fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       681a6b8 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       ed358a8 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       66e6e90 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       3e29640 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar5.cedar.computecanada.ca>      5 commits

       40cf8d3 fixed a bug in input file list of haplotypecaller gatk
       5c903d4 fixed a bug in input file list of haplotypecaller gatk
       7417a63 fixed a bug in input file list of haplotypecaller gatk
       2aca866 fixed a bug in input file list of haplotypecaller gatk
       9c2c698 fixed a bug in input file list of haplotypecaller gatk

  Pubudu Nawarathna Mudiyanselage <pubudu@narval1.narval.calcul.quebec>      52 commits

       1510c0f changed job name for merge_gatk
       0442a69 removed duplicated haplotype calling functions
       08af8f1 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       ec3bd09 Modified interval_padding in dnaseq, exome-seq and chipseq
       9f0b3d3 fix padding issue
       3699643 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding resolved merge conflicts
       a40ed3a Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d111f5f fixed a bug in input file list of haplotypecaller gatk
       4c686e8 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       f595d71 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       eeda83b Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       0b0f2a2 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       d84db85 Added variant calling using GATK4 to chipseq pipeline
       abdc216 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       127c4b4 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       7a650a7 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e74ac07 fixed a bug in input file list of haplotypecaller gatk
       af8bb91 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       3e1af78 Added variant calling using GATK4 to chipseq pipeline
       896157b Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       f523e62 Added variant calling using GATK4 to chipseq pipeline
       67f4788 corrected diffBind path
       fc76ea5 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       95592a4 Modified interval_padding in dnaseq, exome-seq and chipseq
       7cee101 fix padding issue
       1e71673 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       102b90a fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       4ca359b Added variant calling using GATK4 to chipseq pipeline
       588b3fb fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       89e4ed4 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       56bcbc6 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       1879f7e fixed a bug in input file list of haplotypecaller gatk
       a33712d fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       faec78c Added variant calling using GATK4 to chipseq pipeline
       d5189ed Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       a120b8f Modified interval_padding in dnaseq, exome-seq and chipseq
       13b93e9 fix padding issue
       6f6e56c Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       f2e5299 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       f5d7e5e fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       1ba29fb Added variant calling using GATK4 to chipseq pipeline
       b2132fa fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       6d18c04 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d4616a3 changed mugqic_tools version to version 2.10
       9e33161 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       7a7a77f Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       960a2c9 fixed a bug in input file list of haplotypecaller gatk
       7143c0d fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       c79776a Added variant calling using GATK4 to chipseq pipeline
       b31f874 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       1331877 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e051b83 changed default python version to python3

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      4 commits

       6174fcb Merged in chipseq_snp_hotfix (pull request #335)
       9213853 Merged in chip_snp_pubudu (pull request #324)
       72ef484 Merged dev into chip_snp_pubudu
       df9593d Merged in chip_seq_issue_150_fix (pull request #327)

  Robert Eveleigh <eveleigh@narval1.narval.calcul.quebec>      12 commits

       b43c90d snpEff cancer pair file fix
       7c4ce71 panel snpeff fix
       853d0b2 merge varscan2 fixes
       88c2b0e add ram options to gatk_interval_list2bed
       5757d2c dependency fixes at various steps
       f3d82fc bwa to indel realignment dependency fix
       85f21e4 changed mugqic tools version
       e5a1235 use purple annotated strelka2 output for somatic ensemble calling
       99b328c update module verions - issues with bcftools
       41b57ce cit mugqic_tools version fix
       87768ac resource fixes: tumor_pair.exome.ini
       8124f35 resource fixes: dnaseq.base.ini, tumor_pair.extras.ini and cit.ini

  Robert Eveleigh <eveleigh@narval2.narval.calcul.quebec>      9 commits

       2c41ddc bcftool -i/-e fix for tumor pair fastpass varscan2 merge
       0c23107 abacus resource fixes for merging and filtering vcf steps
       8a4df12 memory fixes to merge vcf steps
       85860f8 vardict and varscan2 germline filter fixes
       00dfdf3 fixes to conpair for multiqc intergration
       ff85c18 tumor_pair.dev.ini fix
       a5c24ac change file name : tumor.dev.ini to tumor_pair.dev.ini
       75bc9dc reverted back to python2 for strelka2
       35e5c27 adding cpsr/pcgr reporting with addition of manta, cnvkit to ensemble protocol.  Step clean up as well

  Robert Eveleigh <eveleigh@narval3.narval.calcul.quebec>      7 commits

       653c5c6 more germline filter fixes
       dc40100 resource fix to manta and strelka2, and germline filter improvement
       0a3aef5 add snpEff QCs and fixed conpair concordance output
       f0749eb resource fixes
       9632453 add python2 for varscan2 mugqic tool
       6733ac4 varscan2 one job fix - snp/indel merged vcf now in right place
       12726b3 cpsr/pcgr dependency fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       8948532 python3 fixes for support scripts
       f2ea81d conpair and ensemble filtering fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      2 commits

       54b97ee resource fixes
       889346c python3 fixes when running on abacus

4.1.2        Thu Feb 17 17:40:12 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       7997ec5 Merge remote-tracking branch 'origin/dev' into release_4.1.2
       c83a283 Version bump to 4.1.2
       dd9f185 Updating pipeline READMEs before release
       4bcda97 Version bump to 4.1.1

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       e5a64db Merged in release_4.1.2 (pull request #313)
       99b40bf Merged in release_4.1.2 (pull request #312)
       dc6f51e Merged in release_4.1.2 (pull request #311)
       7dd399f GenPipes - Scheduler : Correcting cluster_cpu parsing for PBS
       7a2aa6d Merged in release_4.1.1 (pull request #306)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       c07382b Merged in ont_covseq_report_ini (pull request #308)
       46e25d1 Adjustment to the `prepare_report` step resources in the base ini file, to try to prevent this job getting stuck in the queue.

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       4d43f87 Merged in patch_4.1.2 (pull request #309)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       9895146 remove cluster_submit_cmd from slurm and pbs ini
       df8c9ba fix hic ini
       f45fdba remode empty quotes
       a955e74 cleanup queue
       4417552 PBS gets nodes and ppn together

4.1.1        Thu Feb 10 15:12:29 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      10 commits

       dc3ff3c Merge remote-tracking branch 'origin/dev' into release_4.1.1
       d2a5407 minor README uedit, again...
       3cd7d34 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       21ee999 Another minor correction in the README...
       549f6dc Merge remote-tracking branch 'origin/dev' into release_4.1.1
       bbc6156 Updating READMEs before release
       150f4a0 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       77e0698 Adding missing space in container help message
       be59560 GenPipes - prep for minor release
       11b8ec2 Version bump to Release 4.1.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       aeccdc3 Merged in release_4.1.1 (pull request #305)
       c756a50 GenPipes - Tumor Pair : Replacing "xrange" calls by "range" calls to resolve Python3-compatibility issues
       95f8ce8 README.md edited online with Bitbucket
       691b2fd Merged in release_4.1 (pull request #304)
       0204b7c README.md edited online with Bitbucket
       ea152c3 corrected typo in README.md for Nanopore_covseq
       ae1c04a Merged in release_4.1 (pull request #303)

4.1.0        Mon Feb 7 21:39:25 2022 +0000        77 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      11 commits

       d00f24a GenPipes - Release : creation of the release branch
       03e6919 Merge remote-tracking branch 'origin/dev' into release_4.1
       185575e Correcting typo in READMEs
       ebc18e3 GenPipes - Release : updating READMEs and VERSION
       30e05e5 GenPipes - Nanopore : correcting minimap2 step
       7fb2d7c GenPipes - BFX : cleaning bash_cmp.py a bit
       f7538c5 GenPipes - Nanopore : fixing minimap2 and dependencies + allowing duplicate samples (if no duplicate readset) in the readset file
       0728713 GenPipes - Utils : renaming of  to
       cf7b0c6 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1da114b adding .gitignore to .gitignore...
       880fbf4 Version bump to 4.0.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       e503967 Merged in release_4.1 (pull request #302)
       8fe012c GenPipes - AmpiconSeq : updating resources for dada2 in cedar.ini
       8a0ea55 GenPipes - BFX : fixing "ln" in bash_cmd.py
       44a9dec Merged in fix_nanopore_eh (pull request #298)
       44ddf6b Merged in ehenrion/genpipes-chipseq-bashini-edited-with-u-1642626170885 (pull request #292)
       0f84b78 GenPipes - ChIP-Seq : bash.ini edited with updated versions of software, fixing Issue #127
       6c8219d Merged in release_4.0.0 (pull request #286)

  jgalvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       e0c9f5f Full squash of covseq_ont branch

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      5 commits

       2841e22 Final correction to the reference genome symlink
       323b5d2 Corrected issue with the reference genome link command
       35ad8b7 Added compatibility to ARTIC primer versions 4 and 4.1
       3f3fd3f Update version of ncov_tools to 1.8
       8bb76cd Added changes to the ini files to allow for ARTIC V4 schemes (and any potential future schemes).

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       fa5b71b Merged in covseq_v4 (pull request #287)
       ce25055 Added force option to reference symlink creation to avoid crash when re-running report.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      10 commits

       8915d42 Merged in covseq_nanopore (pull request #300)
       a242cd1 nanopore - Cit update - Updating cit ini
       6ab0488 nanopore - pycoqc - Fixing step to match with new bfx
       1f4b630 covseq_nanopore - Cit update - Updating cit ini
       803b08a covseq_nanopore - Cit update - Updating cit ini + removing module load in report
       34e981c Merged in covseq_nanopore (pull request #299)
       eaa85d8 covseq_nanopore - Cit - Code cleaning and cit.ini addition
       9786971 covseq_nanopore - Cit update - Updating cit ini
       9726784 covseq_nanopore - General - Renaming covseq and nanopore_covseq class for consistency
       fff0c8b covseq_nanopore - General - Fixing argparse type

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       6eddc14 Using cit magic for walltime
       c84a3b4 Merged in cpu_scheduler_agnostic (pull request #294)
       f73f008 Full squash of covseq_ont branch
       0670278 Merged in config_formater (pull request #289)
       f0047ac Merged in watchdog (pull request #290)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      30 commits

       2c5f270 Add mkdir in FFPE steps
       2860a3b cpu_str to node_str in scheduler class
       4205247 update cit.ini
       d67a799 update nanopore cit.ini
       2c679b9 make cluster_cpu scheduler agnostic
       5694384 x on tumor_pair.py
       13bd411 get symlink back in tumor pair
       b919a3e get bash_cmd back to old self
       9d722c2 fix regressions
       7d068c2 ad tmp to gitignore
       5c19819 update perl in tumor pair
       19aaf71 cluster_memmory -> cluster_mem
       88a70bc tweak tumor pair extra ini
       3fe6b0b missing int casting
       ce9a2e8 pmem for pbs/torque
       f3c406d cleanup sanitycheck log
       4357a19 SanitycheckError args and kargs
       fe6b300 add .message to SanitycheckError
       ddc1401 add cluster_mem to ini option
       e0aec1b force int on walltime read
       87f84cd force int on walltime read
       220b6e0 slurm to time_to_datetime
       853ee6f fix hicseq ini
       07f2b77 walltime method from slurm
       ad88e69 time delta for pbs walltime
       1227182 update utile time delta
       2f679b7 update utile time delta
       2fe4ab3 update utile for torque
       88abd41 add walltime format for pbs
       333c32a rename monitor.sh to watchdog

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       80e3e4d chipseq diffbind pca plot update
       dc0680a methylseq methylkit min cpg issue fix
       f87dde2 testing updates to DiffBind.R in mugqic_tools
       acaaa65 modified differential expression variable names to something meaningful
       da9bf66 fixed issue #129 added differential binding to the atacseq protocol

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       5919bcf Merged in chipseq_atac_diffbind (pull request #288)

4.0.0        Thu Dec 9 19:13:55 2021 +0000        232 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      14 commits

       a34c5a1 Merge remote-tracking branch 'origin/dev' into release_4.0.0
       d2f25eb updated .gitignore
       c02c96d Merge remote-tracking branch 'origin/dev' into release_4.0.0
       9485b93 GenPipes - README : update epiqc REAMDE
       18dc471 Merge remote-tracking branch 'origin/dev' into release_4.0
       8ac6847 GenPipes - README : Updated READMEs before new release
       b1e2df0 GenPipes - config : correcting parsing of walltime from config files
       b83fe9f Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0865f20 GenPipes - DNASeq : correction of sym_link_final_bam dependencies - fixes Issue #119 and Issue #125
       85f8054 GenPipes - Tumor Pair : fixed some enumerate loop...
       5199ee7 GenPipes - Tumor Pair : fixing sym_links steps to avoid duplicate job names
       0e63d7f GenPipes - DNASeq : fixing .bai dependencies
       3453cf0 Version bump to 3.6.3-beta
       d2361a5 Version bump to 3.6.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       963a49a GenPipes - python3 : fixes following py2to3 recommandations
       a227471 GenPipes - Python3 : fixes after last rebase with dev
       4b5f9bd Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       5ca19c5 GenPipes - start switching to python3
       bfa2a48 GenPipes - start switching to python3
       e3f7e92 GenPipes - Python3 : fixing file opening in bfx/readset.py
       7789689 GenPipes - Python3 : fixing conflicts with bfx/readset.py in dev
       7c41602 GenPipes - start switching to python3
       978a291 GenPipes - start switching to python3
       233c108 GenPipes - start switching to python3

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      4 commits

       2601b7f GenPipes - Python3 : removing the shebang from all the bfx scripts
       a99a2e4 GenPipes - Python3 : removing the shebang from all the bfx scripts
       c43eb20 Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       025c3c2 Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh

  Édouard Henrion <henrione@beluga5.int.ets1.calculquebec.ca>      2 commits

       4463bac GenPipes - Config : removed useless debugging messages
       d47adcf GenPipes - Resources : updated install scripts for fgbio and mugqic_tools

  ehenrion <edouard.henrion@mcgill.ca>      8 commits

       baf52f9 Merged in release_4.0.0 (pull request #285)
       540e030 Merged in release_4.0 (pull request #284)
       fe49a83 GenPipes - BFX : corrected typo in gatk4.py
       bdea57f Merged in tumor_pair_sym_link_fix (pull request #282)
       f594754 Merged in hotfix_dnaseq (pull request #280)
       34e54c9 GenPipes - RNASeq : minor typo fixes in README.md
       eca4475 Merged in genpipes_python3_eh (pull request #270)
       34608bc Merged in release_3.6.2 (pull request #269)

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       2b23422 Merged in dnaseq_fixRecal (pull request #276)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       b28b4a8 remove recalibrated bam name error from some input jobs

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      8 commits

       73443c7 Making pandoc report working with pandoc version 2.16.1
       60093e6 Making pandoc report working with pandoc version 2.16.1
       c6ab31b Merged in Paul-Stretenowich/log_reportpy-edited-online-with-bitbucke-1637772592796 (pull request #277)
       e869259 log_report.py switching back to py3 and adding Narval as remote
       0ef01e9 EpiQC - First commit after Rami Coles internship
       4e55e0f General - Fixing genome installation grep issue
       06260f0 EpiQC - First commit after Rami Coles internship
       d458b42 Merge branch 'dev' into epiqc

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      9 commits

       4f89f91 copy contents of readset into inputinfo file(have errors)
       d51b1ce copy contents of readset into inputinfo file(have errors)
       97bc3da creating inputinfo file after copying original file
       45c917c epiqc - completed chromimpute development - right before remove design file and modify readset file
       9500856 epiqc - completed chromimpute development - right before remove design file and modify readset file
       2bb45e0 [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       3a86f94 copy contents of readset into inputinfo file(have errors)
       acfcbc0 copy contents of readset into inputinfo file(have errors)
       7b5b552 creating inputinfo file after copying original file

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       43f4575 update GIAContainer version
       9fec333 changing error log to debug in cit config
       9329e2e add cit options for epiqc
       10199cc send telemetry downloded file to /dev/null
       e942d07 fix warning
       99ec4c9 make log_report.py more robust
       c30189c update R bioconductor in dnaseq_high_coverage

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      47 commits

       4bee1ca bug fix #121. missing module load mugqic_tools 2.8.2
       5e8f3da fixed file path issue in seq2fun corrected the mugqic_tools version in epiQC
       a317b6b corrected a typo
       7f78005 added the reason to use os.remove in a comment
       c4721bc addressed Hector's comments - changed seq2fun db path
       ffd0abe addressed Hector's comments
       6e2ce1f fixed Issue #118: hicseq hic error in interaction_matrices_Chr (mugqic/genpipes) by adding floor division to places which return floating points where intgers needed
       c333a3e modified incorrect seq2fun cvmfs db path
       5e928cf added seq2fun cvmfs db paths
       0c9418a modified seq2fun db paths
       1569435 updated readme file
       22d90bf modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       dc4c608 seq2fun pathway_analysis completed
       de5f87b started pathway analysis-seq2fun
       e162747 improved seq2fun processing
       5d6df18 added seq2fun to the function as the final step. working well but needs improvements
       bdd622d fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       54bf3aa added seq2fun protocol completed merge_fastq partially completed seq2fun function
       0c2d569 addressed Paul's comments
       ef8ae36 fixed a bug
       2e3c742 changed mugqic_tools version to 2.8.3
       06e458e fixed bugs when transferring to python3
       fb9093e fixed bug in chromimpute convert corrected uuid issue in job2jason.py
       0bd6773 fixed an issue for python 3
       0d1622c updated readme file
       9e62546 testing whether env varibales are retrievable dynamically
       ff73039 fixed input files not found error when running the pipeline outside of the project directory
       78af37b added cit.ini
       a4dc68c fixed issues in final report added comments
       ca3dc8a fixed bugs in -o option added IHEC data path
       0bcde68 modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       affd3b0 seq2fun pathway_analysis completed
       129b56b started pathway analysis-seq2fun
       7ef4a6f improved seq2fun processing
       92ac712 added seq2fun to the function as the final step. working well but needs improvements
       97bdedb fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       2d0359b added seq2fun protocol completed merge_fastq partially completed seq2fun function
       93811b3 testing whether env varibales are retrievable dynamically
       72d5f9b fixed input files not found error when running the pipeline outside of the project directory
       b97891f added cit.ini
       b29e717 fixed issues in final report added comments
       cea9f45 fixed bugs in -o option added IHEC data path
       15c5695 completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       b86e6e5 completed fixing errors-working pipeline
       d07286d corrected up to epiqc final report
       9304cba corrected up to epigeek
       562607a corrected up to global dist after chipseq design changes

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      59 commits

       e837166 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       5683bc9 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       faab0ee Addressed Paul's comments
       b3e8185 addressed Paul's comments
       4130054 Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       0225132 space changes in epiqc.py
       ca91e09 completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       dbf8254 completed fixing errors-working pipeline
       fc17d37 General - Fixing genome installation grep issue
       26837e3 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       e3ca2e6 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       9bda353 EpiQC - Fixed bugs
       e4faee4 EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       1dcd8dd EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       851b26d added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       29139a9 General - Fixing genome installation grep issue
       b95fc32 resolve merge conflicts
       701f109 resolved merge conflicts
       d03737e resolve merge conflicts
       8852d61 resolved merge conflicts
       5463486 resolved merge conflicts
       82ff100 resolved merge conflicts
       d833412 resolve merge coonflicts
       e047a3e resolved merge conflicts
       dba6b33 space changes in epiqc.py
       2fe7046 Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       429caf1 General - Fixing genome installation grep issue
       06113d1 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       34d3cc9 EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       c00875a EpiQC - Metrics thresholds can be modified in base.ini
       35a6396 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       9e3e1f4 \Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       698c46f EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       2086d5a EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c2af489 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       1089984 EpiQC - Fixed bugs
       bf60df9 EpiQC - Added signal to noise and epigeec steps
       bced2bd EpiQC - Parallelized chromimpute
       c1b462a EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       4d12401 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       0e763c8 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7f8337e [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       2b769bf resolve merge conflicts
       05b7cd3 resolved merge conflicts
       2e219fa resolve merge conflicts
       aaed941 resolved merge conflicts
       cc0bccf resolved merge conflicts
       231f626 resolved merge conflicts
       a7e5a61 resolve merge coonflicts
       42a741b resolved merge conflicts
       1e5e8af Merge branch 'epiqc' of https://bitbucket.org/mugqic/genpipes into epiqc
       6eef5c8 resolve merge conflicts
       d679a97 resolved merge conflicts
       eb880f5 resolve merge conflicts
       06fe31e resolved merge conflicts
       2296539 resolved merge conflicts
       0e2c9cd resolved merge conflicts
       301d2ee resolve merge coonflicts
       f285a7c resolved merge conflicts

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       b3194fa Merged in seq2fun_fix_bug121 (pull request #281)
       ffc9c01 Merged in epiqc_mugqic_tools (pull request #279)
       a9ac4d6 Merged in epiqc_comments (pull request #278)
       4bd2b26 Merged in epiqc (pull request #271)
       c38208b Merged in seq2fun_denovo_rnaseq (pull request #272)
       e12b111 Merged dev into seq2fun_denovo_rnaseq
       0df8920 Merged in hic_python3_issue_fix (pull request #275)
       1d5f2b0 Merged epiqc into epiqc_ss

  rami.coles@mail.mcgill.ca <rcoles@abacus1.ferrier.genome.mcgill.ca>      2 commits

       70f8ec0 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       5feada1 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec

  rami.coles@mail.mcgill.ca <rcoles@abacus2.ferrier.genome.mcgill.ca>      11 commits

       f20a000 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       74e2a32 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       e911f95 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7a5da00 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       bc66beb Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       d0a7b9d added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6695197 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       a040d0b added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6a5e995 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       53c747f Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       f56d086 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py

  rami.coles@mail.mcgill.ca <rcoles@abacus3.ferrier.genome.mcgill.ca>      21 commits

       6901237 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       88c422d EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       f85d37a EpiQC - Metrics thresholds can be modified in base.ini
       b146548 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       51b6bdb Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       9e118e5 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       a68aad3 EpiQC - Added signal to noise and epigeec steps
       fe957a6 EpiQC - Parallelized chromimpute
       ab1425d Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       1cbf6d0 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       86c0461 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       ef3acde EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       6425046 EpiQC - Metrics thresholds can be modified in base.ini
       e9e9f2f EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       3ecaebb Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       2a773e7 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c800a67 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       d3bd3e1 EpiQC - Fixed bugs
       402b34e EpiQC - Added signal to noise and epigeec steps
       a664314 EpiQC - Parallelized chromimpute
       4bf48b9 EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.

  Robert Syme <rob.syme@gmail.com>      1 commits

       f8dfff5 Merged in servername-fix-dev (pull request #274)

  Rob Syme <rob.syme@gmail.com>      1 commits

       a76c958 Fix example server name in dev

  Shaloo Shalini <shaloo.shalini@gmail.com>      4 commits

       db1afc3 Merged in ss_mermaid_96 (pull request #247)
       dab640b Merged dev into ss_mermaid_96
       6193cda Merged in dev_covseqdoc_ss (pull request #258)
       6dcc0cb Merged in epiqc_ss (pull request #253)

  shaloo <shalz@hotmail.com>      14 commits

       9e54762 Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       4ec66ff Refs #102 removed ChIP-seq pptx link
       7837688 Fixes #102 issues in README.md for epiqc pipeline
       e08c73d Refs #103 Paul's feedback in workflow incorporated
       52ca4d1 Fixes #103 covseq workflow added
       f356ffe Merge remote-tracking branch 'refs/remotes/origin/dev_covseqdoc_ss' into dev_covseqdoc_ss
       878fcf5 Refs #103 covseq doc freebayes update
       8f5625c Typo fix
       81ee462 Fixes typo Refs #103
       a62e20b Refs #103 covseq doc freebayes update
       fb32de5 Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       86e5a38 Refs #102 removed ChIP-seq pptx link
       29ae522 Fixes #102 issues in README.md for epiqc pipeline
       43501a4 Fixes #96 all GenPipes workflows are now coded as mermaid flowcharts and png generated from them

3.6.2        Thu Nov 11 14:14:18 2021 +0000        28 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       3787dea Cleaning before release
       c1aef1d GenPipes - Readset : improved readset parser so that it creates a temptative readset file with unique IDs when the readset file provided has duplicate readert IDs
       cead55d GenPipes - MethylSeq : correted bed2interval_list calls + some minor code reformating
       499a0f6 Version bump to 3.6.1

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      3 commits

       55074b1 GenPipes - README : updating READMEs and install scripts before release
       448d7c8 GenPipes - updating version of mugqic_tools in all the pipelines config files
       1c2d322 GenPipes - Resources : adding some new genome and software install scripts

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       eb01979 Merged in release_3.6.2 (pull request #267)
       2fb1dc1 GenPipes - Readset : correcting the parsing of readset files to stop allowing duplicates headsets in file - fixes Issue #113
       72170d9 GenPipes - RnaDeq denovo Assembly.py : updated merge_trimmomatic_stats outputs to fix insilico_read_normalization_all_report dependencies
       50519a2 GenPipes - RnaDeq denovo Sssembly.py : fixed insilico_read_normalization_all_report dependencies
       2e429b1 Merged in release_3.6.1 (pull request #262)

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       2eca450 Merged in Pierre-Olivier-Quirion/ampliconseqbaseini-edited-online-with-bi-1635366894762 (pull request #266)
       3d108f9 ampliconseq.base.ini edited online with Bitbucket

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       cc3f279 update all pipelines to pandoc 2.16.1
       cd9c041 new pandoc

  Robert Syme <rob.syme@gmail.com>      3 commits

       4ab8e17 Merged in wget_warn_on_fail (pull request #265)
       d20df3e Merged dev into wget_warn_on_fail
       b809a05 Merged in rnaseq_protocol_switch_fix (pull request #263)

  Rob Syme <rob.syme@gmail.com>      9 commits

       20154ee Fix server name
       31fd5fe Test with incorrect server name
       17aea66 switch quote style
       8c4410e Intentially introduce an error in the reporting server.
       beb177e Warn when wget command fails.
       ca54581 Another tiny whitespace fix
       cc33046 Tiny whitespace fix
       b3d620c Switch at correct location
       2cbdcee Fix protocol mixup for rnaseq pipeline

3.6.1        Wed Sep 29 20:53:55 2021 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      6 commits

       27ccb3b correct log_report.py before release
       39e1b94 Merge branch 'master' of bitbucket.org:mugqic/genpipes into release_3.6.1
       0fae24c Updating GenPipes README before release
       df8f9f4 Merge branch 'dev' into release_3.6.1
       96ce81f Re-creating READMEs before release
       530576e Version bump to 3.6.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       1cda65c Merged in release_3.6.1 (pull request #261)
       b472e72 Merged issue_105_fix_eh into dev
       7cd2f09 GenPipes - Tumo Pair: correcting job name in sym_link_fastq_pair
       4df1042 GenPipes - Tumor Pair.py : fixing typo in sym_link_fastq_pair
       d4613be GenPipes - Tumor Pair : fixing issue #106, no more job name overwriting at gym_link_fastq_pair step
       3079e60 GenPipes - DNASeq : fixed bed2interval_list call in gatk_haplotype_caller step
       03c03bf Merged in release_3.6 (pull request #259)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       7947042 Merged in Paul-Stretenowich/sambambapy-edited-online-with-bitbucket-1630531961392 (pull request #260)
       69e1a12 Fixing sambamba sort output files: there were a bai set as output that was blocking the "resume" for covseq pipeline.

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       a026207 more robust log_report
       b941202 remove echo debbug in get wrapper

3.6.0        Mon Aug 30 17:55:36 2021 +0000        370 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       cdaafcd Merge branch 'dev' into release_3.6
       bea48bd GenPipes - updating CHANGELOG and VERSION files (after release)
       2ea1988 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       4554eaf GenPipes - README : updating version to 3.5.0 in all pipeline READMEs
       1f71288 GenPipes - Release 3.5.0 : updating CHANGELOG up to tag 3.5.0

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      16 commits

       9eda417 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       79bc0fd GenPipes - Tumor Pair : no more attemps to write where the reference files are (in case of readonly FS, e.g. CVMFS...)
       79aa7dc Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       86c40fd GenPipes - Tumor Pair : fixed Strelka2 jobs with mkdir
       ed205e1 GenPipes - Tumor Pair : fixed coverage_bed in strelka steps
       b7873b9 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0eafe86 GenPipes - BFX : fixed GATK4 cluster_crosscheck_metrics command
       a9eb687 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       ae4cc08 GenPipes - Tumor Pair ; deleted dev.ini from dev branch
       b6fd739 GenPipes - Tumor Pair : fixed symlink_fastq_pair command
       41ba984 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       84f350f GenPipes - Install script : version change
       7e07739 GenPipes - RNASeq light : fixing typo
       e1c1b50 RNASeq light : importing bash_cmd fix
       3e29063 GenPipes - Release 3.6 : updating the READMEs in dev
       16be170 GenPipes - Pipelines : Fixed dnaseq_high_coverage inheritance with dnaseq and rnaseq_light & rnaseq_denovo_assembly inheritance with rnaseq

  ehenrion <edouard.henrion@mcgill.ca>      13 commits

       677e9cf Merged in release_3.6 (pull request #250)
       3591ec9 GenPipes - ChipSeq : updated mugqic_tools with latest version
       da2b410 GenPipes - Ampliconseq: updated mugqic_tools version in base.ini
       40cb3b4 GenPipes - AmpliconSeq : updated mugqic_tools version in base.ini
       92e4919 Rnaseq_light.py: kallisto_count_matrix dependency
       18d37ad fixing picard2 add_read_groups call
       952ac0a GenPipes - MethylSeq.py : Fixed pipeline inheritance with DNASeq
       de571d6 GenPipes - CoVSeq : setting the gatk4_java_options in base.ini
       11755c1 GenPipes - BFX : fixing add_read_groups in picard2.py
       fa40033 GenPipes - MethlySeq : fixing sambamba_merge_sam_files dependency
       2f63183 Merged in pipeline_inherit_fix_eh (pull request #249)
       aa3ca79 Merged dev into pipeline_inherit_fix_eh
       273fe54 Merged in release_3.5 (pull request #243)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      24 commits

       a1f3040 Merged in covseq_stretenp (pull request #248)
       8374fdb Threshold rename consensus change for LSPQ
       d3c3b5f Code cleaning
       7bc5189 Merged dev into covseq_stretenp
       376368f Upgrading covseq_tools version
       e473abe Changing regex for finding negative controls for LSPQ
       d0291fe covseq - prepare_table - debug
       0ed94bb Merged dev into covseq_stretenp
       b768d92 Adding bfx for covseq_tools and ncov-tools
       0c8407e Upgrading cvoseq_tools version
       6b3f8de Debug
       7dcbf5d Degrouping rename_consensus step for ivar and freebayes
       d476c11 Debug
       00fafc3 Un-grouping freebayes report and ivar report
       1fe319e Un-grouping freebayes report and ivar report
       6b32016 Debug
       3023357 Changing resources for ncovtools
       58d0497 Debug
       29d800f Debug
       8a773f9 Debug
       5ede73b Debug
       ce60e36 Debug
       84cda2b Adding freebayes metrics and reporting
       d0592cc Renaming prefix parameter in ini_section

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       eeee6b1 Merged in tp_debug (pull request #257)
       fae1349 Merged in monitor_limit (pull request #256)
       6283fd9 Merged in kallisto_rnaseq_light (pull request #254)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      13 commits

       8d09333 Bump GiaC to v2.0.3
       55dd941 patch stolen from Ed to fix recalibration step
       a89de00 RETRY ret_code corrected
       64b9f79 RETRY limit on monitor
       c638ab3 RETRY limit on monitor
       bf87297 RETRY limit on monitor
       2c1c895 RETRY limit on monitor
       207a069 fix ini for graham covid
       203ea90 update tp ini for cedar
       64db31b reducing constraint on fit for integration testing
       a17345d minor fixes before release
       248078c fix kallisto dependency
       5583f3f giac updtated to v2.0.2

  P-O Quirion <pioliqui@gmail.com>      2 commits

       fb635da more robust log repport
       55b928d make log report more robust went jobs are failing

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       e56c39d corrected a typo
       9110f86 modified readme files for chipseq differential binding step

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       fbdbe38 Merged in chipseq_readme_fix (pull request #245)

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      40 commits

       ac5a5c5 Homo_sapiens.GRCh38.ini fixes
       675d5dd samtools mpileup spacing fix
       f182779 cit fixes after dev.ini removal
       15e9882 multiqc ini fix
       1cdfe28 samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       5977acd add portal dir to dnaseq
       3f93719 hs37d5 and cit strelka2 fixes
       5a20aa9 removal of cit_gatk4.ini, contents moved to cit.ini
       59473ce tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       d98dd0c merge vardict fix
       a64d34c job.sample fixes for jsonator
       63b37d8 vardict exome fixes
       13b7d0f germline variant annotator fixes
       d366e78 strelka2 germline and sequenza fixes
       19b187c cit adjustments to tumor_pair
       b999375 beluga.ini fixes
       46f7416 strelka2 germline fixes
       3d5693d removal of samtools, substitute samtools germline for strelka2
       ec213ce conflict fixes to dnaseq and tumor pair after dev merge
       8cad167 fix mpileup dependency chain
       e6d6aa5 corrections to mpileup bcftools merge and tumor_pair dependencies
       1dc88b9 samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       9e3ac20 add portal dir to dnaseq
       4a59db3 portal_dir fix
       4dfa21b added portal_dir for cit testing
       528a77b hs37d5 and cit strelka2 fixes
       9c66ec2 removal of cit_gatk4.ini, contents moved to cit.ini
       03a602b tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       30d7b07 merge vardict fix
       83b75bf job.sample fixes for jsonator
       ec14547 vardict exome fixes
       e5a6145 germline variant annotator fixes
       fdd20bf strelka2 germline and sequenza fixes
       a4f8581 cit adjustments to tumor_pair
       7240e3e beluga.ini fixes
       849625a strelka2 germline fixes
       80f5548 removal of samtools, substitute samtools germline for strelka2
       e92837b conflict fixes to dnaseq and tumor pair after dev merge
       9aedd69 fix mpileup dependency chain
       e55e1ff corrections to mpileup bcftools merge and tumor_pair dependencies

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      45 commits

       6ec0367 collect metrics fix - dnaseq
       fd8ee63 PR fixes to sambamba and multqc
       afc4b8f PR conflict fixes
       04c160e updated dev.ini
       53cde24 cit fixes
       359c1ae purple sanity check fixes
       6d12435 cit fixes to purple and vardict exome
       5da7760 fix for multiple tumors with the same control
       d7724ad deliverable fixes - sym link final bams
       63c8bc8 non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       39a2651 bcftools update and filter argument fixes
       6741048 strelka2 germline and ensemble fixes
       b1f7396 sequenza and cit fixes
       c9cef16 fixes to gatk4 for tumor pair
       bc784a6 gatk4 vsqr cit fix and baf plot for b38
       e148bd4 argument fixes for picard imported functions in gatk4 and vqsr fixes
       036809d fixes to multiqc
       7b0712b exome specific fixes
       059c613 updates and fixes for cit
       dbe0ff4 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       3c52e66 remove briaree ini and update dnaseq base
       9343e04 updating beluga ini
       d79e994 updates to beluga.ini and base.ini for dnaseq
       97793f0 ini updates
       f716b0e updated dev.ini
       0f7cccc cit fixes
       45c9d0b purple sanity check fixes
       882096d cit fixes to purple and vardict exome
       97af93a fix for multiple tumors with the same control
       49f4f38 deliverable fixes - sym link final bams
       17f001e non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       6729217 bcftools update and filter argument fixes
       afba10f strelka2 germline and ensemble fixes
       6326449 sequenza and cit fixes
       e8a4ac9 fixes to gatk4 for tumor pair
       3348617 gatk4 vsqr cit fix and baf plot for b38
       b8ec9bd argument fixes for picard imported functions in gatk4 and vqsr fixes
       9b714d0 fixes to multiqc
       ca97ea7 exome specific fixes
       2e20be2 updates and fixes for cit
       6793348 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       8b46abe remove briaree ini and update dnaseq base
       4795a96 updating beluga ini
       04d851f updates to beluga.ini and base.ini for dnaseq
       ab85a84 ini updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      37 commits

       6cec0ef module fix
       7b697cc removal of dev files for CVMFS version
       46ef162 genome ini fixes
       466eb23 dev and hs37d5 fixes
       367666b tumor pair - sambamba markdup fixes
       d5b857f module fix
       e7a2849 sambamba markdup fix
       f7e4fdb strelka cit fixes
       c320d5a added strelka2 input dependency for purple purity
       d3c3292 cit fixes to purple and vardict exome
       89c80a9 muliple normal fix for fastqc
       00531ff dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       d1840a0 fixes and strelka conversion addition to purple
       09d9135 added purity estimate program PURPLE
       aaaff33 strelka bed manipulation fix
       f257d83 vardict exome fix
       f74092d conpair and collectHS metric fixes
       2822b5f tumor_pair qualimap part 2
       51e88ab sym link dnaseq.base into tumor pair
       03e567c updates to b38 variant recal files
       4757879 fixes to tumor_pair on beluga
       52e48ce cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       febb823 strelka cit fixes
       42fd97c added strelka2 input dependency for purple purity
       54fd0cd cit fixes to purple and vardict exome
       858af14 muliple normal fix for fastqc
       bd1f034 dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       1b2c03f fixes and strelka conversion addition to purple
       0957242 added purity estimate program PURPLE
       37529e6 strelka bed manipulation fix
       cc31781 vardict exome fix
       2091afb conpair and collectHS metric fixes
       8689d63 tumor_pair qualimap part 2
       a9888f7 sym link dnaseq.base into tumor pair
       e340eb0 updates to b38 variant recal files
       7a7b46d fixes to tumor_pair on beluga
       a449ed5 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      28 commits

       58152ae varscan version downgrade
       7fd72e1 cit fixes
       3a63b4a fixes to bcftools mpileup for dnaseq mpileup protocol
       b92ebd8 hs37d5/GRCh37 fixes for purple
       51c9658 fastpass varscan merge fix
       2cc7f00 cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       93ed2f3 cit fix for conpair
       9672af2 sym link fixes
       c4d0581 cit fixes to tumor pair
       41ccead fixes to dnaseq
       b07afef updates to GRCh38 annotation file, module updates, cit fixes
       066fb49 fixes to deliverable and b38 ini
       034447e updates to cit and fixes to one job mpileup steps
       17d4c79 updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       4668a0b major fixes to deliverables and completion of beluga test
       9ffad22 fixes to bcftools mpileup for dnaseq mpileup protocol
       180ae35 hs37d5/GRCh37 fixes for purple
       87dcc11 fastpass varscan merge fix
       f662715 cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       233bc03 cit fix for conpair
       8084b8f sym link fixes
       c5840c6 cit fixes to tumor pair
       e4193ca fixes to dnaseq
       14aa4a3 updates to GRCh38 annotation file, module updates, cit fixes
       de7ed14 fixes to deliverable and b38 ini
       9510061 updates to cit and fixes to one job mpileup steps
       510efd1 updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       cad6781 major fixes to deliverables and completion of beluga test

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      43 commits

       05d65f7 dnaseq sym_link fix
       c49e469 dnaseq ini fixes
       6b63aea dnaseq beluga and cedar ini fixes
       0400745 sequenza bug fix
       362539f revert back to old sequenza-utils
       3714cd0 fixes to sequenza for exomes, and minor SV fixes
       444e6c2 exome cit fix to sequenza
       9b97f7e converted -n to -c
       bf70c51 cit fixes - runtimes
       bfd86d9 config file fixes
       75fd219 further sanity-check fixes
       30c2012 sanity check fixes
       adb0724 multiple pair same control fixes and vardict exome fixes
       f0c3f89 cit fixes to vardict
       34b052d cit fixes for tumor exome
       21e17a2 better cram compression of base recalibrated bams with variantBam
       baf4196 sequenza and germline ensemble fixes
       b736a4c sambamba mark_dup fix
       667a4fe module fixes
       cee1774 sym link dnaseq cit to tumor_pair cit
       1e7e2e3 tumor pair cit updates
       8d35d7f dnaseq sym_link fix
       7da9cee dnaseq ini fixes
       3d91118 dnaseq beluga and cedar ini fixes
       e1a83a3 sequenza bug fix
       f494d73 revert back to old sequenza-utils
       ee70ffc fixes to sequenza for exomes, and minor SV fixes
       9ffeafa exome cit fix to sequenza
       cb1a59c converted -n to -c
       00943b3 cit fixes - runtimes
       29c3394 config file fixes
       d3eedbb tumor base fix
       bed7a1d further sanity-check fixes
       d36c1e8 sanity check fixes
       e901171 multiple pair same control fixes and vardict exome fixes
       165e341 cit fixes to vardict
       10aa2b1 cit fixes for tumor exome
       ceef654 better cram compression of base recalibrated bams with variantBam
       3d476d2 sequenza and germline ensemble fixes
       1edc6d7 sambamba mark_dup fix
       a3b1208 module fixes
       caabe63 sym link dnaseq cit to tumor_pair cit
       b9ed61e tumor pair cit updates

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      8 commits

       e7bd745 code cleaning and fixes to exome interval list
       a7d6407 fixes to cedar ini
       cca1a4c fixes to symlinks for paired indel realignment
       211a29e cedar ini and exome update
       f65b070 code cleaning and fixes to exome interval list
       febded3 fixes to cedar ini
       43d4044 fixes to symlinks for paired indel realignment
       09c93a5 cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      4 commits

       c7f3fa9 cedar fixes and GRCh38 fixes
       2baaf56 cedar germline sv updates
       7367777 cedar fixes and GRCh38 fixes
       bd82ef1 cedar germline sv updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      46 commits

       21b8896 fixes to sambamba for lumpy - query sorting instead of corrdinate
       dc9750d updates to metasv - somatic
       c19ad15 Fixes and updates to reference files
       08f0027 remove testing steps
       aad782e updates to SV germline and reference genome tweaks
       491914d somatic sv fixes: lumpy and svaba annotations
       d71a64d json fix
       03bcbc6 fix to germline SV: breakseq2
       97a4bc1 deleting hidden files
       b24f8f3 merge fixes
       9052826 GATK4 fixes - bam indexing and markDupSpark
       fc33ed9 bcftools fixes for tumor pair
       0deb8a9 fingerprint and bug fixes
       57f218a dnaseq - vcftools qc addition: --missing_indv and --depth
       5c82925 GRCh38 fixes
       bb43292 select input for variant caller and fixes to one job calling
       b87e9f5 json and folder updates
       f15fd74 fixes to sCNAphase
       34c9e47 Bug fixes prior to json additions
       ae851f4 merging snv and sv, adding protocols
       c82ce68 Updates and debug
       5788360 Add set somatic and actionable mutations
       bc0ca0a added multiqc and other tweaks
       c51db53 fixes to sambamba for lumpy - query sorting instead of corrdinate
       8da9446 updates to metasv - somatic
       be214eb Fixes and updates to reference files
       f717bb7 remove testing steps
       4c64223 updates to SV germline and reference genome tweaks
       541f671 somatic sv fixes: lumpy and svaba annotations
       c795a24 json fix
       8cc8e79 fix to germline SV: breakseq2
       628ce3d deleting hidden files
       2a83ec4 merge fixes
       90ce1f4 GATK4 fixes - bam indexing and markDupSpark
       4e8c750 bcftools fixes for tumor pair
       4c874cc fingerprint and bug fixes
       3cec8c2 dnaseq - vcftools qc addition: --missing_indv and --depth
       6f3ae37 GRCh38 fixes
       ebfc764 select input for variant caller and fixes to one job calling
       4e0b29b json and folder updates
       fc71742 fixes to sCNAphase
       a8034e2 Bug fixes prior to json additions
       2b23267 merging snv and sv, adding protocols
       3c8b861 Updates and debug
       92db229 Add set somatic and actionable mutations
       39425ef added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      24 commits

       3e96d81 dnaseq - fixes for VariantBam
       2cd2eb5 fixes to purple and mkdir vardict fix
       3ba332b sequenza WGS fixes
       86f1fa3 fixes to metasv for tumor pair
       3d446d9 fixes to sv germline calls
       6f26955 Single job bug fixes
       893626b manta I/O fix and other bug fixes
       9f4109b config updates and b38DH added
       9098c40 dnaseq germline SV updates
       8d737df gatk4 updates and bug fixes
       843ef92 Json related bug fixes
       54fdadc Bug fixes and modification derived from initial PROFYLE benchmarking
       5cd6733 dnaseq - fixes for VariantBam
       047e961 fixes to purple and mkdir vardict fix
       776c53a sequenza WGS fixes
       dc8f820 fixes to metasv for tumor pair
       81b508b fixes to sv germline calls
       46528cd Single job bug fixes
       be5f8d2 manta I/O fix and other bug fixes
       fb73873 config updates and b38DH added
       9c2a8b0 dnaseq germline SV updates
       e86e74a gatk4 updates and bug fixes
       e9b3a2f Json related bug fixes
       0201d38 Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       e9cc405 seqeuenza fixes and sv fixes
       3786a03 gatk4 mutect2 updates
       4ae76b7 cit-based fixes to NGScheckmate
       58c3db3 dnaseq qc additions: NGScheckmate and peddy
       9ae4210 seqeuenza fixes and sv fixes
       dc79320 gatk4 mutect2 updates
       0f9e09a cit-based fixes to NGScheckmate
       45837b4 dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      6 commits

       f00e0b3 Merged in dev_reveleig (pull request #252)
       00fb9e9 Merged in rebasing_tp (pull request #246)
       08361bc Debugging and Guillimin specfic fixes
       f469997 updates to config
       07d625a Debugging and Guillimin specfic fixes
       87f1ef4 updates to config

  Shaloo Shalini <shaloo.shalini@gmail.com>      1 commits

       78f59d5 Merged in ss_wf_chipseq_93 (pull request #244)

  shaloo <shalz@hotmail.com>      1 commits

       8c742ce Updates workflow for chipseq -t chipseq case with differential binding step and dependency update Fixes #93

3.5.0        Mon Jul 12 16:41:20 2021 +0000        1071 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      8 commits

       fd8c1cd GenPipes - Release 3.5.0 : correcting typo in chipseq README
       038dcb2 GenPipes - Release 3.5.0 : correcting typos in READMES
       6e4f4de Version bump to 3.5.0
       48f97ec GenPipes - Release 3.5.0 : correcting typos in pipeline README files
       f892b2e GenPipes - Release 3.5.0 : updating VERSION and README-release files
       82b28e3 GenPipes - Release 3.5.0 : corrected typo in README
       0cf907c GenPipes - Release 3.5.0 : updating version in the pipeline README files
       ddb421a Merge branch 'dev' into release_3.5

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      6 commits

       08bd50a GenPipes - Resources : updated picard index command in genome installation script to run command on a compute node instead of on the login node
       02124fc correcting typo genome install script
       ae646a7 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6981819 GenPipes - Resources : updates in software and reference install scripts
       0659d35 GenPipes - Resources : adding software installation scripts
       53dcb83 GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      47 commits

       339172c Merged in release_3.5 (pull request #242)
       525bfe1 Merged in release_3.5 (pull request #241)
       0d0f7b4 Merged in bash_cmd_for_covseq_update_eh (pull request #236)
       5a54cf9 GenPipes - CoVSeq pipeline : updated call to "zcat" in covseq.py
       2a056bd GenPipes - BFX : updated bash_cmd.py
       2809e5e GenPipes - Install script : added Nanopore pipeline to PATH, fixing Issue #74
       daebeb7 Merged in ehenrion/changelogmd-edited-online-with-bitbucket-1620161608327 (pull request #220)
       210d516 Version bump to 3.4.0
       a36b2cb GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       7835760 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8e2018d GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       b4cc95e GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       f4534a4 GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       e2193e6 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       4b416a0 GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b34a1a5 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       05a90a6 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       5748a6b GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       eb4a375 GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       c1a8b6b VERSION edited online with Bitbucket
       4a0537a GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       2fd4f5e GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       40125ca GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       30917fe GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       810feb8 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       8911666 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       3d76f2c GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       b2b9c45 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       929424c GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       879f109 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       98c7730 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ad82c94 GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       31968ee GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       8dff24a GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       dbbaa11 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       4b617ba GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       04f7714 GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       f1a49af GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       f9c2aef GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       72cd3e9 GernPipes - DNASeq : corrected merge_filter_bcf outputs
       3b635b2 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       04b24d5 GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       55ab017 GenPipes - DNASeq SV : fixing delly input error
       09c1bad GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       a756e77 dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       b514d9f GenPipes - Bug fix : corrected dnaseq.cedar.ini
       26b421e GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       3297e97 Fix Issue #86
       ff3110d Quick correction to address misleading headcrop parameter in the RNAseq base ini.
       e3d7b04 rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      344 commits

       ee6fd30 Merged in covseq_stretenp (pull request #239)
       76bfba9 Fixing sambamba version
       4bc5b99 hic code improve
       2f1529e Debug pdf rendering
       b76a994 Fixing ini
       7a6e9c1 Merge branch 'dev' into covseq_stretenp
       985a5c5 Increasing resources for ncovtools and adding run_name for cit
       83caadf Upgrading bcftools
       cdbb414 Fixing merge_hicrep_scores
       b7ecfc4 Switching to samtools 1.12
       8b29ef3 Fixing merge_hicrep_scores
       beab866 Samtools version update in all pipelines
       00509e8 Merged in covseq_stretenp (pull request #237)
       74a1203 Merged dev into covseq_stretenp
       63cb7a2 Fixing indentation error
       22c619a Merged in covseq_stretenp (pull request #235)
       16854eb Switching to released version of covseqtools
       a69ddeb Cleaning comments
       4d85dff Removing guillimin and mammouth inis
       f25cf78 Merged dev into covseq_stretenp
       e955aa2 Merged in chipseq_design_change (pull request #234)
       5ca220d Merged dev into chipseq_design_change
       4c4678b Removing name in report/metrics job name
       0b1c562 Debug
       8c65979 Debug
       13f9971 Debug
       9cb7669 Merged dev into chipseq_design_change
       bb6e421 Merged dev into covseq_stretenp
       a1650eb Forcing ncovtools execution
       5f908a7 Fixing metadaya ncovtools file creation
       debf790 Adding filtered bam to ncovtools
       13b80b0 Fixing ncovtools
       296dfea Fixing peak file naming
       0133eb0 Forcing snakemake execution
       9e3647e Forcing snakemake execution
       0b32f4b Forcing snakemake execution
       c85a9fe Debug
       a562dd2 Changing resources for ncovtools
       59fdae5 Debug
       2e4dcbe Adding dependencies for reporting job
       558cebc Purging modules
       2147b85 Adding more resources for ncovtools
       2366c90 Debug
       355958d Debug
       4c6955b Debug
       dff338f Debug
       acca409 Debug
       e18a0b7 Debug
       c97470b Debug
       7f82b37 Debug
       4ee9190 Debug
       a50112b Debug
       23233fb Debug
       e143ff1 Debug
       9c6db14 Debug
       98c0231 Debug
       4e8b8e6 Debug
       290e9ab Debug
       032d425 Debug
       72cc786 Adding GenPipes version accseible inside GenPipes
       49210d2 Degrouping ncvotools from report generation because of cd for ncvotools
       f884631 Adding modules reading for reporting
       a3209fd Debug
       4717b02 Debug
       bb6b07c Debug
       8587d12 Debug
       612afbf Debug
       f7e17ec Debug
       b9de420 Debug
       2fde5e2 Debug
       7a124f1 Debug
       feb2364 Debug
       d15188b Debug regex prepare ncovtools
       de76075 Debug
       4b83dab Debug
       8b1e591 Reporting init ncovtools yaml
       b226f89 Debug finding files for report
       bc78d4a Debug test dependencies reporting
       d519ac7 Adding missing sections to beluga ini
       60638b9 Debug
       f110918 Debug
       e482afa Debug
       d77e316 Debug
       11ff7bc First commit on gathering metrics for report
       51df068 Merged dev into covseq_stretenp
       58fd24f Merge
       cc0c761 Fixing minor bug and typo
       2376a97 Update READMEs
       8301531 Updating inis and fixing homer output folder location
       706ccb3 Debug
       c59a207 Debug
       bf45b1d Debug
       7b30adf Debug
       e79e74a Debug
       c78d0ea started to edit the hicseq.py added new file for hicrep
       090de77 Fixing minor bug and typo
       01d52a2 Fixing multiqc dependencies
       8d17dce Debug
       f0a2d01 Debug
       9dc813e Debug
       a3f06f2 Debug
       f8dd651 Debug
       d831843 Debug
       9f8d02b Debug
       39463cd Debug
       daebd07 Debug
       898a523 Debug
       3bdc004 Debug
       333524b Debug
       c017225 Debug
       a75dc5a Changing name of header for a report table
       1eb26a7 Debug
       9eb8c95 Typo
       cbd8f1e Addinx x permission to job2json.py
       6f6a381 Fixing report naming too long
       9c91a7e Fixing report naming too long
       d95ab92 Fixing report naming too long
       ce26275 Fixing mugqicValidator.py for new readset format (chipseq only)
       27598c0 Switching to latest mugqic_tools release
       fe2ab50 Changing Ressources
       7401023 Changing Ressources
       7e11212 Changing Ressources
       0a52ab6 Changing module versions
       3b027e9 Debug
       e5de8a7 Fix
       d1304d0 Increasing picard collect multiple metrics RAM
       fff6438 Debug
       8a0463b Debug
       cb4a7e3 Fix ressources
       268550b Fixing annotation peaks
       305bec9 Fixing annotation peaks
       7c56538 Fixing annotation peaks
       9409017 Adding annotation peaks merge for Narrow peaks only
       8420a20 Iterating only on Narrow peaks for annotation
       2d4e2e7 Iterating only on Narrow peaks for annotation
       366142d Debug report
       9253d10 Debug
       aac7c9c Debug
       9c6dc24 Debug
       c9b7191 Debug
       4ac6ca3 Renaming IHEC metrics section
       36880b9 Debug
       c3858f7 Debug
       b1b68b9 Debug Report
       08b77c9 Debug Report
       916d53f Debug
       745fbda Debug
       0c8164f Fixing report
       c918cb1 Debug
       f578ce1 Debug
       ec566be Debug
       a3b196d Debug
       7768187 Debug
       ddf7588 Debug
       4944d36 Debug
       e7d8d28 Debug
       b9f35c3 Debug
       847dddc Debug
       66a5ca1 Debug
       93ad22e Debug
       9d3cdd8 Changing Output metrics
       723f9ed Debug
       fbf10b9 Debug
       a68f0a6 Debug
       68afdd3 Debug
       1ac59b2 Debug
       edcfbbf Debug
       d5548d1 Debug
       a3be70d Debug
       3b9d879 Debug
       3b5d819 Debug
       9215f33 Debug
       1d7732f Debug
       096571a Debug
       1e22dc8 Debug
       28bbc7f Changing macs2 atacseq param & add macs2 param
       982fb12 Debug
       a86d776 Debug
       9a86cbd Debug
       5c6e6b4 Debug
       f778cc3 Debug
       05187de Debug
       2732d26 Debug
       54f47a2 Debug
       4692237 Debug
       187edbe Debug
       6392d86 Debug
       493624f Debug
       553b9ba Debug
       eca0813 Debug
       721c8d1 Debug
       619d6e4 Debug
       773569d Debug
       7836d83 Debug
       8ab1578 Debug
       db9c786 Debug
       2e8f5d2 Changing R version for mugqic_tools R scripts
       c90cec9 homer_annotate_peaks
       d458666 qc_metrics report
       a4c0fbf Fix ihec_metrics outputs
       66af528 Fix ihec_metrics outputs
       2bca481 Fix ihec_metrics outputs
       d33b9e2 Fixing MultiQC report
       cab6ca6 Fixing MultiQC report
       fe600e9 Fixing MultiQC report
       5754975 Fixing MultiQC report
       f31bf71 Fixing MultiQC report
       fcdc378 Fixing MultiQC report
       d3559e6 Improving MultiQC report
       f424ffa Debug
       a1dcb44 Debug
       44b2647 Debug
       fc18819 Debug
       777c716 Debug
       1217e0b Debug
       cbfcd26 Debug
       a5a3ca3 Debug
       479504f Debug
       44a9b80 Debug
       425dc3c Debug
       a853292 Debug
       eaec07b Debug
       f80b4f2 Debug
       81f990a Debug
       ecc3428 Debug
       89b57b8 Debug
       36e4a2a Debug
       bbf0fdc Debug
       d4d8d64 Debug
       3352180 Debug
       9850046 Debug
       1f291ca Major changes IHEC metrics test
       047fe5f Major changes IHEC metrics test
       9240746 Major changes test
       23c9256 Macs2 report changes
       acfc782 Major change test
       76bdd83 Major change test
       392234f Major change test
       a27e5b6 debug
       d54a590 debug
       4ab48db debug
       cfaad21 debug
       5daadd7 debug
       4ebada3 debug
       5a16691 debug
       807e264 debug
       fd604bd debug
       c6347e7 debug
       aa9e88a debug
       440c290 debug
       96d01f1 debug
       e525356 debug
       7785d11 debug
       dc09b15 debug
       a5c60c6 debug
       7be46d4 debug
       0faaea9 debug
       4db5248 debug
       74b73f3 debug
       2df9a1c debug
       6b583d5 debug
       e3961cb debug
       622832e debug
       8b2d57a debug
       fa973d0 debug
       64fd4ed debug
       6a9071f debug
       ab1f024 debug
       40cb701 Fix test
       f1dce58 Major readset change test
       349a1a3 Fix
       095b13e Fix
       17519c7 Fix
       0baff78 Filtering after merging and changing naming to fit with other pipelines
       cab4a18 Increasing default resources and adding options for markdup
       ffae951 Fixing beluga ini
       33c1ab9 Switching to sambamba markdup for IHEC
       b689ca7 Fix
       9a59883 Fix
       c23b193 Fix
       9d1aa0a Fix
       fb2ad2a Fix
       879f850 Fix
       8cc4c62 Options becomes optional
       7cd7a81 Fix
       804de5c Fixing typo
       9795af9 Fix
       11cd56d Fixing sambamba merge
       fb4bbcf Typo
       823693c Adding mkdir
       11e06a2 Fixing temp typo to tmp
       2738637 Fixing job
       6fe4446 Fixing bash import
       8fbe44a Fixing minor bug and typo
       d9dbf29 Adding missing sections to beluga ini
       5f74c2f Renaming freebayes consensus and ncovtools quick_align
       4e6a1bd Degrouping rename consensus jobs because of quast issue creating infinite dependency loop
       25f9e36 Debug
       a0ed01a Debug
       07ab3d1 Adding missing sections to beluga ini
       cc21ed9 Debug
       a113c66 Merged dev into mgi_stretenp
       490857b hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       c844cea hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       4a793a0 started to edit the hicseq.py added new file for hicrep
       0e991a7 started to edit the hicseq.py added new file for hicrep
       72c9387 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9458e13 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ac314f6 hicseq completed adding basic features of the hicrep analysis.
       cd8fbb9 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       086a8cb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       f9893e5 hicseq completed adding basic features of the hicrep analysis.
       4c811af Cleaning
       e87bf28 Debug
       a360050 Debug
       d955925 Changing freebayes calling to skip bgziping
       cdee352 Reducing Ressources
       8856ee9 Debug
       e77bd4b Debug
       30e04a4 Debug
       670c71e Debug
       292c073 Fixing Jared's script
       d799b7a Reducing Ressources
       e2bb1db Fixing Jared's script
       e6f4df8 Fixing freebayes variant calling
       a7e2b43 Fixing Jared's script
       7a6ef99 Fixing bcftools consensus
       5a88dcd Fixing freebayes
       114d9bb Fixing freebayes
       1a3db5e Fixing freebayes
       46d2839 Fixing freebayes/bcftools
       612e506 Adding bcftools consensus creation following Jared's script
       92df26f Adding bcftools consensus creation following Jared's script
       e73f67d Adding freebayes variant calling
       5ce0e18 Grouping rename jobs into 1 job
       ccae339 Grouping rename jobs into 1 job
       f3daf10 Grouping rename jobs into 1 job
       6d850fb Adding freebayes variant calling
       edf5a76 Grouping rename jobs into 1 job
       76fbe6b Adding freebayes variant calling
       7ab5da2 Adding freebayes variant calling
       00acf27 Adding freebayes variant calling
       c700579 Adding freebayes variant calling
       eefed60 Fixing flagging bug and updatinh year to 2021

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       1719a1f Merged in update_container (pull request #238)
       906d0fd Merged in slurm_cgroup_problem (pull request #233)
       1dd042c Merged in temp (pull request #232)

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      387 commits

       773c8d4 hicseq completed adding basic features of the hicrep analysis.
       ffe7bec Completed developing hicrep and quasar analysis
       9e790f6 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6e6b365 Completed hicrep analysis except the creation of the graph
       dbb237f hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       72bde2f hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       167673f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2471a46 hicseq completed adding basic features of the hicrep analysis.
       68aece4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       5903eca Added pairwise combination for samples
       e732259 started to edit the hicseq.py added new file for hicrep
       af7bfd4 hicseq completed adding basic features of the hicrep analysis.
       ad7a803 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       23746b5 Added pairwise combination for samples
       e50669f started to edit the hicseq.py added new file for hicrep
       96ea65e hicseq [hicrep.py] - Corrected typo
       a952c1f hicseq [hicrep.py] - corrected R_TOOLS path
       3dee5c6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2e52c82 hicseq [hicrep.py] - Corrected typo
       fc9ccc4 hicseq [hicrep.py] - corrected R_TOOLS path
       e93c722 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5076045 hicseq [hicrep.py] - Corrected typo
       1404bb4 hicseq [hicrep.py] - corrected R_TOOLS path
       70a9c00 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5e3e5ff hicseq [hicseq.py] - corrected file after rebase
       b96c0e7 hicseq [hicrep.py] - Corrected typo
       dc6f85b hicseq [hicrep.py] - corrected R_TOOLS path
       750cc6c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       92263db hicseq [hicseq.py] - corrected file after rebase
       b708ea3 hicseq [hicrep.py] - Corrected typo
       b209f6a hicseq [hicrep.py] - corrected R_TOOLS path
       3ba0747 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2bf8e95 hicseq [hicrep.py] - Corrected typo
       8936ce6 hicseq [hicrep.py] - corrected R_TOOLS path
       9ac5f10 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1a5bc63 hicseq [hicrep.py] - Corrected typo
       c01ee52 hicseq [hicrep.py] - corrected R_TOOLS path
       1805da8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       19df2ba hicseq [hicrep.py] - Corrected typo
       33e3949 hicseq [hicrep.py] - corrected R_TOOLS path
       77ebe87 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       85294d6 hicseq [hicrep.py] - Corrected typo
       64840a7 hicseq [hicrep.py] - corrected R_TOOLS path
       3013d61 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       734c279 hicseq [hicrep.py] - Corrected typo
       9a33ab8 hicseq [hicrep.py] - corrected R_TOOLS path
       8330611 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       94b562b hicseq [hicrep.py] - Corrected typo
       005b6cc hicseq [hicrep.py] - corrected R_TOOLS path
       28076a5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       858b057 Added pairwise combination for samples
       e76217b started to edit the hicseq.py added new file for hicrep
       8e835eb hicseq [hicseq.py, readme.md] - modified readmefile
       6bf04f8 hicseq [hicseq.py] - corrected file after rebase
       589ae72 hicseq [hicrep.py] - Corrected typo
       3312996 hicseq [hicrep.py] - corrected R_TOOLS path
       2c6c211 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       54108f4 hicseq [hicrep.py] - Corrected typo
       ca68cf4 hicseq [hicrep.py] - corrected R_TOOLS path
       c1d6c9c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       01817b6 hicseq [hicrep.py] - Corrected typo
       a68d82b hicseq [hicrep.py] - corrected R_TOOLS path
       310b5ed [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       785bd9d hicseq [hicrep.py] - Corrected typo
       4db60e2 hicseq [hicrep.py] - corrected R_TOOLS path
       0baffa5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fc107cc hicseq [hicrep.py] - Corrected typo
       8dbe281 hicseq [hicrep.py] - corrected R_TOOLS path
       f272165 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f57e58 hicseq [hicrep.py] - Corrected typo
       11c5af4 hicseq [hicrep.py] - corrected R_TOOLS path
       e544d4f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0c90edc hicseq [hicrep.py] - Corrected typo
       ec0feb8 hicseq [hicrep.py] - corrected R_TOOLS path
       165aeeb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d4d21a3 hicseq [hicrep.py] - Corrected typo
       f665976 hicseq [hicrep.py] - corrected R_TOOLS path
       9a1295e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9778a4a hicseq [hicseq.py] - corrected file after rebase
       10365ba hicseq [hicrep.py] - Corrected typo
       f470412 hicseq [hicrep.py] - corrected R_TOOLS path
       a7f0ba3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0e1ddc6 hicseq [hicrep.py] - Corrected typo
       956a72d hicseq [hicrep.py] - corrected R_TOOLS path
       1cc5498 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9b7e8b3 hicseq [hicrep.py] - Corrected typo
       06cdfcf hicseq [hicrep.py] - corrected R_TOOLS path
       d4c14d2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0ce13df hicseq [hicrep.py] - Corrected typo
       279dc52 hicseq [hicrep.py] - corrected R_TOOLS path
       7c304c7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       81ad33d hicseq [hicrep.py] - Corrected typo
       87bbe66 hicseq [hicrep.py] - corrected R_TOOLS path
       0c5127b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b67b6b0 hicseq [hicrep.py] - Corrected typo
       95a4a97 hicseq [hicrep.py] - corrected R_TOOLS path
       4c25fc2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6d7a744 hicseq [hicrep.py] - Corrected typo
       d96ad6f hicseq [hicrep.py] - corrected R_TOOLS path
       a1768c5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c52093 Added pairwise combination for samples
       5aff10c started to edit the hicseq.py added new file for hicrep
       55987dd hicseq pipeline [changed the step order]
       525c7d3 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       f1811aa hicseq [hicseq.py] - corrected output -o issue fastq_readset
       5fe3a49 hicseq [hicseq.py, readme.md] - modified readmefile
       df1ac4b hicseq [hicseq.py] - corrected file after rebase
       cfd6f52 hicseq [hicrep.py] - Corrected typo
       2b5daac hicseq [hicrep.py] - corrected R_TOOLS path
       796fe27 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e9d2d2d hicseq [hicrep.py] - Corrected typo
       34b9629 hicseq [hicrep.py] - corrected R_TOOLS path
       beb3b07 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       af6d8ac hicseq [hicrep.py] - Corrected typo
       c97468a hicseq [hicrep.py] - corrected R_TOOLS path
       658e574 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5b2fedf hicseq [hicrep.py] - Corrected typo
       f3c653f hicseq [hicrep.py] - corrected R_TOOLS path
       1690ce4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       89371a2 hicseq [hicrep.py] - Corrected typo
       ef56594 hicseq [hicrep.py] - corrected R_TOOLS path
       ff4b8ce [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       161296e hicseq [hicrep.py] - Corrected typo
       61e8b11 hicseq [hicrep.py] - corrected R_TOOLS path
       a147aa9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       db542ae hicseq [hicrep.py] - Corrected typo
       229f065 hicseq [hicrep.py] - corrected R_TOOLS path
       919463a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       077e2ba hicseq [hicrep.py] - Corrected typo
       8c2d82b hicseq [hicrep.py] - corrected R_TOOLS path
       7dd995d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       78ee1f7 hicseq [hicseq.py] - corrected file after rebase
       42120ff hicseq [hicrep.py] - Corrected typo
       e21766c hicseq [hicrep.py] - corrected R_TOOLS path
       2aebe64 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       51e0f2a hicseq [hicrep.py] - Corrected typo
       3d8ad79 hicseq [hicrep.py] - corrected R_TOOLS path
       59ea5d0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       20ea3a8 hicseq [hicrep.py] - Corrected typo
       df2b34a hicseq [hicrep.py] - corrected R_TOOLS path
       6bdaef1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c4727f5 hicseq [hicrep.py] - Corrected typo
       bd80a99 hicseq [hicrep.py] - corrected R_TOOLS path
       7d00494 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2c05619 hicseq [hicrep.py] - Corrected typo
       593a82d hicseq [hicrep.py] - corrected R_TOOLS path
       60581eb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2f0362a hicseq [hicrep.py] - Corrected typo
       a2dfa4f hicseq [hicrep.py] - corrected R_TOOLS path
       04faafc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8556962 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       57dba1b hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       61fce92 hicseq [quasar_qc - corrected module loading in matrix restructuring
       6bc2c0c hicseq [hicrep.py] - Corrected typo
       23217ca hicseq [hicrep.py] - corrected R_TOOLS path
       fc77bfe [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8190866 hicseq [hicrep.py] - Corrected typo
       7089794 hicseq [hicrep.py] - corrected R_TOOLS path
       e5be8fd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       70754df hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       a83bbe3 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5b5ac2f [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b9a9df [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       29fefa5 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9c80cb0 Completed developing hicrep and quasar analysis
       14a6698 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d4f176e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       69310a4 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       38d65fd Completed hicrep analysis except the creation of the graph
       97d7381 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       a7272dd hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       99142a9 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3c16aa5 hicseq completed adding basic features of the hicrep analysis.
       366611d hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       a26866b Added pairwise combination for samples
       54e99f8 started to edit the hicseq.py added new file for hicrep
       ad05137 hicseq [hicseq.py, readme.md] - modified readmefile
       57f13df hicseq [hicseq.py] - corrected file after rebase
       727d579 hicseq [hicrep.py] - Corrected typo
       cabd6b1 hicseq [hicrep.py] - corrected R_TOOLS path
       4f95216 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d9d53d5 hicseq [hicrep.py] - Corrected typo
       2d50a47 hicseq [hicrep.py] - corrected R_TOOLS path
       0b6541e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b56a96a hicseq [hicrep.py] - Corrected typo
       a41dc59 hicseq [hicrep.py] - corrected R_TOOLS path
       d6a751f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       146c892 hicseq [hicrep.py] - Corrected typo
       552d8f9 hicseq [hicrep.py] - corrected R_TOOLS path
       e529d26 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       743948f hicseq [hicrep.py] - Corrected typo
       bed0154 hicseq [hicrep.py] - corrected R_TOOLS path
       c334f10 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b286fcd hicseq [hicrep.py] - Corrected typo
       548cfc0 hicseq [hicrep.py] - corrected R_TOOLS path
       316c44f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       95bb7ed hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2cd7fd6 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       dad2a18 hicseq [quasar_qc - corrected module loading in matrix restructuring
       912cca5 hicseq [hicrep.py] - Corrected typo
       6f1f520 hicseq [hicrep.py] - corrected R_TOOLS path
       8120389 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3816804 hicseq [hicrep.py] - Corrected typo
       8581ac0 hicseq [hicrep.py] - corrected R_TOOLS path
       5371ed0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       226a4f7 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       be5df36 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       f5876b2 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4465ebe [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c5c922f [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       24c34f0 Completed developing hicrep and quasar analysis
       c7d0aa8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7228208 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d1c440a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       b4bc9bf Completed hicrep analysis except the creation of the graph
       324fbc2 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       ff16ded hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       3d54c9f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1e3be79 hicseq completed adding basic features of the hicrep analysis.
       e35ac3d hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ac2c208 Added pairwise combination for samples
       61bd4a7 started to edit the hicseq.py added new file for hicrep
       c950a43 hicseq [hicseq.py] - corrected file after rebase
       9d755fe hicseq [hicrep.py] - Corrected typo
       9cfb4ee hicseq [hicrep.py] - corrected R_TOOLS path
       9bb9705 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b744ced hicseq [hicrep.py] - Corrected typo
       4477d46 hicseq [hicrep.py] - corrected R_TOOLS path
       74fafa0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2a04493 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       acb090f hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       a7578f9 hicseq [quasar_qc - corrected module loading in matrix restructuring
       0e8398a hicseq [hicrep.py] - Corrected typo
       9d27afa hicseq [hicrep.py] - corrected R_TOOLS path
       7fdc12c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8448782 hicseq [hicrep.py] - Corrected typo
       09b7420 hicseq [hicrep.py] - corrected R_TOOLS path
       e5f6849 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cfdb9d4 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5af4d80 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       1ae8080 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e48da77 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       46ba5a8 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       8a20dba Completed developing hicrep and quasar analysis
       3e6f763 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b9ae848 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       60e5df8 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3c2f2 Completed hicrep analysis except the creation of the graph
       e1dea78 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b086ca5 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       114fbd1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       edd76c8 hicseq completed adding basic features of the hicrep analysis.
       3549fd0 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e372026 Added pairwise combination for samples
       8d80a5b started to edit the hicseq.py added new file for hicrep
       18f4e15 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4fbc1fd hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       c12f62a hicseq [quasar_qc - corrected module loading in matrix restructuring
       5e9ddcd hicseq [hicrep.py] - Corrected typo
       1947930 hicseq [hicrep.py] - corrected R_TOOLS path
       c71ed68 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       157c1b5 hicseq [hicrep.py] - Corrected typo
       5937a02 hicseq [hicrep.py] - corrected R_TOOLS path
       b99ff65 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71c1160 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       03f2537 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c940340 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a8e16f9 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       298874c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       bf24afb Completed developing hicrep and quasar analysis
       e5fe17f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e2fc360 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f67be0b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f1b5f5a Completed hicrep analysis except the creation of the graph
       619bbf1 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       c74201f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0b72f78 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f998fc hicseq completed adding basic features of the hicrep analysis.
       a06d726 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       009d286 Added pairwise combination for samples
       9bc074e started to edit the hicseq.py added new file for hicrep
       bcc23d6 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       20949c2 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       81afdcc hicseq [quasar_qc - corrected module loading in matrix restructuring
       08c6757 hicseq [hicrep.py] - Corrected typo
       51c8f66 hicseq [hicrep.py] - corrected R_TOOLS path
       9c2ee30 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       adc92e7 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b020b24 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       623d8bb [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       7e1ecf6 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       6e22be6 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f2815e4 Completed developing hicrep and quasar analysis
       ce39693 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f2e21fc [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4895f68 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       7a13a03 Completed hicrep analysis except the creation of the graph
       0288fb9 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       708eeec hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e97faa5 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1ae2768 hicseq completed adding basic features of the hicrep analysis.
       646858d hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dab3e7d Added pairwise combination for samples
       1e9e150 started to edit the hicseq.py added new file for hicrep
       dd82144 hicseq [hicrep.py] - Corrected typo
       0c8400b hicseq [quasar_qc.py] - Corrected deletion by mistake
       18a704d hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       18747b8 hicseq [hicrep.py] - corrected R_TOOLS path
       5759e58 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cc9a318 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       dcfaa23 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       df458a5 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       6718520 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       ab0c23a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       63ae5fa Completed developing hicrep and quasar analysis
       506c95b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       78c57d8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b87ff3f created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       1825841 Completed hicrep analysis except the creation of the graph
       643905b hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       aa8395f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c285161 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f10a6b hicseq completed adding basic features of the hicrep analysis.
       662a52c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d43a3c1 Added pairwise combination for samples
       3472e7a started to edit the hicseq.py added new file for hicrep
       4a08538 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       e2ede56 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97d248a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       bb6887e [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       815e22a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       854c1c0 Completed developing hicrep and quasar analysis
       c648da3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0f85058 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54cefb7 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3912e Completed hicrep analysis except the creation of the graph
       c1e1035 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       dc414bd hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       943a7c1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       6be727c hicseq completed adding basic features of the hicrep analysis.
       1e5d444 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       3af40a8 Added pairwise combination for samples
       3180796 started to edit the hicseq.py added new file for hicrep
       e67b079 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6e35955 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f0dd719 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       9556bef [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f17594e Completed developing hicrep and quasar analysis
       a384d9a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       339a801 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       aefd4bb created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       98a1cac Completed hicrep analysis except the creation of the graph
       2ae3543 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       255526b hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ced9d13 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       48a9647 hicseq completed adding basic features of the hicrep analysis.
       3b2f387 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       77e52e2 Added pairwise combination for samples
       0894881 started to edit the hicseq.py added new file for hicrep
       ffb9b96 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7507b5a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       21f3a9c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       2c4d8bb [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       72686d1 Completed developing hicrep and quasar analysis
       7cb650c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       abd1855 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       045f823 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       70c513a Completed hicrep analysis except the creation of the graph
       1cf3bb2 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       06f58f8 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       adef366 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       42888b6 hicseq completed adding basic features of the hicrep analysis.
       bf251fd hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d156fd8 Added pairwise combination for samples
       0ff8252 started to edit the hicseq.py added new file for hicrep
       dd17ce7 Completed developing hicrep and quasar analysis
       db31186 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       467e682 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1e11460 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       e9fac87 Completed hicrep analysis except the creation of the graph
       a2afec7 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       234855f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       db9095f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b8f0c55 hicseq completed adding basic features of the hicrep analysis.
       a96d6c1 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       66d63f7 Added pairwise combination for samples
       a69557d started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      55 commits

       b0eb60b update metylseq ini
       7ae859d crosscheck_fingerprint with EXIT_CODE_WHEN_MISMATCH 0
       f2171d2 add new line on genpipes_file.write()
       74130cd line explicite new lines at the end of fp.write()
       ea3c24b fix for pbs
       f8645cc ajust graham and cedar  ini
       129a427 ajust graham ini
       c5a6875 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       b4cad78 cedar and graham ini maintenance
       8c6d99c tweak monitoring
       b81ace0 force -c on tumor pair ini files
       d5e8b07 remove unused import from hicseq
       8a46c26 log_report formating
       b94022d remove the verbose stuff that I am not monitoring anymore.
       6b44921 add missing args* in super call
       06ff02b typo
       6d6b28b update get_wrapper to v2.0.1
       b65dd6e remove decomissioned hpc ini files
       8fcf649 from output_file to genpipes_file
       ba89a19 --command_file options
       99292d8 remove -o option for output file
       96ee305 used image v2.0.0.rc in for container
       3ecb2fb force TMPDIR to default to /tmp if not set
       1e3cc77 make sure output file is executable
       a6e1730 added option --output for command file
       66cf0f5 update form -n to -c
       e10f090 removing useless shebang
       77d1a0d extra space in samtools sort option
       10a9ffc force pace between stringtie options
       047f332 make sure gemini annotation use the ini defined tmp dir
       fbf5552 force log report to python 2 :
       80d9508 error in path for chunk_genpipes csplit call
       b7eb42a force space before bash line continuation
       164e6f8 missing space in varscan script
       1d3bf91 update form -n to -c
       a40dce4 removing useless shebang
       350f898 extra space in samtools sort option
       5e67fa4 force pace between stringtie options
       ba056af make sure gemini annotation use the ini defined tmp dir
       944fcc7 force log report to python 2 :
       a47fe86 error in path for chunk_genpipes csplit call
       9b19e7a force space before bash line continuation
       f50600c missing space in varscan script
       a087236 update form -n to -c
       5a8d512 excluding wget call from monitoring loop
       1dd1bdb force loading mugqic python 2.7 on cuffmerge
       26183ab force loading mugqic python 2.7 on multiqc 1.7
       29735e8 more verbose went -d/--design is needed for contrast
       0a35601 revert on indel aligner
       c47384b cleanup ini for beluga
       c65298e revert on indel aligner
       b18ca36 cleanup ini for beluga
       8514716 remove dbSNP for oxog ini block fix gatk 4 module name
       d514058 fix picard_collect_oxog_metrics dbSNP vcf
       d630a60 remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       019c257 Corrected hicseq.py file deletions
       0c780b0 Corrected hicseq.py file deletions
       9a30832 Corrected hicseq.py file deletions
       8f1db9a Corrected hicseq.py file deletions
       621cdda Corrected hicseq.py file deletions
       a9dce2d Corrected hicseq.py file deletions
       94bb5a8 Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      14 commits

       cabe4d0 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       772bb68 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ef45d3a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       9bda9db hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       93eb03e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a7f50bf hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       f79cacf hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86ddffc hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       7e5a282 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4de276e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3c940bb hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       757f80f hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       56c8875 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86e066c hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      82 commits

       982fdd8 addressed ED and PaulS's feedback on PR- step 3
       10920cb addressed ED and PaulS's feedback on PR- step 2
       498544e addressed ED and PaulS's feedback on PR- step 1
       d18311c bugfix - unsuccesful
       804bc6c added step information to the md file and chipseq.py
       f7365e1 bugfix - unsuccesful
       f5decb5 chipseq.py - remove space in first line
       b85ff24 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       ed99e0a passed more parameters to R script
       a669637 chipseq.py - changed output directory
       06a3f01 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       0893090 adjusted alignment in functions
       b67ba87 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       48e42ca added Rmarkdown rendering added more paramters
       775730d chipseq.py - changed output directory
       1bd300e chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       101e2f6 chipseq.py - remove space in first line
       0aff721 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b84e237 fixed results files creating at custom output folder
       5830eae added Rmarkdown rendering added more paramters
       cd4a5f1 chipseq.py - changed output directory
       8f70abf chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       17b3232 bugfix - unsuccesful
       1d3b463 chipseq.py - remove space in first line
       07ef003 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       80f6cc0 fixed results files creating at custom output folder
       0d13c4f added Rmarkdown rendering added more paramters
       6f50be2 chipseq.py - changed output directory
       81695c6 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       e3132b2 add diff_bind.py in bfx/
       47d33c2 added step information to the md file and chipseq.py
       53dc923 bugfix - unsuccesful
       c2efc7d chipseq.py - remove space in first line
       b6416d6 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       bb677ee passed more parameters to R script
       15f2313 chipseq.py - changed output directory
       cc24da1 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       4855969 adjusted alignment in functions
       196706e skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0843c3a added Rmarkdown rendering added more paramters
       d25e5d0 chipseq.py - changed output directory
       8c17263 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bf900c7 chipseq.py - remove space in first line
       62bd5a5 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b88a2e3 fixed results files creating at custom output folder
       8497911 added Rmarkdown rendering added more paramters
       6c780b2 chipseq.py - changed output directory
       ed75a1e chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bc9d0f2 added step information to the md file and chipseq.py
       43ab77e bugfix - unsuccesful
       994e168 chipseq.py - remove space in first line
       7413ea5 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       241404f fixed results files creating at custom output folder
       a908dc3 added Rmarkdown rendering added more paramters
       26497f5 chipseq.py - changed output directory
       8e8d8e9 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       f71777e add diff_bind.py in bfx/
       1d1d959 bugfix - unsuccesful
       528bd38 chipseq.py - remove space in first line
       7805b91 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       1a7e957 fixed results files creating at custom output folder
       3f58f29 added Rmarkdown rendering added more paramters
       3973234 chipseq.py - changed output directory
       414155c chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       9f12db9 add diff_bind.py in bfx/
       56028a8 chipseq.py - remove space in first line
       662aefa chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       a8d2fcf differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       456d114 fixed results files creating at custom output folder
       92e6fd7 adjusted alignment in functions
       4f66532 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       e09c78f added Rmarkdown rendering added more paramters
       b65cf7e chipseq.py - changed output directory
       303efcd added paramters to ini files
       f927e56 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       abd9dd7 add diff_bind.py in bfx/
       96170df passed more parameters to R script
       28ceb0e chipseq.py - changed output directory
       696560f added paramters to ini files
       f46eb29 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       b926248 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       641d393 add diff_bind.py in bfx/

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      14 commits

       eabb5bb refer previous commit
       5f8957f refer previous commit
       c225c6d refer previous commit
       1a091f5 refer previous commit
       540862d refer previous commit
       818f474 refer previous commit
       21ce258 refer previous commit
       94281fc refer previous commit
       c7774f2 refer previous commit
       d005007 refer previous commit
       a1cf1c0 refer previous commit
       4adc9ca refer previous commit
       3d25f60 refer previous commit
       25cddb0 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      63 commits

       e34a5f3 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       94b62b9 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       304b534 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       580338d fixed results files creating at custom output folder
       8fbd60e added paramters to ini files
       06b9ddf started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       ec6762a passed more parameters to R script
       58aa030 added paramters to ini files
       99101ea started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       3cc8622 rebasing
       1ed1079 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       6187041 adjusted alignment in functions
       c6c31f2 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       038c730 passed more parameters to R script
       5502008 rebasing and resolving conflicts
       bb6343d added paramters to ini files
       ee7af7d started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       f1498d8 rebasing
       745719b chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       30ab714 adjusted alignment in functions
       d86864f skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       1699bae passed more parameters to R script
       9ba71cf rebasing and resolving conflicts
       ebc7a88 added paramters to ini files
       3055bd9 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       bfe5a10 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       6373577 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       044565e fixed results files creating at custom output folder
       ae1abe7 added paramters to ini files
       ad27f00 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       6878212 passed more parameters to R script
       41ba25a added paramters to ini files
       b350cf2 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2010b54 rebasing
       993dfa2 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ffd8154 adjusted alignment in functions
       80dcb4c skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       186b7c7 passed more parameters to R script
       ac4e38d rebasing and resolving conflicts
       d412a8a added paramters to ini files
       428296b started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2d681ff Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       36aa869 rebasing
       20ea21c chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       8949c9e adjusted alignment in functions
       f9c99d1 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       ac7dbf0 passed more parameters to R script
       3c25313 rebasing and resolving conflicts
       e4d75d2 added paramters to ini files
       4c3b940 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       d8c2155 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       0d24cc4 rebasing
       35fff1f chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ca86b46 adjusted alignment in functions
       1076487 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0bb499b passed more parameters to R script
       716493d rebasing and resolving conflicts
       1e38ab5 added paramters to ini files
       4e22bde started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       78edf35 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis resolved merge conflicts after pull
       0ff1cfd passed more parameters to R script
       3195dec started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       9cceef8 corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       6199e34 Merged in chipseq_diff_analysis (pull request #240)
       ce59006 Merged dev into chipseq_diff_analysis

  Shaloo Shalini <shaloo.shalini@gmail.com>      9 commits

       b8f6ff3 Merged in ss_wf_83 (pull request #229)
       63c5ba1 Merged in ss_wf_82 (pull request #228)
       e38d014 Merged in ss_wf_81 (pull request #227)
       0e0b2a2 Merged in ss_wf_80 (pull request #226)
       afdf55a Merged in ss_wf_79 (pull request #225)
       5526239 Merged in ss_wf_78 (pull request #224)
       cd6f047 Merged in ss_wf_77 (pull request #223)
       02a08a4 Merged in ss_wflow_76 (pull request #222)
       eff450c Merged in ss_wflow_75 (pull request #221)

  shaloo <shalz@hotmail.com>      27 commits

       0c5149a Fixes #83 rnaseq denovo workflow diagram updated for v3.4.0
       116e101 Fixes #82 rnaseq_light workflow updated for v3.4.0
       0e00fcd Fixes #81 rnaseq workflows are now current to 3.4.0
       94b163c Fixes #80 methylseq pipeline workflow is now current with v3.4.0
       2811acf Fixes #79 hicseq pipeline update for v3.4.0
       7679be5 Fixes #78 dnaseq workflow for -t mpileup updated v3.4.0
       37c5ca7 Fixes #77 dnaseq pipeline -t mugqic workflow update for genpipes v3.4.0
       531bbf5 Fixes #76 dnaseq_highcov pipeline workflow updated v3.4.0
       9431515 incorrect update should be for #76 cleaning up
       b5dd068 Fixes #76 dnaseq_higcov pipeline workflow updated in sync gpv3.4.0
       a9998d4 Fixes #75 updated amplicon sequence qiime and dada2 workflows for v3.4.0
       3bd1caa Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       8cbae5e Refs #67 dnaseq -t light option color feedback addressed
       1befe53 Fixes #67 Refs #54 dnaseq -t light workflow schema added
       da22575 Refs #65 chipseq workflow color updated, report links added as per feedback
       0f63c1e Fixes #65 Refs #54 chipseq pipeline updated
       c329847 Refs #68 color added as per feedback
       78895de Fixes #68 Refs #54 nanopore pipeline workflow schema added
       1103637 Refs #66 cleanup after color update
       742de76 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       18010fb Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       dcbf14d Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       97aa943 Refs #53 covseq workflow arrow step 6 -> 9 removed
       2d4a60d Refs #53 Paul's review inputs incorporated
       08e8068 Fixes #53 covseq.py workflow diagram added
       8911c4c Refs #53 added covseq pipeline schema workflow draft under review by Paul
       fb8f6f6 Fixes #51 update covseq pipeline readme to v3.3.0

3.4.0        Thu Apr 29 19:24:01 2021 +0000        784 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       ff50d95 GenPipes - RNASeq : corrected genome_index_folder refencing in star_align

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       2b72afb GenPipes - Release : removing Tumor pair prior release 3.4
       9de7a77 Merge remote-tracking branch 'origin/dev' into release_3.4
       ed27f33 GenPipes - Resources : adding software installation scripts
       8c9f613 GenPipes - Resources : adding software installation scipts
       c2d2418 GenPipes - Resources : adding reference genome installation scripts
       3bc6a80 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       67cc499 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       823f626 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       e42f490 GenPipes - Resources : adding Sus_scrofa genome and ivar & kent installation scripts
       d76314e GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      79 commits

       ccc89c1 Merged in release_3.4 (pull request #219)
       9667e55 Resolve Issue #31
       41c771b Merged dev into eh_cit_correction
       5fd600c GenPipes - BFX : corrected ini section for star_index --outTmpDir
       982a305 GenPipes - BFX : update STAR calls to use --outTmpDir
       7a6549d Merged dev into chipseq_design_change
       0ca1695 Merged dev into chipseq_design_change
       23f6e6d Merged in eh_cit_correction (pull request #207)
       daf9eb8 GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       2e1b75e GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       0db8b44 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       9d5c8af Merged dev into chipseq_design_change
       09b50e3 Merged dev into eh_cit_correction
       af89867 GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       024b251 Merged eh_RNAseq_star_correct into dev
       81e394b GenPipes - RNASeq : corrected genome path in star_align
       0dbd895 Merged eh_fix_callhome_fail_exit_code into dev
       a8bfffd GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b652f60 Merged eh_RNAseq_star_correct into dev
       14f5ca5 VERSION edited online with Bitbucket
       4a2a45d GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       ef26b52 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       01d0090 GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       73f22e3 GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       6ae1139 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       a3781df GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       fde1dbc GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       533b3fa GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       4333752 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       9cda20c GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       ceb9594 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ed1a35b GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       6c61ca5 GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       d250236 GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       a5dc00c GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       1551157 GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       404c931 GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       ba277e4 GernPipes - DNASeq : corrected merge_filter_bcf outputs
       1028ccc VERSION edited online with Bitbucket
       82e6690 Merged in ehenrion/version-edited-online-with-bitbucket-1617908341194 (pull request #205)
       e2249aa VERSION edited online with Bitbucket
       0f75fb6 Merged eh_samtools_cram_output_ini_fix into dev
       74d7c56 GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       a45ed72 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       0da38f2 GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       5ded776 GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       b643725 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       33a4b96 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       053b208 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       da787e9 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       713bd64 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       716110d GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       68c1294 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       26ca812 GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       065df8d GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       98f5391 GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       2d7892d GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       cbe399c GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       b74a476 GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       938318a Merged eh_cit_correction into dev
       4d07678 GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       c095c6c GernPipes - DNASeq : corrected merge_filter_bcf outputs
       c5d36c9 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       0d9ab3e GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       138081b GenPipes - DNASeq SV : fixing delly input error
       79f6e94 Merged in ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004 (pull request #202)
       ce2014c Merged eh_RNAseq_star_correct into ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004
       f1a6ffd GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       f53266b GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8a51f4b Merged in ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543 (pull request #200)
       e0cf3d8 Merged eh_RNAseq_star_correct into ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543
       66c603f GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       77ced18 GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       176669b Merged in eh_fix_delly_issue52 (pull request #199)
       2fde5cb GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       277904f dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       fdbef3a GenPipes - Bug fix : corrected dnaseq.cedar.ini
       5824422 GenPipes - Bug fix : correcting indentation in illumina_run_processing.py
       9cfee2e dnaseq.py edited online with Bitbucket : corrected protocol assgignation

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       a0a8e94 rnaseq.base.ini updated to a newer version of STAR
       22f3dfa Merged in rnaseq_star_update (pull request #195)
       06a8d47 rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      219 commits

       d0835ca Merged in chipseq_design_change (pull request #210)
       f7655cb Merged dev into chipseq_design_change
       ce5c007 Fixing atacseq
       2e0ee42 Merged in chipseq_design_change (pull request #206)
       ef03626 Update READMEs
       376858b Merged dev into chipseq_design_change
       182612e Updating inis and fixing homer output folder location
       f27a211 Merged dev into chipseq_design_change
       0ce5991 Debug
       4c0b11f Debug
       3f6fcc6 Debug
       9cbe87c Debug
       586708f Merge branch 'dev' into chipseq_design_change
       5cb6a40 Debug
       1980df2 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       318284c started to edit the hicseq.py added new file for hicrep
       03e1ebc Fixing minor bug and typo
       4829007 Fixing multiqc dependencies
       9acc4eb Debug
       9b72ff0 Debug
       03c9cb6 Debug
       a159e49 Debug
       2c9cf56 Debug
       76664e8 Debug
       7e90a86 Debug
       4a4ae35 Debug
       5718604 Debug
       4c01435 Debug
       bb6ee3f Debug
       aaf4809 Merged dev into chipseq_design_change
       ae98ed2 Debug
       d105f71 Debug
       9c42079 Changing name of header for a report table
       786ac81 Debug
       d98c7ba Typo
       8b1b565 Addinx x permission to job2json.py
       2de4ae6 Fixing report naming too long
       d53e173 Fixing report naming too long
       a3d412e Merged dev into chipseq_design_change
       40bcb0f Fixing report naming too long
       a399b3d Fixing mugqicValidator.py for new readset format (chipseq only)
       2249d6a Switching to latest mugqic_tools release
       3b77185 Changing Ressources
       34524b8 Changing Ressources
       af360e3 Changing Ressources
       4738bda Changing module versions
       6d5e26a Debug
       860b44d Fix
       ef12823 Merged dev into chipseq_design_change
       6512f2f Changing qualimap bamqc section naming to a parameter in bfx
       9b4d61c Increasing picard collect multiple metrics RAM
       b264e7d Debug
       030fb31 Debug
       8de52dd Fix ressources
       dacaec8 Fixing annotation peaks
       12e5d86 Fixing annotation peaks
       9d9a50f Fixing annotation peaks
       351f5c3 Adding annotation peaks merge for Narrow peaks only
       0a2e098 Iterating only on Narrow peaks for annotation
       fd3bd74 Iterating only on Narrow peaks for annotation
       58da5b0 Debug report
       687969c Debug
       955fe4b Debug
       d7560d9 Debug
       e9e09de Debug
       4b8ab24 Renaming IHEC metrics section
       8f9ee56 Debug
       b73ce49 Debug
       43bebdf Debug Report
       cbf2b31 Debug Report
       d894324 Debug
       5dc81e8 Debug
       98d61cf Fixing report
       3b179a4 Debug
       74cd55d Debug
       d3e03b2 Debug
       6545cbd Debug
       9f90aac Debug
       b94c8a1 Debug
       03279f6 Debug
       830e456 Debug
       79cb654 Debug
       19a0024 Debug
       cf9d252 Debug
       505b1ca Debug
       bbf2852 Changing Output metrics
       e654a97 Debug
       43105ff Debug
       2092de5 Debug
       3c6c7b4 Debug
       0ab7349 Debug
       2ee5418 Debug
       a30cf4a Debug
       ec6d257 Debug
       610181d Debug
       a85c974 Debug
       bfd2e2f Debug
       f21e83e Debug
       3fb5418 Debug
       d77307e Debug
       ee397ef Changing macs2 atacseq param & add macs2 param
       8c4274b Debug
       a28d275 Debug
       a64ce53 Debug
       a689539 Debug
       9587ca8 Debug
       e370f82 Debug
       0d3741d Debug
       db78c58 Debug
       f5d866d Debug
       7965625 Debug
       fdf3d9b Debug
       7847ac0 Debug
       47da9ef Debug
       6a8b1ba Debug
       dddcf88 Debug
       161353f Debug
       73790f9 Debug
       188496e Debug
       361b175 Debug
       f5940ed Debug
       e502b6f Changing R version for mugqic_tools R scripts
       45ef768 homer_annotate_peaks
       6a70faf qc_metrics report
       49e8bfe Fix ihec_metrics outputs
       cc35a5d Fix ihec_metrics outputs
       29c1362 Fix ihec_metrics outputs
       7c68071 Fixing MultiQC report
       12212bd Fixing MultiQC report
       c46e42d Fixing MultiQC report
       73ced5d Fixing MultiQC report
       08c1f43 Fixing MultiQC report
       b403bfb Fixing MultiQC report
       bdc3e46 Improving MultiQC report
       cde658a Debug
       528f047 Debug
       cccc2e9 Debug
       a545b4a Debug
       2b79de3 Debug
       acd95d2 Debug
       c0d1a39 Debug
       92d584e Debug
       0c54b4e Debug
       c919c90 Debug
       e9c5955 Debug
       64657e3 Debug
       69857a4 Debug
       2683f88 Debug
       e8a45bf Debug
       2570791 Debug
       6bf3d10 Debug
       524a5c0 Debug
       985fa27 Debug
       3b6ce4a Debug
       ee40fde Debug
       6eedc38 Debug
       34383e7 Major changes IHEC metrics test
       2385c51 Major changes IHEC metrics test
       72f30ff Major changes test
       9dfe33e Macs2 report changes
       705e66a Major change test
       e23426a Major change test
       7d787c1 Major change test
       58b992c debug
       b972151 debug
       9f4c020 debug
       4a14c42 debug
       e69cf40 debug
       20a8959 debug
       d800207 debug
       d0cebc1 debug
       38616a4 debug
       3f0fea1 debug
       fcca285 debug
       61ff7d3 debug
       d147d39 debug
       10e3030 debug
       70c0c5d debug
       acc9d2c debug
       68ea5f4 debug
       a19d64f debug
       76f5155 debug
       8a92259 debug
       9e6c534 debug
       3b9bf7a debug
       1e35666 debug
       b2f165b debug
       bf55791 debug
       4677380 debug
       5cc5350 debug
       f0d9f38 debug
       b306bbd debug
       1d7d871 debug
       8912342 Fix test
       166b57b Major readset change test
       8db202b Fix
       abfcb3d Fix
       ae7964a Fix
       c9bc750 Filtering after merging and changing naming to fit with other pipelines
       aae4609 Increasing default resources and adding options for markdup
       9f540c4 Fixing beluga ini
       dd795e6 Switching to sambamba markdup for IHEC
       27e1a08 Fix
       ecc127d Fix
       e01ccf2 Fix
       b299062 Fix
       f8cf122 Fix
       264059c Fix
       df2604b Options becomes optional
       1748eb6 Fix
       c09bf2f Fixing typo
       adf2e7f Fix
       6c4cbd3 Fixing sambamba merge
       9cf5988 Typo
       fc6f4cb Adding mkdir
       f060825 Fixing temp typo to tmp
       aca35ab Fixing job
       95e5661 Fixing bash import
       f550811 Fixing minor bug and typo

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      385 commits

       03dbf55 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2bd9808 hicseq completed adding basic features of the hicrep analysis.
       7c19ba7 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3b8a19 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       32ffa56 hicseq [hicrep.py] - Corrected typo
       5321e34 hicseq [hicrep.py] - corrected R_TOOLS path
       7d56a21 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fe8103a started to edit the hicseq.py added new file for hicrep
       9f24b59 hicseq [hicrep.py] - Corrected typo
       78d9e59 hicseq [hicrep.py] - corrected R_TOOLS path
       ac8395d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       56f9aa3 hicseq [hicrep.py] - Corrected typo
       f2f843b hicseq [hicrep.py] - corrected R_TOOLS path
       2e3a2e8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b242698 hicseq [hicseq.py] - corrected file after rebase
       b9c2899 hicseq [hicrep.py] - Corrected typo
       56fb9eb hicseq [hicrep.py] - corrected R_TOOLS path
       bd160f6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3726aea Completed developing hicrep and quasar analysis
       ef6f77b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e342 Completed hicrep analysis except the creation of the graph
       e28dd55 started to edit the hicseq.py added new file for hicrep
       5d82058 hicseq [hicseq.py] - corrected file after rebase
       9f494ca hicseq [hicrep.py] - Corrected typo
       a5f7a90 hicseq [hicrep.py] - corrected R_TOOLS path
       e098a61 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3defd1e hicseq [hicrep.py] - Corrected typo
       3acb37d hicseq [hicrep.py] - corrected R_TOOLS path
       f18061a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ec44a6 hicseq [hicrep.py] - Corrected typo
       ac36fca hicseq [hicrep.py] - corrected R_TOOLS path
       911d669 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8b851f3 hicseq [hicrep.py] - Corrected typo
       aee4c68 hicseq [hicrep.py] - corrected R_TOOLS path
       a4ac85c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       766f8db hicseq [hicrep.py] - Corrected typo
       8f9a52d hicseq [hicrep.py] - corrected R_TOOLS path
       e2a85e8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eca4b35 hicseq [hicrep.py] - Corrected typo
       cae9d46 hicseq [hicrep.py] - corrected R_TOOLS path
       c742e6c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9bd1c6e hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       ea5095f hicseq [hicrep.py] - Corrected typo
       0d7d6c0 hicseq [hicrep.py] - corrected R_TOOLS path
       99e1052 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5deb559 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       31d1b1b hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dc822e8 hicseq [hicseq.py, readme.md] - modified readmefile
       b828c72 hicseq [hicseq.py] - corrected file after rebase
       142dfc3 hicseq [hicrep.py] - Corrected typo
       1fd8238 hicseq [hicrep.py] - corrected R_TOOLS path
       b0d0045 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ffaba4c hicseq [hicrep.py] - Corrected typo
       16866fb hicseq [hicrep.py] - corrected R_TOOLS path
       5c52f5a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ca8542 hicseq [hicrep.py] - Corrected typo
       07ce5ba hicseq [hicrep.py] - corrected R_TOOLS path
       9886c9e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       83e1314 hicseq [hicrep.py] - Corrected typo
       4ee15eb hicseq [hicrep.py] - corrected R_TOOLS path
       8ae9494 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7de7736 hicseq [hicrep.py] - Corrected typo
       cc6564d hicseq [hicrep.py] - corrected R_TOOLS path
       ebcc1a1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1eafa4 hicseq [hicrep.py] - Corrected typo
       380dfd4 hicseq [hicrep.py] - corrected R_TOOLS path
       4d9f853 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c2f186d hicseq [hicrep.py] - Corrected typo
       654198f hicseq [hicrep.py] - corrected R_TOOLS path
       c96c475 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bf2449c hicseq [hicrep.py] - Corrected typo
       342b6f5 hicseq [hicrep.py] - corrected R_TOOLS path
       5831718 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f6e3c5 hicseq [hicseq.py] - corrected file after rebase
       bac4d8d hicseq [hicrep.py] - Corrected typo
       8168840 hicseq [hicrep.py] - corrected R_TOOLS path
       e4fdcd0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7250a4f hicseq [hicrep.py] - Corrected typo
       3942f56 hicseq [hicrep.py] - corrected R_TOOLS path
       10b31a8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b6805e6 hicseq [hicrep.py] - Corrected typo
       03ed962 hicseq [hicrep.py] - corrected R_TOOLS path
       fa3d60b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9c4734f hicseq [hicrep.py] - Corrected typo
       dccfad0 hicseq [hicrep.py] - corrected R_TOOLS path
       907e932 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b3e8d9a hicseq [hicrep.py] - Corrected typo
       cbc0e99 hicseq [hicrep.py] - corrected R_TOOLS path
       a15ba0c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3a6a135 hicseq [hicrep.py] - Corrected typo
       7d17401 hicseq [hicrep.py] - corrected R_TOOLS path
       3c495bf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e833365 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       5dd2d25 hicseq [hicrep.py] - Corrected typo
       62b3943 hicseq [hicrep.py] - corrected R_TOOLS path
       9a8f22a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d95fe9f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff71f6f hicseq completed adding basic features of the hicrep analysis.
       35c3269 Added pairwise combination for samples
       85571a9 started to edit the hicseq.py added new file for hicrep
       72be006 hicseq pipeline [changed the step order]
       0800b92 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       3363d07 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       e5b3070 hicseq [hicseq.py, readme.md] - modified readmefile
       d4dbe19 hicseq [hicseq.py] - corrected file after rebase
       adf7ceb hicseq [hicrep.py] - Corrected typo
       9f9e071 hicseq [hicrep.py] - corrected R_TOOLS path
       8481504 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d286d5d hicseq [hicrep.py] - Corrected typo
       8f8c612 hicseq [hicrep.py] - corrected R_TOOLS path
       85a1cf9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8eca7d4 hicseq [hicrep.py] - Corrected typo
       9229a93 hicseq [hicrep.py] - corrected R_TOOLS path
       5a60eff [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef77769 hicseq [hicrep.py] - Corrected typo
       467b8e9 hicseq [hicrep.py] - corrected R_TOOLS path
       c0559eb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c8ff004 hicseq [hicrep.py] - Corrected typo
       6e1367f hicseq [hicrep.py] - corrected R_TOOLS path
       2266dd9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       999455d hicseq [hicrep.py] - Corrected typo
       240eaf3 hicseq [hicrep.py] - corrected R_TOOLS path
       f2ab622 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d075535 hicseq [hicrep.py] - Corrected typo
       28e544c hicseq [hicrep.py] - corrected R_TOOLS path
       5312b56 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4161b94 hicseq [hicrep.py] - Corrected typo
       7ddd445 hicseq [hicrep.py] - corrected R_TOOLS path
       5a24278 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1dd26f4 hicseq [hicseq.py] - corrected file after rebase
       3d526b7 hicseq [hicrep.py] - Corrected typo
       10196c8 hicseq [hicrep.py] - corrected R_TOOLS path
       66cf3f7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2b72cc9 hicseq [hicrep.py] - Corrected typo
       0a1ae4b hicseq [hicrep.py] - corrected R_TOOLS path
       ec8ec06 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       52ceb83 hicseq [hicrep.py] - Corrected typo
       5691c94 hicseq [hicrep.py] - corrected R_TOOLS path
       3390762 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       071ed29 hicseq [hicrep.py] - Corrected typo
       6f23e5f hicseq [hicrep.py] - corrected R_TOOLS path
       29dd059 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bdd9b87 hicseq [hicrep.py] - Corrected typo
       ee032b0 hicseq [hicrep.py] - corrected R_TOOLS path
       0325e11 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f9cb9ba hicseq [hicrep.py] - Corrected typo
       e2f2d6b hicseq [hicrep.py] - corrected R_TOOLS path
       5c92379 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0f4de0b hicseq [base.ini] - updated mugqic tools version to 2.3.1
       f1361d6 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       fe82858 hicseq [quasar_qc - corrected module loading in matrix restructuring
       589286e hicseq [hicrep.py] - Corrected typo
       9dcbc83 hicseq [hicrep.py] - corrected R_TOOLS path
       85b4bbd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4cc5884 hicseq [hicrep.py] - Corrected typo
       26d8fde hicseq [hicrep.py] - corrected R_TOOLS path
       689a76b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5529f4d hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       2d0af16 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c60cd15 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a5cfa02 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c34a220 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       89a55c1 Completed developing hicrep and quasar analysis
       0b74a0a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54c0ea8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       45c545a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       3ea3286 Completed hicrep analysis except the creation of the graph
       f6f2a51 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       034732a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       777005f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       36ab2e9 hicseq completed adding basic features of the hicrep analysis.
       f9a4ef6 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7fff19e Added pairwise combination for samples
       5ae1358 started to edit the hicseq.py added new file for hicrep
       875a2af hicseq [hicseq.py, readme.md] - modified readmefile
       afb029f hicseq [hicseq.py] - corrected file after rebase
       ae2799d hicseq [hicrep.py] - Corrected typo
       22c4aea hicseq [hicrep.py] - corrected R_TOOLS path
       a8b525d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5a84543 hicseq [hicrep.py] - Corrected typo
       29ef577 hicseq [hicrep.py] - corrected R_TOOLS path
       8e99a72 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0af5ff hicseq [hicrep.py] - Corrected typo
       2c514bf hicseq [hicrep.py] - corrected R_TOOLS path
       9c984a0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65e4495 hicseq [hicrep.py] - Corrected typo
       d6b6bda hicseq [hicrep.py] - corrected R_TOOLS path
       15c7831 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a98b326 hicseq [hicrep.py] - Corrected typo
       bfc0ee4 hicseq [hicrep.py] - corrected R_TOOLS path
       2dd0caa [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6070378 hicseq [hicrep.py] - Corrected typo
       b7f35ae hicseq [hicrep.py] - corrected R_TOOLS path
       c739410 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       473304e hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b536a37 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2e9b266 hicseq [quasar_qc - corrected module loading in matrix restructuring
       e40bf89 hicseq [hicrep.py] - Corrected typo
       ec4b642 hicseq [hicrep.py] - corrected R_TOOLS path
       6ebce20 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8a7e4f1 hicseq [hicrep.py] - Corrected typo
       4ec40e0 hicseq [hicrep.py] - corrected R_TOOLS path
       d7e39f9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1cc96d3 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5e52fc3 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       59c65ac [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fdfbe1 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       e5eb097 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e95cd4d Completed developing hicrep and quasar analysis
       c5c7908 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       08ff8a7 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d6abcd4 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a6156a1 Completed hicrep analysis except the creation of the graph
       4eea55c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       43501ba hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       25c34b4 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a164d8b hicseq completed adding basic features of the hicrep analysis.
       e8aa04c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       491fa1c Added pairwise combination for samples
       b42241a started to edit the hicseq.py added new file for hicrep
       20627d6 hicseq [hicseq.py] - corrected file after rebase
       884d4e4 hicseq [hicrep.py] - Corrected typo
       622d5f3 hicseq [hicrep.py] - corrected R_TOOLS path
       c4a9c60 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fb0ec53 hicseq [hicrep.py] - Corrected typo
       7082fd6 hicseq [hicrep.py] - corrected R_TOOLS path
       f90b5be [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       02600a9 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       456e649 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d96bfb1 hicseq [quasar_qc - corrected module loading in matrix restructuring
       86c1235 hicseq [hicrep.py] - Corrected typo
       79fca77 hicseq [hicrep.py] - corrected R_TOOLS path
       c28c3f8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef796f9 hicseq [hicrep.py] - Corrected typo
       3037557 hicseq [hicrep.py] - corrected R_TOOLS path
       b5b40fe [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e538fb0 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       28dee2f [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d7e7c27 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       82aa89a [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       0032036 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db02527 Completed developing hicrep and quasar analysis
       edf03c1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7b1b760 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2520bb9 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       314a505 Completed hicrep analysis except the creation of the graph
       71f75b6 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       723cecf hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       38ec796 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       fb02b11 hicseq completed adding basic features of the hicrep analysis.
       5664879 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       696a52b Added pairwise combination for samples
       886ea4e started to edit the hicseq.py added new file for hicrep
       10d4d24 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       9d7b865 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2c3b54d hicseq [quasar_qc - corrected module loading in matrix restructuring
       dbdf0ca hicseq [hicrep.py] - Corrected typo
       e70f0ee hicseq [hicrep.py] - corrected R_TOOLS path
       92a5667 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c05b038 hicseq [hicrep.py] - Corrected typo
       1d57da7 hicseq [hicrep.py] - corrected R_TOOLS path
       09b27e0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       565b4f4 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       f9f5282 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       79da62c [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4518e45 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       11e583d [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f7d773d Completed developing hicrep and quasar analysis
       b7f60de [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b3af696 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0034cea created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2a5492d Completed hicrep analysis except the creation of the graph
       d6e8b9b hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       f1f7d98 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ba1d617 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       76e10ba hicseq completed adding basic features of the hicrep analysis.
       72a3137 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c2af6c6 Added pairwise combination for samples
       6f0965b started to edit the hicseq.py added new file for hicrep
       e9f3fe4 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2c206ac hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3e2a53e hicseq [quasar_qc - corrected module loading in matrix restructuring
       f880a49 hicseq [hicrep.py] - Corrected typo
       809f95e hicseq [hicrep.py] - corrected R_TOOLS path
       aaddec5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       251c0c6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       68889f9 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       68003e0 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fef220 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8782a74 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e2caf44 Completed developing hicrep and quasar analysis
       30008d6 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7a64058 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       02e39e7 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ea8610d Completed hicrep analysis except the creation of the graph
       48173c5 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       419d60e hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       b26f623 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f2eba5 hicseq completed adding basic features of the hicrep analysis.
       92b291a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3dcb3e Added pairwise combination for samples
       67488a0 started to edit the hicseq.py added new file for hicrep
       a22ca56 hicseq [hicrep.py] - Corrected typo
       e28aede hicseq [quasar_qc.py] - Corrected deletion by mistake
       5a1f502 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       5b6fff8 hicseq [hicrep.py] - corrected R_TOOLS path
       b02d228 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4ed3b59 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       c361bc6 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5622389 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e4758a1 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c564149 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       760bfd6 Completed developing hicrep and quasar analysis
       275a224 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       627a501 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       33ffee0 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       df084a8 Completed hicrep analysis except the creation of the graph
       f3bebd7 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b7fe382 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       805c4a9 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       dbfde76 hicseq completed adding basic features of the hicrep analysis.
       10e2dfa hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7c2e89b Added pairwise combination for samples
       07099a0 started to edit the hicseq.py added new file for hicrep
       300ba1b hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       688a495 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9aa6c61 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       0f393c6 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       b64480f [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       353a55d Completed developing hicrep and quasar analysis
       7aa1f48 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       008904a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b27275b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       66c49d4 Completed hicrep analysis except the creation of the graph
       a09f038 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5c370b5 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       de693ae hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a504b17 hicseq completed adding basic features of the hicrep analysis.
       fe551f3 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       5220816 Added pairwise combination for samples
       17e1b5c started to edit the hicseq.py added new file for hicrep
       8d592db [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c37d365 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8c15244 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       7f55aeb [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f371fb9 Completed developing hicrep and quasar analysis
       14b332d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       66ded90 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c8cabe6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f066550 Completed hicrep analysis except the creation of the graph
       0a733f0 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       445a468 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ebcf685 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8e3602c hicseq completed adding basic features of the hicrep analysis.
       21b08c3 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       59771c3 Added pairwise combination for samples
       9819482 started to edit the hicseq.py added new file for hicrep
       8fc473c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       673bad3 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       979085b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       5392e00 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       96b5c41 Completed developing hicrep and quasar analysis
       8d38fc1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e517 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dade053 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       14247d2 Completed hicrep analysis except the creation of the graph
       79e675e hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5912593 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c19d1be hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       cc2f34c hicseq completed adding basic features of the hicrep analysis.
       e7e729a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8495df6 Added pairwise combination for samples
       d789632 started to edit the hicseq.py added new file for hicrep
       d8a1537 Completed developing hicrep and quasar analysis
       c13cb74 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       316db2c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a382fb created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       68fc66e Completed hicrep analysis except the creation of the graph
       f2cb6a0 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       702fadb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4fefe1e hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       294e3aa hicseq completed adding basic features of the hicrep analysis.
       5edc80e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       fc95853 Added pairwise combination for samples
       1592510 started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      10 commits

       4b08c93 force loading mugqic python 2.7 on cuffmerge
       ac7cd86 force loading mugqic python 2.7 on multiqc 1.7
       b0875ae more verbose went -d/--design is needed for contrast
       0dc71a9 revert on indel aligner
       37af570 cleanup ini for beluga
       3c0fad0 revert on indel aligner
       5463567 cleanup ini for beluga
       2080dfd remove dbSNP for oxog ini block fix gatk 4 module name
       3a1ab22 fix picard_collect_oxog_metrics dbSNP vcf
       1c5e2ed remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       7bbe4fa Corrected hicseq.py file deletions
       0ca8a12 Corrected hicseq.py file deletions
       e2e57af Corrected hicseq.py file deletions
       9d90aef Corrected hicseq.py file deletions
       312d3db Corrected hicseq.py file deletions
       a95d8b8 Corrected hicseq.py file deletions
       52c9240 Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      11 commits

       39db7bb hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4943274 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       452015d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       61a86ec hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       e798549 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       185c81b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21929bf hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       fe03ec6 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cdca1f9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a146b25 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4025c71 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      11 commits

       c646e67 refer previous commit
       30b0838 refer previous commit
       0e4a6cb refer previous commit
       138fafc refer previous commit
       0fd5742 refer previous commit
       962cc3e refer previous commit
       bb6b1b7 refer previous commit
       a1f0ba2 refer previous commit
       2b79a1b refer previous commit
       ae7e5c8 refer previous commit
       3b6c393 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       30b2d60 corrected md file

  Shaloo Shalini <shaloo.shalini@gmail.com>      17 commits

       f0d1045 Merged in ss_issu_69_rnaseq_strintie_wf (pull request #214)
       6935d81 Merged in ss_issu_67_wf_dnaseq_light (pull request #216)
       2b8fe8c Merged in ss_issu_66_dnaseq_sv_schema (pull request #215)
       17c8fcc Merged dev into ss_issu_69_rnaseq_strintie_wf
       0b1fe1e Merged dev into ss_issu_66_dnaseq_sv_schema
       0f9def7 Merged dev into ss_issu_67_wf_dnaseq_light
       37504fc Merged in ss_issu_67_wf_dnaseq_light (pull request #212)
       b217c76 Merged in ss_issu_65_chipseq_workflow (pull request #211)
       890dc40 Merged in ss_issu_68_nano_wf (pull request #213)
       1c7a288 Merged dev into ss_issu_67_wf_dnaseq_light
       0a6fd27 Merged dev into ss_issu_69_rnaseq_strintie_wf
       c7272d5 Merged dev into ss_issu_68_nano_wf
       0a1ffe3 Merged dev into ss_issu_65_chipseq_workflow
       3f9add6 Merged in ss_issu_66_dnaseq_sv_schema (pull request #209)
       2bddc59 Merged in ss_issu_60_wf_ampllconseq (pull request #203)
       7146106 Merged in ss_covseq_wflow (pull request #197)
       192f7a3 Merged in covseq_readme (pull request #196)

  shaloo <shalz@hotmail.com>      30 commits

       14cb58c Refs #69 Hector and Ed's feedback incorporated
       ea4ea09 Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       de861a5 Refs #66 cleanup after color update
       1b68619 Refs #67 cleanup after color update
       26733d9 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       75101db Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       47c8705 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       dc52d26 Refs #67 dnaseq -t light option color feedback addressed
       0c95d46 Refs #68 color updated as per feedback
       796cf88 Refs #68 color added as per feedback
       1ad517e Refs #65 chipseq workflow color updated, report links added as per feedback
       179955f Fixes #69 Refs #54 rnaseq -t stringtie workflow schema added
       c1a8861 Merge branch 'ss_issu_65_chipseq_workflow' of bitbucket.org:mugqic/genpipes into ss_issu_65_chipseq_workflow
       ce49d39 Fixes #65 Refs #54 chipseq pipeline updated
       8beba66 Fixes #68 Refs #54 nanopore pipeline workflow schema added
       cc91732 Fixes #67 Refs #54 dnaseq -t light workflow schema added
       9f4b8ce Update chipseq pipeline schema for -t chipseq and -t atacseq options to reflect latest dev branch pipeline code
       6f88456 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       53cc1d3 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       e4fd2a2 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       647d39a Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       7b041d9 cleanup DS_Store file
       44db707 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       48626df Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       52a2286 Refs #53 covseq workflow arrow step 6 -> 9 removed
       faf9449 Refs #53 Paul's review inputs incorporated
       5187035 Fixes #53 covseq.py workflow diagram added
       bacd6ed Refs #53 added covseq pipeline schema workflow draft under review by Paul
       0f7b09a Fixes #51 update covseq pipeline readme to v3.3.0
       2fed25b Fixes #51 update covseq pipeline readme to v3.3.0

3.3.0        Fri Feb 19 16:37:40 2021 -0500        641 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       a4d9a14 GenPipes - removing tumor_apir before release
       1a28942 Version bump to 3.2.1-beta
       63d211d Version bump to 3.2.0

  ehenrion <edouard.henrion@mcgill.ca>      4 commits

       0756af3 Merged in release_3.3 (pull request #194)
       7b17e70 GenPipes - Bug fix : corrected dnaseq.cedar.ini
       d62b0c0 GenPipes - Bug fix : removed buggy line in dnaseq.cedar.ini
       30bdd0a GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       ebef4a2 Release 3.2

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      572 commits

       e6d1d8a hicseq pipeline [changed the step order]
       165e060 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       498209b hicseq [hicseq.py] - corrected output -o issue fastq_readset
       8d382dc hicseq [hicseq.py, readme.md] - modified readmefile
       638e16e hicseq [hicseq.py] - corrected file after rebase
       7fe877f hicseq [hicrep.py] - Corrected typo
       2c26336 hicseq [hicrep.py] - corrected R_TOOLS path
       91dd29c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f44f83d hicseq [hicrep.py] - Corrected typo
       124ce0e hicseq [hicrep.py] - corrected R_TOOLS path
       678b9cf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bb10a4e hicseq [hicrep.py] - Corrected typo
       3e2b550 hicseq [hicrep.py] - corrected R_TOOLS path
       8772abd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       77cbfb9 hicseq [hicrep.py] - Corrected typo
       6f820f1 hicseq [hicrep.py] - corrected R_TOOLS path
       2528718 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65191ab hicseq [hicrep.py] - Corrected typo
       28b4377 hicseq [hicrep.py] - corrected R_TOOLS path
       c82045f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6b1b9d6 hicseq [hicrep.py] - Corrected typo
       be461ea hicseq [hicrep.py] - corrected R_TOOLS path
       26cf925 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b06d2ea hicseq [hicrep.py] - Corrected typo
       eff718e hicseq [hicrep.py] - corrected R_TOOLS path
       8c7b0eb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d1f43e7 hicseq [hicrep.py] - Corrected typo
       97174e2 hicseq [hicrep.py] - corrected R_TOOLS path
       d82f677 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       aa9e921 hicseq [hicseq.py] - corrected file after rebase
       96a9420 hicseq [hicrep.py] - Corrected typo
       68f348b hicseq [hicrep.py] - corrected R_TOOLS path
       284f0e1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac094d hicseq [hicrep.py] - Corrected typo
       ec553bd hicseq [hicrep.py] - corrected R_TOOLS path
       e32d3bb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1f326e hicseq [hicrep.py] - Corrected typo
       4b99cbb hicseq [hicrep.py] - corrected R_TOOLS path
       da58a44 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8985879 hicseq [hicrep.py] - Corrected typo
       0d1f31b hicseq [hicrep.py] - corrected R_TOOLS path
       d5aab52 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1d3a959 hicseq [hicrep.py] - Corrected typo
       2ec58eb hicseq [hicrep.py] - corrected R_TOOLS path
       aaf3aee [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2dd3de8 hicseq [hicrep.py] - Corrected typo
       0eb1956 hicseq [hicrep.py] - corrected R_TOOLS path
       c262754 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       944ba8e hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4021b7d hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       e856761 hicseq [quasar_qc - corrected module loading in matrix restructuring
       fc92143 hicseq [hicrep.py] - Corrected typo
       a13ed40 hicseq [hicrep.py] - corrected R_TOOLS path
       ee6bc35 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3b82c00 hicseq [hicrep.py] - Corrected typo
       9a509dc hicseq [hicrep.py] - corrected R_TOOLS path
       edb308a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ca9fa5d hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       d9b7dc4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       aa8821b [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       9419efd [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8d61964 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       2e7870d Completed developing hicrep and quasar analysis
       d54c390 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       3ca0a12 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       845d473 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       42f947a Completed hicrep analysis except the creation of the graph
       0274f32 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       1c1f819 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e27ded7 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5340844 hicseq completed adding basic features of the hicrep analysis.
       5b078cf hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       94230e0 Added pairwise combination for samples
       39b34b0 started to edit the hicseq.py added new file for hicrep
       9e2040f hicseq [hicseq.py, readme.md] - modified readmefile
       4465952 hicseq [hicseq.py] - corrected file after rebase
       ccb344f hicseq [hicrep.py] - Corrected typo
       35df970 hicseq [hicrep.py] - corrected R_TOOLS path
       240f9c5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eecc3b2 hicseq [hicrep.py] - Corrected typo
       78e6c8a hicseq [hicrep.py] - corrected R_TOOLS path
       616e8d3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       84e8ec9 hicseq [hicrep.py] - Corrected typo
       b142de6 hicseq [hicrep.py] - corrected R_TOOLS path
       dd6cb48 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5bad885 hicseq [hicrep.py] - Corrected typo
       030d558 hicseq [hicrep.py] - corrected R_TOOLS path
       e29ae96 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       97622c6 hicseq [hicrep.py] - Corrected typo
       a477d71 hicseq [hicrep.py] - corrected R_TOOLS path
       d05a6b6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       12a4e2f hicseq [hicrep.py] - Corrected typo
       45695cb hicseq [hicrep.py] - corrected R_TOOLS path
       e7041d3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       302c4e9 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       154f9b2 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       8e55fb5 hicseq [quasar_qc - corrected module loading in matrix restructuring
       f1534e3 hicseq [hicrep.py] - Corrected typo
       4557237 hicseq [hicrep.py] - corrected R_TOOLS path
       c2063ea [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e41bb32 hicseq [hicrep.py] - Corrected typo
       ebbbf87 hicseq [hicrep.py] - corrected R_TOOLS path
       4406314 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bd36df0 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       4076547 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6f40129 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       789fb24 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       83ea31b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       0eafabd Completed developing hicrep and quasar analysis
       dc1cf6a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c0974eb [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       82c023b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       0ea49fd Completed hicrep analysis except the creation of the graph
       55f2b85 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bb22368 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       19aec50 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       7cd0785 hicseq completed adding basic features of the hicrep analysis.
       7bdf5b4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ae313c6 Added pairwise combination for samples
       08dae47 started to edit the hicseq.py added new file for hicrep
       eff7114 hicseq [hicseq.py] - corrected file after rebase
       b693105 hicseq [hicrep.py] - Corrected typo
       386ca83 hicseq [hicrep.py] - corrected R_TOOLS path
       5a412f7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0b7fe2a hicseq [hicrep.py] - Corrected typo
       1765bef hicseq [hicrep.py] - corrected R_TOOLS path
       5851499 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1555fed hicseq [base.ini] - updated mugqic tools version to 2.3.1
       d3199ec hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       01d239a hicseq [quasar_qc - corrected module loading in matrix restructuring
       46e7d3a hicseq [hicrep.py] - Corrected typo
       5e92f8e hicseq [hicrep.py] - corrected R_TOOLS path
       f2dbfd7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4fa81e9 hicseq [hicrep.py] - Corrected typo
       b284455 hicseq [hicrep.py] - corrected R_TOOLS path
       b12fdbe [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a78a623 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       125a6e3 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       67272ac [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       19f3edc [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       a8b8d85 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9d37f63 Completed developing hicrep and quasar analysis
       acce383 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       92bfb41 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a3bcdf created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       812ba6d Completed hicrep analysis except the creation of the graph
       27467f3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3cc9f39 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2684371 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       92744be hicseq completed adding basic features of the hicrep analysis.
       06d8e53 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       31b0157 Added pairwise combination for samples
       6775791 started to edit the hicseq.py added new file for hicrep
       38a937e hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4e9dea9 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       6efe674 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b918766 hicseq [hicrep.py] - Corrected typo
       0219397 hicseq [hicrep.py] - corrected R_TOOLS path
       187cba3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac69b2 hicseq [hicrep.py] - Corrected typo
       a80a6bc hicseq [hicrep.py] - corrected R_TOOLS path
       06e2631 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f903c4b hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ebc75da [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c9acbc7 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       dcbe50b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bd63e8e [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4f9a308 Completed developing hicrep and quasar analysis
       6f0d1a8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1052c8d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       704d7f1 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a7f8b2a Completed hicrep analysis except the creation of the graph
       f1fd9af hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       712cb46 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       979429b hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3930cd3 hicseq completed adding basic features of the hicrep analysis.
       b3c1160 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       f493eae Added pairwise combination for samples
       6f98cb3 started to edit the hicseq.py added new file for hicrep
       f58f27f hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b576afd hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       0eb37c3 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b718a5b hicseq [hicrep.py] - Corrected typo
       5f8cf47 hicseq [hicrep.py] - corrected R_TOOLS path
       1fde108 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6478d07 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b262df5 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b340c0d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b6b357 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       d740ed8 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c1b8eef Completed developing hicrep and quasar analysis
       6bc89a8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       95c477e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       53fc03c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       cf6a3af Completed hicrep analysis except the creation of the graph
       d854224 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b96177 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       8438ff6 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b18b191 hicseq completed adding basic features of the hicrep analysis.
       c4ba094 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bd650a8 Added pairwise combination for samples
       6f4136c started to edit the hicseq.py added new file for hicrep
       9364f1e hicseq [hicrep.py] - Corrected typo
       301afce hicseq [quasar_qc.py] - Corrected deletion by mistake
       bd2a8e1 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       9676178 hicseq [hicrep.py] - corrected R_TOOLS path
       a5ca0ba [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       432f586 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       3e68e18 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5c841f7 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       ab730be [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       3bce823 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e1cee2e Completed developing hicrep and quasar analysis
       45b526c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6cf8450 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       21e2e6c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       3441689 Completed hicrep analysis except the creation of the graph
       6570bc9 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       46a28b4 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       645c2f8 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       0189b80 hicseq completed adding basic features of the hicrep analysis.
       2606951 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       6fd22e8 Added pairwise combination for samples
       d664e0e started to edit the hicseq.py added new file for hicrep
       1e51cf6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       4cca061 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2f8c42d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       1c49e45 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cfa04f3 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c0bd365 Completed developing hicrep and quasar analysis
       70d8bdf [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6711a68 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e8effea created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       5fe98a0 Completed hicrep analysis except the creation of the graph
       0bc9a5d hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bea29b6 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       eb89bc0 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       d7bb65d hicseq completed adding basic features of the hicrep analysis.
       ccc376e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       93f957d Added pairwise combination for samples
       1c14ddf started to edit the hicseq.py added new file for hicrep
       0da3f4c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9ae9d5e [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       30a4e55 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8a5ac66 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4fb194f Completed developing hicrep and quasar analysis
       7f2917d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2492118 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9a8cfcd created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       54c7c37 Completed hicrep analysis except the creation of the graph
       c35566c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fe1171a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e9db27d hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8584e29 hicseq completed adding basic features of the hicrep analysis.
       dccac15 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       5689051 Added pairwise combination for samples
       f609b2f started to edit the hicseq.py added new file for hicrep
       31ba309 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       55f7758 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f794061 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       84e93b1 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       915d9ac Completed developing hicrep and quasar analysis
       769c97b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       01ed537 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b8a79f4 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       6163240 Completed hicrep analysis except the creation of the graph
       7efc2ba hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bf7412f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       14256a1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       042082e hicseq completed adding basic features of the hicrep analysis.
       fa1e7d0 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       eba9e46 Added pairwise combination for samples
       9bd17ff started to edit the hicseq.py added new file for hicrep
       99d6222 Completed developing hicrep and quasar analysis
       719115c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7c3222a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       966e814 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ece06e3 Completed hicrep analysis except the creation of the graph
       4ceda72 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a65032 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       91f3b4b hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       15e34d0 hicseq completed adding basic features of the hicrep analysis.
       92baadc hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d9d98a5 Added pairwise combination for samples
       c4853c3 started to edit the hicseq.py added new file for hicrep
       78ab010 hicseq pipeline [changed the step order]
       c2ac146 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       8bf6a4d hicseq [hicseq.py] - corrected output -o issue fastq_readset
       a016c1e hicseq [hicseq.py, readme.md] - modified readmefile
       7fbb757 hicseq [hicseq.py] - corrected file after rebase
       f8fc53c hicseq [hicrep.py] - Corrected typo
       46c6016 hicseq [hicrep.py] - corrected R_TOOLS path
       fb83a92 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       00a5b87 hicseq [hicrep.py] - Corrected typo
       3e840af hicseq [hicrep.py] - corrected R_TOOLS path
       c61122b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d3659cb hicseq [hicrep.py] - Corrected typo
       230b7b6 hicseq [hicrep.py] - corrected R_TOOLS path
       9a454ab [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       857f8c9 hicseq [hicrep.py] - Corrected typo
       5c40245 hicseq [hicrep.py] - corrected R_TOOLS path
       ba545bd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3ff1e50 hicseq [hicrep.py] - Corrected typo
       a0098d0 hicseq [hicrep.py] - corrected R_TOOLS path
       2c820d8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d23933c hicseq [hicrep.py] - Corrected typo
       b064564 hicseq [hicrep.py] - corrected R_TOOLS path
       54bb894 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b81d755 hicseq [hicrep.py] - Corrected typo
       74e48a2 hicseq [hicrep.py] - corrected R_TOOLS path
       40fe0cc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a3680d6 hicseq [hicrep.py] - Corrected typo
       798c19f hicseq [hicrep.py] - corrected R_TOOLS path
       06ec788 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2fd36b1 hicseq [hicseq.py] - corrected file after rebase
       be817b4 hicseq [hicrep.py] - Corrected typo
       5248d32 hicseq [hicrep.py] - corrected R_TOOLS path
       39617b2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec97960 hicseq [hicrep.py] - Corrected typo
       24e7c7c hicseq [hicrep.py] - corrected R_TOOLS path
       4d6d3f8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c17d5f1 hicseq [hicrep.py] - Corrected typo
       aeae41c hicseq [hicrep.py] - corrected R_TOOLS path
       77dabaf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       032cfe9 hicseq [hicrep.py] - Corrected typo
       07a72bd hicseq [hicrep.py] - corrected R_TOOLS path
       1039f67 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       590a600 hicseq [hicrep.py] - Corrected typo
       8f9b546 hicseq [hicrep.py] - corrected R_TOOLS path
       574336c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6c20f49 hicseq [hicrep.py] - Corrected typo
       ff45f2d hicseq [hicrep.py] - corrected R_TOOLS path
       85d661e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec67234 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       463da45 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       7701580 hicseq [quasar_qc - corrected module loading in matrix restructuring
       0c5ee2f hicseq [hicrep.py] - Corrected typo
       1a221b4 hicseq [hicrep.py] - corrected R_TOOLS path
       8966cfb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b02d1a7 hicseq [hicrep.py] - Corrected typo
       4d37159 hicseq [hicrep.py] - corrected R_TOOLS path
       22fbd60 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5541fbc hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       fdeb31d [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d58d40f [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       fa1f4b4 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cdcc790 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7b9f646 Completed developing hicrep and quasar analysis
       31c4df5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c1cc394 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       80dc0ab created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       72cbd00 Completed hicrep analysis except the creation of the graph
       0a93826 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       85cc99e hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2e45578 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       46f2fa3 hicseq completed adding basic features of the hicrep analysis.
       8c3e8d9 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8e5c654 Added pairwise combination for samples
       55862cd started to edit the hicseq.py added new file for hicrep
       e79b132 hicseq [hicseq.py, readme.md] - modified readmefile
       ffea320 hicseq [hicseq.py] - corrected file after rebase
       fed6253 hicseq [hicrep.py] - Corrected typo
       d285ff0 hicseq [hicrep.py] - corrected R_TOOLS path
       e04cb39 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       de2d8e1 hicseq [hicrep.py] - Corrected typo
       9c6f2d2 hicseq [hicrep.py] - corrected R_TOOLS path
       a8ae3ab [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       866b11d hicseq [hicrep.py] - Corrected typo
       82db122 hicseq [hicrep.py] - corrected R_TOOLS path
       4b9e5e6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       dda5352 hicseq [hicrep.py] - Corrected typo
       452f1db hicseq [hicrep.py] - corrected R_TOOLS path
       db683cc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f2e3612 hicseq [hicrep.py] - Corrected typo
       2accb63 hicseq [hicrep.py] - corrected R_TOOLS path
       f8413f0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7938420 hicseq [hicrep.py] - Corrected typo
       e78e318 hicseq [hicrep.py] - corrected R_TOOLS path
       80eff20 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0b0062 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       1160fe2 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3ba1174 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b27f402 hicseq [hicrep.py] - Corrected typo
       01ac15d hicseq [hicrep.py] - corrected R_TOOLS path
       80ad708 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef6d034 hicseq [hicrep.py] - Corrected typo
       045b540 hicseq [hicrep.py] - corrected R_TOOLS path
       f6422cd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8808f05 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       65f84c4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7adcfa1 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       14c0a78 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       935862a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       fd6eef8 Completed developing hicrep and quasar analysis
       8c8248e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       85f0d9e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dc225b2 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       dfe7725 Completed hicrep analysis except the creation of the graph
       db50a5d hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       10222fe hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       5dcda07 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       062551b hicseq completed adding basic features of the hicrep analysis.
       766c1b4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e979eee Added pairwise combination for samples
       cfff916 started to edit the hicseq.py added new file for hicrep
       4443e3d hicseq [hicseq.py] - corrected file after rebase
       ca8b09d hicseq [hicrep.py] - Corrected typo
       584dd18 hicseq [hicrep.py] - corrected R_TOOLS path
       afa77dc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c2535f hicseq [hicrep.py] - Corrected typo
       8db0da7 hicseq [hicrep.py] - corrected R_TOOLS path
       6f264bd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5960be4 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       537be82 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       1313e2b hicseq [quasar_qc - corrected module loading in matrix restructuring
       9bcc430 hicseq [hicrep.py] - Corrected typo
       8cacfa8 hicseq [hicrep.py] - corrected R_TOOLS path
       c8c3cef [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71b02a3 hicseq [hicrep.py] - Corrected typo
       b98db26 hicseq [hicrep.py] - corrected R_TOOLS path
       3698d95 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d461e1e hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5c728c6 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7a9cfa6 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2b66572 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       312c08e [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7e6f9aa Completed developing hicrep and quasar analysis
       334abaf [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       87973bb [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       03aab02 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8359e46 Completed hicrep analysis except the creation of the graph
       8a845b7 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       eee8a21 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       9a55b5c hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2755cba hicseq completed adding basic features of the hicrep analysis.
       3c517c4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8852df1 Added pairwise combination for samples
       54e5f9f started to edit the hicseq.py added new file for hicrep
       9d54100 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       95e4b22 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d63f7b4 hicseq [quasar_qc - corrected module loading in matrix restructuring
       56c0d8c hicseq [hicrep.py] - Corrected typo
       80acf52 hicseq [hicrep.py] - corrected R_TOOLS path
       7cf9e94 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65bafda hicseq [hicrep.py] - Corrected typo
       990b8fa hicseq [hicrep.py] - corrected R_TOOLS path
       b86b171 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e5fdc3a hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff49511 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5da9af4 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8113e90 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       845ec62 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       09a2750 Completed developing hicrep and quasar analysis
       1e1ab4c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4f0cba1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6ba3f9d created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       be91a40 Completed hicrep analysis except the creation of the graph
       e2189f7 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fc54a3a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       09ec984 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       c58406d hicseq completed adding basic features of the hicrep analysis.
       36dd6fa hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       1861910 Added pairwise combination for samples
       ef297a3 started to edit the hicseq.py added new file for hicrep
       5d7c80f hicseq [base.ini] - updated mugqic tools version to 2.3.1
       868a711 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9685043 hicseq [quasar_qc - corrected module loading in matrix restructuring
       f5f88d6 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu Corrected conflicts Conflicts: 	bfx/quasar_qc.py 	pipelines/hicseq/hicseq.base.ini 	pipelines/hicseq/hicseq.py
       5a83860 hicseq [hicrep.py] - Corrected typo
       b9e6c3b hicseq [hicrep.py] - corrected R_TOOLS path
       6d4ce4f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b89f8e7 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       8beb590 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b78f92c [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       18523b8 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       0ff0ff7 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c6f98d6 Completed developing hicrep and quasar analysis
       0b10bc0 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       8b60a2b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6f50b29 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8faa8eb Completed hicrep analysis except the creation of the graph
       0bd1d0a hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b6ab90 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4a7bbca hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       74e0bcb hicseq completed adding basic features of the hicrep analysis.
       9cce006 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       183d5a9 Added pairwise combination for samples
       7f7301f started to edit the hicseq.py added new file for hicrep
       0f33a58 hicseq [hicrep.py] - Corrected typo
       c47a6ee hicseq [quasar_qc.py] - Corrected deletion by mistake
       8b33c6f hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       61be1c4 hicseq [hicrep.py] - corrected R_TOOLS path
       20c3b4a Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       edda3a1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       18506b4 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       75cbaaa [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2a0e1a8 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       15bbbb3 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       44735d2 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       eccd7c7 Completed developing hicrep and quasar analysis
       e9078d8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2ac9239 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       542bad0 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       418c515 Completed hicrep analysis except the creation of the graph
       f587c74 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a9980d hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2d87020 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8df5abc hicseq completed adding basic features of the hicrep analysis.
       f55e196 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       534d63f Added pairwise combination for samples
       756599e started to edit the hicseq.py added new file for hicrep
       b45e44d hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       cf5a2f4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       56985bb [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       b48157c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bdbc7de [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       d1c1ae9 Completed developing hicrep and quasar analysis
       e68d307 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c832c49 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1b63533 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2c192d1 Completed hicrep analysis except the creation of the graph
       e081f17 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       6f3560e hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       35efcde hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5d6fcc1 hicseq completed adding basic features of the hicrep analysis.
       ca0d163 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7216531 Added pairwise combination for samples
       396a06d started to edit the hicseq.py added new file for hicrep
       35412d6 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97ce474 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       38b339b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       53603e1 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db3da4b Completed developing hicrep and quasar analysis
       ba7bbb8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e6625b9 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       a797856 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       720a626 Completed hicrep analysis except the creation of the graph
       b5cbb3f hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bfbb896 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       1099ed0 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5f5f168 hicseq completed adding basic features of the hicrep analysis.
       f303163 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bbbcf3e Added pairwise combination for samples
       a5c29b5 started to edit the hicseq.py added new file for hicrep
       292a0b2 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       cd9d16a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       83f719b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       4c23bf1 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       5f42b4c Completed developing hicrep and quasar analysis
       849b382 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d9ed648 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1767a32 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       068daa9 Completed hicrep analysis except the creation of the graph
       210ad29 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       411d8e1 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0f31f81 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       79aff17 hicseq completed adding basic features of the hicrep analysis.
       b5075cb hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       57e8bae Added pairwise combination for samples
       4a57c1a started to edit the hicseq.py added new file for hicrep
       89006b6 Completed developing hicrep and quasar analysis
       87eb470 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       477281a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       10c5eb1 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       335f4c5 Completed hicrep analysis except the creation of the graph
       10daca1 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       4aeab54 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       509fef1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f096d1 hicseq completed adding basic features of the hicrep analysis.
       b623db5 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c605742 Added pairwise combination for samples
       9ad63fd started to edit the hicseq.py added new file for hicrep

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      16 commits

       a1ddafb Corrected hicseq.py file deletions
       e011fbe Corrected hicseq.py file deletions
       9c40c03 Corrected hicseq.py file deletions
       d8ded44 Corrected hicseq.py file deletions
       0a0903c Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       1fdc640 Corrected hicseq.py file deletions
       cbf6c69 Corrected hicseq.py file deletions
       96e5d3b Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       41a5bc2 Corrected hicseq.py file deletions
       3ee012f Corrected hicseq.py file deletions
       228a9e4 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       9211f3c Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       3f91825 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected merge conflict after rebase Conflicts: 	pipelines/hicseq/hicseq.py
       6291f89 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected module loading
       296b86b Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a757b8d Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      20 commits

       925b082 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3d8131d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ae8f296 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21e7665 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       938c27c hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       38516f3 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       08b1381 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4f835f9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       402a4d1 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       6b60d12 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       dd98639 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2d3a556 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cd42f4e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2631348 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       996ba98 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       d384dc9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       5ce6c34 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       bfd60bc hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       c743b44 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2305332 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudu.nawarathna@mail.mcgill.ca <pnawarat@abacus3.ferrier.genome.mcgill.ca>      1 commits

       e5fa345 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      20 commits

       1608e2b refer previous commit
       a749a35 refer previous commit
       ce31a37 refer previous commit
       4e8a71e refer previous commit
       2ba4345 refer previous commit
       5c06e43 refer previous commit
       bcac098 refer previous commit
       1f85535 refer previous commit
       0496744 refer previous commit
       f0eb790 refer previous commit
       afce945 refer previous commit
       fb33ca9 refer previous commit
       4598445 refer previous commit
       ba1f6b0 refer previous commit
       9227c2a refer previous commit
       3b5869f refer previous commit
       daccd07 refer previous commit
       34642e8 refer previous commit
       4bd83ce refer previous commit
       4ccc106 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       3c7dbdd Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a9fe3dd corrected md file
       fc8ac1a corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       1114ea8 Merged in hicseq_hicrep_pubudu (pull request #170)

3.2.0        Mon Jan 25 12:47:42 2021 -0500        371 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       ed04f33 testing end-of-line character
       957c11d Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       057a8d1 GenPipes - Genomes : updated R installation script to ease installation in dev, also corrected the hicup module called in install_genome.sh
       690f56e Version bump to 3.1.6-beta
       369cad4 Version bump to 3.1.5

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       1ce85b3 GenPipes : illumina_run_processing.py : correcting indentation
       f2d7e72 GenPipes - Tumor Pair pipeline : removed dev vawk module in tumor_pair.base.ini
       62ac679 GenPipes - DNASeq pipeline - removing mugqic_dev modules in gatk4.ini
       39ceaa2 GenPipes - Tumor Pair pipeline : removing mugqic_dev modules
       a84f5aa Merged in genpipes-3.1.5_release (pull request #166)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      4 commits

       42981f1 Merged in rnaseq_light_docs (pull request #183)
       a190faa rnaseq_light.py edited to adjust docstrings to address issue raised by Shaloo here : https://github.com/c3g/GenPipes/issues/63
       ceb0c8b Merged in JoseHectorGalvezLopez/nanoporebaseini-edited-online-with-bitbu-1596652406491 (pull request #179)
       35bc91b Edited the nanopore ini file to address the mugqic_tools error.

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       52f9015 modify cedar ini
       f38492f a bit ugly resolution from argparse overriding issue of the type argument

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      154 commits

       41671c6 Merged in mgi_stretenp (pull request #192)
       63a5468 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       a3cffe0 Changing qualimap bamqc section naming to a parameter in bfx
       6e33047 Changing qualimap bamqc section naming to a parameter in bfx
       ad67e66 Merged in mgi_stretenp (pull request #191)
       356bd20 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       dcbf891 igvtools ressources change
       f5026a8 Reducing ressources
       f15a4f8 Fixing awk
       c74a4c6 Changing default ressources
       fede180 igvtools ressources change
       58c298c Reducing ressources
       29ec475 Fixing awk
       cfca926 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       6618992 Changing default ressources
       9f663b0 Changing default ressources
       e74c221 Merged in mgi_stretenp (pull request #188)
       21f213f Fix
       cec6282 Fix
       907ab2a Fix interval_list checking
       08c1c54 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       9ba2711 Checking if interval_list is already on genome folder
       b6f2b91 Merged dev into mgi_stretenp
       247697b samtools bam2fq typo
       92d206f Adding kraken to beluga ini
       bb375b9 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       43d126e Switching to ivar 1.3 and NML dehosting filters
       7b86a54 Updating ivar trim usage
       d04bd81 Adding htslib to ivar module
       3ea8a2e Fixing kraken output and adding warning on rnaseq htseq-count
       9f4fde9 Switching to latest kraken2 db build
       cad6ed3 Forcing pigz to overwrite
       522905b Switching to latest kraken2 release
       de8e702 Fixing kraken
       e2249cc Fixing kraken bfx
       28dcdf3 Fixing kraken bfx
       63339cd Fixing kraken bfx
       0f47d1e Fixing kraken bfx
       338b71d Fix kraken
       3f5c373 Fix kraken
       d13c20d Addinf kraken analysis for metrics
       272131b Update sambamba sort ini
       6945104 Fix rename ln test
       940b62c Fix
       3ac107c Fixing tee
       624b2c2 Test
       c469a36 Fix
       9993fca Adding metrics on dehosted bam
       ff29ecf renaming cit covseq file
       01e6e6c cit ini for cit test
       67b78d6 Fixing hybrid genome path and output selection
       5d8d2f7 renaming cit covseq file
       15ba0d2 cit ini for cit test
       f6f9782 Switching to ivar 1.3 and NML dehosting filters
       4e16c72 Updating ivar trim usage
       4a07fed Adding htslib to ivar module
       05475b4 Fixing kraken output and adding warning on rnaseq htseq-count
       82ba084 Switching to latest kraken2 db build
       05e7d80 Forcing pigz to overwrite
       a11ec96 Switching to latest kraken2 release
       e43b7fb Fixing kraken
       2f5805b Fixing kraken bfx
       879ee55 Fixing kraken bfx
       19b050c Fixing kraken bfx
       2cb8818 Fixing kraken bfx
       eb614ee Fix kraken
       74227cc Fix kraken
       c797c2a Addinf kraken analysis for metrics
       52b1110 Update sambamba sort ini
       a06f588 Fix rename ln test
       398fcef Fix
       012b267 Fixing tee
       54b5fa5 Test
       33201b8 Fix
       516add1 Adding metrics on dehosted bam
       d5b321a Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       0b9b56a Reducing sambamba filtering default cpus
       77cddeb Fixing rename consensus
       916b091 renaming cit covseq file
       a550953 cit ini for cit test
       67a8174 Changing ini sections to use Illumina beds files as default
       fcfce99 Changing default seq_method
       84afdaa Cleaning and switching to cvmfs genome
       ed0fc19 Fix cit
       409ffcd Cit fix
       ccd50ff Cit test
       c999d20 Fix
       24fd48c Fix
       a11172c Fix
       03180fb Fix rename consensus symlink
       7533de7 Fixing rename consensus
       d33ab92 Fixing tsv
       b09b78e Fixing tsv
       2ab5bd4 Fixing tsv renaming
       5395489 Fixing output rename header + tsv for ncov-tools
       b0eb526 Fixing picard metrics
       13a5d31 Collecting picard metrics on raw AND filtered bam
       95dadb9 Fixing select output file
       8cfc8ea Fixing hybrid genome path and output selection
       35c8d49 Fixing merging step for 1 sample with multiple readsets
       210ebdd Update inis
       fe84c0b quast -> Quast
       13c60bf renaming cit covseq file
       9704c19 cit ini for cit test
       11fd809 Reducing sambamba filtering default cpus
       5843d7b Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       1fccae9 Changing ini sections to use Illumina beds files as default
       4c321ca Changing default seq_method
       ad84e07 Cleaning and switching to cvmfs genome
       5b66a1b Fix cit
       8ac2d78 Cit fix
       aa5609b Cit test
       a5a46a5 Fix
       360e7a1 Fix
       46bc396 Fix
       93fa452 Fix rename consensus symlink
       c48a298 Fixing rename consensus
       816eb57 Fixing tsv
       24e9734 Fixing tsv
       679da0c Fixing tsv renaming
       52b8e1b Fixing output rename header + tsv for ncov-tools
       8727d67 Fixing picard metrics
       c6a1a2d Collecting picard metrics on raw AND filtered bam
       860d44b Fixing select output file
       26be158 Fixing hybrid genome path and output selection
       facc8cb Fixing merging step for 1 sample with multiple readsets
       1442f6a Update inis
       126d53d quast -> Quast
       3a7504d renaming cit covseq file
       efac275 cit ini for cit test
       412da13 Changing ini sections to use Illumina beds files as default
       24ccf1f Changing default seq_method
       f4d97e9 Cleaning and switching to cvmfs genome
       3c8daf0 Fix cit
       e6e9738 Cit fix
       3ba059d Cit test
       7a29d01 Fix
       e6a8686 Fix
       e179616 Fix
       e528c79 Fix rename consensus symlink
       4de7076 Fixing rename consensus
       5430169 Fixing tsv
       2952b1a Fixing tsv
       4ca82a1 Fixing tsv renaming
       f7df72a Fixing output rename header + tsv for ncov-tools
       00b6052 Fixing picard metrics
       c009672 Collecting picard metrics on raw AND filtered bam
       fc6322f Fixing select output file
       a33b34c Fixing hybrid genome path and output selection
       75332a9 Fixing merging step for 1 sample with multiple readsets
       4a4b67e Update inis
       527a801 quast -> Quast
       6bc79c0 renaming cit covseq file
       d467187 cit ini for cit test

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       406e38b Merged in fix_monitor (pull request #190)
       e885338 Merged in monitor_bug (pull request #189)
       af56bf4 Merged in chunk_slurm_submit (pull request #185)
       bfe072a Merged in chunk_slurm_submit (pull request #182)
       b10eac6 Merged in rsync_in_chipseq (pull request #180)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      24 commits

       67c0af4 Merge branch 'release_3.2' of bitbucket.org:mugqic/genpipes into release_3.2
       90052a6 update readme for container install
       1294b96 Remove tumor_pair pipeline from release
       04c0c5a revert on indel aligner
       d75d458 cleanup ini for beluga
       cd12551 remove dbSNP for oxog ini block fix gatk 4 module name
       7d84e33 fix picard_collect_oxog_metrics dbSNP vcf
       97d6a47 remove HOME_DEV paths
       70d1310 make job2json more robust
       a05a243 remove mugqic_dev vawk
       7baea32 Fix README conflict
       3bda6d7 Fix out.out in monitor
       f0b99b8 rename to monitor.sh
       1ec8ba8 update control batch description
       bd24c9a rename deley_sbatch script
       c0f8012 fix script usage
       c579431 add export to sourced files
       f00a2d1 make sure already submited jobs id are sourced
       0a4444d if error
       2858608 remove debbug line
       5f6e768 add delay and chunking script
       88026d5 update for cit testing
       92edc02 replace -a by -r in rsync
       2551c78 fix master conflic with deleted readme file

  P-O Quirion <pioliqui@gmail.com>      2 commits

       0e48637 tweek monitor and chunk
       4f50d9e Fix autorestart and cleanup on failes interrupts

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      26 commits

       4f04bf4 covseq mpileup command fix
       acfe7b9 covseq qualimap fix, and high coverage ini fix
       9e0c8d7 add covseq - dnaseq consistency
       c6e5099 dnaseq high coverage trimmomatic to skewer
       8409952 tumor_pair fixes
       0c79ba9 import fixes
       5fb88cb bam.bai to .bai fix
       8894fc9 cit fixes to b38 - samtools and other b38 annotations
       87b0da8 conflict fixes to dnaseq and tumor pair after dev merge
       5a4dc7c fix mpileup dependency chain
       f8cda00 corrections to mpileup bcftools merge and tumor_pair dependencies
       b8ee480 bam.bai to .bai fix
       4c7b3f7 fixes to indentation in sambamba merge
       cb807fc fixes to sambamba merge
       3bc8f0e further gatk4 hc fixes
       0cdd2df update gatk4 hc arguments with gatk4 suite update
       8d7b9b8 fixes to gatk4 mark dup
       e7029c6 gatk4 mark dup walltime fix
       1d49954 variant recal fix
       9095ef1 variant recal dependency fix
       927011a cit fixes to b38 - samtools and other b38 annotations
       03ee662 minor dnaseq.py fixes
       4ffeb59 conflict fixes to dnaseq and tumor pair after dev merge
       cd935c7 fix mpileup dependency chain
       ef39fb7 corrections to mpileup bcftools merge and tumor_pair dependencies
       4d72c00 conforming deliverables to cit conventions

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      29 commits

       91f7655 exome specific fixes
       e6d4caf updates and fixes for cit
       72d072d remove briaree ini and update dnaseq base
       1547aa7 updates to beluga.ini and base.ini for dnaseq
       8e64f79 ini updates
       c925a61 gatk4 vsqr cit fix and baf plot for b38
       dd8c844 add cram to input selection
       8b2d9c3 multiqc fix
       55d8a64 argument fixes for picard imported functions in gatk4 and vqsr fixes
       12f55ac indel realignment consistency issues between dnaseq and tumor_pair
       1c6f1bc add mark dup metric file to multiqc
       4964f16 cit b38 fixes
       01f6a48 fix to bash.ln for cit
       2e1a585 fixes to multiqc
       d5cbe66 exome specific fixes
       e407f1d updates and fixes for cit
       a3bd732 Updates to bcftools/samtools for dnaseq_mpileup
       d962723 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       fc43094 remove briaree ini and update dnaseq base
       72e8439 adding 1 job processing to specific steps for cit.ini
       f4911b3 fix of dev genome reference
       2447075 issues between dnaseq.base.ini and dnaseq.beluga.ini
       cc7c806 updates to modules for beluga
       d5f0301 updating beluga ini
       6775896 added tumor pair beluga ini
       16bcac2 updates to beluga.ini and base.ini for dnaseq
       8e8f22e updates to beluga.ini
       f8ea263 ini updates
       de84a17 module updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      15 commits

       369bd8a picard2 and high coverage fixes
       84885da updates to b38 variant recal files
       3db1c5a fixes to tumor_pair on beluga
       a40514c cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       e04c3f8 fixes to fixmate input file name
       9fe39b0 picard2 and high coverage fixes
       30cf939 tumor_pair beluga ini fix
       c77c821 tumor_pair qualimap part 2
       5c039c6 qualimap tumor_pair fix
       b980f2a fixes to vsqr gatk4
       ed99603 gatk4 fixes to callable loci and DoC
       6162e8e sym link dnaseq.base into tumor pair
       950a285 updates to b38 variant recal files
       f6c42cb fixes to tumor_pair on beluga
       d8eaf58 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      24 commits

       2f926d6 fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       800ecb2 cit fixes after rebasing
       a37298c dependency fixes
       ba400c9 updates to GRCh38 annotation file, module updates, cit fixes
       01f0249 updates to cit and fixes to one job mpileup steps
       3e213bd updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8b745c1 major fixes to deliverables and completion of beluga test
       3078a9a fix modules dev to cvmfs
       2345caf fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       3c28403 fixes to merge_filter_bcf
       662e007 gatk4 bsqr fixes and mpileup cat
       0f3aefc mpileup protocol fix
       aefbb37 cit fixes after rebasing
       6f47ad3 update to gatk4 mutect2 filtering procedures
       419f376 create cit for gatk4 due to deviation from argument usage
       04fcaac cit fixes to gatk4 + sym links for recalibration
       802cc65 dependency fixes
       7c3b8d3 updates to GRCh38 annotation file, module updates, cit fixes
       2882826 fixes to deliverable and b38 ini
       81204f6 updates to cit and fixes to one job mpileup steps
       dfa3288 updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8bc536d fixes to metasv annotations
       556e48b major fixes to deliverables and completion of beluga test
       973087a fix modules dev to cvmfs

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      1 commits

       b5c1486 final PR fixes

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      7 commits

       9ac20da code cleaning and fixes to exome interval list
       f55a738 fixes to cedar ini
       deca50f fixes to symlinks for paired indel realignment
       31525ee code cleaning and fixes to exome interval list
       fa7328f fixes to cedar ini
       2d3d367 fixes to symlinks for paired indel realignment
       5475ca4 cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      5 commits

       7fbd4f9 fixes to metasv, adding metasv germline
       5a6a9d7 cedar fixes and GRCh38 fixes
       b42ed16 cedar dnaseq updates and svaba germline added
       2e389e9 cedar germline sv updates
       67a2e47 sequence dictionary and module updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      35 commits

       7d71caf updates to metasv - somatic
       eada99c updates to SV germline and reference genome tweaks
       0aa9d42 somatic sv fixes: lumpy and svaba annotations
       0f27790 fix to germline SV: breakseq2
       df8096a merge fixes
       a96622a json and folder updates
       ee184e2 fixes to sCNAphase
       fa8f99f merging snv and sv, adding protocols
       6cae335 Fixes to indel realigner
       1609905 Updates and debug
       04229ba Add set somatic and actionable mutations
       937e543 added multiqc and other tweaks
       2be1a97 add metrics step for metasv
       a9697c7 updates to metasv - somatic
       f622b92 Fixes and updates to reference files
       909c2ef remove testing steps
       0cb7a4b updates to SV germline and reference genome tweaks
       de92d57 somatic sv fixes: lumpy and svaba annotations
       573c45d fix to germline SV: breakseq2
       f2fc9d6 merge fixes
       77bc119 GATK4 fixes - bam indexing and markDupSpark
       21051bc bcftools fixes for tumor pair
       44f3d04 fingerprint and bug fixes
       c9c2049 dnaseq - vcftools qc addition: --missing_indv and --depth
       9428baa select input for variant caller and fixes to one job calling
       1c78b28 json and folder updates
       3731f20 fixes to sCNAphase
       402bccb Added json info
       036606a Bug fixes prior to json additions
       bd6e6cb merging snv and sv, adding protocols
       d441d47 Fixes to indel realigner
       c07293c Add deliverables module
       f674ec9 Updates and debug
       3161cf6 Add set somatic and actionable mutations
       e8d90db added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      14 commits

       23f5be7 remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       614ee8f fixes to metasv for tumor pair
       9f1b500 Bug fixes and modification derived from initial PROFYLE benchmarking
       245df59 remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       079fd7e fixes to metasv for tumor pair
       48bb3c5 Single job bug fixes
       b77cfaa manta I/O fix and other bug fixes
       d14ba6d config updates and b38DH added
       6b52879 dnaseq germline SV updates
       b82f0f6 gatk4 updates and bug fixes
       14310b4 gatk4 updates and bug fixes
       77126e8 Fix line break types
       9fe9453 Json related bug fixes
       dcf0b8a Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       04b5490 fixes to breakseq2 and metasv
       b0626ec cit-based fix - verifyBAMid
       b897ef1 dnaseq qc additions: NGScheckmate and peddy
       0cbfb0d fixes to breakseq2 and metasv
       612c446 cit-based fix - verifyBAMid
       52023e1 cit-based fixes to NGScheckmate
       4dec63d dnaseq - multiqc fix
       29e26f7 dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      4 commits

       f3c972e Merged in rebasing_tp (pull request #186)
       f6e74bd Debugging and Guillimin specfic fixes
       4810e9a Debugging and Guillimin specfic fixes
       0d16e4d updates to config

  Rom Grk <romgrk.cc@gmail.com>      2 commits

       1d5912c Merged in fix-watch-portal-folder (pull request #187)
       c22599c fix: use of undeclared variables

covid_1.0        Mon Aug 3 19:56:47 2020 -0400        781 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      19 commits

       a24b8db GenPipes Covid Release - Renamed README.md as INFO.md and added GenPipes-like README.md for CovSeq pipeline
       73a694f Version bump to Covid Release 1.0
       eda88d0 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       3555e1b Version bump to Covid Release 1.0
       a2f2178 GenPipes - Covid Release : adding resources
       9e07fc4 GenPipes - HiC-Seq pipeline : corrected fastq_readName_edit input path
       dd58212 GenPipes - DNA-Seq pipeline : correcting symlink creation in sambamba_merge_sam_file
       7e00f1d Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1799cea updated version of MUGQIC_TOOLS in installation script
       fb8ad6e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7a9b872 GenPipes bug correction - corrected ln() function in bfx/bash_cmd.py
       62402ad GenPipes Update - debugging use portal_output_dir variable : check for both undef and empty value
       2830d23 GenPipes Genome - added Sus_scrofa.sh (Pig genome) installation script
       81d580a GenPipes Soft - added kent.sh installation script
       756d9a0 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       076e200 GenPipes Update - DNASeq - updated sym_link for better handling of path reconstruction
       19cc2c0 GenPipes Update : fixing path of config files when passed to job2json script
       1d204f9 DNASeq - removed use of 'os.path.abspath' in call of 'ln()'
       399aac0 DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable

  Édouard Henrion <henrione@beluga3.int.ets1.calculquebec.ca>      4 commits

       4640569 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       8c46541 GenPipes - ChIPSeq update : a bit of code cleaning and simplifying
       b03ee2f Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       72c7538 GENOME INSTALLATION - updated install genome.sh : added bismark genome preparation + refined genome digest command

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      2 commits

       0dc068b Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       c744938 Merged in eh_quick_fixes (pull request #144)

  Édouard Henrion <henrione@ip18.m>      6 commits

       b190926 GenPIpes Update - corrected one problematic sym_link call...
       ed444a7 GenPipes Update - corrected pipeline behavior regarding PORTAL_OUTPUT_DIR environment variable : if te variable is empty or not set, then no JSON file at all will be generated
       c540751 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6331eed C3G Software - added demuxlet.sh installation script
       d63b356 Genome Update - install_genome.sh
       a899066 Genome Update - some updates to Homo_sapiens.GRCh38.sh

  ehenrion <edouard.henrion@mcgill.ca>      15 commits

       7d9f811 BFX- Software - iVar installation script update with latest version
       3c3a80c common.py : genpipes replacing mugqic_pipeline....
       8658f7b GenPipes RNA-SEq - calling DESeq2 instead of DESeq for differential_expression
       6a0c825 Coorected typo in README.md
       d59cb59 GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       af972e1 GenPipes - pipelines/dnaseq.py  : corrected prefix generation in SymLinkFastq step
       8ebbd14 GenPipes - pipelines/common.py corrected outputs name generation patterns for SamToFastq & Trimmomatic steps
       65fe80a Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918129493 (pull request #142)
       97afb33 GenPipes - dnaseq.py : bug correction - typo removed
       4f9a8dd GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       fc67a7b GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       9fbd085 GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       16a5635 GenPipes - nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       853c806 updated methylseq.base.ini, useless comments removed
       be965bd updated methylseq.base.ini, useless comments removed

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      27 commits

       d485728 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       7cb4024 final edits to the nanopore pipeline
       f1db8a9 Added support for gzipped fastq
       8da4a1c Added the nanopore CIT ini file
       86988df Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       71f1658 Final corrections before merge to dev
       6994d69 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       83a6f5f More corrections to the INI files for other servers to allow proper running with SLURM
       caf6cb2 Corrected an error in the mp2b ini file
       e2637c0 Made corrections to the nanopore ini files for all servers that were causing the pipeline to break when running with SLURM
       a10fe05 Merge conflict resolutions
       e08cd6d Corrected problem with nanopore readset parsing that caused problem with the paths
       5d884c1 Added full documentation for Nanopore pipeline, as well as the Graham config file.
       e2edb70 Included commands necessary to add readgroup tag to alignments in minimap2
       278feb7 Final corrections before testing on other servers
       cc0f765 Corrected merge error related to .gitignore
       cb83788 Fixed bug caused by missing module import in nanopore.py
       649e04b Added minimap2 script that was missing from previous commit
       0f2a89d First working version of the nanopore pipeline
       7acad99 Added full documentation for Nanopore pipeline, as well as the Graham config file.
       7c528bd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       56699ff Included commands necessary to add readgroup tag to alignments in minimap2
       3d53188 Final corrections before testing on other servers
       558f3e9 More bug corrections for nanopore pipeline after initial testing. Switched to only one protocol, with an optional first step (guppy)
       5107ed7 Fixed bug caused by missing module import in nanopore.py
       653c4e0 Added minimap2 script that was missing from previous commit
       a0fe5eb First working version of the nanopore pipeline

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      8 commits

       50bc9f3 Merged in Jose-Hector-Galvez/rnaseq_lightbaseini-edited-online-with-b-1580485656011 (pull request #171)
       e2c91c8 Added module_perl to the rnaseq_light ini file.
       007835b Merged in Jose-Hector-Galvez/jobpy-edited-online-with-bitbucket-1576266027164 (pull request #154)
       c500f6b Suggestion, add `module purge` to all jobs that load modules, to avoid conflicts between modules.
       4c35a7c Merged in Jose-Hector-Galvez/gatk4py-edited-online-with-bitbucket-1575403540009 (pull request #135)
       89e70fe gatk4.py edited to correct for inconsistencies with configuration parameters within functions.
       eefc95a Merged in Jose-Hector-Galvez/found-a-bug-in-the-schedulerpy-script-i--1575322351581 (pull request #127)
       afae6ee Found a bug in the scheduler.py script. I am adding a line to correct it.

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      1 commits

       3ec3f1c Added CoVSeQ ini files for Graham and Cedar. Corrected a few errors on the Beluga ini

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       10e762b Merged in mgi_stretenp (pull request #177)
       3960717 Merged in nanopore_jhg (pull request #173)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      644 commits

       41514c0 Merged in mgi_stretenp (pull request #175)
       d6ae129 Adding WARN for not changing R ver in rnaseq denovo
       f51f693 Fixing rnaseq denovo R issue by creating deseq and a deseq2
       4a1a5bf Fixing back rnaseq denovo
       56fb3fb Fixing rnaseq denovo
       1d29eb6 Fixing rnaseq denovo assembly R versions
       50fb470 Fixing trinity report on cit
       1ae8b17 deliverables fix
       c51f7e1 Switching to 10% as minor variants threshold
       dca8fea Removing kraken module & Changing R version for rnaseq
       c49a073 Including sambamba changes in ini
       29b5aec Including sambamba changes into other ini
       c4a7ef3 Fixing gatk_indel_realigner if 1 job
       8fd747f Including sambamba modifs in ini
       d826d81 sambamba merge realigned high coverage
       f135baf Switching to sambamba merge sam
       c55297c Including bwa sambamba into high coverage
       466906d Inluding a with sambamba bam.bai to picard mark duplicates
       858c1f5 Merge branch 'dnaseq_picard_fix' into mgi_stretenp
       a3354a5 Typo
       dd39ded Fixing for cit run
       1b4d014 Merge branch 'dev' into mgi_stretenp
       d9566ff Using Illumina as default
       34cd422 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       2502aba Renaming outputs and default genome
       8ccb242 Fixing consensus renaming
       075ab5a Fixing consensus renaming
       9dec668 Renaming + flaging consensus seq
       621449a Fixing quast options
       da26384 Fixing qualimap output
       cb98c5e Fixing outputs quast and snpeff
       6d22060 Fixing mkdir
       d15425c Fixing typo
       921d565 Fixing typo
       53c68fd Changing to samtools bam2fq and fixing options ini
       4f5c5ce Fixing quast
       f7042ec Fixing latest commit
       04f2916 Adding intermediate bam file on host removal step
       d931a89 Forcing pigz to overwritte gz
       b0711a8 fixing host removal
       4cfa142 Adding pigz within bam2fq to be able to skip step
       7916935 quast fix
       4788bae quast + snpeff + host removal fix
       2e3bc59 Fixing cutadapt input files
       75d6382 Fixing type
       308ace1 Removing indexing after name sorting
       f76fdad Fixing path creation at host removal step
       d5aa6f0 Fixing snpeff
       9202478 Removing print files
       afa134e Fixing input choosing
       979eafe Fixing choosing inout files mapping
       9b6ae43 Fixing host removal
       f7ba725 Fixing host removal
       f0ac000 Fixing quast step
       8ce9bf2 Fixing quast
       c07ffae Fixing quast step
       5415ab7 Fixing param requirements
       9218a8e Fixing typo
       52999ae Fixing pigz
       fec2cea Not using kraken anymore
       3363c61 Switching steps order
       7e51be5 Adding 3 steps
       863d61b Fixing picard multiple metrics raw
       4d8cc43 Insert Size metrics on raw bam
       f2e8b42 sambamba flagstat fix
       0a8c640 Renaming snpeff step/job
       d473404 Fixing flagstat
       8163c9a Fixing flagstat
       8461e19 flagstat on all the bams
       a380f68 Zipping output of snpeff
       28969ba Fixing snpeff
       a899bdd Fixing cutadapt
       2f86db5 Flagstat on raw bam
       f2b8a9f Fix
       f4fba5a Fixing renaming step
       b0fcd1c Switching to fgbio trim primer
       e9c450f Updating metrics
       6582bc9 Addition of consensus step
       319e50d MGI init commit
       6df0cca Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       350131a DNA-Seq - Fix update
       5cce818 DNA-Seq - Fix update
       838cc92 DNA-Seq - Fix update
       990fc56 DNA-Seq - Fix update
       7e7b810 DNA-Seq - Fix update
       0a8d7ee DNA-Seq - Fix update
       fe22929 DNA-Seq - Fixmate with samtools and sorting with sambamba
       3c84dd1 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       e0fa9ff DNA-Seq - Fix update including input_file_revision_eh changes
       d03669f DNA-Seq - Fix update
       5f22222 DNA-Seq - Fix update
       18dd106 DNA-Seq - Fix update
       08c2864 DNA-Seq - Fix update
       117181c DNA-Seq - Fixmate with samtools and sorting with sambamba
       813c02e DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c193068 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       0e8d422 DNA-Seq - Fix update including input_file_revision_eh changes
       80f732a Zipping and indexing vcfs
       6b1179b Adding parameter to ivar consensus
       a463f5b Changing caller
       75a071a Addition of consensus step
       5639985 MGI init commit
       9811842 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       34dac0c cutadapt bfx
       9ccadd3 DNA-Seq - Fix update
       8733abc DNA-Seq - Fix update
       e007f9a DNA-Seq - Fix update
       93e342c DNA-Seq - Fix update
       0160ab5 DNA-Seq - Fix update
       567bd8a DNA-Seq - Fixmate with samtools and sorting with sambamba
       9e4119b DNA-Seq - Fix update
       7ed58f8 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       aefdd2c DNA-Seq - Fix update including input_file_revision_eh changes
       8a94a76 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       36fbd80 DNA-Seq - Fix update
       cb43b48 DNA-Seq - Fix update
       069de1c DNA-Seq - Fix update
       2775e6f DNA-Seq - Fix update
       348c2c6 DNA-Seq - Fix update
       809b43d DNA-Seq - Fix update
       08c4183 DNA-Seq - Fix update
       7314301 DNA-Seq - Fixmate with samtools and sorting with sambamba
       b4728a8 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       d22a3b1 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7a608be DNA-Seq - Fix update including input_file_revision_eh changes
       96dcc43 MGI-Seq - First commit
       c6e80af Adding snpEff step for mgi pipeline
       0e998c0 Fixing sambamba flagstat
       79e8114 Fixing sambamba module
       3bd9321 Adding module sambamba to ini
       812e317 Fixing sambamba indexing
       b0cbd93 Switching from picard to sambamba withing methylseq
       3f5d5a3 Fix
       271b81d fix
       6b0c09f Fixing input output files
       79e49f4 Fixing filtering bam
       c4dee0f Fixing renaming step
       f131124 Fix renaming step
       c765138 fix
       c04a496 Fix
       eea7c7b Fixing filtering
       c0bc132 Adding filtering step
       7d34c4e Switch to fgbio
       f38cd7d ivar triming switch
       15b7438 Switching to fgbio trim primer
       cec53cb Updating metrics
       b8de4e3 MGI init commit
       11b01b0 cutadapt bfx
       e3858aa DNA-Seq - Fix update
       1eaa3dc DNA-Seq - Fix update
       d23743f DNA-Seq - Fix update
       bd59ddb DNA-Seq - Fix update
       4eb3268 DNA-Seq - Fix update
       9d8e58b DNA-Seq - Fix update
       c93e233 DNA-Seq - Fix update
       2934b70 DNA-Seq - Fixmate with samtools and sorting with sambamba
       ff53bb3 DNA-Seq - Fix update
       5d3fb26 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       bb731b1 DNA-Seq - Fix update
       1b68132 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       3a6e72d DNA-Seq - Fix update including input_file_revision_eh changes
       01551e5 DNA-Seq - Fix update
       dcaedc8 DNA-Seq - Fix update
       a01251a DNA-Seq - Fix update
       56a50d0 DNA-Seq - Fix update
       32e070f DNA-Seq - Fix update
       e4778b0 DNA-Seq - Fix update
       9456045 DNA-Seq - Fix update
       083f6a6 DNA-Seq - Fixmate with samtools and sorting with sambamba
       67603dd DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       a59d994 DNA-Seq - Fix update including input_file_revision_eh changes
       c203912 MGI-Seq - First commit
       c5a1baa Switch to fgbio
       b4ca5f7 Fixing alignment
       f3378e1 ivar triming switch
       8a0045d Switching to fgbio trim primer
       d479b85 Fixing ivar trim
       80c578a Filtering reads
       30f51c7 Updating metrics
       ea1149b fgbio
       22daf31 Fixing ivar trim
       6462908 Using ivar trim instead of fgbio
       fc6de5f Adding ivar primer trimming
       7049f31 Default bwa parameters to include pairs
       0e834bc Zipping and indexing vcfs
       0935479 Adding parameter to ivar consensus
       d8751cc Changing caller
       1e6dd74 ivar bfx
       3399cbb Addition of consensus step
       f46ccb5 ivar module
       52f8a2b MGI init commit
       3bd0932 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       2fbd1e3 Adding trim_primers to fgbio bfx
       327b963 cutadapt bfx
       7652e5b DNA-Seq - Fix update
       0ab9060 DNA-Seq - Fix update
       66b240b DNA-Seq - Fix update
       472c6b0 DNA-Seq - Fix update
       e1e7318 DNA-Seq - Fix update
       6954eaf DNA-Seq - Fix update
       db7bab2 DNA-Seq - Fix update
       2a0bb00 DNA-Seq - Fix update
       60e7d9a DNA-Seq - Fixmate with samtools and sorting with sambamba
       7be45c1 DNA-Seq - Fix update
       1dddc4f DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       27ac18b DNA-Seq - Fix update including input_file_revision_eh changes
       748e668 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       0b53786 DNA-Seq - Fix update
       6218cdd DNA-Seq - Fix update
       56547aa DNA-Seq - Fix update
       05baa00 DNA-Seq - Fix update
       ff5a871 DNA-Seq - Fix update
       591339e DNA-Seq - Fix update
       de100ee DNA-Seq - Fix update
       5e3e598 DNA-Seq - Fix update
       584181f DNA-Seq - Fix update
       02439f0 DNA-Seq - Fixmate with samtools and sorting with sambamba
       90e9ee9 DNA-Seq - Fix update
       ac56691 DNA-Seq - Fix update
       f73bce0 DNA-Seq - Fix update
       22e030d DNA-Seq - Fix update
       6928926 DNA-Seq - Fix update
       0facc34 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c6bf187 DNA-Seq - Fix update
       111cbcb DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ef0a25e DNA-Seq - Fix update ini files
       2eed96b DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       417ff67 DNA-Seq - Fix update bai output
       114d004 DNA-Seq - Fix update indexing
       fc9901c DNA-Seq - Fix update including input_file_revision_eh changes
       0016ba2 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       77e5041 MGI-Seq - First commit
       9447bed Renaming outputs and default genome
       6f45655 Fixing consensus renaming
       87a6264 Fixing consensus renaming
       cb95534 Renaming + flaging consensus seq
       6c7058d Fixing quast options
       1ec14ff Fixing qualimap output
       90eb6b0 Fixing outputs quast and snpeff
       3dbf0a1 Fixing mkdir
       8cd2ba3 Fixing typo
       200cb4b Fixing typo
       2b01f81 Changing to samtools bam2fq and fixing options ini
       d5dc5b7 Fixing quast
       f32dd81 Fixing latest commit
       d58a7b2 Adding intermediate bam file on host removal step
       dbd486b Forcing pigz to overwritte gz
       95505ba fixing host removal
       a119100 Adding pigz within bam2fq to be able to skip step
       de7633c quast fix
       d0383f5 quast + snpeff + host removal fix
       f165e56 Fixing cutadapt input files
       2dfd2e6 Fixing type
       6e67b6d Removing indexing after name sorting
       4b0740e Fixing path creation at host removal step
       7020fe0 Fixing snpeff
       5e6b7a9 Removing print files
       031d9df Fixing input choosing
       256cdc7 Fixing choosing inout files mapping
       cd7e52e Fixing host removal
       c11d220 Fixing host removal
       a3d58c9 Fixing quast step
       8abe6c0 Fixing quast
       5a68e07 Fixing quast step
       b141fa8 Fixing param requirements
       a7318fd Fixing typo
       aa32b62 Fixing pigz
       617f817 Not using kraken anymore
       df8f32c Switching steps order
       10c3586 Adding 3 steps
       eb7d612 Fixing picard multiple metrics raw
       8d852a7 Insert Size metrics on raw bam
       67dcc5f sambamba flagstat fix
       08fb926 Renaming snpeff step/job
       068d072 Fixing flagstat
       770c233 Fixing flagstat
       4ff0fb6 flagstat on all the bams
       1e3ef28 Zipping output of snpeff
       5c219e6 Fixing snpeff
       c0faa37 Fixing cutadapt
       effce8a Flagstat on raw bam
       8c7dd5a Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       221dea2 Adding snpEff step for mgi pipeline
       86a0acf Fixing sambamba flagstat
       dd8163f Fixing sambamba module
       609e6ca Adding module sambamba to ini
       71316e9 Fixing sambamba indexing
       9040722 Switching from picard to sambamba withing methylseq
       9d3bdfd Fix
       e5a7f40 fix
       59aeaff Fixing input output files
       da3b1cd Fixing filtering bam
       cb775da Fixing renaming step
       f78909c Fix renaming step
       49dba7d fix
       50ac6ee Fix
       860591e Fixing filtering
       cc0f319 Adding filtering step
       4d36888 Switch to fgbio
       6066f52 ivar triming switch
       5790ea2 Switching to fgbio trim primer
       45c2449 Updating metrics
       ffc8fd1 MGI init commit
       13d4f10 cutadapt bfx
       2a61965 DNA-Seq - Fix update
       7c8d9e3 DNA-Seq - Fix update
       aad95a2 DNA-Seq - Fix update
       db01cdb DNA-Seq - Fix update
       c8fd1f4 DNA-Seq - Fix update
       19b7a74 DNA-Seq - Fix update
       386a604 DNA-Seq - Fix update
       6063ddd DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e1f0b5 DNA-Seq - Fix update
       30121a7 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       4b9d8e1 DNA-Seq - Fix update
       a83ee80 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       033dd3e DNA-Seq - Fix update including input_file_revision_eh changes
       a9102d3 DNA-Seq - Fix update
       a4ee47f DNA-Seq - Fix update
       18d832e DNA-Seq - Fix update
       6623f8c DNA-Seq - Fix update
       1d002d1 DNA-Seq - Fix update
       fe5c405 DNA-Seq - Fix update
       dbd8034 DNA-Seq - Fix update
       74d2bdd DNA-Seq - Fixmate with samtools and sorting with sambamba
       1678639 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       5f66227 DNA-Seq - Fix update including input_file_revision_eh changes
       562180b MGI-Seq - First commit
       a16695b Switch to fgbio
       1d92909 Fixing alignment
       01f5df7 ivar triming switch
       e49e632 Switching to fgbio trim primer
       057e2e8 Fixing ivar trim
       d9638b5 Filtering reads
       c3ab6f1 Updating metrics
       dc9a43e fgbio
       579e0dd Fixing ivar trim
       603bb4a Using ivar trim instead of fgbio
       26a2da8 Adding ivar primer trimming
       f25c75d Default bwa parameters to include pairs
       6a8c301 Zipping and indexing vcfs
       cf6c19f Adding parameter to ivar consensus
       c9231e5 Changing caller
       b9bfed0 ivar bfx
       70cc54a Addition of consensus step
       5a990e2 ivar module
       cd17c0d MGI init commit
       5bef158 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       f7a155d Adding trim_primers to fgbio bfx
       8177037 cutadapt bfx
       a702d27 DNA-Seq - Fix update
       f6fe0e2 DNA-Seq - Fix update
       4a42439 DNA-Seq - Fix update
       a880e3b DNA-Seq - Fix update
       2cd5d01 DNA-Seq - Fix update
       9f7e048 DNA-Seq - Fix update
       75b18ca DNA-Seq - Fix update
       9a0adea DNA-Seq - Fix update
       bb259fe DNA-Seq - Fixmate with samtools and sorting with sambamba
       0c96937 DNA-Seq - Fix update
       75e589b DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       dd3999c DNA-Seq - Fix update including input_file_revision_eh changes
       2a2a43c DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       cb09a7b DNA-Seq - Fix update
       31653bd DNA-Seq - Fix update
       513818c DNA-Seq - Fix update
       b202a5b DNA-Seq - Fix update
       c87c2cd DNA-Seq - Fix update
       a4926cb DNA-Seq - Fix update
       c06ee52 DNA-Seq - Fix update
       6f894d2 DNA-Seq - Fix update
       3d5daeb DNA-Seq - Fix update
       ed9e7fe DNA-Seq - Fixmate with samtools and sorting with sambamba
       4215929 DNA-Seq - Fix update
       3dca79c DNA-Seq - Fix update
       7bf5141 DNA-Seq - Fix update
       ff24ef9 DNA-Seq - Fix update
       f535202 DNA-Seq - Fix update
       a6b428c DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       7e35168 DNA-Seq - Fix update
       302e6c9 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       03fc8c2 DNA-Seq - Fix update ini files
       0cd4885 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       ade8279 DNA-Seq - Fix update bai output
       8ee3267 DNA-Seq - Fix update indexing
       4f25a98 DNA-Seq - Fix update including input_file_revision_eh changes
       314107e DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       6d82d3f MGI-Seq - First commit
       2cbfa5c Adding snpEff step for mgi pipeline
       8f5d839 Fixing sambamba flagstat
       4415a7a Fixing sambamba module
       644bc8d Adding module sambamba to ini
       38cd45a Fixing sambamba indexing
       08dc201 Switching from picard to sambamba withing methylseq
       e89813f Fix
       9ef7ecb fix
       0d6651a Fixing input output files
       192bb65 Fixing filtering bam
       030b759 Fixing renaming step
       1348f00 Fix renaming step
       bc42658 fix
       49a16fc Fix
       c69973d Fixing filtering
       76219a8 Adding filtering step
       7ed54e8 Switch to fgbio
       b6cfddf Fixing alignment
       4796cff ivar triming switch
       3ec5f7f Switching to fgbio trim primer
       3a84e9e Fixing ivar trim
       b673f47 Filtering reads
       9fdafbd Updating metrics
       3f63d5c fgbio
       56bec97 Fixing ivar trim
       30475ca Using ivar trim instead of fgbio
       d5b4054 Adding ivar primer trimming
       925645e Default bwa parameters to include pairs
       314a869 Zipping and indexing vcfs
       191bfb3 Adding parameter to ivar consensus
       db7c70a Changing caller
       d9e302b ivar bfx
       09d9f6b Addition of consensus step
       d9c3f32 ivar module
       f6a8e17 MGI init commit
       49c7f1d Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       eae7037 Adding trim_primers to fgbio bfx
       bfc47be cutadapt bfx
       ffcfe26 DNA-Seq - Fix update
       0f17426 DNA-Seq - Fix update
       34b0965 DNA-Seq - Fix update
       b7323de DNA-Seq - Fix update
       e0909f1 DNA-Seq - Fix update
       48e203e DNA-Seq - Fix update
       0cd7e5e DNA-Seq - Fix update
       b3d7401 DNA-Seq - Fix update
       6ff06b0 DNA-Seq - Fixmate with samtools and sorting with sambamba
       b000f97 DNA-Seq - Fix update
       52deaf3 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       2acf3c0 DNA-Seq - Fix update including input_file_revision_eh changes
       0fe44df DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       3772b9b DNA-Seq - Fix update
       479e9ed DNA-Seq - Fix update
       5f25120 DNA-Seq - Fix update
       aa1687c DNA-Seq - Fix update
       8632426 DNA-Seq - Fix update
       c5d384b DNA-Seq - Fix update
       746d2f0 DNA-Seq - Fix update
       eef3b27 DNA-Seq - Fix update
       f1d1b55 DNA-Seq - Fix update
       d2d4430 DNA-Seq - Fixmate with samtools and sorting with sambamba
       5683297 DNA-Seq - Fix update
       a7e450c DNA-Seq - Fix update
       5cb9574 DNA-Seq - Fix update
       e888e09 DNA-Seq - Fix update
       d366f67 DNA-Seq - Fix update
       ba6aa73 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       9dc3bf6 DNA-Seq - Fix update
       c01dcf5 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ca81d44 DNA-Seq - Fix update ini files
       78519f2 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       0c922cc DNA-Seq - Fix update bai output
       6117bc9 DNA-Seq - Fix update indexing
       5e9b268 DNA-Seq - Fix update including input_file_revision_eh changes
       519beea DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       43bb771 MGI-Seq - First commit
       f2cd11a Switch to fgbio
       c3aa00d Fixing alignment
       8d157fd ivar triming switch
       42a5e2d Switching to fgbio trim primer
       6030fa7 Fixing ivar trim
       6bb7e08 Filtering reads
       330e345 Updating metrics
       b1153cd fgbio
       4450c70 Fixing ivar trim
       a68648b Using ivar trim instead of fgbio
       8a2bdf3 Adding ivar primer trimming
       8856367 Default bwa parameters to include pairs
       4e96864 Zipping and indexing vcfs
       6e9529f Adding parameter to ivar consensus
       4b826ef Changing caller
       b8bc79a ivar bfx
       dab5579 Addition of consensus step
       c91d067 ivar module
       ec2a176 MGI init commit
       260059a Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       7760a79 Adding trim_primers to fgbio bfx
       3180972 Fixing haplotype_caller
       66d8cdd Fixing haplotype_caller
       3bdb0aa cutadapt bfx
       2d274eb Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       f083b66 MGI-Seq - First commit
       538907e Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into dnaseq_picard_fix
       b0197fe DNA-Seq - Fix update
       54b4741 DNA-Seq - Fix update
       452f869 DNA-Seq - Fix update
       d0cb631 DNA-Seq - Fix update
       6e3580e DNA-Seq - Fix update
       222b79b DNA-Seq - Fix update
       5cf6269 DNA-Seq - Fix update
       de8d901 DNA-Seq - Fix update
       d87cd96 DNA-Seq - Fix update
       b3b9b1b DNA-Seq - Fixmate with samtools and sorting with sambamba
       5e55bfc DNA-Seq - Fix update
       9f20fe7 DNA-Seq - Fix update
       123b730 DNA-Seq - Fix update
       5ae600a DNA-Seq - Fix update
       8119523 DNA-Seq - Fix update
       d6007d2 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       1ff8cac DNA-Seq - Fix update
       58d662a DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       83ba60d DNA-Seq - Fix update ini files
       f7bef7d DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       10d50aa DNA-Seq - Fix update bai output
       1319c26 DNA-Seq - Fix update indexing
       fe4e7b9 DNA-Seq - Fix update including input_file_revision_eh changes
       80def2a DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       894566c DNA-Seq - Fix update
       1413d45 DNA-Seq - Fix update
       ab945dd DNA-Seq - Fix update
       97186dd DNA-Seq - Fix update
       fc158a4 DNA-Seq - Fix update
       86b2e80 DNA-Seq - Fix update
       d2ecd4b DNA-Seq - Fix update
       98928a5 DNA-Seq - Fix update
       bdee1f8 DNA-Seq - Fix update
       8afb500 DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e8699b DNA-Seq - Fix update
       1e47ab1 DNA-Seq - Fix update
       bea5e94 DNA-Seq - Fix update
       eae982e DNA-Seq - Fix update
       15d745e DNA-Seq - Fix update
       0553a5d DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       ef79ced DNA-Seq - Fix update
       afaa484 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7dfa730 DNA-Seq - Fix update ini files
       3d22e67 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       6dd6e0c DNA-Seq - Fix update bai output
       54a5fc0 DNA-Seq - Fix update indexing
       20cee1f DNA-Seq - Fix update including input_file_revision_eh changes
       a1579f4 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       c5e9966 Merged in ihec_metrics (pull request #126)
       7bceefb Merge branch 'ihec_metrics' of bitbucket.org:mugqic/genpipes into ihec_metrics
       388c9bf ChIP-Seq - Fix IHEC metrics
       b5683f2 ChIP-Seq - Fix IHEC metrics
       7c0f68f Merged in chipseq_atacseq_mode (pull request #121)
       7a101a7 Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       e83caef MethylSeq - Trimmomatic resources revisited for methylseq
       f9bbdb7 General - Genome installation with cvmfs grep
       a8c4d93 ChIP-Seq - Fixing GenPipes metrics
       d82ed89 ChIP-Seq - Sambamba loading properly in mito calcul
       e6c4274 ChIP-Seq - Reducing ihec report ressources
       f3755e9 ChIP-Seq - Fixing typo
       3153017 ChIP-Seq - Increasing IHEC preprocess ressources
       a02fd32 ChIP-Seq - Fix IHEC metrics
       7f1f4a9 ChIP-Seq - Fix IHEC report template md file
       805c58b ChIP-Seq - Fix IHEC report template md file
       19b0b62 ChIP-Seq - Fix IHEC report template md file
       3c573f8 ChIP-Seq - Fixing merge metrics
       1bb96cb ChIP-Seq - Fixing merge metrics
       9714313 ChIP-Seq - Fixing merge metrics
       ad5e531 ChIP-Seq - Fixing merge metrics
       f3ebb74 ChIP-Seq - Fixing merge metrics
       b5658b8 ChIP-Seq - Adding IHEC report template md file
       3ad67e0 ChIP-Seq - Fixing merge metrics
       bef0f18 ChIP-Seq - Fixing merge metrics
       ec3dea0 ChIP-Seq - Fixing merge metrics
       4883d4c ChIP-Seq - Fixing merge metrics
       51bd5f2 ChIP-Seq - Adding merge IHEC metrics
       c310d1e ChIP-Seq - Fixing metrics
       b1ec49d ChIP-Seq - Fixing metrics
       5084171 ChIP-Seq - Fixing metrics
       fa035a1 ChIP-Seq - Fixing metrics
       85c8c53 ChIP-Seq - Fixing metrics
       41081d6 ChIP-Seq - Fixing metrics
       1acdfa5 ChIP-Seq - Fixing metrics
       cfe4059 ChIP-Seq - Fixing metrics & report
       d29881f ChIP-Seq - Fixing metrics
       d1d4385 ChIP-Seq - Fixing metrics
       5e801f9 ChIP-Seq - Fixing metrics
       a1a4784 ChIP-Seq - Fixing metrics
       0f847f1 ChIP-Seq - Fixing metrics
       2f9084a ChIP-Seq - Fixing metrics
       e545fb4 ChIP-Seq - Fixing metrics
       56f6af0 ChIP-Seq - Fixing metrics
       878335c ChIP-Seq - Fixing metrics
       380e122 ChIP-Seq - Fixing metrics
       f6f26f4 ChIP-Seq - Fixing metrics
       343c9ff ChIP-Seq - Fixing metrics
       3da997d ChIP-Seq - Fixing metrics
       db54da1 ChIP-Seq - Fixing metrics
       6a91f0d ChIP-Seq - Adding metrics
       b6a1e60 ChIP-Seq - Fixing bwa missing import
       319655b ChIP-Seq - Fixing chipseq pipeline
       1a47e31 ChIP-Seq - Adding ATAC-Seq protocol
       5f113fc Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       a440d4e General - Genome installation with cvmfs grep
       d7cd595 ChIP-Seq - Fixing GenPipes metrics
       bbfb6bb ChIP-Seq - Sambamba loading properly in mito calcul
       466041c Merge branch 'dev' into chipseq_atacseq_mode
       15611bf ChIP-Seq - Reducing ihec report ressources
       56c39f6 ChIP-Seq - Fixing typo
       f88ee59 ChIP-Seq - Increasing IHEC preprocess ressources
       0bd494d Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       1d727d0 ChIP-Seq - Fix IHEC metrics
       09d7280 Merged dev into chipseq_atacseq_mode
       364e381 ChIP-Seq - Fix IHEC report template md file
       24b64ed ChIP-Seq - Fix IHEC report template md file
       9f8c10e ChIP-Seq - Fix IHEC report template md file
       0d811bd ChIP-Seq - Fixing merge metrics
       70a5ed0 ChIP-Seq - Fixing merge metrics
       b1da600 ChIP-Seq - Fixing merge metrics
       aee5b14 ChIP-Seq - Fixing merge metrics
       74f44e2 ChIP-Seq - Fixing merge metrics
       b096b63 ChIP-Seq - Adding IHEC report template md file
       12b2f31 ChIP-Seq - Fixing merge metrics
       a807a5e ChIP-Seq - Fixing merge metrics
       5131a32 ChIP-Seq - Fixing merge metrics
       60b0f6f ChIP-Seq - Fixing merge metrics
       05c2623 ChIP-Seq - Adding merge IHEC metrics
       9664631 ChIP-Seq - Fixing metrics
       1e2bc38 ChIP-Seq - Fixing metrics
       c6a08a6 ChIP-Seq - Fixing metrics
       e9e96da ChIP-Seq - Fixing metrics
       809cb52 ChIP-Seq - Fixing metrics
       7873647 ChIP-Seq - Fixing metrics
       4b1bb18 ChIP-Seq - Fixing metrics
       e806abc ChIP-Seq - Fixing metrics & report
       ca01c9c ChIP-Seq - Fixing metrics
       3c72d0f ChIP-Seq - Fixing metrics
       1199e5f ChIP-Seq - Fixing metrics
       1ee3d36 ChIP-Seq - Fixing metrics
       59a2828 ChIP-Seq - Fixing metrics
       5157069 ChIP-Seq - Fixing metrics
       4a8a555 ChIP-Seq - Fixing metrics
       2cb7a7b ChIP-Seq - Fixing metrics
       a48c88e ChIP-Seq - Fixing metrics
       235552a ChIP-Seq - Fixing metrics
       5b47faf ChIP-Seq - Fixing metrics
       31535af ChIP-Seq - Fixing metrics
       2d751ea ChIP-Seq - Fixing metrics
       85ea7e4 ChIP-Seq - Fixing metrics
       b62f519 ChIP-Seq - Adding metrics
       bbe8b2c ChIP-Seq - Fixing bwa missing import
       5bdb8aa ChIP-Seq - Fixing chipseq pipeline
       2ae5803 ChIP-Seq - Adding ATAC-Seq protocol

  Paul Stretenowich <pstretenowich@CYPRUS.local>      4 commits

       20aa88f Cleaning mgi.py
       cc25234 Cleaning mgi.py
       397efae Cleaning mgi.py
       97bd776 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      10 commits

       31826de Merged in remove_pacbio (pull request #176)
       68d5195 Merged in poq/fast_module_check (pull request #164)
       652e207 Merged poq/graham_ini into dev
       9c6da9a Merged poq/graham_ini into dev
       1f66c6c Merged update_beluga_ini into dev
       70167f1 Merged update_beluga_ini into dev
       e7c7f34 Merged update_beluga_ini into dev
       a54d824 Merged in add_container (pull request #167)
       8a076c7 Merged in poq/graham_ini (pull request #168)
       f18307f Merged in master (pull request #161)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      33 commits

       241f706 remove pacbio from the repo/release
       fd02e50 no mail in cit.ini
       68d8a66 add sh script for steps in pbs too
       d79b273 Add SARS-CoV2 genome file
       b4c4740 tweek memory usage beluga denovo
       7fcd09a update cluster for rnaseq star index
       473294d use 1.1.0 genpipes_in_container release
       0d591bb make module show sure it raise with older version
       1aa57cb chipseq cedar and graham ini
       1329a13 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6f7be7f more graham ini
       b455c4c more graham ini
       8267603 copy dnaseq ini from cedar to graham
       b448597 fix fastq symlink on graham
       153aa9e less log
       40ab0ff cleanup
       0ac3121 speedup module check
       20ff188 one module call
       1fb3ae7 add mem to star_* on rnaseq
       4d682b6 update again
       19a81cf mem-per-cpu
       e5ce420 Add generic ini for chipseq update README
       d103761 create graham ini file
       05ef09e log when all sbatch submit is dont
       c1f2c53 feedback at submit time
       0ac30b3 feedback at submit time
       fd755d8 fix ini typo
       77b2c9b wrapper for slum, pbs and batch
       bc80a4d update wrap options.
       304c7eb add default wrapper
       22d9efa remove docker
       004a9a4 put --wrap import at the top
       3516da4 add wrapper to all pipelines

  P-O Quirion <pioliqui@gmail.com>      2 commits

       328a00e make nanopore executable
       ab1de9a Add automatic wrapper option

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       123f6c5 Merged in fix-watch-folder (pull request #165)

  Rom Grk <romgrk.cc@gmail.com>      1 commits

       1baeef5 watch_portal_folder.py: fix undefined variable

  ufgauthi <ulysse.fortiergauthier@mcgill.ca>      1 commits

       0a60ef2 Bug Fix by replacing sacct delimiter | by ^

  Ulysse Fortier Gauthier <ulysse.fortiergauthier@mcgill.ca>      1 commits

       bf54dbf Merged in ufg_log_report_fix (pull request #157)

3.1.5        Wed Jan 15 11:58:16 2020 -0500        424 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      78 commits

       e0844c3 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       986f2a3 HICSeq pipeline - improving trimmomatic resources for SULRM + added graham config file for hicseq pipeline
       6714a5d Genome installation - added America Mink genome (Neovison_vison.NNQGG.v01) installation script
       d60d136 Software installation - added 'demuxlet' installation script
       6822cbf Software installation - updated ucsc.sh with lat4est version i.e. v387
       b350d2d Software installation - replaced call of lsb_release by universal commands (avoid 'lsb_release command not found' error)
       5c87fd7 Software installation - corrected regular expression within genome_digest function
       23e9315 removed some test code introduced by latest commit...
       9445ee9 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dev
       00577d4 GenPipes RNA-Seq de novo - corrected bugs introduced in commit 4a72735
       f3b9a92 GenPipes BFX - added bash_cmd python wrapper to wrap basic bash commands
       4a72735 GenPipes RNA-Seq de novo - updated pipeline with better sample assignemnt to jobs for better JSON building
       d89d4fc Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1b48215 GenPipes JSON - updated pipeline.py and scheduler.py so that the copy of the original JSONs to the PORTAL is made before building the pipeline steps
       6360349 GenPipes Trimmomatic - corrected assignment of samples to the job to correct analysis JSON generation per sample
       fbaa84a MethylSeq - corrected assignment of samples to the jobs to correct analysis JSON generation per sample
       d1cd79b GenPipes analysis JSON - updated jsonator.py in reagrds of the updated .base.ini files, now containing 'source' & 'version' to ease the process + minor updates to core/pipeline.py
       fc1b6c0 GenPipes config files - updated most of the .base.ini files with missing 'source' & 'version' to avoid issues with the 'jsonator'
       f44ccb1 GenPipes Anlalysis JSON file - corrected core/pipeline.py regarding the use of PORTAL_OUTPUT_DIR
       9751513 GenPipes Analysis JSON file - added the system to update the submission_date
       a9ba832 GenPipes Anlalysis JSON file - pipelines now create JSON analysis file as a default behavior
       7b42d15 GenPipes Analysis JSON file - added the submission_date to the JSON
       7e04404 Analysis JSON - updated jsonator.py with project_name and submission_date. Also updated version to 1.0.1
       f422649 GenPipes utils - minor updates : updated some comments
       38937e6 Software install - updated install_module.sh with finer LIBDIR for patching and better patching to avoid overwritting potential pre-existing RPATH
       fb2bffc Software install - updated R_Bioconductor.sh with new R packages and finer LIBDIR for patching
       721d941 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       074045a GenPipes Sanity Check - Refining the report
       0d0e0bd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       1478455 GenPipes Sanity Check - Adjusting the log level
       cf3ec9e GenPIpes Sanity Check - set log level to 'warning' instead of info when running sanity check, thus removing some useless if statements and refining display of the messages
       6cf36ee GenPipes Sanity Check mode : removed some useless comments in core/pipeline.py
       b46df3f GenPIpes Sanity Check mode - created the SanitycheckError class in config.py to use instead of Error
       f4d6323 Merge branch 'master' of bitbucket.org:mugqic/genpipes into rnaseq_denovo_jhg
       989996a GenPipes MethylSeq - updated base.ini with [samtools_cram_output] section
       066c24e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       433950d GenPipes Sanity Check - DNASeq high coverage, RNASeq denovo assembly & RNASeq light pipelines are updated regarding the sanity-check mode
       0fbecd2 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       6b00ce5 GenPipes Sanity Check - updated pipelines for sanity-check mode responsiveness
       afa406a GenPipes Sanity Check - updated common.py to be sanity-check responsive
       4f25a0f GenPipes Sanity Check - update _raise() function in config.py & refined code of sanity-check mode in pipeline.py
       6a457d5 GenPipes Sanity Check - updated readset.py design.py to be sanity-check responsive & tools.py with proper import statements
       6f11825 MethylSeq pipeline - updated pipeline to make SINGLE-END mode actually work
       b3ac5a9 GenPipes Sanity Check mode : updated PacBio Assembly pipeline as the first try to test sanity check mode
       f288544 GenPipes Sanity Check mode : updated pipeline.py with the sanity check mode fully functionnal
       b5a78dd GenPipes Sanity Check mode : updated sample.py to reflect the updates that have been done in config.py
       3203ce9 GenPipes Sanity Check mode : updated config.py to avoid raising errors when sanity-check mode is on, logging them instead
       daabbf4 GenPipes code cleaning - cleaned some 'import' calls : stop using 'import *' and specify which modules/classes/functions to import instead
       2699bc9 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       975e0ff GenPipes Core - updated pipeline.py to catch , relative to changes in config.py. Added the sanity-check mode : the pipeline now does not stop at fist met error/exection : wait untill the end to print the errors
       1f628d3 GenPipes Core - updated config.py wth PO's code : catches  instead of
       0f4c860 Software upadte - updated mugqic_tools installation script with version 2.2.3
       68f32a9 Software update - updated R_Bioconductor.sh : added binless and DoubletFinder packages to the installation, updated the installer depending on the version of R
       2bca735 Saoftware update : update cellranger-atac.sh with the latest versoin 1.1.0
       1478321 Software update - added th shebang to the C3G wrappers for installed binaries
       697eaec Genome installation - corrected typo in 'lambda_page'
       e7010eb Software update - added script to install gemBS in C3G softwatre stack
       9c4f295 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       93f92cd Software update - updated installation scripts for LUMPY-SV, mugqic_tools & SMRTLink
       39b691a Merge branch 'master' of bitbucket.org:mugqic/genpipes
       06fdb58 Software update - updated MultiQC installation script : shebang of MultiQC scripts now uses #!/usr/bin/env python
       e66f717 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1e4c96f Software update : added the installation of methylKit packages in R_Bioconductor.sh - also added the generation of a file listing all the packages installed for all the new installed R versions
       251f418 GenPipes utils - updated nanuq2mugqic_pipelines script : added support for iSeq projects
       742aa73 GenPipes utils - updated nanuq2mugqic_pipelines script : added -rn/--run parameter, standing for Nanuq run ID, to fetch only readsets procesed in specified run(s)
       5555171 Software update - updated install_module.sh so that it does not wrap nor patch the executable binaries when installing on DEV space
       d4fb87c Software update - added MiXCR v3.0.5 installation script - MiXCR: a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data
       f020d92 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3e7d563 Software update - updated some installation scripts with the latest version of software : bowtie v1.2.2 - bowtie2 v2.3.5 - MultiQC v1.7
       09335ff install_genome.sh - added the generation of the 100-bin GC content file for bisulfite genome references
       4e1331f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       04047df Genome installation : added Bisulfite Genome Reference creation and indexing with Bismark
       8916af0 updated LIBDIR in R_Bioconductor.sh
       0ad8456 update Bismark version to 0.21.0
       eff8002 updated python installation script with version 3.7.3 and added pip as a symling to pip3
       db27d82 Added fastq software installation script
       e1f627f Version bump to 3.1.5-beta
       1786fb3 Version bump to 3.1.4

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      1 commits

       650797a GenPipes - DEBUGGING - DNASeq SamToFastq & SymLink steps corrected + working bfx/bash_cmd.py

  Édouard Henrion <henrione@gra-login2.graham.sharcnet>      4 commits

       9e9adc7 GenPipes JSON - debugged call to job2json when output_dir is different that '.'
       b7eb0ed Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       754a7d9 GenPipes software update - updated patching in install_module.sh
       5bf7663 GenPipes software update - updated R_Bioconductor.sh with the latest requested libraries & updated patching

  Édouard Henrion <henrione@ip18.m>      3 commits

       140b411 GenPipes Scheduler - Corrected bug in job2json call
       c77a535 Merge branch 'eh_methylseq_single_end' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       a592ddf GenPipes MethylSeq - updates for single-mode

  Édouard Henrion <henrione@ip20.m>      1 commits

       41ed7a2 Genome Installation - updated install_genome.sh with common grep version

  ehenrion <edouard.henrion@mcgill.ca>      83 commits

       25f72f0 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1578931756922 (pull request #163)
       64b7414 DNASeq.py - removed the use of 'os/path.abspath' in call of 'ln()'
       3c750f0 bash_cmd.py - added import os
       7a05d36 DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable
       e9213a6 bash_cmd.py - remove call to "sleep" in the ln() function. Replaced it by a call to "ls" in order to flush the cache
       0115491 GenPipes - BUG correction - trailing space removed in output files of DNASeq skewer trimming step
       1fdf579 Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578602990982 (pull request #162)
       b74f38d BUG correction : corrected bash_cmd.py "ln" function with correct format keys when buiding command
       07d826f Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578588554527 (pull request #158)
       828d573 Bug Correction : corrected 'ln' function in bash_cmd.py to avoid a python "TypeError: cannot concatenate 'str' and 'int' objects" exception
       4205582 Merged in ehenrion/dnaseqbelugaini-edited-online-with-bitbu-1576097030188 (pull request #147)
       fcaef48 GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       5202100 Merged in eh_quick_fixes (pull request #143)
       815632f Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918024470 (pull request #141)
       92a7c11 GenPIpes - dnaseq.py : bug correction - typo removed
       44d9eb5 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575399438745 (pull request #132)
       1526642 Merged in ehenrion/commonpy-edited-online-with-bitbucket-1575398053532 (pull request #131)
       482a0bb GenPipes - bug correction in pipelines/common.py : corrected name for bam and fastq files created by the pipeline : don't depend and raw file names anymore, but built from sample & readset name given in the readset file
       2a25960 GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       72e7f1e GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       c946bb0 Merged in ehenrion/schedulerpy-edited-online-with-bitbucket-1575392160612 (pull request #129)
       9bf49cb Merged in ehenrion/nanuq2mugqic_pipelinespy-edited-online-w-1575392552042 (pull request #130)
       3893e10 GenPipes : nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       2f38d5d GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       a876d02 GenPipes - corrected scheduler.py with missing argument line
       8b8bca8 Merged in dev (pull request #125)
       7cd7b1c RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - dowgraded trinity version to 2.0.4_patch
       656db12 RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - corrected typo in trinity version...
       01c2bee RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - updated trinity version to 2.2.0_patch which contains a patch (from C3G developppers) to avoid  printing buffer 'hiccup'
       b9325d1 CIT - RNASeq Denovo Assembly - update default cluster_walltime to 4:00:00 in cit.ini
       8235190 Merged ehenrion/cit-dnaseq-updated-cit_cluster_walltim-1574180387297 into dev
       238abec CIT - DNASeq - updated 'cit_cluster_walltime' for gatk_callable_loci step in cit.ini
       6a891c2 DNASeq - gatk_callable_loci - adjusted memory and cpu allocation in dnaseq.beluga.ini
       dec27c8 CIT - DNASeq - Adjusted trimmomatic resource allocation through java_other_options, threads settings and mem-per-cpu use
       3d34243 CIT - RNASeq - corrected 'cluster_walltime' for wiggle  step in cit.ini
       a24e7b1 RNASeq - wiggle step - adjusted memory allocation using 'mem-per-cpu' instead of 'mem'
       3e3ee1d CIT - RNASeq - updated java threads to 5 through 'java_other_options' in rnaseq.beluga.ini
       1a1d0f1 CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       bae00cd CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       cc5df68 CIT - HiCSeq - redefined hicup_align walltime and mem-per-cpu in cit.ini
       5cd0d02 CIT - Pacbio assembly - set specific walltime for smrtanalysis_summarize_polishing in cit.ini
       1c19791 HICSeq pipeline - trimmomatic resources revisited in hicseq.beluga.ini : removed buffer_size
       29024f2 CIT - RNASeq denovo Assembly - set localfit to true in cit.ini for differential_expression_deseq
       622b41a CIT - Pacbio Assembly - redefined walltime for pacbio_tools_split_reads in cit.ini
       09e2235 RNASeq denove Assembly - edited rnaseq_denovo_assembly.beluga.ini adjusted memory assignment for insilico_read_normalisation_readsets
       3a25fef RNASeq denovo Assembly : edited rnaseq_denovo_assembly.beluga.ini with beter resources assignement
       b4a3e98 Software update - trimmomatic.py - updated trimmomatic command with the use of the  java_other_options parameter provided by the ini files
       47e5e41 Merged in ehenrion/rnaseq-metricspy-ihec_metrics_rnaseq--1573489397505 (pull request #124)
       87c3190 RNASeq - metrics.py - ihec_metrics_rnaseq : added file in the input_file list to correct job dependency
       fc680b1 RNASeq cufflinks - correcting dependencies for ihec_metrics : needed rnaseqc report file to be added to the output_file list
       25f6cc3 CIT - Pacbio assembly - set threads in beluga.ini
       aa3d185 CIT - RNASeq_light - set threads for trimmomatic in beluga.ini
       b199242 CIT - RNASeq_light - redefined walltime for callisto_count_matrix in cit.ini
       24d79ff CIT - Pacbio Assembly - redefined walltime for preassembly in cit.ini
       3c8f736 CIT - RNASeq denovo Assembly - redefined walltime for align_estimate_abundance in cit.ini
       3134f05 CIT - RNASeq denovo Assembly - updates cit_cluster_walltime to 24h
       24351cb CIT - RNASeq denovo Assembly - redefined walltime for transdecoder and align_estimate_abundance in cit.ini
       c17d791 CIT - Pacbio Assembly - redefined walltime for smrtanalysis_filtering in cit.ini
       19ed125 CIT - DNASeq - redefined walltime for sambamba_merge_realigned in cit.ini
       d92aa27 CIT - DNASeq High Coverage - redefined trimmomatic walltime in cit.ini
       74979fa CIT - HiCSeq - redefined trimmomatic walltime in cit.ini
       ed4fb13 CIT - RNASeq_light - redefined walltime for trimmomatic in cit.ini
       0cf351d RNASeq denovo Assembly pipeline - rnaseq_denovo_assembly.py - corrected command quoting in trinity step
       8386936 RNASeq denovo Assembly Pipeline - rnaseq_denovo_assembly.beluga.ini - corrected max_memory parameter for trinity step : using Gb instead of Mb
       373bd60 RNA-Seq de novo assembly pipeline - rnaseq_denovo_assembly.beluga.ini - corrected Jellyfish memory parameter : using Gb instead of Mb unit
       d2e4c46 HiC-Seq pipeline - hicseq.base.ini - corrected typo in HiCUP module name
       ceaae53 HiC-Seq pipeline - hicseq.base.ini -  update HiCUP version to v0.7.2
       7860ea8 dnaseq.py - corrected file names dependencies
       c198c42 dnaseq.py - corrected wrong file extension in metrics_vcf_stats step
       bc952b6 HiC-Seq pipeline - corrected typo in cram_output options definition - hicseq.base.ini
       d8c5f58 rnaseq_light.cedar.ini : corrected typo in cluster_server assigned value
       19bf83c rnaseq_light.mp2b.ini corrected minor typo on cluster_server
       27f7080 rnaseq_denovo_assembly.py : corrected differential expression jobs with samples
       5010b02 corrected typo in ampliconseq.beluga.ini
       1f16017 DNASeq High Coverage README.md updated with `picard_fixmate` step documentation
       eda0617 GenPipes - Readme updated for `picard_fixmate` step in dnaseq_high_coverage.py
       a05f490 Merged in eh_methylseq_single_end (pull request #115)
       6c88ebf Merged in sanity_check_mode (pull request #114)
       fa25623 rnaseq_denovo_assembly.beluga.ini :  updated [insilico_read_normalization_all] with missing core allocation
       305d96e README.md updated : removed Guillimin settings section from README
       56d668c README.md updated : corrected the setting of MUGQIC_INSTALL_HOME_DEV for mp2b
       8d4f1b4 rnaseq_denovo_assembly.mp2b.ini - corrected resources requirments : all the steps now run on one single node !
       33909cf chipseq.base.ini - updated MultiQC version to 1.6 to support Homer

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       6c76ddc Merged in lefebvref/rnaseq_denovo_assemblybaseini-edited-onl-1557891680322 (pull request #110)
       53bd32a Previous default expression in base.ini would end up basically only retaining transmembrane proteins. Not good.

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      2 commits

       5bea83c Corrected error with beluga ini for RNAseq denovo
       02b131c Corrected minor error in the help message that said that the default protocol was cuflinks. Starting from version 3.1.4, stringtie is the default protocol

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      10 commits

       153d4c4 Merged in Jose-Hector-Galvez/update-to-the-rnaseq_light-readmemd-file-1569265173884 (pull request #120)
       e24a066 Update to the rnaseq_light README.md file to address issue brought up by GSoD collaborator.
       5f4086f Merged in Jose-Hector-Galvez/readmemd-edited-online-with-bitbucket-1565184718156 (pull request #118)
       ed95482 README.md edited to remove warning labels.
       341cab2 Merged in rnaseq_denovo_jhg (pull request #103)
       05e7dfe Merged in rnaseq_jhg (pull request #112)
       9a24bea Merged dev into rnaseq_jhg
       6f98c76 Merged master into rnaseq_denovo_jhg
       5505147 Merged master into rnaseq_jhg
       160079d Merged master into rnaseq_denovo_jhg

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      3 commits

       e618184 Corrected error in the beluga ini for RNAseq de novo
       e9610a9 Corrected error in the beluga ini for RNAseq de novo
       ceade42 Corrected error in the beluga ini for RNAseq de novo

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       cc5eac7 add header in the cram file (all base.ini)
       92b5c4f add cram to tumor_pair pipeline
       7822413 implement the launch of the entire pipeline if no step argument is given
       689732a add cram to dnaseq_high_coverage pipeline
       71569ce add cram to methylseq pipeline
       c242a05 add cram to rnaseq pipeline
       70ccb03 add cram to hicseq pipeline
       0e78873 add cram to chipseq ini files
       830cd38 add cram to dnaseq ini files
       f9d780a add cram to ChipSeq
       5f1d736 correct bugs in dnaseq light
       0313f9c add cram creation to dnaseq pipeline
       bfde214 add generic function to create CRAM from BAM
       1956514 make samtools view's output not mandatory removable

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      2 commits

       52cd608 temporary fix for tbi missing output in haplotypecaller when running only 1 job
       6a98813 Merged in cram (pull request #111)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      9 commits

       ee86893 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0fd7a88 Nextflow install script
       cc8694c Merged in rnaseq_rnaseqc (pull request #116)
       d59ea5e General - Correcting genome ini path
       caa989b Merge branch 'dev' of bitbucket.org:mugqic/genpipes into rnaseq_rnaseqc
       8e90fcf RNA-Seq Pipeline - Typo update
       6bfe800 RNA-Seq Pipeline - Typo update
       bb8a0cd RNA-Seq Pipeline - Removing verify_bam_id step
       7eb13c4 Chip-Seq Pipeline - Adding homer_make_ucsc_file to mp2b.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      12 commits

       a8d1aae Merged in add_container (pull request #122)
       9b20dde Merged in cit_fix_cedar_2019-04 (pull request #109)
       4eb85d8 gatk4.py added 2 spaces in multi line string
       ef09195 ampliconseq.py edited online with Bitbucket
       9e04a50 pacbio_assembly.cedar.ini edited online with Bitbucket
       2b02306 Merged in log_repport_ci (pull request #108)
       04e9b6b Merged in log_repport_ci (pull request #107)
       f171ba1 Merged in log_repport_ci (pull request #106)
       0da3f3e Merged in log_repport_ci (pull request #105)
       f9458d8 Merged in log_repport_ci (pull request #104)
       1b7e74f Merged in fix_multiqc (pull request #89)
       0e66455 Merged in fix_multiqc (pull request #77)

  Pierre-Olivier Quirion <poq@beluga3.int.ets1.calculquebec.ca>      2 commits

       feeaf91 Log when loading Job. More robust output parser.
       5181b5c Log when loading Job. More robust output parser.

  P-O Quirion <pioliqui@gmail.com>      158 commits

       161d1d1 fixed motifMaker ini: used default memory
       ba32a47 typo in waltime kallisto_count_matrix
       50b824b trimmomatic thread ini
       20c0f56 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       795edd9 string tie ini
       d4e255d sync thread in java for trimomatic
       9bd93fa differential_expression_and_goseq_rsem return list list
       42bf1d8  fix scheduler logging; reduce sleep to 0.1 sec
       0902f42 ignore .pyc
       bca87cc remove the methylkit_differential_analysis step form the pipeline
       f9cef1b add gitignore
       dd0267e missing import
       3d6b4e7 log regression scheduler
       f4ab24a fix beluga ini
       6d5d769 fix beluga ini
       ef68329 Merge branch 'add_container' into merge
       a872d63 Merge branch 'master' into dev
       d31246f fix batch system account for container
       313fcb2 remove commented code
       281e3cd Merge branch 'add_container' of bitbucket.org:mugqic/genpipes into add_container
       ebb67ba Add default config for container ini
       116615f fix regression
       7542a09 fix scheduler when container not in use
       8c9d6e3 use mugqic stack java in rnaseq on cedar.ini
       bb6fc69 simple default ini for container
       b82cace default ini for container
       331715d Basic container ini file
       b787ff7 Working singularity version tested on graham
       2fe8f47 WIP exec line prototype
       b2b4d97 add container option
       215d352 Merge branch 'cit_fix_cedar_2019-04' of bitbucket.org:mugqic/genpipes into cit_fix_cedar_2019-04
       4092e77 amplicon R tweak
       7b389fc missing quotes in pacbio
       138f335 amplicon memory
       32627c8 amplicon memory
       141554d sleuth again wrong name
       8c34b2b renaming sleuth to differential_expression_sleuth
       963d105 fix typo in ampliconseq ini
       ff92d21 fix amplicon method typo
       da2cdf8 tweak default mem for rnaseq denovo and light
       75e2f69 more ampliconseq deadlock: otu_assigning
       1cfac73 parallel_assign_taxonomy_uclust.py deadlock
       6ecfef8 rna fix
       cecb514 update cit
       b4997b3 packbio tmpdir again
       c392576 chr19 option for cit chicagorun
       46d9149 fix log report on output id
       5cfc33c fix log report on output id
       0cc3e17 update waltime in amplicon and regex in log
       5cc0ce5 migrate form testdata ini to cit
       2c02e92 mem teak
       51bf13d tweak to get in bynode queue
       3f05dbe tweak ini for cit
       078d0e0 dnaseq ini tweak for slurm; sh_create_baitmap \\
       cf79f46 timout amplicon
       410185e remove mpileup_metrics_vcf_stats step from mpileup
       9d100d8 pbutgcns_wf need full path an tmp need to exist
       ad31e6e ini file tweak
       ac6bf87 tweak mem usage
       b4f2cf8 ini tweak plus perl module load in mummer
       3ae1b09 tweak mem usage per-cpu for ampliconseq
       bb8fc31 Ensembl v 90 for GRCh38
       890cd48 mummer using deprecated perl syntax
       16326cf more default time for rnaseq cufflink
       3d1aa24 mugqic_tools/2.2.4 in rnaseq_light
       ffcf62f samll sample for both gatk call in variant_recalibrator
       7b25d7b change R_Bioconductor version/cleanup code
       bc97eae update rnaseq_light references
       505ffeb tweek memory usage in dnaseq_high_coverage
       3ce4bb6 add time to ampliconseq
       c1b5a6a reformat errors in the core.config module
       05aa94b redirect small_sample_option to the right function
       3d51765 tweak cit for dnaseq
       720c3fa fix ampliconseq dada2
       390ceee tweak VariantRecalibrator for testing
       f74afdf fix slurm_time_to_datetime
       733f258 set cit ini file and config
       7d9c898 add runtime to csv output
       86e84f3 Add GENPIPES_CIT env var to run in testing mode
       4cb6e24 fix all slumtmpdir sting renaming problem
       8f3652f update mugqic_tools in hicseq
       2107daa second run of pipeline fixing
       458e919 remove extra gz extention from anotation
       690acf4 flash_log is a list of one
       b829c31 update hicseq.base
       81a4690 add interactive option in smrtanalysis bash cmd
       7855217 fix regression
       18aa48b fix scheduler when container not in use
       ae7824c use mugqic stack java in rnaseq on cedar.ini
       b8eee7f simple default ini for container
       2522d29 default ini for container
       60de9f1 Basic container ini file
       a5b2901 Working singularity version tested on graham
       5229f0e WIP exec line prototype
       56e0a9b add container option
       27a2cb3 amplicon R tweak
       738c59d missing quotes in pacbio
       d73e165 amplicon memory
       a770708 amplicon memory
       80498da sleuth again wrong name
       53c5205 renaming sleuth to differential_expression_sleuth
       1e5d20e fix typo in ampliconseq ini
       7bb900d fix amplicon method typo
       85b51a0 tweak default mem for rnaseq denovo and light
       7df9cae more ampliconseq deadlock: otu_assigning
       777bb33 parallel_assign_taxonomy_uclust.py deadlock
       24bd06e rna fix
       a9e2370 update cit
       b120e4c packbio tmpdir again
       8e1835f chr19 option for cit chicagorun
       eec3f95 fix log report on output id
       5488adc fix log report on output id
       34960c0 update waltime in amplicon and regex in log
       ca1fe4d migrate form testdata ini to cit
       eca1836 mem teak
       058f071 tweak to get in bynode queue
       10ff44a tweak ini for cit
       c7568c6 dnaseq ini tweak for slurm; sh_create_baitmap \\
       1c50e41 timout amplicon
       ff1a9c7 remove mpileup_metrics_vcf_stats step from mpileup
       d8d8a3d pbutgcns_wf need full path an tmp need to exist
       48ff75b ini file tweak
       68d6f46 tweak mem usage
       160a123 ini tweak plus perl module load in mummer
       ee6a7aa tweak mem usage per-cpu for ampliconseq
       f0acafe Ensembl v 90 for GRCh38
       23788c3 mummer using deprecated perl syntax
       14d0500 more default time for rnaseq cufflink
       4e3a8b0 mugqic_tools/2.2.4 in rnaseq_light
       624781c samll sample for both gatk call in variant_recalibrator
       8307e59 change R_Bioconductor version/cleanup code
       069dd4a update rnaseq_light references
       3e01c15 tweek memory usage in dnaseq_high_coverage
       c0c5905 add time to ampliconseq
       d29d926 reformat errors in the core.config module
       43fdeda redirect small_sample_option to the right function
       cfefb02 tweak cit for dnaseq
       509d333 fix ampliconseq dada2
       927a732 tweak VariantRecalibrator for testing
       0a2659d fix slurm_time_to_datetime
       0f48c2d set cit ini file and config
       20c8521 add runtime to csv output
       02ce535 Add GENPIPES_CIT env var to run in testing mode
       dbbb0cd fix all slumtmpdir sting renaming problem
       19c6c30 update mugqic_tools in hicseq
       27e5164 second run of pipeline fixing
       32a1a40 remove extra gz extention from anotation
       4edbfc5 flash_log is a list of one
       a396c5f update hicseq.base
       d806fce add interactive option in smrtanalysis bash cmd
       d9567ac add space
       31ddd6c add check for output job id
       63431c4 get real path in log_report
       bec1796 file exist file missing mixup
       4120c3a bug fixes
       3e94acd get pipeline info in tsv file can also mute stdout
       29f5670 run dnaseq multiqc on all samples at once
       8d73217 add indentation

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       43f70c6 fix to preprocess step

  Robert Syme <robsyme@beluga4.int.ets1.calculquebec.ca>      1 commits

       8581436 We don't need to track compiled python .pyc in git

  Robert Syme <rob.syme@gmail.com>      7 commits

       30ea46d Merged in spacesfix (pull request #117)
       0c4d6ce Remove spaces from Rscript commands.
       1f3c0f6 More trailing space cleanup.
       534124c Cleanup trailing spaces
       16b7d01 Merged in amplicon-pathfix (pull request #113)
       63d804a Add support for --output-dir argument
       eb8389c Trailing whitespace cleanup

  Rola Dali <rola.dali@mail.mcgill.ca>      31 commits

       92202e8 dnaseq_high_coverage.base.ini edited Varscan version to 2.4.3
       3f0b732 dnaseq.mp2b.ini edited online with Bitbucket
       07b6376 dnaseq.beluga.ini edited online with Bitbucket
       a0bad07 dnaseq.cedar.ini edited online with Bitbucket
       f379676 dnaseq.beluga.ini edited online with Bitbucket
       2805dc6 dnaseq.cedar.ini edited online with Bitbucket
       a4c3418 pacbio_assembly.beluga.ini edited online with Bitbucket
       bdd8b72 dnaseq.beluga.ini edited online with Bitbucket
       b648fb5 dnaseq.cedar.ini edited online with Bitbucket
       287e82b dnaseq.beluga.ini edited online with Bitbucket
       a760a50 dnaseq.cedar.ini edited online with Bitbucket
       62de3fb hicseq.beluga.ini edited online with Bitbucket
       1fd7787 hicseq.cedar.ini edited online with Bitbucket
       42eb77a ampliconseq.cedar.ini edited online with Bitbucket
       fabb3d5 ampliconseq.mp2b.ini edited online with Bitbucket
       193c97b ampliconseq.beluga.ini edited online with Bitbucket
       5c4c908 ampliconseq.cedar.ini edited online with Bitbucket
       59491ff dnaseq.cedar.ini edited online with Bitbucket
       9c9a75f dnaseq.beluga.ini edited online with Bitbucket
       a2e7e80 rnaseq.beluga.ini edited online with Bitbucket
       d1777dd hicseq.beluga.ini edited online with Bitbucket
       5a6a198 hicseq.cedar.ini edited online with Bitbucket
       855a02c hicseq.cedar.ini edited online with Bitbucket
       a4cb4c1 hicseq.beluga.ini edited online with Bitbucket
       efe78e4 README.md edited online with Bitbucket
       c0b67c1 README.md edited online with Bitbucket
       5c78251 README.md edited online with Bitbucket
       982f934 README.md edited online with Bitbucket
       4402aef README.md edited online with Bitbucket
       ed0fd9e README.md edited online with Bitbucket
       b97304e README.md edited online with Bitbucket

3.1.4        Tue Mar 26 14:03:32 2019 -0400        198 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      69 commits

       cea5856 Updated all the GenPipes .base.ini files with the latest verison of mugqic_R_packages i.e. 1.0.6
       892cef5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       cc01f4e MethylSeq - added reference to Bismark in markdown config_and_references.md file for a better report document - also updated the metrics table assignement in MethylSeq.ihec_sample_metrics_report.md and bfx/report/MethylSeq.ihec_sample_metrics_targeted_report.md
       540ee0e Merge branch 'master' of bitbucket.org:mugqic/genpipes into slurm_report
       3ab4a6f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f8c58f3 updated mugqic_tools version to 2.2.2 in ampliconseq.base.ini
       4a100ec updated R_mugqic_packages.sh with version 1.0.6 of the mugqic R packages
       1d490a7 updated pipeline inis file with the latest version of mugqic_tools i.e. v2.2.2
       d1ad834 MethylSeq - MethylKit DMR - to fit with C3G code standards, getting rid of bfx/methylkit.py and use bfx/tools.py instead to call methylkit.R, which is part of mugqic_tools
       d5e1184 MethylSeq - methylKit DMR - correted typo in bfx/methylkit.py
       a989966 MethylSeq - methylKit DMR - correted typo in methylseq.base.ini
       ba9c8b1 MethylSeq - MethylKit DMR - change the call to methylKit.R : now using R instead of Rscript
       c6c82bf Pac Bio Assembly pipeline - updated the bfx wrapper circlator.py : make sure th eoutput is removed before launching circlator, otherwise it will throw an error
       a763f6d Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       9650b32 PacBio Assembly pipeline - minor update in pacbio_assembly.py for code reading purposes
       e41a72e Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       0ddffdd updated svaba.sh installation script with latest version of SvABA : 1.1.0
       43e7527 MethylSeq pipeline - DMR analysis - updated R module version in the base.ini file to make sure methylKit is available
       3576f8d Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       90d9206 updated mugqic_tools version to 2.2.2 in methylseq.base.ini
       34f4ccb update mugqic_tools.sh with the latest version : 2.2.2
       f1a34d0 updated install_module.sh to better handle LD_LIBARAY_PATH in the LIBDIR
       0fb5b2c Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0aa9f92 corrected Octopus installation script
       a5a7bd9 Pacbio Assembly pipeline - corrected basemodification and MotifMaker bfx subroutine wrappers : now sourcing /etc/setup.sh
       bfefcf7 modified R module in install_genome.sh : now using mugqic/R_Bioconductor/3.5.1_3.7
       805760f updated mugqic_tools.sh with version 2.2.1
       02d4ef8 Adding installation scripts for CMake and Octopus
       10a9042 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       193fe28 MethylSeq pipeline - DMR analysis - change file namings in prepare_methylkit step
       f6f936f MethylSeq pipeline - DMR analysis - adjusting walltime for filter_snp_cpg step
       aca1e34 MethylSeq pipieline - DMR analysis - clean code in bfx/tools.py
       732a39a updated wrapping and patching in install_module.sh
       56891fa updating software install scripts with newer versions
       82723b7 bedtools.py - removed unnecessary argument to bedtools.coverage subroutine
       d126c1b MethylSeq pipeline - methylseq.py - fixing inputs assignment in IHEC metrics step
       4eead82 MethylSeq pipeline - methylseq.py - correcting sort_sam call and affecting job name to bismark methylation call job
       27f85ae Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e23a982 BFXDEV-429 - updated syntax of temporary diretory after sourcing of /etc/setup.sh for smrtanalysis tools
       53e5131 updated resources/modules/R_Bioconductor.sh script : refine LIBDIR definition - revised wrappers creation
       3f20912 updated install_module.sh : refine create_c3g_wrappers subroutine to catch executable binaries - updated LIBDIR definition for centOS cases
       58286be updated delly.sh script : lasted version 0.8.1 & added DELLY_PATH in the modulefile
       8a1084f udpated lumpy-sv.sh script : added a sed command to make lumpyexpress.config compliant with our environment
       1c7b2c4 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       58bde8f PacBio Assembly - updated bfx/smrtanalysis.py as part of the solution of BFXDEV-429
       7ee6531 Updated syntax to tmp_dir to finish fixing BFXDEV-429 - first usefull related update was made in commit ad3caf1, when bumping to version 3.1.3
       807bee4 minor - updated syntax
       e93e06f Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       dafb3c0 removing occurences of DEV paths within the ini files of GenPipes
       3d256c1 Adjusting .ini section for some htslbfx subroutines which were pointing to [DEFAULT] section of the .ini file
       3fe0d06 BFXDEV-752
       25bc30d Adjusting .ini section for some qualimap ray sambamba and other subroutines which were pointing to [DEFAULT] section of the pipeline .ini files
       6d4dd16 corrected some typo in jsonator.py to have fastq & bam file paths properly handled
       630de05 Adjusting .ini section for some htslib subroutines which were pointing to [DEFAULT] section of the .ini file
       5b03880 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       891780f comment line removed in install.module.sh
       f6720fb Added ctpr.sh script to install CTPR package in C3G software stack
       0a3e8ab Adjusting .ini section for some blast * blat subroutines which were pointing to [DEFAULT] section of the .ini file
       3ebb9cc BFXEDV-529 - edited RNAseq to consider star_cycle_number instead of cycle_number for star.index step - also ajusted read_size value for sjdbOverhang value completion
       54fd040 Homo_sapiens.GRCh38.sh - Added the creation of a Ribosomal RNA interval list file, for use of Picard Tools CollectRnaSeqMetrics
       fe136f7 Adjusting .ini section for some bcftools subroutines which still were pointing to [DEFAULT] section
       709347d Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       4344c94 changing samples assignment in trimmomatic jobs for a better analysis JSON creation
       321c65d Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       24e90fa Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9cc90a8 added genpipes.sh
       0157c92 updated GATK installation script
       ea08443 Version bump to 3.1.4-beta
       7e318dc Version bump to 3.1.3

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       9cc62b8 MethylSeq Pipiline - MethylKit DMR analysis - updated Cedar .ini file with relevant cluster resources, changed filter_snp_cpg subroutine int bfx/tools.py to use bedops instead of bedtools, added bedops module in the .base.ini file
       e26c2af MethylSeq Pipeline - updated filter_snp_cpg call with more cluster resources and removed pipes to avoid 'Broken pipe' error...
       7d77cd8 MethylSeq pipeline - DMR analysis - fixing cedar ini file + cleaning code before pull request

  ehenrion <edouard.henrion@mcgill.ca>      14 commits

       3ef3a04 pipeline.py  : edited the help so that it actually shows that SLURM is the default scheduler used by GenPipes
       d3edaa6 methylseq.py corrected report files in metrics step
       86646c4 Merged in slurm_report (pull request #73)
       b2dbb46 job.py : corrected typo
       0961c1e job.py : added missing setter class for name and multiqc_files attributes
       4f48cc6 Merged in methylseq (pull request #67)
       10bcac6 dnaseq.cedar.ini : resolved conflicts
       a18369a rnaseq_denovo_assembly.base.ini : corrected missing variables
       d80c178 tumor_pair.base.exome.ini - updated COSMIC path
       1a66a54 tumor_pair.base.ini - updated COSMIC path
       a3c68a3 rnaseq_light.base.ini : commetn adapter_fasta because it was pointing to some old place
       0f799c3 BFXDEV-737 - updated ucsc.py for a better handling of Danio Rerio GRCz11 genome build
       b768d04 BFXDEV-734 - updated ucsc.py for a better handling of Mus Musculus GRCm38 genome build
       54dab24 AmpliconSeq - reassign silva_db to correct path for dada2 analysis - ampliconseq.base.ini

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      28 commits

       7da61f0 Corrected slurm parameters for stringtie jobs, including stringtie_abund which was causing errors
       fb2d51d Removed DEV from sleuth configuration
       68878b1 Added note to RNAseq README regarding switch to stringtie protocol by default
       94342a5 Merge branch 'rnaseq_light_jhg' of bitbucket.org:mugqic/genpipes into rnaseq_light_jhg
       68bf886 Specified genome support for RNAseq_light in the pipeline's README
       5088cbd Corrected INI file to inlcude latest version of Mugqic_tools
       e455bbc Modified INI file to include latest version of MUGQIC_tools as well as more memory for kallisto steps by default
       e974d45 Merge branch 'rnaseq_jhg' of https://bitbucket.org/mugqic/genpipes into rnaseq_jhg
       49bbed5 Corrected minor parameters for abacus on the RNAseq ini
       330c861 corrected merge issue with the rnaseq_ligth cedar INI
       0d7a7d7 Corrected minor issues with rnaseq_light pipeline, and added documentation to the README
       160fd60 Corrected some issues with the UCSC script that was causing issues with non GRCh38 genomes
       52daaea Added documentation for stringtie protocol. Corrected some problems with the base inif file. Minor corrections to the stringtie script and rnaseq.py script.
       214930b Fixed merge issues
       497139c merged to master and resolved conficts
       95be45a Merged to latest version of master branch
       90e40ba Corrected a mistake on the methylseq ini files for Cedar and Mp2b
       5b68417 Corrected CVMFS Stringtie module in INIs
       615ad85 Fixed merge conflict with the cedar ini
       548a27f Merge branch 'master' into hgalvez
       b6a9e4e Second testing version stringtie/ballgown pipeline. Includes ballgown now, as well as corrections for stringtie. Additionally, stringtie is now default protocol for RNA-seq.
       08bbb3c Fixed a small typo with the last commit
       1abef45 Fixed bug for the order or arguments passed to kallisto
       abbbfb2 Updated version of Bioconductor used by default in the rnaseq_light pipeline (modified base.ini file).
       a97c715 Merge branch 'hector-cedar' of https://bitbucket.org/mugqic/genpipes into hector-cedar
       7a877b3 Minimally functional rnaseq_light pipeline with the added sleuth step. Still requires further testing on clusters beyond abacus. Also, some harmonization of the required genome reference files would be advisable.
       3c142b8 Merge branch 'master' into hector-cedar
       2ccccde Eliminated warning message from rnaseq_light pipeline for single read samples

  jose.hector.galvez@computationalgenomics.ca <hgalvez@abacus2.ferrier.genome.mcgill.ca>      2 commits

       15b6ecb Fixed issues with stringtie merge and abund. Ready to test in other servers
       4f3e5e8 Corrected job.py to allow for the definition of the multiqc_files parameter, should fix error with multiqc and Kallisto

  Jose Hector Galvez <Hector@rhododendron.genome.mcgill.ca>      10 commits

       4a56148 latest modifications to the tools.py script
       9363e82 Fixed minor typo on rnaseq.py
       397a903 Fixed indent issue in the stringtie job creation
       35f2904 Fixed minor bug in the stringtie job creation
       ebe3d34 Fixed minor bug in the stringtie job creation
       ffd2a51 Fixed stringtie module location in the rnaseq.base.ini file
       a76c8e4 Fixed sample name issue within stringtie.py
       c6b938f corrected minor typo
       18c6967 correct minor typo
       01cc53d First commit adding stringtie functionality, still testing

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      22 commits

       447c061 Merged in rnaseq_light_jhg (pull request #75)
       daabf76 Merged in rnaseq_jhg (pull request #76)
       de0d01e Merged master into rnaseq_jhg
       7293ec0 Merged master into rnaseq_light_jhg
       e032d8f Merged in rnaseq_jhg (pull request #69)
       fecf7c7 Merged master into rnaseq_jhg
       d687b0f Merged in rnaseq_light_jhg (pull request #68)
       e78a608 Merged master into rnaseq_light_jhg
       9225ebb Merged master into rnaseq_jhg
       3868e89 Merged master into rnaseq_jhg
       573a0e0 Merged master into rnaseq_light_jhg
       640d57c Merged master into rnaseq_jhg
       77c6491 Merged master into rnaseq_jhg
       16c4eb7 Merged in missing-cedar-ini-files (pull request #64)
       c957398 Merged master into hgalvez
       072c23d Merged hector-multiqc into hector-cedar
       753edc2 Merged master into hector-cedar
       4569486 Merged hector-cedar into hector-multiqc
       0585acc Merged master into hector-multiqc
       ff9fcfb Merged master into hector-cedar
       05daa1d Merged in hector-multiqc (pull request #38)
       eb7c1fa Merged master into hector-cedar

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      3 commits

       b897897 Added cedar and mp2b ini files for four pipelines: ampliconseq, dnaseq_high_coverage, pacbio_assembly, rnaseq_light
       289f74a modified cedar ini to support stringtie
       ac01fb0 Commit of my ini file for rnaseq_light in cedar and modifications to the pipeline to allow for bootstraps in kallisto

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      36 commits

       eb48294 update to last version of mugqic_tools 2.2.1
       4a71fcf add dnaseq_high_coverage.mp2b.ini file
       af718ac correct select_input_file call
       ea4250e correct select_input_file call
       5fbb456 correct select_input_file call
       c1fcfd2 fix haplotype list parrelization
       84001d2 fix haplotype list parrelization
       260a336 correct select_input_file call
       950c980 fix haplotype list parrelization
       bd92b74 correct select_input_file call
       6866c66 fix haplotype list parrelization
       a4427e4 fix number of haplotype and bed compatibilty
       32f78b6 fix number of haplotype and bed compatibilty
       c9bc68a correct select_input_file call
       9a4a376 remove file format error
       35b3dcb correct select_input_file call
       e6ef881 include both gatk and gatk4 lib import
       0a34bc7 correct bgzip_tabix libnrary call
       068983e correct flash log input file
       6a2b070 correct Mutect2 new arguments
       4da201a correct realigned bam name changed in inherited function
       bccc5cb add sequence_dictionary_variant proprety in dnaseq
       dbcd5a4 correct typo in DNAseq base ini
       12dbf3e remove typo in deilverable command
       dc11a6f remove verifyBamId for methylSeq pipeline
       637606e correct typo in the tabix_bgzip call
       2667523 correct a typo
       de9c8b3 add md5sum generation
       d905e59 add required=F for markDup other_option
       908ffc3 change htslib.bgzip_tabix_vcf to htslib.bgzip_tabix as thereis no _vcf function in the lib
       c98c362 make the add UMI steps skip automatically if no UMI in  the readset file
       043ba65 add known_mills option
       40e6536 add required=F for markDup other_option
       5142880 primers addition bug correction
       c859b1a Merge branch 'master' of bitbucket.org:mugqic/genpipes
       79ce228 add other_option paramter to compair in order to support other reference than GRCh37

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       c084a8e Merged in slurm_report (pull request #72)
       4258a45 Merged in new_style_class (pull request #63)

  P-O Quirion <pioliqui@gmail.com>      6 commits

       899d587 parse exit status, and only 50 lines after pro/epilogue
       a50700c log error when output logs have the wrong job number
       6b28d34 log report for slurm
       a6c5553 id setter for job
       01c96cf All class are new style, add setter to some getter
       8bf27ae working on Job

  Rola Dali <rola.dali@mail.mcgill.ca>      3 commits

       b72b06b Merged in beluga_inis (pull request #74)
       6077a58 adding draft beluga inis
       258876e README.md edited online with Bitbucket

3.1.3        Tue Dec 18 16:37:25 2018 -0500        108 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      52 commits

       e239ca0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       6cd0351 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       43c5f25 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1b3f70c Software installation scripts : new ones (mainly for Tumor SV pipeline) and updated ones Genome installation scripts : updated Homo_sapiens.GRCh37.sh (uncommented some useful commands) & Homo_sapiens.hg19.sh (updated vcf annotation process)
       8d3a636 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0a312f3 updated version of mugqic_tools
       395c7dc Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       1506eac Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       84fa517 MethylSeq - DMR - adusted 'other+options' for dmr step
       e6a306d new addings to AMP_Scanner installation script
       5e5094f Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       fd32bb5 MethylSeq pipeline - updated base.ini file with MethylKit differential analysis step requirements
       26d842c MethylSeq pipeline - added MethylKit differential analysis steps & cleaned the code
       e25b313 added methylkit.py as a bfx wrapper to call methylKit.R from mugqic_tools
       5869184 debugged 'filter_snp_cpg' & 'prepare_methylkit' subroutines used by MethylSeq pipeline
       0c798db Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       f6047c4 updated VarScan installation script with lates VarScan version v2.4.3
       efa28c8 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       869f99c updated Conpair installation script to latest Conpair version v0.2
       7e29421 added Manta installation script
       2233a27 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ad29c15 added Ruby, LAST & Picky installation scripts
       f085968 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       0683e04 updated AMP_Scanner.sh installation script
       498f190 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       41a042a add AMP_Scanner installation script
       2de2537 updated regexp in install_module.sh
       934fdec update supernova.sh with latest version 2.1.1
       6a82f17 update cellranger.sh with latest version 3.0.0
       82783db Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       36e6b1c Merge branch 'master' of bitbucket.org:mugqic/genpipes
       bad5792 CCCG-1143 - added Caenorhabditis_elegans.ce11 genome installation script
       7033838 BFXDEV-578 - updated FasTree installatino script
       c35f7b2 BFXDEV-756 - changed the '-j' pipeline parameter default value to 'slurm' instead of 'pbs'
       4847e41 added StringTie installation script
       6cc4a37 added RNA-Seq Light section in the README
       90ea223 Version bump to 3.1.3-beta
       0b4dcef Version bump to 3.1.2
       2e259ec MethylSeq - updating ini files after testings on guillimin
       19568c7 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       cbefa41 MethylSeq - cleaned methylkit_dmr subroutine
       8be6e37 MethylSeq - added MethylKit DMR steps
       2ef5d07 MethylSeq - adding the bfx wrappers for DMR steps
       35aac40 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e6d9edf Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       685e8a5 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e7f4335 fixing ini for trimmomatic16S and updating mammouth ini for dada2 protocol steps
       b716290 AmpliconSeq : resolving asva step dependencies
       57f08b8 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e8af44b AmpliconSeq - updates to Dada2 protocol
       523d1da AmpliconSeq - updated last step of dada2 protocol for better dependencies and coding standards
       8a9605e BFXDEV-674 - updates : flagstats calls are now made after alignment and deduplication instead of during the metrics step

  Édouard Henrion <henrione@cedar1.cedar.computecanada.ca>      2 commits

       7670b7e ChIP-Seq pipeline - correcting cedar.ini file for missing 'homer_make_ucsc_file' step requirements
       b66e935 AmpliconSeq - call 'zless' instead of 'less' to avoid issues on Graham and Cedar systems

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       fcafa26 Merge branch 'dada2' of bitbucket.org:mugqic/genpipes into dada2
       fa5e298 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       a92edf3 AmpliconSeq - adding cedar ini file

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      6 commits

       d6515ac Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       c3180f2 update database references in ampliconseq.base.ini
       37685bd asva.py pipeline merged to ampliconseq as another protocol asva.R deleted because moved to mugqic_tools new trimmomatic16S function added in bfx/trimmomatic.py
       5181050 update nanuq2mugqic_pipelines.py to also fetch the primer sequences ; needed by the new AmpliconSeq protocol (dada2)
       0694899 removed cutPrimer and set sys.path correctly
       dfe1513 Creating dada2 branch content

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      4 commits

       c54fe2b added pool parameter to bfx/dada2.py & pipelines/ampliconseq/ampliconseq.base.ini
       94f9800 removing deprecated pipelines/ASVA/asva.py
       9b09c66 first working verson of dada2 protocol for ampliconseq pipeline
       5acf527 minor updates (line break & indentation)

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       5627b8e GATK4 Updates - gatk4.ini - removed mugqic_dev module
       0f32dbf Updates GATK4 - dnaseq.cedar.ini - removed mugqic_dev module
       75b9153 Updatres GATK4 - dnaseq.base.ini - remove mugqic_dev module
       dfda31e MethylSeq UMI - methylseq.base.ini - updated mugqic_tools to mugqic/mugqic_tools/2.2.0 for fgbio tools
       7c9c703 MethylSeq UMI - methylseq.base.ini - updated module_fgbio to cvmfs version of the module
       b29b117 methylseq.py - commented subroutine all_sample_metrics_report as it has been remove from the pipeline (because useless)
       b8f3fb6 fgbio.py - correted typo in addumi surbroutine

  Emmanuel Gonzalez <emmanuel.gonzalez@mcgill.ca>      1 commits

       223052d Merged in dada2 (pull request #49)

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       e1c7ef8 README.md edited online with Bitbucket
       6d1bd5c README.md edited online with Bitbucket

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a6fe363 remove all_sample_metrics_report as it does the same as ihec_metrics
       b446c86 change the ini section for module in tools.methylseq_ihec_metrics_report
       598028b Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       cee0464 adding umi and on-target metrics
       96f4dec adding the missing count parameter -c to samtools.count
       ba6545f correct typos
       0d03902 correct small bugs for UMI integration
       317fe6a Merge branch 'methylSeq_UMI' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       acb5536 corect typo

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       7f369e1 Merged in methylSeq_UMI (pull request #52)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      11 commits

       1d2a7d3 add on target metrics and UMI
       86cea70 add target_cpg_profile function in metrics lib
       ad2dd57 add samtools mapped_count function
       b322c94 add samtools mapped_count function
       9d49c2a correst typo
       6ccf4ee integrate UMI step in ini file
       cf4c008 include bam UMI annotation step
       f74f1ae include UMI annotated bam to the select input file of the picard_merge_sam_files step
       5c52b3d Add other_option  system to MarkDuplicate
       2efad22 add the UMI field in the readset file
       9ef8610 add the UMI field in the readset file

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      1 commits

       cf47de1 updates to cedar ini

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      2 commits

       5b912aa Cedar resource and gatk4 metric fixes
       672ad33 Updates to cedar.ini

  Robert Eveleigh <eveleigh@ip16.m>      2 commits

       4177822 GATK4 mp2b file added and improvements to alt contig exclusions
       91fb349 Added mp2b ini and exome improvements

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       e9c66c6 Updates GATK4 and annotations

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       3c0ac31 Resolve conflicts
       71b7496 improvements to exome handling

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       ad3caf1 Merged in dnaseq_gatk4 (pull request #47)

  Rola Dali <rola.dali@mail.mcgill.ca>      1 commits

       dc7996d README.md edited online with Bitbucket

3.1.2        Wed Nov 21 15:05:01 2018 -0500        30 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      17 commits

       be396f0 updated mugqic_tools version 2.1.12
       f187cdf removed _dev modules from ini files
       2bfc307 Updated R_Bioconductor.sh : pointing to system libraries improved
       e80f729 PacBio Assembly pipeline - corrected incompatibility bug with Guillimin
       1be0124 updated install_module.sh to accomodate both CentOS7 & Ubuntu16.04 system libraries
       ff709f1 update HiCUP install script with latest version 0.7.0
       659fbc9 corrected typo in ampliconseq.py causing pipeline crach at plot heatmap
       c8b5634 updated README-release.txt with correct URL for GenPipes download page
       09fe369 Changing resources requirements for [gatk_merge_and_call_individual_gvcfs] in DNA=-Seq pipeline
       30a9d33 corrected reference to kallisto tx2gene file in the RNA-Seq light pipeline
       0e7dae6 corrected reference to kallisto index in the RNA-Seq light pipeline
       e8ec0ac Corrected bug in HicSeq pipeline when running in capture mode
       58841e5 updated RNASeq light base ini file with kallisto version on CVMFS
       78a4d4b corrected core/pipeline.py to avoid pipeline erroring when using '--json' parameter...
       1e2ef94 corrected typo in R_Bioconductor.sh
       e4f493f Version bump to 3.1.2-beta
       c6b48be Version bump to 3.1.1

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      4 commits

       929e37c Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f4479be freebayes install script
       9fadba9 platypus install script
       7909a4e vcfanno fgbio & delly install scripts

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       362f101 Merged in hicup_arima (pull request #51)
       f3593a6 tumor_pair.base.ini : changed gemini version to 0.20.1
       31d2d00 tumor_pair.guillimin.ini : removed some more _dev modules
       b689763 tumor_pair.guillimin.ini : removed mugqic_dev module
       6e9b3f9 hicseq.base.ini : updated HiCUP version to 0.7.0
       07b9a21 dnaseq_high_coverage.base.ini : setting ram parameter witin section igvtools_compute_tdf
       d210893 Tumor_pair pipeline : bug fix in bfx/bcbio_variation_recall.py Correted typo in the executable call

  Rola Dali <rola.dali@mail.mcgill.ca>      2 commits

       6d94474 editing hicseq.py for Arima compatibility
       33d597d adding HiC Arima digest to install_genome

3.1.1        Thu Nov 1 15:32:25 2018 -0400        161 commits

  David Bujold <david.bujold@mail.mcgill.ca>      1 commits

       b0adf94 Merged in pipeline_stats (pull request #20)

  dbujold <david.bujold@mail.mcgill.ca>      2 commits

       86c8a72 Display JSON log statistics into tables and figures on the log VM.
       6f716cb Python CGI script to tranform pipelines stats log file into a JSON document.

  Edouard Henrion <edouard.henrion@mcgill.ca>      32 commits

       2495721 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b14189d update install_module.sh script with integration of apt along yum as system libraries resources
       0de17dd update R installation script with integration of apt along yum as system libraries resources
       7a57dfe Merge branch 'master' of bitbucket.org:mugqic/genpipes
       de73f87 update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       cd6f62a update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       052e9fc Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e8af086 updating genome ini file generation with versioning
       e2f53d8 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b929742 update smrtlink.sh with latest SMRTLink version
       8130d05 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       13080e2 updated R_Bioconductor.sh with new packages and corrections
       63e12b5 new version of multiqc.sh to install latest version
       f499c1b aded popoolation2 installation script
       211dc62 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       d58b1c3 updated gatk.sh so it can handle both versions 3 & 4 installation
       4b5565f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       c4004ed updated R_Bioconductor.sh - new libraries & patches
       38e01cd updated python.sh : make -j12 & configure command
       b06c410 resources/modules/install_module.sh
       d1fe8bf update flash.sh with parallele make -j12
       6eca5f3 watch_portal_folder: removed sample_name in filename
       f0e6222 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1ef165c corrected gatk.sh : paths in the modulefile were wrong
       749cd97 updated rnaseq.py : report steps are now included in the analysis JSON file
       8978a59 updated ampiconseq.py : both steps merge_flash_stats & merge_uchime_stats are now inclded in the analysis JSON file
       b4d4b04 updated rnaseq_light.py : samples added to each job including report jobs, reviewed indentations and spacing
       8e75e02 updated rnaseq_denovo_assembly.py : samples added to each job (some were still remaining) including report jobs
       61af5bf updated common.py : added samples to the call of rmarkdown.render
       cc974a3 updated bfx rmarkdown wrapper : JSON file generation added
       1d9b834 Version bump to 3.1.1-beta
       bd721f1 Version bump to 3.1.0

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      1 commits

       b3e4560 removing all the remaining MUGQIC_INSTALL_HOME_DEV in all the ini files

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      6 commits

       4a4d0e5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ca805a0 updated install_module.sh : new installation procedure integration, adding the patching the C3G executatbles with patchelf : making sure all system libraries are now searched in /cvmfs/soft.mugqic/yum/centos7/1.0
       f771e9b updated python.sh with some new package : umap-learn
       63a855e udpated R_Bioconductor.sh with new packages. Also, now integrates the new installation procedures including pathcing the C3G executables (i.e. use of patchelf)
       41fcad8 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       4a5958c updated archive URL in flash.sh install script

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      7 commits

       c243ec8 updated R_Bioconductor.sh script : set the PAGER variable to /usr/bin/less
       f315771 updated pipeline READMEs : all the steps are now shown independently of the available pipeline protocols
       5fa7bb4 updated R_Bioconductor.sh with some new packages in the install list
       1c468e0 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b24d7f8 updated silva.sh to install latest vesion of silva DB
       04140e3 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e9fcbff added Glycine_max.v2.sh installatino script for Glycine (Soybean) genome installation

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      16 commits

       bdf93cc updated hicup.sh
       795b9e5 updated smrtlink.sh with the version 6.0.0 of SMRTLink
       0a9fb4b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e141ec7 updated pipeline cedar ini files to avoid using java from CC software stack
       efbf07a modified core/pipeline.py to avoid having to set  when not generating the anaylsis JSON file
       256907f modified DNA-Seq README with better step descriptions
       76a26cb Updated the pipeline workflow diagram download links : path of the full-size picture instead of the resized one
       d0fd874 Added a link to download the pipeline workflow diagram along with the diagram picture itself
       11cbf4a Updated pipeline README.md files, with workflow diagram pictures embeded
       30c53cc deleted pipelines/tumor_pair.base.exome.ini
       bd6f000 moving tumor_pair.base.exome.ini from 'pipelines/' to 'pipelines/tumor_pair/'
       9362d72 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       baf53ba updated ampliconseq pipeline following Emmanuel gonzalez comments
       6a85b8f added source EnsemblPlants to install_genome.sh script
       690f20b added skewer installation script
       59996b5 BFXDEV-591 - updated ampliconseq.base.ini based on Emmanuel Gonzalez feedback

  ehenrion <edouard.henrion@computationalgenomics.ca>      2 commits

       3b3ff68 updated mugqic_pipelines.sh so that it now refers to genpipes repository on bitbucket
       97a4552 corrected typo within gatk.sh installation script

  ehenrion <edouard.henrion@mcgill.ca>      18 commits

       1238882 methylseq.py : corrected dependencies assignment for ihect_sample_metrics_report step
       85bfc18 rnaseq_light.py : typo correction
       2057f1e rnaseq_light.py : corrected typo
       91be7b8 rnaseq_light.py : corrected typo in job parameter assignement
       5de9c4e dnaseq.cedar.ini removed 'module_java=java/1.8.0_121' from cedar.ini
       b4ce9ef methylseq.base.ini : modified [bismark_align] section : maximum insert_size now set to 1000
       146eaa7 methylseq.mammouth.ini : modified bismark align section
       a4eea15 methylseq.cedar.ini : modified bismark_align walltime
       a22490f methylseq.base.ini modified [bismark_methyl_call] section within the base.ini
       a5276a3 methylseq.base.ini - modified bismark_align parameters and resources within the base.ini file
       4fd03a1 dnaseq.py edited online with Bitbucket Correct typo at line 724 : "jobs" replaces "obs"
       425c086 dnaseq.cedar.ini : added variant_recalibrator section to defined resources on cedar
       d55362f Merged in revert-pr-41 (pull request #46)
       17aa18f Revert "Can run on HPC with slurm and  containers (pull request #41)"
       60be1aa rnaseq.py edited online with Bitbucket added missing job (metrics.wigzip) to the json analysis file
       1c5f56d README.md edited online with Bitbucket
       6dccf5b README.md edited online with Bitbucket
       3b3407b README.md edited online with Bitbucket

  Éloi Mercier <emercier@cedar5.cedar.computecanada.ca>      1 commits

       c6593cf in cedar.ini and graham ini of rnaseq, chipseq and dnaseq: change assembly_dir to MUGQIC_INSTALL; in dnaseq.graham.ini: uncomment assembly_dir variable

  emercier <eloi.mercier@mcgill.ca>      10 commits

       ec3183f all pyc files removed
       f0568f4 Update R_module in ini files to mugqic/R_Bioconductor/3.5.0_3.7 (except for illumina_run_processing)
       5a14b4d in ampliconseq, pacbio and rnaseq.guillimin.ini: change lm queue (depreciated) to meta queue
       e13be16 in modules/weblogo.sh: remove whitespaces at the beginning of the echo blocks
       4675ad6 in dnaseq and rnaseq.base.ini: change R_Bioconductor to stable version 3.4.3_3.6
       fcad025 in install_all_genome.sh: add Danio_rerio.GRCz11.sh
       d353ecf add install script for genome zebrafish Danio_rerio.GRCz11
       9416156 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0b53f64 in install_genome.sh: small fix to create_kallisto_index and create_transcripts2genes_file functions
       f5681c4 in install_genome.sh: fix bug in create_transcripts2genes_file to work with recent version of Ensembl

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      4 commits

       53d11f3 rnaseq.base.ini edited online with Bitbucket: commented explicit adapter file parameter
       8653d02 R_Bioconductor.sh edited online with Bitbucket: added PopSV. Currently commented out since not tested
       0b42f96 R_Bioconductor.sh edited online with Bitbucket: Added a few dependencies to the list
       3a268e6 R_Bioconductor.sh edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       a23ba48 Added dev install scripts for delly, lumpy, sv, vcfanno

  José Héctor Gálvez López <hgalvez@ip16.m>      2 commits

       88e4929 Further refinements to Mp2b based on feedback from mammouth admins
       39bf760 Added ini files tailored for Mp2b based on Cedar ini files for the following pipelines: RNA-seq, RNA-seq de novo, DNA-seq, Hi-C seq, Methylseq, and ChIP-seq.

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       1f95d68 Add cedar ini file for methylSeq
       8b07176 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b1ed77c removing pipe empty of jobs in tumor_pair
       7e70dc5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       db20957 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9e0178b correct bug in ucsc bedGraphToBigWig which raise an error for specie with MT chromosome name except for GRCh37 - BFXDEV-737

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      7 commits

       2e021f6 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       85ead2b resolve issue with multiple inputs in DNAseq picard markduplicates
       1e00830 resolve pull conflict
       b0042a0 update dnaseq.cedar.in file
       4074ab5 adjust dnaseqto add select input file to some of the steps - need to be continued
       e38920e modify Slurm scheduler delay (sleep) from 0.5 to 0.2
       b1d3d4b ChipSeq - add mutliqc param in the cedar ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       9c7c815 Merged in add_container (pull request #41)

  P-O Quirion <pioliqui@gmail.com>      5 commits

       cfa63b2 Merge branch 'master' into add_container
       44274bd Basic container ini file
       4296d31 Working singularity version tested on graham
       b14afd9 WIP exec line prototype
       8b92416 add container option

  Rola Dali <rola.dali@mail.mcgill.ca>      8 commits

       795776b methylseq.cedar.ini edited online with Bitbucket: added cluster_walltime to ini to avoid errors
       47addd2 chipseq.base.ini edited online with Bitbucket: bigwig and run_spp resources are not enough; edited them to avoid failure
       50747a7 Merged in IHEC_metrics (pull request #44)
       bdb0b76 dnaseq.base.ini edited online with Bitbucket: change threads for GATK due to errors
       142921b update csvToreadset.R in utils to use column names due to changes in nanuq csv
       40cdff0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       8709b9d README.md edited online with Bitbucket
       13f08a5 README.md edited online with Bitbucket: -j slurm and tutorial

  Romain Grégoire <romgrk.cc@gmail.com>      3 commits

       6f8ea86 Merged in add-json-version (pull request #50)
       1271bc1 Merged in dashboard-display-launching-user (pull request #43)
       73e768a Merged in logs-add-checksum (pull request #42)

  Rom Grk <romgrk.cc@gmail.com>      28 commits

       ff11dbc jsonator.py: add version number
       bb72124 watch_portal_folder: read sample_name from filename but stay backward-compatible
       921e54d job2json: add sample_name to file for genpipes dashboard
       784848f copy sample json files during script run
       b414b04 scheduler.py: clean unused value
       5a9d8af Sample: remove .json_dump property
       7d45ce1 lint pipeline.py
       337a1c8 lint job2json
       014c19d job2json: guard __main__
       26b9ef9 lint job2json
       bf4f976 job2json: use $USER of user running script
       8980701 common.py: remove unused import
       622b611 common.py: fix log issues
       676f7be common.py: add unique md5 checksum to logs
       56545fe Revert "watch_portal_folder: put sample_name in filename"
       a8badcd Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0526ecc watch_portal_folder: put sample_name in filename
       13b5638 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0411ae9 watch_portal_folder: fix memory errors
       e5e3579 watch_portal_folder: add logging
       46d8a35 watch_portal_folder: fix typo
       f54a127 watch_portal_folder: add cache option
       7f89e82 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       2fbdca8 watch_portal_folder: implement update by diff
       45cc189 watch_portal_folder: update script
       67c1471 watch_portal_folder: skip missing files
       8a9dd08 watch_portal_folder: more resilient to network errors
       32bf92a watch_portal_folder.py: dont watch if no interval is provided

3.1.0        Wed Mar 28 15:46:33 2018 -0400        188 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      77 commits

       ff52e7c MethylSeq - added 'ram' parameter to the igvtools_compute_tdf step in methylseq.base.ini
       68661e8 updated jb2json.py with a better locking system : now creates a folder instead of a file in order to create the lock
       55bf1b0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       43c9cc1 updated core/scheduler.py to remove job2json call when --json parameter is omitted
       21b3bdf updated gatk.sh with version gatk-4.0.2.1
       9908731 updated longranger.sh with LongRanger version 2.2.2
       d03abd5 updated hicseq.py : calling of the newly built bfx libraries and minor indentations changes
       c52b582 reviewed topdom.py wrapper for a better handling of input and putput files, so that job dependencies do not break...
       8ae89bd aded locking file system to jsonator.py and job2json.py to avoid multiple and synchronous writing attempts to JSON file
       4d1eaa6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       70ddf87 updated hicseq.base.ini with newer version of mugqic_tools (2.1.11) and with revised resource allocation for hic_align step
       3f79827 removed useless loading of mugqic_tools module when calling job2json.py
       4d7da45 added locking file system to job2json to avoid multiple and synchronous writing attempts to the JSON file
       6aafbf0 removed useless loading of mugqic_tools module when calling job2json.py
       205d577 mugqic_tools.sh : swith to mugqic_tools version 2.1.11
       0ea6336 improved the process regarding the update of the analysis JSON file when resuming a pipeline execution
       62dd71d removed useless import from core/scheduler.py
       ac6e9ce added the locking file system to avoid multiple & simultaneous writing attempts on the same file
       ffd35a6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       e654508 added new bash tools wrappers to be called by hicseq.py pipeline script + some minor indentation & syntaxe updates
       843d657 added new bfx libraries to be called by hicseq.py pipeline script
       1f6be1b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e1dd68f added some new software installation scripts
       e2abf45 updated 'core/job.py' : added 'samples' parameter to 'concat_jobs' and 'pipe_jobs'
       8cf84f4 added the initialisation of 'portal_output_dir' to '/lb/project/mugqic/analyste_dev/portal_output_dir' within all of the pipeline .base.ini files
       9bae17c added the generation of the analysis JSON file to the SLURM scheduler class
       6e6363d Merge branch 'master' of bitbucket.org:mugqic/genpipes into cedar
       87fb014 updated core/scheduler.py : call of job2json.py is now done without /home/ehenrion/work/portal_repo/genpipes
       1681c6f Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       f95e690 corrected trimmomatic input selection to reflect the file names what might be outputed from picard_sam_to_fastq
       f22ee01 update verifyBamID subroutine within common.py so that it now accepts .vcf.gz files (formerly was only accepting .vcf files)
       998e1bf updated version of cellranger to 2.1.1 within the bash installation script
       ba5099f updated pipeline READMEs
       b5cecbd updated core/pipeline.py to make the analysis JSON file creation optional : no JSON file created by default
       3fae6e8 updated bfx/macs2.py so that it follows the C3G coding standards
       3b426ea added installation script for hdf5 and zlib libraries
       2e82f07 corrected typo in install_genome.sh
       ff59b0a updated install_genome.sh script with newer software version as well as minor indentation fixes
       3d26798 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1bd0b98 updated kallisto installation script
       cf3219a updated STAR version in the installation script
       19997fa updated dbSNP version to 150 for Homo_sapiens GRCh37 installation script
       1e6f032 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f89e056 updated kallisto version to 0.44.0 in the genome installation script - corrected typo within create_transcripts2genes_file subroutine
       7bf13a9 updated dbSNP and dbNSFP in Homo_sapiens genome installation scipts
       30299e2 added cDNA fasta to the genome installation scripts
       e17fcb1 added the portal_ouput_dir to in base ini file
       acea9bf Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       fe61627 Merge branch 'portal-integration' of bitbucket.org:mugqic/genpipes into portal-integration
       e0010e2 Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       1f63dc5 Added README-GenAP_coding_standards.txt which gives the coding guidelines for who may want to participate in the C3G/GenAP developments
       ed345d2 BFXDEV-721 - updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       9a59536 updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       0ed0ef6 updated Prokka installation script so that 'prokka --setupdb' is launched afer the installation
       743668f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       22b2601 updated Supernova bash installation script to version 2.0.0
       c62cf42 updated the bash installation script of python : installation of QIIME has been updated so that emperor is well installed before QIIME
       96c98a4 added the initialization of the LFS variable, otherwise it is recognize on guillimin...
       b02491e updated the bash installation script of python : after the installation of the specified version of Python is successfull, the script now also takes care of the installation of all the needed python libraries
       7609147 update type within python_lib.sh script
       1ab7a7e RNASeq pipeline : value of parameter 'localfit' used by DESeq, is set to default i.e. empty (which also means 'false'), so Parametric Dispertion Fit is performed by default instead of Local Fit
       896c3b2 DNASeq pipeline : updated versions of mugqic_tools (to 2.1.10) and of samtools (to 1.4.1) within the base.ini file
       74a5c9a updated picard installation script : now installs version 2.17.3
       8d2dde0 added the installation script of Prokka, a tool for prokaryotic genome annotation
       cbfe9d7 BFXDEV-673 - remove some verbose during the execution of job2json.py
       66893f8 updated mugqic_tools.sh so it now installs version 2.1.10
       170f39b added the bash installation script for the Illumina InterOp parser
       dc0ebcb BFXDEV-674 - MethylSeq pipeline - minor updates within the metrics .md report template file, for standardization purpose
       d9e7b0b updated ortograph within sample metrics .md report files
       e0b41ff BFXDEV-673 - updated scheduler.py to generalize the use of job2json to all schedulers
       5250d80 BFXDEV-674 - added a report .md file for IHEC metrics reports for targeted anaylsis
       cb466e7 BFXDEV-674 - updated .md files for metrics reporting : revised headers & descriptions
       460f24a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       2833f6a updated gemini.sh bash installation script : added PYTHONHOME setting when loading the gemini modules
       f483e37 update release instructions to generate proper README for all the pipelines including HicSeq
       2d70be4 Version bump to 3.0.1-beta
       3cb8610 Version bump to 3.0.0 - updated

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      2 commits

       8983180 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       49c6b32 updated kallisto bash installation script

  ehenrion <edouard.henrion@computationalgenomics.ca>      3 commits

       3f60795 updated R_Bioconductor.sh with good indentation and new packages installation
       b1215c3 corrected wrong environment variable name within blast.sh
       dc4fbb7 updated the version to v359 within ucsc.sh script

  ehenrion <edouard.henrion@mcgill.ca>      6 commits

       39fea26 rnaseq_denovo_assembly.base.ini edited online with Bitbucket
       f9278f9 README.md edited : "MUGQIC pipelines" replaced by "GenPipes"
       e6abfd9 apis_mellifera.sh deleted : was the exact replicate of Apis_mellifera.sh
       8e9ab1c README.md edited Updated links to the pipeline pages
       469b9b4 README.md edited online with Bitbucket updated some links to reflect the repository renaming to genpipes
       87524a9 smrtanalysis.py - standardized command-line format and indentation

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      3 commits

       99b1a5d removing pyc files from rnaseq
       1c7731a in rnaseq.mammouth.ini add section for bed_graph to set up ppn to 1
       6d91088 In nanuq2mugqic: change readset file name to readset_<seq_type>.tsv

  Eric Fournier <ericfournier2@yahoo.ca>      2 commits

       6fafc8d Merged in ericfournier2/genpipes (pull request #35)
       4138a72 Fix list within list bug which breaks chipseq pipeline.

  Gary Leveque <gary.leveque@gmail.com>      1 commits

       6857307 Merged in gary_pacbio_assembly (pull request #29)

  gary.leveque@mail.mcgill.ca <gleveque@abacus2.ferrier.genome.mcgill.ca>      7 commits

       fd1e91c additions made to smrtanalysis.py for basemodification and motifMaker steps
       b88dd9e Addressed issues commented by Edouard; tested on abacus and mammouth
       39f57bc added cluster_server= to pacbio_assembly.mammouth.ini
       b135ae0 revised versions of .base and mammouth.ini files; I was changing between sw and lm nodes
       6bac368 revision of pacbio_assembly.mammouth.ini, back to qwork
       356aadc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into gary_pacbio_assembly
       a24451f Addition of base modification detection and generation of a motif_summary.csv steps to the pacbio HGAP assembly pipeline;  see BFXDEV-703

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       2455b6d Merged in hector-cedar (pull request #37)

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      1 commits

       36b892f Corrected minor bug in the create_scheduler() function that was creating errors when using slurm

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       05bbcec Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3b818dc externalize excluded chormosome to be specified in the genome ini
       a2d807b remove conflict
       3896d77 uniformize verify_bam_id vcf path
       f4bf8cc add a delay (sleep) after the symlink creation to avoid issue (Invalid Job Dependency) during submission
       322348d Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f734dfb test new url
       3752940 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0834b9 add NovaSeq to nanuq2mugqic script
       0dd2ba3 DNAseq - split the pipeline into 2 different pipeline GATK best_practice or old mpileup - BFXDEV-716
       5c586ff Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9d48a58 chnage RNAseq_denovo inheritance to RNAseq
       8a250c8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7716f96 include a test of inputs type for some library function which does not enforce the list type

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       804ae6c Merged in cedar (pull request #34)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      24 commits

       e2b2a7a resolve conflict with master
       82d1bee Merge branch 'cedar' of bitbucket.org:mugqic/genpipes into cedar
       443fcf6 restore low memory ini now the slurm bug is corrected
       78d156b update rnaseq for cedar
       62520da update rnaseqDN on cedar
       9a27402 update hicseq on cedar
       b9842b9 update dnaseq ini on cedar
       b695bcf add more more RAM to compensate slurm bug
       924fff5 add 0.5s sleep to let slurm submiting the job correctly
       282fa9e Increase memory request to avoid I/O buffering hit the memory limit
       5f157cd Changing I/O block size
       3b395c5 Changing I/O block size
       614564d remove confilct pulling master
       be81579 Modifying I/O block size for cedar
       496e10d update ini file
       7f77b19 update log with slurm scheduler
       8180609 add chipseq ini file for cedar
       0e27909 update RNAseq and link to mugqic_dev
       d26fef6 update DNAseq and link to mugqic_dev
       10c2069 remove lattency for jobsubimission to slurm
       d192a91 remove conflict
       77e5d35 update cedar scheduler
       67b8006 update cedar ini
       120f92d DNAseq -ini file for cedar

  Mathieu Bourgey <mbourgey@cedar5.cedar.computecanada.ca>      9 commits

       9a9ab77 DNAseq - update ini file
       910be76 remove scheduler module testing
       9e416b1 change igvtools exe to igvtools jar in order to have control of the ram usage
       9d3c7bd adjust RNA and DNA ini files;  generate fake prologue and epilogue; add 2s delay between each job submission- BFXDEV-683
       95f20c6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       f247a6a create RNAseq ini file - BFXDEV-683
       6f5cca5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       297a2c6 create a DNAseq ini for cedar; remove typo in base DNaseq ini  - BFXDEV-683
       1dce04b create a version of the scheuler class for slurm  - BFXDEV-683

  Mathieu Bourgey <mbourgey@gra-login4.graham.sharcnet>      1 commits

       fea1501 Add Graham ini for RNA and DNA

  Rola Dali <rola.dali@mail.mcgill.ca>      17 commits

       6978e14 Merged in IHEC_metrics (pull request #33)
       d66737b added .hic file generation to capture hic
       3a9b9a3 adding RobusTAD TAD scoring to hicseq.py
       3c79721 adding multiqc to chipseq. should customize yaml file and check homer module
       ad0ff99 changing rnaseq dependencies to ensure no repeats when job is complete
       e581c3a adding mugqicValidator.py to utils to validate basic structure of readset and design file
       8ee1083 removing job outputs from chipseq ihec metric method to accomodate samples running without input
       befab29 Merged in IHEC_metrics (pull request #31)
       2a0d785 fixing the tmp_dir for macs
       6a363de Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       a10ad51 changing genome to fit with merged homer.maketagdir and changing ihec matrics output name
       bba10e0 Merged in IHEC_metrics (pull request #28)
       28850c5 config.py edited online with Bitbucket
       00b9e34 adding query_module=spider to mammouth ini
       00a746f allowing module spider on mammouth to reduce module loading time--committing to test on guillimin and abacus
       e3c9fab Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fea082b fixing MACS bug BFXDEV-712

  romain.gregoire@mcgill.ca <rgregoir@abacus2.ferrier.genome.mcgill.ca>      2 commits

       53fe3b9 pass config files directly to job2json.py
       94d21ed Merge branch 'portal-integration' of https://bitbucket.org/mugqic/genpipes into portal-integration

  Romain Gregoire <romgrk.cc@gmail.com>      1 commits

       655f7c9 watch_portal_folder.py: fix --username argument

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       2bcae74 Merged in portal-integration (pull request #32)

  Rom Grk <romgrk.cc@gmail.com>      15 commits

       c3a0080 watch_portal_folder.py: remove username option
       f3c3f34 send $USER to portal integration
       bdb24f8 watch_portal_folder.py: fix file path
       4981525 watch_portal_folder.py: sort files sent by modification time
       e1f361a export CONFIG_FILES to allow loading config from pipeline
       430d0c4 export CONFIG_FILES to allow loading config from pipeline
       446122c job2json.py: add mugqic dir to python path
       c35b6f9 watch_portal_folder.py: safer handling of response
       13196bf watch_portal_folder.py: fix data posting
       bbbfe8a watch_portal_folder.py: fix script exit
       bff9c64 watch_portal_folder.py: fix argument passing
       47b1073 use uuid to avoid collisions when buffering JSON files
       7f99bdc job2json.py: make a copy of the updated JSON files for the portal
       629cef6 add watch_portal_folder.py
       63e4de6 jsonator.py: make a copy of the JSON files to a buffer folder

3.0.0        Thu Dec 7 14:19:49 2017 -0500        444 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      247 commits

       8c5a5c6 MethylSeq pipeline - BFXDEV-674 - updated wiggle_tracks step with more comprehensive output file names
       2024838 updated ucsc.py with simplified if-else statement for more clarity and corrected the 'chr' prefixing behavior
       bcf3d8a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       8594288 updated Homo_sapiens.GRCh38.sh installation script with vcf indexes
       cb43d49 updated jellyfish installation script with version 2.2.3
       9e3c955 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7e14f09 Updated README for all the pipelines
       a6b244e updated mugqic_pipelines.sh script with the most recent version of the mugqic_pipelines, now called GenAP_Pipes, version 3.0.0
       86092c6 Version bump to 3.0.1-beta
       b8f4310 Version bump to 3.0.0
       3844eee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c9ad368 added tagInfo.txt in homer.py
       3370260 DNA-Seq pipeline - update base.ini file by removing all the paths which were still pointing to a '_dev' location
       7b15520 BFXDEV-674 - MethylSeq pipeline - updated bedtools.py intersect function to carry header from input bed to output bed
       da1ece6 added tagInfo.txt
       1b84ec6 Version bump to 3.1.0-beta
       03f3dda Version bump to 3.0.0
       5d0a774 added the README file for RNA-Seq Light Pipeline
       837b4fa slightly updated release instructions
       1343ab6 Version bump to 3.0.1-beta
       188037f Version bump to 3.1.0-beta
       03638bd Version bump to 3.0.0
       4c027b6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8768a77 updated version to 2.0.2 for cellranger in the bash installation script
       00eb96b RNA-Seq pipeline - updated .base.ini file : removed specific module loading for differential expression step + minor indentation updates within trinity.py for RNA-Seq de-novo Assembly pipeline
       a4ab697 RNA-Seq pipeline - update .base.ini file with newer software versions
       626856f DNA-Seq pipeline - updated base.ini with GATK version 3.7
       a7b56ba updated base.ini with GATK version 3.7
       dde4f14 master VERSION changed to 3.0.0-beta
       bfdcbae PacBio Assembly Pipeline - updated circlator step : created a circlator.py within the bfx folder and review the coding/calling of the circlator step within the pipeline python wrapper
       132f2f2 BFXDEV-674 - updated MethylSeq pipeline to correct wrong dependencies causing some jobs to always be consider as 'NOT up to date'
       8cbbff9 BFXDEV-673 - corrected error in bfx/jsonator.py occuring when modifying the list of softwares from an already existing JSON file
       9428046 some more minor updates on bfx/tools.py regarding standard indentation and parameters naming
       1aaf4f9 MethylSeq - IHEC metric report jobs are now labelled with the sample name
       29e74ac minor updates on bfx/tools.py especially to make indentation uniform across the whole file
       573dba5 updated & added software install scripts
       c8fa5cd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d964b2 commit prior merging to avoid conflict
       2baf982 AmpliconSeq pipeline - adding mammouth walltime for qiime_otu_assigning step
       fb8c75d resolving conflicts in bfx/tools.py
       fdddddb BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq pipeline
       c961e15 corrected typo introduced after resolving conflicts...
       f4116b0 BFXDEV-673 - updated jsonator.py to handle Illumina Run Processing ini file entries when generating JSON analysis file for Illumina Run Processing pipeline
       d573722 BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       fb346c9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       763aac0 BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       8b40129 BFXDEV-673 - adding analysis JSON file generation to the Tumor Pair pipeline
       9653250 BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq De Novo Assembly pipeline
       f961ef1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1dcae3f BFXDEV-673 - added the analysis JSON file generation to the RnaSeq pipeline
       f107e04 BFXDEV-673 - added the analysis JSON file generation to the PacBio Assembly pipeline
       b806ac6 BFXDEV-673 - updated jsonator.py to handle PacBio Assembly readset files when generating JSON analysis file for PacBio Assembly pipeline
       fdd8577 BFXDEV-673 - added JSON analysis file generation to the DnaSeq high Coverage Pipeline
       d48d642 BFXDEV-673 - adding JSON analysis file generation to DnaSeq pipeline
       b8824ef BFXDEV-673 - adding JSON analysis file generation ito ChipSeq pipeline
       d19b304 BFXDEV-673 - updated AmpliconSeq pipeline with analysis JSON file generation
       102f66c BFXDEV-673 - updating jsonator.py to generalized the way dbsnp_version and server entries are handled
       054e410 BFXDEV-674 - MethylSeq pipeline - adding samtools to the loaded modules for methylseq_metrics_report
       0169e50 BFXDEV-673 - minor update to the help content
       18e8e68 BFXDEV-674 - updated wiggle tracks generation tools by spliting bedgraph and wiggle traks into 2 different jobs, also managing .bw file for GRCh37 build to make it UCSC compatible
       cc7af63 updated 'tmp_dir' in guillimin .ini files : now using  space which is automatically cleaned up after the job ends
       04d9fed BFXDEV-674 - updated methylseq_ihec_metrics_report to reflect changes in mugqic_tools IHEC_methylseq_metrics.sh
       6e854c2 Picard mark_duplicate minor update
       8958fae BFXDEV-673 - corrected error in variable assignment within jsonator.py
       4872276 BFXDEV-673 - jasonataor.py now handles the case where user has used the 'force' argument for his pipeline execution in order to entirely re-create the analysis JSON file
       22c477d BFXDEV-673 - review job2json.py passed parameters to handle start & end job date/time
       61225ee BFXDEV-673 - JSON analysis file key 'hpc_center' changed to 'server'
       4b2da7f BFXDEV-674 - MethylSeq pipline - updated mugqic_tools from dev to prod within pipelines/methylseq/methylseq.base.ini
       3e8afc3 BFXDEV-673 - updated core/scheduler.py to handle both job_start_date and jobs_end_date for the json analysis file
       c777a64 BFXDEV-673 - updated core/pipeline.py for a better generation of the json analysis file
       1d735a0 BFXDEV-673 - updated bfx/jsonator.py
       1b0b956 BFXDEV-673 - updated utils/job2json.py
       98b0d7f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       625ca4f added optparse library installation in R_Bioconductor.sh
       f751eb5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c41c514 removing MULTIQC from mugqic tools and updating version to 2.1.9
       00118d6 version 1.0.5 of R_mugqic_packages
       18ab851 BFXDEV-674 - removed useless bedtools parameters from mammouth and guillimin ini files
       a17600f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       b98dcf9 adding some new bash install scripts
       89a48fd add chicago module
       0c6d44a fix conflicts
       9919d1a keep locale module change on libuser version of the repo
       6a3eff8 keep locale module change on libuser version of the repo
       2b55b3f BFXDEV-674 - corrected input flagstat file name for metrics steps
       9ac6e81 BFXDEV-674 - adding on_target context handling
       69eb2bb BFXDEV-673 - corrected log file assignment to the JSON file
       7dd7c0d BFXDEV-673 - BFXDEV-674 - updated syntax and added IHEC metrics report
       5fccb8b importing samtools from common.py
       18e87a5 BFXDEV-668 - added ihec_metrics step to rnaseq.py
       167e9c3 BFXDEV-674 - removed useless parameter settings from ini file
       f7e2c94 BFXDEV-673 - updated cluster_hpc_center parameter within the ini files
       5ca6c52 removed unnecessary trailing comma dnaseq.py
       62c711a removed unnecessary parameters from dnaseq.mammouth.ini
       e3186d4 BFXDEV-674 - corrected call of the bash script which generate the metrisc for IHEC
       10c042d BFXDEV-673 - updated jsonator for a better behavior of json file updates
       bfc0042 updated syntax for unrequired paramaters in bedtools.py
       8e42663 BFXDEV-674 - updated bedtools with proper syntax standards
       f9538e8 BFXDEV-668 - BFXDEV-675 - updated bedtools graph other_options parameter to fit with group syntax standards
       790234e BFXDEV-668 - BFXDEV-675 - corrected typo in ucsc.py and bedtools.py
       af16960 BFXDEV-674 - adjusted walltime for bissnp step
       ba015ef BFXDEV-673 - corrected error in job2json.py
       3f5c34e BFXDEV-674 - corrected typo in ucsc.py
       47cd3f0 BFXDEV-674 - updated guillimin.ini
       fb5d03c BFXDEV674 - updated base ini
       9100666 unset batch.ini file
       4b8594c BFXDEV-674 - briaree ini file
       14ae3c1 BFXDEV-674 - mammouth ini file
       58f6b28 BFXDEV-673 - updated some 'if' statements to avoid syntax errors...
       f036046 BFXDEV-674 - updated call of ucsc.bedGraphToBigWig within the pipeline wrapper
       8e187ab BFXDEV-674 - added walltime to bismark align step
       a2bcb0c BFXDEV-674 - updated bedGraphToBigWig subroutine in ucsc.py, added bedToBigBed to ucsc.py, and updated minor things in bedtools.py
       0506547 BFXDEV-668 - BFXDEV-675 - updated bedtools.py graph subroutine which nows calls ucsc.py to avoid code redundancy of bedGraphToBigWig
       4a6e9cd BFXDEV-668 - BFXDEV-675 - updated ucsc.py to handle cases where bedGraph is in .gz format
       895d820 BFXDEV-674 - type correction in job2json.py
       a322705 BFXDEV-674 - cancel the creation of one folder per sample for the json files as on file per sample is enough
       12a085c BFXDEV-674 - updated json file path in scheduler.py
       c2f2827 BFXDEV-674 - updated bismark align & dedup outputs in order to create good dependencies for all_sample_metrics_report and ihec_sample_metrics_report steps
       607ee22 BFXDEV-674 - changed the json file location to a more simple one and got rid of the resume subroutine since not used anymore
       5b6f982 BFXDEV-674 - corrected error in scheduler when launching job2json command
       25a4ede BFXDEV-674 - another typo correted in core/scheduler.py and removed the dev references from pipelines/methylseq/methylseq.base.ini
       d5aa74d BFXDEV-674 - correcting typo in core/scheduler.py during Json generation
       1dd7bd4 BFXDEV-674 - added some missig report files and updated job2json file command and tool
       c8266f3 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       da1de3c BFXDEV-674 - updated pipeline wrapper with report files, ihec metrics and updated metrics computing
       c59a302 BFXDEV-674 - added HPC center name to the ini files
       a128e70 BFXDEV-674 - updated bedtools when working on MethylSeq pipeline to a better handling of job piping from bamtobed to coverage
       57b0229 MethylSeq pipeline - updated pipeline python wrapper - BFXDEV-674
       807a26b BFXDEV-673 - add a python script in utils to handle the appending of information to the JSON file
       31a092d BFXDEV-673 - add the bfx python wrapper to take care of the creation of the JSON file
       4c575e9 BFXDEV-673 - added a sample list object to the Job class to handle JSON file generation
       c9e342e ucsc bfx tools - bedgraph_to_bigbwig : removal of the temporary sorted bedgraph after bigWig is created
       f7f0113 MethylSeq pipeline - updated version of the pipeline : steps have been condensed, json file for each sample has been added, more metrics have been added to fit with HiHeq requirements, whole genome datasets are now handle correctly
       be4179e MethylSeq pipeline - updated guillimin ini file for bismark_dedup and bissnp requested resources
       3ee60cd MethylSeq pipeline - very minor changes within the base ini file
       34cddb6 DNA-Seq pipeline - added the verifyBamID step to the pipeline (BFXDEV-619) and json file generation add-on
       a638e4c DNA-Seq pipeline - added verifyBamID settings to the base ini file - BFXDEV-619
       52fb55a MethylSeq pipeline - updated common.py so that all the common functions (i.e. sam_to_fastq, timmomatic, merge_trimmomatic_stats, verify_bam_id) now generate a json dump to be append to each sample json file
       1425d42 MethylSeq pipeline - updated scheduler.py to handle the json file generation while the pipeline is running, i.e. adds a json section sample json files as soon as a job successfully ends
       5dae040 MethylSeq pipeline - updated pipeline.py to take care of the creation of the json file for each sample
       8a91523 MethylSeq pipeline - modified sample.py to handle json file creation during pipeline execution
       16c8f93 MethylSeq pipeline - minor update of the verifyBamID python wrapper
       c330a37 MethylSeq pipeline - reviewed picard add_read_group command
       585c92b MethylSeq pipeline - reviewed parameters passed to bismark align
       fc34adf MethylSeq pipeline - removed jsonator import from bfx, waiting for it to be fully implemented
       3674baf MethylSeq pipeline - updated the pipeline with new metrics calculations as well as merged some steps together
       b579365 BFXDEV-619 - added verifyBamID .Rmd template file for report
       91071fe MethylSeq pipeline - added a python wrapper for all the tools related to methylation profiling
       9a2defc BFXDEV-619 - added verifyBamID in pipelines/common.py
       a3b94d3 MethylSeq pipeline - added GC Bias R tools to metrics
       5513358 BFXDEV-619 - added verifyBamID in bfx/tools.py
       168fb69 MethylSeq pipeline - Added ucsc.py, with 'bedgraph_to_bigbwig' which is now called by new bedtools.graph function. Also called from methylseq.py
       c3ce449 MethylSeq pipeline - added bedtools 'coverage' & 'bamtobed' to bfx/bedtools.py, also updated 'graph' & 'intersect'
       946bd70 BFXDEV-661 - DNAseq - added the use module R_Bioconductor within picard_collect_multiple_metrics
       a5b24ac BFXDEV-642 - RNAseq - passed the adjusted p-value as a markdown parameter to the differential expression report
       6e5470a BFXDEV-644 - RNAseq - added the loading of module_python and the use of 'bash pipefail' for step differential_expression_goseq_report
       38db219 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d8a33d9 execute permissions updated
       1d8f3a7 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       dee964a MethylSeq pipeline - create specific ini file to handle capture datasets
       3ce2005 MethylSeq pipeline - create specific ini file to handle capture datasets
       c042cb2 corrected typo in python3.sh install script
       6406b9c added index file as part of the output files in add_read_groups function, in picard.py
       8a940f3 added index file as part of the output files in add_read_groups function
       1c2f36a MethylSeq pipeline - updated wiggle track step with splitted forward and reverse strand reads
       040eefb MethylSeq pipeline - corrected typo in bfx/bedtools.py
       a1c2b6c new CHANGELOG coming along release v2.3.1
       0f801b8 updated Homo_sapiens.GRCh38.sh install script with Ensembl87, dbSNP_149 and dbNSFPv3.4
       ec64899 updated rnammer_transcriptome resources within rnaseq_denovo_assembly.guillimin.ini
       ae6119c updated surpi.sh install script
       ae8ebb5 updated snap.sh install script
       08908b5 updated seqtk.sh install script
       b6d06b8 updated bowtie2.sh install script
       f040ef6 updated RAPSearch2.sh install script
       ad4b019 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       41fdf29 corrected bowtie.sh install script
       e977277 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       79bb77e updated versions of some module install scripts and added some other
       3f90a97 updated versions of some module install scripts
       f6ffda1 MethylSeq pipeline - adjuted resource allocation for bismark align step
       53c3d85 resolving conflicts and merging
       e5e0074 minor changes to resolve conflicts before merging mMaster with MethylSeq
       eb21575 Small changes and updates before merging MethylSeq branch to Master to prevent conflicts
       93d9f01 MethylSeq pipeline - updates done after pur pull request review
       e5eb8fa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       adb1cbe updated version to v346 in UCSC module installation script
       faa8215 AmpliconSeq pipeline - updated version of R_Bioconductor to 3.3.3_3.4
       23d2966 new genome installation script : Apis_mallifera.sh
       9173ee8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       396af25 updated mirbase installation script : removed useless loading of mugqic/R_Bioconductor/3.2.3_3.2
       15d7243 Apis_mellifera.sh
       834765f updated module installation scripts
       996b10f new bash scripts for new modules
       08b52b0 MethylSeq pipeline - few minor corrections before the pull request
       dced67c MethylSeq pipeline - adding of the README.md file
       212451e MethylSeq pipeline - reducing ram to 25G for bissnp step to avoid 'out of memory' issue
       e3d09b5 MethylSeq pipeline - added index creation at then end of picard_add_read_groups step to avoid error in next step (picard_merge_sam_files) when only one readset per sample
       571d0f8 MethylSeq pipeline - corrected ini file for bismark_methyl_call and bismark_coverage2cytosine section where assembly_dir was used instead of bismark_assembly_dir
       f837d91 MethylSeq pipeline - added the creation of the bam index after filtering is done, and updated picard.build_sam_index regarding the ini_section parameter
       dafcd65 MethylSeq pipeline - reducing ram for bissnp step to avoid 'out of memory' issue
       0f4507d MethylSeq pipeline - pipeline with the mapping quality filter step as well as with the picard_calculate_hs_metric step to get ontarget mapping rate
       67e3842 MethylSeq pipeline - updated ini file with newer module versions as well as removing all MUGQIC_INSTALL_HOME_DEV reference and setting them to MUGQIC_INSTALL_HOME
       8f2a7dc MethylSeq pipeline - updates after correction of the section name within the base.ini file
       e6dfe0c MethylSeq pipeline - corrected section name within the base.ini file
       9c5c314 MethylSeq pipeline - corrected typo in the bismark python wrapper
       af1f9f9 MethylSeq pipeline - added tmp_dir parameter to bismark_align step
       ccbf775 MethylSeq pipeilne - addings of .md files for the first steps
       aaf3ba3 MethylSeq pipeline - corrected dependency problem within bismark_align step
       57db20f MethylSeq pipeline - corrected metrics input and step order
       6248dbf MethylSeq pipeline - addings for report generation
       deb45a6 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       02259ac MethylSeq - corrected typo within the picard2_collect_multiple_metrics which caused to always restart the step
       83cecac MethylSeq - Updates and corrections, pipeline is now fully functionnal on test data
       8bd10db MethylSeq - BisSNP step implemented
       77176b4 MethylSeq - debugging bed_graph step
       b58b12d MethyleSeq - correct output files directory for methylation_call step
       d8d7c74 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       1422706 MethylSeq - debugging methylation_profile step
       ae3a9df MethylSeq - updates and corrections, wiggle traks are now generated properly
       4b96e31 MethylSeq - updated pipeline up to 'methylation_profile', still need to work on BisSNP (last step)
       7a2a721 MethylSeq - minor bug corrections
       6b0f820 MethylSeq - First working version of the pipeline : generates bismark alignment and metrics
       89f4a60 MethylSeq - updates of pipeline up to methylation call & also removed/merged some metrics steps
       86a9ae4 MethylSeq - refine metrics calling and preparation, up to pUC19 & lambda reads step
       ce3c86d MethylSeq - updated pipeline until metrics step
       ff5bee3 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       78c8fd0 MethylSeq - updates and corrections up to step 10 flagstat
       5a55472 Methylseq - updates and corrections up to step 10 flagstat
       7b0a9a9 MethylSeq - setps 7 & 8 added
       88dfebd update bedtools version in the /bedtools.sh module installation script
       8f9c198 BEDTools - added intersect function to bedtools.px wrapper
       63b4628 MethylSeq - debugging wrong dependencies du to wrong output file names
       17746b5 MethyleSeq - debugging file names in step 4 & 5 inorder to make step6 picard_merge_sam_files work
       5057bc0 MethyleSeq - defining parameter of step6 picard_merge_sam_files within the .ini file
       b63b02c MethylSeq - preparing step 6 picard_merge_sam_file by reviewing/changing file names upstream (steps 4 & 5)
       c589684 MethyleSeq - bug fix : quoted parameters of add_or_replace_read_groups function in picard(2).py wrappers
       206f128 MethyleSeq - corrected wrong reference to .ini file section regarding step5 picard_add_read_groups
       d0006b1 MethylSeq - missing output file is now passed to the picard.add_or_replace_read_groups function
       cae1cdf MethylSeq - corrected step5 picard_add_read_groups
       6294b65 MethylSeq - updated step5 of the pipeline : AddOrReplaceReadGroups
       0019f59 updated module installation scripts for bowties2 htslib python & samtools
       a8e9c4a PICARD - bug correction in python wrappers (picard.py & picard2.py)
       acd0257 GENAP MODULES - added bismark.sh script for Bismark module installation
       9c70103 MethylSeq - redefined output file for bismark_align step and added AddOrReplaceReadGroups function to picard(2).py wrappers
       222ce02 MethylSeq - python wrappers and scripts updates
       edb1d9b MethySeq - corrected typo in __init__.py
       41f69bd MethyleSeq - creation of the (empty) files as a first commit to the branch

  ehenrion <edouard.henrion@mcgill.ca>      19 commits

       4e7c266 BFXDEV-673 - updated jsonator.py for a better handling of module names & versions
       39aac2d Merged in methylseq (pull request #23)
       11c2e27 chipseq.base.ini : edited module_deeptools to remove reference to mugqic_dev
       bc8a9be README.md updated RAC_ID export line
       a29213c README.md edited : added the export of the $RAC_ID variable which will be used on Cedar for job sumission
       5f66230 README.md edited : added $MUGQIC_INSTALL_HOME_DEV setting for cedar
       27a43ac MethylSeq pipeline - edited guillimin ini file : more walltime for gatk_depth_of_coverage & more cores for bismark_dedup
       0c5f3a0 MethylSeq pipeline - changed flagstat output name so it is more obviously related to ontarget bam
       3c1f730 MethylSeq pipeline - inverted order of input file for methylation_call step : "*.readset_sorted.dedup.bam" is now set before "*.sorted.dedup.bam" to avoid unnecessary sorting of bam files...
       f79719e MethylSeq pipeline - GCbias and mapping_qual_fileter jobs have been added to metrics, their commands were created but not submitted...
       e77d93e MethylSeq pipeline - methylseq.guillimin.ini adjusted some ppn values for guillimin
       793836e MethylSeq pipeline - methylseq.base.ini edited bismark align parameters
       4d6b5d0 MethylSeq pipeline - methylseq.base.ini added module_R
       c6a2446 methylseq.py : added missing variable
       4f247e6 README.md edited online with Bitbucket
       321816e README.md edited online with Bitbucket
       0bb51e1 README.md edited online with Bitbucket
       a42811f bedtools.py edited online with Bitbucket removed unused and unfinished function genomecov...
       e76018c bedtools.py edited online with Bitbucket

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      43 commits

       33f242d removed depreciated reasignment of module_snpeff in snpsift_annotate step in dnaseq.base.ini
       bad3d89 renamed and duplicated snp_effect section in mammouth.ino to mpileup_snp_effect and haplotype_caller_snp_effect in order to correctly set ppn=1
       54251b5 Moved report.kallisto job after kallisto_count_matrix since the it needs the output of kallisto_count_matrix; added dependancies to copy_tx2genes_file so it waits until kallisto_count_matrix is done
       19a3beb Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       73b4c97 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eac4107 fixing local merge conflict resolution issue
       e4a5abe (hopefully) fixed conflict
       1596188 added a few more os.path.join; added an import to tools
       d4e5e0f removing kallisto.sh from module/dev
       5f0a9b8 remove all compiled files
       220ed53 adding os.path.join where possible
       030d09d remove unused sections in mammouth and guillimin ini files
       87289ad resolve merge conflicts
       46cb645 remove commented lines
       2941bfa move rmarkdown calls to appropriate step
       b25aef3 fix ini file to correctly set walltime
       76e4c73 fix an issue with dash in sample names in report
       fa04e0a minor improvments to kallisto report
       b55eb98 change path to file in kallisto.Rmd
       6d3e20f change job names to fit sections in ini file; replace mention of samples by readsets in md files
       4fd75b2 Changed R code in kallisto.Rmd; add more columns and fix error in the stat table; changed kallisto method text
       1728a1c add step to generate transcript count matrix
       341b0fd Change text RNAseq report for differential expression
       1431e2e update report files
       9bb9dbf added option for single reads, added new parameters in ini file
       f9310f0 added trimmomatic stats to report
       f81641c added first version of RNAseq Light report
       7e53c0f added a mkdir command to the merge script
       efb008f Change path of merged abudance file
       d8dedd1 fix job dependancies
       1639f59 fix exploratory function
       9d875eb change path for call to rnaseq_light_kallisto.sh
       4fc4f3b adding new step for exploratory analysis
       3831ba6 added a step for merging individual abundance files
       d720876 add call to module_tools
       d2996e1 fix path to abundanceTranscript2geneLevel function
       5dc4760 adding configuration ini files for mammouth and guillimin
       6c2900a make it clear transcrptome must end by idx
       feca8f7 fix path, disable exploratory
       aed4e03 new RNAseq_light pipeline with kallisto, dev version
       053dbf0 added RSEM 1.3.0
       2c86a7f updated kallisto
       74e5c0c update sailfish to 0.9.2

  eloi.mercier@mcgill.ca <eloi.mercier@mcgill.ca>      1 commits

       c82520a Approval granted by Rola. Merged in RNAseq_light_dev (pull request #24)

  eloi.mercier@mcgill.ca <emercier@abacus1.ferrier.genome.mcgill.ca>      3 commits

       22084d1 revert sailfish version change
       c8ecd2c ajout Salmon 0.8.0
       e2986c5 mise a jour kallisto 0.43.0

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      21 commits

       6150a8a remove bad output in ihec_metrics_report
       07709bf Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cc3d29d MethylSeq - add MethylSeq.picard_remove_duplicates.md report file
       f8a97db fixing conflict in resources/genomes/install_genome.sh
       2a9190b change Dedup from bismark to remove duplicat from picard
       9a9b70f change ihec methylseq metrics to work on single sample
       7206787 methylseq - generate ihexc report per sample and move methyl_profile lib to tools lib
       900f112 add specific temp dir to the sort step while generating bigwig
       a9b1ed3 increase general recalibration walltime in  Dnaseq to 96h
       764bb1d Allow bedtools.graph to support not having the other_options set in the ini
       961e0e1 removing .DS_store file
       fc1defd ChIPseq - address reviewer coments - BFXDEV-675
       f8bf304 Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       cb48b45 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       1162758 RNAseq - add concat job  - BFXDEV-668
       c859120  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ee2eab0  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       188462d  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ba1f73c  RNAseq & ChIPseq-   Update (chip) and debug (Rna) IHEC metrics steps - BFXDEV-668  - BFXDEV-675
       9e71317  RNAseq -  implement IHEC RNA metrics step - BFXDEV-668
       7cfbb6d  RNAseq - BFX - implement mugqic_tools module for the IHEC RNA metrics generation script - BFXDEV-668

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      10 commits

       34124b6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a9f69a4 add protocole compatibility to rnaseq_light
       cd1caeb Make other pipeline supporting several prtocols - BFXDEV-692
       ffdf4f6 Create 2 protocols with different steps for hicseq - BFXDEV-692
       554f8ad Make other pipeline supporting several prtocols - BFXDEV-692
       5ac65fc allow several prtocols with different step list - BFXDEV-692
       e634cd2 Create 2 protocols with different steps for hicseq - BFXDEV-692
       ece6ab0 test multi protocole pipeline
       b109814  RNAseq & ChIPseq-   Update ini file for the new release of mugqic_tools 2.1.9 - BFXDEV-668  - BFXDEV-675
       e057114 ChipSeq - finish ihec metrics, preprocess and reformat -  BFXDEV-675

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       f19e749 Merged in IHEC_metrics (pull request #21)

  pascale.marquis2@mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      2 commits

       852ef1c /#cluster_queue
       eb94727 python/2.7.12

  Pascale Marquis <pmarquis@lg-1r17-n03.guillimin.clumeq.ca>      1 commits

       1052ba1 update tumor_pair.base.ini

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       df577fb Minor bug fixes and addition of base exome ini

  Rola Dali <rdali@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       46371a9 starting the hicup_align step

  Rola Dali <rola.dali@mail.mcgill.ca>      88 commits

       9be3474 Merged in IHEC_metrics (pull request #27)
       130d0b5 changes to mammouth.ini to set all ppn=1; changed module spider back to module show since it is incompatible with abacus
       71cbe96 Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       f93c39e module show changed to module spider in config.py to accelerate module checking
       a64c34f fixing homer dependencies
       d274293 editing homer tag directory output back to folder
       78dc2fc Merged in IHEC_metrics (pull request #26)
       6b97871 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       ac73327 homer edits to generalize methods
       649af97 resolving merge conflicts
       fd5cb38 dependencies in rnaseq metrics
       23921c9 rnaseq metrics
       62da3e6 fix ihec_metric job names
       fb83174 adding chip_type to ini to avoid pipeline crash
       b7e2090 run_spp fuctional
       85dd69a added run_spp to calculate dnsc and rsc. TO TEST on mp2
       1459b13 adding TMPDIR to bigWig creation
       a058aa1 sample name change in ihec_metrics
       31c83eb fixing TMPDir issues on guillimin-to test
       80e315f pipeline to run without input & chnage in ihec_metrics format
       14fde80 moving homer code to homer.py moduke
       d33094e Merged in hicseq (pull request #25)
       fe6d261 resolving juicer conflict
       3cd0685 resolve merge conflict
       eff83a2 resolving merge conflict
       574b511 adding C3G logo
       400b590 commit for pull request
       18c1112 commit for pull request
       27b87a3 adding bedops to annotate capture file
       198416b adding runchicago to capture hic
       742376e Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c26ea1d Resolved merge conflict
       cc289e6 commit before Mathieu's pull
       64cd090 adding capture seq
       46b4f6b merging capture hicseq with hicseq
       4b96552 chicseq
       b19d167 Merged in hicseq (pull request #22)
       a25f08b added new chromosome contigs based on ensemble and NCBI genomes BFXDEV-670
       11b2bb9 changing perl to perl/env  BFXDEV-670
       a790fa1 added genome.py for genome related methods BFXDEV-670
       6a63124 changes for pull request BFXDEV-670
       f9a2552 changes for pull request BFXDEV-670
       0d65329 changes for pull request BFXDEV-670
       aa2e874 edits for merge request:wrappers. BFXDEV-670
       de38b5a BFXDEV-670 hicseq merge edits
       d69045f creat_hic_file edits BFXDEV-670
       1649642 adding juicer.sh installation script BFXDEV-670
       665ea38 Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c4adfa9 testing create_hic_file BFXDEV-670
       4a67afe modified hicseq.briaree.ini and batch BFXDEV-670
       a37fe8a deleted files moved to mugqic_tools BFXDEV-670
       bcf5b1a moved genome_digest before ncrna method which is failing in install_genome.sh BFXDEV-670
       861cc79 added module install files and genome digest BFXDEV-670
       1523060 changes for pull request
       ad4a141 commit before merge changes
       98d7947 added samtools_bam_sort BFXDEV-670
       6fd3b93 edited ini files BFXDEV-670
       18e0774 samtools_bam_sort and multiqc_report in testing
       3c9468b Fixed hicup/TAD/HiCPlotter restrart bugs BFXDEV-670
       6a58aa4 split chr interaction matrices from plotting BFXDEV-670
       dca827d HiC v1.0 is ready for testing on guillimin and abacus BFXDEV-670
       2b40a8c bam merging in testing BFXDEV-670
       65763b1 samtools_merge_bams in testing BFXDEV-670
       d2ff711 added petagDistDistribution to HomerQc plots BFXDEV-670
       e8edda6 added output_dir property to reduce code redundancy
       ded82f6 fixed path in identify_compartments
       1b31f40 added identify_peaks BFXDEV-670
       d212eba compartments and TAD identification now working BFXDEV-670
       384d836 chr and genome interaction system now working BFXDEV-670
       6714aab interaction matrix plots testing BFXDEV-670
       edadc74 homer archiving and Qc plotting now working
       8c94b93 first 6 steps working BFXDEV-670
       5c3eb28 resolving git merge issues
       e1d46d6 Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       399cc72 syntax changes
       c545fc6 added fastq_readName_Edit() BFXDEV-670
       2b7ee9e mammouth.ini cpu set to 1
       4afc530 genome digests on mammouth
       e4e52d3 update genome digest files
       9277398 homer_tag_directory archiving in testing BFXDEV-670
       acb862f basic formatting changes
       53257d6 make_tag_directory in testing
       148eb65 works to produce hicup bam and library Qc. ITS ALIVE :)
       958242d ITS ALIVE
       c446402 need to expand  variables
       099e5ef pipeline now accepts enzyme
       ffbdf83 hicup_align producing script; to test tomo
       884af26 initialising hicseq analysis pipeline

  Xiaojian SHAO <xshao@lg-1r14-n04.guillimin.clumeq.ca>      2 commits

       6ddab51 add ppn to wiggle_tracks step
       06bd6a7 add ppn to wiggle_tracks step

  Xiaojian SHAO <xshao@lg-1r17-n03.guillimin.clumeq.ca>      3 commits

       6fd6ec5 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cd938a3 add walltime to bismark aligner
       8a5b826 ppn_Changes_in_Guillimin.ini

  Xiaojian SHAO <xshao@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       db2901f methylseq: edit on ppn setting. -xiaojian

2.3.0        Mon Feb 27 13:40:01 2017 -0500        82 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      13 commits

       7717c35 Pre-release - adding the installation scripts for all the new software used by tumor_pair pipeline
       8eada67 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8e17bca Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       59cf940 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       026ccb1 BFXDEV-602 - PacBio Assembly pipeline now contains a new optional step (step 10) to circularize the successful/clean assembly contigs with circlator - also comes with some minor unrelated updates
       06e095f updates brought to many module install scripts, and adding of new modules
       9d565ba RNASeq - corrected a typo inserted after correting bedtools bug...
       982a014 RNASeq - corrected a bug bedtools.graph function : samtools_options now handles reverse strand specific parameters, avoiding an empty begGraph for reverse strand
       72d1a3f updating python & python libraries installation bash scripts
       6242c97 BUG correction within mpileup function : parameters were shift after introduction of 'ini_section' parameter for tumor_pair prpose
       1ca38fa DNASeq - bug correction after merging tumor_pair branch to master
       e6b482b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d8876d FastQC - updated bash install script to match our standards

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       0e875a5 README.md edited online with Bitbucket

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      39 commits

       d6be75a remove .gitignore
       3821384 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8e31719 rnaseq - increase abacus ressource for STAR
       6f5a17b remove conflict with tumor_pair changes
       07cab4f remove bfx/samtools.py conflict with tumor_pair changes
       abecc6c remove bfx/picard.py conflict with tumor_pair changes
       4368b7c merge and remove conflicts
       2d5b7c5 merge and remove conflicts
       dfcbf62 tumor_pair- update germline_loh ensemble
       376ff92 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       07b9e8b add encode  sample+readaset+pipelien to log system
       f4c9b60 Tumor pair - add comment on the bed system
       ba47702 Tumor pair - change vardict input output to correct deepndency
       0e20097 Tumor pair - change varscan input output to correct deepndency
       4ab037e sequence_dictionary.py - split_by_size - correct bug create a first job with everybody when the nb of split was too high
       b5a7b33 tumor_pair - add output bam index file when only 1 realignment is produced
       a5b4491 tumor_pair - remove file name bug in indel realignement mv
       4429b19 tumor_pair - rebuild gatk indel recalibration input/output scheme in tumor
       b9030e5 bfx/gatk - indel realigner add target list file as input for dependency
       440001c tumor_pair - remove duplicated lines in guillimin.ini file (from base.ini)
       6c1e606 tumor_pair - correct symlink creation
       b15f118 tumor_pair - remove few issues (paths, dependency, input/output)
       329a56c tumor_pair - put WGS as default (instead of WES)
       271ba29 remove issue with uncorrect interval list when no bed is attached to the project
       517e8c3 adjusting tumor_pair
       c049fe2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       c2c0986 tumor_pair - debug bfx/bcbio_variation.py
       5c6034d tumor_pair - remove conflicts
       d2c0015 rtumor_pair - emove conflicts
       7c4bd8a remove conflict
       c67b837 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       9cbd222 tumor_pair - correct typo in scapel download adress - BFXDEV-477
       760faaa Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       8d2a243 tumor_pair - add scalpel install script - BFXDEV-477
       e6beef4 tumor-pair - mv tools.py from highcov branch to tumor_pair - BFXDEV-475
       1bc346f tumor-pair - adding bedfile spliting process - BFXDEV-476
       5acad27 tumor_pair - extract tumor_pair code for the high coverage branch - BFXDEV-475
       6c73c2b tumor-pair - adding bedfile spliting process & start implementation in tumor_pair.py - BFXDEV-476
       905c1ac tumor_pair - extract tumor_pair code for the high coverage branch

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       b5b1d51 tumor_pair - add comments to dict2beds function - BFXDEV-521
       e0590a4 tumor_pair - add feature to dict2beds function - BFXDEV-521
       7c0efb4 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       bd7e169 tumor_pair - modify ini file to take into account the new version of scalpel (CVMFS) - BFXDEV-477
       8e735ac tumor_pair - add space charter before scalpel option - BFXDEV-478
       33a9347 pull origin Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       5ed5783 tumor_pair - add the two pass option to scalpel - BFXDEV-478
       bf963a3 tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing - BFXDEV-476
       18110f5 tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      5 commits

       192f4d3 Bug squashes and speed improvements to ensemble processes
       40df5d9 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       93a2142 Modified to remove analyses of alt contigs + split Varscan2 by chromosome
       8d3f03f Beta version: module files created and code tested on wes and wgs on abacus and guillimin
       755b0c0 update to Mutect2, added Vardict, re-added samtools and added ensemble approach

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      4 commits

       e94352f Fixes to samtools modules calling for pairs
       9729909 Numerous speed improvements and addition of fast variant calling
       3f1aceb Additional fixes to resource allocation issues
       0f04f2a Dealt with comments, added paired indel realignment, varscan2, and seperate somatic and germline call files

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      11 commits

       67b5c53 Debugging samtools pair bcftools
       3a4996e Fixes to samtools due to dnaseq changes
       a920ef8 Update README.md
       2c7ce40 Bug squashes and preformance improvements to ensemble process
       3349a81 Bug squashing and speed improvements to ensemble process
       15a7f87 Dependency and other bug fixes. Guillimin ini and install scripts updated
       a21682a tumor_pair - fixes to tumor_pair README.md
       1966434 tumor_pair - fixes to README.md and resolving of conflicts
       c53612c tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       dd6df55 tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       c708f88 Updates/fixes from guillimin test

2.2.1        Mon Dec 19 10:57:33 2016 -0500        212 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       cd7c8e7 Added proper error message when running script with too old Python version.

  Edouard Henrion <edouard.henrion@mcgill.ca>      138 commits

       68f0366 GenAP Pipelines - updated tmp_dir variable within all the .base.ini files for a better use of the memory in the compute nodes on abacus
       c364c3e modules - updated version of VSEARCH (2.3.4) within the module installation script
       ed99e0e GENAP PIPELINES - updated all the .base.ini files to set tmp_dir to /lb/scratch/ehenrion instead of /lb/scratch/ on abacus
       d0faf2e RNASeq - added __init__.py
       a8d1101 PICARD - bug correction in python wrappers (picard.py & picard2.py)
       7d2e9b2 install_genome.sh - corrected MUGQIC_INSTALL_HOME_DEV link
       b2ec14e install_genome.sh - set abacus2 variable to manage mammouth execution properly
       c3549ed install_genome.sh - updated for a better handling of mammouth cases
       4726932 Homo_sapiens.hg19.sh install script - updated URL to retrieve dbSNP vcf
       b0dfe82 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       465b913 RNASeq denovo Assembly - updated mugqic_tools version to 2.1.6
       55e6456 MUGQIC MODULES - updated bedtools version to 2.26.0 within the bash install script
       797fe54 MUGQIC MODULES - updated samtools version to 1.3.1 within the bash install script
       f846746 MUGQIC MODULES - updated bowtie2 version to 2.2.9 within the bash install script
       7ef8596 MUGQIC MODULES - updated bowtie version to 1.1.2 within the bash install script
       eca9f4d MUGQIC MODULES - updated bismark version to 0.16.3 within the bash install script
       b92f63f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2b0f2e0 RNASeq_denovo_assembly - updated guillimin specific .ini file for trinity & insilico_read_normalization steps
       5c11609 DNASeq - updated guillimin specific .ini file for snp_effect step
       639668a DNASeq - updated .ini file for snp_effect step
       14f4c0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b3a054e PICARD 2.0.1 - updated picard2.py wrapper with proper syntax
       7a20e2f briaree pipeline .ini files
       b4b5e8c RNASeq_denovo_assembly - updated mammouth .ini file
       b3ffcc7 python.sh - updated python version (2.7.12) within the bash installation script & added all the python library installation steps : no need for python_lib.sh anymore
       731b029 Homo_sapiens.GRCh38.sh - updated Ensembl version (85) & dbNSFP version (3.2c)
       ba15075 install_genome.sh - update cmd_or_job function to automatically handle mammouth environment cases
       5d07d31 exonerate.sh - update version to 2.4.0
       c1ac0e4 RNASeq_denovo_assembly - back to hmmer version 2.3.2 (from mugqic modules)
       1f64757 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       16c3464 exonerate - within the bash install script : corrected broken link to the archive
       875e3eb Genome reference .ini files updated for each reference installed on CVMFS
       8c4c2ed module install script - corrected emboss.sh : unresoved conflicts are now resolved
       21a06d7 module install script - corrected emboss.sh : unresoved conflicts are now resolved
       5f46cf6 Python packages - added TEToolkit package to the installation scripts
       f4a2adb FASTQC - updated installation script with newer version and integration of /lb/project/mugqic/analyste_dev as a possible installation path
       a396762 RNASeq_deNovo_assembly - updated ini file with integration of picard v2.0.1 and working version of trinity (2.0.4)
       916b8d6 PICARD 2.0.1 - added picard2.py to the bfx tools
       ef03d6b commiting minor updates, prior the new release
       aa0f55f RNASeq_denovo - briaree specific config for insilico_read_normalization_all
       29fc52a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd77e98 corrected version of mugqic_tools in rnaseq.base.ini
       daad645 corrected version of mugqic_tools in rnaseq.base.ini
       e2d140e RNASeq de-novo assembly - bug correction when calling trinity with new default parameters
       87d6833 RNASeq de-novo assembly - updated ini file for briaree
       1605b8f modules - new STAR version handled in the installation script
       8f88daa RNASeq de-novo assembly - updated trinity call with newer version of trinity and new default parameter
       73ca888 RNASeq - update differenial expression to handle local fit option if needed (i.e. if parametric dispersion fit fails)
       ec71e94 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d0d765 new reference genome installation scripts (new species)
       0029c2a metaMarkerSeq - guillimin .ini file adjustments regarding qiime
       e4c1a5b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b5cc05d added butter installation script
       8ec7a74 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       16f2ae7 qiime db installation scripts and some new module scipts
       8fe19f8 chimera_unite.sh added
       ad09a1e new module installation scripts
       a4f4b05 some updated & new module installation scripts
       8fd6a71 updated genome installation scripts
       083ff4e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       846d931 metaMarkerSeq - removed bug with amplicon_type within the .ini file
       2ba2d4c metaMarkerSeq - updated python & bash script as well as .ini with the newer silva db version
       4504c46 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e68ffa3 RNASeq - fixed bug in gq_seq_utils_exploratory_analysis_rnaseq
       903ef8d mammouth ini files added for dnaseq_high_coverage
       def0206 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c822e0d mammouth ini files updated for pacbio & dnaseq
       887e518 adding pipeline ini files specific to brairee environment
       e5e3585 DNASeq - review some parameters within ini file during the pre-release pipeline testing
       cece203 RNASeq - allow sample names to have '.' within their name without crashing during gq_seq_utils_exploratory_analysis_rnaseq
       46d069f fixed typo
       54f260b RNASeq - fix mammouth setting for htseq_count step in rnaseq.mammouth.ini
       ada6738 RNASeq - restoring bfx/report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd
       35c5cad ChIPSeq - metaMarkerSeq - adding briaree specific ini file
       f0112e1 RNASeq - correct default genome in the ini file
       4fb234d pushing ampliconseq.guillimin.ini after testing on guillimin
       90ff7dc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4ec8c82 pushing ampliconseq.mammouth.ini after testing on mammouth
       4d8e583 Corrected typo in dnaseq.base.ini
       6e9c3b2 correct mammouth master branch divergence
       8ab80b3 updating module and genome scripts from mammouth
       7e28490 remove conflict
       2c439fe conflicts resolution
       4db2a9e some residual diffs to commit
       de4c043 redo the call of bcftools module within samtools.py
       b69fa6e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ad24f0c resolving conflicts after merging...
       6330a7f commits before merging
       5e4eac4 DNASeq/RNASeq - code correction after testing on briaree
       580589a metaMarkerSeq - Merging branch to master git add resources/genomes/Bos_taurus.UMD3.1.sh resources/modules/perl.sh resources/modules/star.sh resources/modules/vcftools.sh
       0fae673 metaMarkerSeq - reformat the report by removing some redundancies in the paragraphs and ensuring correct links to the tables/figures - BFXTD-26
       8873109 BFXTD-26 - corrected the template and links used for report generation - MetaMarkerSeq
       a0c0188 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e18bb24 correct a bug in wiggle step which causes *forward.bedGraph and *reverse.bedGraph files to be identical
       850b148 metaMarkerSeq - commits after reviewing pull request comments sent by Marc Michaud - BFXDEV-453
       2390160 threshold instead of treshold
       51848d9 metaMarkerSeq - debug Krona input generation and updated md files for report - BFXDEV-453
       475035b meatMarkerSeq - get rid of remaining  variables - BFXDEV-453
       0fd827d AmpliconSeq - debugging report template files
       076e6f7 bug fixes and file names correction
       0441415 ampliconseq - tuning of the ini file - BFXDEV-434 BFXDEV-448
       feac3c0 ampliconseq - new/updated module install scripts related to the pipeline needs - BFXDEV-434 BFXDEV-448
       3865cfa ampliconseq - code review & debug prior to relase - BFXDEV-434 BFXDEV-448
       664910a AmpliconSeq - modified some parameters in .ini file & other minor syntax changes within the wrapper - BFXDEV-434
       b22e6e2 ampliconseq - conflict resolution - BFXDEV-434
       28788ce ampliconseq - code review of ampliconseq.py - BFXDEV-434
       d07b3d8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ba4bca9 BFXDEV-523 - error correction & version update in the gemini installation script
       4e7e0a6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f22df95 BFXDEV-516 - Macaque genome installation bash scripts + other minor updates
       90d34cd BFX-4532 - added prepend-path entry for dnaseq_high_coverage in mugqic_pipelines.sh
       7f41180 RNA-Seq - remove useless parameter from rnaseq.base.ini
       7442c01 RNA-Seq - increase default walltime for star_align & star_index steps
       c9d3694 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a9e67a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       839b257 genomes - updated versions of genome installation scripts, essentially fixing STAR indexes installation - BFXDEV-494
       a3da825 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       192a971 STAR - updated STAR version (to 2.5.1b) in the installation script (star.sh) as well as in the RNA-Seq pipeline .ini files - BFXDEV-514
       646ad6a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       394913d Version bump to 2.1.0-beta
       0dac262 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7a78685 BFXDEV-490 - recommit after conflict resolution
       8fc43e3 BFXDEV-490 - updating chipseq.base.ini, dnaseq.base.ini, pacbio_assembly.base.ini, & rnaseq.base.ini
       cfd56b9 ampliconseq - module - adding vsearch.sh in the module script folder - BFXDEV-434
       fb1f4d2 ampliconseq - base ini file updated - BFXDEV-434
       ed7611f ampliconseq - updated ampliconseq ini with new module versions - BFXDEV-434
       8c047cc ampliconseq - updated module versions - BFXDEV-434
       565ab99 Merge branch 'ampliconseq' of bitbucket.org:mugqic/mugqic_pipelines into ampliconseq
       d2ef4f4 AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEV-434
       b15af32 AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEX-434
       ace14e4 AmpliconSeq - revised and standardized code of ampliconseq.py - BFXDEV-434
       76c924e AmpliconSeq - updated ampliconseq.base.ini - BFXDEV-434
       83d3b0e AmpliconSeq - new VSearch library added to bfx - BFXDEV-434
       96ce607 AmpliconSeq - updated tools.py with AmpliconSeq functions in bfx - BFXDEV-434
       b886764 AmpliconSeq - new Qiime library added to bfx - BFXDEV-434
       bb008af AmpliconSeq - new Krona library added to bfx - BFXDEV-434
       d9016c8 update of resources/modules/mugqic_tools.sh - BFXDEV-490
       0611584 ampliconseq - updated syntax and some corrections

  Edouard Henrion <henrione@briaree2.rqchp.qc.ca>      2 commits

       91a47a5 RNASeq - adding .ini file for briaree
       ce7f2ae small adjustments after cloning/before testing on briaree

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      1 commits

       5806c4a ampliconseq - bump README.md to the latest version

  ehenrion <edouard.henrion@mcgill.ca>      3 commits

       62133d0 README.md edited online with Bitbucket
       d988c47 README.md edited online with Bitbucket
       5edde68 README.md edited online with Bitbucket

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      1 commits

       0c86240 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       7227342 mini and spades modules
       be71d1c nxtrim and quest modules updates
       5ec026e dev install scripts for MAKER

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a157437 AmpliconSeq - remove ini recursive interpolation
       84a3e9b remove conflict
       0e5bd2a High coverage - debug vt and gemini step - BFXDEV-542 BFXDEV-541
       28ce1e1 ampliconseq - remove dependencies issues and remove try instances - BFXDEV-531
       342169d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ed67908 DNASEQ - bwa always set LB in RG tags: when library barecode not given used sample name instead - BFXDEV-354
       e28eff7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       73bce2e ressource - python lib: add futures for multithearding
       a7c0b13 update install module general script

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       2f9db74 HighCoverage - add missing README.md file
       21f73de Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5e66e13 Removing unmaintained pipelines (PUURE & rRNATAGGER) from master
       435f7d5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3366710 RNAseq synchromize mammouth ini with the base ini (missing tuxedo_hard_clip) - BFXDEV-515

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      5 commits

       8a0bb76 dnaseq.py  - gatk DoC job name wasn't following the rule: <stepName>.<sampleName>.<OtherInfo> . corrected
       5056c01 Merged in irp_modifications (pull request #17)
       4991ce6 change mail for bug & issue in  README.md edited online with Bitbucket
       ebe9e98 README.md edited online with Bitbucket
       0cdc302 README.md edited online with Bitbucket

  mmichaud <marc.michaud@mail.mcgill.ca>      20 commits

       ed1aece Increase max jobs per step to 24 (cover most case without overloading the cluster) (cherry picked from commit 4b06703)
       5836649 BFXDEV-559 Error when a lane contains only double-index librairies on a single-index run (cherry picked from commit 19bd532)
       b2f9c7e Increase even more STAR walltime to allow fetching the reference index from cvmfs.
       a8b9ffe Don't mark the blast job as failed if there is no blast results. If grep can't find any results, it exits with a code 1. After, when we test the $PIPESTATUS, we intercept that error code and mark the job as failed. By adding a " && true", we preserve exit codes but we reset the $PIPESTATUS.
       5c5d933 Increase STAR walltime to allow fetching the reference index from cvmfs.
       61b608b BFXDEV-538 Add perl for bed2interval_list
       4d72471 BFXDEV-538 Update modules version
       e39b202 BFX-4696 Use a newer version of python having all the needed packages.
       d141b2f BFXDEV-544 VerifyBamId: Don't filter the variant file with the genomic target file.
       0237e8a Merge branch 'master' into irp_modifications
       5849c79 Illumina Run Processing: Accept "." in asembly name in the genomic database. (fix regex)
       e020c64 Illumina Run Processing: Accept "." in asembly name in the genomic database.
       8806df0 Illumina Run Processing: Use cvmfs version of verify_bam_id.
       6735cae BFXDEV-512 Illumina Run Processing: Use "Parc" silva database.
       8aa74f2 Illumina Run Processing: Add "nanuq_environment" configuration variable to be able to ease testing on the QC environment.
       3d2bec9 BFXDEV-533 Illumina Run Processing: The exclude bam option doesn't exclude the sorted bam.
       776850b BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       7d1694d BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       edfb2e4 BFXDEV-527 Illumina Run Processing: Merge jobs of a task when their count exceed a configurable threshold
       cd6cfde Tweak cluster parameters : Use parallelGCThread for all available cores. BvaTools run faster with only 4 threads.

  ptranvan <patrick.tranvan@mail.mcgill.ca>      24 commits

       46ee612 Eulidean distance option for PCOA
       8ca5e0e README link correction
       18b9818 Module perl correction
       d63b313 Module perl correction
       b4fc989 Format edition
       39c2882 . file deletion
       c2cd2c3 BFXDEV-434 README update
       e0d27cc BFXDEV-434 Create inis for other clusters
       f4475b9 BFXDEV-434 Deploy db files in $MUGQIC_INSTALL_HOME
       dc7f9fa - Other pipelines have been added to the branch. - Create inis for other clusters, leverage overlaading of parameters. - put db parameters in DEFAULT section,. - QIIME section has been exploded in the config file. - Tutorial in README file.
       27f3fa5 README and configuration file correction
       5a9fdb7 Configuration file modification.
       3419fac Adding alternative normalization method
       a7ade41 Dependance correction
       d4e383d add a step (close ref)
       6b63f3f CLuster algorithm change
       a5d8876 tutorial.txt modification
       72ef558 Minor modifications (report and option)
       9cf47dc Amplicon-Seq pipeline + report upload for test.
       c9c6125 Change name: metagenomic pipeline to Amplicon-Seq pipeline Adding new features (alpha, beta plots) Final step implemented
       d7db137 Adding metagenomics pipeline (1st part)
       075ec74 Adding metagenomics (amplicon) pipeline.
       990ba4b Merge branch 'metagenomics' of https://bitbucket.org/mugqic/mugqic_pipelines into metagenomics
       f779d06 test

2.2.0        Mon Feb 8 12:04:44 2016 -0500        405 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       dc03d86 Added link to the GenAP project in front page.

  Edouard Henrion <edouard.henrion@mcgill.ca>      63 commits

       477b332 Version bump to 2.2.0
       caa197d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d87ea7 BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       371dd66 add report section in DNA-Seq High Coverage pipeline ini file
       514d5b6 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       be6d52a updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       61ac85c ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       042bb03 modules - updated vcftools VERSION in vcftools.sh -  BFXDEV-490
       f71fc01 minor update in module file created by perl.sh
       d3c2bbf add report section in DNA-Seq High Coverage pipeline ini file
       cd5d58c added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       7d0e6f0 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       7bf0a39 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       349ae58 modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       ab4b16d module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       9e04f6a modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       e93f812 minor update in module file created by perl.sh
       fe2a9cb rnaseq & rnaseqdn ini files updated after testings on guillimin - BFXDEV-490
       39e91f1 add report section in DNA-Seq High Coverage pipeline ini file
       e3bce29 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       ecc918b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a4f7977 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       442830d updated ini files - BFXDEV-490
       175b162 BFXDEV-490 - resolving conflict on guillimin
       36b6f5d BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       d1b51cb Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       a3e2ab8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ec80af3 modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       b6e687a minor update in module file created by perl.sh
       97bd155 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       682edbe mugqic_tools.sh VERSION=2.1.5, again...
       fd3f0e3 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       424768a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4b4c08a ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       f6415e1 modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       045c237 module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       d2072f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c22a596 modules - corrected typo in perl.sh when creating module file - BFXDEV-490
       f89dbb3 modules/genomes - updated module version calls as well as database releases - BFXDEV-490
       3970c21 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       17531c7 module - debugged scalpel.sh script - BFXDEV-490
       a05c71b module - updated scalpel version - BFXDEV-490
       1c9d868 module - debug  in snap.sh  - BFXDEV-490
       eee7e56 module - updated version of Python (to 2.7.11) in gemini.sh - BFXDEV-490
       8018094 module - removed module load calls in shortstack.sh - BFXDEV-490
       8c69953 module - updated star.sh - BFXDEV
       071a17b module - updated shortstack.sh snap.sh - BFXDEV
       3fa5330 module - corrected typo in sphinx.sh - BFXDEV-490
       976b856 module - vcftools updated version to 0.1.14 - BFXDEV-490
       0d5709b module - debugged gemini.sh - BFXDEV-490
       add9770 module - samtools updated to 1.3 - BFXDEV-490
       15428f8 module - debug the name of the archive - BFXDEV-490
       7ab3bbf module - UCSC version set to v326 instead of 'latest' - BFXDEV-490
       16571fc update version of bwa module - BFXDEV-490
       bf6422a some more updated modules - BFXDEV-490
       0cf2a4d module & genome updates - BFXDEV-490
       39c1758 module updates for to the release - BFXDEV-490
       bb09f04 resources/modules/dev/epacts.sh has been removed (really)
       0815aec resources/modules/dev/epacts.sh has been removed (really)
       27058bc BFXDEV-490 - update of resources/modules/mugqic_tools.sh
       0021aaa Merge branch 'gatk_variants' of bitbucket.org:mugqic/mugqic_pipelines into gatk_variants
       a41229e gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted - BFXDEV-481
       b3ae471 gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted

  Edouard Henrion <henrione@ip03.m>      2 commits

       839c464 updated chipseq.base.ini after mugqic_tools update - BFXDEV-501
       476279f updated chipseq.base.ini after mugqic_tools update - BFXDEV-501

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       565ea8b BFXDEV-491, missing argument to re.sub for the SINGLE_END STAR commends case
       bc9c0df Reverte back to older scheme. Was assuming bacteria only.

  Francois Lefebvre <lefebvrf@gmail.com>      37 commits

       ae9667f Mammouth rnaseq ini required cluster_cpu=-l nodes=1:ppn=1 for new steps related to rRNA estimation
       e21b18c -S flag was missing from the last samtools view command in the hard clip command
       f5734af Install scripts
       6eabc35 minor mod to R_Bioconductor (removed tabs in HERE DOCS)
       432b3ab sleuth R package
       93c3b1d Removed old R installation scripts
       7758ea3 Added slash to URL to be able to retrieve latest Bioc version
       b80611a sspace-longread dev install script
       b626deb Kallisto install script (abacus only)
       0f9449a popoolation install scripts modifications
       dd9e031 more notes on pacts
       1e9442f Added notes to EPACTS install script
       ee112de Fastx and EPACTS install scripts
       a1b5289 install scripts for ray, proved, sortmerna
       a086744 Added celera and loRDEC install scripts
       c56d297 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01d05df Merge branch 'knitr_test'
       62efc5f Added PBJelly install script + networkx python module
       1691002 Fixed typos
       394fca5 Replaced library calls in rmarkdown.py with direct ::render call. Added missing section header to exploratory report.
       7d61414 no message
       cf13caf Modified exploratory analysis report text
       a8fd54f import rmarkdown in wrapper
       abd9f89 damn dots
       336b156 Testing out markdown::render paradigm
       95df301 row.names=F for exploratory summary table
       0f11b55 Working directory madness fixed
       4de25fc Trying out rmarkdown::render() instead
       22d5cc6 Trying html output
       2c3caa8 knitr would set the wd to input document location...
       c624897 dir.create problems
       960c5e5 quotes missing
       01e5ea5 no message
       94b1e9b EOF is simpler
       5742f00 Had forgotten the module calls
       65e698e knit for exploratory + other changes
       39f86fe Created Platanus install script

  Gary Leveque <gary.leveque@gmail.com>      2 commits

       248c455 pacbio_assembly.base.ini edited online with Bitbucket
       93cfad1 pacbio_assembly.base.ini edited online with Bitbucket --changed smrtanalysis_version = 2.3.0.140936.p2 to: smrtanalysis_version = 2.3.0.140936.p4

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      1 commits

       2d582b4 changes necessary for bacterial RNAseq using STAR --see BFXDEV-449

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      52 commits

       f7e38cb BFXDEV-405 correct single library issue in insilico read normalization
       a745799 Merged in rnaseq_denovo_assembly_new_features (pull request #10)
       b1ecdcc BFXDEV-397 resolved rebase conflict resources/modules/verifyBamID.sh
       23a6e08 BFXDEV-397 PRJBFX-1187 added genome and software installs for verifyBamID, resolved issues from pull request #10
       296aab1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       91707a5 BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features removed dev modules
       f80579a BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features
       8ebbedb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d492790 BFXDEV-397 added spaces to table captions to goseq and dge filtered outputs
       c03656d BFXDEV-397 added DGE and GOseq using filtereg contigs (based on trinotate annotations)
       c22ce84 README.md edited online with Bitbucket
       706303f BFXDEV-439 force installation or rmarkdown and knitr
       ab8526c BFXDEV-439 resolved rebase conflicts rnaseq_denovo_assembly_new_features
       716c646 BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       028430f BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       1905e38 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, commit new libraries
       13cefa9 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report, regenerated rnaseq_de_novo pipeline README markdown page
       8d5737a BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report
       b40353b BFXDEV-397 correct dependency problem in exploratory analysis using filtered isoforms
       59fdc4e BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       107fe2a BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       b85f60b BFXDEV-396 added exploratory analysis using filtered transcripts
       17e2143 BFXDEV-423 BFXDEV-399 BFXDEV-397 corrected blastx, abundance estimates to matrix and generated tabular and fasta filtered assembly using python SeqIO
       de22cfa BFXDEV-432 define module picard rnaseq_denovo_assembly base ini
       3820738 Merge branch 'rnaseq_denovo_assembly_new_features' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       f249f9c BFXDEV-397 rnaseq_denovo_assembly.py and rnaseq_denovo_assembly.base.ini rebased from master
       9639c9f BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       1922f92 BFXDEV-397 tested filtered trinity output using real data (1 sample)
       afad437 BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d223c74 detected bugs during tests
       6d0e428  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       f91e427 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       c4597fa BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       b2afc9e detected bugs during tests
       bd91a58  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       2939c0e BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       5ef76ae BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       641a653 BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       cb1e022 BFXDEV-397 tested filtered trinity output using real data (1 sample)
       d9a9c46 BFXDEV-397 merging conflicting versions of rnaseq_denovo_assembly.base.ini and rnaseq_denovo_assembly.py files
       5c5b83e Merge branch 'rnaseq_denovo_assembly_new_features', remote branch 'origin' into rnaseq_denovo_assembly_new_features
       62acde9 BFXDEV-397 fixing differences when rebasing from master
       33d0dd0 BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d3e9ff6 detected bugs during tests
       73ec341  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       8a0bb01 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       4d099e5 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       126102e detected bugs during tests
       f9d902a  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       0cb2d53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       5579bc4 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       ff1f6da BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo

  lletourn <louis.letourneau@mail.mcgill.ca>      18 commits

       74503b2 Merge branch 'master' into highCoverageVariants
       bfdf023 fixed io_buffer default
       8b2cfcc BFXDEV-392 First implementation of high depth calling
       d124958 BFXDEV-370 Fixed merging and output naming bugs
       7381a03 Merge branch 'tumorPair' into highCoverageVariants
       da6e1fb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fffb12b BFXDEV-392 need varscan for high coverage
       4bdd553 Merge branch 'master' into tumorPair
       8ae0dd3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       23c3e59 set default ppn for igv to 2
       bba9c56 BFXDEV-379 remove -gpfs from the ini files
       6cf6ac4 Version bump to 2.2.0-beta
       8e1345d BFXDEV-370 added merging step
       27983c2 Merge branch 'master' into tumorPair
       61e35f3 BFXDEV-370 Added indels and COSMIC
       ffe5e6d BFXDEV-370 fixed scalpel script, added LD_LIB_PATH
       8c0cfc7 BFXDEV-370 Added the scalpel module
       7452be9 bvatools version bump

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      105 commits

       db54b51 resolving conflict
       e7dc251 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       3a2c22f resolving conflict
       c554b5d resolving conflict
       93a55f8 RNAseqDN - update module and remove trinity version check - BFXDEV-490
       48fafb9 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       2be21e1 Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       4103b9f DNASEQ - update SnpEff command for versuion 4.2 - BFXDEV-490
       8e1bbc2 DNASEQ - support single end libray for metrics and dnasample metrics steps & allow filter_nstrech to use non recalibrated compress data - BFXDEV-503 BFXDEV-505
       4e0cd6c resolving conflict
       228c23a resolving conflict
       6e42b3d resolving conflict
       fea21f4 RNAseqDN - update module and remove trinity version check - BFXDEV-490
       9cccf23 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       50f9346 Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       71c233c DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       3a3c119 resolving merging conflict
       a0491a4 RNAseqDN - update module and remove trinity version check - BFXDEV-490
       ca27f3f RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       5b7e920 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7481dc7 Rnaseq_denovo_assembly - correct single library type bugs - BFXDEV-405
       b041650 DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       132e479 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ffcee7d DNASEQ - support single end libray for metrics and dnasample metrics steps - BFXDEV-503
       3f22b9d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       10e4a0a updating modules (test on abacus human done) - BFXDEV-490
       5de58e2 rnaseq - support fr-stranded single end library in wiggle tracks and metrics - BFXDEV-499
       106788a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9b97e88 rnaseq - make  wiggle tracks working for non UCSC genome ref - BFXDEV-498
       9bb4b27 correct the version of pandoc - BFXDEV-490
       e3f47cd update chipseq ini to latest version of modules - BFXDEV-490
       ae112c0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f7b90b5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0e0ee5 update human genome sh files - BFXDEV490
       e3c6325 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       acc1d01 module - remove type in macs2 that make the module file invisible to the module system - BFXDEV-490
       4c413a2 module - add gemini and bwakit modules - BFXDEV-490
       c9a2f31 module - add star indexes for read length 75 and 150 bp - BFXDEV-490
       14464b1 module - adding a step print for the module creation - BFXDEV-490
       39f89bc module - add gemini and bwakit ; update bcftools htslib java picard prinseq-lite rnaseqc samtools - BFXDEV-490
       b50d111 BFXDEV-490 - updating modules gnuplot ; hmmer; igvtools; macs2; pandoc; python_lib
       61770a4 fixing conflict between master and highCoverageVariants branches
       a3491eb BFXDEV-490 - updating modules bedtools; blast; gatk; ucsc
       bd023c2 bfx/gatk.py - removing conflict to allow gatk_variants branch to be merged
       292d51c Ressource - bump vcftools version to 0.1.14 - BFXDEV-489
       13aa925 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5b8a705 genomes - add rat Rnor_6 - BFXDEV-487
       9ef56ce genomes - install genome by default in dev and in prod if argument MUGQIC_INSTALL_HOME is given - BFXDEV-485
       4ac216d ressource - fix version of matplot lib to 1.4.3 for Qiime compatibility - BFXDEV-483
       28175da gatk_variants - danseq -  add mark_duplicates cleaning files - BFXDEV-471
       1c06f0c gatk_variant - module - bump mugqic_tools install script to 2.1.3 - BFXDEV-484
       62c067e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       da7245b ressource - update python to 2.7.10 and add scikit-bio lib - BFXDEV-483
       7c774e6 remove tumor_pair from highCoverage
       5a20f1f DNAseq - add cleaning list to the new steps - BFXDEV-471
       9266e03 DNAseq - change destionation path in the copy of SNV metrics zip file - BFXDEV-473
       71b439e Pacbio assembly - modify whitelist option usage to do not generate additional ini and xml - BFXDEV-456
       1f4eb9d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be62216 Bump  module to the new smartanalysis patch 5
       15b99f4 gatk_variants - add baits intervals file for picard_hs_metrics - BFXDEV-467
       ad0a5c9 ingstats.py - Adding Array callrate check (Missing or Low) and add support for different manifest version through the -m option - BFXDEV-445 - BFXDEV-466 - BFXDEV-451
       75379ef DNAseq - synchronize code and ini and correct runnning issue - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       7a52905 DNAseq - remove some smalls dev typos and bugs - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       a19894c PacBio - add require=False to the whitelist param - BFXDEV-456
       e8223e6 PacBio - incorpore whitelist option during intial filtering - BFXDEV-456
       53b2fbf PacBio - add addtionnal filtering config xml file for  whitelisted assembly - BFXDEV-456
       064ae2e PacBio - add specific ini file for whitelisted assembly - BFXDEV-456
       8e298ef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7769b77 update bwa install script to the newer version
       7261373 change star index place in resources/genomes/install_genome.sh
       0904455 update install_genome.sh to use the coirrect version of R_packages
       a3e0340 bump bos_taurus install script to ensembl 80
       c57be72 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       94223ff bump bos_taurus install script to ensembl 81
       17bea25 make ignstats poiting on nanuq (not QC nanuq)
       382c952 Merge remote-tracking branch 'origin/processingReport'
       051e378 RNAseq - repair missing dependency - BFXDEV-427
       26d5d41 debug paired_tumor scalpel vcf merging
       d2c011b removed merge conflicts
       29a319b DNAseq - regenarate updated README.md
       85017c1 DNAseq - change fixemate description text: picard -> bvatools
       ba184cc High_coverage - needs dict2BEDs.py which is only in dev version of mugqic_tools - modifiy the ini file
       3d88798 DNAseq - change post fixmate sorting to generate and indexed bam - BFXDEV-421
       83e7e16 modify igvtools lib to allow modifiy paramter through the ini file - BFXDEV-419
       d5d86bf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       71c3d26 DNAseq - change coverage metrics report to remove the CCDS part of the table - BFXDEV-381
       0de4e65 rnaseq add bam_hard_clip step : Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature - BFXDEV-377 + change rRNA metrics output to avoid name conflict - BFXDEV-401 - BFXDEV-407
       8a63681 update python and python lib to include Qiime ligth version installation - BFXDEV-412
       8bb30ee bump trinity install script to version 2.0.6
       1b2a34d rnaseq - remove tophat/botwie commented section and modules
       f53ed32 rnaseq_denovo_assembly - correct single library issue and add in comment the correponding change in the base.ini - BFXDEV-405
       d37b2c8 nanuq2mugqic_pipelines.py - support nanuq group info for Miseq/Hiseq/Pacbio technology - BFXDEV-418
       b03537a update install module general script
       2d1dbe4 bump smrtanalysis module to 2.3.0 patch 4
       27bad46 update install module general script
       006f59c update install module general script
       08acf20 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcb9589 bump smrtanalysis module to 2.3.0 patch 4
       ba4b9f1 update bvatools module for release 1.6
       c72acd4 core/pipeline.py remove duplicates sections in reports - BFXDEV-388
       4cb0ba7 ChipSeq - remove expected input files for broad peaks when generating homer annotation report - BFXDEV-393
       797c968 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3a54c85 modify resources/modules/smrtanalysis.sh for patch 3
       ce104cd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d309b7d RESSOURCES - R_bioconductor - Addthe package 'circlize' to the list of R package to automatically install

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      32 commits

       a57fd6d RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       129d564 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       4695881 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       226c59b RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       1cbb0f1 chipseq - update module -  BFXDEV-490
       9e24486 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       501c3ac DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       bc5920a RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       849d732 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       db6c64a chipseq - update module -  BFXDEV-490
       65e3861 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       e1637bf DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       c9e8db0 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       0d077aa RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       32d1d5b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4c43bb2 chipseq - update module -  BFXDEV-490
       7a68160 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       7e2f61b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a95aa1 DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       58f9477 BFXDEV-465 - correct GATK DoC small bug to support when no bed fileis given in the readset file
       269b60f DNAseq - Adding variant recalibration BFXDEV-436
       557f75a DNAseq - gatk DoC will use the bed file as intervalls if the bed file is in the readset shett  - BFXDEV-465
       ce3014a DNAseq - implement haplotype caller and mpileup annotation and filtering using old foinction as background - BFXDEV-463
       2aa0008 DNAseq - create new pipeline steps - BFXDEV-463
       f46876a DNAseq - Add gvcf Combining the set of sample and genotyping - BFXDEV-440
       157b7e9 DNAseq - starting to implement GATK gvcf merging - BFXDEV-440
       3c407b2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       74c6c2a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66cfdc6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcfdd2c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c83d71c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1b17411 bump pacbio module to patch 4 -  BFXDEV-415

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      6 commits

       5db43a8 Merged in irp_end_fastq_notification (pull request #14)
       de7d7f8 Merged in irp_bcl2fastq2 (pull request #13)
       d3acddd Merged in illumina_run_processing_sprint (pull request #11)
       82747c8 Merged in pacBio_whitelist (pull request #7)
       e16e90d Merged in ignstats (pull request #8)
       8e2cf4c rnaseq.base.ini increase star io limit to 4G by default

  mmichaud <marc.michaud@mail.mcgill.ca>      62 commits

       4d820a7 BFXDEV-504 - Notify nanuq when fastqs are completed for a lane (fix url) - illuminaRunProcessing
       7ab569e BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       cdafb23 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       fd8579e BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       7282b30 BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       2d9a896 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       e0be331 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       6e04c62 BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       34ca19a BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       09ff700 Increase resources for qc.
       3dfce03 BFXDEV-488 Illumina Run Processing: Use bcl2fastq2. Fix for dual-indexing.
       5e8133e Use 2 core for rsync.
       48a4bd8 Merge branch 'master' into irp_bcl2fastq2
       a9560d4 Fine tune cluster settings according to the history of the jobs of the last three months.
       b917ee5 BFXDEV-413 Fix code according to code-review.
       a8fd8bf BFXDEV-488 Illumina Run Processing: Use bcl2fastq2.
       be6ec3f Merge branch 'master' into illumina_run_processing_sprint
       e0c3e20 BFXDEV-468 VerifyBamId: Fix job name when not using bed file.
       ec19be9 BFXDEV-468 VerifyBamId: Run even we don't have a bed file.
       1ffafbb BFXDEV-468 VerifyBamId: Use a version supporting "chr" in chromosome name.
       dec09d8 BFXDEV-468 VerifyBamId: Fix output for nanuq.
       eed1cc1 BFXDEV-468 VerifyBamId: Shorter job name.
       750495f BFXDEV-468 VerifyBamId: Properly concat jobs
       ce01a76 BFXDEV-468 VerifyBamId: Nanuq friendly output.
       9d5fa0e BFXDEV-468 VerifyBamId: Don't use a compressed VCF + Nanuq friendly output.
       c4bd81b BFXDEV-468 VerifyBamId: Add missing configuration.
       72c77d2 BFXDEV-468 VerifyBamId: Fix annotation file name.
       0001099 BFXDEV-468 VerifyBamId: Fix wrong genome for sample.
       15a3a2b BFXDEV-468 VerifyBamId: Changes to the annotation filename.
       7fc91bd BFXDEV-468 Don't run verifyBamId when there is no "dbnfp_af_field" on the genome.
       d747ef6 Illumina Run Processing: Fix barcode counter path with cvmfs
       68c0089 BFXDEV-468 IRP - VerifyBamId: Optional dbsnp_version and dbnsfp_af_field in genome ini file.
       5a3fd82 Use cvmfs blast and bwa.
       75af739 BFXDEV-468 Add verifyBamID in Illumina Run Processing
       94ccde0 Merge branch 'master' into illumina_run_processing_sprint
       7d6cf51 BFXDEV-447 Add rRNA estimate using silva blast database. (Fix change output file name)
       dd9f3f7 BFXDEV-464 Don't use an hardcoded genome list for the alignment: parse the value from the sample sheet and validate that we have the genome on disk.
       015435e BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name (fix redo of pipeline always restarting the STAR align).
       f0ecd63 BFXDEV-447 Add rRNA estimate using silva blast database. (Fix error when creating result file)
       2d39be9 BFXDEV-413 Add Estimated Genome Size from the nanuq CSV.
       c73b3fc BFXDEV-447 Add rRNA estimate using silva blast database.
       fdcb45e BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name.
       1819bf2 Add sample sheets examples and add description about "Genomic Database"
       be10d09 BFXDEV-457 BFXDEV-459 Update documentation.
       b337a06 Don't warn when skipping alignment on a sample without "Genomic Database"
       a73ca70 BFXDEV-459 Merge start_copy and copy steps.
       583237b BFXDEV-457 Add the possibility to regroup all the md5 jobs into one (add "one_job=1" in "[md5]" config section).
       56cdf2d BFXDEV-458 Fix mask calculation problem on lanes with nextera/truseq samples.
       041d4e8 More useful usage message.
       3cd065c Merge branch 'master' into irp_genomic_database
       d90c137 Merge branch 'master' into irp_genomic_database
       2f0ef65 Illumina Run Processing: Increase memory for BVATool DoC and Readsqc
       a864f31 BFXDEV-369 Illumina Run Processing: Fix GRCm38 genome detection.
       7eae59d BFXDEV-369 Illumina Run Processing: Use the reference specified in the request form
       64fb218 Increase processors per node for Fastq and bwa alignment. Ignore setfacl exit code for the copy step.
       fe44577 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0) - Fix passFail threshold according to the new value.
       cafa62c BFXDEV-417 IGNStats: Output Array Call Rate as ratio, not percentage (0.95 not 95.0%)
       eb810e0 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0)
       3724583 BFXDEV-417 Modifiy IGN stats parser to allow the upload of data to a nanuq server.
       ce52854 Use 13 core for BWA job to avoid excessive memory usage.
       8505fab BFXDEV-385 Fix bed files handling in Illumina Run Processing.
       29d5c0f BFXDEV-380 Change permissions of the run files on the destination.

  noreply@clumeq.ca <libuser@lg-1r14-n01.guillimin.clumeq.ca>      15 commits

       ad55029 genomes - updated bash script for human genome GRCh37 and GRCh38  - BFXDEV-490
       4d653bf Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       d2fd863 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       e4b6339 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       1c79f7e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       36b260f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       5655027 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       7e8870a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0a73d79 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eb30107 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8c32446 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0b9b70c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0233d70 mugqic_tools is now up to date && resources/modules/dev/epacts.sh has been removed
       64c821d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a65de9e temporary modify install_genome.sh to run the job in batch (job submission is blocked

  ptranvan <patrick.tranvan@mail.mcgill.ca>      1 commits

       2a1917b Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      4 commits

       f184966 Addition of htslib module and resolved pull request comments
       457c369 addition of vcf preprocessing, snpeff, gemini and ini adjustments
       09ad7f1 addition of vcf preprocessing, snpeff, gemini and ini adjustments
       05f05b2 Merge branch 'highCoverageVariants' of https://bitbucket.org/mugqic/mugqic_pipelines

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      2 commits

       60fd230 Merge branch 'highCoverageVariants' of bitbucket.org:mugqic/mugqic_pipelines into highCoverageVariants
       76802ee Updates to pairedTumor: addition of CombineVariants

2.1.1        Mon Apr 13 22:23:46 2015 -0400        172 commits

  Francois Lefebvre <lefebvrf@gmail.com>      4 commits

       65710dc QUAST and Minia dev install scripts
       5b404f3 Updated kmergenie version in install script and moved to dev
       6db368e Added NxTrim mugqic_dev install script
       45b9fb0 Updated a bunch of module dev install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      109 commits

       4db053f Added Gorilla_gorilla.gorGor3.sh in resources/genomes/old/
       deb0b82 Changed report job names with '_' instead of '.' to avoid scheduler cpu setting > 1
       4d8b38a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b471f6f Fixed bug in RNA-Seq De Novo Assembly : use RSEM plugin in Trinity instead of external one
       6c79dd4 README.md edited online with Bitbucket
       2e53a76 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       455a50c Fixed missing pandoc module in rnaseq
       0ffb7ea Added variables total and average read length in pacbio assembly stats report + changed locale from fr_FR to en_CA
       70150e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       121faa0 Added sample_cutoff_mer_size to pacbio report job names
       fcb6b89 README.md edited online with Bitbucket
       272830b README.md edited online with Bitbucket
       89c55cb README.md edited online with Bitbucket
       8f6eeb7 Fixed cpu bug for blastp_transdecoder_uniprot job in RNA-Seq De Novo assembly
       cae10e5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       28eb954 Minor comment updates
       2c4a463 README.md edited online with Bitbucket
       c4262fa README.md edited online with Bitbucket
       3f664a5 README.md edited online with Bitbucket
       dc8ec12 README.md edited online with Bitbucket
       c60d742 README.md edited online with Bitbucket
       b38e195 README.md edited online with Bitbucket
       6d4206c README.md edited online with Bitbucket
       4a6982c README.md edited online with Bitbucket
       39ee0b9 README.md edited online with Bitbucket
       1f13501 README.md edited online with Bitbucket
       0376786 README.md edited online with Bitbucket
       43d0ac6 README.md edited online with Bitbucket
       be50cd0 README.md edited online with Bitbucket
       f85506e README.md edited online with Bitbucket
       5e335b4 Added GNU Lesser General Public License (LGPL) to MUGQIC Pipelines
       0b79b46  BFXDEV-59 Updated ChIP-Seq report redesign with sample metrics without trimming
       ec88ea5 Fixed bug cluster_cpu for blastx_trinity_uniprot
       5962a9f Removed UniRef BLAST in RNA-Seq De Novo Assembly since it is too long; factorized differential expression code; renamed design variable into contrast
       611223e Adjusted Trinity butterfly cluster resources for RNA-Seq De Novo Assembly on abacus
       371d60c  BFXDEV-59 Completed ChIP-Seq report redesign
       7eedd05 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       01a41aa  BFXDEV-59 ChIP-Seq report redesign
       f8bf89f BFXDEV-59 Completed PacBio Assembly report redesign
       479ea0b BFXDEV-59 Differential expression RNA-Seq De Novo Assembly report redesign commit
       76c9bc3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       c7a9036 BFXDEV-59 More RNA-Seq De Novo Assembly report redesign commit
       ac2daf3 README.md edited online with Bitbucket
       e7f156a README.md edited online with Bitbucket
       3caf504 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       6771f03 BFXDEV-59 First RNA-Seq De Novo Assembly report redesign commit
       5c858ad Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5d82ab0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       94e9930 Decreased rnaseq cufflinks default pmem to 2700 for guillimin
       79ea01b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       80cfde8 Increased rnaseq cufflinks default pmem to 3700 for guillimin
       9e655b9 Fixed BED file abspath in bvatools
       c183fad BFXDEV-59 Completed RNA-Seq report redesign
       61676eb BFXDEV-59 Added metrics steps for RNA-Seq report
       e760ac7 BFXDEV-59 Fixed merge conflict
       fbedae5 BFXDEV-59 More and more commit for partial HTML report
       7a95c84 Increased dnaseq compute_effects ram
       c6c7ead BFXDEV-59 Even more commit for RNA-Seq report
       75b2a7e BFXDEV-59 Even more commit for RNA-Seq report
       5da9885 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       0463d00 BFXDEV-348 Fixed with readset name = <sample>.<library_barcode>.<run>.<lane>
       3072da1 BFXDEV-59 Even more commit for RNA-Seq report
       57f2aef BFXDEV-59 More commit for RNA-Seq report
       ad7e670 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be92407 Updated mugqic_tools to version 2.1.0
       672b681 BFXDEV-59 First commit for RNA-Seq report
       f194ec3 BFXDEV-59 DNA-Seq report minor fix
       7f98cbf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       40a04d0 BFXDEV-59 DNA-Seq report redesign done + fix
       b4b408f BFXDEV-59 DNA-Seq report redesign done
       426e878 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       a0b4223 Minor doc fix
       bc1d662 BFXDEV-59 Even more report redesign commit
       09de8ac BFXDEV-59 More report redesign commit
       91e8e95 BFXDEV-59 First report redesign commit
       77f66c4 Added UCSC genomes in install_all_genomes.sh
       aa6a8e0 Create genome rrna.fa with grep -i 'rrna'... + remove variation sequences containing '.' in Ensembl vcf.gz which make GATK crash
       525faa8 BFXDEV-295 Updated mugqic_R_packages to 1.0.3 for RNA-Seq De Novo pipeline
       da8253e BFXDEV-295 minor fix for RNA-Seq De Novo pipeline
       e161ad9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4d177be Updated Python pysam to version 0.8.2
       e7ca27f README.md edited online with Bitbucket
       48073ca README.md edited online with Bitbucket
       70a0304 README.md edited online with Bitbucket
       1e24257 Increased cores for homer_annotate_peaks in chipseq
       3c2c24d Fix new section names for blast on uniprot in RNA-Seq De Novo pipeline
       e699c64 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       aeedee4 Added chipseq in PATH of module mugqic_pipelines
       7ce201c Added ccds filepath in Rattus_norvegicus.Rnor_5.0.ini
       7b08143 BFXDEV-295 Moved trinotate Pfam + UniProt DBs in /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod/genomes/[blast|pfam]_db/ for RNA-Seq De Novo pipeline
       c47a519 BFXDEV-295 Update RNA-Seq De Novo pipeline with modules in prod and trinotate updated
       33812be Minor mammouth adjustment regarding increased ram and core in snp_effect job for DNA-Seq Pipeline
       9689e6d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c345049 Increased ram and core in snp_effect job for DNA-Seq Pipeline
       8743bbb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f5b422 BFXDEV-295 Updated RNA-Seq De Novo assembly pipeline with Trinity 2.0 and Trinotate 2.0
       b484ff7 README.md edited online with Bitbucket
       4b19c20 README.md edited online with Bitbucket
       1b0f09c Minor fix in python lib install
       9180612 Separated python install and python lib install + minor weblogo install update
       e06560d Minor fix in Perl lib install
       b202849 Fixed missing bvatools.depth_of_coverage other_options + snpsift_annotate module_snpeff=mugqic/snpEff/4.0 + minor uppercased '.insert_size_Histogram.pdf' for picard.collect_multiple_metrics output file in dnaseq pipeline
       f8b35eb Fixed missing 'nodes=1' in pacbio_assembly.guillimin.ini smrtanalysis_run_ca
       e7d4a0e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       62a0e00 Updated bedtools to version 2.22.1
       9328bc9 Standardized perl module install and moved CPAN libs in a different script
       0fb4c64 BFXDEV-335 Removed adapter FASTA files except adapters-fluidigm.fa
       4aba5c1 BFXDEV-335 Create adapter FASTA from readset file for Trimmomatic, if not defined in config file
       ba78748 Version bump to 2.1.1-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      14 commits

       1792f2c Version bump to 2.1.1
       8810822 BFXDEV-375 Fixed ram sorting issues when using star
       e9fcf54 Merge branch 'master' into rna_metrics
       a7e3465 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ba6662 BFXDEV-368 IGN script to extract stats
       8f69140 Fixed code for tools that need to be downloaded manually like gatk
       f6ee04e Updated gatk
       4622df3 Updated samtools, bcftools, htslib
       6eac315 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a44561f BFXDEV-351 removed md5 from markdup and added it to recal
       8abc6cc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1db78f9 BFXDEV-355 Removed CCDS from the options BFXDEV-356 added compression to haplotype caller output
       9a0a698 BFXDEV-351 Removed MD5 from markdup, added it to recalibration
       b33195c BFXDEV-346 Split jobs in a more uniform way

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       05937d6 RNA-seq - rna_metrics : test are ok; new files annoation files are created ; point to the production assembly folder -  BFXDEV-345
       b2c948b ressource -Genome install script add the creation of the ref_flat file format of the annotation from the gtf and correct a path in the GRCh38 file (All.vcf.gz) - BFXDEV-374
       7d8f72f RNAseq- update module version not found in CVMFS (bowtie, bvatools, mugqic_tools) - BFXDEV-373
       e57ed5e RNAseq- remove redundant step in the step initialization - BFXDEV-371
       3138867 release new version of resources/modules/mugqic_tools.sh
       9aa58f7 RNAseq - rna_metrics - finish the rRNA and rna metrcis modifications - BFXDEV-345
       cd5f239 RNAseq -rna_metrics -fit the new bam2fq from bvatools_dev
       d2048eb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rna_metrics
       d651fee RNAseq - rRNA metrics add bvatools | bwa | picard && rrnaBMAcount.py

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       70099f3 rnaseq - include correlation matrix manual estimation using an utils.py new function
       fd53801 RNAmetrics - include ini_section argument to some bfs picard function to allow to use them several time with different parameter in the same ini
       e7c5a01 rnaMetrics - remove confilct pipelines/rnaseq/rnaseq.base.ini
       9ad7fbb rnaMetrics - update pipelines/rnaseq/rnaseq.py pipelines/rnaseq/rnaseq.base.ini
       9d3b87d RNAseq - metrics update ini
       1913bc1 RNAseq - remove conflict
       deca058 RNAseq - metrics RNA - update base ini
       2fc2006 RNAseq - remove rnaseqc; add picard_rna_metrics ; partial add estimate_ribosomal_rna - BFXDEV-345

  mmichaud <marc.michaud@mail.mcgill.ca>      27 commits

       2a821a4 bcl2fastq module name change.
       00b481f Merge branch 'master' into irp_rna_metrics
       c2c33a3 MPS-1740 Use updated production genomes.
       5f3c74e MPS-1740 Use released version 1.5 of bvaTools
       d191969 Increase STAR sort memory
       dd936d0 BFXDEV-363 IRP: Don't copy phasing and matrix files.
       ae9713b BFXDEV-353 IRP: Standardize job name
       e83e25f Merge branch 'master' into irp_rna_metrics
       87aaf0d Merging rna_metrics on irp-rna_metrics
       b0dff3b BFXDEV-353 Use dev version of mugqic tools
       f94c135 Merge branch 'master' into irp_rna_metrics
       75398a3 BFXDEV-353 Use new version of rRNABAMcounter.
       5f443d4 BFXDEV-353 Use a nanuq friendly name for the rRNA metrics file
       ea661af BFXDEV-353 Add the rnaseqc '-rRNA' option to set an optional ribosomal rna interval list. Setting an empty file will skip the rRNA count.
       91c0d87 Simplify illumina run processing RG tag logic. As library, run and lane are mandatory there is no need to validate that they exists.
       4abaa7b BFXDEV-353 Fix bwa other_options
       81e7423 BFXDEV-353 Use DEV versions of bvatools and genome. Update mugqic_tools to 2.1.0.
       79a02a7 BFXDEV-353 Update rna-seq metrics according to rna-seq pipeline
       081c6b8 Merge branch 'master' into irp_rna_metrics
       fcc8392 Code format
       d62f855 Merge branch 'master' into irp_rna_metrics
       0bc4cbe Add missing thread parameter for bvatools_depth_of_coverage.
       221eaae BFXDEV-339 Use rrna file in rnaseqc
       8e95c99 BFXDEV-339 Use rrna file as ribosomal annotation (instead of ncrna)
       5efa104 BFXDEV-338 New Nanuq MPS run association
       7aab26f BFXDEV-338 New Nanuq MPS run association
       6e6de07 BFXDEV-338 Add run_id in all commands available parameters

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       559973d correction other_options

2.1.0        Wed Feb 4 14:40:16 2015 -0500        123 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       eefcc47 Added qualimap installa script + updated population
       d9397b0 Changes to sailfisj , R install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      104 commits

       09344ee Version bump to 2.1.0
       b3e8561 Changed procs= to nodes=1:ppn= in rnaseq and rnaseq de novo guillimin config
       a63500d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a195837 Added libgd in PacBio mammouth config
       a863cdd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e13ac41 Added pacbio_assembly.mammouth.ini + removed unnecessary wgs module
       1cb2ea5 Added pacbio_assembly.guillimin.ini + adjusted guillimin cluster settings
       b9f1264 README.md edited online with Bitbucket
       ff58b26 Removed 'daemon' from scheduler options
       2d4459d BFXDEV-292 Removed optional trimming dependencies for report in chipseq pipeline
       6b9cdad Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f40f39a BFXDEV-292 Fix metrics and report dependencies in chipseq pipeline
       c1fb324 README.md edited online with Bitbucket
       d8ade06 Removed tmp_dir validation since some directories are available on exec nodes but not on login nodes
       44c931e Removed memtime from PacBio pipeline + updated smrtanalysis to version 2.3.0 patch 2
       698d3f7 BFXDEV-292 Added cleaning in chipseq pipeline
       cac6305 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       086ee05 Updated smrtanalysis module 2.3.0 with patch 2
       0c6d970 README.md edited online with Bitbucket
       f808abe Updated ChIP-Seq README + main READMEs
       e789ffd Added BiSNP module install script
       26a0629 Update mugqic_tools to 2.0.3 in rnaseq config
       7af099d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d10a058 BFXDEV-292 Fix metrics tmp file preventing job to be up to date + module ImageMagick to use convert command on mammouth
       367f86e Removed '\' before /ltmp/fillon in mammouth config files
       0c1988f Added picard tmp_dir type validation
       cf82065 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0884bf5 BFXDEV-292 Added config files for all clusters in chipseq pipeline
       2a57215 Updated mugqic_tools install to version 2.0.3
       387c225 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       5017938 BFXDEV-292 Last steps implementation in chipseq pipeline
       52d1f78 BFXDEV-292 Added UCSC genomes hg19, mm9, mm10, rn5
       751e701 BFXDEV-319 Adjusted cluster settings in rnaseq.mammouth.ini + .ini absolute path for report
       edc2651 BFXDEV-35 Standardized memtime module install
       90f4aa5 Minor update in README-release.txt
       ab65f5d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       f37bed0 BFXDEV-292 Set peak annotation with internal homer genome in chipseq pipeline (incomplete)
       9b35602 BFXDEV-35 Moved wgs-assembler module install in dev
       77859cc BFXDEV-35 Standardized vcftool module install
       a07ea0b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       2ee4315 BFXDEV-35 Standardized ucsc module install
       284138b BFXDEV-35 Standardized trinotate module install
       9f31173 BFXDEV-35 Standardized trimmomatic module install
       1fab82a BFXDEV-35 Standardized tophat module install
       00b9355 BFXDEV-35 Standardized tmhmm module install
       2917504 BFXDEV-35 Standardized tabix module install
       6c1ec5b BFXDEV-35 Standardized star module install
       5554746 BFXDEV-35 Standardized snpEff module install
       25f53a1 BFXDEV-35 Fix smrtanalysis archive exec permission
       275894f BFXDEV-35 Standardized smrtanalysis module install
       6f3b27a BFXDEV-35 Standardized signalp module install
       a335bfa BFXDEV-35 Standardized samtools module install
       7cfda90 BFXDEV-319 Added rnaseq.mammouth.ini
       dac40de BFXDEV-319 Removed java option -Dsamjdk.use_async_io=true in all .ini files
       e394bf7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       705f2c6 BFXDEV-292 Adjusted cluster settings in chipseq pipeline
       a87203d README.md edited online with Bitbucket
       2524cb9 README.md edited online with Bitbucket
       af39fb4 README.md edited online with Bitbucket
       e180a83 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0daae6e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       93c1a28 Fixed bam sub arguments in dnaseq bwa_mem_picard_sort_sam
       934b42f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       a5931e4 BFXDEV-35 Standardized rsem module install
       dbaebd1 Moved repeatmasker module install script in dev
       0c30ff4 BFXDEV-292 Minor docstring change in chipseq pipeline
       dc0b988 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ef6c09e BFXDEV-292 Added GTF for homer_annotate_peaks in chipseq pipeline
       26894dc README.md edited online with Bitbucket
       c7cdd1a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       4f0b085 BFXDEV-292 Fixed input/output files bugs in chipseq pipeline
       4c1cf7d BFXDEV-35 Standardized weblogo module install
       bf36cc0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       7952a5f BFXDEV-292 homer_find_motifs_genome and annotation_graphs steps in chipseq pipeline
       b29cb3d BFXDEV-35 Standardized rnaseqc module install
       05208ed BFXDEV-35 Standardized rnammer module install
       e33776a BFXDEV-35 Standardized more prinseq-lite module install
       fcc8bc9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d15d9b0 BFXDEV-35 Minor fix in AMOS module
       07bba1e BFXDEV-35 Fixed path bugs in MUMmer and AMOS module install
       99363ee BFXDEV-35 Standardized MUMmer module install
       423c3c7 BFXDEV-35 Minor fix in module install template
       c00bdfc BFXDEV-35 Standardized java module install
       d185a50 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       8bed5d9 BFXDEV-292 homer_annotate_peaks step in chipseq pipeline
       0bca196 BFXDEV-35 Standardized igvtools module install
       1581fc9 BFXDEV-35 Standardized hmmer module install
       2acdc49 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       cddbce7 BFXDEV-292 More macs callpeak in chipseq pipeline
       ed7cf4e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b2d02b1 Fixed macs2 with generic shebang
       cf61d7e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       bbe5703 BFXDEV-292 Fixed missing samtools module for MACS2
       ff1a4ec Standardized gnuplot module install
       dcdc211 Standardized exonerate module install
       5c47746 Minor change in README-release.txt
       076c02b Version bump to 2.1.0-beta
       953894c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       1090d38 More chipseq deelopment
       0841b02 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ace0e30 BFXDEV-292 Beginning of macs2_callpeak step in chipseq
       b404e22 Added qc_plots_R step in chipseq
       43e70dc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       51f3444 BFXDEV-292 First draft of chipseq pipeline

  lletourn <louis.letourneau@mail.mcgill.ca>      3 commits

       9f9f60c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d2d65f0 BFXDEV-327 Used only one thread for haplotypecaller because of a race condition
       639f650 BFXDEV-327 Used only one thread for haplotypecaller because of a race condition

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       7ec0344 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1905eb6 COMMON - add fastq2 = None in sam_to_fastq pipeline step wehen the read are single - BFXDEV-321

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       12ff20e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       233244b RNAseq - correct stdin input issue of htseq-count - BFXDEV-318
       fa33528 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       47c4cc8 RNAseq - htse-count: pipe samtools view -F 4 output in htseq-count instead of using the bam to remove error due to unmapped reads - BFXDEV-318

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       953de3d BFXDEV-332 Fix depth_of_coverage 'other_options' by removing extra quotes
       abca180 BFXDEV-332 Base the jobs configuration on the perl pipeline
       619c076 BFXDEV-329 Add run and lane number to the alignment and metrics jobs name
       6f4440f BFXDEV-202 Increase RAM for barcode counting
       eac645f BFXDEV-328 Increase mem walltime to 48h (and set ram for SortSam)
       ccb73f3 BFXDEV-331 Run Processing: Use a RG id that is unique across multiple lanes
       c2218da BFXDEV-317 Run processing: Use the old name for the coverage and onTarget metrics
       11ec9e7 BFXDEV-316 Fix errors when using a custom Illumina sheet file

2.0.2        Mon Jan 12 16:56:16 2015 -0500        80 commits

  Joël Fillon <joel.fillon@mcgill.ca>      54 commits

       104778a Version bump to 2.0.2
       bbdf438 Version bump to 2.0.2
       9d984ce Fixed cluster resources in rnaseq_denovo_assembly.mammouth.ini
       ee96fe3 Minor fix in pacbio_tools_split_reads config section name
       8831a3a Even more standardization of module install
       a2c5824 More standardization of module install
       2b334e0 Standardized mugqic_pipelines module install
       395d3ec Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66230bd Standardized mugqic_tools module install
       d9a56ab Standardized trinity module install
       0549be9 Standardized cufflinks module install
       ddbd8f2 Standardized MACS2 module install
       386921e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9ba8d39 Standardized picard module install
       55a5ccb Standardized cd-hit module install
       d693b9a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f05fe1 BFXDEV-310 Updated picard to version 1.123
       a118e0d Standardized bwa module install
       4dc7bb7 Updated bvatools to version 1.4 for dnaseq and puure
       cca258b Standardized bvatools module install
       94910f9 Standardized bowtie2 module install
       3e8f201 Standardized bowtie module install
       1b70a25 In picard_sam_to_fastq, updated skip test if FASTQ column present and BAM missing
       345909a BFXDEV-290 Added job_input_files, job_output_files in JSON export
       f6e623e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       39e9ebf More JSON export for Daemon Scheduler
       7295671 Standardized blast module install
       947111a Standardized bedtools module install
       9325049 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c52d277 Standardized GATK module install
       a69bd82 More module install generalization; first test with prinseq-lite
       d116ebd Minor wget output file fix in prinseq-lite.sh and module install template
       c30ad66 Added prinseq-lite module install script
       b2059f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ca9017 Fixed pb wget --spider in install_genome.sh + various up-to-date checks in Homo sapiens assemblies install
       21e59af BFXDEV-284 In nanuq2mugqic_pipelines.py, use native httplib, retrieve BED files and create adapters file
       e9620ff BFXDEV-290 First draft of daemon scheduler
       b778e5d BFXDEV-74 Finished cleaning in rnaseq denovo assmbly pipeline
       a171581 BFXDEV-74 Started cleaning in rnaseq denovo assmbly pipeline
       c1cee29 Reorganised lib functions more compact
       1cc8e33 BFXDEV-74 Added cleaning in rnaseq pipeline
       13475e2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       84954f3 BFXDEV-74 Added cleaning for dnaseq pipeline + reformatted various bfx modules
       1948d9d README.md edited online with Bitbucket
       848a0ab Fixed bug use lstat instead of stat to check job up-to-date status without following symlinks
       fdfcbdd Minor aesthetic updates on pipelines doctrings + README.md
       da66a41 Generated pipelines README.md from docstrings using --help
       89718b0 Added steps docstrings in RNA-Seq De Novo Assembly Pipeline
       ab5e595 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       232bed5 BFXDEV-296 Added steps docstring for RNA-Seq pipeline
       8d6f9a0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b665c4c Minor release doc update
       bf9afb6 Update mugqic_pipelines module to 2.0.1
       530edcf Version bump to 2.1.0-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       ccb812e Fixed picard installer, reverted back to 123

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       693ea1d PacBIo - add input files to the first pacBio job file (cp job) - BFXDEV-305
       ce7f0b3 Revert "PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305"
       fc6eb26 PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305
       162af5e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2ff273b RNAseq - modify docstrings - BFXDEV-303
       c58b6d3 DNAseq - correct report  dependency (typo) - BFXDEV-304
       ab81bc3 DNAseq - correct report  dependency - BFXDEV-304
       12830d3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b6f4e8f PacBio de novo - to allows running multiples sampes in parralele in the same analysis - BFXDEV-301

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       fcd7134 pull before pushing Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       05b352d  PacbioQC - Fix code vs ini section name - BFXDEV-315
       d085269 RNAseq - correct discrepency in hsteq_count ini calls - link to BFXDEV-312
       6773ef7 pull before pushing
       d5d82f0 RNAseq - correct errounous section header in ini files - BFXDEV-312
       98a456a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fdd92aa RNAseq -update guillimin cluster ini file - BFXDEV-307
       b4a12a5 RNAseq - fix report dependencies - BFXDEV-306

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c14a281 BFXDEV-314 Illumina Run Processing: Output nanuq friendly metrics files
       3141ccc BFXDEV-313 Run processing: Don't depend on the start fastq notification
       5956b71 BFXDEV-310 Use Picard 1.123 with a bug fix in MarkDuplicate
       d018879 BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       793e91c BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       000d043 BFXDEV-309 Use the hiseq qos by default on abacus
       91cfb26 Run processing: Check the library type when the library source is 'library' to determine if the sample is from RNA
       94f196e Run processing: Fix configuration for miseq (copy's destination folder)

2.0.1        Wed Dec 17 09:56:23 2014 -0500        33 commits

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       76a44bc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fe96c58 no message
       0f6a3ec Modified R_Bioconductor script to accommodate older gcc's

  Joël Fillon <joel.fillon@mcgill.ca>      22 commits

       bc72fc8 Version bump to 2.0.1
       5868b85 Updated mugqic_R_packages to 1.0.1 and mugqic_tools to 2.0.2
       3230102 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd49fed Updated /mugqic_tools.sh to 2.0.2
       2b979a2 Fixed dnaseq species_vcf_format_descriptor with absolute path (no /sb/programs/analyste)
       51fef9a BFXDEV-299 Fixed bug GATK realign with unmapped parameter not skipping other chromosomes
       6890d45 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9a53a11 BFXDEV-296 Added DNA-Seq step docstring documentation
       a820bbc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4475a6b BFXDEV-296 Added detailed --help option with output in markdown format
       cb49e7b Added detailed --help option with output in markdown format
       bf7bb0a More PacBio README.md update
       cd04651 README.md edited online with Bitbucket
       93e28c7 README.md edited online with Bitbucket
       48619e0 README.md edited online with Bitbucket
       2e79a55 README.md edited online with Bitbucket
       3a17272 PacBio Assembly README.md generated by --help
       93e8124 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       69073ee Trinity + normalization settings for guillimin in rnaseq_denovo_assembly pipeline
       8d8571f BFXDEV-26 Added cluster_max_jobs param in config files + added warning in PBS scheduler if this limit is reached
       14194c6 Minor update in README-release.txt
       38f24d7 Version bump to 2.1.0-beta

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c661248 Run processing: Copy the rnaseqc metrics file to the bam directory to be compatible with nanuq
       e9ecd7a BFXDEV-165 Add missing Perl module for bcltoFastq
       44c06a0 BFXDEV-296 Changes to run processing documentation to use the generic help
       ecb2447 Run processing: Put rnaseqc results in a specific folder by librairy
       336a265 Run processing config: Only send email on abort and increase cluster wall-time to 48h for the fastq job
       3e8733a Run processing: Add symbolic link to STAR created bam to fix output file check on re-run. The original file is renammed to follow nanuq naming conventions
       7e8dddb Run processing: Fix rnaseqc sample list
       7e38acf BFXDEV-297 Add RNA-SeQC in run processing pipeline

2.0.0        Thu Dec 11 18:12:54 2014 -0500        669 commits

  Eric Audemard <audemard@ip03.m>      1 commits

       dc40fe2 update mammouth .ini

  Eric Audemard <eaudemard@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       62e2464 fix bug to find lib in pipeline python

  Eric Audemard <eaudemard@lg-1r17-n02.guillimin.clumeq.ca>      1 commits

       02b7ea6 fixe bug on puure before Ray execution

  Eric Audemard <eaudemard@lg-1r17-n03.guillimin.clumeq.ca>      2 commits

       3165f9a Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       bbe9d01 Bug fixed. Pipeline tested and validated on guillimin (step 1 to 21)

  Eric Audemard <eaudemard@lg-1r17-n04.guillimin.clumeq.ca>      2 commits

       cd0dfbc Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       ad3c17e bug fixed : 1) unmapped read in realign 2) small script error in puure

  eric audemard <eaudemar@imac6-ub.(none)>      15 commits

       3b2258d bug puure
       c0acd51 bug puure
       0a33731 bug puure
       c4ca954 update base.ini for guillimin
       b47c7c9 add base.ini for guillimin
       ab7535a 3nd script of PUURe done + some bug
       5e7c1fe Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a0b5418 2nd script of PUURe done + some bug
       9576d84 add file for puure pipeline
       844021a bug correction on 1st  script of PUURe
       c49f558 first script of PUURe done
       8979304 add art software (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) : Simulation Tools add perl lib path for SVDetect
       315bfc1 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_resources
       2b250c7 add script for install: SVMerge, RetroSeq, TIGRA, cmake
       f89034f add biopython

  eric.audemard@mail.mcgill.ca <eaudemar@abacus2.ferrier.genome.mcgill.ca>      4 commits

       caa90df add ini files. bug fixed on mkdir cov directory
       4771482 add ini files
       46f547e debuging step: fixe some bug
       00d5fb8 bug fixed on puure + add all output in dnaseq metrics (gatk, picard)

  Francois Lefebvre <lefebvr3@ip03.m>      1 commits

       e84abf3 picard version 1.125 module

  François Lefebvre <lefebvrf@gmail.com>      3 commits

       4ae7aff tophat and bowtie2 according to template install script
       35e4743 Minor changes to deploy script
       cc01370 no message

  Francois Lefebvre <lefebvrf@gmail.com>      54 commits

       22fd75c Updated R install scripts to reflect mugqic_pipelineS new name and deprecation of mugqic_resources
       ccb88f5 fix to deploy script
       3e061c3 update R packages list in resources
       22d0545 Fixed bug when -v is specified
       cd63425 more cleanup
       a3e5bbe more cleanup
       2be18df no message
       8536bf3 Fixed bug in prefix mecanism
       2e83482 Fixed default behaviour of update mode
       2bb2812 R_Bioconductor.sh: update mode is default
       7b04a5f First commit of the new R install scripts
       6e908ac STAR install script update to 2.4.0e, accounting for new folder structure in the archive (bin/,source/, etc)
       39922d9 shortStack install script
       b032d1d tcltk is part of base install
       c80bb97 Updated star install dev script
       46888c0 magittr and other packages added to install list
       d8773a3 Guillimin before mammouth
       b819fb6 cufflinks version change
       a62b519 Updated dependencies list (gqMicroarrays)
       0a5aeba Added sailfish module install script
       dfb486d bowtie 2.2.2
       97999d4 minor changes to top hat and vienna install script
       00d65aa Newer org.MeSH.* packages are too numerous and their installation take forever. Excluding them from the org.* packages list.
       7d51b6a Created gmap-gsnap install script (dev)
       af8b6c0 Added HTSFilter to del list
       21cc201 more package dependencies
       f34a34d Fixes to Trinotate related install scripts
       b937e6f Removed G phase1 from R deploy script
       8ae1b01 Multiple install scripts related to trinotate
       d3a8fe2 module install script for BEERS, RNA-seq Illlumina reads simulator
       7aca566 no message
       ba527a2 GCC module call for phase1 only: before compilation and called by R module for run time.
       915b729 modified example R.sh calls; need sh -l for module call to work on guillimin phase 1 ....
       2402dee bash synthax fix
       34d5126 entrezdirect installation script
       2de4f30 Added -p option to R.sh: name of an env variable that specifies a prefix to -i and -m.
       ba1638a no message
       9903c2a cron modified, next clip module
       c071c1e no message
       6c8c26f Nextclip install script
       99a9084 wrong file name in cron script
       c4ed552 no message
       d21308e temporary commits.. sourcing from Dropbox because problems with abacus
       42a0781 Cleaning and and chomping module files too
       21ef6d8 Changed chmod commands
       0e7cacb Problem with roxygen2 on CRAN. install.packages does not find it for R<3.0.2. Msg posted to r-help. R.sh will not roxygenize until this is fixed.
       67164c5 Last commit to R install scripts. Three files added: - R.sh is the workhorse which checks installs R/Bioc if necessary and performs updates. List of dependencies is hard-coded within this script. By default it will install/update the latest R version to $MUGQIC_INSTALL_HOME_DEV, hacking build.R to insure any subsequent package installation is performed with umask 0002. It can also install specifyc versions with -v, install to a specific location with -i and create the module file at a specific location with -m. The -f option forces rm -rf on any previous installation. The latest R version number is obtained by parsing the VERISON file from a R-latest.tar.gz download.
       3c6a6c6 tar is silent, chmod done outside R
       baacdad Combined installation/update in one script. package list and old install script will be removed in later commits
       d186f78 Added R/Bioc Update script
       453fdbd added varscan install script
       8a2f543 Updated R install script
       a3a16f3 list of R packages deps now linking to mugqic_ressources
       22b0242 no message

  Joël Fillon <joel.fillon@mcgill.ca>      363 commits

       bd175de Updated README-release.txt with more info
       21cbf65 Version bump to 2.0.0
       c3aea54 Added utils/ to PATH in module mugqic_pipelines
       ae6ae81 Updated module mugqic_R_packages to version 1.0.0
       b91ac7b Standardized mugqic_pipelines module install script
       493da6e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       6eb4df2 Updated PacBio pipeline steps as docstrings
       2ebef9d Updated module mugqic_tools to version 2.0.1
       718961a Fixed bug metrics.matrix sort tmpMatrix.txt before join
       ab9e001 Fixed bug Picard sort_sam: add BAM index as output file only if sort_order='coordinate'
       6610814 Added trimmomatic adapters FASTA files in bfx
       975beb6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01f101c Minor update in rnaseq de novo pipeline
       4b11101 README.md edited online with Bitbucket
       d5b5827 README.md edited online with Bitbucket
       e897da5 Updated SnpEff to 3.6 and snpeff_genome=GRCh37.75 by default in dnaseq pipeline
       078f739 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d0e8ec0 Separate bed2interval_list job from picard.calculate_hs_metrics in dnaseq pipeline
       347f888 Updated mugqic/tools/... to mugqic/mugqic_tools/... in config files
       5b93fc3 Updated and standardized module install mugqic_tools-2.0.0
       b7626d2 README.md edited online with Bitbucket
       fd69eeb Renamed mugqic_pipeline into mugqic_pipelines
       36aed54 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       66a477b Added multiple input files for star step in rnaseq
       35ac8b0 Added star memory/cpu settings in rnaseq.batch.ini
       9990a91 Minor output changes in star.py
       dc44a1d Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f7667c9 Added step multiple candidate input file support + updated trimmomatic and dnaseq bwa_mem_sort_sam, picard_merge_sam_files accordingly
       a6dda26 Added BAM index file as output of Picard merge_sam_files and sort_sam
       8c3c87f Removed unused trimmomatic skip option
       bc4bf73 Adjusted picard_merge_sam_files ppn value in dnaseq.guillimin.ini
       11c164c Added symlink to BAM index (*.bai) in nanuq2mugqic_pipelines.py
       ae0a63f Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6cca375 Fixed BED Files split bug in nanuq2mugqic_pipelines.py
       5edef09 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       98f2dd85 Updated dnaseq dbNSFP2.0.txt to dbNSFP2.4.txt.gz
       a5926c5 Replaced dnaseq dbsnp and known_sites parameters by known_variants + updated genome config .ini files with dbsnp_version + minor fixes
       11f6245 Added vcf.gz.tbi in genome install
       4cf7e43 Added log_report.pl and nanuq2mugqic_pipelines.py
       86167d7 Reverted to previous pacbio_assembly name
       fd9fd0f README.md edited online with Bitbucket
       3e0de92 Added RNA-Seq report step (needs update for cuffdiff section)
       165073d Removed differential_expression.goseq for cuffdiff and added self.gq_seq_utils_exploratory_analysis_rnaseq jobs in RNA-Seq pipeline
       3ea16d3 Merged in jfillon/readmemd-edited-online-with-bitbucket-1417115922753 (pull request #6)
       7dbb188 README.md edited online with Bitbucket
       8c70917 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a3f02c3 Recovered previous gq_seq_utils_exploratory_rnaseq modifs
       11863cb Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       62964d7 Update rnaseq_denovo_assembly.guillimin.ini with generic metaq, proc, pmem cluster settings
       29a8087 Upgraded smrtanalysis module 2.3.0 with patch 1
       08552cb Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f604b0f Recovered rnaseq.base.ini with new differential_expression and goseq settings
       6fb2716 Updated all R modules with mugqic/R_Bioconductor/3.1.2_3.0 and mugqic/mugqic_R_packages/0.1 if necessary
       4bca911 Fixed merge conflict
       8aab691 Added exploratory rnaseq step (partial)
       fe1bcae Added genes file in rnaseq config
       1101652 Updated Star module to 2.4.0f1
       f64a60c README.md edited online with Bitbucket
       9c0818a README.md edited online with Bitbucket
       c0fc855 README.md edited online with Bitbucket
       45c78c0 Merged in jfillon/readmemd-edited-online-with-bitbucket-1416841132369 (pull request #5)
       a7dfa97 README.md edited online with Bitbucket
       6583d2f Completed implementation of cleaning feature
       cb4c0c6 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       54378d5 More on pipeline cleaning feature
       8113c51 Fixed merge conflicts + first implementation of pipeline cleaning feature
       8c3c120 First implementation of exploratory step in rnaseq pipeline
       4ccebd1 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       9c33f31 Added Gene Ontology files in genomes + added goseq in rnaseq pipeline
       5b2a7a9 First attempt to install genome Gene Ontologies annotations
       1a34018 Fixed bug forgot super() call in Illumina and PacBioAssembly pipelines
       a46d772 Updated all config base.ini with module_mugqic_tools=mugqic/tools/1.10.5
       97243c7 Tagged mugqic_tools 1.10.5
       10f05ae README.md edited online with Bitbucket
       d624f80 README.md edited online with Bitbucket
       b52b592 Added design file description + minor changes in README.md
       c79b9c9 Updated rnaseq_denovo_assembly differential expression cluster settings
       e8d968c Updated rnammer cluster resources
       b882786 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e703a73 Added Python interpreter version check
       9fb0472 README.md edited online with Bitbucket
       fe5d954 Updated ini files with RAP_ID comment and rnammer settings
       6fd646d Removed openmpi settings on abacus
       1486f61 Design parsing support for empty or '0' fields
       171bf2c Removed command dump in .done file due to bug
       f4226eb Renamed readset and design file 'SampleID' column into 'Sample'
       74bef6d Updated cluster queue settings for rnammer_transcriptome in rnaseq_denovo_transcriptome
       a8ba35e Fixed version _argparser bug
       f9796a5 Dumped command in .done file + moved version in MUGQICPipeline
       a3d02bb Added differential_expression in RnaSeq pipeline
       3873197 Fixed syntax errors in rnaseq
       b0494b9 Added readset path normalization and variable expansion in readset file parsing
       531b02b Fixed PacBio readset parsing bug
       f17853b Added guillimin config file for rnaseq_denovo_assembly pipeline
       bf1e01f Moved genome_configs into resources/genomes/config; updated README.md accordingly.
       7eb3ada Imported mugqic_resources repository as 'resources' subtree
       be12079 In PacBio summarizePolishing, removed redundant cmph5 sort and used a copy of cmph5 file to prevent Job from being always out of date
       a3cf844 Updated PacBio pipeline with smrtanalysis 2.3.0
       76f9f33 Improved Job debug log
       d658c84 Minor cosmetics changes in README.md and PacBio readset file
       7ce1851 Added Bitbucket URL in command line help
       aab8a99 README.md edited online with Bitbucket
       0c66e6f Added pipelines/pacbio_assembly/README.md; removed pipelines/pacbio_assembly/README.txt
       49b0f0c Added
       48eb836 README.md edited online with Bitbucket
       e0ad931 Added pipelines/rnaseq/README.md
       329a627 README.md edited online with Bitbucket
       744fe94 Added pipelines/dnaseq/README.md
       e76ca73 README.md edited online with Bitbucket
       3426e21 README.md edited online with Bitbucket
       0e4226f README.md edited online with Bitbucket
       24a3bf3 README.md edited online with Bitbucket
       415922f README.md edited online with Bitbucket
       8bd144e README.md edited online with Bitbucket
       34db02a README.md edited online with Bitbucket
       44024ca README.md edited online with Bitbucket
       a8a0227 README.md edited online with Bitbucket
       0d57ec0 README.md edited online with Bitbucket
       f6bb2f4 README.md edited online with Bitbucket
       5442dec README.md edited online with Bitbucket
       d738c55 README.md edited online with Bitbucket
       73a1d8c README.md edited online with Bitbucket
       47ae71c More on README.md
       0ae6c2d README.md edited online with Bitbucket
       9a9751d More general documentation on pipelines
       cce8833 Replaced internal RAP ID with generic  variable in dnaseq.[guillimin|mammouth].ini
       49f5e04 Updated igv_genome config param to match generic <genome>.fa.fai index + updated all genome configs with Ensembl 77
       bd09fa8 Module versions update in config files
       ace8e6f Fixed pacbio filtering dependency bug
       c971597 Minor comment update in homer module install
       0ba8de7 Changed DnaSeq snpeff_genome to hg19
       8abed9c Fixed bug readset.fastq[12] attributes should be writable
       dd78287 Reverted SMRT Analysis module install script to version 2.3.0
       9a34f43 Added SMRT Analysis 2.2.0.133377-patch-3 in module install script
       68c77e5 Update SMRT Analysis module instal with 2.3.0 patch 0
       b06999f Fixed dbSNP download_path bug
       4612110 Updated dbSNP to build 142 for Homo sapiens genomes GRCh37 and GRCh38
       b8e26ad In genomes/install_genome.sh, added functions skipping if up to date + added rRNA FASTA creation and BWA indexing if present in ncRNA fasta
       49c8cf3 Moved star_index/ into genome/
       e0fe746 Added install_all_genomes.sh ; fix STAR module version for index creation
       7e92a39 Added STAR index
       d67a0c4 Added version number + updated README.md (incomplete)
       515ab99 Added star module install script
       39d4b84 Updated gqSeqUtils report calls and parameters
       162ec0b Fixed merged conflicts
       93a94e9 Updated RNAseq nozzle reports with Python version
       e4bae26 Updated DNAseq nozzle reports with Python version
       5ad7b16 Renamed bio module into bfx
       ad6c8b1 Added config file list support in gq_seq_utils report  and pacbio_assembly
       83cd491 Added report and compile steps in pacBioAssembly pipeline
       59fafec Added log.debug for job.is_up2date + set config.filepath with config.trace.ini file
       c90118a Added raise NotImplementedError for PUUre and RRNATagger pipelines
       4702215 Updated config genome paths with new genome files organization
       13471eb Fixed bug smrtanalysis.reference_uploader transforms '-' into '_' in fasta filename
       1d948bc Moved /sb/programs/analyste/genomes/blast_db/README.txt on bitbucket genomes/blast.sh + move oryCun2.sh in genomes/old/
       46062c5 Updated default dnaseq genome paths + readset library, run, lane and config sequencing_center are optional for dnaseq bwa mem
       45ec207 Fixed permissions bug for Homo sapiens genomes install
       eb7ea9c Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d052f61 RELEASE -> VERSION; DBSNP from NCBI for Homo sapiens; symlink DBSNP to Ensembl VCF for other species
       f103474 Added pacbio mummer step
       f2e84d7 Added pacbio steps up to blast
       29af2de Added job input/output files debug in pipeline
       3ae18d9 Fix cluster settings from job name instead of step name
       1da5b28 Fixed pacbio pbutgcns cmd missing &&
       412287e Check job input/output files as valid paths instead of regular files
       a072760 Update genome installs with Ensembl release 77
       a38bdfd Fixed rnaseq cluster_submit_cmd_suffix bug
       7e0c6e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       799ccb4 Fixed permissions on genome files
       7b1b697 Use Ensembl genome instead of 1000genome
       c7dd349 RNASeqQC now uses GTF with transcript_id only
       7aa6712 Removed abacus config .ini
       b6689cc Updated ini files with proper genome keys
       2e9c839 Added genome configs
       802c39d More on pacBio assembly steps
       49bdc5c Updated rnaseq genome config names
       e91a300 Increased abacus cores and memory for insilico_read_normalization_readsets
       59e7acc Added DBD::SQLite in perl install dependencies
       3fae17c Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       ca5b125 Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       f551b41 Fixed trimmomatic headcrop length order + insilico_read_normalization_all job name
       eeb1ff6 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       37a1260 More pacbio steps + here-document for bash command to avoid escape characters in output
       4cdd7a9 Added DBI and PDF::API2 as CPAN dependencies in Perl install
       0c4d6cd Fix RSEM reference dependencies
       59f01f9 Standardized perl module install with reduced CPAN dependencies
       f51ae0d Minor trimmomatic thread update in rnaseq_denovo_assembly.base.ini
       c3be9ab Added rnaseq_denovo_assembly.mammouth.ini + minor path fix
       74e88c0 module_tools -> module_mugqic_tools + pacbio filtering step (incomplete)
       1870e1a Added rat,dog,chicken genomes + python fix compiled in ucs4 + matplotlib manual install
       3599716 Standardized chipSeq modules as separate weblogo, homer, macs install scripts
       07f6d3b Tagged mugqic_tools 1.10.2
       f29ef00 Various fixes in rnaseq + pacbio first draft
       ea16a1d Removed picard_index and sam_index subdirectories + error tweak for gunzip human_g1k_v37.fasta.gz
       2cb5f0b Added dnaseq alignment readset subdirectory + mugqic_log not sent if pipeline has no jobs
       94590ba Cleaned genome_configs
       17a3b31 Cleaned genome_configs
       ec62c58 Minor documentation change
       a62ab92 Removed old RRNATagger files
       3d52780 Added MUGQIC parent pipeline with call_home log function
       a409b72 Added config trace creation from merged input config files
       178d1c1 Change readset file column name Sample -> SampleID; remove -V from qsub; updated gtf tophat index path in rnaseq
       0b1bd44 Creation of the major model genomes with new standard organisation
       8a6a093 Fixed puure.abacus.ini filename
       dba20ba Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       cf3c3ba Fixed merge_and_call_gvcf output file name + updated dnaseq with mugqic/GenomeAnalysisTK/3.2-2
       8614cd1 Updated rnaseq config files and parameters + minor dnaseq config update
       ab1c6d1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a1466e7 GENOME_INSTALL_TEMPLATE with Ensembl support
       0b2880b Fixed mkdir metrics in illumina.py
       66b529d Minor fix in dnaseq.batch.ini
       e002351 Minor fix trinotate nb of columns + scheduler date formatting
       bff0301 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       bc2b9fc RNASeq De Novo pipeline implemented except nozzle report
       23a4a72 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e621271 Fixed samtools_sort and baf_plot options for dnaseq on mammouth
       4491132 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       7b211be Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       c91f835 Added fastasplit in RNASeq De Novo pipeline
       dc50442 Minor fix in dnaseq.mammouth.ini tmp_dir and R version
       a6570ef Minor scheduler formatting
       2a9c298 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       edaf6f1 Added normalization steps in RNA-Seq De Novo pipeline + fixed python lib bug with absolute path
       f88a448 Modified -j and -c options descriptions + removed all .pyc files
       47a1c95 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95b82a3 Restandardized mugqic_tools install template
       2ab7d2a Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       10f5ecc Added first draft of Puure pipeline
       82006ac Added dnaseq mammouth config file
       e4b99ee Added dnaseq guillimin config file
       6f1cf24 More on GENOME_INSTALL_TEMPLATE.sh
       5dc0690 First draft of GENOME_INSTALL_TEMPLATE.sh
       4eb9823 Minor fix in Python module install
       c950b25 Fixed python module install + minor fix in MODULE_INSTALL_TEMPLATE.sh
       ddd1e5d Standardized python module install
       8eede7f Moved AMOS module install in dev and standardized it
       e2912e0 More config standardization
       e22b5c4 Minor fix in ucsc module install
       0c9abe7 Standardized UCSC module installation + minor template improvements
       2873b94 Fixed merge conflicts in bio/rrna_amplicons.py
       475e92c Reorganized all config parameters
       b8d0673 Minor fixes
       d609431 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       d73fe49 Fixed minor tophat bugs + scheduler proper exit code
       e19aaae Added rnaseq steps up to cufflinks
       6cf0a06 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b65411f Added cufflinks module
       104113f Added gtf in htseq + argparser factorisation in pipeline classes
       8b68b09 Added bio/htseq.py
       dd38967 Fixed newline and output file bugs in R mugqicPipelineReport + batch scheduler with job.done
       96fcc8e Fixed merge conflict in samtools.py
       631d37f First htseq draft in rnaseq + snv fix in dnaseq
       655c6bb Added multiple ini files argument feature + first draft rnaseq raw_counts
       7f92b49 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a351008 Added rnaseq wiggle step
       3d3d844 Fixed input/output files bug in concat_jobs()
       bdf6c8b Removed all .pyc from git
       f6995a9 Added rnaseqc step in rnaseq + other minor fixes
       cec67b7 Grouped common dnaseq and rnaseq features in a parent pipeline Illumina
       192486b Fixed minor bugs
       4441fb9 First rnaseq draft; dnaseq with os.path.join
       443c9c9 Update input/output files for metrics_snv_graph and deliverable
       f776268 Modified sys.path according to dir changes
       710c0bd Reorganised directories
       55e7904 All dnaseq steps implemented (need to test though)
       7cde825 Added dnaseq steps: rawmpileup, rawmpileup_cat, snp_and_indel_bcf, merge_filter_bcf, filter_nstretches
       af2c35a Added dnaseq steps merge_and_call_gvcf, dna_sample_metrics
       b0ea4b8 Added dnaseq haplotype_caller step
       952ebec Formatted batch scheduler + calculate_hs_metrics, callable_loci, extract_common_snp_freq, baf_plot steps
       2a3b223 More dnaseq steps + job input files validation
       dba5767 Moved python scripts in mugqic_pipeline directory
       c941e47 Removed step loop + standardized step names
       96e6e96 Added indel_realigner step + fixed job is_up2date with dependencies
       6cedd33 Added first GATK functions + extended dnaseq pipeline
       2c2937c Added first gatk draft + parse Nanuq readset file
       397fb36 Added global config object + picard module
       4dd79a1 Added group, pipe, concat job functions
       fef73cd Added DnaSeq trim and mem + force option
       24057c2 Torque scheduler test OK
       b484d11 Added torque submit header
       dc45f44 Added scheduler and log
       b6fada2 Updated mugqic_tools install with version 1.10
       16e9aa3 Added argument parsing and more
       2e6cc0e First Python version with basic args parsing
       3ee3e3e Core classes coded
       97d8940 More on python prototype
       fa8dd9a Added DB_File Perl module + minor permissions fix
       a24fd3a Added prerequisite for modules rnammer and trinotate
       e6d3ded Standardized trinotate and dependencies module install
       af7cea5 Standardized BWA module install
       53b9605 Moved R&D modules in dev/
       fbccba1 Added archive check in MODULE_INSTALL_TEMPLATE.sh to avoid redundant archive download
       1e72c28 Minor fix on template install: archive dir with /
       8898ddf Restored blat module in dev
       54192b2 Removed blat module (already part of ucsc module)
       6a70836 Prefix java module by openjdk
       d5575c0 Updated Java module install with OpenJDK
       3e9546b Updated RSEM and Trinity module install with latest version
       f80d0d2 First Python draft
       65176a9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       28f5ab6 Fixed conflicts in merge master into bam2fastq branch
       26c226a Resolved conflicts with master branch
       1729a4a More on object-oriented dnaSeq pipeline
       2378320 More on object-oriented design
       a3dabbe More on object-oriented redesign
       3d942d8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       455ebed More on object-oriented pipelines
       20bc7af Minor variable name changes
       b2c5d0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       d017c57 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       9b1ae90 Incomplete version of bam2fastq dnaseq pipeline
       e87e5b9 Added Parse::Range in Perl modules
       29bd3d9 Solved merging conflicts in rnaSeqDeNovo.mammouth.ini
       a91dbdc First object-oriented version of RNA-Seq De Novo Pipeline
       2255e9c Minor RSEM install script fix
       a7c1275 Updated RSEM install script by modifyin Perl script shebangs with /usr/bin/env/perl
       a37edd0 Updated mugqic_tools to version 1.8
       063caa2 Added PerlIO::gzip as module dependency
       2fa1750 Added dev/prod/ environment in MODULE_INSTALL_TEMPLATE.sh + modified all shebangs with bash
       ae10289 Merged dnaSeq.pl conflict from master to bamToFastq
       a6d0d76 More on bamToFastq process (uncomplete)
       f911e8a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       89a3066 Merged master in local branch and resolved conflict
       ced737e First bamToFastq attempt in dnaSeq.pl
       338afbe Minor fix in snpEff archive name
       8154c76 Standardized snpEff, VCFtools install scripts; added Java install script; added various permissions
       17a4f95 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       c351f09 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       b64d2e7 Merged with master
       27c29de Even more on bam2fastq
       9c38e40 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       04f66e4 More on bam2fastq
       ca1172c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       e37f19a Some steps toward bam2fastq support
       ab015b8 Standardized tophat and cufflinks module install script
       201cc7d Standardized bowtie, exonerate, samtools module install scripts
       bbc946e Fixed module directory named "pipeline" instead of "mugqic_pipeline"
       25141fc Standardized mugqic_pipeline module install script
       732d095 Added others permissions for MUMmer and wgs-assembler
       493b28b Updated mugqic_tools version number to 1.6
       689141b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ea2d0ae Installed trinity and rsem in production
       9b4f23f Standardized Trimmomatic module install
       ff81b27 Updated RSEM version, RSEM and Trinity permissions
       e51836d Upgraded to mugqic_tools-1.5
       f205a7a Added BLASTDB variable in blast module install script
       b3f1686 Minor permission fix
       9147435 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e894200 Updated Trinity module install with last version
       819415a Updated mugqic_tools install script according to MODULE_INSTALL_TEMPLATE.sh
       5093baf Renamed module install template
       9473a90 Added read/execute permissions for others in module install template
       0e6cf37 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       99c14da Added ape + ctc R dependencies for Trinity differential expression analysis
       f6648d9 Updated Trinity edgeR PATH
       1f140b6 Removed comment in Trinity install script
       0348f39 Updated trinity and rsem module install scripts in dev
       5e0865f Moved some install script to dev/
       91784d0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       29666f4 Clean up obsolete/dev install scripts
       4ec387a Removed Makefile draft
       90d94d3 Minor change
       7711123 Initial import of genome and module install scripts from mugqic_pipeline repository to mugqic_resources repository

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       9f787f3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95ac9e3 BFXDEV-133 detected perl bug on homer and weblogo scripts using our own perl module, updated perl scripts shebangs
       4b5a245 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       72a19c1 BFXDEV-133 Adapt the software/module installation to guillimin phase2 requirements
       90db5b5 BFXDEV-82 version 1.7 of mugqic tools, added bug correction for R-tools related to chipSEQ pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      12 commits

       0418260 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d33f8e2 prepend popoolation path to PERL5LIBS
       c5a79a8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       05d6532 added Picard Mark duplicates no optical v 1.90 to the dev software repository
       5edf0ed pmarquis: Added genome setup for arabidopsis thaliana-TAIR10
       4bc9c17 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a8a9a0f added env-variable POPOOLATION_HOME to the module file
       13c9b20 added an installation script + module for popoolation
       f2046f1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       2c85b4f added plotrix to the list of R packages
       eaab7c1 creating R v.3.0.2 module in dev
       3a6476d installing R on the test environment

  jtrembla <jtremblay514@gmail.com>      24 commits

       4d370a6 rrnatagger.py continuation
       a228ba6 completed bio/rrna_amplicons.py
       a4bbfc7 continued coding of rrna_amplicons.py.
       f82ac52 writing rrnatagger pipeline in python. First step works.
       cf40514 Added two functions (that actually works :-)). I need to convert the rest to python syntax.
       bc82663 added rrna_amplicons
       d045a6a renamed RRNATagger to rrnatagger
       d0a7f34 Added rrnatagger.py and list of steps.
       9f616ee added lib path.
       b0d9af6 Added bamtools.
       907ef6b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       87e713a Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e3ddbbb minor fixes.
       058fa82 Updated bwa version.
       6432b27 New module install scripts.
       cecd97f Added path (prepend-env) to root/R-tools.
       c61c904 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       119c7ad minor modif to module list for perl.
       2fe2654 Added comments to smrtanalysis installation (2.2.0). Tested and installed.
       9858f66 Added muscle aligner.
       9403ab7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources corrected typo, then merges..
       ac2a43e corrected typo.
       5c1e23a Added chmod at the end of the install module "scripts".
       79e2d4a Added list/script to install modules. Could be improved...

  Julien Tremblay <jtrembla@ip03.m>      1 commits

       72347df added module DB_File to perl module installation.

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       7d2d5a3 Added RRNATagger-tools to prepend path in the module file. BFXDEV-31
       b792e6a added python install script.
       1be6f84 Fixed two missing paths for PERL5LIB
       0b26b2d Added perl to mugqic modules.
       82b9d19 Modified install scripts according to our group's template.
       3156d84 Added gnuplot-4.6.4 and memtime-1.3

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      1 commits

       923d661 Removed genelenght.r. Use mugqic_tools/R-tools or gqSeqUtils directly with HERE docs or Rscript -e

  lefebvrf <lefebvrf@gmail.com>      14 commits

       0be2b76 Removed superfluous (already in base) dependencies from list
       4cf1272 violin plots package added to list of dependencees
       16882a5 PerlIO::gzip is a dependency for Trinity clusterProfiler as R package
       0530c57 Guillimin phase 2 added to daily R deploy
       9a37c3c R.sh needed a module call to a more recent gcc than system one.
       e90b9ca devtools::install_local() would install packages to home folder if exists...   .libPaths(.Library) solves the problem
       923b50a —vanilla removed from R.sh too
       a13f582 blurb added related to last commit
       d028e16 Add creation of Rprofile.site to R install script. This will force using cairo X11 backend since cairo is not always set to default when available…
       321311b Added ViennaRNA and mirdeep2 install scripts
       33dfb87 Small hack to tools::build.R when installing R allows umask 002 (!)
       7496779 Corrected exit code, first version of R_deploy
       b8db61a Additional setdif(reps) to avoid duplicate installations -> slow
       2082c05 Fixed leading tabs problem in module file by using <<-EOF here tag. Neutralized R_LIBS to insure installation proceeds correctly when R module already loaded Added vanille biocLite() call

  lletourn <louis.letourneau@mail.mcgill.ca>      45 commits

       e76f595 Merged old perl changed into python
       0f7f0f8 Fixed haplotypeCaller output file name extension
       4098ccc Version bump mugqic pipeline to 1.4
       b9625ab Version bump to 1.5-beta
       7e8d28c Bumped bowtie to 1.1.1
       2528f92 BVATools version bump
       277abf2 FastQC version 0.11.2
       c1d9878 Version bumped Ray to 2.3.1, removed unneccessary fix now
       45b85da PIcard Version bump to 1.118
       9db0c54 Added pysam
       fc0ea4c BWA version bump AND fixed the script, it was a bad merged script
       cc5aa56 Updated GATK and mutect
       2cd15d6 Merged diffs
       b10f91e Added gsalib and added configure switch
       cb1777f Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       3677916 Version bump pipeline to 1.3
       58fff5f Force output file name
       ef35a74 Version bump bvatools to 1.3
       a441cda Version bump of mugqic_tools and snpeff
       9b6fa74 Added ascat
       a669439 BWA Version bump to 0.7.8
       6109739 Added the jellyfish tool
       6108922 Updated gatk
       09e5270 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       44e530a Added Mutect home
       9a96cbf Version bump
       6827d55 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       adda765 Version bumped wgs-assembler, bwa, picard
       674305d Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       9e6a4c6 Version bump bowtie2 bwa picard
       287aa23 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c2aac9a Added aspera instructions
       701df8c Version bumped vcftools to 0.1.11
       6f5a6ed Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c8f949a Version bumped blast to 2.2.29
       d0dab7f Bumped version of picard
       4e6fcc4 Updated snpeff
       ff4f397 Version Bump of the pipeline
       34c15b0 Added variable to access the repo's location
       fa87e19 new Ray version
       a6402ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ac6513f Updated BVATools to 1.1
       ea6aff9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       dca240f Added BVATools
       71bfcf2 Added wgs assemble

  Louis Letourneau <louis.letourneau@mail.mcgill.ca>      1 commits

       e1c5341 BFXDEV-246 Version bump

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      47 commits

       0eff10b Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6272cb0 Python - RNAseq : allows more flexibity in fdr and p-value for goseq AND remove native goseq approach - BFXDEV-289
       e78717a Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8694d0c Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       218dbdc Python - RNAseq: update rnaseq.py
       a133310 Python RNAseq - updates
       6a061db Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       5e8df54 Python RNAseq - updates
       61a12cf Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b1f9bc8 change RNAseq python update
       5c9e82f change RNAseq python cuffnorm output to include sample names
       015fab6 remove RNAseq python step bug exploratory v2
       582e5e3 remove RNAseq python step bug exploratory
       b14fbf9 remove RNAseq python code conflict
       62a0e4f PYTHON -RNASEQ: star, cufflinks, htseq updates
       40970ff Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e600a91 Python: RNA - update star + DGE anbd goseq input oputpuyt file correction
       01fccf9 resolved pipelines/rnaseq/rnaseq.py confilcts
       87f92a3 major RNA implementation: STAR, picard, cufflinks, etc...
       a39a35e Python RNAseq: add cufflinks > 2.2.0 new workflow - cullinks - cuffmerge step done
       d902b4d Python RNAseq : update STAR - add optional read sorting during alignment
       8f572bc Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8305fe7 Python - RNAseq: add STAR alignment 2 pass && add utils folder/function for methods generic unrelated to any software or tools
       c713736 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       4cd30e3 update samtoiols module install to 1.0
       860eba6 python : RNAseq - change ini
       a6ea523 On goinig adding star to RNAseq
       33f80d6 add dnacopy in R.sh edited online with Bitbucket
       3bdf7cb Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       0f1910e correct python.sh script that was not workin for the modules: BFXDEV-46
       27348a0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       793e633 add specific RNAseq packages to R package list
       1016ccd replace correct permissions in mugqic_tools - BFXDEV-55
       fda2802 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       40d5ef5 ngoing python module package
       03dab09 update breakdancer module
       7f60581 update breakdancer module
       494455b pdate breakdancer module
       9deedd9 change in samtools modules and add breakdancer module
       02fde14 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       bfdd544 add module igv
       60e6732 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       6615a3b create module install sh for the mugqic pipeline
       30983a9 update nugqic_tools to the new version 1.2
       834b58e Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       19f7508 update mugqic_tools
       5e75b3d update modules/mugqic_tools to point to the new repository mugqic_tools with tag v1.0

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       14060c1 Remove conflicts in modules/mugqic_tools.sh and modules/dev/star_dev.sh
       cb81397 update modules/mugqic_tools.sh to 1.10.4
       b2db257 remove decrepated python module script and add a new one
       16d3068 replace dev module in module/dev/ cufflinks_dev.sh  star_dev.sh
       fd6119e Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       85a4523 add genomes/oryCun2.sh  modules/cufflinks_dev.sh  modules/star_dev.sh
       b35f6de up-date mugqic_tools.sh version to 1.7

  mmichaud <marc.michaud@mail.mcgill.ca>      54 commits

       4617cd8 Run processing: Remove duplicate '.' in dup bam file
       f72f529 Run processing: Run blast on bam file, fallback on fastqs
       ae2fcbd Run processing: Add the suffix '.sorted' in the readset bam name
       7bb7a50 Run processing: Merge md5 jobs (fastq+bam) and run qc graphs on BAM (or fallback to fastq)
       5603e88 Run processing: Run DepthOfCoverage even if there is no BED file
       a81e719 Run processing: Use bcl2fastq 'force' option instead of removing destination folder
       1890439 Run processing: Remove dependencies dependencies from copy job
       6f33eaa Run Processing: Use a kind of factory to manage the different aligners
       e006fa9 Run processing: Remove the dependency to samtools by using picard to generate the STAR BAM index.
       f46022c Run processing: Move md5 step later, to optimize job allocation and to minimize condition when the md5 job is finished before the copy job is submitted
       8d06973 Run processing: Remove unused adaptor columns parsing
       fd4705a Run processing: More documentation (sample sheets)
       60abc0c Run processing: Output the index metrics file in the output folder, not the run folder
       88a0907 Run processing: Use xfer queue for downloading BED files
       ba1b4d6 Run processing: Add samtools module version and increase number of core for STAR (a test job took 71G of ram)
       3186121 Run processing: Manually generate BAM index file when using STAR.
       1b42a02 Run processing: Fix usage in readme
       f4e87a7 Run processing: Add documentation
       a6e1838 Run processing: Add missing '.' in bam name
       7aad80b Run processing: Fix copy step exclusion
       7100d7c Run processing: Copy output file is now in the copy destination folder
       073b8a9 Run processing: Change back to manual copy job depedencies gathering
       ff84ee0 Run processing: Initial RNA-seq support with STAR
       b5969a3 Run processing: Return empty list when there are no input for the copy job
       b5b30a6 Run processing: Fix copy job inputs
       58e4326 Run processing: Fix end copy notification job queue
       43c153a Run processing: Fix qc output file
       9d9e87a Run Processing: Fix race condition in HsMetrics interval file creation
       7dc2b9b Run processing: Fix configuration for blast and qc
       edb3269 Run processing: Generate sample sheet in the submit job method
       f5f23ff Run processing: Fix config for index step, now using 'index' category
       348972b Run processing: BAM metrics are run in markdup output
       176e96f Run processing: Various fixes for the first test run
       63944b9 Run processing: Various fixes for the first test run
       1043bb1 Run processing: Various fixes for the first test run
       c8ee9b0 Run processing: Add basic support for different aligners
       aa16b2d Run processing: Get copy step dependencies by introspection
       c74a06b Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       3af310b Run processing: Delete existing 'Unaligned' folder when running fastqs
       583c80f Put Illumina configure bcl2fastq in a job
       aa93387 Run processing: Use . as job name separator
       ddafcfe Use $MUGQIC_INSTALL_HOME variable in config file to replace hardcoded paths
       79d09cc Fix arguments for IlluminaRunProcessing by removing readset argument from MUGQICPipeline
       cbe8ddc Run processing: Use run_dir instead of run_directory to follow output_dir convention.
       aa213fe Run processing: Improvement to the sample sheet generator to loop only one time
       7d16bcf Run processing: Seperate config files for MiSeq and HiSeq. The base config file still hold for a HiSeq
       df783bc Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       2b53d11 Run processing: Add wget download of sample sheets and bed files. Fix Casava sample sheet generation.
       07aacb6 Add some imports from __future__ to ease the transition to python 3
       7af4b28 Add copy job dependencies
       82841bf BFXDEV-283 Fix reference for coverage calculation, add copy step.
       213e716 Change illumina_run_processing name to follow code conventions
       8d8e1af Code convention changes
       f51f5e0 First draft of the illumina run processing pipeline in python.

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.ferrier.genome.mcgill.ca>      1 commits

       acde237 adding sh tools for riboRNA

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      1 commits

       f6ddfbe BFXDEV-112 install genomes from IGENOMES

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       2fc3f9d usearch.sh

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      2 commits

       481bc98 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       7ff159c gene_length

1.4        Mon Nov 17 13:15:48 2014 -0500        139 commits

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       667c6f7 README.md edited online with Bitbucket
       ee9efd5 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvr3@ip03.m>      6 commits

       3355905 Previous minor bug fix (dots in sample names) actually introduced major bug.
       0bc3e58 temp files remove
       93c2905 rnaseqc other options possible
       c045541 rnaseq .ini file more tophat options. BFXDEV-215
       8270318 Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.
       0ef6afc Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.

  Francois Lefebvre <lefebvrf@gmail.com>      7 commits

       242f2bc Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c004e51 Passing projectName to SubmitTocluster can create invalid job names. Replaced with string
       e578fb4 variable name was inappropriate
       eab0310 Updated rnaSeq mammouth template .ini file.
       202215e Defined a job name prefix for wig zip.  'metrics'  as a job name was not enough information
       59f1528 cuffRescolumns and dgeRescolumns now in goseq param section. Also adjusted those values for UCSC genomes hg19 and mm10  templates, original value did not work.
       524dd82 overwrite.sheets=TRUE to avoid common problem of updated projects

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      4 commits

       257c44d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4b8a23b patch3 applied; see BFXDEV-260
       3c6e850 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ff5fe3d change R module due to crashing nozzle report generation, see BFXDEV-255

  Joël Fillon <joel.fillon@mcgill.ca>      34 commits

       0860442 Updated pacBio .ini config file with new smrtanalysis module name: 2.2.0.133377-patch-3
       c2d5000 Updated config genome paths according to new genome organization
       030e601 Updated chipSeq mammouth config with new Homer 4.7
       be23877 Added explicit Python module in rnaseq cuffcompare step
       9d16e5a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c256913 Added newline after mugqicLog command
       8f3c5f3 Updated rnaSeq.mammouth_hg19.ini with generic /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod path
       1fd8a63 Removed explicit RAP ID in RNASeq De Novo config files
       0efe22d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6aac202 Removed moduleVersion.htseq ini ini config files
       ded1a46 Fixed missing java module in igvtools
       2c37b3e Minor fix in dnaSeq step_range help
       cf08611 Added default RAP ID in RNA-Seq De Novo guillimin config file
       e67468a Added --exclude <fastq> option in illuminaRunProcessing rsync for samples having BAM files
       98420e4 Uncommented phone home in dnaseq
       0e6e8a8 BFXDEV-203 Create one .done with checksum instead of one per output file + update config files with default adapters path + update perlpods removing -e option
       9886061 Fixed bug missing Library Source -> column not mandatory + updated pipelines/rnaseq/rnaSeq.guillimin.ini with accurate module versions
       e07ffd8 README.md edited online with Bitbucket
       d02ae67 Added comment to update Resource Allocation Project ID
       13bf5c9 BFXDEV-221 Migrated abacus/guillimin config files from msub to qsub
       d8a411e MUGQIC call home is now run inside bash script after job submissions instead of bash creation.
       5c2627a Updated chipseq mammouth config with python 2.7.6
       29163dd README.md edited online with Bitbucket
       bad3648 Added call home feature notice in README.md
       9413014 Fixed RNA-seq de novo resume jobs; added env variable WORK_DIR in pipeline; formatted bash output
       b2f3fd5 Another glob fix for prefix path check in LoadConfig::getParam
       65b9923 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       768184e Minor fix for config param prefix path
       74d4328 Added prefixpath check type in LoadConfig::getParam
       249a5cb Fixed transdecoder bug ln: creating symbolic link 'Trinity.fasta.transdecoder.pfam.dat.domtbl': File exists
       59ea4b6 BFXDEV-32 Fixed wrong transdecoder file path for missing PFAM 'cds.' prefix
       9c6fe6f BFXDEV-32 Fixed pfam missingcds. ID prefix + blast clusterCPU tag for guillimin and abacus
       73ce180 Added chipSeq pipeline cleaning
       9f6e203 Fixed rnammer missing modules hmmer 2.3.2 and trinity

  jtrembla <jtremblay514@gmail.com>      14 commits

       6c0720b --arg for low abundant clusters. after 99% clustering. BFXDEV-31
       d09db89 Added low abundance cluster cutoff argument. (After 99% ID clustering step). BFXDEV-31
       9166850 Put more lenient parameters for itags QC to make it more 'universal' for most projects, especially those for which quality of reads 2 is low. BFXDEV-31
       05e7b99 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f60fa07 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a996d3f Indentation correction. BFXDEV-31
       38d8da8 Fixed parsing for blastdb command. BFXDEV-30
       1b20686 Updated README for 454 data processing instructions. BFXDEV-31
       8bafb19 Added even more description. BFXDEV-31
       c003951 Added details to description of output. BFXDEV-31
       fcd21e5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9dcc9b Added step to filter number of blast results. BFXDEV-30
       5d0c803 Removed --vanilla. BFXDEV-30
       309ef5e replaced -num_target_seqs with -num_alignments. BFXDEV-30

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      1 commits

       a16e5e9 Added --sampleSheet argument to getMiSeqBarcodes.pl BFXDEV-31

  lefebvrf <lefebvrf@gmail.com>      2 commits

       b7f7d66 Changed TopHatBowtie.pm to put an end to  the .fa.fa.fasta.fa symlinks madness when installing genomes. Parameter is now the bowtie index basename, consistent with the tool's documentation.
       d29a172 Fixed dots in sample names bug BFXDEV-51

  lletourn <louis.letourneau@mail.mcgill.ca>      31 commits

       3f785de Version bump to 1.4
       38d0a8f BFXDEV-39 Fixed realigner when only one job is given
       f728ac9 Updated bvatools
       5e168e6 Tweak guillimin parameters
       ea83ecd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1748e1a Updated parameters
       69b8f2e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5da052a Adjusted cluster requirements
       4cc7adf Removed useless param
       4dace3e BFXDEV-256 Added step range to paired variants
       cdbd4a0 BFXDEV-256 Added step range to paired variants
       b535570 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46d016c BFXDEV-254 Fixed GATK 3.2 change on CatVariants
       10fb56b BFXDEV-252 Removed flagstat
       125945f Fixed BAQ from pileup and ram from fixmate
       71366af Changed picard to version 1.118 to fix the freeze when an exception occurs in asyncWriter
       e9112eb BFXDEV-216 Removed per lane metrics
       0568788 BFXDEV-245 Fixed uninitialized error when no bams are present
       682a997 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6d2b309 Fixed CCDS location
       85578b8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       54fde94 Merged master
       407b9fe BFXDEV-236 optimized settings and split human builds
       419e0a0 Added missing perl module
       976ba1f Fixed missing HS metric, add 2 cores to bwa
       9cf6dc5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c3fab1c Missing validation silent
       5cfb61f Add genome versions of ini files
       de4a383 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       dd39c54 Fixed typo in mergeAndCallGVCF section
       3366970 Version bump to 1.4-beta

  Marc Michaud <marc.michaud@mail.mcgill.ca>      1 commits

       6b8f9f7 Add missing use

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      13 commits

       3bd1f63 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0e991b0 Update pairedVariant
       469415f RNAseq replace headcrop at the good position in the trimmmomatic command; remove by default headcrop from the ini files - BFXDEV-267
       534615d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       402db54 PairedfVariant.pl: change where pindel get the insert size info - BFXDEV-266
       87a9da6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       192600e PairedVariant: allow mutec to run without cosmic file - BFXDEV-263
       62e335d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9696dbc update dnaseq and paired variant ini files
       37ad08b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       21e3c79 remove conflicys in pipelines/dnaseq/pairedVariants.abacus.ini
       a156ae5 pairedVariant.pl formAt SV and CNV to new standard - part of BFXDEV-41
       cb5b91e cuffdiff now should not be relaunch in a restart if it exit correctly during the previous analysis BFXDEV-212

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       3c0fb6e update ini to fit the new mugqic_tools tag 1.10.4 - BFXDEV-275
       e5bd476 correct dnaseq bwa samnse dependency bug - BFXDEV-253
       6fdf802 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f3abd8c correct cuffdiff input double array issue when checking the job object is up to date

  mmichaud <marc.michaud@mail.mcgill.ca>      14 commits

       b9ef532 Fix genome path
       ec31045 BFXDEV-269 Change genome files hierarchy
       0234cd1 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       5eb75aa Fix usage to reflect new -s argument
       f16ff09 Allow more ram for DepthOfCoverage
       4033170 Fix CPU limit usage error by using less thread for the GC
       88e5096 More RAM for QC (to avoid java heap space errors) + More threads and more walltime for mem (to avoid wall-time exceeded errors)
       d946aa7 Run processing: Add option to force download of sample sheets
       7feb923 Fix quote escaping in filter
       8883201 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       0e60d3d BFXDEV-171 Trim index in the generated sample sheet according to the first/last index specified as parameter
       c5a35cb BFXDEV-210 Use Parse::Range for steps to run, as all other pipelines
       598ff88 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8e8ca8e BFXDEV-211 Don't send email when jobs are successfully completed

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      4 commits

       0cb571a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f934df3 MiSeq ini
       1f34db5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31ae450 add ini MiSeq

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      1 commits

       446b4d8 ppn=16 pour mem

  Pascale Marquis <pmarquis@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       90e7a61 rm illuminaRunProcessingMiSeq_PM.ini

1.3        Mon Jun 2 10:02:07 2014 -0400        109 commits

  Joël Fillon <joel.fillon@mcgill.ca>      17 commits

       9931258 Partial MUGQIC remote log for RRNATagger
       da59f18 Remote MUGQIC Log Report in chipSeq, dnaSeq, pacBioAssembly, rnaSeq, rnaSeqDeNovoAssembly
       521765d Updated RNA-Seq De Novo config files with latest module R/3.1.0
       f59efb9 Updated rnaseq_denovo default config ini files with trinotate steps
       0194584 More Trinotate steps in RNASeq De Novo
       9a5cd2d RNASeq De Novo config files conflict solved + openjdk
       1cfd288 Added file cleaning for RNA-Seq De Novo pipeline
       4382d80 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1092478 Solved conflicts when merging master
       42e7aca Updated RNASeq De Novo pipeline with Trinity version 20140413p1
       240d576 Beginning cleaning
       e1c9021 Cleaning of cleaning...
       ad6ed7a Fix on samToFastq/trimming dependencies in chipSeq pipeline
       02d3cbe Added BAM file check in samToFastq
       5e1aeb2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq_basic
       8ea8c98 Basic bamToFastq support for all pieplines
       9afcd9e Version bump

  jtrembla <jtremblay514@gmail.com>      32 commits

       6e3d8ed Added ini file for production purposes (QC assembly among others). BFXDEV-30
       462ecaa fixed end step for BB. BFXDEV-31
       06a4423 Added BigBrother modifs to RRNATagger pipelines. BFXDEV-31
       e0c59ac Added sample counting in pipeline loop. BFXDEV-31
       ebd465f Added a description of output files. BFXDEV-31
       e0272b5 Fixes to itags_QC. decision if primers are present or not. BFXDEV-31
       b199a3d mooooooore fixes. BFXDEV-31
       d67905f more fixes to ini. BFXDEV-31
       7b8a64c updated ini files. BFXDEV-31
       609c204 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       befd271 Changed arg to adapt from percent to X hard cutoff.
       16248cd Replaced percent cutoff by X cov cutoff. BFXDEV-30
       d79cc38 added / updated ini files.
       65ed0f3 ini files of RRNATagger changes. BFXDEV-31 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9e30f84 updated ini files for RRNATagger. BFXDEV-31
       37c1b0c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4cb24e1 ini file modif. BFXDEV-30
       de66a9b update celera specs for hgap3. BFXDEV-30
       0faaea6 Updated file check for restart at blasr ste. BFXDEV-30
       83c0cbf Modified way seeds.m4.filtered is handled. BFXDEV-30
       00482bd Updated ini file for pacbio assembly on guillimin. BFXDEV-30
       9174116 fixed lib path. BFXDEV-30
       f823b93 Fixed tab indentations. BFXDEV-30
       5f8123d Upgrade pipeline from HGAP2 to HGAP3. BFXDEV-30
       98d20c4 Removed unused SMRTpipe .xml files. Only keep filtering xml file. BFXDEV-30
       468aacb Updated mugqic tools module.BFXDEV-30
       a1ae253 Fixed output .mugqic.done file for pacbio dcmegablast. BFXDEV-30
       e7426a8 Added module load perl in subroutines. BFXDEV-31
       056cae4 added perl in ini files. BFXDEV-30
       cf512b2 pacbio ini file for abacus. BFXDEV-31
       85ebc8a Added blast parameters to ini files for guillimin. BFXDEV-31
       52a9f7f Added blast step for pacbio rRNA tags data type. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       cc1fab0 updated ini file abacus.
       9363d1f Updated ini file for abacus. BFXDEV-30
       c3a1189 replaced compareSequences with pbalign for numthreads. BFXDEV-30
       7189270 ini file for hgap3 on abacus. BFXDEV-30
       37df4a8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1bf18fd Fixes to HGAP3. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      38 commits

       8f38558 Version bump to 1.3
       c051f9a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       42735dd BFXDEV-207 Fixed when interval lists are created
       102ee1d Fixed gvcf bugs
       613a990 Added emtpy quotes to empty keys so they don't become ARRAYs
       b594056 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d39b5d9 Updated BVATools version to 1.3
       9eae448 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9cc1832 BFXDEV-204 removed varfilter
       82592ff Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8b4f2bd BFXDEV-196 Added onTarget metric
       2f57d4f Completed POD documentation
       8176ab0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1d6793d BFXDEV-161 Added last steps to Haplotyper caller
       240ca52 BFXDEV-198 Added simpleChrName in the ini since it was removed from the pm
       7f0d539 Fixed params
       39a7b38 BFXDEV-196 Added onTarget metric
       8e7a04e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5825461 Ram was too close to max
       db33ed9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89c5216 BFXDEV-193 Use new fixmate from bvatools. Fixed bad step dependency
       ff9d44c BFXDEV-192 Added possibility to have extra flags in indel realigner
       f30f24e BFXDEV-182 CollectMetrics sometimes needs to GC so 2 cores are needed
       20768a4 Merged changes
       de76e02 BFXDEV-181 Updated java vm version
       d946871 BFXDEV-153 Added a way to ignore readset status
       4ad3df1 BFXDEV-176 create MD5 on alignement
       52b96be BFXDEV-174 Test that the data is valid
       25f2a40 Merge branch 'master' into haplotypeCaller
       d4a588d BFXDEV-168 added depth ratio plots BFXDEV-161 added haplotype caller
       6596bf3 BFXDEV-166 fixed no index read, but one is associated
       7b4bdc2 BFXDEV-162 fixed hiseq recognition
       3ad3319 BFXDEV-157 Add callable region generation stats and BED in DNASeq
       9f0e145 Added BVATools depth of coverage as a transition phase
       17c7a8f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c775e18 Don't use BAQ since we use recalibration, this was decided awhile ago
       1c6b7a4 BFXDEV-156 fixed csv encoding issues
       46573b7 Updated settings for phase1 + phase2 merge of hardware

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       35060f6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7faafe0 Create pipeline cleaning lib; partially implemented with RNA cleaning sub only BFXDEV-74
       201124c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5177736 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       27c987a change how exit status is catched and exit the correct status in case of pipe discrepency BFXDEV-140

  mmichaud <marc.michaud@mail.mcgill.ca>      11 commits

       8c059d4 BFXDEV-155 Add rat and mouse alignment
       4e9b1e1 BFXDEV-185 Separate config file (HiSeq, MiSeq)
       637526f BFXDEV-164 Fetch sample sheets when they aren't specified and they aren't on disk
       49e0b8c BFXDEV-183 Download each bed file only once
       f87afda BFXDEV-184 Don't rsync Images folder. Was used on miSeq for debuging purpose
       ae9782f BFXDEV-186 Don't use the BAM generated by markdup
       c5d36af BFXDEV-180 Use processingSheetId as dependency id instead of sample name which isn't unique
       6f3d923 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       ec78976 BFXDEV-171 Add option to specify first and last index
       e3ad967 Add option to specify first and last index
       d8eacb4 Run Processing: Add dependency on the metrics in the copy step, in case the markdup & BAM MD5 was already done in a previous pipeline, but not other metrics

1.2        Fri Mar 28 16:11:44 2014 -0400        200 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       c18dfb6 RSEM more cpus on mammouth
       9fff58f rnaseqdeno mammouth ini tweaks for trimming and RSEM

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.(none)>      3 commits

       0072121 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8bfdcff module load python in rnaseq.guillimin.ini --for htseq-count
       bb66dcc rnaSeq.guillimin.ini changed default tmpDir  --BFXDEV-144

  Joël Fillon <joel.fillon@mcgill.ca>      46 commits

       2ca8f73 Version bump
       0ce7d1a Solved conflict for merge master and bam2fastq branches
       8696c9b Changed gzip to zip to compress rnaseq de novo outputs
       d86c302 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c09c292 Changed archive command from gzip to basic zip
       63656f1 Minor fix: removed semicolumns in chipseq default config file
       f383320 Removed unused variables and functions in all pipelines
       3924dbe Added Version number in RNA-Seq De Novo Usage
       3662f2e Minor mugqic/tools version update in rnaseq de novo ini files
       9c9b747 Minor ini param adjustment + comment fix
       be76b26 Removed hard-coded path to modules.sh: not required when invoked with Bash shell
       56881db Updated default project paths with new /gs/ partition
       4c24b2b Removed deprecated GetFastaAlias.pm
       2d92ee8 Replaced shebang #!/usr/bin/perl with #!/usr/bin/env perl in all Perl scripts
       6b7b376 Fixed bug SequenceDictionaryParser filepath with environment variables + set param [annotateDbNSFP] dbNSFP not mandatory
       842ca5f Removed redundant file existence test in SequenceDictionaryParser.pm
       851df47 Check all config modules only once when config hash is built, to reduce runtime
       b5f3c9f Minor fix for adapters path in RNASeq De Novo config files
       a86bd79 Updated adapters paths in RNASeq De Novo ,ini files
       9d73b18 Added  [Warning/Error] in msg display
       d11cb76 Removed old lib/LoadModules.pm
       1126cd9 More minor bug fixes for getParam validation
       670c892 Fixed minor getParam validation bugs
       e002a7a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       872d712 Added validation in all getParam functions
       d9621c0 Fixed merge conflict in lib/Picard.pm
       88411ce Added moduleLoad validation + major style reformatting
       c130dc9 Added && between job commands to get right exit status
       080410a Added param definition validation and module availability in LoadConfig
       2ee8ea0 Added raw read dir and adapter file path validation in Trimmomatic lib
       a461511 And more and more about parallel normalization
       c32e316 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into sample_norm
       8b316d0 Fixed metrics parent step
       f6e131b Even more parallel normalization
       da10fb1 More parallel normalization
       5b5522c Normalization parallelized by sample
       7ef8faa First draft of resume-jobs
       39c49a5 Use module mugqic_dev/R until properly deployed
       813d197 Fixed typo
       3d51ae6 Fixed normalization stats file name
       b667299 First stable RNA-Seq De Novo pipeline version
       de50a89 Updated blast results and DGE matrices file names
       04e1d5b Merged conflicted rnaSeqDeNovoAssembly.pl
       d8717e4 Added POD documentation + fixed bug blast-longest-transcript file
       cdc15b6 Added BLAST for longest transcript only, with results header
       df03b3b Added metrics and deliverables steps in RNA-Seq De Novo pipeline

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       90ac5a4 BFXDEV-133 incompatibility between /usr/bin/perl and mugqic/perl/5.18.2. Added perl HOMER_HOME/bin/ to the program execution in peak annotations and motifs usign Homer
       5e62ee0 BFXDEV-133 detected bug dependencies between trimming and alignment chipseq pipeline
       972c98b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0fe7a1c BFXDEV-133 updated software, parameters, corrected bugs in chipseq pipeline for guillimin phase2
       2401cf3 BFXDEV-133 adapt chipseq pipeline and configuration file to guillimin phase2

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      1 commits

       6a2eca8 README.md edited online with Bitbucket

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      17 commits

       7f02064 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ba4cfb0 bug in ini files: the following variables are not defined : genomeSize, annotation distances, markdup, report variables
       942ef2a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       050cc3c BFXDEV-123 add design example files to chipseq and rnaseq pipelines, update the user manual
       599e532 BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       4006f00 BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       8290e24 BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper - typo
       b79e9f7 BFXDEV-36 Generated PODs documentation for rnaSeq.pl wrapper
       c1873ce BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper
       5bfe14e changed my email by johanna.sandoval@mail.mcgill.ca in standard ini files
       9d2232b BFXDEV-85 added flagstats/ generated a file with number of reads after filtering/mark as duplicates
       1d069aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6680dcf BFXDEV-84 correct bug for restart when merge lanes step failed
       394b3ae BFX-1799 wrong variable initialization for genomeSize, detected when genome is not human or mouse
       5ac34dc comment skip trimming step from standard ini files, added imagemagick to mammouth ini
       e2b57e1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8c0e000 BFXDEV-77 bug in merge metrics, instructions to run Chipseq in the pipeline directory

  jtrembla <jtremblay514@gmail.com>      15 commits

       8d91c82 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline Merge latest change prior to modifs to PyNAST (unset LD_LIB).
       05f5335 Added unset LD_LIBRARY_PATH to PyNAST step. BFXDEV-31
       8b79904 Updates to pacbio stats step. BFXDEV-30
       956d5e3 Fixed pdf report for nc1. BFXDEV-31
       c578e83 updated ini files for RRNATagger. BFXDEV-31
       abb4b98 Added module load openmpi to PyNAST. BFXDEV-31
       0eeea2e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cd179fc Added nozzle report for RRNATagger_diversity. BFXDEV-31
       22b7fea Removed Iterator::Fastx from the wrapper as it is not even used here. BFXDEV-31
       884627f fixed file checking for restart mechanism for sub loadPulses. BFXDEV-30
       d275e03 fixed file to check for input in sub referenceUploader. BFXDEV-30
       cab6153 Minor fixes to restart mechanisms. Removed inputFofn for filtering step and loadPulses step. BFXDEV-30
       e9c1d81 Corrected some parameters for celera assembly step. now on lm2 by default. BFXDEV-31
       752dd0a Implement module load and getParam checks. BFXDEV-31
       8afc93f added missing path for abacus ini files. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      33 commits

       3cf8824 fixed path to primer files. BFXDEV-31
       60e8f92 fwd and rev primers options now optional. BFXDEV-31
       7f86972 Added modules for R and mugqic_tools to rarefactionPlots.R . BFXDEV-31
       add6ebd mugqic.done fixes to rarefaction subroutines. Added mugqic tools module to appropriate subroutines. BFXDEV-31
       72338c7 Updated help screen and removed appending ./scripts/ to PATH in curr dir. BFXDEV-31
       6a982c8 Removed scripts/ dir and moved it to mugqic_tools
       c99e920 minor fix for tmpdir definition. BFXDEV-31
       8e3af42 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f806248 minor fix to ini file (for mummer). BFXDEV-30
       dabb656 Modifications to README and added mapping file example. BFXDEV-30
       d3473a4 README for RRNATagger. BFXDEV-31
       82ad834 Added RRNATagger (16S/18S/ITS rRNA amplicon) pipeline.  BFXDEV-31
       a9eeb68 fixes to ini files (pacbio pipeline). BFXDEV-30
       a6bc0cc put only 1 mersize (14) in the gullimin ini files. BFXDEV-30
       7bcfc65 Fixed bug with fofns. Minor modif to main loop. BFXDEV-30
       126618d Changed /dev/null in the order of commands so no empty consensus.fasta anymore. BFXDEV-30
       e9d2882 Forgot a && before gunzip in variantCaller  step! . BFXDEV-30
       907be7d Fix for uncompressed consensus.fasta.gz. BFXDEV-30
       d5e5250 Fixes for compatibility with latest version of pacbio assembly pipeline. BFXDEV-30
       fea30cc Fix for .bas.h5 files. BFXDEV-30
       935aea0 Updates to PacBio assembly pipeline which now support multiple polishing rounds.
       35d5972 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       28c7b82 Compute assembly stats from polished assembly and not unpolised ones. BFXDEV-30
       ebdbae0 Fixed and updated pacbio ini file for abacus. BFXDEV-30
       920b75f changed callVariants for variantCaller . BFXDEV-30
       faf6c43 typo corr. BFXDEV-30
       0b00031 Added additional instructions in the README. and changes relative lib path. BFXDEV-30
       468184f Forgot to update this library for pacBio pipeline. BFXDEV-30
       eda60c3 updated README. BFXDEV-30
       6217fda corrected for typo gnuplot . BFXDEV-30
       89f4c68 corrected relative path of lib folder
       fdf9305 BFXDEV-30 modified parameters so they are more generic.
       27ab913 Loop/dependency fixes to PacBio pipeline. Added compile stats step at the end of pipeline.

  lefebvrf <lefebvrf@gmail.com>      2 commits

       5e46e90 necessary to honour cairo as X11 backend for R graphics with current module installation
       8707d78 vanilla will hinder reading Rprofile.site, which in  our modules will not be used to force cairo as X11 backend when available

  lletourn <louis.letourneau@mail.mcgill.ca>      27 commits

       8c107ea BFXDEV-149 Fixed the way BVATools is called for coverage.
       b839f5a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c7f0a75 Removed path test since they contain env variables
       31322fb BFXDEV-124 fixed center when using mem
       c2625cc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e420186 BFXDEV-116 Fixed reverse adapter
       4f5f9d2 BFXDEV-111 Fixed many module versions
       973c7b0 BFXDEV-114 Added R loading
       4be23c3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       419f6fe BFXDEV-109 Configured java version
       b9989e5 Updated guillimin pipeline
       7034346 Updated depth of coverage ini
       db2b3a8 Fixed case when there are no BED files
       b79a3a0 Fixed module typo
       be34921 tmp hack so nanuq takes coverage graphs
       8cf40e6 BFXDEV-89 BFXDEV-88 Change GATK to BVATools for DepthOfCoverage and support multi bed in project
       7f45e96 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b0a885e Added the missing mappability key
       5b568ba Added line to keep overlapping pairs
       73d736b Added umask to dnaSeq jobs
       7e1f584 Updated module versions
       67b9a78 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       472a3cd BFXDEV-73 Fixed undef jobId if step is up 2 date
       a6e3c7c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2d9d66a Added details to the release guide
       70f59da Version bump 1.2-beta
       daabc8d Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       3cd0dee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22a87a0 DNAseq: correct SNVmetrics dependency BFXDEV-136; RNAseq: remove hstseq dummy/null module usage BFXDEV-137
       cbdca41 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3558370 DNASEQ: correct depthofCoverage missing bed file field in sample sheet BFXDEV-135
       38db1aa ask bash to source /etc/profile.d/modules.sh BFXDEV-125
       184fd53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3cbe61a add the -T option to cuffcompare BFXDEV-93 and validate previous commit for BFXDEV-47
       a671e06 README.md edited online with Bitbucket
       f9a7616 README.md edited online with Bitbucket
       8f494a8 README.md add the bioinnfo email adress
       6cb9d24 correct the 1st line typo in chipSeq.pl
       f1df7ba Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31f3438 replace ; by && in pipelines except for the rm of .done
       8620232 replace a default ppn=20 by ppn=12 otherwise the job will never be launched

  mmichaud <marc.michaud@mail.mcgill.ca>      33 commits

       e242aa1 BFXDEV-147 Fix index mask of a single-index lane in a dual-index run
       3c527b2 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       bce0666 Use new version of BVATools with simpleChrName support
       0cf3980 RunProcessing: gentle perlcritic compliance
       6e3b82e Fix BED file list from SampleSheet
       4753244 Fix again the readsqc BVATools path
       becd68c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c8fc20f MD5 job doesn't need 5 cores, only one is ok
       7b03a8f Fix readsqc using current session bvatools jar file instead of the one loaded on the compute node
       230a617 Fix run processing when there is no target file
       cb1be91 Merge branch 'master' into runProcessing
       59829b2 Use new module loader for BVATools qc
       61c25ba Add cranR module version
       7e38002 Merge master into runProcessing branch
       1261ab2 Support spaces in bed files
       2b5b54e Don't align mice yet
       9bb906f Support spaces in bed files. Fix bvatool module version
       3e2571b Fix usage message (to show optional parameters) and print error message on dying
       f549388 Die when there is a barcode collision. Fix bwa rgCenter for mem
       172452b Add a parameter for the number of mismatches
       da5e47a Default values for both sample sheets path
       214c3d9 Merge branch 'master' into runProcessing
       9584091 Update BWA module version
       1a6989d Fix rsync command: quote character were not escaped correctly
       a8e58fd Add missing depedencies for the copy job (metrics)
       f19991c Don't create qc folder when the qc job will not run
       7385c8e Depth of coverage: add reference parameter
       7508bd4 Change '_' to '.' as seperators in the metrics job name
       2ea6103 Enforce processingSheetId column in the sample sheet only when processing a run
       f79fd8e BFXDEV-76 IlluminaRunProcessing Add target coverage metrics & change job ids to support multiple sample with the same name in the same lane
       8e90e62 Merging upstream changes of the barcode metrics jobs. Less core and memory used, skip last read
       5a2b132 Merging master into runProcessing branch
       70383b3 BFXDEV-76 Illumina Run Processing Pipeline

  Pascale Marquis <pmarquis@lg-1r17-n02.guillimin.clumeq.ca>      2 commits

       2172223 rnaSeq.guillimin_hg19.ini
       d75c72d rnaSeq.guillimin_mm10.ini

1.1        Mon Dec 23 14:23:34 2013 -0500        137 commits

  Joël Fillon <joel.fillon@mcgill.ca>      69 commits

       3051392 Commented code
       78e706c Check rawReadDir value in configuration file
       3e137a0 Minor R version fix
       2eec3e9 Added trimMetrics step, fixed blast outfmt single quote bug
       758a68b Added option file check, module availability check
       7fabc9b Added genes/isoforms length values in edgeR results
       faefc93 Added blast description in edgeR results
       c54b336 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e61eb0f Cleanup of rnaseq_denovo files
       02ed5e6 Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       378d4d9 Updated ocnfig .ini files with trim parameters
       2da3da4 Fixed missing escape $R_TOOLS
       470df00 Fixed differential expression bug
       36711d4 Added differentialExpression step
       467076b Added multi-step dependency support
       eba1405 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       003c41b Fixed Trinity.pm merge conflicts
       6046b1b Added trim step + various fixes
       efc7b74 Minor fix
       5c140d0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e31fb3c Added guillimin config file + minor modifs
       faf0d71 Updated .ini config files with cluster-dependant processing values
       19f6b6d Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       dabafcf Added blastCPUperJob config tag
       0eb0a4a Added blast step
       01b495f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       885d5f1 Updated default rnaseq tmpdir on guillimin: /sb/scratch/jfillon
       7768316 Added first version of BLAST step
       03e5423 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       0a05c4e Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       fcb3c5d Minor fix
       dc1e136 Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       a80a6c2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       29f179f Updated config ini files
       63a2a5d Minor fix bowtie CPU
       b3cf0fc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       27fd66f Massive cleanup
       f2eeb87 Another minor fix in mammouth ini file
       60fac50 Minor fix mammouth ini file
       2822aab Updated mammouth ini file
       a7f4778 Added TrinityQC step + code reorganization
       340a63d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       9d7d75e Moved normalization parameters in config file
       950b713 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       3fd9531 Fixed bug edgeR semicolumn
       f850faf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       c9e35cf Added edgeR step
       42d400d Added rnaseq_denovo.mamouth.ini
       a20391f Fixed bug unsorted list of fastq.gz from find
       4a35ad6 Extended wall time for normalization and trinity
       73bc9e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       db70b57 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       02f65b8 rsem prepare reference separately + job log relative path
       ef38d5f First version of RNA-Seq de novo pipeline with normalization, trinity, rsem
       39a95e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       74d4ac0 Further development of RNA-Seq de novo assembly
       92aa251 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e43a3e2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       8946fc1 Further development of rnaseq_denovo pipeline
       06fb9bc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       d2a0332 Import of old deNovoAssembly in rnaseq_denovo branch
       fc8082b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       052b804 Minor fix
       eddaed2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       ab7f621 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       9eb1383 First draft of de novo RNA-Seq normalization
       f4c998a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       5873297 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       23f34b2 Fixed bug missing '\' before

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      12 commits

       29ff238 Readme for chipseq pipeline
       6f007a9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       062d698 document mugqic pipelines setup using md - correcting bugs
       30364fc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       fad747c document mugqic pipelines setup using md - correcting bugs
       5a0a009 document mugqic pipelines setup using md - correcting bugs
       f771f1a document mugqic pipelines setup using md - correcting bugs
       848249a document mugqic pipelines setup using md - correcting bugs
       58dca89 document mugqic pipelines setup using md - correcting bugs
       9652290 document mugqic pipelines setup using md - correcting bugs
       4cf001d document mugqic pipelines setup using md
       4e9a82a Miscelaneous graphs chipseq pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      21 commits

       a5f8c70 corrected bug in qcTagsjobid
       d3b16fd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ecd9434 Added PODs documentation to chipseq pipeline wrapper
       3555ab4 correcting links to project's directories on README.md
       1b94392 link to wiki on general README file
       31db440 BFXDEV-63 Added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       589689e added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       eb5a068 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       4e1e7df added variables for chipseq peaks annotation plots
       ea839b3 corrected bug loading ReadMetrics and read metrics sample/run separator
       9298aae Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       0d69c52 commit chipseq report branch
       95aab1a added annotation statistics + graphs to the chipseq pipeline
       94bcb8b adapt chipseq pipeline to the new resume jobs functionality
       73ee7a7 adapt chipseq pipeline to the new resume jobs functionality- config file
       f1dabdd adapt chipseq pipeline to the new resume jobs functionality
       4fe9953 Adapted chipseq pipeline to resume jobs changes
       9d6f3d8 Transfer chipseq pipeline metrics libraries to the pipeline space
       df472f0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       2e5f613 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       67ae379 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report

  Johanna Sandoval <sandoval@ip03.m>      2 commits

       2a5fbe0 document mugqic pipelines setup using md
       6e68220 adding mammouth configuration file for chipseq pipeline

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      8 commits

       1971ef1 Added memtime to dcmegablast and blastdbcmd. BFXDEV-30
       9f0b66a Fixed unwanted mofifs to BLAST.pm. BFXDEV-30
       b430c4c Modifications to perl packages for Pacbiopipeline. BFXDEV-30
       bb5a19e Major modifications to the pacbio pipeline. Functional version tested on abacus and guillimin. BFXDEV-30
       a17cd5a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       97cd35a Added readme for pacbio assembly. BFXDEV-30
       6432aad Remove my info fields in ini file. BFXDEV-30
       67c11c1 Added PacBio assembly pipeline libs/wrapper/files, etc. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      19 commits

       4a6d8af Version bump 1.1
       ef0d743 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2541c18 BFXDEV-68 Added mutect to paired pipeline
       d3998b6 Removed gpfs for guillimin
       ed3de9b Fixed conflicts
       11d0fe5 Fixed undef on steps that are up2date
       bc77881 BFXDEV-67 use a true csv parser
       bdbb7ce Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline
       4d7c7fe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1a8cf59 BFXDEV-45 Fixed single read handling, fixed undefs
       5d28000 BFXDEV-15 Changed the name of lane level bam metrics
       ed46494 Fixed starting dependency and final job output
       4555ab0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9252050 Added pre dependency and final sample output
       aa452ff Changed BWA version for run processing and added some metrics
       a29bac1 Added runprocessing parameters
       ab65932 BFXDEV-52 DNAseq realignment restart generates duplicate unmapped bam
       05272a8 BFXDEV-48 added missing close BFXDEV-49 added support for alignment+QC only RNASeqQC
       e1ecac6 Version bump

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       da9919d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fa47c13 replace ; by && - BFXDEV-47 and update the rnaseq ini files
       e0ccc02 dnaSeq.pl correct typo (realign step) : BFXDEV-54
       371c03f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       799fca2 change/test replacing ; by && un command line
       0a07178 Trimmomatic.pm remove single trimm bug

1.0        Fri Nov 8 15:03:24 2013 -0500        794 commits

  David  Morais <dmorais@ccs.usherbrooke.ca>      2 commits

       2620120 Merged in daveM_denovo (pull request #2)
       0764a25 Merged in daveM_denovo (pull request #1)

  David Morais <moraisd2@ip03.m>      94 commits

       4c13dd1 remove entry
       f458e98 added new entries
       d8e671b split butterfly command lines in chunks
       a963e53 creates BLAST best hits reports
       61d9bf3 Fixed bug
       5b9c160 Modify how the butterfly files are splitted
       0234e16 fixed bug
       5b2bbbc modify command line
       f29fe2e fixed bug
       f2bf4ab Added full path to groupInfo left and right samples
       e8794a0 Added option for de Novo Assembly
       3490e43 Implemented and tested Merge Trimming Remove duplicates de Novo assembly Blast BWA
       8094aca fix paths
       191bdcc fixed picard version and env variable
       17934d1 change module add trinity
       891bcf9 change module add blast
       60bfa9f change module add blast
       78f4ba6 change module add trinity
       59905eb fixed typo
       aacc3ec fixed typo
       e78b5b6 Modified HtseqCount::matrixMake call. Now group is the last parameter so it can be left out.
       c13eb91 Fix $group problem. Now if group is not specified the variable $group takes $sampleName value.
       dcabd4b Fixed blast db error name
       592659c perldoc
       a94ec1d adding scripts needed by de Novo Assembly
       a64a2ce add comments
       289ec2f De novo Assembly config file sample. It needs to be modified according to your each project.
       ba64039 mammouth_deNovoAssembly.ini
       befb424 This is the main de novo RNA assembly pipeline. This is the first commit.
       774254d modified config variable
       5be92fd tidying up
       7d190c5 tidying up
       8d9704b tidying up
       03fc63e tidying up
       9b75e54 tidying up
       24e31b0 tidying up
       2cb0b70 tidying up
       94a0d32 tidying up
       6a8f8fc tidying up
       24aa7af tidying up
       42608a8 tidying up
       2d3f8a9 tidying up
       23e5b5e Coded full command line. Added sub contigStats
       d0f849d Modified $laneDirectory
       88bb2c7 Modified $laneDirectory
       db5f506 Modified $laneDirectory in the singleCommand
       bc774be Modified $laneDirectory to assembly/ and added $group variable to deal with groups on deNovoTranscriptome assembly.
       55eb46d Modified perldoc
       0ad8b58 Modified $laneDirectory
       c21cab5 Modified $laneDirectory from alignment to assembly
       2037ad6 Modified $laneDirectory
       6bb3964 Modified $laneDirectory
       229bba4 Assingned reads/ to $defaultDir
       f0b1c61 Assigned  $laneDirectory to 'reads/'
       eb94920 Modified $laneDirectory = reads/
       921adb0 new line
       b187d45 Output bast results to /db directory
       7c1754a change from own branch
       c58b61f First commit. Library to create differencial expression analysis.
       4086a94 Modified Doc
       cc24e6b First Commit. Library to generate basic statistics on raw read count. Not Yet Ready for use.
       1df2ca7 Library to create generate basic statistics of each sample. First commit. Not ready for use yet.
       edf8193 Added the possibility of looping through more tha one DB. In This case the DB must be passed as an argument to the function.
       3dac1f1 added comment
       fccb0b1 removed quality offset from file name
       db155a0 removed quality offset from file name
       301427d added indexing function
       e373c99 add function
       34fd4a4 add function
       986c00c add function
       7183abe modify hash groupInfo
       a77802b fixed input values
       6bae877 BLAST lib, first commit
       78008bd modify looping by group option
       998b38a tiding up
       6e243f5 fist commit Trinity
       fb5182c change split setting
       719e16b added new entries to sampleInfo hash
       4bb98a4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       644271f  a simple multiformat file Spliter. First commit
       8d3ccb3 Remove duplicate reads. First commit
       cb8cfdf modify to work with nanuq sample sheet
       56d16c8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       80dd2fd First commit
       7fbacf8 first attempt on STD
       736365b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2c33bfc Reads the libraries form script_dir/lib
       c3d669c added function that allows the script to read libs from the script dir
       d9a378d add header
       e32573a remove package folder
       754f187 first commit. Lib to read module list
       16372cb first Commit. lib to read config files
       10e24cd change repo name
       c56a0af adding space.

  eric audemard <eaudemar@imac6-ub.(none)>      5 commits

       f536084 add software install script : igvtools virusFinder SVDetect add genome install script : virus
       2bd3239 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f55f728 add install canFam3 and new bwa
       7bdcac7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       270bac9 adding deNovoSV pipeline from the mugqic_eaudem

  Francois Lefebvre <flefebvr@abacus2.(none)>      4 commits

       6f6477c febvre <ddcccccZZZZlefebvr@abacus2.(none)>
       a4f6031 Integrated exploratory analysis step
       1ea865b Changed Rscript -e template
       a3468dc Fixed perl formatGtfCufflinks.pl call

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       d312e26 fixed problem when R packages list has duplicates
       3fae3d0 Merged in module_install_scripts (pull request #3)

  Francois Lefebvre <lefebvrf@gmail.com>      53 commits

       2ec42e9 duplicate section names
       2a175b0 duplicate sampleOutputRoot
       e9e50a1 duplciate param in ini
       2abd12c duplicate sortSam sections  in ini file removed
       e8f8792 fixed bwa module name dnaseq.abacus.ini
       5b3c1a9 kmergenie install script
       810d081 updated rnaseq.abacus.ini  wiggle chrsize path
       8321d0f Updated R install script
       d3087f7 R install script now installing Vennerable, gqSeqUtils, gqUtils from remotes.
       9caf202 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       b6591f3 Exploratory analysis step, still to test
       ac2361f Added RSEM install script  + updated Trinity one
       789cb0f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       6e37c42 no message
       b014a95 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       7d5e0aa Changed bwa to mugqic/bwa/0.7.5a, problem between RNA-seQC and 0.6.2
       165e72d Updates R/Bioc script, now also installing all org packages
       aa3dcc8 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       2e48fa5 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       92e196a no message
       5d60cfa Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8d4ef87 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       5124369 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       34000ad hg19 install script" add gtf mask for chrM, tRNA, rRNA for cufflinks
       77f6acd Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8165a03 Walltimes for htseq, cufflinks, cuffdiff too small  for typical data
       d99c877 htseq walltime
       03aa1e1 Added 72h wall time align step
       f70ac2f Updated bwa install script
       5556091 Added blast module install script
       a4cd056 Added few more R packages for install list
       baaacb2 Removed MANPATH prepend from mugqic/R defnition
       8eb8a25 chmod in in mugqic_tools.sh install script.
       42d4ef9 <tools,perl-tools,R-tools,java-tools> from the svn now in tool_shed/
       9b773f4 R packages lists updated
       0725e01 Corrected bash_profiel for mammoth and changed R install script to 3.0.0
       83cc87f Added R package "outliers" to R packages dependencies master list
       5c1d4a6 Updated mammouth $MUGQIC_INSTALL_HOME value in guessing script.
       bea7789 Few change to abacus wall time rnaseq
       0e30bd7 Previous fix to Metrics.pm did not work
       30f0bfb minor fix to Metrics.pm (will be depr. anyway at some point). chrome size file for wiggle left to default in default abacus .ini
       16132ed Added -T sort unix sort in metrics.pm to correct problem on guillimin. Added job name suffixes in rnaseq.pl to make job names unique on Guillimin (dependencies)
       60c84cc Merge remote-tracking branch 'origin/rnaseq_test'
       ca2a4ac ..
       59e116e Added hg19 installation
       c83481d Added mm10 installation
       c0dabe4 no message
       c99255e updated to top hat 2.08 (bug in 2.07) + more genome scripts
       57248b5 Added Eric's abyss install script
       b20104a drafted genome installation script
       e82f783 Corrected htseq module by adding module load python. Also started genome setup scripts
       a4c4232 Finished python/numpy script. This script will not be 100% portable, need to set locate BLAS and LAPACK
       6da084c Added modules/ directory to hold modules related content

  Joel Fillon <fillon@ip03.m>      4 commits

       c474058 Added java module in BWA lib + dos2unix rnaSeq.guillimin.ini
       5bf363d Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       929fb10 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       5f3e86b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output

  Joel Fillon <jfillon@abacus1.(none)>      4 commits

       a051fa5 Added readRestartFile  function
       a7a1404 Print out MUGQIC command exit status.
       f55bec0 Missing ";"
       a7a2f6d Missing ";"

  Joel Fillon <jfillon@abacus2.(none)>      1 commits

       16641e8 Minor misspelling

  Joël Fillon <joel.fillon@free.fr>      27 commits

       67995c4 Added simple Perl script tool_shed/getLogTextReport.pl to create log reports
       f7ddba6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       3d61e71 Simplified getLogTextReport parameters
       e12feca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       aac70e2 Print out full job log path + 2 exit code outputs
       d9e2d1a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       75de931 Fixed lib/SubmitToCluster.pm conflicts
       253bdfb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       e7d552e Added number of jobs in log output
       020e4e9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       8ac8484 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       4991e38 Added timestamp in job log filenames .o
       ae0b820 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       b42133b Reorganisation of job logs in an array
       2db5f57 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       5d39c04 Added a few comments
       bc9df34 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       f3998be Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       30cf022 Fixed exit status + walltime
       39079a5 Added getLogTextReport function
       3c10f40 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       964e92d Fixed package name
       4c6b06b Merged master into logReport
       5e38b43 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       70eff6e Added readRestartFile log function + Exit Status code in .o
       d7f0f17 Fixed deliverable typo
       5fd5004 added python tools path in mugqic_tools.sh

  Joel Fillon <joel.fillon@mcgill.ca>      8 commits

       81a3ec9 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       e19c7bd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2128dd2 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2235fc5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       dab8f29 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       4cf3157 Check well-defined variables
       7b03705 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c4eaaeb dos2unix rnaSeq.guillimin.ini

  Joël Fillon <joel.fillon@mcgill.ca>      95 commits

       a3fddd2 Removed Log lib (now merged in mugqic_tools/getLogReport.pl)
       32f25ea A bit of a cleanup
       8dbabd9 Synchronized rnaseq .ini files
       c8fe2ba Fixed bug set @INC with relative path to mugqic pipeline lib
       06ca9c5 Moved module and genome files to mugqic_resources repository
       3f11508 Fixed bug missing '\' before
       b163546 Standardization of rnaseq .ini file for the 3 clusters
       6e7375f Removed prereq in module install script template
       d76adbf Removed prereq + fix R module load to compile kmergenie
       784d0c7 Updated Picard module install script now based on template
       750caca Fixed inline comment bug in module install script template + fixed permissions in ea-utils module install script
       6fc665c Fixed ea-utils compilation error on guillimin
       24ab3fe Minor comment changes
       dd1f730 Added ea-utils module install script + minor README change
       a8533c2 Renamed MODULE_INSTALL_SCRIPT_TEMPLATE.sh + deleted old archive + minor fix
       e6a1552 Added module file permissions
       7e85ddc Created a module install script template; modified tabix install script
       1b2f9bc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       40b98e7 Deleted redundant modules/add_me_to_bashprofile.sh; added README.txt instead
       e00d040 dos2unix all pipeline .ini files
       1467393 Updated chipSeq pipeline .ini files with latest module versions
       0133053 Updated rnaseq .ini files with latest module versions
       ce59f4f Merged
       e187b2d Minor change in R module install script
       88ebc06 Changed permissions in R module install script
       d5cfef1 Added vcftools module install script
       8a5b369 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a788a15 Added archive storage + permissions + cleanup in snpEff module install script
       127d327 Removed unnecessary lib BiocInstaller for R install
       f1bc1a0 Reorder chipseq modules in .ini files
       ecb3c6b Removed duplicated modules in chipSeq .ini files
       4280f36 Added archive storage + permissions + cleanup in Trinity module install script
       1976a51 Added archive storage in UCSC module install script + update dnaseq ini file with bwa/0.6.2-ptx renaming
       dfc61b0 Added archive storage in Trimmomatic module install script
       64e9a66 Added archive storage in tophat module install script
       04ba4a9 Added zip archive storage in picard module install script
       d1c2d54 Added archive storage msg in GATK module install script
       eed094d Added archive local storage in igvtools module install script
       8bdeb5b Removed bedGraphToBigWig module install script (now part of ucsc module)
       ff546a7 Added permissions in Tophat module install script
       9c25d16 Updated dnaseq/validation.pl and dnaseq/pairedVariants.pl to comply with new version of SubmitToCluster
       cbaf966 Added permissions in Picard module install script
       7a726e6 Added permissions and cleanup in IGVTools module install script + minor changes
       3709582 Added permissions and cleanup in hmmer module install script
       34df850 Added permissions and cleanup in gemLibrary module install script + minor aesthetic changes
       89ac0f9 Added permissions + minor fix in gatk module install script
       ed0aeca Added permissions and cleanup in fastx module install script
       d7dbb81 Added permissions and cleanup in fastqc module install script
       bd89096 Added permissions and cleanup in exonerate module install script
       ff61379 Fixed permission bug in UCSC module install script
       568ef4e Fixed another permission bug in module install scripts
       daa7cf2 Fixed bedtools install script bug for archives with different naming system
       c3a3067 Added permissions and cleanup in cufflinks install script
       1cb9cc5 Added permissions and cleanup in bwa install script
       3ef584d Updated default pipeline .ini files with new mugqic/ucsc/20130924 for bedGraphToBigWig
       942bc47 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0eac0a1 Removed -j option in make; added kentUtils.Documentation.txt in bin/
       e2253bc New install script for UCSC genome browser 'kent' bioinformatic utilities
       7e1b165 Reorganized install scripts for a5, bedtools, blast, blat, bowtie, bowtie2
       d097e82 Reorganized BEDtools install script
       99f688c Removed last workdir parameter in chipSeq, rnaSeq; removed commented code; fixed blast+ bug in blast install module
       a985b24 Added shebang in pipeline bash
       27c480f Removed leading ':' in chipSeq JOB_IDS lists
       669a97d Clarified bash by using job variables + header
       6e72b67 Removed sampleOutputRoot tag in all .ini
       ce93e25 Removed 'mkdir output_jobs' commented lines
       89d60fd Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       68e948a Removed output_jobs subdirs in chipSeq pipeline + added job dependencies in log
       c47fd5d Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       454e4ab Fancy install script from Johanna
       e9f61b1 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       054fbfb Removed final ';' in MACS2 cmd + aesthetic change
       6e049a9 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8471e38 Removed ';' at the end of Homer cmds
       ae79f34 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       66a9193 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2445956 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       310035f formatting changes + chipSeq redirect job output to specific location
       830d57b  ->
       19f116e Minor formatting changes
       890c915 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       995a8b4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       a05a090 Added java module in BWA lib
       622217d Reorganized job_output for DNAseq pipeline
       405043d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8f6b4e5 Added tophat 2.0.9
       6ef6c34 Minor formatting changes
       e01c182 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       9e220c5 Minor job parsing change
       4af616a Added Job ID number in log
       7ab1321 Redirect job list output into a specific file
       7438d64 Job outputs go into jobs_output directory
       967b888 Minor coding standardization changes
       472a9b2 Minor error msg edit
       abf1067 Moved getLogTextReport.pl in tool_shed/perl-tools/

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      49 commits

       8e63db4 changed config files to fit the new SubmitToCluster:initSubmit structure
       98b7290 changed config files to fit the new SubmitToCluster:initSubmit structure
       f2c2aa0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b56095e adapters small rna from Trimmomatic-0.22
       61c25df Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       85c99b1 added module and installation script for gemLibrary
       3998d3d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       15e6a92 added chipseq.guillimin.ini parameters for hg19
       152e3a7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7559eaf changed deseq.R : added flag comment.char=" " when reading edger_results.csv
       b641ee2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e46c9c1 update install genomes script in tool_shed directory
       9849932 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5bde35 prepend perl-tools directory to PERL5LIB
       e819469 adapted ini file to changes in trimming parameters (headcrop)
       af929ce Changed readstats module, added ReadMetrics.pm library
       ef32d66 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       db61b39 Parsers or trimmomatic and flagstats reads statisttics, create a readstats.csv file
       a8d8917 remove unused function (old tag directory generator)
       dfc0b2e compute read stats using the ReadMetrics library
       7faa225 prepend perl-tools to the PERL5LIB variable
       b4c6417 corrected freeBayes install on guillimin
       79ee5c6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       e3e8549 Adapt install freebayes for Guillimin
       f46883d Generate QC stats using R
       3fa1222 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       09f7ec2 tag directory generation removed tmp sam file creation - added samtools
       c3ca36b tag directory generation removed tmp sam file
       d045ac1 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8288285 Freebayes install and create module
       237ed03 filter aligned unique reads using samtools
       176e710 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       df89fcb added qc plots, corrected motif weblogo bug: seqlogo available on v 2.8.2
       51077d3 generate chIPSeq QC metrics R script
       6087d82 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a99e139 Increase walltime for qctags generation
       8e47a39 Initialize BWA_JOB_IDS per sample, create design 0 directories for MACS and annotations
       143492d generate wiggle tracks: remove -i parameter
       9666a94 Validate differences between design file and nanuq project file
       c606dbd Validate missing values in design file
       0f2715e adapt read statistics to changes in Metrics library
       494d67e pbam flag for paired end reads on MACs
       d0cfe12 chipseq pipeline configuration file for abacus
       ed0bf0b chipseq pipeline wrapper
       2665867 chipseq pipeline HOMER tags, annotations and motifs
       dde6603 chipseq pipeline MACS2 peak calls
       69e8682 uncommited changes, will control variable names from my scripts
       aa273e9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       70afa15 transform jobIdPrefix to avoid illegal bash variable names

  Julien Tremblay <jtrembla@abacus2.(none)>      2 commits

       fde3058 modifications to skip unnecessery bam merges.
       cbdfb71 Added dnaclust install trace

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      3 commits

       4edd2fd Fixed bowtietophat module, added cpu param to align, drafted template guillimin
       136c03f test
       d75425a Added Config-Simple and guillimin ini wc

  lefebvrf <lefebvrf@gmail.com>      5 commits

       c7b8ab6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       63e674f quiet option  is hard coded, makes it hard to diagnose cufflinks… -q could be enabled in previously added option parameter
       bfc8e65 Cufflinks otherOoption parameter, to allow for instance for -M rRNA, MT, etc mask
       62b1928 abacus htseq wall time increased from 3 to 12h
       4ed3cba unused raw_count folder was created

  lletourn <louis.letourneau@mail.mcgill.ca>      174 commits

       b07568d Added restart implementation
       6633b73 Fixed pairedVariants with new structure
       50e3e85 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       29fc49a Added release instructions
       5833330 Removed deprecated code.\nAdded a version
       6294737 Fixed step order bug
       3a50453 Fixed region bug from previous commit
       490b11e Added a way to skip trimming. Added a mem example
       a385b07 Added a way to skip trimming
       c60ac77 Merged master
       a535288 Added lib barcode to lane data
       441b96b Added missing ;
       8fa3ded Updated old code
       399a2d9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9e7894 Fixed mem bugs
       1c2a330 Fixed badly named CFG key
       df22d10 Fixed paths for install. Crucial for these to work since they hardcode the paths
       fbe736a Merged Master, some fixes and new dnaSeq report
       9bcc238 Updated dbsnp
       623de3e Amos install
       0ffb61e Fixed region for torque
       abd0940 Amos install
       0370f0e Added split for stats
       508f191 Fixed if test
       2202926 Fixed dependant files
       57f85c2 Fixed some bugs found while testing public datasets
       e9d5b4e Fixed dir search order
       dab8a76 Put appending and cleaning of dones in SubmitToCluster. More central, easier to maintain
       2cd0d6b Fixed to support any year 2XXX
       b9aafe9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1eb4c13 Merged master with outputDir fixes
       ac70626 Merged HEAD, removal of job output dir
       a4c2162 Change delimiter for regions, msub doesn't like colon
       d1c131e Don't generate the pileup, it takes way too much space
       9805b64 Added hg1k install and GNU parallel
       2be9b17 Forced default picard buffer to the abacus optimal one
       6230215 Finished updated the rnaSeq.pl pipeline
       431b192 Merged master, added java module, dos2unix rnaSeq.ini
       b2898cb Merged Master
       cf0fca3 Finished implementing reset in dnaSeq
       a2fc84a merged
       fcd41d2 added missing calls for output dir
       708959e added sampleOutput param to inis
       ef0f696 Use the reight blat tools
       5fdb7f3 installed v4 preview
       8336865 Added beagle 3.3.2
       ba10cb4 Refactored in the new Job object and input output file tests for restart
       ea2786e Updated the way snp calling works to support, nb jobs instead of window.
       f28a83a Removed dup index generation since we added recalibration
       4a6c701 Fixed the way output directories are initialized
       10f160f Updated version
       0f80b86 Added pigz
       b9e96e5 Merged master
       9845505 Added PipelineUtil pacakge to all.\nCompleted BLASt implementation
       d218799 Implemented per-job indel realigner Implemented multi sample snp calling Implemented per-job/per chromosome snp calling
       3359fa3 Use ini for igvtools version Updated igvtools version for correct exit code
       a3860f4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       608ae55 Fixed already processed lanes.
       66430db Version bump fixed execute bit
       4ed8abc Trial implementation of new job restart flags
       2e09503 Implemented restart job handling
       9912c46 Added link to ptx patch
       c984117 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       daa8f5a Added and updated modules
       5448d09 Added tabix
       35fa5c7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       75c9afd Fixed many bugs
       ae12147 Fixed many bugs
       fb3af3a Fixed hash generation
       42f88d8 Changed trimming
       e7eb18f Fixed index vs ref problem
       1dbc2b8 Fixed index vs ref problem
       84c1895 Split metrics...we probably should split it further.
       7747bcd Added TRINITY_HOME
       70cd9b8 Added a parsable file containing trimming stats
       d032c5a Fixed default params
       9ea1211  Changed metrics position in steps
       96ca681 Fixed parameter passing.
       ea8cecc Fixed problems in RNASeq denovo pipeline
       5cf8898 Added support for non-multiplex lanes....
       1d8a190 Fixed trimmomatic bug Added usage to sampleSetup Fixed some ini params in rnaseq
       a405492 Added options for skipping samples
       4d7e1fa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d06a44f Fixed many issues and params for the RNASeq denovo
       162142a Format stats from new format
       7b61961 Version bump
       972fa53 Added missing file warning
       42b94f2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c347b1f Made many fixes to the denovo RNA Assembly pipeline
       a995968 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       05042bf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7ec9500 Modified permissions
       ceea8d2 Fixed mpileup bam location
       58152cf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       49b1409 Fixed parallele threads in single reads
       e311d8d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d7e10f5 Fixed bwa mem RG bug.
       8fa5378 Updated version
       01e0dbd Added gatk to module list
       8ff53e4 Fixed special single issue
       1c6edb2 Added recalibration
       36cbabd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f027679 Fixed generation for single ends
       660ceda Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8ef1373 Fixed validation for new project layout Added BWA mem support
       89b8a33 Added snpEff and mutect Changed mode
       070bcee Fixed output redirection
       add31c7 Fixed multi runId resolution
       351be7c Cleaned up downloaded genomes
       40edd6e Sample setup scripts that works with miseq and hiseq and fetches the nanuq project by ID
       69f8ad4 Added gallus gallus 3
       a408a59 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b3ee340 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f461b67 Fixed default ref.
       210079b Use parallele bwa by default
       1d973c6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e83d068 Fixed libs and dna pipeline to use official project hierarchy
       b0a33ed Fixed cleanup code
       2af9a96 Updated the ray module
       643a97b Added the exonerate module Fixed the ray module
       5c1c619 Fixed append bug
       41a888f Change file mode
       12bdf75 Fixed configurable module load to trimmomatic Fixed pass output dir, not hard coded in trimmomatic
       3f20c53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46ac74f Added missing rawReadsDir
       6f4161d Partially fixed raw_read support for dnaSeq Fixed markdup call in dnaSeq
       f861134 fixed conflicts
       c973658 Fixed whole genome coverage
       32b1e16 Added module versioning, fixed paired pileup bug
       94c6f97 Added SV to the paired pipeline
       2e99fd6 Added snv calling
       2bfad55 Rolled back changes because it broke existing pipelines
       0bc8182 Added SV to the paired pipeline
       d67cedd Added DNAC support
       79e909f Added DNAC support
       b9b11aa Merge branch 'mugqic_lletourn'
       0cceb49 Use job ids Added variant calling Added more flexibility when calling modules
       18a60ab Added thresholds to genomeCoverage
       a70bf88 Merge branch 'master' into mugqic_lletourn
       275d08a Fixed usage
       cee2b18 Added for snp calling
       4e93043 Added metrics and steps to the dnaSeq pipeline
       8222088 Change job dependencies from names to ids.
       77a4c61 Fixed 2 paths bugs
       30396f0 Don't use alternate contigs
       0c43b15 Changed executable attribute
       7b04c50 Keep dictionary ordering
       f08a3ea Fixed usage message
       4117265 Added paired snp variant calling
       a5d2e16 Added more flexibility with tmp.
       67b2b05 Completed first half of the validation pipeline Added adapters file for common paired adapters
       d0df7f4 Added metrics to validation Added IGV tootls Added targetCoverage
       d66f466 Fixed header
       466071a Fixed issues with multi single+pair
       df02765 Added Validation pipeline
       b48e4b5 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       01ce3d4 Fixed index generation
       0f46ad3 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       0067ebf Fixed timestamp checking
       08d14e0 Added trim counts statistics
       0923113 Fixed Trimmomatic MAJOR bug Added more settings to trimmomatic
       275cd84 Use the right path for trimmomatic
       f652072 Fixed output job directory naming
       cc07363 Fixed bad realigner output file name
       7c15b68 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       992ef66 Fixed bad argument
       c4479ab Added guillimin conf file
       5411b17 Fixed job dependency code.
       53a75c0 Finished bam generation pipeline
       faf8ed1 Merge branch 'master' of ssh://bitbucket.org/mugqic/mugqic_pipeline
       4e30674 Added realigner
       eda6d12 Fixed uninitialized bugs
       2b17806 First pass for dnaSeq pipeline
       5b1fc8b Fixed typo

  Louis Letourneau <lletourn@abacus1.(none)>      1 commits

       9286b83 Fixed many bugs

  Mathieu Bourgey <bourgey1@ip03.m>      2 commits

       3b2af18 RNASEQDN: modify trinity to check for previous assembly + variable name change
       d49a53a TOOLS: chenge permission of non-executable scripts

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      129 commits

       a51199f samtools allow pileup with nb de region = 1 => pas de region
       ea59d1f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3e902cc dnaseq.pl && samtools allow pileup with nb de region = 1 => pas de region
       0ff4f48 upadte rnaseq.pl: missing convertion to new job object for stranded wiggle printToSubmit calls
       ae5a4a8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f79aadb dnaseq annotation bug correct and ini update
       d4e4586 rnaseq.pl: correct exploratoiy dependency if start at step 9
       1b1abe4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       33df28e dnaseq correct metrics restart bug & update dnaseq ini
       6db6f70 remove tools_shed for the new tool repo
       c653f58 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e1bc285 metrics/SNVGraphs add input output test for restart
       72b31c8 update Tools
       0e1b10e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b758404 update cuffdiff stranded
       af5dba2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       371205c switch ToolShed to Tools
       c0c1d65 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       665a098 RTOOLS update
       3887f71 snvGraphMetrics.R  update
       38b94df DNAseq.pl update for metrics
       4732109 DNAseq.pl update for metrics
       daceb8a Merge branch 'dnaSeq_Report'
       c82147a DNAseq.pl update for metrics and report
       b31654c DNAseq.pl update for metrics and report
       75f28e5 DNAseq.pl update for metrics and report
       425b21e add tool_shed/python-tools/vcfStats.py
       ac96abb remove mater to branch conflict
       d133749 rnaSeq.pl edited online with Bitbucket
       23e60fc remove conflict between danSeq_report and master
       e67ddd8 add tool_shed/tools/splitSnpEffStat.awk
       e4cfaee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5ce57f GENOMES: add Tetrahymena_thermophila.sh
       5f468e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4bd157a sampleSetup.pl patch
       a2dbacc samplesetup correction
       f4ba095 Trimmomatic.pm edited online with Bitbucket
       00d0af0 merging matser in dnaseq_report
       7ff833e rnaSeq.pl edited online with Bitbucket
       4556bda DNArepport update
       91d8b3c goseq.R edited online with Bitbucket
       e2148e7 goseq.R edited online with Bitbucket
       5a202eb DNAreport update
       a7a6d13 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cc7af8d TOOLs gtf2geneSize.awk protability
       05fdcce RNAseq update
       8b603f9 RNAseq update
       58cc9c6 RNAseq update
       dcf225d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b69a50b RNAseq correct htseqcount for stranded RNA
       f452010 RNAseq Strand-specificity correction
       7754be2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       bc5c5a5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a5450ec TOOLS rpkmSautrationadd more flexibilty on the file format
       4960f3e RNASEQ: add output directory creation in fpkm folder
       904595f General : Get back the change lost in commit 0007ea2
       84818ac Revert "RNASEQ: correct dependency issue"
       d45c5d7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d3b3bb1 Revert "RNASEQ: correct dependency issue"
       7aabc00 save before revert
       0007ea2 RNASEQ: correct dependency issue
       a3c5a79 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       eb6d023 RNASEQ : add deliverable + test exploratory
       f62067d Merged in rnaseq_lef (pull request #4)
       63a0163 RNASEQ - GqSeqUtils.pm update
       23099ae RNASEQ modify exploratory for the pull request #4
       80da056 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       af49cae RNAseq update
       c258e2c REDO lost commit: rnaSeq.pl use absolute path of design file
       7e28edd SubmitToCluster.pm edited online with Bitbucket
       79e8b60 add log.pm lib
       19289aa add logfile for cluster output file path
       ccf1fa8 RNASEQ: rnaseq.pl update
       5e8fcfc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2896458 RNASEQ: make design file path an absolute path to conflict if the working directory does match the relative path of the file given in argument
       363892b rnaSeq.pl edited online with Bitbucket
       6435ebe RNASEQ: rnaseq.pl && cufflinks.pm update
       3e54a64 SubmitToCluster.pm edited online with Bitbucket
       9ac95ff RNASEQ: rawCountMetrics update
       509b645 RNASEQ: rawCountMetrics update
       c484eed Merge branch 'logReport' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       6e725ad Merge branch 'master' into logReport
       a7faedc RNASEQ: mv dgeMetrics (step 10) as rawCountMetrics (step 8)
       04752e2 add log.pm lib
       3fcc283 add logfile for cluster output file path
       c527ebc RNAseq: change wingzip call in the metrics
       cbbd97a RNAseq update rnaseq.pl
       f4e3215 RNAseq: update rnaseq.pl
       2a8c65c RNAseq: update rnaseq.pl change matrixreadcount form step 10 to 7
       b43088a RNAseq update rnaseq.pl
       f672d80 RNAseq: update rnaseq.pl
       5535e50 TOOLS: add formatDenovoCombinedGTF.py & RNAseq: update cuffcompare
       a348a62 RNAseq: update rnaseq add cuffcompare
       cb485b8 rnaSeq.pl edited online with Bitbucket
       177aa34 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0c0f4a1 TOOL: Adding thepython tools folder and the getRefBinedGC.py script in it
       77b49e2 Merge branch 'mugqic_mathieu'
       1ebe60b Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       22ad9e5 update RNA
       36d735d RNASEQ: add correct cpu request for alignment
       e171861 correcting sv code small bugs
       553a3b0 finish breakdancer filter and  add pindel to dnaSeq
       ce1cd91 start adding breakdancer filter and pindel to dnaSeq
       59fa224 resolve Htseqconfilct
       0ad6350 merge with mugqic_mathieu
       c1a278a debug the code
       1498b14 Merge branch 'master' into mugqic_mathieu
       c2b1cc4 resolve lib/SubmitToCluster.pm conflict
       4df938a resolve HtseqCount conflict
       3268ff4 try to update local master to real master
       97b4432 everyhting except delivrable are done; debugging
       a0a5544  goseq done; 1st test on going
       eb5785e metrics: done; General module nomeclature conversion
       5efae81  metrics: updates; SAMtools: add viewFilter sub function ; tophatBowtie: bug correction and simplification
       2e48636  metrics updates
       be331d2 matrix done; DGE done
       4223cc7 fpkm done; DTE done ; DGE on going ; metrics on going
       08778dc correcting wiggle; fpkm and rawCount done ; DTE on going
       1ed2d7f wiggle done
       cd318c1 update metrics ini.file and few correction tophat/botwite
       052889d merging ok ; metrics started; new functions and changes in the picard lib
       2f6709a merging done ; metrics started
       66ea6a4 update to the master branch
       7166f74 Tophat-botwie done; merging started
       f9662b8 uppdating TopHat lib and rnaSeq.abacus.ini
       634a238 uppdating my branch
       9f81ed0 Starting bowtie/tophat lib
       43c3e7a Global dependency and RNAseq Triming and submit with working directory argument
       d6a6324 starting the rnaseq pipeline

  mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      96 commits

       ad59572 RNAseq update cuffdiff new dependency update
       278c9d2 RNASEQ: change variable nemae of bowtie fasta
       483e00d RNAseq update metrics stats
       8b02e13 RNAseq update trimming metrics
       2149362 RNAseq correct DESseq wrong dependency
       398a121 RNAseq upstae trimming merge stats; update edgeR
       b07481e rna update
       bed41b6 rna update
       4c45b56 remove tmp test for maty in dnaseq.pl
       45c8ba2 rna update
       ab2a149 R-tool: change saturation graph format
       7d052f1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d0bd499 RNASEQ: add format protability to goseq
       1cd3ae0 RNASEQ: add format protability to goseq
       2c7b901 RNASEQ: add format protability to goseq
       be3833a RNASEQ: add format protability to goseq
       3289aac RNASEQ: add format protability to goseq
       971938c RNASEQ: add format protability to goseq
       6cc1155 RNASEQ: add format protability to goseq
       9f84207 RNASEQ: rnaseq.pl update
       e31e700 RNASEQ: rnaseq.pl update
       fc1faad RNASEQ: add format protability to goseq
       18ffd9a RNASEQ: add format protability to goseq
       e725e01 RNASEQ: edger update
       7fbbb5a RNASEQ: goseq update
       c6dbccb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       73a94bf RNASEQ: goseq update
       06b2e4a RNASEQ: goseq update
       f93dba0 RNASEQ: allow non native GOseq analysis
       656c63c RNAseq: cuffdiff result fillter now include in the merge with fpkm
       02e9ef8 RNAseq: cuffdiff result fillter now include in the merge with fpkm
       bc787bd RNAseq update
       49ab91f RNAseq update
       e6515d3 change mugiqc_tools.install script to avoir facl conflict
       9676b39 update tools changes
       287a8e2 update tools changes
       e35a443 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       58ad96b RNASEQ: resolve splited metrics dependency
       4008a3d Rtools edger old index name  not suported by the actual edger version now use
       483781a RNAseq remove fpkm stats as they are done also in rna-seqc
       c47230d RNAseq correct bugs - see BFXDEV-20 for details
       8b614c0 RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence 2
       b5eee70 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       79804b1 RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence
       3b0ce98 RNASEQ: bowtie reference the new code should now be portable
       4b6bf8e RNASEQ: allow readstat on single library
       d079f59 RNASEQ: correct readstat output format
       a61e005 RNASEQ: metrics changes bug correction & Cufflinks correction
       28285f9 Metrics add  argument
       a1032f1 RNASEQ - correct typo in the mergebam step
       68b309e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       318f27b TOOLS remove type blatbestHit.awk
       30a3793 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5c8851c TOOLS adding blastbestHit.awk blatbestHit.awk
       e8481fd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       26152f7 RNASEQ: rnaSeq.abacus.ini correct rnaQC fasta variable name
       9bcdcf9 TOOLS: gtf2tmpMatrix.awk - adapt for line witout gene name
       57441b6 RNASEQ : adjust saturation plot title
       f05e6a7  RNASEQ : add saturation thread Number in the abacus .ini file
       1d072b2 Tools : make gtf2tmpMatrix.awk portable for various type of GTF file
       fbf1e9a tools : adding gtf2geneSize.awk and gtf2tmpMatrix.awk which were initially in tools but not anymore in the module
       2de9a90 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3ec125d RNASEQ : allow running without reference GTF .2
       a19c305 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a48116d RNASEQ : allow running without reference GTF
       4c9f803 pairedVariant - finish controlFREEC allow either bam or mpileup input everthing is specified in the ini file -  BFXDEV-7
       4a2439b R-tools change rpkmSaturation the count file need to be tab separeted
       93eda90 R-tools correct the polting bug of the rpkmSaturation - BFXDEV-3
       c67cdd0 rnaseq - replace my email by the  MAIL variable
       b1039a2 pairedVariant - add Control FREEC lib
       75351f1 update rnaseq.pl for single mode
       d155782 update rpkmSaturation to the last version
       6714873 correct formatGTFCufflinks.pl
       2220b26 Merge branch 'master' into mugqic_mathieu
       8aef01e Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       7ff4088 correct rnaseq ; remove special caharacter in subimti to cluster job ID varable
       39e7fd1 correct cuffdiff deNovo
       9db7cdf correct conflict btw mugqic_mathieu and master
       314eadc correct small bugs
       0307f25 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       17dadbe correcting pindel in lib/pindel.pm and in pairedVariant.pl form dnaseq pipeline
       3927743 remove python module call in htseqCount::readCountPortable
       ed644b5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fcc0e41 Merge branch 'master' into mugqic_mathieu
       d7d4284 correcting last bug in rnaseq
       4cb464e debug rnaSeq.pl
       7725d78 trimmomatic add java module loading
       1c9c334 submitToCluster add default  value =  if undef
       70afb54 rnaSeq debug
       2a78c37 resolve conflicts
       3c840f2 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       208d0cc debug
       cc57bb2 update
       8889bf1 rnaseq test debug; submitcluster modification; trimmomatic modifications
       bd4afe1 debugging RNAseq
       249ded6 debugging RNAseq

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      26 commits

       58770fe RNASEQ: change merge test
       73407b0 RNASEQ: change merge test
       54a772a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22f64a3 RNASEQ: add mkdir metrics in step 2
       be482de Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ee774f9 goseq.R: check and remove results that are not reported in the GO.db R package
       220e9fb lib GqSeqUtils report call change
       aa69e7e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       08a72c9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       20ce1c9 formatDenovoCombinedGTF.py update
       3183fd0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f6ce3df RNAseq update
       c833eee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d1e8dcc RNASEQ update unstrand wiggle bug correct
       49d7bbc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e7dbfbc RNASEQ correct stranded wiggle array assignation
       352f1ed add the java module call at Metrics::rnaseqQC
       7825d84 add the java module call before at each picard function
       aeb1c37 RNASEQ: add mamouth ini file
       2f858ec MODULES - add temporary download folder in several module install scripts
       e63d649 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       44263b2 rnaseq: resolve dependency pb
       9873fc2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89fbe54 RNAseq modify some output location for the reporting
       27bbc10 STATS: correct metrics:readstat for using output trimmomatic
       772a5e5 RNASEQDN: old mamouth ini file by the new one that fit changes

  Mathieu Bourgey <mbourgey@abacus2.(none)>      1 commits

       f36f0e2 rnaseq.pl debug; Metrics debug; Picard debug; SaMtools debug; SubmitToCluster debug; TophatBowtie debug

  Maxime Caron <mcaron@abacus1.(none)>      2 commits

       4ad23e1 test
       0b8c7fb test

  Pascale Marquis <pmarquis@abacus2.(none)>      5 commits

       168ea11 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       221a12e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b91d4da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0b27236 zorro.sh
       9f170dc Tetrahymena_thermophila.sh

