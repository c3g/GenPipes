38 tags, 9840 commits

HEAD        Tue Mar 14 17:48:02 2023 +0000        0 commits

4.4.1        Tue Mar 14 17:52:02 2023 +0000        10 commits

  ehenrion <edouard.henrion@mcgill.ca>      10 commits

       0d2e78fc2 Merged in release_4.4.1 (pull request #411)
       412c2e3d6 GenPipes - in prep for bug-fix release 4.4.1
       c27edf154 Resources - update mosdepth version to 0.3.3 in install script
       6a5222f2d EpiQC - correcting chromimpute_preprocess command : do not use "rm" and force symlink creation with the "-f" flag
       087fd9b3f Resources - adding / updating install scripts
       252b6ee50 RNASeq Light - Added mkdir to the kalliso command
       07f45a23f Resources - adding 'gget' install script and updating others
       ed8e359e6 Version bump to 4.4.1-beta
       20c470b78 Merged in release_4.4.0 (pull request #410)
       8cabdc1c6 Version bump to 4.4.0

4.4.0        Thu Mar 9 18:24:15 2023 +0000        156 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      24 commits

       afb8ae79f RNASeq - Fix arriba inputs
       18fb04116 RNASeq - fixing Arriba inputs and dependencies
       f4913ea1a RNASeq - Set the default build back to GRCh37
       658004c6e EpiQC - Fixing `chromimpute_convert` dependency with `chromimpute_preprocess`
       222c9d9f9 GenPipes - Resources : updated java and nextflow installation scripts
       be0334b98 GenPipes - Resources : some fixes in install_genome.sh
       77071b2ab GenPipes - Linting + multiqc dependency stengthening
       efeb90696 Minor update : just removing trailing comma
       039fef269 GenPipes - MultiQC : updated all the multiqc call with latest MultiQC version and removed the loading of Python module
       28943fc62 HiCSeq - corrected typo
       8095f0fec RNASeq - Adding 'cancer' and 'variants' protocols + removing 'cufflinks' protocol
       0ba1053e4 BFX - Linx : minor formating updates
       d33db1318 Tumor Pair - SV : correcting purple.py after rebase
       7e462a315 Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...
       5efad869d BFX - Conpair : corrected rm command
       d9b76059f Tumor Pair - fixed typo in filename
       7074e1008 Tumor Pair - fixes qualimap inputs to multiqc
       483a1c3b1 Tumor Pair - updating multiqc to 1.14
       26ac4f1c4 Tumor Pair - update inputs to multiqc
       aad32e3cb Tumor Pair - Strengthening the depedencies for the multiqc reports
       232129009 BFX - MultiQC : removed python module from command
       ada37b0c9 Tumor Pair - corrected multiqc dependencies
       ec0038784 Tumor PAir - Adding Compair and Purple to multiqc call
       fe29ad4c8 Version bump to 4.3.2

  ehenrion <edouard.henrion@mcgill.ca>      49 commits

       1c2947681 Merged in release_4.4.0 (pull request #409)
       bc806684a Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       468975cb8 updating README
       b0db9b96e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       67fe5bbe1 RNASEq Light - Fixed Kallisto step
       9e0141702 Resources - added HMMcopy package in defaults R libraries
       ffd223d03 RNASeq Light - Adding kallisto wrapper in bfx
       91bca69e8 RNASeq Light - revamping the call to Kallisto to have per gene counts instead of per readset
       d83071de7 Resources - Adding install script for hmmcopy-utils and ichorCNA
       79b74d038 Resources - Adding AGeNT install script and updating script for VarDict
       d2a8fef7b RNASeq Denovo Assembly - fixed trinity job dependencies
       9629f9a34 RNASeq Light - removing kallisto output folder before runnig kallisto to avoid mixed or bad results in case of pipeline restart
       eb7df8b18 EpiQC - fixed signal_to_noise dependencies
       48f936c11 DNASeq High Coverage - increased default walltime for gemini_annotations
       43121ea93 HiCSeq - fixing "homer_tag_directory" outputs
       4f7bfea47 Merge remote-tracking branch 'origin/dev' into release_4.4.0
       cf2828915 GenPipes - updating READMEs before release
       70d27f7c1 Resources - Adding some new install scripts
       699683cae Resources - updated version in some install scripts
       3d0c96d9c EpiQC - improved restart mechanisms
       dc3822fba RNASeq - Arriba : updated path of 'blacklist', 'known_fusions' and 'protein_domains' annotation files for GRCh37 in rnaseq.base.ini
       1eef454ce Merged dev into cit_fixes_mcj
       925fa891a Nanopore CoVSeq - minor linting update
       6f4b9a9bc EpiQC - fixed file path : distinction between files for job dependency and files for direct access
       d69696e65 GenPipes - resources : added ceratits_capitata.EGII-3.2 install script + version update in java.sh and picard.sh
       f9ca3c1c4 Tumor_pair - corrected java19 module
       07f2bfde6 PCGR report - added 'ls' command on the html report output to confirm job success
       991a65dbd Fixed mark_duplicate calls in all pipelines being from picard or gatk + Fixed pcgr report argument order in tumor_pair
       b6a29235c GenPipes - RNASeq : correctied calls to gatk_mark_dupliates + fixed arriba annotations for b38
       691ba5b64 GenPipes - Tumor Pair : updating cit.ini for GRCh37 default
       582055e4f GenPiopes - RNASeq : fixed picard/gatk_mark_duplicate
       362b72432 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       c5c436a03 RNASeq - Stringtie : updating resources for picrad mark duplicates
       0cec83cd8 RNASeq - Fixing sym link in star_align
       7faa99ae7 GenPipes - resources : update of GRCh38 install scripts
       9edf7b594 Resources - modules : adding pandas to the default python libraries
       bccb70316 Merged in tpsv (pull request #401)
       598e4c3f4 Resources - version updates for bcftools, htslib, samtools, python3 and multiqc
       fc5d59ff4 Adding 'core/__pycache/' to .gitignore
       317283ef0 EpiQC - minor change in chromimpute
       834752e20 Resources - Genome : updated 'create_bwa_index' with sym links to .fai and .dict
       aad898f31 Merged in tp_fix_multiqc (pull request #399)
       74f7b4512 Merged dev into tp_fix_multiqc
       6017946ad RNASEQ - Wiggle : corrected deeptools call
       1f444bca3 Merged dev into Gerardo_rnaseq_bigwig
       fd2b64af1 Merged in hotfix_dev (pull request #391)
       f64dc2b0c DNASEQ - ensemble_metasv dependency fix
       94efcdb75 GenPipes - Resources : Adding sv-prep installation script + updating gridss, gripss, linx, and purple scripts
       aa742e321 Merged in release_4.3.2 (pull request #389)

  Gerardo Zapata Abogado <gerardo_za_94@hotmail.com>      52 commits

       7b7706d4a Merged in Gerardo_rnaseq_bigwig (pull request #390)
       ba3ce5c71 rnaseq_denovo [in silico] small fix of the base.ini
       36d1ac84a rnaseq.py [wiggle] _32 addition of strand2
       841e7af7f rnaseq.py [wiggle] _30 addition of strand
       2003ed750 rnaseq.py [wiggle] _30 bamcoverage to wiggle
       e7f2e4916 rnaseq.py [wiggle] _29 added the addition of merged + forward + reverse - missing "s"
       354ffb02d rnaseq.py [wiggle] _29 added the addition of merged + forward + reverse
       99a979a4a readed deeptools.py
       eb2834dc0 removed Deeptool.py
       4bf5befae rnaseq.py [wiggle] _28 small changes to .ini and .py
       287500833 rnaseq.py [wiggle] _27 small changes to .ini and .py
       36a123e2b rnaseq.py [wiggle] _26 small changes to .ini and .py
       a09aabd21 rnaseq.py [wiggle] _24 small changes to .ini and .py
       24279a3bc rnaseq.py [wiggle] _24 small changes to .ini and .py
       12f3d6a72 rnaseq.py [wiggle] _22 small changes to .ini and .py
       39c127570 rnaseq.py [wiggle] _22 small changes to .ini and .py
       296ed2b98 rnaseq.py [wiggle] _21 small changes to .ini and .py
       4602dab69 rnaseq.py [wiggle] _20 changes to add option for for & rev.
       8a1490107 rnaseq.py [wiggle] _19 need to implement changes to add option for for & rev.
       b56bbc0c7 rnaseq.py [wiggle] _18 location of ifelse loop
       3834821fa rnaseq.py [wiggle] _17 set up option to run separate .bw??
       ecf6046c8 rnaseq.py [wiggle] Deeptools _16 add forward and reverse options
       7265f011a rnaseq.py [wiggle] Deeptools & Ini changes _15 small changes
       e2f2d5582 rnaseq.py [wiggle] Deeptools - _14 string errors
       b46b26b12 rnaseq.py [raw_countMatrix] -  _12 for loop removal
       56ee87422 rnaseq.py [raw_countMatrix] -  _11 bigwigs .zip
       d53d4f4ba rnaseq.py [wiggle] - Deeptools _11 "missing comma" :|
       112882f70 rnaseq.py [wiggle] - Deeptools _11 "output_location error" + "archive"
       6cdd381a0 rnaseq.py [wiggle] - Deeptools _10 "name" "samplename" location
       4877bab61 rnaseq.py [wiggle] - Deeptools _9 "name" "samplename" location
       340669c27 rnaseq.py [wiggle] - Deeptools _8 missing comma
       35c25218f rnaseq.py [wiggle] - Deeptools _t "name" "samplename"
       218a97899 rnaseq.py [wiggle] - Deeptools _6 missing comma
       a361eb147 rnaseq.py [wiggle] - Deeptools _5 "name" "samples"
       89a3e4788 rnaseq.py - [wiggle] - concate_jobs
       ca916a9d9 rnaseq.py - fix Deeptools.bamCoverage()
       26b586954 Deeptools.py - fix bamCoverage
       97733826e rnaseq.py [wiggle] - Deeptools
       08e56f631 rnaseq.py [wiggle] - fix "tracks_directory"_4
       f5b8293b9 rnaseq.py [wiggle] - fix "tracks_directory"_3
       c2b03d99c rnaseq.py [wiggle] - fix "tracks_directory"_2
       91073b16c rnaseq.py [wiggle] - fix "tracks_directory"
       ed2f1c48c ranseq.py - added new wiggle step, removed previous
       0ae71bb29 rnaseq.py Modify [wiggle] ini
       bb9b9120a Added - Deeptols.py
       fe825b146 rnaseq.py - readcounts - ".csv" to ".tsv" - ALL
       0bb883c02 rnaseq.py - rawCountMatrix - ".csv" to ".tsv" - ALL
       1b247225c rnaseq.py - readCounts - ".csv" to ".tsv" 2
       90bb9e292 rnaseq.py - readCounts - ".csv" to ".tsv"
       355a51915 rnaseq.py - rawCountMatrix - ".csv" to ".tsv"
       0ed06e096 test commit Gerardo
       2b424a138 test commit

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      4 commits

       2121b4f25 methylseq : do not require methylation_protocol and mapping_implementation parameters to be set for hybrid protocol
       791aeb882 GenPipes methylseq : make single-pass dragen align option compatible with bismark in hybrid protocol
       b0507cf39 GenPipes - methylseq : fix cp from dragen
       2d69b61fc issue caused by qiime_catenate loading two python modules

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      23 commits

       33a2685f2 Merged in cit_fixes_mcj (pull request #408)
       9e6f6e809 Rnaseq cancer : update genome version used for cpsr and pcgr
       a48b0ed58 Merged dev into cit_fixes_mcj
       a95676086 Merged in cit_fixes_mcj (pull request #407)
       4cab1006f gatk_indel_realigner : set -fixMisencodedQuals flag only when QualityOffsets are 64 (solexa/illumina1.5+)
       0e459caab skewer trimming : adjust trim settings based on quality offset in readset file
       b28269b72 Merge branch 'cit_fixes_mcj' of bitbucket.org:mugqic/genpipes into cit_fixes_mcj merge cit updates with local changes
       575d3052e skewer trimming : adjust trimming settings based on quality offset provided in readset files
       88128229b Merged in EpiBrain_methylseq (pull request #406)
       2ab0e92b7 Merged in cit_fixes_mcj (pull request #405)
       1b9d47c7e Merged in EpiBrain_methylseq (pull request #404)
       ab58ff992 GenPipes - rnaseq : adding cd command to run_arriba
       5801352eb removing outfile options from arriba
       349822756 Merged in cit_fixes_mcj (pull request #403)
       25f08c9b5 GenPipes arriba : added specific paths for arriba output files to run_arriba command
       33515ecd3 GenPipes - rnaseq : removed -fixMisencodedQuals from ini for gatk_indel_realigner
       1ffffc2e7 GenPipes - rnaseq : fixed file extension in ini for run_arriba
       83f0d22ff GenPipes - epiqc : add mkdir to chromimpute_convert
       2c6d4c719 Merged in ampliconseq_abacus_mcj (pull request #398)
       ce9f2755b Merged in rnaseq_star_1pass (pull request #394)
       fae754eea Merged dev into rnaseq_star_1pass
       cbfc7659b Merged in tumor_pair_mcj (pull request #393)
       c9343cf6e Merged dev into rnaseq_star_1pass

  Mareike Janiak <mcj43@narval1.narval.calcul.quebec>      1 commits

       0522542f4 ability to skip 2-pass star alignment and only do 1 pass

  Mareike Janiak <mcj43@narval2.narval.calcul.quebec>      1 commits

       cd4a4556c gatk_indel_realigner : added two cd commands

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      1 commits

       388584362 test paul

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       638b81a54 Merged in rebase_rnaseq_variant (pull request #370)

4.3.2        Thu Dec 8 17:21:23 2022 +0000        132 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      54 commits

       ca99865e0 Merge remote-tracking branch 'origin/dev' into release_4.3.2
       118e25a36 GenPipes - Updating READMEs for release
       fdb289a0d Tumor Pair SV - properly setting the svaba indels vcf file in the config
       f568fe4df Tumor Pair - rearrangement of the ini for cit
       475e2b944 Tumor Pair SV - fixing metasv.ensemble call + fixed scones step for GRCh38
       eb7112e07 Tumor Pair SV - fixing step order
       16cd02b99 minor - replacing 'CentOS6' by 'root'
       53e3e05fc Some more fixes induced by cit
       29e469727 RNASeq denovo asembly - fixed filter_annotated_components_exploratory
       c4d242b27 RNASeq denovo assembly - fixed differential expression steps
       77f4ca476 AmpliconSeq - enuring the use of Python2 for all the qiime commands
       143bc39ac Tumor Pair Ensemble - fixed strelka dependencies to manta_indels
       ab5692610 RNASeq - updated mugqic_tols to later dev version
       0aa6a1b23 DNASeq SV - updating cit.ini
       5c6358fc5 DNASeq - corrected typo in gatk.py
       ba6439b8b RNASeq denov assembly - updated mugqic_tools to later dev version
       73a4da529 AmpliconSeq - setting mugqic_tools to later dev version
       7227a9409 AmpliconSeq - some more formating
       93f057f7d RNASeq - correction of rnaseqc and htseq_count after cit
       2fdb99b52 AmpliconSeq - correctig job names
       152103b7a DNASeq SV - fix breaseq2 call
       8c7383221 RNASeq - test batch effect corecction in DESeq2
       06a8eeb7b RNASeq - fixed typo
       b367867ff RNASeq - fix after cit
       74bc42022 HiCSeq - fix quality_scores afeter cit
       6ef0a2190 DNASeq SV - fix svtyper after cit
       6f4cb2deb pipelines/hicseq/hicseq.base.ini
       8264c71f4 DNASeq High Coverage - fix preprocess_vcf after cit
       7174d35fe Ampliconseq - fix flash step after cit
       a008986e5 Tumor Pair - Fixing strlka dependencies to manta_sv
       d60941184 RNASEQ + RNASeq Denovo : fixed batch corretion for RNASeq and added it to RNASeq denovo
       40cf4033b Corrections after cit + formating
       185ed8310 GenPipes - Fixes after cit on Béluga
       2bebb47ef Fixes after cit + code standardization
       c9b9e2bb8 DNASeq - SV : updating metasv module to 0.5.5
       b7d0ce4d2 DNASeq - SV : re-using delly VCF in meta_sv ensemble analysis, depending on meta_sv version
       3004136df DNASeq - SV : corrected manta_sv call by removing delly_vcf
       7c0444012  DNASeq - correcting typo
       efa9aa387 DNASeq - correcting typo
       341de84dd DNASeq - generalize use of self.output_dirs
       6775aad45 DNASeq - specify self.output_dir when setting files and folders path
       0284d9aae DNAseq - SV : correcting ini sectino definition
       313dc0028 BFX - correcting typo in `cpsr.py` and `htslib.py`
       224f3eaf1 DNASeq - bug fixes for SV protocol
       48c925f95 DNASeq High Coverage - update sambamba_merge_sam_files in config
       ed1b33a66 Tumor Pair - format standardization
       3c82ce6f6 Tumor Pair - minor fix
       b21669d52 Tumor Pair - bfx/pcgr.py updated to fit with new version and keep backward conpatibility.
       e830d5825 Tumor Pair - Report : updated with new version of PCGR (1.0.3) for both 'pcgr' and 'cpsr' reports
       c024bd6bf DNASeq - corrected merge_sam_files + extract unmapped steps
       8279cfa40 Core - Scheduler : corredcted PBS mem-per-cpu settings when --force_mem_per_cpu
       a34fda403 Core - Scheduler : corrected `memory` call in `cpu` from dev
       4a045d5a4 Core - Scheduler : corrected `memory` call in `cpu`
       80afe6dfd Version bump to 4.3.1

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      45 commits

       88e0ba82e CovSeq - Improving restart mechanisms
       a1ef93526 ChIP-Seq - improved restart mechanisms
       962c90d93 AmpliconSeq - Improving the restart mechanisms
       aea76ebf9 AmpliconSeq - Improving restart/resume mechanisms
       4e401e9fa Tumor Pair SV - correcting svaba ref configuration
       d63c56d58 Tunmor Pair SV - CNV detection now done with cnvkit : removing SCoNes, not as well maintained as cnvkit.
       952693d8f Resources - corrected gridss install script
       3031bd63a DNASeq - returning to compressed filtered vcf for cnvkit
       705ace9ac BFX - cnvkit : updated inputs in scatter and segment
       025165cc6 CORE - removing useless log.debug() in core/pipeline.py
       29e81cc36 DNASeq SV - testing with no compression during cnvkit_batch.vcf_flt
       cefc1dbb4 Tunmor Pair - Fix cit ini
       b411c0c92 RNASeq light - correcting R version for kalisto differential expression
       b5e1d2ebe Tumor Pair - fixing cit ini
       5b440da92 Resources - put back DEV instead of wrongly modified GENFS
       ca5745802 Modules - adding install script for tools needed in Tumor Pair SV
       7bc7fc395 DNASeq SV - corrected bcftools filtering before cnvkit + corrected metasv_ensemble command
       b67c4aaa3 DNASeq SV - improved bcftools options for vcf filtering before cnvkit_batch
       8536530ce DNASeq SV - corrected filtered vcf outputed by cnvkit_batch.vcf_flt
       ae84cc421 DNASeq SV - improved filetering in cnvkit_batch.vcf_flt
       56f73fd1b DNASeq + Tumor Pair - fixed breakseq2 wrapper
       4dd168e04 Core - Pipeline : corrected log.error to log.debug
       f8eaf1dd5 Resources - updated mugqic_tools instalation script, adding java-tools in the module file
       be5efebd9 RNASeq denovo assembly - corrected outputs of edger and deseq2 jobs
       f6a7bb8d9 INI - updating mugqic_tools to latest production version for AmpliconSeq, RNASeq and RNASeq denovo assembly pipelines
       2b12dea55 DNASeq - correcting input dependencies for manta_sv jobs for better restart/resume of the pipeline
       c9493affc RNASeq - refined ballgown job outputs to improved pipeline resume/restart mechanisms
       6300bafff RNASeq - imporoved restart/resume mechanisms
       a5df32f91 DNASeq + Tumor Pair - imporoved restart/resume mechanisms
       1c7a63f4c AmpliconSeq - fixes  on python modules and restart mechanisms
       9b21440d8 RNASeq denovo assembly - correcting arguments in calls of py_parseMergeCsv
       fc88f6a24 DNASeq - bug fix
       1417bd8e9 Tools - correcting typo in bfx/tools.py
       096206be7 Config - corrected pipelines/common_ini/Homo_sapiens.GRCh38.ini
       4deb294b7 GenPipes - Config : reorganisation of the ini files contained within the repository
       a7883c1a5 GenPipes - config : use common_ini for GRCh38 specific config
       e784b3595 GenPipes - DNASeq & Tumor Pair : cleaning ini files
       723b201ae DNASeq SV - cleaning resources/genomes/config/Homo_sapiens.GRCh38.ini a bit...
       0b1e36a3a Resources - Adding some genome installation scripts
       b316c215c RNASeq denovo assembly - Corre ting calls to py_parseMergeCsv
       cadad2206 AmpliconSeq - Correcting python call in qiime_otu_table
       9f3f3e75d Resources - adding some more module installation scripts
       20c0951e3 Install script - updates in module install scripts
       3341cbc88 Tumor Pair - Improving the restart mecanisms
       df5aa276c GenPipes - DEV : updating CHANGELOG and setting beta version

  ehenrion <edouard.henrion@mcgill.ca>      17 commits

       ed7c18351 Merged in release_4.3.2 (pull request #388)
       23f764a4b GenPipes - Resources : updated mugqic_tools install script with latest version
       15091a007 DNASeq SV - updated pair_diretory path with non-absolute path
       b18b9d3af EPIQC - worked on file path to avoid depencencies errors and missing files
       5678259df HicSeq - stop uding mugqic_dev
       552c0fba4 RNASeq Light - fixing kallisto_count_matrix dependencies
       ee70b175c Tumor Pair SV - fix delly sv annotations dependencies
       71f958976 Tumor Pair SV - fixing SVannot call
       3e68d78ac GenPipes - more formatting and linting
       90b953538 GenPipes - code cleaning and formating
       2dcd2e6e3 CoVSeQ - using report directory to output freebayes and ivar reports
       ef89a8f43 GenPipes - improved averall restart/resume mechanisims with relative path
       10727fed4 HicSeq - working on restart/resume mechanisms + doing some reformating
       ce5e36f27 Merged in tp_hotfix (pull request #385)
       03d0ab285 DNASeq - gnomad : corrected link to vcf
       118309c84 Merged in release_4.3.1 (pull request #383)
       baeee47d9 Merged in release_4.3.1 (pull request #382)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      16 commits

       0148c73f3 General - Linting
       dacf35a07 Merged in HotFix_dev (pull request #387)
       b9419c5b8 tumour_pair - artifacts_metrics - Fixing output file name for GenPipes to be able to skip step
       e799a375b covseq - sambamba_filtering - Adding SE options in ini file
       27b4f5506 covseq - sambamba_filtering - Debug for SE mode
       cafd9749c rnaseq_light - base ini - Updating mugqic_R_packages
       7ee0df69d Merged in bwa_sambamba_splitting (pull request #386)
       0182fc0f0 dnaseq_high_coverage - bwa_mem_sambamba_sort_sam - Adding bwa_mem_sambamba_sort_sam in ini
       5c2dccb9c dnaseq/tumour_pair - sambamba_sort - Minor improvement
       8b57df09e dnaseq/tumour_pair - sambamba_sort - deleting bai if existing and chmod the output bai
       1050404cc dnaseq/tumour_pair - sambamba_sort - Removing index from sorting as the bai is generated by sorting
       506605f43 Merged in bwa_sambamba_splitting (pull request #384)
       656678db0 dnaseq/tumour_pair - base ini - Cleaning commented out code
       dd4ca6901 dnaseq/tumour_pair - sambamba_sort_index - Fixing index input
       bef6394ea dnaseq/tumour_pair - alignment - Switching in protocols
       2fceb1040 dnaseq/tumour_pair - alignment - First commit

4.3.1        Tue Oct 4 14:11:09 2022 +0000        85 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      7 commits

       a647b011d Merge remote-tracking branch 'origin/dev' into release_4.3.1
       759bbfb9f GenPipes : updating READMEs for release
       241077185 GenPipes - PBS Scheduler : stop  using -mem/-pmem parameters anymore in PBS job submission commands
       57667889e Version bump to 4.3.0 - for real
       937a5ca4b Version bump to 4.3.0
       69d74b05f Merge remote-tracking branch 'origin/dev' into release_4.3.0.1
       ee071fd8d Updating READMEs

  ehenrion <edouard.henrion@mcgill.ca>      9 commits

       f0feb9adf Merged in release_4.3.1 (pull request #381)
       94d7a6cbb GenPipes - ChIPSeq : Updated mugqic_tools to 2.10.10 in chipseq.base.ini to prevent crash at differential_expression step
       6da854f3c GenPipes - ChIPSeq : updated mugqic_tools to 2.10.9 in chipseq.base.ini
       db7a44852 tumor_pair.extras.ini edited online with Bitbucket
       285568f96 GenPipes - README : replacing former "CentOS6" by "root" in cvmfs path and variables
       4d824c56f Merged in ehenrion/tumor_pairpy-edited-online-with-bitbucke-1655753579034 (pull request #369)
       f322ebd79 tumor_pair.py edited online with Bitbucket
       e5980e368 Merged in release_4.3.0.1 (pull request #367)
       110f8ad7e Merged in release_4.3.0.1 (pull request #366)

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      1 commits

       c4a59f503 Changed the core pipeline script to modify the trace config ini file to include the timestamp. Additionally, the full command is now added to the ini header as a comment for future reference.

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      7 commits

       468bbfec3 Merged in config_timestamp (pull request #374)
       ed5738de5 Merged in covseq_nanopore_GPUqueue_bugfix (pull request #373)
       e1a911218 nanopore_covseq.base.ini edited to re-add the GPU_QUEUE value. Without it the pipeline crashes on Abacus.
       bd66f86dd Merged in MetaSV-Delly-bug (pull request #372)
       693ee1d97 Removed Delly input file requirement from dependency list in metasv.py
       bf99d700b Removed delly referneces from the metasv step in the tumor pair pipeline.
       4183eb1f4 MetaSV has a delly argument that should not be there since it does not support Delly (see https://github.com/bioinform/metasv/issues/110). I removed the delly parameter.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      30 commits

       e2f0de69c Merged in HotFix_dev (pull request #378)
       56b2fcc5e rnaseq_light - base ini - Improving syntax
       511f59a26 methylseq - General - updating way to mkdir via bfx
       dec5d5069 methylseq - filter_snp_cpg - Debug
       908b158f3 methylseq - filter_snp_cpg - Allowing to change cpg threshold for cit minimal
       1b7086496 methylseq - methylkit - Changing the way Rscript is called
       5d2395f51 rnaseq - gq_seq_utils_exploratory_analysis_rnaseq - Upgrading to latest R_packages + fixes
       210a47372 rnaseq - gq_seq_utils_exploratory_analysis_rnaseq - Upgrading to latest R_packages
       2fe1ac6ed chipseq - differential_binding - Adding parameter contrastnb for cit atacseq only (set it to "cit" for cit minimal)
       794a06639 hicseq - homer_tag_directory - Adding chrList from ini file
       1619b24b1 hicseq - quasar_qc - Allowing more than 1 "_" in chr name
       9117fe423 hicseq - General - Upgrading mugqic_tools and trimmomatic
       5269243a6 hicseq - identify_TADs_RobusTAD - Adding min and max window parameters
       cf6418827 hicseq - multiqc - Upgrading to latest multiqc version
       b1dead1d1 hicseq - hicup - Moving remove before hicup conf
       0a72d80f2 hicseq - hicup - remove recursive
       28b162909 hicseq - hicup - Re-adding removing of results otherwise the step is not able to restart
       d837f8dae hicseq - General - Improving the way genome assembly name is handled
       6e2d1aa0a hicseq - hic - Improving juicer call
       9e539f6b2 hicseq - hicup - Debug
       0c1b223c2 hicseq - hicrep - Code cleaning
       3018b34ff hicseq - hicup - No longer creating folder and removing things in the step
       3ccf80788 hicseq - hicup - Debug
       42245c79e hicseq - General - Improving hicup job and removing path hardcoded in ini
       2424c6aa2 hicseq - General - Upgrading default genome to hg38
       c58d79530 epiqc - General - Improving code
       ce9973551 Merged in dragen_update1 (pull request #365)
       0c39884cb rnaseq - General - Updating doc.
       307995288 General - General - Formalizing README doc.
       e6776b455 epiqc - General - Improving code

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       eb2067d41 Merged in force_mem_per_cpu (pull request #380)
       9822d934a Merged in fail_on_log_pattern (pull request #376)
       e1eeb1773 Merged in sample_readset (pull request #338)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       c3e418bf7 support G and M only
       13c7c51f9 Remove explicit request for memory in pbs call when --force_mem_per_cpu option is used
       a74d9547c add a fail_on_pattern ini option
       4b796dea7 mv sample_tumor_pairs.py to core
       92049a777 move sample readset and design to core from bfx

  pubudumanoj <pubudumanoj@gmail.com>      2 commits

       d2990b248 modified inputinfofile and chr_sizes path
       847d47659 modified inputinfofile and chr_sizes path

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      18 commits

       c6bcd4239 remove comments in tumor pair sv
       7a85f1d85 fix more bugs in tumor pair sv protocol
       2f76d8fc4 fixed more bugs
       ca6c2826e fixed major issues in sv protocol
       3908f5e0f modified the methylseq/cit.ini
       663b11f4b fixed a regression
       3166cdf1b chnaged resources in cit.ini
       d5aa23879 changed mugqic_tools version to 2.10.9
       ba2956457 fixing issues
       a93935fcf added methylkit differential binding back to the analysis
       830b31751 modified base.ini
       6e7dd92d3 Merge branch 'epiqc_code_improvement' of bitbucket.org:mugqic/genpipes into epiqc_code_improvement
       0c1426bc1 Corrected the inputinfo file paths in ini
       f787ac781 expanded env_var in inputinfo file path and chr_sizes path
       4e7e8e7a1 Corrected the inputinfo file paths in ini
       9684c5e2f fixed some issues in epiQC md file
       b772f84a6 Added multiline comments to methylseq pipeline
       191ed379c expanded env_var in inputinfo file path and chr_sizes path

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       12be5b9e1 Merged in tumor_pair_hotfix_metsv (pull request #375)
       ceadb4e8b Merged in methylseq_cit_update (pull request #371)
       f2678bed9 Merged in epiqc_code_improvement (pull request #368)

4.3.0        Tue Jun 14 17:50:22 2022 +0000        118 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      21 commits

       3f33a35a6 Version bump to 4.3.0
       76b0b6264 Merge remote-tracking branch 'origin/dev' into release_4.3.0
       ec90b7abe GenPipes - prep for new release
       863677573 Updating CHANGELOG from master after release
       1e06b04a9 Version bump to 4.2.1
       4e9074168 GenPipes - DNASeq : fixing dependencies in recalibratino step for unmapped_reads
       46042e1cc GenPipes - remove useless sections in ini
       40d6fe2f4 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7dff261f5 GenPipes - move 'sambamba_merge_sam_files' into common.py + renamed 'sambamba_merge_sam_files' to 'sambamba_merge_sam_extract_unmapped' in DNASeq
       34ea28a1d Cleaning bfx/sambamba from samtools dependencies
       adf5e7ff7 GenPipes - DNASeq : revisited recal jobs to avoid useless restart of the step
       ed49153b2 GenPipes - DNASeq : Fixing dependencies for recal .bam.bai
       b7ef560e3 GenPipes - DNASeq : Recreating of the index of the recal bam file after merging of the unmapped reads
       373fb134e GenPipes - DNASeq : correcting the job which merges the unmapped reads to the recalibrated fastqs
       d30127dda GenPipes - DNASeq - Fix job queueing for merge unmapped to recalibrated
       1fc52437f fix typo
       afffee2b2 GenPipes - DNASeq : corrected unmapped_reads job queueing
       379dfdfb4 GenPipes - DNASeq : removing useless part of code...
       4b897f01c GenPipes - DNASeq : extracting unmapped reads from the merged bam, then re-instering them into the recalibrated bam
       460b6a3c8 Version bump to 4.2.0
       5cc640ceb Version bump to 4.1.4-beta

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       5690df407 GenPipes - C3G Software Stack : showing only the 3 top versions of the softwares and adding '...' when necessary to indicate more older versions are also available
       34b06e68d GenPipes - C3G Software Stack : reverse-sorting the software versions when combining the metadata JSONs
       416fe4a48 GenPipes - C3G Software Stack : improved sorting of JSON output
       864465e62 GenPipes - C3G Software Stack : sort JSON output
       d3f03f6db GenPipes - C3G Software Stack : parametrized the redirection of the output to a file or to stdout
       701298b12 GenPipes - C3G Software Stack : removed the sending of the JSON
       85b700910 GenPipes - C3G Software Stack : corrected typo
       833d3caf3 GenPipes - C3G Software Stack : corrected message when file is sent
       3303a59fc GenPipes - C3G Software Stack : Added sending of the combined JSON to url + reformating the module helpers
       0b8b7d59c GenPipes - module helper : updating the pythonSearcher class

  ehenrion <edouard.henrion@mcgill.ca>      9 commits

       2e0a7ecc9 Merged in release_4.3.0 (pull request #364)
       390169f30 Merged in release_4.3.0 (pull request #363)
       47523b322 Merged in release_4.2.1 (pull request #359)
       5f690ffe9 Merged in soft_jsondb_gsoc2020_eh (pull request #181)
       ba84d8c17 dnaseq.base.ini edited online with Bitbucket
       943e103cd Merged in unmapped_reads (pull request #355)
       ef817ea69 fixed typo
       ecdcb0ae3 dnaseq.base.ini edited online with Bitbucket
       9f80a8c9f Merged in release_4.2.0 (pull request #354)

  Moonshroom <yatharthrai16@ducic.ac.in>      4 commits

       a51bba0f1 Some documentation edits, README updates, fixed bad variable naming
       8d223caba Migration to argparse
       c6bde0f2f Added -h/--help flag
       810321f07 1. Pushing pre-final code 2. Added CLI args 3. Added CLI Documentation

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      20 commits

       cd2845a88 Merged in HotFix_dev (pull request #357)
       a5caedfd6 General - General - Fixing UnboundLocalError: local variable 'job_name_prefix' referenced before assignment
       100567c62 Merged in HotFix_dev (pull request #356)
       a35809f59 General - General - Fixing typo
       edc01ae2d General - General - Fixing typo
       20801c451 General - General - Adding step_wrapper to slurm and batch mode (fixing regression)
       e7869220f General - general - Standardizing inis + fixing methylseq.base.ini
       343ec56db Merged in HotFix_dev (pull request #350)
       401be7909 ampliconseq - General - Upgrading to latest release of mugqic_tools
       a52c6d626 covseq - General - Upgrading to latest release of covseq_tools
       92e0a5f74 ampliconseq - flash - Debug previous flash commit
       d73c35dfb ampliconseq - ampliconLengthParser - Fixing head with pipefail
       c74fe1ca5 ampliconseq - flash - Fixing None appearing in cmd
       720962b8a covseq - qualimap - Fixing qualimap parameter for threads
       92081f8e4 covseq - prepare_report - Fixing pipefail issue for ncovtools: we always want ncvotools error to be ignored and continue the job
       e40e7cc16 covseq - prepare_report - Fixing pipefail issue when grep doesn't find anything
       2c6b46c58 covseq - prepare_table - Adding threads param
       5387d2981 covseq - prepare_table - Adding samtools as module required
       b8cfc5eaa methylseq - metrics - Changing ini path of file
       f20a986a4 General - general - Standardizing inis + fixing methylseq.base.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       cf3d7277c Merged in release_4.2.1 (pull request #358)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       c7b3aa6dc Merge remote-tracking branch 'origin/dev' into release_4.2.1
       b09200275 Update readme and version for release
       5469a9037 remove -nt completely from indel realigner
       ccf852903 abacus ini typo
       720fd57f3 typo in tabix option, and no threads in gatk_indel_realigner
       68183a842 always force overwrite on tabix
       ef425e7e2 tweak and fix base ini

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      39 commits

       79ef5180f addressed comments from Paul
       2011bfb5c removed methylseqraw class
       a726b9778 Fix issues in input and output dependencies
       66e55e3bf resolved differences in dev branch files
       9c8dd415d Modified dragen.ini file
       c5acaceb5 fixed documentation typo
       5949d7bc5 fixed issues in output file tracking added documentation
       fea7a92cd fixed issue in dragen hybrid protocol in single-pass
       8fca6d1c6 added dragen_bedgraph step
       454abe680 added split_dragen_methylation_report
       313ce780c fixed issues in dragen_methylation_call step working pipeline
       e2a952e2d trying to fix dependencies
       0c481543c added missing job to the dragen_aling
       15dda5161 fixed a bug
       d3e987499 Add changes in common and methylseq
       c85318991 Implement jsonator disabling at job level
       4d64d7340 Debug dragen protocol
       d1c7bf2a5 update dragen protocol
       8ed0ab89b Implement Modified dependency of dragen command to allo mv commands Added dragen function
       2ce9490bc Implement Modified dependency of dragen command to allo mv commands Added dragen function
       dfa349d63 fixing issues with the done file tracking system
       861201122 Trying to fix dependency issue in dragen
       61945dc12 resolved merge conflicts
       f22f0f663 resolved merge conflicts
       6664a8571 added dragen mark_duplication
       427d7a02c fixed a bug
       33caf9049 update dragen protocol
       3cb1965f4 added missing job to the dragen_aling
       e0e21a72a Resolved incorrect tool link
       93a084b3a fixed a bug
       30ef76c40 fixed a bug
       fa534de9e Add changes in common and methylseq
       e8833f502 Removing debug outputdir log info
       8485ae7dd Implement jsonator disabling at job level
       98cc9725e Debug dragen protocol
       5ec4bfc87 update dragen protocol
       c32117b8d pipelines/methylseq/methylseq.dragen.ini
       5e71366ec Implement Modified dependency of dragen command to allo mv commands Added dragen function
       179cabf84 Implement Modified dependency of dragen command to allo mv commands Modify dependency of bash commands Added bfx/dragen.py

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       486907ab6 modified readme file
       be1755ab7 modified readme file with dragen update
       d82bab012 added comments and fix hybrid mode single pass directional-complement
       865a65f5e fixed ihec metric by adding a job to create an empty metric file for estimated_library_size
       cd94f5562 changed dragen structure and added a new protocol

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       42112207e Merged in dragen_update2 (pull request #362)
       7faffc517 Merged in dragen_update_version2 (pull request #346)

4.2.0        Wed Jun 1 20:01:30 2022 +0000        230 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      16 commits

       695798b17 Generating READMEs for 4.2.0
       9947d057e Version bump to 4.1.3
       d7f194248 Merge remote-tracking branch 'origin/dev' into release_4.1.3
       e12f9ed23 in prep for a release
       e135fce3c GenPipes - Tumor Pair : corrected location of the coverage_bed in strelka jobs + prepend strelka job with 'rm' command to ensure folder is clean
       b7ed6142e GenPipes - Tumor Pair : some code cleaning before correcting strelka errors
       1d34130e9 GenPipes - DNASseq : updating base.ini
       20453e2b3 GenPipes - DNASEQ : fixing base.ini for SV protocol
       fe3ccfcd8 GenPipes - Scheduler : Re-using the default python, from ini, to call job2json.py
       61218acd8 GenPipes - DNASeq SV : fix sambamba_merge_realigned
       fa7cc656b GenPipes - DNASeq SV : using Python2 for breakseq2
       b790c9135 GenPipes - BFX : CNVkit now uses the muqgic/CNVkit module, instead of a pyhton module which would
       9db4af258 GenPipes - DNASeq SV : set use of Python2 for MantaSV in base.ini
       97513bd2b GenPipes - DNASeq SV : fixing resource dependencies
       5685c7b89 GenPipes - DNASeq SV : fixed delly call
       1ef5a9ecf Version bump to 4.1.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      4 commits

       da1d31b48 Resources - install scripts : update to mugqic_tools 2.10.4 + install_modules updates for improved installation in dev space
       5616bb75a GenPipes - Resources : ensuring the use of rpath when patching binaries
       b15769cb6 GenPipes - Resources : adding some new installation scripts
       802c02d2e GenPipes - Resources : Updating some install scripts with newer version + adeind bioinfokit and lofreq install scripts

  Édouard Henrion <henrione@narval1.narval.calcul.quebec>      3 commits

       13b4a69ab GenPipes - Tumor Pair : setting Python2 for Strelka2 in the config file
       69201a362 GenPipes - Tumor Pair : set the output_dirs global variable to better handle outputs, thus dependencies
       a9bf30ce7 GenPipes - Resourcess : adding pycoqc install script + updating other scripts

  Édouard Henrion <henrione@narval2.narval.calcul.quebec>      4 commits

       70c9b89af GenPipes - Tumor Pair : fixed dependency issues by using relative path everywhere in the pipeline, no more absolute path. Also setting 'samples' for all the jobs.
       1156e33a8 GenPipes - Tumor Pair : updating mugqic_tools to latest version 2.10.2
       2a541dece GenPipes - BFX : updated GATK wrappers to remove print_reads outputs before running
       576724b0a GenPipes - Tumor Pair : fixing Strelka & Manta SV step

  Édouard Henrion <henrione@narval3.narval.calcul.quebec>      1 commits

       bab3a852f GenPipes - Tumor Pair : Fixing conpair with specific call to Python2

  ehenrion <edouard.henrion@mcgill.ca>      16 commits

       927e40e45 Merged in release_4.2.0 (pull request #353)
       a546428e1 Merged in release_4.1.3 (pull request #352)
       f1e6716e0 Merged in release_4.1.3 (pull request #351)
       43e87b722 dnaseq.base.ini edited online with Bitbucket
       cdcebcb1c Merged in tumor_pair_hotfix_dev (pull request #339)
       020505790 chipseq.base.ini edited online with Bitbucket
       f52b3e24e Merged in tumorpair_hotfix_strelka (pull request #331)
       4db4e81b7 Merged in dnaseqsv_fix_eh (pull request #329)
       66eca370d Merged in fix_dnaseqsv_eh (pull request #320)
       e634cddfc Merged in fix_delly_call_dnaseqsv_eh (pull request #317)
       2c20036d7 GenPipes - Scheduler : stop using python from ini, use default mugqic python3 instead to ensure job2json is working in al cases
       2ee5d49ee GenPipes - COVSEQ : corrected covseq.graham.ini for dna_sample_qualimap
       0b5b8967c GenPipes - COVSEQ : corrected covseq.cedar.ini for dna_sample_qualimap
       41b755e33 GenPipes - COVSEQ : corrected covseq.beluga.ini for dna_sample_qualimap
       2364a6c8b GenPipes - Readset : corrected "writer" call in checkDuplicateReadets
       a1be2b96b Merged in release_4.1.2 (pull request #314)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       97a0701d5 Merged in JoseHector-Galvez-Lopez/nanopore_covseqbaseini-edited-online-wit-1645559267177 (pull request #316)
       7c1b1f530 Final adjustment to the nanopore_covseq base ini to optimize job submission on Abacus.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      54 commits

       9c2409954 General - general - Adding cpulimit to be used by default on Abacus for all jobs (discussed with PO)
       65971dc02 Merged in covseq_stretenp (pull request #341)
       a03cb2501 Merged in HotFix_dev (pull request #342)
       10b13c0a1 covseq - prepare_report - Adding pandoc from cvmfs
       9f12d3ec9 General - ini - Updating rnaseq and dnaseq ini
       058e9d043 rnaseq - star_index - Auto-Calculating star genomeSAindexNbases based on genome length
       8e4d4d5ae rnaseq - bedtools - Not normalizing if input reads too low
       d1261c0c8 Merged in HotFix_dev (pull request #340)
       129c70907 methylseq - bismark - Hardcoding the multicore value for now
       4c787bde6 General - ini - Standardizing ini base section
       5c0240031 Merged in Paul-Stretenowich/covseqbaseini-edited-online-with-bitbuck-1649436403200 (pull request #336)
       e29820b13 covseq - ini - Fixing bam2fq number of threads hard coded
       75863a056 Merged in HotFix_dev (pull request #334)
       15069c5e9 chipseq - macs2_callpeak + homer_annotate - making genomre_size not required in ini
       2456dec9d Merged in HotFix_dev (pull request #330)
       3f46d5c3f chipseq - homer_annotate_peaks - Adding way to customize genome size in ini for cit
       47ab87a4f chipseq - macs2_callpeak - Adding way to customize genome size in ini for cit
       f9878828f General - scheduler - Addressing PO's comment
       835b0f56f chipseq - homer_annotate_peaks - Debug genome_size
       2942015ca chipseq - Misc - Adding genome size + log to .o file in batch mode
       a57654c5d chipseq - homer_annotate_peaks + homer_find_motifs_genome - Adding way of using genome fasta and not only uscs naming
       ee79ed560 chipseq - homer_make_ucsc_file_bigWig - Adding ini_section variable + IMPORTANT: adding way to skip exit code 141 due to (z)cat | head sending SIGPIPE despite not being an error
       4c9fbc228 chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       4fbf2740e chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       7c972feae General - scheduler - Fixing batch mode job2json
       81ceffe2e General - scheduler - Fixing batch mode job2json
       484353ca6 General - scheduler - Fixing batch mode
       a9d439cc5 General - scheduler - Fixing batch mode
       81187562f chipseq - base.ini - Fixing cluster_cpu old way
       3d5c6a37a Merged in covseq_nanopore (pull request #325)
       5d8cadb6e nanopore_covseq - General - Code cleaning
       564fc0357 nanopore_covseq - General - Code cleaning
       06c695b6c covseq - rename_consensus_header - Updating date from 2021 to 2022 by default
       20df473f0 nanopore_covseq - General - Adding step_wrapper to slurm and batch mode + Adding pipefail in jobs to catch errors in pipes
       5231d853b nanopore_covseq - General - Changing default resources
       9db1946fb nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       cba39ac51 nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       f91a244b1 nanopore_covseq - jsonator - Adding NanoporeCoVSeq for json Nanopore
       e818ccf6d nanopore_covseq - General - Upgrading to python3 ini by default
       27453fbe6 nanopore_covseq - prepare_report - Re-implementing cpulimit to be set at ini level and not automatically
       73b5b5bda Merged in HotFix_dev (pull request #328)
       2d75e4388 dnaseq - General - HotFix
       01c2b49b3 Merged in chipseq_urgent_fix (pull request #326)
       96a5f32f8 chipseq - General - Upgrading to python3
       17f337f20 chipseq - differential_binding - Adding pandoc to module load
       da53a9de5 Merged in chipseq_urgent_fix (pull request #322)
       97dfd745a Merged in covseq_nanopore (pull request #321)
       0215ccf67 chipseq - macs2_callpeak - Fixing typo error
       32d6175b1 nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       c4e2ad5d5 nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       a8017e374 nanopore_covseq - prepare_report - Adding cpulimit for Abacus
       20f047ed2 Merged in Fixing-Issue-#144 (pull request #315)
       2b3c2b552 chipseq - macs2_callpeak - Fixing other_options to be used via ini
       38d0c9a48 chipseq - macs2_callpeak - Fixing other_options to be used via ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       868b49de0 Merged in tp_reports (pull request #348)
       364981844 Merged in tp_reports (pull request #347)
       ea6260fe3 Merged in new_ini (pull request #319)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      20 commits

       98bc236cf tabix ovewright existing index in tumor pair
       166b252fc fix ram in depth coverage gatk
       433f7a32d fix format2pcgr config
       423079142 up mugqic tool version in tumor pair
       5a9f72c7f fix compute_effect with cancer
       e42f7782e revert R for gq_seq_utils_exploratory_analysis_rnaseq
       4eb6fd659 python 3.10.4 for filter ensemble
       edc55691a remove sed from manta_sv
       de2d4252a text file in not compress
       7fabc99bf gatk_indel_realigner does not support multithread
       34d2e11f0 Cleanup dev files
       8e0e9fe34 typo
       cee03217a typo ini files
       d40a229d3 let cpulimit do its stuff on abacus
       4c6514fac let cpulimit do its stuff on abacus
       0e25a16da fix rnaseq denovo crashes
       3a884782a disable reporting to mgcill.hpc with NO_MUGQIC_REPORT variable set
       2924b41be fix nodes set in pbs
       6abbcc4e7 more robust log_report.py
       d93e09478 new ini file setup with common ini for clusters

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       9ea6e292e fixed a bug in input file list of haplotypecaller
       81cce4e44 Added new parameters to change the FDR and P-value of the differential binding analysis (chip-seq)
       c7ff99763 fixed a bug in input file list of haplotypecaller
       900278313 fixed a bug in input file list of haplotypecaller
       75c58e42f fixed a bug in input file list of haplotypecaller
       4d7ba2cde fixed a bug in input file list of haplotypecaller
       c3b27d8b5 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       193ef5a5e Added variant calling using GATK4 to chipseq pipeline

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar1.cedar.computecanada.ca>      6 commits

       c5d8d3e92 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       7605e7cee fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       681a6b864 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       ed358a8d0 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       66e6e90c6 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       3e29640a1 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar5.cedar.computecanada.ca>      5 commits

       40cf8d3fb fixed a bug in input file list of haplotypecaller gatk
       5c903d45b fixed a bug in input file list of haplotypecaller gatk
       7417a6374 fixed a bug in input file list of haplotypecaller gatk
       2aca86606 fixed a bug in input file list of haplotypecaller gatk
       9c2c69883 fixed a bug in input file list of haplotypecaller gatk

  Pubudu Nawarathna Mudiyanselage <pubudu@narval1.narval.calcul.quebec>      52 commits

       1510c0f70 changed job name for merge_gatk
       0442a69aa removed duplicated haplotype calling functions
       08af8f18b Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       ec3bd09ee Modified interval_padding in dnaseq, exome-seq and chipseq
       9f0b3d372 fix padding issue
       36996431f fixed a bug in input file list of haplotypecaller extended the function to change inter-padding resolved merge conflicts
       a40ed3afa Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d111f5f75 fixed a bug in input file list of haplotypecaller gatk
       4c686e83e fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       f595d7196 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       eeda83bf4 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       0b0f2a23a fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       d84db8521 Added variant calling using GATK4 to chipseq pipeline
       abdc21623 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       127c4b4d8 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       7a650a7b2 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e74ac0767 fixed a bug in input file list of haplotypecaller gatk
       af8bb9100 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       3e1af7800 Added variant calling using GATK4 to chipseq pipeline
       896157b1f Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       f523e62dd Added variant calling using GATK4 to chipseq pipeline
       67f4788c5 corrected diffBind path
       fc76ea5f3 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       95592a42c Modified interval_padding in dnaseq, exome-seq and chipseq
       7cee1016f fix padding issue
       1e716736f Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       102b90a3e fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       4ca359b1f Added variant calling using GATK4 to chipseq pipeline
       588b3fb75 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       89e4ed424 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       56bcbc6b5 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       1879f7eb9 fixed a bug in input file list of haplotypecaller gatk
       a33712db1 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       faec78c34 Added variant calling using GATK4 to chipseq pipeline
       d5189ed7e Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       a120b8f6f Modified interval_padding in dnaseq, exome-seq and chipseq
       13b93e9f8 fix padding issue
       6f6e56c4a Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       f2e529915 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       f5d7e5ecb fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       1ba29fb52 Added variant calling using GATK4 to chipseq pipeline
       b2132fa70 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       6d18c042c Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d4616a36e changed mugqic_tools version to version 2.10
       9e33161b9 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       7a7a77f75 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       960a2c9d0 fixed a bug in input file list of haplotypecaller gatk
       7143c0ddf fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       c79776a42 Added variant calling using GATK4 to chipseq pipeline
       b31f874a6 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       133187797 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e051b83d5 changed default python version to python3

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      4 commits

       6174fcb8a Merged in chipseq_snp_hotfix (pull request #335)
       9213853ec Merged in chip_snp_pubudu (pull request #324)
       72ef48442 Merged dev into chip_snp_pubudu
       df9593d9c Merged in chip_seq_issue_150_fix (pull request #327)

  Robert Eveleigh <eveleigh@narval1.narval.calcul.quebec>      12 commits

       b43c90d48 snpEff cancer pair file fix
       7c4ce7115 panel snpeff fix
       853d0b2a5 merge varscan2 fixes
       88c2b0e25 add ram options to gatk_interval_list2bed
       5757d2cac dependency fixes at various steps
       f3d82fcc6 bwa to indel realignment dependency fix
       85f21e45a changed mugqic tools version
       e5a123555 use purple annotated strelka2 output for somatic ensemble calling
       99b328caa update module verions - issues with bcftools
       41b57cea7 cit mugqic_tools version fix
       87768ac52 resource fixes: tumor_pair.exome.ini
       8124f3507 resource fixes: dnaseq.base.ini, tumor_pair.extras.ini and cit.ini

  Robert Eveleigh <eveleigh@narval2.narval.calcul.quebec>      9 commits

       2c41ddc9a bcftool -i/-e fix for tumor pair fastpass varscan2 merge
       0c23107b9 abacus resource fixes for merging and filtering vcf steps
       8a4df12fe memory fixes to merge vcf steps
       85860f894 vardict and varscan2 germline filter fixes
       00dfdf326 fixes to conpair for multiqc intergration
       ff85c1894 tumor_pair.dev.ini fix
       a5c24ac53 change file name : tumor.dev.ini to tumor_pair.dev.ini
       75bc9dcd2 reverted back to python2 for strelka2
       35e5c27fc adding cpsr/pcgr reporting with addition of manta, cnvkit to ensemble protocol.  Step clean up as well

  Robert Eveleigh <eveleigh@narval3.narval.calcul.quebec>      7 commits

       653c5c63b more germline filter fixes
       dc4010010 resource fix to manta and strelka2, and germline filter improvement
       0a3aef5d9 add snpEff QCs and fixed conpair concordance output
       f0749eb87 resource fixes
       9632453d7 add python2 for varscan2 mugqic tool
       6733ac4fc varscan2 one job fix - snp/indel merged vcf now in right place
       12726b387 cpsr/pcgr dependency fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       8948532d8 python3 fixes for support scripts
       f2ea81d3d conpair and ensemble filtering fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      2 commits

       54b97eedb resource fixes
       889346c3c python3 fixes when running on abacus

4.1.2        Thu Feb 17 17:40:12 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       7997ec593 Merge remote-tracking branch 'origin/dev' into release_4.1.2
       c83a283fb Version bump to 4.1.2
       dd9f1853f Updating pipeline READMEs before release
       4bcda9716 Version bump to 4.1.1

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       e5a64db91 Merged in release_4.1.2 (pull request #313)
       99b40bf6a Merged in release_4.1.2 (pull request #312)
       dc6f51e64 Merged in release_4.1.2 (pull request #311)
       7dd399ffe GenPipes - Scheduler : Correcting cluster_cpu parsing for PBS
       7a2aa6db9 Merged in release_4.1.1 (pull request #306)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       c07382b9b Merged in ont_covseq_report_ini (pull request #308)
       46e25d146 Adjustment to the `prepare_report` step resources in the base ini file, to try to prevent this job getting stuck in the queue.

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       4d43f87f4 Merged in patch_4.1.2 (pull request #309)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       989514651 remove cluster_submit_cmd from slurm and pbs ini
       df8c9ba5c fix hic ini
       f45fdba41 remode empty quotes
       a955e743a cleanup queue
       441755268 PBS gets nodes and ppn together

4.1.1        Thu Feb 10 15:12:29 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      10 commits

       dc3ff3c20 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       d2a5407b1 minor README uedit, again...
       3cd7d347e Merge remote-tracking branch 'origin/dev' into release_4.1.1
       21ee999d7 Another minor correction in the README...
       549f6dc25 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       bbc6156de Updating READMEs before release
       150f4a010 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       77e0698d8 Adding missing space in container help message
       be5956049 GenPipes - prep for minor release
       11b8ec2fc Version bump to Release 4.1.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       aeccdc342 Merged in release_4.1.1 (pull request #305)
       c756a500f GenPipes - Tumor Pair : Replacing "xrange" calls by "range" calls to resolve Python3-compatibility issues
       95f8ce822 README.md edited online with Bitbucket
       691b2fdc6 Merged in release_4.1 (pull request #304)
       0204b7cb9 README.md edited online with Bitbucket
       ea152c3e2 corrected typo in README.md for Nanopore_covseq
       ae1c04a9d Merged in release_4.1 (pull request #303)

4.1.0        Mon Feb 7 21:39:25 2022 +0000        77 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      11 commits

       d00f24a31 GenPipes - Release : creation of the release branch
       03e691923 Merge remote-tracking branch 'origin/dev' into release_4.1
       185575e8e Correcting typo in READMEs
       ebc18e3f9 GenPipes - Release : updating READMEs and VERSION
       30e05e5a4 GenPipes - Nanopore : correcting minimap2 step
       7fb2d7c02 GenPipes - BFX : cleaning bash_cmp.py a bit
       f7538c5ae GenPipes - Nanopore : fixing minimap2 and dependencies + allowing duplicate samples (if no duplicate readset) in the readset file
       072871379 GenPipes - Utils : renaming of  to
       cf7b0c65b Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1da114b56 adding .gitignore to .gitignore...
       880fbf481 Version bump to 4.0.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       e503967d9 Merged in release_4.1 (pull request #302)
       8fe012c75 GenPipes - AmpiconSeq : updating resources for dada2 in cedar.ini
       8a0ea558c GenPipes - BFX : fixing "ln" in bash_cmd.py
       44a9dec21 Merged in fix_nanopore_eh (pull request #298)
       44ddf6bc3 Merged in ehenrion/genpipes-chipseq-bashini-edited-with-u-1642626170885 (pull request #292)
       0f84b7827 GenPipes - ChIP-Seq : bash.ini edited with updated versions of software, fixing Issue #127
       6c8219d2b Merged in release_4.0.0 (pull request #286)

  jgalvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       e0c9f5f78 Full squash of covseq_ont branch

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      5 commits

       2841e22b1 Final correction to the reference genome symlink
       323b5d202 Corrected issue with the reference genome link command
       35ad8b757 Added compatibility to ARTIC primer versions 4 and 4.1
       3f3fd3f17 Update version of ncov_tools to 1.8
       8bb76cd84 Added changes to the ini files to allow for ARTIC V4 schemes (and any potential future schemes).

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       fa5b71bc3 Merged in covseq_v4 (pull request #287)
       ce2505546 Added force option to reference symlink creation to avoid crash when re-running report.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      10 commits

       8915d42ac Merged in covseq_nanopore (pull request #300)
       a242cd1b1 nanopore - Cit update - Updating cit ini
       6ab048812 nanopore - pycoqc - Fixing step to match with new bfx
       1f4b6307a covseq_nanopore - Cit update - Updating cit ini
       803b08afc covseq_nanopore - Cit update - Updating cit ini + removing module load in report
       34e981c10 Merged in covseq_nanopore (pull request #299)
       eaa85d8ef covseq_nanopore - Cit - Code cleaning and cit.ini addition
       978697104 covseq_nanopore - Cit update - Updating cit ini
       972678409 covseq_nanopore - General - Renaming covseq and nanopore_covseq class for consistency
       fff0c8b15 covseq_nanopore - General - Fixing argparse type

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       6eddc14a4 Using cit magic for walltime
       c84a3b4cc Merged in cpu_scheduler_agnostic (pull request #294)
       f73f008cc Full squash of covseq_ont branch
       067027810 Merged in config_formater (pull request #289)
       f0047ac38 Merged in watchdog (pull request #290)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      30 commits

       2c5f270bf Add mkdir in FFPE steps
       2860a3b7a cpu_str to node_str in scheduler class
       42052473b update cit.ini
       d67a799e3 update nanopore cit.ini
       2c679b9ce make cluster_cpu scheduler agnostic
       5694384d1 x on tumor_pair.py
       13bd411ee get symlink back in tumor pair
       b919a3e3f get bash_cmd back to old self
       9d722c259 fix regressions
       7d068c2f4 ad tmp to gitignore
       5c198199c update perl in tumor pair
       19aaf7113 cluster_memmory -> cluster_mem
       88a70bc76 tweak tumor pair extra ini
       3fe6b0bb6 missing int casting
       ce9a2e814 pmem for pbs/torque
       f3c406ddb cleanup sanitycheck log
       4357a19c6 SanitycheckError args and kargs
       fe6b300ba add .message to SanitycheckError
       ddc14010f add cluster_mem to ini option
       e0aec1b83 force int on walltime read
       87f84cd2d force int on walltime read
       220b6e0b0 slurm to time_to_datetime
       853ee6ff1 fix hicseq ini
       07f2b77c2 walltime method from slurm
       ad88e69da time delta for pbs walltime
       1227182e0 update utile time delta
       2f679b7b7 update utile time delta
       2fe4ab3e6 update utile for torque
       88abd413d add walltime format for pbs
       333c32a5d rename monitor.sh to watchdog

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       80e3e4d39 chipseq diffbind pca plot update
       dc0680a25 methylseq methylkit min cpg issue fix
       f87dde231 testing updates to DiffBind.R in mugqic_tools
       acaaa6519 modified differential expression variable names to something meaningful
       da9bf6679 fixed issue #129 added differential binding to the atacseq protocol

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       5919bcf98 Merged in chipseq_atac_diffbind (pull request #288)

4.0.0        Thu Dec 9 19:13:55 2021 +0000        232 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      14 commits

       a34c5a193 Merge remote-tracking branch 'origin/dev' into release_4.0.0
       d2f25eb31 updated .gitignore
       c02c96d70 Merge remote-tracking branch 'origin/dev' into release_4.0.0
       9485b9388 GenPipes - README : update epiqc REAMDE
       18dc47137 Merge remote-tracking branch 'origin/dev' into release_4.0
       8ac684765 GenPipes - README : Updated READMEs before new release
       b1e2df04b GenPipes - config : correcting parsing of walltime from config files
       b83fe9fcd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0865f206c GenPipes - DNASeq : correction of sym_link_final_bam dependencies - fixes Issue #119 and Issue #125
       85f8054fe GenPipes - Tumor Pair : fixed some enumerate loop...
       5199ee7f7 GenPipes - Tumor Pair : fixing sym_links steps to avoid duplicate job names
       0e63d7fa3 GenPipes - DNASeq : fixing .bai dependencies
       3453cf00f Version bump to 3.6.3-beta
       d2361a56b Version bump to 3.6.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       963a49a7d GenPipes - python3 : fixes following py2to3 recommandations
       a2274713b GenPipes - Python3 : fixes after last rebase with dev
       4b5f9bda1 Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       5ca19c5f9 GenPipes - start switching to python3
       bfa2a48f8 GenPipes - start switching to python3
       e3f7e9281 GenPipes - Python3 : fixing file opening in bfx/readset.py
       77896893f GenPipes - Python3 : fixing conflicts with bfx/readset.py in dev
       7c41602df GenPipes - start switching to python3
       978a291c2 GenPipes - start switching to python3
       233c1086c GenPipes - start switching to python3

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      4 commits

       2601b7f78 GenPipes - Python3 : removing the shebang from all the bfx scripts
       a99a2e427 GenPipes - Python3 : removing the shebang from all the bfx scripts
       c43eb204c Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       025c3c2fb Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh

  Édouard Henrion <henrione@beluga5.int.ets1.calculquebec.ca>      2 commits

       4463baced GenPipes - Config : removed useless debugging messages
       d47adcfe2 GenPipes - Resources : updated install scripts for fgbio and mugqic_tools

  ehenrion <edouard.henrion@mcgill.ca>      8 commits

       baf52f94d Merged in release_4.0.0 (pull request #285)
       540e030df Merged in release_4.0 (pull request #284)
       fe49a832d GenPipes - BFX : corrected typo in gatk4.py
       bdea57fe8 Merged in tumor_pair_sym_link_fix (pull request #282)
       f5947541e Merged in hotfix_dnaseq (pull request #280)
       34e54c904 GenPipes - RNASeq : minor typo fixes in README.md
       eca4475dd Merged in genpipes_python3_eh (pull request #270)
       34608bc70 Merged in release_3.6.2 (pull request #269)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       b28b4a8c1 remove recalibrated bam name error from some input jobs

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       2b234225e Merged in dnaseq_fixRecal (pull request #276)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      8 commits

       73443c7af Making pandoc report working with pandoc version 2.16.1
       60093e6b2 Making pandoc report working with pandoc version 2.16.1
       c6ab31bdb Merged in Paul-Stretenowich/log_reportpy-edited-online-with-bitbucke-1637772592796 (pull request #277)
       e86925998 log_report.py switching back to py3 and adding Narval as remote
       0ef01e960 EpiQC - First commit after Rami Coles internship
       4e55e0fb5 General - Fixing genome installation grep issue
       06260f0da EpiQC - First commit after Rami Coles internship
       d458b4228 Merge branch 'dev' into epiqc

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      9 commits

       4f89f915d copy contents of readset into inputinfo file(have errors)
       d51b1ceaa copy contents of readset into inputinfo file(have errors)
       97bc3da35 creating inputinfo file after copying original file
       45c917cf9 epiqc - completed chromimpute development - right before remove design file and modify readset file
       95008563e epiqc - completed chromimpute development - right before remove design file and modify readset file
       2bb45e0aa [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       3a86f9488 copy contents of readset into inputinfo file(have errors)
       acfcbc04f copy contents of readset into inputinfo file(have errors)
       7b5b552e5 creating inputinfo file after copying original file

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       43f45759b update GIAContainer version
       9fec33373 changing error log to debug in cit config
       9329e2e85 add cit options for epiqc
       10199ccd0 send telemetry downloded file to /dev/null
       e942d077c fix warning
       99ec4c929 make log_report.py more robust
       c30189c93 update R bioconductor in dnaseq_high_coverage

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      47 commits

       4bee1cab7 bug fix #121. missing module load mugqic_tools 2.8.2
       5e8f3dabf fixed file path issue in seq2fun corrected the mugqic_tools version in epiQC
       a317b6bd6 corrected a typo
       7f7800533 added the reason to use os.remove in a comment
       c4721bc94 addressed Hector's comments - changed seq2fun db path
       ffd0abe5b addressed Hector's comments
       6e2ce1f95 fixed Issue #118: hicseq hic error in interaction_matrices_Chr (mugqic/genpipes) by adding floor division to places which return floating points where intgers needed
       c333a3e82 modified incorrect seq2fun cvmfs db path
       5e928cfc5 added seq2fun cvmfs db paths
       0c9418a91 modified seq2fun db paths
       1569435d2 updated readme file
       22d90bf78 modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       dc4c60880 seq2fun pathway_analysis completed
       de5f87b32 started pathway analysis-seq2fun
       e16274751 improved seq2fun processing
       5d6df1896 added seq2fun to the function as the final step. working well but needs improvements
       bdd622d0e fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       54bf3aa8a added seq2fun protocol completed merge_fastq partially completed seq2fun function
       0c2d56930 addressed Paul's comments
       ef8ae3602 fixed a bug
       2e3c7428f changed mugqic_tools version to 2.8.3
       06e458ef5 fixed bugs when transferring to python3
       fb9093e1b fixed bug in chromimpute convert corrected uuid issue in job2jason.py
       0bd677308 fixed an issue for python 3
       0d1622c57 updated readme file
       9e6254603 testing whether env varibales are retrievable dynamically
       ff73039b9 fixed input files not found error when running the pipeline outside of the project directory
       78af37b95 added cit.ini
       a4dc68c94 fixed issues in final report added comments
       ca3dc8aaa fixed bugs in -o option added IHEC data path
       0bcde6879 modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       affd3b060 seq2fun pathway_analysis completed
       129b56b5a started pathway analysis-seq2fun
       7ef4a6fd9 improved seq2fun processing
       92ac7120f added seq2fun to the function as the final step. working well but needs improvements
       97bdedbfa fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       2d0359b47 added seq2fun protocol completed merge_fastq partially completed seq2fun function
       93811b3fe testing whether env varibales are retrievable dynamically
       72d5f9bc1 fixed input files not found error when running the pipeline outside of the project directory
       b97891f2e added cit.ini
       b29e71777 fixed issues in final report added comments
       cea9f45f1 fixed bugs in -o option added IHEC data path
       15c569593 completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       b86e6e560 completed fixing errors-working pipeline
       d07286d8e corrected up to epiqc final report
       9304cba63 corrected up to epigeek
       562607a92 corrected up to global dist after chipseq design changes

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      59 commits

       e83716689 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       5683bc997 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       faab0eedf Addressed Paul's comments
       b3e81850e addressed Paul's comments
       4130054c3 Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       0225132cb space changes in epiqc.py
       ca91e09bc completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       dbf825497 completed fixing errors-working pipeline
       fc17d37a1 General - Fixing genome installation grep issue
       26837e3cc EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       e3ca2e657 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       9bda3533d EpiQC - Fixed bugs
       e4faee46c EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       1dcd8dd64 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       851b26da0 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       29139a981 General - Fixing genome installation grep issue
       b95fc32e8 resolve merge conflicts
       701f1099f resolved merge conflicts
       d03737e54 resolve merge conflicts
       8852d61cb resolved merge conflicts
       546348680 resolved merge conflicts
       82ff10026 resolved merge conflicts
       d833412de resolve merge coonflicts
       e047a3e35 resolved merge conflicts
       dba6b339b space changes in epiqc.py
       2fe70467f Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       429caf1e5 General - Fixing genome installation grep issue
       06113d178 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       34d3cc942 EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       c00875a3c EpiQC - Metrics thresholds can be modified in base.ini
       35a639627 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       9e3e1f4ec \Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       698c46ff9 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       2086d5af7 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c2af489dd EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       1089984b8 EpiQC - Fixed bugs
       bf60df96c EpiQC - Added signal to noise and epigeec steps
       bced2bd57 EpiQC - Parallelized chromimpute
       c1b462aba EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       4d12401f8 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       0e763c872 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7f8337e77 [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       2b769bff3 resolve merge conflicts
       05b7cd3a8 resolved merge conflicts
       2e219fa5a resolve merge conflicts
       aaed941ee resolved merge conflicts
       cc0bccf13 resolved merge conflicts
       231f6264f resolved merge conflicts
       a7e5a61e9 resolve merge coonflicts
       42a741bbc resolved merge conflicts
       1e5e8af1f Merge branch 'epiqc' of https://bitbucket.org/mugqic/genpipes into epiqc
       6eef5c8d2 resolve merge conflicts
       d679a978c resolved merge conflicts
       eb880f5df resolve merge conflicts
       06fe31e23 resolved merge conflicts
       229653917 resolved merge conflicts
       0e2c9cdbc resolved merge conflicts
       301d2eecb resolve merge coonflicts
       f285a7c14 resolved merge conflicts

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       b3194fa8e Merged in seq2fun_fix_bug121 (pull request #281)
       ffc9c01f8 Merged in epiqc_mugqic_tools (pull request #279)
       a9ac4d618 Merged in epiqc_comments (pull request #278)
       4bd2b26d8 Merged in epiqc (pull request #271)
       c38208b83 Merged in seq2fun_denovo_rnaseq (pull request #272)
       e12b11181 Merged dev into seq2fun_denovo_rnaseq
       0df8920de Merged in hic_python3_issue_fix (pull request #275)
       1d5f2b078 Merged epiqc into epiqc_ss

  rami.coles@mail.mcgill.ca <rcoles@abacus1.ferrier.genome.mcgill.ca>      2 commits

       70f8ec037 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       5feada103 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec

  rami.coles@mail.mcgill.ca <rcoles@abacus2.ferrier.genome.mcgill.ca>      11 commits

       f20a00000 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       74e2a32ac Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       e911f9554 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7a5da0064 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       bc66beb50 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       d0a7b9d81 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6695197a5 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       a040d0b2b added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6a5e99583 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       53c747f5e Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       f56d086f3 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py

  rami.coles@mail.mcgill.ca <rcoles@abacus3.ferrier.genome.mcgill.ca>      21 commits

       690123790 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       88c422d4a EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       f85d37a70 EpiQC - Metrics thresholds can be modified in base.ini
       b14654865 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       51b6bdb3f Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       9e118e512 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       a68aad386 EpiQC - Added signal to noise and epigeec steps
       fe957a6b3 EpiQC - Parallelized chromimpute
       ab1425d29 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       1cbf6d003 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       86c04611b EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       ef3acdeaf EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       642504639 EpiQC - Metrics thresholds can be modified in base.ini
       e9e9f2f8a EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       3ecaebb63 Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       2a773e754 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c800a679c EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       d3bd3e12f EpiQC - Fixed bugs
       402b34ed2 EpiQC - Added signal to noise and epigeec steps
       a6643145d EpiQC - Parallelized chromimpute
       4bf48b977 EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.

  Robert Syme <rob.syme@gmail.com>      1 commits

       f8dfff54f Merged in servername-fix-dev (pull request #274)

  Rob Syme <rob.syme@gmail.com>      1 commits

       a76c958a7 Fix example server name in dev

  Shaloo Shalini <shaloo.shalini@gmail.com>      4 commits

       db1afc378 Merged in ss_mermaid_96 (pull request #247)
       dab640b6f Merged dev into ss_mermaid_96
       6193cdac7 Merged in dev_covseqdoc_ss (pull request #258)
       6dcc0cb7b Merged in epiqc_ss (pull request #253)

  shaloo <shalz@hotmail.com>      14 commits

       9e5476210 Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       4ec66ffc1 Refs #102 removed ChIP-seq pptx link
       783768871 Fixes #102 issues in README.md for epiqc pipeline
       e08c73d81 Refs #103 Paul's feedback in workflow incorporated
       52ca4d125 Fixes #103 covseq workflow added
       f356ffe4b Merge remote-tracking branch 'refs/remotes/origin/dev_covseqdoc_ss' into dev_covseqdoc_ss
       878fcf5e1 Refs #103 covseq doc freebayes update
       8f5625c8d Typo fix
       81ee462de Fixes typo Refs #103
       a62e20bf3 Refs #103 covseq doc freebayes update
       fb32de59c Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       86e5a38fd Refs #102 removed ChIP-seq pptx link
       29ae522d6 Fixes #102 issues in README.md for epiqc pipeline
       43501a4ba Fixes #96 all GenPipes workflows are now coded as mermaid flowcharts and png generated from them

3.6.2        Thu Nov 11 14:14:18 2021 +0000        28 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       3787dea82 Cleaning before release
       c1aef1d93 GenPipes - Readset : improved readset parser so that it creates a temptative readset file with unique IDs when the readset file provided has duplicate readert IDs
       cead55d8f GenPipes - MethylSeq : correted bed2interval_list calls + some minor code reformating
       499a0f63c Version bump to 3.6.1

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      3 commits

       55074b1d6 GenPipes - README : updating READMEs and install scripts before release
       448d7c85a GenPipes - updating version of mugqic_tools in all the pipelines config files
       1c2d32237 GenPipes - Resources : adding some new genome and software install scripts

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       eb019794f Merged in release_3.6.2 (pull request #267)
       2fb1dc1f1 GenPipes - Readset : correcting the parsing of readset files to stop allowing duplicates headsets in file - fixes Issue #113
       72170d962 GenPipes - RnaDeq denovo Assembly.py : updated merge_trimmomatic_stats outputs to fix insilico_read_normalization_all_report dependencies
       50519a286 GenPipes - RnaDeq denovo Sssembly.py : fixed insilico_read_normalization_all_report dependencies
       2e429b1ab Merged in release_3.6.1 (pull request #262)

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       2eca45036 Merged in Pierre-Olivier-Quirion/ampliconseqbaseini-edited-online-with-bi-1635366894762 (pull request #266)
       3d108f933 ampliconseq.base.ini edited online with Bitbucket

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       cc3f279da update all pipelines to pandoc 2.16.1
       cd9c041f0 new pandoc

  Robert Syme <rob.syme@gmail.com>      3 commits

       4ab8e1799 Merged in wget_warn_on_fail (pull request #265)
       d20df3ea1 Merged dev into wget_warn_on_fail
       b809a05cf Merged in rnaseq_protocol_switch_fix (pull request #263)

  Rob Syme <rob.syme@gmail.com>      9 commits

       20154ee66 Fix server name
       31fd5fe57 Test with incorrect server name
       17aea667e switch quote style
       8c4410e9b Intentially introduce an error in the reporting server.
       beb177eeb Warn when wget command fails.
       ca54581c7 Another tiny whitespace fix
       cc33046a7 Tiny whitespace fix
       b3d620c83 Switch at correct location
       2cbdcee19 Fix protocol mixup for rnaseq pipeline

3.6.1        Wed Sep 29 20:53:55 2021 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      6 commits

       27ccb3b1e correct log_report.py before release
       39e1b94d3 Merge branch 'master' of bitbucket.org:mugqic/genpipes into release_3.6.1
       0fae24c17 Updating GenPipes README before release
       df8f9f45e Merge branch 'dev' into release_3.6.1
       96ce81fbb Re-creating READMEs before release
       530576e10 Version bump to 3.6.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       1cda65c20 Merged in release_3.6.1 (pull request #261)
       b472e728d Merged issue_105_fix_eh into dev
       7cd2f0943 GenPipes - Tumo Pair: correcting job name in sym_link_fastq_pair
       4df104246 GenPipes - Tumor Pair.py : fixing typo in sym_link_fastq_pair
       d4613be2a GenPipes - Tumor Pair : fixing issue #106, no more job name overwriting at gym_link_fastq_pair step
       3079e608f GenPipes - DNASeq : fixed bed2interval_list call in gatk_haplotype_caller step
       03c03bf87 Merged in release_3.6 (pull request #259)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       794704250 Merged in Paul-Stretenowich/sambambapy-edited-online-with-bitbucket-1630531961392 (pull request #260)
       69e1a121c Fixing sambamba sort output files: there were a bai set as output that was blocking the "resume" for covseq pipeline.

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       a02620715 more robust log_report
       b941202d4 remove echo debbug in get wrapper

3.6.0        Mon Aug 30 17:55:36 2021 +0000        370 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       cdaafcd7a Merge branch 'dev' into release_3.6
       bea48bd42 GenPipes - updating CHANGELOG and VERSION files (after release)
       2ea1988fb Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       4554eaff2 GenPipes - README : updating version to 3.5.0 in all pipeline READMEs
       1f71288e1 GenPipes - Release 3.5.0 : updating CHANGELOG up to tag 3.5.0

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      16 commits

       9eda417ac Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       79bc0fdcc GenPipes - Tumor Pair : no more attemps to write where the reference files are (in case of readonly FS, e.g. CVMFS...)
       79aa7dc27 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       86c40fdb5 GenPipes - Tumor Pair : fixed Strelka2 jobs with mkdir
       ed205e15f GenPipes - Tumor Pair : fixed coverage_bed in strelka steps
       b7873b93c Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0eafe86bd GenPipes - BFX : fixed GATK4 cluster_crosscheck_metrics command
       a9eb687fd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       ae4cc0840 GenPipes - Tumor Pair ; deleted dev.ini from dev branch
       b6fd7399b GenPipes - Tumor Pair : fixed symlink_fastq_pair command
       41ba98458 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       84f350f37 GenPipes - Install script : version change
       7e07739d0 GenPipes - RNASeq light : fixing typo
       e1c1b5040 RNASeq light : importing bash_cmd fix
       3e29063be GenPipes - Release 3.6 : updating the READMEs in dev
       16be17064 GenPipes - Pipelines : Fixed dnaseq_high_coverage inheritance with dnaseq and rnaseq_light & rnaseq_denovo_assembly inheritance with rnaseq

  ehenrion <edouard.henrion@mcgill.ca>      13 commits

       677e9cf86 Merged in release_3.6 (pull request #250)
       3591ec925 GenPipes - ChipSeq : updated mugqic_tools with latest version
       da2b4102d GenPipes - Ampliconseq: updated mugqic_tools version in base.ini
       40cb3b420 GenPipes - AmpliconSeq : updated mugqic_tools version in base.ini
       92e491944 Rnaseq_light.py: kallisto_count_matrix dependency
       18d37ad17 fixing picard2 add_read_groups call
       952ac0aef GenPipes - MethylSeq.py : Fixed pipeline inheritance with DNASeq
       de571d6bd GenPipes - CoVSeq : setting the gatk4_java_options in base.ini
       11755c1ff GenPipes - BFX : fixing add_read_groups in picard2.py
       fa40033ed GenPipes - MethlySeq : fixing sambamba_merge_sam_files dependency
       2f631833e Merged in pipeline_inherit_fix_eh (pull request #249)
       aa3ca7950 Merged dev into pipeline_inherit_fix_eh
       273fe5485 Merged in release_3.5 (pull request #243)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      24 commits

       a1f3040b5 Merged in covseq_stretenp (pull request #248)
       8374fdb95 Threshold rename consensus change for LSPQ
       d3c3b5fee Code cleaning
       7bc51893b Merged dev into covseq_stretenp
       376368f64 Upgrading covseq_tools version
       e473abe5b Changing regex for finding negative controls for LSPQ
       d0291febd covseq - prepare_table - debug
       0ed94bb62 Merged dev into covseq_stretenp
       b768d9211 Adding bfx for covseq_tools and ncov-tools
       0c8407e61 Upgrading cvoseq_tools version
       6b3f8dec6 Debug
       7dcbf5d28 Degrouping rename_consensus step for ivar and freebayes
       d476c1185 Debug
       00fafc332 Un-grouping freebayes report and ivar report
       1fe319eb4 Un-grouping freebayes report and ivar report
       6b32016a8 Debug
       3023357a6 Changing resources for ncovtools
       58d049722 Debug
       29d800f58 Debug
       8a773f9cd Debug
       5ede73b0e Debug
       ce60e3659 Debug
       84cda2bef Adding freebayes metrics and reporting
       d0592cc6c Renaming prefix parameter in ini_section

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       eeee6b146 Merged in tp_debug (pull request #257)
       fae134938 Merged in monitor_limit (pull request #256)
       6283fd90d Merged in kallisto_rnaseq_light (pull request #254)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      13 commits

       8d09333a3 Bump GiaC to v2.0.3
       55dd94150 patch stolen from Ed to fix recalibration step
       a89de0038 RETRY ret_code corrected
       64b9f79a4 RETRY limit on monitor
       c638ab3d7 RETRY limit on monitor
       bf872974e RETRY limit on monitor
       2c1c89525 RETRY limit on monitor
       207a06943 fix ini for graham covid
       203ea908d update tp ini for cedar
       64db31b2e reducing constraint on fit for integration testing
       a17345dc6 minor fixes before release
       248078c2a fix kallisto dependency
       5583f3fe0 giac updtated to v2.0.2

  P-O Quirion <pioliqui@gmail.com>      2 commits

       fb635da80 more robust log repport
       55b928dc1 make log report more robust went jobs are failing

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       e56c39d85 corrected a typo
       9110f8673 modified readme files for chipseq differential binding step

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       fbdbe38a2 Merged in chipseq_readme_fix (pull request #245)

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      40 commits

       ac5a5c590 Homo_sapiens.GRCh38.ini fixes
       675d5ddcc samtools mpileup spacing fix
       f1827795a cit fixes after dev.ini removal
       15e98828f multiqc ini fix
       1cdfe282f samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       5977acdc2 add portal dir to dnaseq
       3f93719a4 hs37d5 and cit strelka2 fixes
       5a20aa9f6 removal of cit_gatk4.ini, contents moved to cit.ini
       59473cebe tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       d98dd0ce9 merge vardict fix
       a64d34c20 job.sample fixes for jsonator
       63b37d854 vardict exome fixes
       13b7d0f7e germline variant annotator fixes
       d366e782b strelka2 germline and sequenza fixes
       19b187c92 cit adjustments to tumor_pair
       b99937561 beluga.ini fixes
       46f7416ab strelka2 germline fixes
       3d5693d2b removal of samtools, substitute samtools germline for strelka2
       ec213cecd conflict fixes to dnaseq and tumor pair after dev merge
       8cad167d2 fix mpileup dependency chain
       e6d6aa5bb corrections to mpileup bcftools merge and tumor_pair dependencies
       1dc88b918 samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       9e3ac20d4 add portal dir to dnaseq
       4a59db33f portal_dir fix
       4dfa21b9c added portal_dir for cit testing
       528a77be2 hs37d5 and cit strelka2 fixes
       9c66ec284 removal of cit_gatk4.ini, contents moved to cit.ini
       03a602b26 tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       30d7b0756 merge vardict fix
       83b75bf7a job.sample fixes for jsonator
       ec1454724 vardict exome fixes
       e5a614571 germline variant annotator fixes
       fdd20bfbc strelka2 germline and sequenza fixes
       a4f8581e5 cit adjustments to tumor_pair
       7240e3e90 beluga.ini fixes
       849625a49 strelka2 germline fixes
       80f55487d removal of samtools, substitute samtools germline for strelka2
       e92837b1d conflict fixes to dnaseq and tumor pair after dev merge
       9aedd698c fix mpileup dependency chain
       e55e1ffda corrections to mpileup bcftools merge and tumor_pair dependencies

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      45 commits

       6ec0367c8 collect metrics fix - dnaseq
       fd8ee6388 PR fixes to sambamba and multqc
       afc4b8f4b PR conflict fixes
       04c160e45 updated dev.ini
       53cde245a cit fixes
       359c1ae04 purple sanity check fixes
       6d124355c cit fixes to purple and vardict exome
       5da776052 fix for multiple tumors with the same control
       d7724ad5f deliverable fixes - sym link final bams
       63c8bc872 non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       39a26519f bcftools update and filter argument fixes
       67410487c strelka2 germline and ensemble fixes
       b1f739682 sequenza and cit fixes
       c9cef1698 fixes to gatk4 for tumor pair
       bc784a604 gatk4 vsqr cit fix and baf plot for b38
       e148bd4f9 argument fixes for picard imported functions in gatk4 and vqsr fixes
       036809d40 fixes to multiqc
       7b0712be4 exome specific fixes
       059c613a5 updates and fixes for cit
       dbe0ff4f8 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       3c52e6669 remove briaree ini and update dnaseq base
       9343e040c updating beluga ini
       d79e994ce updates to beluga.ini and base.ini for dnaseq
       97793f0a9 ini updates
       f716b0e90 updated dev.ini
       0f7cccc37 cit fixes
       45c9d0bc5 purple sanity check fixes
       882096d50 cit fixes to purple and vardict exome
       97af93a5c fix for multiple tumors with the same control
       49f4f38c4 deliverable fixes - sym link final bams
       17f001eb3 non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       6729217af bcftools update and filter argument fixes
       afba10fe8 strelka2 germline and ensemble fixes
       63264495b sequenza and cit fixes
       e8a4ac9dc fixes to gatk4 for tumor pair
       3348617d0 gatk4 vsqr cit fix and baf plot for b38
       b8ec9bd93 argument fixes for picard imported functions in gatk4 and vqsr fixes
       9b714d09d fixes to multiqc
       ca97ea737 exome specific fixes
       2e20be2fd updates and fixes for cit
       67933483b cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       8b46abeb3 remove briaree ini and update dnaseq base
       4795a96fc updating beluga ini
       04d851fe9 updates to beluga.ini and base.ini for dnaseq
       ab85a8468 ini updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      37 commits

       6cec0eff2 module fix
       7b697cc34 removal of dev files for CVMFS version
       46ef162ca genome ini fixes
       466eb2371 dev and hs37d5 fixes
       367666b7d tumor pair - sambamba markdup fixes
       d5b857fba module fix
       e7a28494f sambamba markdup fix
       f7e4fdbcc strelka cit fixes
       c320d5a4e added strelka2 input dependency for purple purity
       d3c32925a cit fixes to purple and vardict exome
       89c80a9f9 muliple normal fix for fastqc
       00531ff39 dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       d1840a09a fixes and strelka conversion addition to purple
       09d913527 added purity estimate program PURPLE
       aaaff33b1 strelka bed manipulation fix
       f257d83bf vardict exome fix
       f74092d96 conpair and collectHS metric fixes
       2822b5f7f tumor_pair qualimap part 2
       51e88ab54 sym link dnaseq.base into tumor pair
       03e567cd0 updates to b38 variant recal files
       4757879e5 fixes to tumor_pair on beluga
       52e48cefa cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       febb82345 strelka cit fixes
       42fd97c64 added strelka2 input dependency for purple purity
       54fd0cd3c cit fixes to purple and vardict exome
       858af1468 muliple normal fix for fastqc
       bd1f03468 dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       1b2c03fe0 fixes and strelka conversion addition to purple
       0957242ae added purity estimate program PURPLE
       37529e6c5 strelka bed manipulation fix
       cc3178199 vardict exome fix
       2091afb77 conpair and collectHS metric fixes
       8689d63c8 tumor_pair qualimap part 2
       a9888f745 sym link dnaseq.base into tumor pair
       e340eb047 updates to b38 variant recal files
       7a7b46da8 fixes to tumor_pair on beluga
       a449ed5d1 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      28 commits

       58152ae75 varscan version downgrade
       7fd72e181 cit fixes
       3a63b4ad9 fixes to bcftools mpileup for dnaseq mpileup protocol
       b92ebd8a1 hs37d5/GRCh37 fixes for purple
       51c9658c4 fastpass varscan merge fix
       2cc7f0025 cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       93ed2f309 cit fix for conpair
       9672af284 sym link fixes
       c4d0581ef cit fixes to tumor pair
       41ccead32 fixes to dnaseq
       b07afef93 updates to GRCh38 annotation file, module updates, cit fixes
       066fb4997 fixes to deliverable and b38 ini
       034447e9e updates to cit and fixes to one job mpileup steps
       17d4c791e updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       4668a0bcd major fixes to deliverables and completion of beluga test
       9ffad226c fixes to bcftools mpileup for dnaseq mpileup protocol
       180ae35cd hs37d5/GRCh37 fixes for purple
       87dcc11f3 fastpass varscan merge fix
       f6627158b cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       233bc03be cit fix for conpair
       8084b8f71 sym link fixes
       c5840c6b3 cit fixes to tumor pair
       e4193cae4 fixes to dnaseq
       14aa4a300 updates to GRCh38 annotation file, module updates, cit fixes
       de7ed14d7 fixes to deliverable and b38 ini
       951006179 updates to cit and fixes to one job mpileup steps
       510efd1fb updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       cad6781f2 major fixes to deliverables and completion of beluga test

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      43 commits

       05d65f792 dnaseq sym_link fix
       c49e4692e dnaseq ini fixes
       6b63aea53 dnaseq beluga and cedar ini fixes
       040074579 sequenza bug fix
       362539f0a revert back to old sequenza-utils
       3714cd06e fixes to sequenza for exomes, and minor SV fixes
       444e6c266 exome cit fix to sequenza
       9b97f7e87 converted -n to -c
       bf70c513c cit fixes - runtimes
       bfd86d95f config file fixes
       75fd21985 further sanity-check fixes
       30c201265 sanity check fixes
       adb0724e2 multiple pair same control fixes and vardict exome fixes
       f0c3f899f cit fixes to vardict
       34b052def cit fixes for tumor exome
       21e17a243 better cram compression of base recalibrated bams with variantBam
       baf4196f0 sequenza and germline ensemble fixes
       b736a4cf3 sambamba mark_dup fix
       667a4fe7d module fixes
       cee177488 sym link dnaseq cit to tumor_pair cit
       1e7e2e35e tumor pair cit updates
       8d35d7f60 dnaseq sym_link fix
       7da9cee84 dnaseq ini fixes
       3d9111871 dnaseq beluga and cedar ini fixes
       e1a83a360 sequenza bug fix
       f494d73c0 revert back to old sequenza-utils
       ee70ffc89 fixes to sequenza for exomes, and minor SV fixes
       9ffeafa3e exome cit fix to sequenza
       cb1a59cfd converted -n to -c
       00943b326 cit fixes - runtimes
       29c33947d config file fixes
       d3eedbbf6 tumor base fix
       bed7a1d4a further sanity-check fixes
       d36c1e8ea sanity check fixes
       e901171a1 multiple pair same control fixes and vardict exome fixes
       165e341b0 cit fixes to vardict
       10aa2b106 cit fixes for tumor exome
       ceef654e6 better cram compression of base recalibrated bams with variantBam
       3d476d226 sequenza and germline ensemble fixes
       1edc6d7f8 sambamba mark_dup fix
       a3b1208ff module fixes
       caabe6305 sym link dnaseq cit to tumor_pair cit
       b9ed61eba tumor pair cit updates

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      8 commits

       e7bd7451d code cleaning and fixes to exome interval list
       a7d64071d fixes to cedar ini
       cca1a4c6b fixes to symlinks for paired indel realignment
       211a29e2f cedar ini and exome update
       f65b070de code cleaning and fixes to exome interval list
       febded3f1 fixes to cedar ini
       43d40444f fixes to symlinks for paired indel realignment
       09c93a523 cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      4 commits

       c7f3fa9e3 cedar fixes and GRCh38 fixes
       2baaf5641 cedar germline sv updates
       7367777a5 cedar fixes and GRCh38 fixes
       bd82ef133 cedar germline sv updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      46 commits

       21b88961b fixes to sambamba for lumpy - query sorting instead of corrdinate
       dc9750dfb updates to metasv - somatic
       c19ad156a Fixes and updates to reference files
       08f002729 remove testing steps
       aad782e6f updates to SV germline and reference genome tweaks
       491914df7 somatic sv fixes: lumpy and svaba annotations
       d71a64de1 json fix
       03bcbc695 fix to germline SV: breakseq2
       97a4bc18e deleting hidden files
       b24f8f3de merge fixes
       9052826ca GATK4 fixes - bam indexing and markDupSpark
       fc33ed9bd bcftools fixes for tumor pair
       0deb8a949 fingerprint and bug fixes
       57f218ae8 dnaseq - vcftools qc addition: --missing_indv and --depth
       5c82925c9 GRCh38 fixes
       bb43292ff select input for variant caller and fixes to one job calling
       b87e9f501 json and folder updates
       f15fd7409 fixes to sCNAphase
       34c9e47a8 Bug fixes prior to json additions
       ae851f40a merging snv and sv, adding protocols
       c82ce6851 Updates and debug
       57883607b Add set somatic and actionable mutations
       bc0ca0ac6 added multiqc and other tweaks
       c51db5364 fixes to sambamba for lumpy - query sorting instead of corrdinate
       8da9446d3 updates to metasv - somatic
       be214eb30 Fixes and updates to reference files
       f717bb797 remove testing steps
       4c6422379 updates to SV germline and reference genome tweaks
       541f67155 somatic sv fixes: lumpy and svaba annotations
       c795a24e8 json fix
       8cc8e79bc fix to germline SV: breakseq2
       628ce3d2a deleting hidden files
       2a83ec4bd merge fixes
       90ce1f49e GATK4 fixes - bam indexing and markDupSpark
       4e8c75011 bcftools fixes for tumor pair
       4c874cc4c fingerprint and bug fixes
       3cec8c22c dnaseq - vcftools qc addition: --missing_indv and --depth
       6f3ae37a4 GRCh38 fixes
       ebfc76432 select input for variant caller and fixes to one job calling
       4e0b29b56 json and folder updates
       fc717428e fixes to sCNAphase
       a8034e242 Bug fixes prior to json additions
       2b2326710 merging snv and sv, adding protocols
       3c8b861ca Updates and debug
       92db2299e Add set somatic and actionable mutations
       39425ef16 added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      24 commits

       3e96d8182 dnaseq - fixes for VariantBam
       2cd2eb5dc fixes to purple and mkdir vardict fix
       3ba332b82 sequenza WGS fixes
       86f1fa323 fixes to metasv for tumor pair
       3d446d9c2 fixes to sv germline calls
       6f26955cb Single job bug fixes
       893626b5a manta I/O fix and other bug fixes
       9f4109b58 config updates and b38DH added
       9098c4094 dnaseq germline SV updates
       8d737df1e gatk4 updates and bug fixes
       843ef92bf Json related bug fixes
       54fdadc70 Bug fixes and modification derived from initial PROFYLE benchmarking
       5cd673341 dnaseq - fixes for VariantBam
       047e9613f fixes to purple and mkdir vardict fix
       776c53a9b sequenza WGS fixes
       dc8f8202b fixes to metasv for tumor pair
       81b508baf fixes to sv germline calls
       46528cd52 Single job bug fixes
       be5f8d27e manta I/O fix and other bug fixes
       fb73873fe config updates and b38DH added
       9c2a8b0b8 dnaseq germline SV updates
       e86e74ae6 gatk4 updates and bug fixes
       e9b3a2f8d Json related bug fixes
       0201d3861 Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       e9cc4056d seqeuenza fixes and sv fixes
       3786a0369 gatk4 mutect2 updates
       4ae76b754 cit-based fixes to NGScheckmate
       58c3db328 dnaseq qc additions: NGScheckmate and peddy
       9ae421034 seqeuenza fixes and sv fixes
       dc793204b gatk4 mutect2 updates
       0f9e09a59 cit-based fixes to NGScheckmate
       45837b45e dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      6 commits

       f00e0b393 Merged in dev_reveleig (pull request #252)
       00fb9e974 Merged in rebasing_tp (pull request #246)
       08361bc03 Debugging and Guillimin specfic fixes
       f4699979f updates to config
       07d625a8a Debugging and Guillimin specfic fixes
       87f1ef417 updates to config

  Shaloo Shalini <shaloo.shalini@gmail.com>      1 commits

       78f59d52a Merged in ss_wf_chipseq_93 (pull request #244)

  shaloo <shalz@hotmail.com>      1 commits

       8c742ce24 Updates workflow for chipseq -t chipseq case with differential binding step and dependency update Fixes #93

3.5.0        Mon Jul 12 16:41:20 2021 +0000        1071 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      8 commits

       fd8c1cd54 GenPipes - Release 3.5.0 : correcting typo in chipseq README
       038dcb2c9 GenPipes - Release 3.5.0 : correcting typos in READMES
       6e4f4de1c Version bump to 3.5.0
       48f97ec18 GenPipes - Release 3.5.0 : correcting typos in pipeline README files
       f892b2e97 GenPipes - Release 3.5.0 : updating VERSION and README-release files
       82b28e314 GenPipes - Release 3.5.0 : corrected typo in README
       0cf907cd6 GenPipes - Release 3.5.0 : updating version in the pipeline README files
       ddb421a45 Merge branch 'dev' into release_3.5

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      6 commits

       08bd50a8e GenPipes - Resources : updated picard index command in genome installation script to run command on a compute node instead of on the login node
       02124fc01 correcting typo genome install script
       ae646a798 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       69818199e GenPipes - Resources : updates in software and reference install scripts
       0659d35bd GenPipes - Resources : adding software installation scripts
       53dcb8302 GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      47 commits

       339172cf8 Merged in release_3.5 (pull request #242)
       525bfe10e Merged in release_3.5 (pull request #241)
       0d0f7b403 Merged in bash_cmd_for_covseq_update_eh (pull request #236)
       5a54cf99f GenPipes - CoVSeq pipeline : updated call to "zcat" in covseq.py
       2a056bd5a GenPipes - BFX : updated bash_cmd.py
       2809e5efb GenPipes - Install script : added Nanopore pipeline to PATH, fixing Issue #74
       daebeb7a4 Merged in ehenrion/changelogmd-edited-online-with-bitbucket-1620161608327 (pull request #220)
       210d516af Version bump to 3.4.0
       a36b2cbb4 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       7835760cd GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8e2018d45 GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       b4cc95ecc GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       f4534a4a3 GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       e2193e687 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       4b416a0d3 GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b34a1a562 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       05a90a65a GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       5748a6be5 GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       eb4a375b7 GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       c1a8b6bbc VERSION edited online with Bitbucket
       4a0537a89 GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       2fd4f5ea8 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       40125ca00 GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       30917fea0 GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       810feb8c7 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       8911666de GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       3d76f2c16 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       b2b9c458c GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       929424c54 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       879f109cc GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       98c773043 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ad82c94eb GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       31968eea1 GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       8dff24afb GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       dbbaa1126 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       4b617ba2c GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       04f7714a4 GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       f1a49af84 GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       f9c2aefcb GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       72cd3e9d2 GernPipes - DNASeq : corrected merge_filter_bcf outputs
       3b635b239 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       04b24d588 GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       55ab0170e GenPipes - DNASeq SV : fixing delly input error
       09c1bad32 GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       a756e77b0 dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       b514d9f04 GenPipes - Bug fix : corrected dnaseq.cedar.ini
       26b421e52 GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       3297e974b Fix Issue #86
       ff3110dc2 Quick correction to address misleading headcrop parameter in the RNAseq base ini.
       e3d7b04fd rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      344 commits

       ee6fd3087 Merged in covseq_stretenp (pull request #239)
       76bfba93a Fixing sambamba version
       4bc5b9971 hic code improve
       2f1529eed Debug pdf rendering
       b76a994cd Fixing ini
       7a6e9c1d8 Merge branch 'dev' into covseq_stretenp
       985a5c5f2 Increasing resources for ncovtools and adding run_name for cit
       83caadfed Upgrading bcftools
       cdbb41411 Fixing merge_hicrep_scores
       b7ecfc427 Switching to samtools 1.12
       8b29ef3a6 Fixing merge_hicrep_scores
       beab86625 Samtools version update in all pipelines
       00509e830 Merged in covseq_stretenp (pull request #237)
       74a12033d Merged dev into covseq_stretenp
       63cb7a2f9 Fixing indentation error
       22c619af0 Merged in covseq_stretenp (pull request #235)
       16854eb78 Switching to released version of covseqtools
       a69ddeb3e Cleaning comments
       4d85dff49 Removing guillimin and mammouth inis
       f25cf7864 Merged dev into covseq_stretenp
       e955aa2dd Merged in chipseq_design_change (pull request #234)
       5ca220d8d Merged dev into chipseq_design_change
       4c4678bbb Removing name in report/metrics job name
       0b1c5622c Debug
       8c6597933 Debug
       13f9971dc Debug
       9cb766937 Merged dev into chipseq_design_change
       bb6e42107 Merged dev into covseq_stretenp
       a1650eb4f Forcing ncovtools execution
       5f908a77c Fixing metadaya ncovtools file creation
       debf790cd Adding filtered bam to ncovtools
       13b80b056 Fixing ncovtools
       296dfea85 Fixing peak file naming
       0133eb089 Forcing snakemake execution
       9e3647ed2 Forcing snakemake execution
       0b32f4bd7 Forcing snakemake execution
       c85a9fe47 Debug
       a562dd2c6 Changing resources for ncovtools
       59fdae5b4 Debug
       2e4dcbeff Adding dependencies for reporting job
       558cebc4b Purging modules
       2147b8595 Adding more resources for ncovtools
       2366c9035 Debug
       355958da9 Debug
       4c6955bba Debug
       dff338fea Debug
       acca40963 Debug
       e18a0b7de Debug
       c97470b55 Debug
       7f82b371a Debug
       4ee919029 Debug
       a50112bd1 Debug
       23233fb2d Debug
       e143ff192 Debug
       9c6db14ce Debug
       98c023176 Debug
       4e8b8e6a2 Debug
       290e9abd5 Debug
       032d4258a Debug
       72cc78608 Adding GenPipes version accseible inside GenPipes
       49210d2ae Degrouping ncvotools from report generation because of cd for ncvotools
       f8846315d Adding modules reading for reporting
       a3209fda2 Debug
       4717b023b Debug
       bb6b07c74 Debug
       8587d1256 Debug
       612afbf79 Debug
       f7e17ec2d Debug
       b9de420c3 Debug
       2fde5e2f5 Debug
       7a124f18e Debug
       feb236412 Debug
       d15188b22 Debug regex prepare ncovtools
       de760753f Debug
       4b83dab55 Debug
       8b1e5917b Reporting init ncovtools yaml
       b226f89c6 Debug finding files for report
       bc78d4ae9 Debug test dependencies reporting
       d519ac787 Adding missing sections to beluga ini
       60638b969 Debug
       f11091870 Debug
       e482afa1d Debug
       d77e3166b Debug
       11ff7bce4 First commit on gathering metrics for report
       51df06845 Merged dev into covseq_stretenp
       58fd24f4b Merge
       cc0c76153 Fixing minor bug and typo
       2376a9736 Update READMEs
       830153118 Updating inis and fixing homer output folder location
       706ccb32b Debug
       c59a207a2 Debug
       bf45b1d17 Debug
       7b30adfbd Debug
       e79e74ac4 Debug
       c78d0ea57 started to edit the hicseq.py added new file for hicrep
       090de7723 Fixing minor bug and typo
       01d52a21f Fixing multiqc dependencies
       8d17dce06 Debug
       f0a2d0170 Debug
       9dc813e03 Debug
       a3f06f2b0 Debug
       f8dd65165 Debug
       d831843c5 Debug
       9f8d02b3d Debug
       39463cdbf Debug
       daebd07f4 Debug
       898a52308 Debug
       3bdc004a6 Debug
       333524bd8 Debug
       c017225df Debug
       a75dc5abc Changing name of header for a report table
       1eb26a735 Debug
       9eb8c9570 Typo
       cbd8f1e0a Addinx x permission to job2json.py
       6f6a3817e Fixing report naming too long
       9c91a7e3b Fixing report naming too long
       d95ab9261 Fixing report naming too long
       ce26275fc Fixing mugqicValidator.py for new readset format (chipseq only)
       27598c0d6 Switching to latest mugqic_tools release
       fe2ab5094 Changing Ressources
       740102314 Changing Ressources
       7e112120c Changing Ressources
       0a52ab63d Changing module versions
       3b027e9ac Debug
       e5de8a76d Fix
       d1304d013 Increasing picard collect multiple metrics RAM
       fff643846 Debug
       8a0463be6 Debug
       cb4a7e3e7 Fix ressources
       268550baa Fixing annotation peaks
       305bec917 Fixing annotation peaks
       7c56538ab Fixing annotation peaks
       940901769 Adding annotation peaks merge for Narrow peaks only
       8420a20fd Iterating only on Narrow peaks for annotation
       2d4e2e775 Iterating only on Narrow peaks for annotation
       366142d45 Debug report
       9253d109d Debug
       aac7c9cc8 Debug
       9c6dc244f Debug
       c9b7191b2 Debug
       4ac6ca3fc Renaming IHEC metrics section
       36880b945 Debug
       c3858f71c Debug
       b1b68b91f Debug Report
       08b77c98e Debug Report
       916d53f9f Debug
       745fbda19 Debug
       0c8164f05 Fixing report
       c918cb1e6 Debug
       f578ce168 Debug
       ec566be32 Debug
       a3b196de5 Debug
       7768187c0 Debug
       ddf75886f Debug
       4944d36da Debug
       e7d8d28f3 Debug
       b9f35c3bb Debug
       847dddcde Debug
       66a5ca174 Debug
       93ad22ea4 Debug
       9d3cdd81d Changing Output metrics
       723f9ed06 Debug
       fbf10b929 Debug
       a68f0a650 Debug
       68afdd3e3 Debug
       1ac59b201 Debug
       edcfbbfde Debug
       d5548d1b9 Debug
       a3be70d05 Debug
       3b9d879fe Debug
       3b5d819b6 Debug
       9215f3389 Debug
       1d7732f9d Debug
       096571afa Debug
       1e22dc820 Debug
       28bbc7f26 Changing macs2 atacseq param & add macs2 param
       982fb12be Debug
       a86d776b8 Debug
       9a86cbd76 Debug
       5c6e6b4db Debug
       f778cc3b2 Debug
       05187de21 Debug
       2732d265d Debug
       54f47a2dd Debug
       469223788 Debug
       187edbe26 Debug
       6392d86bb Debug
       493624f21 Debug
       553b9baf1 Debug
       eca0813ff Debug
       721c8d102 Debug
       619d6e4ce Debug
       773569d70 Debug
       7836d837a Debug
       8ab157880 Debug
       db9c7868a Debug
       2e8f5d2fa Changing R version for mugqic_tools R scripts
       c90cec996 homer_annotate_peaks
       d4586662e qc_metrics report
       a4c0fbfab Fix ihec_metrics outputs
       66af5287d Fix ihec_metrics outputs
       2bca4810b Fix ihec_metrics outputs
       d33b9e236 Fixing MultiQC report
       cab6ca635 Fixing MultiQC report
       fe600e9aa Fixing MultiQC report
       5754975f6 Fixing MultiQC report
       f31bf7114 Fixing MultiQC report
       fcdc378b2 Fixing MultiQC report
       d3559e6ff Improving MultiQC report
       f424ffac3 Debug
       a1dcb44a9 Debug
       44b2647fb Debug
       fc1881914 Debug
       777c716c4 Debug
       1217e0b63 Debug
       cbfcd261f Debug
       a5a3ca354 Debug
       479504f54 Debug
       44a9b8082 Debug
       425dc3cd3 Debug
       a853292f2 Debug
       eaec07b0a Debug
       f80b4f2e6 Debug
       81f990a44 Debug
       ecc34289e Debug
       89b57b8ac Debug
       36e4a2af5 Debug
       bbf0fdc1a Debug
       d4d8d6482 Debug
       335218079 Debug
       985004658 Debug
       1f291cafd Major changes IHEC metrics test
       047fe5f38 Major changes IHEC metrics test
       924074699 Major changes test
       23c925645 Macs2 report changes
       acfc78227 Major change test
       76bdd83e0 Major change test
       392234ff9 Major change test
       a27e5b69a debug
       d54a5903a debug
       4ab48dbcf debug
       cfaad2128 debug
       5daadd72f debug
       4ebada3fd debug
       5a16691c1 debug
       807e2645b debug
       fd604bd84 debug
       c6347e79d debug
       aa9e88a5d debug
       440c29048 debug
       96d01f102 debug
       e52535672 debug
       7785d116f debug
       dc09b1545 debug
       a5c60c6cb debug
       7be46d4dc debug
       0faaea9a4 debug
       4db524801 debug
       74b73f352 debug
       2df9a1c52 debug
       6b583d5a8 debug
       e3961cb82 debug
       622832e3a debug
       8b2d57a1b debug
       fa973d0ec debug
       64fd4ed91 debug
       6a9071f54 debug
       ab1f02402 debug
       40cb701b1 Fix test
       f1dce58b9 Major readset change test
       349a1a3d7 Fix
       095b13e27 Fix
       17519c7fc Fix
       0baff781c Filtering after merging and changing naming to fit with other pipelines
       cab4a18b8 Increasing default resources and adding options for markdup
       ffae951c3 Fixing beluga ini
       33c1ab908 Switching to sambamba markdup for IHEC
       b689ca7fa Fix
       9a59883e5 Fix
       c23b19389 Fix
       9d1aa0a62 Fix
       fb2ad2a9a Fix
       879f850af Fix
       8cc4c624c Options becomes optional
       7cd7a814f Fix
       804de5cc7 Fixing typo
       9795af9dc Fix
       11cd56d48 Fixing sambamba merge
       fb4bbcfa7 Typo
       823693c4f Adding mkdir
       11e06a2fe Fixing temp typo to tmp
       273863710 Fixing job
       6fe444698 Fixing bash import
       8fbe44ae1 Fixing minor bug and typo
       d9dbf29ac Adding missing sections to beluga ini
       5f74c2f69 Renaming freebayes consensus and ncovtools quick_align
       4e6a1bdc6 Degrouping rename consensus jobs because of quast issue creating infinite dependency loop
       25f9e362d Debug
       a0ed01a75 Debug
       07ab3d110 Adding missing sections to beluga ini
       cc21ed983 Debug
       a113c663c Merged dev into mgi_stretenp
       490857bf2 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       c844ceaf2 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       4a793a05a started to edit the hicseq.py added new file for hicrep
       0e991a7d8 started to edit the hicseq.py added new file for hicrep
       72c9387ca hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9458e13dd hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ac314f69e hicseq completed adding basic features of the hicrep analysis.
       cd8fbb9e9 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       086a8cbde hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       f9893e590 hicseq completed adding basic features of the hicrep analysis.
       4c811af21 Cleaning
       e87bf2899 Debug
       a360050de Debug
       d95592529 Changing freebayes calling to skip bgziping
       cdee35236 Reducing Ressources
       8856ee9b3 Debug
       e77bd4ba5 Debug
       30e04a4e2 Debug
       670c71e56 Debug
       292c07341 Fixing Jared's script
       d799b7a28 Reducing Ressources
       e2bb1dbd3 Fixing Jared's script
       e6f4df817 Fixing freebayes variant calling
       a7e2b43eb Fixing Jared's script
       7a6ef9994 Fixing bcftools consensus
       5a88dcd53 Fixing freebayes
       114d9bbaa Fixing freebayes
       1a3db5e3a Fixing freebayes
       46d283916 Fixing freebayes/bcftools
       612e50688 Adding bcftools consensus creation following Jared's script
       92df26f28 Adding bcftools consensus creation following Jared's script
       e73f67df7 Adding freebayes variant calling
       5ce0e1800 Grouping rename jobs into 1 job
       ccae3393b Grouping rename jobs into 1 job
       f3daf1022 Grouping rename jobs into 1 job
       6d850fb5e Adding freebayes variant calling
       edf5a762c Grouping rename jobs into 1 job
       76fbe6bf4 Adding freebayes variant calling
       7ab5da228 Adding freebayes variant calling
       00acf27ba Adding freebayes variant calling
       c700579bb Adding freebayes variant calling
       eefed60de Fixing flagging bug and updatinh year to 2021

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       1719a1f08 Merged in update_container (pull request #238)
       906d0fd3e Merged in slurm_cgroup_problem (pull request #233)
       1dd042cd3 Merged in temp (pull request #232)

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      387 commits

       773c8d4b9 hicseq completed adding basic features of the hicrep analysis.
       ffe7becc2 Completed developing hicrep and quasar analysis
       9e790f672 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6e6b3659c Completed hicrep analysis except the creation of the graph
       dbb237f85 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       72bde2f12 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       167673f59 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2471a46c7 hicseq completed adding basic features of the hicrep analysis.
       68aece4fe hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       5903eca7a Added pairwise combination for samples
       e7322592d started to edit the hicseq.py added new file for hicrep
       af7bfd464 hicseq completed adding basic features of the hicrep analysis.
       ad7a803ee hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       23746b5f8 Added pairwise combination for samples
       e50669f32 started to edit the hicseq.py added new file for hicrep
       96ea65eb5 hicseq [hicrep.py] - Corrected typo
       a952c1f01 hicseq [hicrep.py] - corrected R_TOOLS path
       3dee5c653 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2e52c82a4 hicseq [hicrep.py] - Corrected typo
       fc9ccc44a hicseq [hicrep.py] - corrected R_TOOLS path
       e93c722d5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       507604578 hicseq [hicrep.py] - Corrected typo
       1404bb4b6 hicseq [hicrep.py] - corrected R_TOOLS path
       70a9c0069 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5e3e5ffc4 hicseq [hicseq.py] - corrected file after rebase
       b96c0e715 hicseq [hicrep.py] - Corrected typo
       dc6f85b8c hicseq [hicrep.py] - corrected R_TOOLS path
       750cc6c1d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       92263db44 hicseq [hicseq.py] - corrected file after rebase
       b708ea3c9 hicseq [hicrep.py] - Corrected typo
       b209f6a17 hicseq [hicrep.py] - corrected R_TOOLS path
       3ba0747ff [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2bf8e9502 hicseq [hicrep.py] - Corrected typo
       8936ce6d2 hicseq [hicrep.py] - corrected R_TOOLS path
       9ac5f101b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1a5bc63c5 hicseq [hicrep.py] - Corrected typo
       c01ee521a hicseq [hicrep.py] - corrected R_TOOLS path
       1805da837 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       19df2bab0 hicseq [hicrep.py] - Corrected typo
       33e3949d6 hicseq [hicrep.py] - corrected R_TOOLS path
       77ebe87e3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       85294d650 hicseq [hicrep.py] - Corrected typo
       64840a7c9 hicseq [hicrep.py] - corrected R_TOOLS path
       3013d61df [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       734c2793f hicseq [hicrep.py] - Corrected typo
       9a33ab8ad hicseq [hicrep.py] - corrected R_TOOLS path
       833061195 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       94b562bfc hicseq [hicrep.py] - Corrected typo
       005b6cc7f hicseq [hicrep.py] - corrected R_TOOLS path
       28076a5cd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       858b05706 Added pairwise combination for samples
       e76217b20 started to edit the hicseq.py added new file for hicrep
       8e835eb85 hicseq [hicseq.py, readme.md] - modified readmefile
       6bf04f8dd hicseq [hicseq.py] - corrected file after rebase
       589ae72f2 hicseq [hicrep.py] - Corrected typo
       331299614 hicseq [hicrep.py] - corrected R_TOOLS path
       2c6c21167 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       54108f490 hicseq [hicrep.py] - Corrected typo
       ca68cf408 hicseq [hicrep.py] - corrected R_TOOLS path
       c1d6c9c0e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       01817b641 hicseq [hicrep.py] - Corrected typo
       a68d82b4f hicseq [hicrep.py] - corrected R_TOOLS path
       310b5eda9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       785bd9dbe hicseq [hicrep.py] - Corrected typo
       4db60e27a hicseq [hicrep.py] - corrected R_TOOLS path
       0baffa5c2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fc107ccc2 hicseq [hicrep.py] - Corrected typo
       8dbe281d0 hicseq [hicrep.py] - corrected R_TOOLS path
       f2721650d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f57e58ab hicseq [hicrep.py] - Corrected typo
       11c5af450 hicseq [hicrep.py] - corrected R_TOOLS path
       e544d4f3a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0c90edc12 hicseq [hicrep.py] - Corrected typo
       ec0feb848 hicseq [hicrep.py] - corrected R_TOOLS path
       165aeeb4f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d4d21a3f1 hicseq [hicrep.py] - Corrected typo
       f665976e6 hicseq [hicrep.py] - corrected R_TOOLS path
       9a1295e56 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9778a4a8e hicseq [hicseq.py] - corrected file after rebase
       10365ba01 hicseq [hicrep.py] - Corrected typo
       f470412cc hicseq [hicrep.py] - corrected R_TOOLS path
       a7f0ba323 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0e1ddc65b hicseq [hicrep.py] - Corrected typo
       956a72de4 hicseq [hicrep.py] - corrected R_TOOLS path
       1cc549834 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9b7e8b321 hicseq [hicrep.py] - Corrected typo
       06cdfcf40 hicseq [hicrep.py] - corrected R_TOOLS path
       d4c14d238 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0ce13dfe7 hicseq [hicrep.py] - Corrected typo
       279dc5298 hicseq [hicrep.py] - corrected R_TOOLS path
       7c304c72a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       81ad33d28 hicseq [hicrep.py] - Corrected typo
       87bbe6682 hicseq [hicrep.py] - corrected R_TOOLS path
       0c5127b14 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b67b6b046 hicseq [hicrep.py] - Corrected typo
       95a4a9767 hicseq [hicrep.py] - corrected R_TOOLS path
       4c25fc252 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6d7a744cf hicseq [hicrep.py] - Corrected typo
       d96ad6fdd hicseq [hicrep.py] - corrected R_TOOLS path
       a1768c57f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c52093ed Added pairwise combination for samples
       5aff10cf1 started to edit the hicseq.py added new file for hicrep
       55987dd49 hicseq pipeline [changed the step order]
       525c7d397 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       f1811aa31 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       5fe3a4931 hicseq [hicseq.py, readme.md] - modified readmefile
       df1ac4b22 hicseq [hicseq.py] - corrected file after rebase
       cfd6f524b hicseq [hicrep.py] - Corrected typo
       2b5daacaf hicseq [hicrep.py] - corrected R_TOOLS path
       796fe27bf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e9d2d2db3 hicseq [hicrep.py] - Corrected typo
       34b96297f hicseq [hicrep.py] - corrected R_TOOLS path
       beb3b076d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       af6d8aca5 hicseq [hicrep.py] - Corrected typo
       c97468afc hicseq [hicrep.py] - corrected R_TOOLS path
       658e574d6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5b2fedf9a hicseq [hicrep.py] - Corrected typo
       f3c653f27 hicseq [hicrep.py] - corrected R_TOOLS path
       1690ce438 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       89371a2d4 hicseq [hicrep.py] - Corrected typo
       ef5659485 hicseq [hicrep.py] - corrected R_TOOLS path
       ff4b8ceb7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       161296e2e hicseq [hicrep.py] - Corrected typo
       61e8b115c hicseq [hicrep.py] - corrected R_TOOLS path
       a147aa988 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       db542ae8b hicseq [hicrep.py] - Corrected typo
       229f0656c hicseq [hicrep.py] - corrected R_TOOLS path
       919463a22 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       077e2ba07 hicseq [hicrep.py] - Corrected typo
       8c2d82bec hicseq [hicrep.py] - corrected R_TOOLS path
       7dd995d6b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       78ee1f785 hicseq [hicseq.py] - corrected file after rebase
       42120ffde hicseq [hicrep.py] - Corrected typo
       e21766cda hicseq [hicrep.py] - corrected R_TOOLS path
       2aebe64d8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       51e0f2aaf hicseq [hicrep.py] - Corrected typo
       3d8ad7977 hicseq [hicrep.py] - corrected R_TOOLS path
       59ea5d009 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       20ea3a824 hicseq [hicrep.py] - Corrected typo
       df2b34a36 hicseq [hicrep.py] - corrected R_TOOLS path
       6bdaef1f1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c4727f54c hicseq [hicrep.py] - Corrected typo
       bd80a99e2 hicseq [hicrep.py] - corrected R_TOOLS path
       7d004943f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2c056190f hicseq [hicrep.py] - Corrected typo
       593a82da5 hicseq [hicrep.py] - corrected R_TOOLS path
       60581eb38 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2f0362a1d hicseq [hicrep.py] - Corrected typo
       a2dfa4f4a hicseq [hicrep.py] - corrected R_TOOLS path
       04faafc29 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       85569620b hicseq [base.ini] - updated mugqic tools version to 2.3.1
       57dba1bd2 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       61fce92a4 hicseq [quasar_qc - corrected module loading in matrix restructuring
       6bc2c0cea hicseq [hicrep.py] - Corrected typo
       23217cace hicseq [hicrep.py] - corrected R_TOOLS path
       fc77bfe88 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8190866b4 hicseq [hicrep.py] - Corrected typo
       708979437 hicseq [hicrep.py] - corrected R_TOOLS path
       e5be8fdd7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       70754df08 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       a83bbe3f4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5b5ac2f31 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b9a9dfdc [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       29fefa5ea [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9c80cb01e Completed developing hicrep and quasar analysis
       14a6698e7 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d4f176efe [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       69310a4f6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       38d65fda9 Completed hicrep analysis except the creation of the graph
       97d738117 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       a7272ddb4 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       99142a9a9 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3c16aa5ca hicseq completed adding basic features of the hicrep analysis.
       366611dce hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       a26866b4a Added pairwise combination for samples
       54e99f84f started to edit the hicseq.py added new file for hicrep
       ad051374d hicseq [hicseq.py, readme.md] - modified readmefile
       57f13dfb6 hicseq [hicseq.py] - corrected file after rebase
       727d579ff hicseq [hicrep.py] - Corrected typo
       cabd6b11f hicseq [hicrep.py] - corrected R_TOOLS path
       4f95216fb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d9d53d590 hicseq [hicrep.py] - Corrected typo
       2d50a471f hicseq [hicrep.py] - corrected R_TOOLS path
       0b6541ef4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b56a96a8e hicseq [hicrep.py] - Corrected typo
       a41dc59de hicseq [hicrep.py] - corrected R_TOOLS path
       d6a751fcd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       146c8921e hicseq [hicrep.py] - Corrected typo
       552d8f973 hicseq [hicrep.py] - corrected R_TOOLS path
       e529d26bd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       743948f1d hicseq [hicrep.py] - Corrected typo
       bed015461 hicseq [hicrep.py] - corrected R_TOOLS path
       c334f10bc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b286fcd50 hicseq [hicrep.py] - Corrected typo
       548cfc05d hicseq [hicrep.py] - corrected R_TOOLS path
       316c44f76 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       95bb7ed9f hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2cd7fd621 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       dad2a18fa hicseq [quasar_qc - corrected module loading in matrix restructuring
       912cca583 hicseq [hicrep.py] - Corrected typo
       6f1f5207d hicseq [hicrep.py] - corrected R_TOOLS path
       8120389bd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3816804a8 hicseq [hicrep.py] - Corrected typo
       8581ac0b9 hicseq [hicrep.py] - corrected R_TOOLS path
       5371ed0ed [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       226a4f75f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       be5df368f [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       f5876b20a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4465ebe0e [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c5c922fe3 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       24c34f0fb Completed developing hicrep and quasar analysis
       c7d0aa89a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       722820889 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d1c440a69 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       b4bc9bfe6 Completed hicrep analysis except the creation of the graph
       324fbc2f7 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       ff16ded05 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       3d54c9fd8 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1e3be791f hicseq completed adding basic features of the hicrep analysis.
       e35ac3d87 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ac2c2082b Added pairwise combination for samples
       61bd4a7a5 started to edit the hicseq.py added new file for hicrep
       c950a435b hicseq [hicseq.py] - corrected file after rebase
       9d755feb5 hicseq [hicrep.py] - Corrected typo
       9cfb4eef9 hicseq [hicrep.py] - corrected R_TOOLS path
       9bb97058c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b744cedee hicseq [hicrep.py] - Corrected typo
       4477d4678 hicseq [hicrep.py] - corrected R_TOOLS path
       74fafa09a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2a0449336 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       acb090f56 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       a7578f998 hicseq [quasar_qc - corrected module loading in matrix restructuring
       0e8398afc hicseq [hicrep.py] - Corrected typo
       9d27afa6c hicseq [hicrep.py] - corrected R_TOOLS path
       7fdc12cfc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8448782a4 hicseq [hicrep.py] - Corrected typo
       09b742080 hicseq [hicrep.py] - corrected R_TOOLS path
       e5f684992 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cfdb9d463 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5af4d8054 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       1ae80800a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e48da77f1 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       46ba5a853 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       8a20dbad6 Completed developing hicrep and quasar analysis
       3e6f7634f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b9ae84853 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       60e5df89d created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3c2f201 Completed hicrep analysis except the creation of the graph
       e1dea78c5 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b086ca533 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       114fbd179 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       edd76c8c9 hicseq completed adding basic features of the hicrep analysis.
       3549fd07a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e3720267a Added pairwise combination for samples
       8d80a5b04 started to edit the hicseq.py added new file for hicrep
       18f4e1531 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4fbc1fdcc hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       c12f62a87 hicseq [quasar_qc - corrected module loading in matrix restructuring
       5e9ddcd45 hicseq [hicrep.py] - Corrected typo
       194793082 hicseq [hicrep.py] - corrected R_TOOLS path
       c71ed68bb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       157c1b51a hicseq [hicrep.py] - Corrected typo
       5937a022a hicseq [hicrep.py] - corrected R_TOOLS path
       b99ff65f0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71c1160f4 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       03f2537a2 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c94034087 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a8e16f945 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       298874c44 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       bf24afb65 Completed developing hicrep and quasar analysis
       e5fe17ff4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e2fc36049 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f67be0b69 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f1b5f5ade Completed hicrep analysis except the creation of the graph
       619bbf169 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       c74201f59 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0b72f780b hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f998fcf7 hicseq completed adding basic features of the hicrep analysis.
       a06d72680 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       009d28680 Added pairwise combination for samples
       9bc074eb0 started to edit the hicseq.py added new file for hicrep
       bcc23d6bf hicseq [base.ini] - updated mugqic tools version to 2.3.1
       20949c263 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       81afdcce5 hicseq [quasar_qc - corrected module loading in matrix restructuring
       08c675799 hicseq [hicrep.py] - Corrected typo
       51c8f6644 hicseq [hicrep.py] - corrected R_TOOLS path
       9c2ee30da [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       adc92e7a8 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b020b249e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       623d8bbc0 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       7e1ecf64c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       6e22be6d1 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f2815e41b Completed developing hicrep and quasar analysis
       ce39693c5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f2e21fcf4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4895f68ca created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       7a13a0385 Completed hicrep analysis except the creation of the graph
       0288fb9a3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       708eeec30 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e97faa549 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1ae276853 hicseq completed adding basic features of the hicrep analysis.
       646858dbf hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dab3e7d4b Added pairwise combination for samples
       1e9e150b3 started to edit the hicseq.py added new file for hicrep
       dd821443c hicseq [hicrep.py] - Corrected typo
       0c8400be0 hicseq [quasar_qc.py] - Corrected deletion by mistake
       18a704dfd hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       18747b859 hicseq [hicrep.py] - corrected R_TOOLS path
       5759e5820 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cc9a31847 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       dcfaa23c6 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       df458a53d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       6718520df [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       ab0c23aaf [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       63ae5faa9 Completed developing hicrep and quasar analysis
       506c95bb0 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       78c57d891 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b87ff3f90 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       18258411f Completed hicrep analysis except the creation of the graph
       643905bd2 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       aa8395fb9 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c28516168 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f10a6b30 hicseq completed adding basic features of the hicrep analysis.
       662a52c7c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d43a3c190 Added pairwise combination for samples
       3472e7a58 started to edit the hicseq.py added new file for hicrep
       4a08538b4 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       e2ede562b [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97d248a14 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       bb6887e8c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       815e22a64 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       854c1c0d6 Completed developing hicrep and quasar analysis
       c648da3c9 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0f8505886 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54cefb76b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3912eb1 Completed hicrep analysis except the creation of the graph
       c1e103535 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       dc414bd2a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       943a7c188 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       6be727c2a hicseq completed adding basic features of the hicrep analysis.
       1e5d44462 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       3af40a81c Added pairwise combination for samples
       3180796fc started to edit the hicseq.py added new file for hicrep
       e67b07979 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6e359557d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f0dd71996 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       9556bef7c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f17594ed1 Completed developing hicrep and quasar analysis
       a384d9a55 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       339a8017c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       aefd4bbd0 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       98a1cac68 Completed hicrep analysis except the creation of the graph
       2ae354314 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       255526b81 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ced9d1389 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       48a964715 hicseq completed adding basic features of the hicrep analysis.
       3b2f387e8 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       77e52e21a Added pairwise combination for samples
       0894881c7 started to edit the hicseq.py added new file for hicrep
       ffb9b96af [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7507b5a86 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       21f3a9cb5 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       2c4d8bbe6 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       72686d186 Completed developing hicrep and quasar analysis
       7cb650caa [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       abd1855f3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       045f82334 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       70c513a42 Completed hicrep analysis except the creation of the graph
       1cf3bb2e5 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       06f58f8a2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       adef366f6 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       42888b6fb hicseq completed adding basic features of the hicrep analysis.
       bf251fd69 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d156fd81b Added pairwise combination for samples
       0ff825209 started to edit the hicseq.py added new file for hicrep
       dd17ce73b Completed developing hicrep and quasar analysis
       db3118620 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       467e68295 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1e11460c9 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       e9fac8721 Completed hicrep analysis except the creation of the graph
       a2afec73b hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       234855f29 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       db9095ff1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b8f0c5535 hicseq completed adding basic features of the hicrep analysis.
       a96d6c1d8 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       66d63f73b Added pairwise combination for samples
       a69557d9b started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      55 commits

       b0eb60b87 update metylseq ini
       7ae859df6 crosscheck_fingerprint with EXIT_CODE_WHEN_MISMATCH 0
       f2171d231 add new line on genpipes_file.write()
       74130cd09 line explicite new lines at the end of fp.write()
       ea3c24b79 fix for pbs
       f8645cc5e ajust graham and cedar  ini
       129a42731 ajust graham ini
       c5a687597 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       b4cad788e cedar and graham ini maintenance
       8c6d99c44 tweak monitoring
       b81ace0dc force -c on tumor pair ini files
       d5e8b0781 remove unused import from hicseq
       8a46c262a log_report formating
       b94022d29 remove the verbose stuff that I am not monitoring anymore.
       6b44921e4 add missing args* in super call
       06ff02bac typo
       6d6b28bc6 update get_wrapper to v2.0.1
       b65dd6e14 remove decomissioned hpc ini files
       8fcf649e2 from output_file to genpipes_file
       ba89a19e9 --command_file options
       99292d871 remove -o option for output file
       96ee30565 used image v2.0.0.rc in for container
       3ecb2fb34 force TMPDIR to default to /tmp if not set
       1e3cc77a5 make sure output file is executable
       a6e173092 added option --output for command file
       66cf0f53c update form -n to -c
       e10f09038 removing useless shebang
       77d1a0dcf extra space in samtools sort option
       10a9ffcf7 force pace between stringtie options
       047f332d6 make sure gemini annotation use the ini defined tmp dir
       fbf55522d force log report to python 2 :
       80d950858 error in path for chunk_genpipes csplit call
       b7eb42a04 force space before bash line continuation
       164e6f8a8 missing space in varscan script
       1d3bf9179 update form -n to -c
       a40dce4f7 removing useless shebang
       350f89880 extra space in samtools sort option
       5e67fa4c9 force pace between stringtie options
       ba056af73 make sure gemini annotation use the ini defined tmp dir
       944fcc7d0 force log report to python 2 :
       a47fe8604 error in path for chunk_genpipes csplit call
       9b19e7aa0 force space before bash line continuation
       f50600c74 missing space in varscan script
       a087236b9 update form -n to -c
       5a8d51212 excluding wget call from monitoring loop
       1dd1bdbfc force loading mugqic python 2.7 on cuffmerge
       26183abe8 force loading mugqic python 2.7 on multiqc 1.7
       29735e863 more verbose went -d/--design is needed for contrast
       0a3560173 revert on indel aligner
       c47384b08 cleanup ini for beluga
       c65298ead revert on indel aligner
       b18ca362e cleanup ini for beluga
       85147169b remove dbSNP for oxog ini block fix gatk 4 module name
       d5140585a fix picard_collect_oxog_metrics dbSNP vcf
       d630a6033 remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       019c25726 Corrected hicseq.py file deletions
       0c780b0da Corrected hicseq.py file deletions
       9a30832a9 Corrected hicseq.py file deletions
       8f1db9ace Corrected hicseq.py file deletions
       621cddab9 Corrected hicseq.py file deletions
       a9dce2dee Corrected hicseq.py file deletions
       94bb5a848 Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      14 commits

       cabe4d025 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       772bb6865 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ef45d3a70 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       9bda9db81 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       93eb03e9d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a7f50bfdb hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       f79cacf71 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86ddffc25 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       7e5a282ab hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4de276e89 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3c940bb4e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       757f80fce hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       56c88759a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86e066cb7 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      82 commits

       982fdd833 addressed ED and PaulS's feedback on PR- step 3
       10920cb2c addressed ED and PaulS's feedback on PR- step 2
       498544e2c addressed ED and PaulS's feedback on PR- step 1
       d18311ca5 bugfix - unsuccesful
       804bc6c29 added step information to the md file and chipseq.py
       f7365e14f bugfix - unsuccesful
       f5decb532 chipseq.py - remove space in first line
       b85ff2417 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       ed99e0a82 passed more parameters to R script
       a66963774 chipseq.py - changed output directory
       06a3f01bb chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       0893090c6 adjusted alignment in functions
       b67ba87b5 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       48e42ca17 added Rmarkdown rendering added more paramters
       775730d41 chipseq.py - changed output directory
       1bd300e1f chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       101e2f61f chipseq.py - remove space in first line
       0aff721bc differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b84e23751 fixed results files creating at custom output folder
       5830eaeef added Rmarkdown rendering added more paramters
       cd4a5f106 chipseq.py - changed output directory
       8f70abf35 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       17b3232d7 bugfix - unsuccesful
       1d3b46394 chipseq.py - remove space in first line
       07ef003a1 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       80f6cc0d6 fixed results files creating at custom output folder
       0d13c4f9d added Rmarkdown rendering added more paramters
       6f50be2c9 chipseq.py - changed output directory
       81695c6f4 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       e3132b272 add diff_bind.py in bfx/
       47d33c2f7 added step information to the md file and chipseq.py
       53dc92357 bugfix - unsuccesful
       c2efc7d24 chipseq.py - remove space in first line
       b6416d605 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       bb677ee68 passed more parameters to R script
       15f231386 chipseq.py - changed output directory
       cc24da1ae chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       485596995 adjusted alignment in functions
       196706e5d skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0843c3a74 added Rmarkdown rendering added more paramters
       d25e5d0d9 chipseq.py - changed output directory
       8c17263d1 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bf900c736 chipseq.py - remove space in first line
       62bd5a5e4 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b88a2e3b4 fixed results files creating at custom output folder
       84979112b added Rmarkdown rendering added more paramters
       6c780b250 chipseq.py - changed output directory
       ed75a1eeb chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bc9d0f2ce added step information to the md file and chipseq.py
       43ab77e19 bugfix - unsuccesful
       994e168ad chipseq.py - remove space in first line
       7413ea5b3 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       241404f3a fixed results files creating at custom output folder
       a908dc37d added Rmarkdown rendering added more paramters
       26497f526 chipseq.py - changed output directory
       8e8d8e96b chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       f71777e0b add diff_bind.py in bfx/
       1d1d959a6 bugfix - unsuccesful
       528bd3873 chipseq.py - remove space in first line
       7805b91a0 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       1a7e95766 fixed results files creating at custom output folder
       3f58f2966 added Rmarkdown rendering added more paramters
       39732342c chipseq.py - changed output directory
       414155ce6 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       9f12db9b6 add diff_bind.py in bfx/
       56028a81e chipseq.py - remove space in first line
       662aefaa1 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       a8d2fcf0a differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       456d11434 fixed results files creating at custom output folder
       92e6fd779 adjusted alignment in functions
       4f6653275 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       e09c78fb9 added Rmarkdown rendering added more paramters
       b65cf7e72 chipseq.py - changed output directory
       303efcda1 added paramters to ini files
       f927e563c chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       abd9dd769 add diff_bind.py in bfx/
       96170dff7 passed more parameters to R script
       28ceb0ef0 chipseq.py - changed output directory
       696560f71 added paramters to ini files
       f46eb293d chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       b92624883 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       641d393a4 add diff_bind.py in bfx/

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      14 commits

       eabb5bb7a refer previous commit
       5f8957f33 refer previous commit
       c225c6d10 refer previous commit
       1a091f54a refer previous commit
       540862d7b refer previous commit
       818f474ad refer previous commit
       21ce258d2 refer previous commit
       94281fc6d refer previous commit
       c7774f2f4 refer previous commit
       d005007ef refer previous commit
       a1cf1c0d9 refer previous commit
       4adc9ca6d refer previous commit
       3d25f6095 refer previous commit
       25cddb0ab refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      63 commits

       e34a5f380 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       94b62b910 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       304b53471 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       580338df6 fixed results files creating at custom output folder
       8fbd60e44 added paramters to ini files
       06b9ddf5c started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       ec6762a25 passed more parameters to R script
       58aa030ef added paramters to ini files
       99101ea02 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       3cc862239 rebasing
       1ed107919 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       618704193 adjusted alignment in functions
       c6c31f2ae skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       038c730d5 passed more parameters to R script
       5502008fe rebasing and resolving conflicts
       bb6343d52 added paramters to ini files
       ee7af7d2d started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       f1498d8cb rebasing
       745719ba4 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       30ab71425 adjusted alignment in functions
       d86864ff5 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       1699baef2 passed more parameters to R script
       9ba71cf08 rebasing and resolving conflicts
       ebc7a88bb added paramters to ini files
       3055bd9ea started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       bfe5a10d0 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       637357744 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       044565e53 fixed results files creating at custom output folder
       ae1abe78c added paramters to ini files
       ad27f00ba started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       687821214 passed more parameters to R script
       41ba25a3f added paramters to ini files
       b350cf2ba started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2010b544f rebasing
       993dfa20f chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ffd815434 adjusted alignment in functions
       80dcb4c7d skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       186b7c7c8 passed more parameters to R script
       ac4e38d59 rebasing and resolving conflicts
       d412a8a73 added paramters to ini files
       428296b86 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2d681ff4c Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       36aa86938 rebasing
       20ea21c84 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       8949c9eb5 adjusted alignment in functions
       f9c99d134 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       ac7dbf095 passed more parameters to R script
       3c253135d rebasing and resolving conflicts
       e4d75d2ce added paramters to ini files
       4c3b94017 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       d8c2155e3 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       0d24cc410 rebasing
       35fff1fa4 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ca86b46cb adjusted alignment in functions
       1076487b1 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0bb499ba7 passed more parameters to R script
       716493df7 rebasing and resolving conflicts
       1e38ab5ce added paramters to ini files
       4e22bde84 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       78edf353f Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis resolved merge conflicts after pull
       0ff1cfd67 passed more parameters to R script
       3195decef started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       9cceef8dc corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       6199e349c Merged in chipseq_diff_analysis (pull request #240)
       ce590061e Merged dev into chipseq_diff_analysis

  Shaloo Shalini <shaloo.shalini@gmail.com>      9 commits

       b8f6ff31f Merged in ss_wf_83 (pull request #229)
       63c5ba1f9 Merged in ss_wf_82 (pull request #228)
       e38d01422 Merged in ss_wf_81 (pull request #227)
       0e0b2a232 Merged in ss_wf_80 (pull request #226)
       afdf55ab2 Merged in ss_wf_79 (pull request #225)
       55262398e Merged in ss_wf_78 (pull request #224)
       cd6f047aa Merged in ss_wf_77 (pull request #223)
       02a08a421 Merged in ss_wflow_76 (pull request #222)
       eff450c0d Merged in ss_wflow_75 (pull request #221)

  shaloo <shalz@hotmail.com>      27 commits

       0c5149a0a Fixes #83 rnaseq denovo workflow diagram updated for v3.4.0
       116e1018f Fixes #82 rnaseq_light workflow updated for v3.4.0
       0e00fcdb1 Fixes #81 rnaseq workflows are now current to 3.4.0
       94b163c4a Fixes #80 methylseq pipeline workflow is now current with v3.4.0
       2811acfed Fixes #79 hicseq pipeline update for v3.4.0
       7679be511 Fixes #78 dnaseq workflow for -t mpileup updated v3.4.0
       37c5ca7ff Fixes #77 dnaseq pipeline -t mugqic workflow update for genpipes v3.4.0
       531bbf5fd Fixes #76 dnaseq_highcov pipeline workflow updated v3.4.0
       943151547 incorrect update should be for #76 cleaning up
       b5dd068c1 Fixes #76 dnaseq_higcov pipeline workflow updated in sync gpv3.4.0
       a9998d404 Fixes #75 updated amplicon sequence qiime and dada2 workflows for v3.4.0
       3bd1caa34 Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       8cbae5ef5 Refs #67 dnaseq -t light option color feedback addressed
       1befe5328 Fixes #67 Refs #54 dnaseq -t light workflow schema added
       da2257566 Refs #65 chipseq workflow color updated, report links added as per feedback
       0f63c1ef2 Fixes #65 Refs #54 chipseq pipeline updated
       c32984777 Refs #68 color added as per feedback
       78895de9d Fixes #68 Refs #54 nanopore pipeline workflow schema added
       1103637f3 Refs #66 cleanup after color update
       742de7661 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       18010fb42 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       dcbf14d8e Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       97aa943ae Refs #53 covseq workflow arrow step 6 -> 9 removed
       2d4a60d37 Refs #53 Paul's review inputs incorporated
       08e8068b6 Fixes #53 covseq.py workflow diagram added
       8911c4cbe Refs #53 added covseq pipeline schema workflow draft under review by Paul
       fb8f6f6eb Fixes #51 update covseq pipeline readme to v3.3.0

3.4.0        Thu Apr 29 19:24:01 2021 +0000        784 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       ff50d9505 GenPipes - RNASeq : corrected genome_index_folder refencing in star_align

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       2b72afb3a GenPipes - Release : removing Tumor pair prior release 3.4
       9de7a7797 Merge remote-tracking branch 'origin/dev' into release_3.4
       ed27f3314 GenPipes - Resources : adding software installation scripts
       8c9f61390 GenPipes - Resources : adding software installation scipts
       c2d24183e GenPipes - Resources : adding reference genome installation scripts
       3bc6a8013 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       67cc49944 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       823f62647 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       e42f490a6 GenPipes - Resources : adding Sus_scrofa genome and ivar & kent installation scripts
       d76314e02 GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      79 commits

       ccc89c112 Merged in release_3.4 (pull request #219)
       9667e558b Resolve Issue #31
       41c771b57 Merged dev into eh_cit_correction
       5fd600c39 GenPipes - BFX : corrected ini section for star_index --outTmpDir
       982a3050b GenPipes - BFX : update STAR calls to use --outTmpDir
       7a6549de1 Merged dev into chipseq_design_change
       0ca16958e Merged dev into chipseq_design_change
       23f6e6d47 Merged in eh_cit_correction (pull request #207)
       daf9eb8a6 GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       2e1b75e92 GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       0db8b4449 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       9d5c8af87 Merged dev into chipseq_design_change
       09b50e334 Merged dev into eh_cit_correction
       af89867e5 GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       024b25173 Merged eh_RNAseq_star_correct into dev
       81e394bfe GenPipes - RNASeq : corrected genome path in star_align
       0dbd8958a Merged eh_fix_callhome_fail_exit_code into dev
       a8bfffdad GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b652f60de Merged eh_RNAseq_star_correct into dev
       14f5ca516 VERSION edited online with Bitbucket
       4a2a45dff GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       ef26b52d5 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       01d00902e GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       73f22e32c GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       6ae113957 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       a3781df20 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       fde1dbc93 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       533b3fa8a GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       433375287 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       9cda20cd2 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       ceb95945b GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ed1a35b40 GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       6c61ca59c GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       d250236cd GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       a5dc00cd8 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       155115709 GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       404c9312f GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       ba277e429 GernPipes - DNASeq : corrected merge_filter_bcf outputs
       1028ccc34 VERSION edited online with Bitbucket
       82e6690ed Merged in ehenrion/version-edited-online-with-bitbucket-1617908341194 (pull request #205)
       e2249aad4 VERSION edited online with Bitbucket
       0f75fb63e Merged eh_samtools_cram_output_ini_fix into dev
       74d7c56c4 GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       a45ed7252 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       0da38f2cd GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       5ded7767e GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       b64372520 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       33a4b9636 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       053b208b9 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       da787e924 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       713bd648a GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       716110de3 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       68c129412 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       26ca812c7 GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       065df8d3c GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       98f539131 GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       2d7892d04 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       cbe399cdd GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       b74a4763d GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       938318a81 Merged eh_cit_correction into dev
       4d076784a GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       c095c6cae GernPipes - DNASeq : corrected merge_filter_bcf outputs
       c5d36c968 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       0d9ab3eac GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       138081b8f GenPipes - DNASeq SV : fixing delly input error
       79f6e9495 Merged in ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004 (pull request #202)
       ce2014cc7 Merged eh_RNAseq_star_correct into ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004
       f1a6ffd19 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       f53266bd4 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8a51f4b39 Merged in ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543 (pull request #200)
       e0cf3d830 Merged eh_RNAseq_star_correct into ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543
       66c603f10 GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       77ced18c9 GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       176669b42 Merged in eh_fix_delly_issue52 (pull request #199)
       2fde5cb4a GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       277904f51 dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       fdbef3adc GenPipes - Bug fix : corrected dnaseq.cedar.ini
       5824422cd GenPipes - Bug fix : correcting indentation in illumina_run_processing.py
       9cfee2e7a dnaseq.py edited online with Bitbucket : corrected protocol assgignation

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       a0a8e94f7 rnaseq.base.ini updated to a newer version of STAR
       22f3dfad2 Merged in rnaseq_star_update (pull request #195)
       06a8d473d rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      219 commits

       d0835ca2b Merged in chipseq_design_change (pull request #210)
       f7655cbd8 Merged dev into chipseq_design_change
       ce5c0079b Fixing atacseq
       2e0ee42ae Merged in chipseq_design_change (pull request #206)
       ef03626c8 Update READMEs
       376858baf Merged dev into chipseq_design_change
       182612e76 Updating inis and fixing homer output folder location
       f27a211b9 Merged dev into chipseq_design_change
       0ce599179 Debug
       4c0b11f13 Debug
       3f6fcc6de Debug
       9cbe87c7c Debug
       586708f07 Merge branch 'dev' into chipseq_design_change
       5cb6a403c Debug
       1980df200 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       318284c34 started to edit the hicseq.py added new file for hicrep
       03e1ebc64 Fixing minor bug and typo
       48290070d Fixing multiqc dependencies
       9acc4eb7c Debug
       9b72ff083 Debug
       03c9cb675 Debug
       a159e49db Debug
       2c9cf5666 Debug
       76664e891 Debug
       7e90a8669 Debug
       4a4ae351d Debug
       5718604fa Debug
       4c01435ca Debug
       bb6ee3f15 Debug
       aaf4809ba Merged dev into chipseq_design_change
       ae98ed2e2 Debug
       d105f7199 Debug
       9c42079da Changing name of header for a report table
       786ac8147 Debug
       d98c7ba81 Typo
       8b1b565fe Addinx x permission to job2json.py
       2de4ae65c Fixing report naming too long
       d53e173d7 Fixing report naming too long
       a3d412e0c Merged dev into chipseq_design_change
       40bcb0f26 Fixing report naming too long
       a399b3d6c Fixing mugqicValidator.py for new readset format (chipseq only)
       2249d6a0d Switching to latest mugqic_tools release
       3b77185d5 Changing Ressources
       34524b8c3 Changing Ressources
       af360e36e Changing Ressources
       4738bda02 Changing module versions
       6d5e26a19 Debug
       860b44d8e Fix
       ef12823e7 Merged dev into chipseq_design_change
       6512f2f14 Changing qualimap bamqc section naming to a parameter in bfx
       9b4d61cee Increasing picard collect multiple metrics RAM
       b264e7dca Debug
       030fb31a5 Debug
       8de52dd6b Fix ressources
       dacaec867 Fixing annotation peaks
       12e5d8669 Fixing annotation peaks
       9d9a50f59 Fixing annotation peaks
       351f5c312 Adding annotation peaks merge for Narrow peaks only
       0a2e098cd Iterating only on Narrow peaks for annotation
       fd3bd7464 Iterating only on Narrow peaks for annotation
       58da5b0be Debug report
       687969cd4 Debug
       955fe4bb3 Debug
       d7560d905 Debug
       e9e09de8f Debug
       4b8ab24c6 Renaming IHEC metrics section
       8f9ee567d Debug
       b73ce4987 Debug
       43bebdfb5 Debug Report
       cbf2b319b Debug Report
       d89432452 Debug
       5dc81e888 Debug
       98d61cfc5 Fixing report
       3b179a45f Debug
       74cd55db9 Debug
       d3e03b287 Debug
       6545cbd6d Debug
       9f90aacf0 Debug
       b94c8a1e9 Debug
       03279f65f Debug
       830e45650 Debug
       79cb6542a Debug
       19a0024e3 Debug
       cf9d2527a Debug
       505b1caea Debug
       bbf285204 Changing Output metrics
       e654a974e Debug
       43105ff5c Debug
       2092de56d Debug
       3c6c7b45a Debug
       0ab7349ab Debug
       2ee5418ca Debug
       a30cf4a7d Debug
       ec6d257a0 Debug
       610181d7f Debug
       a85c9742e Debug
       bfd2e2fe4 Debug
       f21e83eaa Debug
       3fb54184b Debug
       d77307e13 Debug
       ee397ef87 Changing macs2 atacseq param & add macs2 param
       8c4274b4a Debug
       a28d27589 Debug
       a64ce5328 Debug
       a689539ed Debug
       9587ca818 Debug
       e370f8201 Debug
       0d3741d76 Debug
       db78c589b Debug
       f5d866d1c Debug
       796562560 Debug
       fdf3d9bf5 Debug
       7847ac055 Debug
       47da9ef24 Debug
       6a8b1ba9c Debug
       dddcf8840 Debug
       161353f83 Debug
       73790f902 Debug
       188496e06 Debug
       361b17542 Debug
       f5940edaf Debug
       e502b6fba Changing R version for mugqic_tools R scripts
       45ef768f5 homer_annotate_peaks
       6a70faf51 qc_metrics report
       49e8bfe86 Fix ihec_metrics outputs
       cc35a5dc0 Fix ihec_metrics outputs
       29c13621c Fix ihec_metrics outputs
       7c6807182 Fixing MultiQC report
       12212bd6c Fixing MultiQC report
       c46e42dc3 Fixing MultiQC report
       73ced5d58 Fixing MultiQC report
       08c1f4356 Fixing MultiQC report
       b403bfbdf Fixing MultiQC report
       bdc3e4679 Improving MultiQC report
       cde658a17 Debug
       528f0476d Debug
       cccc2e92f Debug
       a545b4ab8 Debug
       2b79de344 Debug
       acd95d270 Debug
       c0d1a393f Debug
       92d584e25 Debug
       0c54b4ed2 Debug
       c919c9001 Debug
       e9c595591 Debug
       64657e3be Debug
       69857a490 Debug
       2683f8859 Debug
       e8a45bf5e Debug
       2570791b7 Debug
       6bf3d107d Debug
       524a5c0c1 Debug
       985fa272e Debug
       3b6ce4a52 Debug
       ee40fde5f Debug
       6eedc38f4 Debug
       34383e759 Major changes IHEC metrics test
       2385c51a1 Major changes IHEC metrics test
       72f30ff75 Major changes test
       9dfe33e83 Macs2 report changes
       705e66a70 Major change test
       e23426a5a Major change test
       7d787c142 Major change test
       58b992c58 debug
       b97215113 debug
       9f4c020bc debug
       4a14c4261 debug
       e69cf403a debug
       20a8959ac debug
       d800207ac debug
       d0cebc148 debug
       38616a43f debug
       3f0fea1a6 debug
       fcca28563 debug
       61ff7d36b debug
       d147d39a3 debug
       10e30302d debug
       70c0c5d54 debug
       acc9d2c8c debug
       68ea5f437 debug
       a19d64f93 debug
       76f515579 debug
       8a9225951 debug
       9e6c5349e debug
       3b9bf7a40 debug
       1e35666d1 debug
       b2f165b7d debug
       bf5579128 debug
       4677380f2 debug
       5cc5350d8 debug
       f0d9f3857 debug
       b306bbd46 debug
       1d7d87168 debug
       891234246 Fix test
       166b57b56 Major readset change test
       8db202b74 Fix
       abfcb3d8c Fix
       ae7964acc Fix
       c9bc75067 Filtering after merging and changing naming to fit with other pipelines
       aae460942 Increasing default resources and adding options for markdup
       9f540c4d5 Fixing beluga ini
       dd795e633 Switching to sambamba markdup for IHEC
       27e1a0877 Fix
       ecc127db5 Fix
       e01ccf2a8 Fix
       b299062aa Fix
       f8cf12259 Fix
       264059ca7 Fix
       df2604b8d Options becomes optional
       1748eb6c5 Fix
       c09bf2f24 Fixing typo
       adf2e7fa1 Fix
       6c4cbd38c Fixing sambamba merge
       9cf59889e Typo
       fc6f4cb8c Adding mkdir
       f060825c4 Fixing temp typo to tmp
       aca35ab97 Fixing job
       95e56619f Fixing bash import
       f5508117c Fixing minor bug and typo

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      385 commits

       03dbf551a hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2bd980875 hicseq completed adding basic features of the hicrep analysis.
       7c19ba7af hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3b8a1916 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       32ffa5633 hicseq [hicrep.py] - Corrected typo
       5321e34a1 hicseq [hicrep.py] - corrected R_TOOLS path
       7d56a2117 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fe8103ab5 started to edit the hicseq.py added new file for hicrep
       9f24b59a1 hicseq [hicrep.py] - Corrected typo
       78d9e5965 hicseq [hicrep.py] - corrected R_TOOLS path
       ac8395d56 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       56f9aa365 hicseq [hicrep.py] - Corrected typo
       f2f843b59 hicseq [hicrep.py] - corrected R_TOOLS path
       2e3a2e80a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b242698fb hicseq [hicseq.py] - corrected file after rebase
       b9c2899fd hicseq [hicrep.py] - Corrected typo
       56fb9ebb8 hicseq [hicrep.py] - corrected R_TOOLS path
       bd160f60d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3726aea98 Completed developing hicrep and quasar analysis
       ef6f77b4f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e34261 Completed hicrep analysis except the creation of the graph
       e28dd55a9 started to edit the hicseq.py added new file for hicrep
       5d82058d0 hicseq [hicseq.py] - corrected file after rebase
       9f494ca9a hicseq [hicrep.py] - Corrected typo
       a5f7a9039 hicseq [hicrep.py] - corrected R_TOOLS path
       e098a617d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3defd1e3d hicseq [hicrep.py] - Corrected typo
       3acb37df9 hicseq [hicrep.py] - corrected R_TOOLS path
       f18061a56 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ec44a6bd hicseq [hicrep.py] - Corrected typo
       ac36fca92 hicseq [hicrep.py] - corrected R_TOOLS path
       911d66932 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8b851f363 hicseq [hicrep.py] - Corrected typo
       aee4c68a5 hicseq [hicrep.py] - corrected R_TOOLS path
       a4ac85cef [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       766f8db63 hicseq [hicrep.py] - Corrected typo
       8f9a52d3d hicseq [hicrep.py] - corrected R_TOOLS path
       e2a85e877 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eca4b3568 hicseq [hicrep.py] - Corrected typo
       cae9d46a6 hicseq [hicrep.py] - corrected R_TOOLS path
       c742e6c9d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9bd1c6ee0 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       ea5095ff4 hicseq [hicrep.py] - Corrected typo
       0d7d6c085 hicseq [hicrep.py] - corrected R_TOOLS path
       99e1052f1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5deb5597a hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       31d1b1bbe hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dc822e8de hicseq [hicseq.py, readme.md] - modified readmefile
       b828c729e hicseq [hicseq.py] - corrected file after rebase
       142dfc305 hicseq [hicrep.py] - Corrected typo
       1fd823887 hicseq [hicrep.py] - corrected R_TOOLS path
       b0d0045a2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ffaba4c66 hicseq [hicrep.py] - Corrected typo
       16866fb8a hicseq [hicrep.py] - corrected R_TOOLS path
       5c52f5ab7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ca854211 hicseq [hicrep.py] - Corrected typo
       07ce5ba99 hicseq [hicrep.py] - corrected R_TOOLS path
       9886c9e82 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       83e131415 hicseq [hicrep.py] - Corrected typo
       4ee15eb65 hicseq [hicrep.py] - corrected R_TOOLS path
       8ae9494cf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7de77365d hicseq [hicrep.py] - Corrected typo
       cc6564d51 hicseq [hicrep.py] - corrected R_TOOLS path
       ebcc1a11f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1eafa496 hicseq [hicrep.py] - Corrected typo
       380dfd472 hicseq [hicrep.py] - corrected R_TOOLS path
       4d9f85314 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c2f186d71 hicseq [hicrep.py] - Corrected typo
       654198f79 hicseq [hicrep.py] - corrected R_TOOLS path
       c96c47572 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bf2449c99 hicseq [hicrep.py] - Corrected typo
       342b6f586 hicseq [hicrep.py] - corrected R_TOOLS path
       5831718ac [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f6e3c5d9 hicseq [hicseq.py] - corrected file after rebase
       bac4d8d15 hicseq [hicrep.py] - Corrected typo
       81688407c hicseq [hicrep.py] - corrected R_TOOLS path
       e4fdcd0ed [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7250a4f29 hicseq [hicrep.py] - Corrected typo
       3942f5664 hicseq [hicrep.py] - corrected R_TOOLS path
       10b31a811 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b6805e643 hicseq [hicrep.py] - Corrected typo
       03ed962c7 hicseq [hicrep.py] - corrected R_TOOLS path
       fa3d60b0b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9c4734f52 hicseq [hicrep.py] - Corrected typo
       dccfad08e hicseq [hicrep.py] - corrected R_TOOLS path
       907e93256 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b3e8d9a94 hicseq [hicrep.py] - Corrected typo
       cbc0e992f hicseq [hicrep.py] - corrected R_TOOLS path
       a15ba0c94 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3a6a1358d hicseq [hicrep.py] - Corrected typo
       7d17401f6 hicseq [hicrep.py] - corrected R_TOOLS path
       3c495bf12 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e83336555 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       5dd2d25e2 hicseq [hicrep.py] - Corrected typo
       62b39434a hicseq [hicrep.py] - corrected R_TOOLS path
       9a8f22a75 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d95fe9f25 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff71f6fdb hicseq completed adding basic features of the hicrep analysis.
       35c32698d Added pairwise combination for samples
       85571a9d8 started to edit the hicseq.py added new file for hicrep
       72be0065c hicseq pipeline [changed the step order]
       0800b9250 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       3363d0752 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       e5b3070e4 hicseq [hicseq.py, readme.md] - modified readmefile
       d4dbe19a9 hicseq [hicseq.py] - corrected file after rebase
       adf7ceb83 hicseq [hicrep.py] - Corrected typo
       9f9e07123 hicseq [hicrep.py] - corrected R_TOOLS path
       848150431 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d286d5dd2 hicseq [hicrep.py] - Corrected typo
       8f8c6128c hicseq [hicrep.py] - corrected R_TOOLS path
       85a1cf964 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8eca7d4ae hicseq [hicrep.py] - Corrected typo
       9229a936f hicseq [hicrep.py] - corrected R_TOOLS path
       5a60eff29 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef77769e4 hicseq [hicrep.py] - Corrected typo
       467b8e941 hicseq [hicrep.py] - corrected R_TOOLS path
       c0559eb92 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c8ff0042f hicseq [hicrep.py] - Corrected typo
       6e1367f12 hicseq [hicrep.py] - corrected R_TOOLS path
       2266dd99f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       999455da7 hicseq [hicrep.py] - Corrected typo
       240eaf33f hicseq [hicrep.py] - corrected R_TOOLS path
       f2ab622da [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d0755352e hicseq [hicrep.py] - Corrected typo
       28e544c5e hicseq [hicrep.py] - corrected R_TOOLS path
       5312b5602 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4161b94ff hicseq [hicrep.py] - Corrected typo
       7ddd445a6 hicseq [hicrep.py] - corrected R_TOOLS path
       5a24278ca [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1dd26f408 hicseq [hicseq.py] - corrected file after rebase
       3d526b79d hicseq [hicrep.py] - Corrected typo
       10196c8f1 hicseq [hicrep.py] - corrected R_TOOLS path
       66cf3f7d4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2b72cc9a7 hicseq [hicrep.py] - Corrected typo
       0a1ae4bd8 hicseq [hicrep.py] - corrected R_TOOLS path
       ec8ec0634 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       52ceb831c hicseq [hicrep.py] - Corrected typo
       5691c9485 hicseq [hicrep.py] - corrected R_TOOLS path
       33907621a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       071ed2979 hicseq [hicrep.py] - Corrected typo
       6f23e5fd4 hicseq [hicrep.py] - corrected R_TOOLS path
       29dd05991 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bdd9b8789 hicseq [hicrep.py] - Corrected typo
       ee032b051 hicseq [hicrep.py] - corrected R_TOOLS path
       0325e110b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f9cb9baff hicseq [hicrep.py] - Corrected typo
       e2f2d6bc9 hicseq [hicrep.py] - corrected R_TOOLS path
       5c92379e3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0f4de0bbb hicseq [base.ini] - updated mugqic tools version to 2.3.1
       f1361d63f hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       fe8285893 hicseq [quasar_qc - corrected module loading in matrix restructuring
       589286e27 hicseq [hicrep.py] - Corrected typo
       9dcbc8338 hicseq [hicrep.py] - corrected R_TOOLS path
       85b4bbd5d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4cc5884ee hicseq [hicrep.py] - Corrected typo
       26d8fde46 hicseq [hicrep.py] - corrected R_TOOLS path
       689a76b66 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5529f4d14 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       2d0af167d [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c60cd1577 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a5cfa02c0 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c34a220a4 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       89a55c1d8 Completed developing hicrep and quasar analysis
       0b74a0aa1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54c0ea89c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       45c545a2b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       3ea32860c Completed hicrep analysis except the creation of the graph
       f6f2a519f hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       034732add hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       777005f77 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       36ab2e9b4 hicseq completed adding basic features of the hicrep analysis.
       f9a4ef6b2 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7fff19eee Added pairwise combination for samples
       5ae1358ef started to edit the hicseq.py added new file for hicrep
       875a2afdb hicseq [hicseq.py, readme.md] - modified readmefile
       afb029fe3 hicseq [hicseq.py] - corrected file after rebase
       ae2799ddb hicseq [hicrep.py] - Corrected typo
       22c4aeaad hicseq [hicrep.py] - corrected R_TOOLS path
       a8b525d61 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5a84543e3 hicseq [hicrep.py] - Corrected typo
       29ef57793 hicseq [hicrep.py] - corrected R_TOOLS path
       8e99a729c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0af5ff67 hicseq [hicrep.py] - Corrected typo
       2c514bf02 hicseq [hicrep.py] - corrected R_TOOLS path
       9c984a06a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65e4495a9 hicseq [hicrep.py] - Corrected typo
       d6b6bdabf hicseq [hicrep.py] - corrected R_TOOLS path
       15c7831ce [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a98b3265c hicseq [hicrep.py] - Corrected typo
       bfc0ee4f6 hicseq [hicrep.py] - corrected R_TOOLS path
       2dd0caa5c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6070378dc hicseq [hicrep.py] - Corrected typo
       b7f35aee9 hicseq [hicrep.py] - corrected R_TOOLS path
       c739410ff [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       473304ea6 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b536a3740 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2e9b26618 hicseq [quasar_qc - corrected module loading in matrix restructuring
       e40bf89f9 hicseq [hicrep.py] - Corrected typo
       ec4b642d0 hicseq [hicrep.py] - corrected R_TOOLS path
       6ebce209e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8a7e4f113 hicseq [hicrep.py] - Corrected typo
       4ec40e07f hicseq [hicrep.py] - corrected R_TOOLS path
       d7e39f9ae [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1cc96d3b7 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5e52fc322 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       59c65ac57 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fdfbe14e [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       e5eb0973d [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e95cd4db1 Completed developing hicrep and quasar analysis
       c5c79082a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       08ff8a785 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d6abcd4d7 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a6156a17e Completed hicrep analysis except the creation of the graph
       4eea55cc9 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       43501ba2a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       25c34b435 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a164d8bf9 hicseq completed adding basic features of the hicrep analysis.
       e8aa04c33 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       491fa1ce3 Added pairwise combination for samples
       b42241ad8 started to edit the hicseq.py added new file for hicrep
       20627d696 hicseq [hicseq.py] - corrected file after rebase
       884d4e4ab hicseq [hicrep.py] - Corrected typo
       622d5f3ac hicseq [hicrep.py] - corrected R_TOOLS path
       c4a9c60ce [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fb0ec53bd hicseq [hicrep.py] - Corrected typo
       7082fd6a7 hicseq [hicrep.py] - corrected R_TOOLS path
       f90b5be46 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       02600a93f hicseq [base.ini] - updated mugqic tools version to 2.3.1
       456e6495b hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d96bfb1cc hicseq [quasar_qc - corrected module loading in matrix restructuring
       86c12350e hicseq [hicrep.py] - Corrected typo
       79fca7718 hicseq [hicrep.py] - corrected R_TOOLS path
       c28c3f88a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef796f988 hicseq [hicrep.py] - Corrected typo
       30375570a hicseq [hicrep.py] - corrected R_TOOLS path
       b5b40fe12 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e538fb0c6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       28dee2f78 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d7e7c279e [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       82aa89a0d [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       00320364a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db025271a Completed developing hicrep and quasar analysis
       edf03c154 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7b1b76054 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2520bb909 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       314a50546 Completed hicrep analysis except the creation of the graph
       71f75b6d2 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       723cecf56 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       38ec796af hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       fb02b1158 hicseq completed adding basic features of the hicrep analysis.
       56648793a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       696a52bdf Added pairwise combination for samples
       886ea4edd started to edit the hicseq.py added new file for hicrep
       10d4d24cb hicseq [base.ini] - updated mugqic tools version to 2.3.1
       9d7b8658b hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2c3b54d68 hicseq [quasar_qc - corrected module loading in matrix restructuring
       dbdf0ca64 hicseq [hicrep.py] - Corrected typo
       e70f0ee0c hicseq [hicrep.py] - corrected R_TOOLS path
       92a5667ae [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c05b038f1 hicseq [hicrep.py] - Corrected typo
       1d57da7e6 hicseq [hicrep.py] - corrected R_TOOLS path
       09b27e056 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       565b4f464 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       f9f528266 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       79da62c2b [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4518e4563 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       11e583d0c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f7d773d06 Completed developing hicrep and quasar analysis
       b7f60de8c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b3af69605 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0034cea7a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2a5492d42 Completed hicrep analysis except the creation of the graph
       d6e8b9b57 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       f1f7d98cb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ba1d617db hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       76e10ba82 hicseq completed adding basic features of the hicrep analysis.
       72a3137a0 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c2af6c676 Added pairwise combination for samples
       6f0965b26 started to edit the hicseq.py added new file for hicrep
       e9f3fe40d hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2c206ace8 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3e2a53ed3 hicseq [quasar_qc - corrected module loading in matrix restructuring
       f880a4977 hicseq [hicrep.py] - Corrected typo
       809f95e67 hicseq [hicrep.py] - corrected R_TOOLS path
       aaddec59e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       251c0c6fd hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       68889f950 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       68003e05c [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fef22083 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8782a74b4 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e2caf4499 Completed developing hicrep and quasar analysis
       30008d63a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7a640583d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       02e39e7c3 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ea8610d24 Completed hicrep analysis except the creation of the graph
       48173c514 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       419d60e00 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       b26f6236f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f2eba582 hicseq completed adding basic features of the hicrep analysis.
       92b291a50 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3dcb3ec6 Added pairwise combination for samples
       67488a04f started to edit the hicseq.py added new file for hicrep
       a22ca5679 hicseq [hicrep.py] - Corrected typo
       e28aedeef hicseq [quasar_qc.py] - Corrected deletion by mistake
       5a1f50254 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       5b6fff8da hicseq [hicrep.py] - corrected R_TOOLS path
       b02d2280c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4ed3b5932 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       c361bc620 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       562238934 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e4758a198 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c56414948 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       760bfd648 Completed developing hicrep and quasar analysis
       275a2243e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       627a50162 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       33ffee0b1 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       df084a84e Completed hicrep analysis except the creation of the graph
       f3bebd7d3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b7fe38284 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       805c4a96d hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       dbfde7643 hicseq completed adding basic features of the hicrep analysis.
       10e2dfa73 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7c2e89b82 Added pairwise combination for samples
       07099a0a3 started to edit the hicseq.py added new file for hicrep
       300ba1bb6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       688a49561 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9aa6c6145 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       0f393c68f [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       b64480f83 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       353a55d11 Completed developing hicrep and quasar analysis
       7aa1f48b8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       008904a01 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b27275b6a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       66c49d4ee Completed hicrep analysis except the creation of the graph
       a09f0387d hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5c370b552 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       de693aec8 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a504b1719 hicseq completed adding basic features of the hicrep analysis.
       fe551f3cf hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       522081626 Added pairwise combination for samples
       17e1b5c09 started to edit the hicseq.py added new file for hicrep
       8d592db3c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c37d365fb [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8c152446b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       7f55aebe3 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f371fb96b Completed developing hicrep and quasar analysis
       14b332d0e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       66ded908e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c8cabe6bd created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f066550e3 Completed hicrep analysis except the creation of the graph
       0a733f012 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       445a46848 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ebcf685ed hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8e3602c47 hicseq completed adding basic features of the hicrep analysis.
       21b08c351 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       59771c353 Added pairwise combination for samples
       981948253 started to edit the hicseq.py added new file for hicrep
       8fc473c6f [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       673bad3e5 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       979085be4 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       5392e0067 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       96b5c41e0 Completed developing hicrep and quasar analysis
       8d38fc18b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e517d3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dade053a8 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       14247d254 Completed hicrep analysis except the creation of the graph
       79e675ea9 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       59125936a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c19d1be79 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       cc2f34c1e hicseq completed adding basic features of the hicrep analysis.
       e7e729af4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8495df617 Added pairwise combination for samples
       d789632af started to edit the hicseq.py added new file for hicrep
       d8a1537b9 Completed developing hicrep and quasar analysis
       c13cb74d1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       316db2c72 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a382fb8e created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       68fc66e6d Completed hicrep analysis except the creation of the graph
       f2cb6a0bc hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       702fadb26 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4fefe1e97 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       294e3aa3f hicseq completed adding basic features of the hicrep analysis.
       5edc80e24 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       fc95853d4 Added pairwise combination for samples
       159251065 started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      10 commits

       4b08c93c3 force loading mugqic python 2.7 on cuffmerge
       ac7cd868c force loading mugqic python 2.7 on multiqc 1.7
       b0875ae99 more verbose went -d/--design is needed for contrast
       0dc71a9f2 revert on indel aligner
       37af57010 cleanup ini for beluga
       3c0fad0e6 revert on indel aligner
       54635678b cleanup ini for beluga
       2080dfdaa remove dbSNP for oxog ini block fix gatk 4 module name
       3a1ab2231 fix picard_collect_oxog_metrics dbSNP vcf
       1c5e2ed4b remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       7bbe4faf2 Corrected hicseq.py file deletions
       0ca8a123a Corrected hicseq.py file deletions
       e2e57af06 Corrected hicseq.py file deletions
       9d90aeffc Corrected hicseq.py file deletions
       312d3dbae Corrected hicseq.py file deletions
       a95d8b8a9 Corrected hicseq.py file deletions
       52c9240cb Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      11 commits

       39db7bb7d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       494327429 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       452015d9e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       61a86ecec hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       e7985498c hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       185c81b68 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21929bf11 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       fe03ec69d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cdca1f90c hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a146b2518 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4025c7165 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      11 commits

       c646e6755 refer previous commit
       30b0838ea refer previous commit
       0e4a6cb1a refer previous commit
       138fafc15 refer previous commit
       0fd5742d7 refer previous commit
       962cc3e50 refer previous commit
       bb6b1b7d3 refer previous commit
       a1f0ba284 refer previous commit
       2b79a1bd4 refer previous commit
       ae7e5c821 refer previous commit
       3b6c393a8 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       30b2d60ad corrected md file

  Shaloo Shalini <shaloo.shalini@gmail.com>      17 commits

       f0d10456f Merged in ss_issu_69_rnaseq_strintie_wf (pull request #214)
       6935d81ce Merged in ss_issu_67_wf_dnaseq_light (pull request #216)
       2b8fe8c2e Merged in ss_issu_66_dnaseq_sv_schema (pull request #215)
       17c8fcc4b Merged dev into ss_issu_69_rnaseq_strintie_wf
       0b1fe1eb3 Merged dev into ss_issu_66_dnaseq_sv_schema
       0f9def7cd Merged dev into ss_issu_67_wf_dnaseq_light
       37504fcd5 Merged in ss_issu_67_wf_dnaseq_light (pull request #212)
       b217c763e Merged in ss_issu_65_chipseq_workflow (pull request #211)
       890dc4078 Merged in ss_issu_68_nano_wf (pull request #213)
       1c7a288ca Merged dev into ss_issu_67_wf_dnaseq_light
       0a6fd2740 Merged dev into ss_issu_69_rnaseq_strintie_wf
       c7272d598 Merged dev into ss_issu_68_nano_wf
       0a1ffe39c Merged dev into ss_issu_65_chipseq_workflow
       3f9add6e5 Merged in ss_issu_66_dnaseq_sv_schema (pull request #209)
       2bddc5939 Merged in ss_issu_60_wf_ampllconseq (pull request #203)
       71461066a Merged in ss_covseq_wflow (pull request #197)
       192f7a305 Merged in covseq_readme (pull request #196)

  shaloo <shalz@hotmail.com>      30 commits

       14cb58cfa Refs #69 Hector and Ed's feedback incorporated
       ea4ea09a5 Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       de861a5d0 Refs #66 cleanup after color update
       1b6861958 Refs #67 cleanup after color update
       26733d9d2 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       75101db13 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       47c870543 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       dc52d268c Refs #67 dnaseq -t light option color feedback addressed
       0c95d465b Refs #68 color updated as per feedback
       796cf88cc Refs #68 color added as per feedback
       1ad517e60 Refs #65 chipseq workflow color updated, report links added as per feedback
       179955f4f Fixes #69 Refs #54 rnaseq -t stringtie workflow schema added
       c1a886138 Merge branch 'ss_issu_65_chipseq_workflow' of bitbucket.org:mugqic/genpipes into ss_issu_65_chipseq_workflow
       ce49d390c Fixes #65 Refs #54 chipseq pipeline updated
       8beba66d0 Fixes #68 Refs #54 nanopore pipeline workflow schema added
       cc91732d6 Fixes #67 Refs #54 dnaseq -t light workflow schema added
       9f4b8ce6e Update chipseq pipeline schema for -t chipseq and -t atacseq options to reflect latest dev branch pipeline code
       6f8845682 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       53cc1d368 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       e4fd2a28a Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       647d39a25 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       7b041d923 cleanup DS_Store file
       44db70722 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       48626df06 Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       52a228660 Refs #53 covseq workflow arrow step 6 -> 9 removed
       faf944996 Refs #53 Paul's review inputs incorporated
       51870352c Fixes #53 covseq.py workflow diagram added
       bacd6ed25 Refs #53 added covseq pipeline schema workflow draft under review by Paul
       0f7b09a65 Fixes #51 update covseq pipeline readme to v3.3.0
       2fed25b4b Fixes #51 update covseq pipeline readme to v3.3.0

3.3.0        Fri Feb 19 16:37:40 2021 -0500        641 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       a4d9a14e2 GenPipes - removing tumor_apir before release
       1a28942ac Version bump to 3.2.1-beta
       63d211df4 Version bump to 3.2.0

  ehenrion <edouard.henrion@mcgill.ca>      4 commits

       0756af3e3 Merged in release_3.3 (pull request #194)
       7b17e70bd GenPipes - Bug fix : corrected dnaseq.cedar.ini
       d62b0c014 GenPipes - Bug fix : removed buggy line in dnaseq.cedar.ini
       30bdd0a45 GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       ebef4a209 Release 3.2

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      572 commits

       e6d1d8aa6 hicseq pipeline [changed the step order]
       165e06006 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       498209baa hicseq [hicseq.py] - corrected output -o issue fastq_readset
       8d382dc57 hicseq [hicseq.py, readme.md] - modified readmefile
       638e16e1b hicseq [hicseq.py] - corrected file after rebase
       7fe877f44 hicseq [hicrep.py] - Corrected typo
       2c26336ab hicseq [hicrep.py] - corrected R_TOOLS path
       91dd29c56 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f44f83d71 hicseq [hicrep.py] - Corrected typo
       124ce0e30 hicseq [hicrep.py] - corrected R_TOOLS path
       678b9cf3a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bb10a4edf hicseq [hicrep.py] - Corrected typo
       3e2b5502c hicseq [hicrep.py] - corrected R_TOOLS path
       8772abdb8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       77cbfb905 hicseq [hicrep.py] - Corrected typo
       6f820f10f hicseq [hicrep.py] - corrected R_TOOLS path
       2528718fd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65191abb4 hicseq [hicrep.py] - Corrected typo
       28b43776a hicseq [hicrep.py] - corrected R_TOOLS path
       c82045f9f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6b1b9d621 hicseq [hicrep.py] - Corrected typo
       be461ea1d hicseq [hicrep.py] - corrected R_TOOLS path
       26cf92542 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b06d2eafa hicseq [hicrep.py] - Corrected typo
       eff718e0d hicseq [hicrep.py] - corrected R_TOOLS path
       8c7b0eb2e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d1f43e783 hicseq [hicrep.py] - Corrected typo
       97174e239 hicseq [hicrep.py] - corrected R_TOOLS path
       d82f67715 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       aa9e9215d hicseq [hicseq.py] - corrected file after rebase
       96a94204f hicseq [hicrep.py] - Corrected typo
       68f348bf5 hicseq [hicrep.py] - corrected R_TOOLS path
       284f0e1ee [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac094de0 hicseq [hicrep.py] - Corrected typo
       ec553bd9e hicseq [hicrep.py] - corrected R_TOOLS path
       e32d3bbf4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1f326ef5 hicseq [hicrep.py] - Corrected typo
       4b99cbbc8 hicseq [hicrep.py] - corrected R_TOOLS path
       da58a44bc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8985879c0 hicseq [hicrep.py] - Corrected typo
       0d1f31b6d hicseq [hicrep.py] - corrected R_TOOLS path
       d5aab5221 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1d3a959b4 hicseq [hicrep.py] - Corrected typo
       2ec58eb59 hicseq [hicrep.py] - corrected R_TOOLS path
       aaf3aee31 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2dd3de83e hicseq [hicrep.py] - Corrected typo
       0eb19568e hicseq [hicrep.py] - corrected R_TOOLS path
       c2627543c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       944ba8e7a hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4021b7d75 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       e8567612b hicseq [quasar_qc - corrected module loading in matrix restructuring
       fc92143eb hicseq [hicrep.py] - Corrected typo
       a13ed40c2 hicseq [hicrep.py] - corrected R_TOOLS path
       ee6bc35d1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3b82c00fe hicseq [hicrep.py] - Corrected typo
       9a509dcc9 hicseq [hicrep.py] - corrected R_TOOLS path
       edb308a67 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ca9fa5df8 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       d9b7dc4e1 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       aa8821b64 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       9419efd75 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8d61964ca [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       2e7870dd3 Completed developing hicrep and quasar analysis
       d54c390ab [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       3ca0a12a4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       845d4734b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       42f947ad3 Completed hicrep analysis except the creation of the graph
       0274f3227 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       1c1f819fc hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e27ded703 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       53408440d hicseq completed adding basic features of the hicrep analysis.
       5b078cfe4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       94230e0e6 Added pairwise combination for samples
       39b34b00e started to edit the hicseq.py added new file for hicrep
       9e2040f47 hicseq [hicseq.py, readme.md] - modified readmefile
       446595214 hicseq [hicseq.py] - corrected file after rebase
       ccb344f89 hicseq [hicrep.py] - Corrected typo
       35df970dd hicseq [hicrep.py] - corrected R_TOOLS path
       240f9c5db [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eecc3b296 hicseq [hicrep.py] - Corrected typo
       78e6c8ab4 hicseq [hicrep.py] - corrected R_TOOLS path
       616e8d371 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       84e8ec919 hicseq [hicrep.py] - Corrected typo
       b142de6e6 hicseq [hicrep.py] - corrected R_TOOLS path
       dd6cb48bb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5bad885b4 hicseq [hicrep.py] - Corrected typo
       030d55842 hicseq [hicrep.py] - corrected R_TOOLS path
       e29ae96c9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       97622c6d8 hicseq [hicrep.py] - Corrected typo
       a477d712d hicseq [hicrep.py] - corrected R_TOOLS path
       d05a6b625 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       12a4e2f4d hicseq [hicrep.py] - Corrected typo
       45695cbdf hicseq [hicrep.py] - corrected R_TOOLS path
       e7041d318 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       302c4e96d hicseq [base.ini] - updated mugqic tools version to 2.3.1
       154f9b278 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       8e55fb576 hicseq [quasar_qc - corrected module loading in matrix restructuring
       f1534e3b6 hicseq [hicrep.py] - Corrected typo
       45572379a hicseq [hicrep.py] - corrected R_TOOLS path
       c2063ea9a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e41bb32fb hicseq [hicrep.py] - Corrected typo
       ebbbf879a hicseq [hicrep.py] - corrected R_TOOLS path
       440631443 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bd36df09f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       40765477e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6f40129e3 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       789fb24b0 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       83ea31bd2 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       0eafabdac Completed developing hicrep and quasar analysis
       dc1cf6aa3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c0974eb5b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       82c023bbb created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       0ea49fdbc Completed hicrep analysis except the creation of the graph
       55f2b8505 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bb223680c hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       19aec5068 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       7cd078591 hicseq completed adding basic features of the hicrep analysis.
       7bdf5b494 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ae313c6d4 Added pairwise combination for samples
       08dae47d0 started to edit the hicseq.py added new file for hicrep
       eff71142f hicseq [hicseq.py] - corrected file after rebase
       b6931059f hicseq [hicrep.py] - Corrected typo
       386ca8339 hicseq [hicrep.py] - corrected R_TOOLS path
       5a412f71d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0b7fe2a23 hicseq [hicrep.py] - Corrected typo
       1765bef41 hicseq [hicrep.py] - corrected R_TOOLS path
       5851499dc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1555fed40 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       d3199ece5 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       01d239a9b hicseq [quasar_qc - corrected module loading in matrix restructuring
       46e7d3a52 hicseq [hicrep.py] - Corrected typo
       5e92f8e64 hicseq [hicrep.py] - corrected R_TOOLS path
       f2dbfd7b5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4fa81e9bb hicseq [hicrep.py] - Corrected typo
       b284455e3 hicseq [hicrep.py] - corrected R_TOOLS path
       b12fdbe51 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a78a623e3 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       125a6e305 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       67272ac2d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       19f3edc7f [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       a8b8d85c8 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9d37f638c Completed developing hicrep and quasar analysis
       acce38352 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       92bfb414c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a3bcdf29 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       812ba6d03 Completed hicrep analysis except the creation of the graph
       27467f368 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3cc9f3933 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2684371d4 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       92744be25 hicseq completed adding basic features of the hicrep analysis.
       06d8e53b0 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       31b01578e Added pairwise combination for samples
       6775791d3 started to edit the hicseq.py added new file for hicrep
       38a937ec2 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4e9dea90b hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       6efe674ba hicseq [quasar_qc - corrected module loading in matrix restructuring
       b918766ee hicseq [hicrep.py] - Corrected typo
       0219397e9 hicseq [hicrep.py] - corrected R_TOOLS path
       187cba374 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac69b20a hicseq [hicrep.py] - Corrected typo
       a80a6bcb3 hicseq [hicrep.py] - corrected R_TOOLS path
       06e2631ed [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f903c4b68 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ebc75da3e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c9acbc74b [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       dcbe50b4f [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bd63e8e6b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4f9a3087c Completed developing hicrep and quasar analysis
       6f0d1a8a5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1052c8dbe [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       704d7f190 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a7f8b2a96 Completed hicrep analysis except the creation of the graph
       f1fd9af6f hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       712cb4621 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       979429b83 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3930cd3fa hicseq completed adding basic features of the hicrep analysis.
       b3c1160aa hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       f493eae58 Added pairwise combination for samples
       6f98cb3a1 started to edit the hicseq.py added new file for hicrep
       f58f27f5d hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b576afd01 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       0eb37c3f2 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b718a5b2c hicseq [hicrep.py] - Corrected typo
       5f8cf4740 hicseq [hicrep.py] - corrected R_TOOLS path
       1fde1089a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6478d0797 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b262df5bf [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b340c0d5c [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b6b35788 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       d740ed8ab [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c1b8eef75 Completed developing hicrep and quasar analysis
       6bc89a861 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       95c477e89 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       53fc03ce8 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       cf6a3af07 Completed hicrep analysis except the creation of the graph
       d8542242e hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b961772e hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       8438ff65d hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b18b1915c hicseq completed adding basic features of the hicrep analysis.
       c4ba0941e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bd650a8bf Added pairwise combination for samples
       6f4136c3d started to edit the hicseq.py added new file for hicrep
       9364f1e4a hicseq [hicrep.py] - Corrected typo
       301afce04 hicseq [quasar_qc.py] - Corrected deletion by mistake
       bd2a8e128 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       96761788d hicseq [hicrep.py] - corrected R_TOOLS path
       a5ca0babb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       432f5862f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       3e68e18c1 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5c841f744 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       ab730be1b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       3bce8239d [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e1cee2e16 Completed developing hicrep and quasar analysis
       45b526c8d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6cf845071 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       21e2e6ce1 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       34416892a Completed hicrep analysis except the creation of the graph
       6570bc973 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       46a28b493 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       645c2f89a hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       0189b8082 hicseq completed adding basic features of the hicrep analysis.
       260695138 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       6fd22e84d Added pairwise combination for samples
       d664e0e54 started to edit the hicseq.py added new file for hicrep
       1e51cf6de hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       4cca0610c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2f8c42d4d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       1c49e4580 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cfa04f301 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c0bd36519 Completed developing hicrep and quasar analysis
       70d8bdfe8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6711a6899 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e8effea3c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       5fe98a017 Completed hicrep analysis except the creation of the graph
       0bc9a5dd2 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bea29b638 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       eb89bc046 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       d7bb65dff hicseq completed adding basic features of the hicrep analysis.
       ccc376e43 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       93f957d35 Added pairwise combination for samples
       1c14ddf21 started to edit the hicseq.py added new file for hicrep
       0da3f4c4c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9ae9d5eda [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       30a4e554b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8a5ac66bb [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4fb194fae Completed developing hicrep and quasar analysis
       7f2917d21 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2492118cf [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9a8cfcd53 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       54c7c37f8 Completed hicrep analysis except the creation of the graph
       c35566ccc hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fe1171a03 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e9db27ddf hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8584e29d6 hicseq completed adding basic features of the hicrep analysis.
       dccac155c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       568905105 Added pairwise combination for samples
       f609b2f9f started to edit the hicseq.py added new file for hicrep
       31ba309e9 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       55f775882 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f794061c5 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       84e93b1bf [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       915d9ace6 Completed developing hicrep and quasar analysis
       769c97be1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       01ed5378a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b8a79f4f1 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       616324014 Completed hicrep analysis except the creation of the graph
       7efc2ba3e hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bf7412fdb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       14256a142 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       042082e58 hicseq completed adding basic features of the hicrep analysis.
       fa1e7d0a2 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       eba9e4687 Added pairwise combination for samples
       9bd17ff30 started to edit the hicseq.py added new file for hicrep
       99d622283 Completed developing hicrep and quasar analysis
       719115c36 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7c3222a92 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       966e814b2 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ece06e332 Completed hicrep analysis except the creation of the graph
       4ceda720c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a650327a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       91f3b4b76 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       15e34d054 hicseq completed adding basic features of the hicrep analysis.
       92baadc30 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d9d98a541 Added pairwise combination for samples
       c4853c36e started to edit the hicseq.py added new file for hicrep
       78ab01055 hicseq pipeline [changed the step order]
       c2ac146d2 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       8bf6a4d7a hicseq [hicseq.py] - corrected output -o issue fastq_readset
       a016c1ebd hicseq [hicseq.py, readme.md] - modified readmefile
       7fbb757c7 hicseq [hicseq.py] - corrected file after rebase
       f8fc53cba hicseq [hicrep.py] - Corrected typo
       46c6016a9 hicseq [hicrep.py] - corrected R_TOOLS path
       fb83a9249 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       00a5b87fa hicseq [hicrep.py] - Corrected typo
       3e840af9e hicseq [hicrep.py] - corrected R_TOOLS path
       c61122bc6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d3659cb4c hicseq [hicrep.py] - Corrected typo
       230b7b68d hicseq [hicrep.py] - corrected R_TOOLS path
       9a454ab89 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       857f8c9a7 hicseq [hicrep.py] - Corrected typo
       5c4024507 hicseq [hicrep.py] - corrected R_TOOLS path
       ba545bd4b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3ff1e5094 hicseq [hicrep.py] - Corrected typo
       a0098d0ef hicseq [hicrep.py] - corrected R_TOOLS path
       2c820d824 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d23933c71 hicseq [hicrep.py] - Corrected typo
       b06456439 hicseq [hicrep.py] - corrected R_TOOLS path
       54bb894c4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b81d75509 hicseq [hicrep.py] - Corrected typo
       74e48a258 hicseq [hicrep.py] - corrected R_TOOLS path
       40fe0ccb4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a3680d650 hicseq [hicrep.py] - Corrected typo
       798c19f6f hicseq [hicrep.py] - corrected R_TOOLS path
       06ec78860 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2fd36b1d5 hicseq [hicseq.py] - corrected file after rebase
       be817b445 hicseq [hicrep.py] - Corrected typo
       5248d3250 hicseq [hicrep.py] - corrected R_TOOLS path
       39617b231 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec9796015 hicseq [hicrep.py] - Corrected typo
       24e7c7c96 hicseq [hicrep.py] - corrected R_TOOLS path
       4d6d3f889 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c17d5f120 hicseq [hicrep.py] - Corrected typo
       aeae41c10 hicseq [hicrep.py] - corrected R_TOOLS path
       77dabaf91 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       032cfe90a hicseq [hicrep.py] - Corrected typo
       07a72bd18 hicseq [hicrep.py] - corrected R_TOOLS path
       1039f672f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       590a6004d hicseq [hicrep.py] - Corrected typo
       8f9b54667 hicseq [hicrep.py] - corrected R_TOOLS path
       574336c63 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6c20f4925 hicseq [hicrep.py] - Corrected typo
       ff45f2df3 hicseq [hicrep.py] - corrected R_TOOLS path
       85d661e34 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec67234d4 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       463da45c8 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       7701580c5 hicseq [quasar_qc - corrected module loading in matrix restructuring
       0c5ee2fb7 hicseq [hicrep.py] - Corrected typo
       1a221b482 hicseq [hicrep.py] - corrected R_TOOLS path
       8966cfbcb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b02d1a7f4 hicseq [hicrep.py] - Corrected typo
       4d37159f4 hicseq [hicrep.py] - corrected R_TOOLS path
       22fbd6065 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5541fbc28 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       fdeb31d16 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d58d40f22 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       fa1f4b49c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cdcc79068 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7b9f6466b Completed developing hicrep and quasar analysis
       31c4df576 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c1cc394cf [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       80dc0ab59 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       72cbd0039 Completed hicrep analysis except the creation of the graph
       0a938267c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       85cc99e50 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2e455782a hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       46f2fa315 hicseq completed adding basic features of the hicrep analysis.
       8c3e8d90d hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8e5c6545d Added pairwise combination for samples
       55862cd88 started to edit the hicseq.py added new file for hicrep
       e79b132de hicseq [hicseq.py, readme.md] - modified readmefile
       ffea320a9 hicseq [hicseq.py] - corrected file after rebase
       fed625353 hicseq [hicrep.py] - Corrected typo
       d285ff0e6 hicseq [hicrep.py] - corrected R_TOOLS path
       e04cb3926 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       de2d8e170 hicseq [hicrep.py] - Corrected typo
       9c6f2d219 hicseq [hicrep.py] - corrected R_TOOLS path
       a8ae3ab42 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       866b11d85 hicseq [hicrep.py] - Corrected typo
       82db1223c hicseq [hicrep.py] - corrected R_TOOLS path
       4b9e5e65c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       dda535240 hicseq [hicrep.py] - Corrected typo
       452f1db6d hicseq [hicrep.py] - corrected R_TOOLS path
       db683cc7a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f2e361292 hicseq [hicrep.py] - Corrected typo
       2accb63c3 hicseq [hicrep.py] - corrected R_TOOLS path
       f8413f07f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7938420d5 hicseq [hicrep.py] - Corrected typo
       e78e318b7 hicseq [hicrep.py] - corrected R_TOOLS path
       80eff2042 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0b006288 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       1160fe259 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3ba117420 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b27f4021a hicseq [hicrep.py] - Corrected typo
       01ac15dab hicseq [hicrep.py] - corrected R_TOOLS path
       80ad7082a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef6d0348e hicseq [hicrep.py] - Corrected typo
       045b540a4 hicseq [hicrep.py] - corrected R_TOOLS path
       f6422cd3b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8808f0558 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       65f84c4ef [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7adcfa1eb [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       14c0a7896 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       935862a90 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       fd6eef8ac Completed developing hicrep and quasar analysis
       8c8248e46 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       85f0d9efc [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dc225b23a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       dfe7725f2 Completed hicrep analysis except the creation of the graph
       db50a5d97 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       10222fe0f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       5dcda07ce hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       062551bf2 hicseq completed adding basic features of the hicrep analysis.
       766c1b493 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e979eee43 Added pairwise combination for samples
       cfff9164f started to edit the hicseq.py added new file for hicrep
       4443e3d2e hicseq [hicseq.py] - corrected file after rebase
       ca8b09dbe hicseq [hicrep.py] - Corrected typo
       584dd1847 hicseq [hicrep.py] - corrected R_TOOLS path
       afa77dc2c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c2535f69 hicseq [hicrep.py] - Corrected typo
       8db0da7e2 hicseq [hicrep.py] - corrected R_TOOLS path
       6f264bd60 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5960be4ca hicseq [base.ini] - updated mugqic tools version to 2.3.1
       537be8296 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       1313e2b78 hicseq [quasar_qc - corrected module loading in matrix restructuring
       9bcc430cc hicseq [hicrep.py] - Corrected typo
       8cacfa8a7 hicseq [hicrep.py] - corrected R_TOOLS path
       c8c3cef02 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71b02a39e hicseq [hicrep.py] - Corrected typo
       b98db26dd hicseq [hicrep.py] - corrected R_TOOLS path
       3698d9506 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d461e1e88 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5c728c649 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7a9cfa6b0 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2b665729d [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       312c08eac [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7e6f9aa02 Completed developing hicrep and quasar analysis
       334abafc6 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       87973bb4b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       03aab022c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8359e46d5 Completed hicrep analysis except the creation of the graph
       8a845b7b4 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       eee8a212c hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       9a55b5c2e hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2755cba76 hicseq completed adding basic features of the hicrep analysis.
       3c517c462 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8852df1ff Added pairwise combination for samples
       54e5f9fdd started to edit the hicseq.py added new file for hicrep
       9d5410099 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       95e4b22fc hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d63f7b4d9 hicseq [quasar_qc - corrected module loading in matrix restructuring
       56c0d8c42 hicseq [hicrep.py] - Corrected typo
       80acf5283 hicseq [hicrep.py] - corrected R_TOOLS path
       7cf9e9469 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65bafda3f hicseq [hicrep.py] - Corrected typo
       990b8fae3 hicseq [hicrep.py] - corrected R_TOOLS path
       b86b171e8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e5fdc3a9b hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff4951167 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5da9af436 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8113e902a [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       845ec6267 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       09a27509a Completed developing hicrep and quasar analysis
       1e1ab4c32 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4f0cba114 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6ba3f9dd0 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       be91a4064 Completed hicrep analysis except the creation of the graph
       e2189f7c6 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fc54a3a23 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       09ec9842f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       c58406d48 hicseq completed adding basic features of the hicrep analysis.
       36dd6fab9 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       1861910ad Added pairwise combination for samples
       ef297a339 started to edit the hicseq.py added new file for hicrep
       5d7c80f9a hicseq [base.ini] - updated mugqic tools version to 2.3.1
       868a71108 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9685043dd hicseq [quasar_qc - corrected module loading in matrix restructuring
       f5f88d63d Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu Corrected conflicts Conflicts: 	bfx/quasar_qc.py 	pipelines/hicseq/hicseq.base.ini 	pipelines/hicseq/hicseq.py
       5a8386028 hicseq [hicrep.py] - Corrected typo
       b9e6c3b87 hicseq [hicrep.py] - corrected R_TOOLS path
       6d4ce4f3c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b89f8e745 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       8beb59079 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b78f92c16 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       18523b876 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       0ff0ff7fa [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c6f98d60f Completed developing hicrep and quasar analysis
       0b10bc0b1 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       8b60a2bc6 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6f50b29bc created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8faa8ebf5 Completed hicrep analysis except the creation of the graph
       0bd1d0a36 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b6ab90b2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4a7bbcaf9 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       74e0bcb9d hicseq completed adding basic features of the hicrep analysis.
       9cce00606 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       183d5a955 Added pairwise combination for samples
       7f7301fdd started to edit the hicseq.py added new file for hicrep
       0f33a580e hicseq [hicrep.py] - Corrected typo
       c47a6eec4 hicseq [quasar_qc.py] - Corrected deletion by mistake
       8b33c6fed hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       61be1c4db hicseq [hicrep.py] - corrected R_TOOLS path
       20c3b4a87 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       edda3a126 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       18506b49e hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       75cbaaa96 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2a0e1a8a1 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       15bbbb377 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       44735d2f4 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       eccd7c756 Completed developing hicrep and quasar analysis
       e9078d817 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2ac9239fe [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       542bad0b2 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       418c51563 Completed hicrep analysis except the creation of the graph
       f587c7427 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a9980d36 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2d8702023 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8df5abcd7 hicseq completed adding basic features of the hicrep analysis.
       f55e1969c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       534d63fdc Added pairwise combination for samples
       756599e17 started to edit the hicseq.py added new file for hicrep
       b45e44d66 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       cf5a2f409 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       56985bb0d [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       b48157cd0 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bdbc7de85 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       d1c1ae9b5 Completed developing hicrep and quasar analysis
       e68d3079b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c832c4911 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1b63533a6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2c192d1d0 Completed hicrep analysis except the creation of the graph
       e081f1727 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       6f3560e48 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       35efcde1e hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5d6fcc133 hicseq completed adding basic features of the hicrep analysis.
       ca0d1635e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       721653183 Added pairwise combination for samples
       396a06d96 started to edit the hicseq.py added new file for hicrep
       35412d6c4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97ce47465 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       38b339b39 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       53603e16c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db3da4bf9 Completed developing hicrep and quasar analysis
       ba7bbb879 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e6625b9bc [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       a7978565e created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       720a6263d Completed hicrep analysis except the creation of the graph
       b5cbb3fe3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bfbb89683 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       1099ed048 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5f5f1685f hicseq completed adding basic features of the hicrep analysis.
       f30316307 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bbbcf3ef4 Added pairwise combination for samples
       a5c29b574 started to edit the hicseq.py added new file for hicrep
       292a0b2b2 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       cd9d16ade [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       83f719b4a [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       4c23bf13b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       5f42b4c2c Completed developing hicrep and quasar analysis
       849b382b5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d9ed64810 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1767a32bb created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       068daa9dc Completed hicrep analysis except the creation of the graph
       210ad2912 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       411d8e10b hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0f31f813c hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       79aff1763 hicseq completed adding basic features of the hicrep analysis.
       b5075cb22 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       57e8baef6 Added pairwise combination for samples
       4a57c1a75 started to edit the hicseq.py added new file for hicrep
       89006b6c3 Completed developing hicrep and quasar analysis
       87eb470f0 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       477281afa [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       10c5eb187 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       335f4c528 Completed hicrep analysis except the creation of the graph
       10daca130 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       4aeab549c hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       509fef1cc hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f096d180 hicseq completed adding basic features of the hicrep analysis.
       b623db54b hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c6057427c Added pairwise combination for samples
       9ad63fde3 started to edit the hicseq.py added new file for hicrep

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      16 commits

       a1ddafb1d Corrected hicseq.py file deletions
       e011fbe49 Corrected hicseq.py file deletions
       9c40c0373 Corrected hicseq.py file deletions
       d8ded4481 Corrected hicseq.py file deletions
       0a0903c52 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       1fdc6405d Corrected hicseq.py file deletions
       cbf6c692e Corrected hicseq.py file deletions
       96e5d3bd6 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       41a5bc296 Corrected hicseq.py file deletions
       3ee012f38 Corrected hicseq.py file deletions
       228a9e46f Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       9211f3c9b Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       3f9182510 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected merge conflict after rebase Conflicts: 	pipelines/hicseq/hicseq.py
       6291f8980 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected module loading
       296b86bf9 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a757b8dfc Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      20 commits

       925b082a7 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3d8131d05 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ae8f296e0 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21e766523 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       938c27c9b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       38516f39a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       08b138116 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4f835f9ba hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       402a4d1ae hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       6b60d1290 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       dd98639a0 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2d3a5561d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cd42f4e54 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       26313486b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       996ba98bb hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       d384dc9b0 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       5ce6c34db hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       bfd60bcda hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       c743b4459 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2305332a1 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudu.nawarathna@mail.mcgill.ca <pnawarat@abacus3.ferrier.genome.mcgill.ca>      1 commits

       e5fa3455e Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      20 commits

       1608e2b11 refer previous commit
       a749a350a refer previous commit
       ce31a37f2 refer previous commit
       4e8a71e44 refer previous commit
       2ba43453b refer previous commit
       5c06e4346 refer previous commit
       bcac09899 refer previous commit
       1f8553540 refer previous commit
       04967444e refer previous commit
       f0eb790eb refer previous commit
       afce94504 refer previous commit
       fb33ca902 refer previous commit
       459844518 refer previous commit
       ba1f6b0f9 refer previous commit
       9227c2a40 refer previous commit
       3b5869f79 refer previous commit
       daccd0792 refer previous commit
       34642e8dd refer previous commit
       4bd83cecb refer previous commit
       4ccc10640 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       3c7dbddcd Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a9fe3dd44 corrected md file
       fc8ac1ade corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       1114ea8a5 Merged in hicseq_hicrep_pubudu (pull request #170)

3.2.0        Mon Jan 25 12:47:42 2021 -0500        1147 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      19 commits

       ed04f3308 testing end-of-line character
       957c11dbf Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       057a8d130 GenPipes - Genomes : updated R installation script to ease installation in dev, also corrected the hicup module called in install_genome.sh
       9e07fc4fa GenPipes - HiC-Seq pipeline : corrected fastq_readName_edit input path
       dd582121b GenPipes - DNA-Seq pipeline : correcting symlink creation in sambamba_merge_sam_file
       7e00f1d3a Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1799cea72 updated version of MUGQIC_TOOLS in installation script
       fb8ad6e09 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7a9b872ee GenPipes bug correction - corrected ln() function in bfx/bash_cmd.py
       62402add6 GenPipes Update - debugging use portal_output_dir variable : check for both undef and empty value
       2830d2351 GenPipes Genome - added Sus_scrofa.sh (Pig genome) installation script
       81d580adb GenPipes Soft - added kent.sh installation script
       756d9a0ef Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       076e2008f GenPipes Update - DNASeq - updated sym_link for better handling of path reconstruction
       19cc2c06a GenPipes Update : fixing path of config files when passed to job2json script
       690f56e12 Version bump to 3.1.6-beta
       369cad4a6 Version bump to 3.1.5
       1d204f91f DNASeq - removed use of 'os.path.abspath' in call of 'ln()'
       399aac0be DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable

  Édouard Henrion <henrione@beluga3.int.ets1.calculquebec.ca>      4 commits

       46405691e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       8c465415c GenPipes - ChIPSeq update : a bit of code cleaning and simplifying
       b03ee2f66 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       72c753810 GENOME INSTALLATION - updated install genome.sh : added bismark genome preparation + refined genome digest command

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      2 commits

       0dc068bb3 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       c7449389d Merged in eh_quick_fixes (pull request #144)

  Édouard Henrion <henrione@ip18.m>      6 commits

       b19092652 GenPIpes Update - corrected one problematic sym_link call...
       ed444a743 GenPipes Update - corrected pipeline behavior regarding PORTAL_OUTPUT_DIR environment variable : if te variable is empty or not set, then no JSON file at all will be generated
       c54075179 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6331eed2b C3G Software - added demuxlet.sh installation script
       d63b35682 Genome Update - install_genome.sh
       a899066e1 Genome Update - some updates to Homo_sapiens.GRCh38.sh

  ehenrion <edouard.henrion@mcgill.ca>      20 commits

       1ce85b320 GenPipes : illumina_run_processing.py : correcting indentation
       f2d7e7240 GenPipes - Tumor Pair pipeline : removed dev vawk module in tumor_pair.base.ini
       62ac67971 GenPipes - DNASeq pipeline - removing mugqic_dev modules in gatk4.ini
       39ceaa259 GenPipes - Tumor Pair pipeline : removing mugqic_dev modules
       7d9f81141 BFX- Software - iVar installation script update with latest version
       3c3a80c29 common.py : genpipes replacing mugqic_pipeline....
       8658f7b5a GenPipes RNA-SEq - calling DESeq2 instead of DESeq for differential_expression
       6a0c82559 Coorected typo in README.md
       a84f5aa2f Merged in genpipes-3.1.5_release (pull request #166)
       d59cb5940 GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       af972e150 GenPipes - pipelines/dnaseq.py  : corrected prefix generation in SymLinkFastq step
       8ebbd14e5 GenPipes - pipelines/common.py corrected outputs name generation patterns for SamToFastq & Trimmomatic steps
       65fe80a3a Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918129493 (pull request #142)
       97afb33b8 GenPipes - dnaseq.py : bug correction - typo removed
       4f9a8dd76 GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       fc67a7bc4 GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       9fbd085d0 GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       16a5635bc GenPipes - nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       853c806dc updated methylseq.base.ini, useless comments removed
       be965bdc1 updated methylseq.base.ini, useless comments removed

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      27 commits

       d48572880 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       7cb402439 final edits to the nanopore pipeline
       f1db8a944 Added support for gzipped fastq
       8da4a1c5d Added the nanopore CIT ini file
       86988dfd1 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       71f1658cc Final corrections before merge to dev
       6994d6932 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       83a6f5f70 More corrections to the INI files for other servers to allow proper running with SLURM
       caf6cb2e4 Corrected an error in the mp2b ini file
       e2637c0d4 Made corrections to the nanopore ini files for all servers that were causing the pipeline to break when running with SLURM
       a10fe05c2 Merge conflict resolutions
       e08cd6d7b Corrected problem with nanopore readset parsing that caused problem with the paths
       5d884c133 Added full documentation for Nanopore pipeline, as well as the Graham config file.
       e2edb7042 Included commands necessary to add readgroup tag to alignments in minimap2
       278feb744 Final corrections before testing on other servers
       cc0f765a1 Corrected merge error related to .gitignore
       cb8378851 Fixed bug caused by missing module import in nanopore.py
       649e04b23 Added minimap2 script that was missing from previous commit
       0f2a89ddd First working version of the nanopore pipeline
       7acad993c Added full documentation for Nanopore pipeline, as well as the Graham config file.
       7c528bd8b Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       56699ff4b Included commands necessary to add readgroup tag to alignments in minimap2
       3d53188da Final corrections before testing on other servers
       558f3e9b3 More bug corrections for nanopore pipeline after initial testing. Switched to only one protocol, with an optional first step (guppy)
       5107ed7c2 Fixed bug caused by missing module import in nanopore.py
       653c4e08c Added minimap2 script that was missing from previous commit
       a0fe5ebd7 First working version of the nanopore pipeline

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      8 commits

       50bc9f3a3 Merged in Jose-Hector-Galvez/rnaseq_lightbaseini-edited-online-with-b-1580485656011 (pull request #171)
       e2c91c8b2 Added module_perl to the rnaseq_light ini file.
       007835bd8 Merged in Jose-Hector-Galvez/jobpy-edited-online-with-bitbucket-1576266027164 (pull request #154)
       c500f6b4e Suggestion, add `module purge` to all jobs that load modules, to avoid conflicts between modules.
       4c35a7c85 Merged in Jose-Hector-Galvez/gatk4py-edited-online-with-bitbucket-1575403540009 (pull request #135)
       89e70fe06 gatk4.py edited to correct for inconsistencies with configuration parameters within functions.
       eefc95aa2 Merged in Jose-Hector-Galvez/found-a-bug-in-the-schedulerpy-script-i--1575322351581 (pull request #127)
       afae6eecf Found a bug in the scheduler.py script. I am adding a line to correct it.

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      1 commits

       3ec3f1c76 Added CoVSeQ ini files for Graham and Cedar. Corrected a few errors on the Beluga ini

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      6 commits

       42981f133 Merged in rnaseq_light_docs (pull request #183)
       a190faaaa rnaseq_light.py edited to adjust docstrings to address issue raised by Shaloo here : https://github.com/c3g/GenPipes/issues/63
       ceb0c8bbb Merged in JoseHectorGalvezLopez/nanoporebaseini-edited-online-with-bitbu-1596652406491 (pull request #179)
       35bc91bb9 Edited the nanopore ini file to address the mugqic_tools error.
       10e762bb8 Merged in mgi_stretenp (pull request #177)
       396071701 Merged in nanopore_jhg (pull request #173)

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       52f9015b8 modify cedar ini
       f38492f7b a bit ugly resolution from argparse overriding issue of the type argument

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      798 commits

       41671c605 Merged in mgi_stretenp (pull request #192)
       63a546832 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       a3cffe0a1 Changing qualimap bamqc section naming to a parameter in bfx
       6e33047a2 Changing qualimap bamqc section naming to a parameter in bfx
       ad67e66e3 Merged in mgi_stretenp (pull request #191)
       356bd208f Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       dcbf891d2 igvtools ressources change
       f5026a89e Reducing ressources
       f15a4f8ef Fixing awk
       c74a4c670 Changing default ressources
       fede180dc igvtools ressources change
       58c298c2a Reducing ressources
       29ec475f6 Fixing awk
       cfca926e9 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       66189928e Changing default ressources
       9f663b00e Changing default ressources
       e74c221f7 Merged in mgi_stretenp (pull request #188)
       21f213f09 Fix
       cec6282b9 Fix
       907ab2a00 Fix interval_list checking
       08c1c5493 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       9ba271179 Checking if interval_list is already on genome folder
       b6f2b9195 Merged dev into mgi_stretenp
       247697b67 samtools bam2fq typo
       92d206f9f Adding kraken to beluga ini
       bb375b9b4 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       43d126ebb Switching to ivar 1.3 and NML dehosting filters
       7b86a544b Updating ivar trim usage
       d04bd81ef Adding htslib to ivar module
       3ea8a2e58 Fixing kraken output and adding warning on rnaseq htseq-count
       9f4fde94c Switching to latest kraken2 db build
       cad6ed3e7 Forcing pigz to overwrite
       522905b7d Switching to latest kraken2 release
       de8e702b1 Fixing kraken
       e2249ccd3 Fixing kraken bfx
       28dcdf305 Fixing kraken bfx
       63339cd2c Fixing kraken bfx
       0f47d1e74 Fixing kraken bfx
       338b71d5d Fix kraken
       3f5c373c8 Fix kraken
       d13c20d0b Addinf kraken analysis for metrics
       272131bb0 Update sambamba sort ini
       694510458 Fix rename ln test
       940b62c64 Fix
       3ac107c79 Fixing tee
       624b2c23e Test
       c469a365e Fix
       9993fca82 Adding metrics on dehosted bam
       ff29ecf68 renaming cit covseq file
       01e6e6c23 cit ini for cit test
       67b78d6ac Fixing hybrid genome path and output selection
       5d8d2f7e3 renaming cit covseq file
       15ba0d218 cit ini for cit test
       f6f9782c2 Switching to ivar 1.3 and NML dehosting filters
       4e16c7252 Updating ivar trim usage
       4a07fed8d Adding htslib to ivar module
       05475b4a0 Fixing kraken output and adding warning on rnaseq htseq-count
       82ba08478 Switching to latest kraken2 db build
       05e7d80cf Forcing pigz to overwrite
       a11ec96c9 Switching to latest kraken2 release
       e43b7fb3c Fixing kraken
       2f5805b76 Fixing kraken bfx
       879ee559a Fixing kraken bfx
       19b050c96 Fixing kraken bfx
       2cb881820 Fixing kraken bfx
       eb614ee56 Fix kraken
       74227ccaa Fix kraken
       c797c2aa6 Addinf kraken analysis for metrics
       52b1110fe Update sambamba sort ini
       a06f58862 Fix rename ln test
       398fcefcf Fix
       012b267d8 Fixing tee
       54b5fa5b3 Test
       33201b803 Fix
       516add157 Adding metrics on dehosted bam
       d5b321a21 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       0b9b56a0d Reducing sambamba filtering default cpus
       77cddeb15 Fixing rename consensus
       916b0912f renaming cit covseq file
       a5509539c cit ini for cit test
       67a817420 Changing ini sections to use Illumina beds files as default
       fcfce9971 Changing default seq_method
       84afdaa0e Cleaning and switching to cvmfs genome
       ed0fc1967 Fix cit
       409ffcd2b Cit fix
       ccd50ffb2 Cit test
       c999d20c2 Fix
       24fd48c37 Fix
       a11172c3c Fix
       03180fbb7 Fix rename consensus symlink
       7533de7b3 Fixing rename consensus
       d33ab92b1 Fixing tsv
       b09b78e0d Fixing tsv
       2ab5bd4a7 Fixing tsv renaming
       539548990 Fixing output rename header + tsv for ncov-tools
       b0eb526d9 Fixing picard metrics
       13a5d3103 Collecting picard metrics on raw AND filtered bam
       95dadb9ba Fixing select output file
       8cfc8ea8b Fixing hybrid genome path and output selection
       35c8d4925 Fixing merging step for 1 sample with multiple readsets
       210ebdd09 Update inis
       fe84c0b02 quast -> Quast
       13c60bf9e renaming cit covseq file
       9704c198c cit ini for cit test
       11fd8090b Reducing sambamba filtering default cpus
       5843d7bba Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       1fccae9af Changing ini sections to use Illumina beds files as default
       4c321cab4 Changing default seq_method
       ad84e0756 Cleaning and switching to cvmfs genome
       5b66a1b33 Fix cit
       8ac2d785f Cit fix
       aa5609b0d Cit test
       a5a46a530 Fix
       360e7a1c3 Fix
       46bc3964a Fix
       93fa45253 Fix rename consensus symlink
       c48a29892 Fixing rename consensus
       816eb5798 Fixing tsv
       24e973434 Fixing tsv
       679da0c2b Fixing tsv renaming
       52b8e1b80 Fixing output rename header + tsv for ncov-tools
       8727d67a5 Fixing picard metrics
       c6a1a2d1b Collecting picard metrics on raw AND filtered bam
       860d44bf9 Fixing select output file
       26be1581a Fixing hybrid genome path and output selection
       facc8cba4 Fixing merging step for 1 sample with multiple readsets
       1442f6a17 Update inis
       126d53dd3 quast -> Quast
       3a7504d56 renaming cit covseq file
       efac275b4 cit ini for cit test
       412da1300 Changing ini sections to use Illumina beds files as default
       24ccf1f98 Changing default seq_method
       f4d97e9c8 Cleaning and switching to cvmfs genome
       3c8daf00f Fix cit
       e6e973853 Cit fix
       3ba059d91 Cit test
       7a29d012c Fix
       e6a868652 Fix
       e179616bc Fix
       e528c792e Fix rename consensus symlink
       4de707684 Fixing rename consensus
       543016938 Fixing tsv
       2952b1ace Fixing tsv
       4ca82a109 Fixing tsv renaming
       f7df72ad7 Fixing output rename header + tsv for ncov-tools
       00b605227 Fixing picard metrics
       c009672a2 Collecting picard metrics on raw AND filtered bam
       fc6322f34 Fixing select output file
       a33b34c53 Fixing hybrid genome path and output selection
       75332a91b Fixing merging step for 1 sample with multiple readsets
       4a4b67eb7 Update inis
       527a801a0 quast -> Quast
       6bc79c0a4 renaming cit covseq file
       d467187ee cit ini for cit test
       41514c00d Merged in mgi_stretenp (pull request #175)
       d6ae1290f Adding WARN for not changing R ver in rnaseq denovo
       f51f69388 Fixing rnaseq denovo R issue by creating deseq and a deseq2
       4a1a5bff7 Fixing back rnaseq denovo
       56fb3fbfb Fixing rnaseq denovo
       1d29eb6d9 Fixing rnaseq denovo assembly R versions
       50fb4704c Fixing trinity report on cit
       1ae8b1773 deliverables fix
       c51f7e18d Switching to 10% as minor variants threshold
       dca8feac8 Removing kraken module & Changing R version for rnaseq
       c49a0733c Including sambamba changes in ini
       29b5aec32 Including sambamba changes into other ini
       c4a7ef31c Fixing gatk_indel_realigner if 1 job
       8fd747f1d Including sambamba modifs in ini
       d826d8198 sambamba merge realigned high coverage
       f135bafcf Switching to sambamba merge sam
       c55297c32 Including bwa sambamba into high coverage
       466906d18 Inluding a with sambamba bam.bai to picard mark duplicates
       858c1f539 Merge branch 'dnaseq_picard_fix' into mgi_stretenp
       a3354a526 Typo
       dd39deddd Fixing for cit run
       1b4d01494 Merge branch 'dev' into mgi_stretenp
       d9566ff33 Using Illumina as default
       34cd42234 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       2502abac3 Renaming outputs and default genome
       8ccb242e5 Fixing consensus renaming
       075ab5a6a Fixing consensus renaming
       9dec66856 Renaming + flaging consensus seq
       621449abf Fixing quast options
       da2638473 Fixing qualimap output
       cb98c5efc Fixing outputs quast and snpeff
       6d220605d Fixing mkdir
       d15425c80 Fixing typo
       921d5657f Fixing typo
       53c68fdda Changing to samtools bam2fq and fixing options ini
       4f5c5ceb1 Fixing quast
       f7042ec7c Fixing latest commit
       04f291650 Adding intermediate bam file on host removal step
       d931a8938 Forcing pigz to overwritte gz
       b0711a8e6 fixing host removal
       4cfa142e3 Adding pigz within bam2fq to be able to skip step
       7916935ba quast fix
       4788baed4 quast + snpeff + host removal fix
       2e3bc597a Fixing cutadapt input files
       75d63826f Fixing type
       308ace13f Removing indexing after name sorting
       f76fdadcf Fixing path creation at host removal step
       d5aa6f0e0 Fixing snpeff
       92024780b Removing print files
       afa134ef8 Fixing input choosing
       979eafe25 Fixing choosing inout files mapping
       9b6ae435c Fixing host removal
       f7ba72545 Fixing host removal
       f0ac000f0 Fixing quast step
       8ce9bf21b Fixing quast
       c07ffae73 Fixing quast step
       5415ab72f Fixing param requirements
       9218a8e71 Fixing typo
       52999aeec Fixing pigz
       fec2ceae5 Not using kraken anymore
       3363c61c7 Switching steps order
       7e51be51c Adding 3 steps
       863d61b11 Fixing picard multiple metrics raw
       4d8cc43bc Insert Size metrics on raw bam
       f2e8b42b1 sambamba flagstat fix
       0a8c640fc Renaming snpeff step/job
       d47340446 Fixing flagstat
       8163c9a31 Fixing flagstat
       8461e1983 flagstat on all the bams
       a380f6816 Zipping output of snpeff
       28969ba9d Fixing snpeff
       a899bdd57 Fixing cutadapt
       2f86db554 Flagstat on raw bam
       f2b8a9f0b Fix
       f4fba5a7f Fixing renaming step
       b0fcd1cb4 Switching to fgbio trim primer
       e9c450f8f Updating metrics
       6582bc9b7 Addition of consensus step
       319e50d62 MGI init commit
       6df0cca45 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       350131aa0 DNA-Seq - Fix update
       5cce81872 DNA-Seq - Fix update
       838cc9208 DNA-Seq - Fix update
       990fc5644 DNA-Seq - Fix update
       7e7b81073 DNA-Seq - Fix update
       0a8d7eee3 DNA-Seq - Fix update
       fe22929a3 DNA-Seq - Fixmate with samtools and sorting with sambamba
       3c84dd1e8 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       e0fa9ff8c DNA-Seq - Fix update including input_file_revision_eh changes
       d03669f11 DNA-Seq - Fix update
       5f22222dc DNA-Seq - Fix update
       18dd10693 DNA-Seq - Fix update
       08c2864ba DNA-Seq - Fix update
       117181c22 DNA-Seq - Fixmate with samtools and sorting with sambamba
       813c02e60 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c193068cb DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       0e8d422c8 DNA-Seq - Fix update including input_file_revision_eh changes
       80f732a49 Zipping and indexing vcfs
       6b1179bb3 Adding parameter to ivar consensus
       a463f5b34 Changing caller
       75a071ac4 Addition of consensus step
       563998512 MGI init commit
       98118422b Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       34dac0c61 cutadapt bfx
       9ccadd399 DNA-Seq - Fix update
       8733abc8e DNA-Seq - Fix update
       e007f9a4d DNA-Seq - Fix update
       93e342ce7 DNA-Seq - Fix update
       0160ab53d DNA-Seq - Fix update
       567bd8a9d DNA-Seq - Fixmate with samtools and sorting with sambamba
       9e4119b0b DNA-Seq - Fix update
       7ed58f895 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       aefdd2c6c DNA-Seq - Fix update including input_file_revision_eh changes
       8a94a76ef DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       36fbd8068 DNA-Seq - Fix update
       cb43b4863 DNA-Seq - Fix update
       069de1c20 DNA-Seq - Fix update
       2775e6f9e DNA-Seq - Fix update
       348c2c614 DNA-Seq - Fix update
       809b43dd2 DNA-Seq - Fix update
       08c418306 DNA-Seq - Fix update
       7314301c9 DNA-Seq - Fixmate with samtools and sorting with sambamba
       b4728a855 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       d22a3b113 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7a608bef4 DNA-Seq - Fix update including input_file_revision_eh changes
       96dcc43ac MGI-Seq - First commit
       c6e80af21 Adding snpEff step for mgi pipeline
       0e998c0a0 Fixing sambamba flagstat
       79e811492 Fixing sambamba module
       3bd9321b8 Adding module sambamba to ini
       812e317b9 Fixing sambamba indexing
       b0cbd9318 Switching from picard to sambamba withing methylseq
       3f5d5a3be Fix
       271b81d57 fix
       6b0c09fa4 Fixing input output files
       79e49f4ba Fixing filtering bam
       c4dee0fba Fixing renaming step
       f13112485 Fix renaming step
       c765138d3 fix
       c04a496ea Fix
       eea7c7ba8 Fixing filtering
       c0bc132f2 Adding filtering step
       7d34c4e44 Switch to fgbio
       f38cd7d91 ivar triming switch
       15b7438f0 Switching to fgbio trim primer
       cec53cbc0 Updating metrics
       b8de4e31b MGI init commit
       11b01b0da cutadapt bfx
       e3858aa55 DNA-Seq - Fix update
       1eaa3dc1e DNA-Seq - Fix update
       d23743fc7 DNA-Seq - Fix update
       bd59ddbbb DNA-Seq - Fix update
       4eb32689b DNA-Seq - Fix update
       9d8e58ba1 DNA-Seq - Fix update
       c93e23392 DNA-Seq - Fix update
       2934b708e DNA-Seq - Fixmate with samtools and sorting with sambamba
       ff53bb30c DNA-Seq - Fix update
       5d3fb26d3 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       bb731b117 DNA-Seq - Fix update
       1b6813225 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       3a6e72dcf DNA-Seq - Fix update including input_file_revision_eh changes
       01551e523 DNA-Seq - Fix update
       dcaedc85a DNA-Seq - Fix update
       a01251ad7 DNA-Seq - Fix update
       56a50d0c8 DNA-Seq - Fix update
       32e070f2f DNA-Seq - Fix update
       e4778b077 DNA-Seq - Fix update
       94560456d DNA-Seq - Fix update
       083f6a698 DNA-Seq - Fixmate with samtools and sorting with sambamba
       67603dda6 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       a59d9940c DNA-Seq - Fix update including input_file_revision_eh changes
       c20391216 MGI-Seq - First commit
       c5a1baa56 Switch to fgbio
       b4ca5f770 Fixing alignment
       f3378e14a ivar triming switch
       8a0045dca Switching to fgbio trim primer
       d479b8565 Fixing ivar trim
       80c578a2e Filtering reads
       30f51c7d4 Updating metrics
       ea1149b20 fgbio
       22daf31fc Fixing ivar trim
       646290822 Using ivar trim instead of fgbio
       fc6de5f64 Adding ivar primer trimming
       7049f3116 Default bwa parameters to include pairs
       0e834bc96 Zipping and indexing vcfs
       0935479dc Adding parameter to ivar consensus
       d8751ccce Changing caller
       1e6dd7499 ivar bfx
       3399cbb15 Addition of consensus step
       f46ccb5c1 ivar module
       52f8a2b15 MGI init commit
       3bd093242 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       2fbd1e343 Adding trim_primers to fgbio bfx
       327b963f8 cutadapt bfx
       7652e5b8e DNA-Seq - Fix update
       0ab9060cd DNA-Seq - Fix update
       66b240b8c DNA-Seq - Fix update
       472c6b0e4 DNA-Seq - Fix update
       e1e7318a8 DNA-Seq - Fix update
       6954eaf75 DNA-Seq - Fix update
       db7bab265 DNA-Seq - Fix update
       2a0bb0045 DNA-Seq - Fix update
       60e7d9a89 DNA-Seq - Fixmate with samtools and sorting with sambamba
       7be45c10f DNA-Seq - Fix update
       1dddc4f40 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       27ac18b64 DNA-Seq - Fix update including input_file_revision_eh changes
       748e668ad DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       0b5378609 DNA-Seq - Fix update
       6218cdd22 DNA-Seq - Fix update
       56547aa0c DNA-Seq - Fix update
       05baa00f0 DNA-Seq - Fix update
       ff5a87127 DNA-Seq - Fix update
       591339e79 DNA-Seq - Fix update
       de100ee11 DNA-Seq - Fix update
       5e3e598c8 DNA-Seq - Fix update
       584181f91 DNA-Seq - Fix update
       02439f06c DNA-Seq - Fixmate with samtools and sorting with sambamba
       90e9ee942 DNA-Seq - Fix update
       ac5669180 DNA-Seq - Fix update
       f73bce0cc DNA-Seq - Fix update
       22e030d5d DNA-Seq - Fix update
       692892673 DNA-Seq - Fix update
       0facc342e DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c6bf187f7 DNA-Seq - Fix update
       111cbcb5d DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ef0a25e29 DNA-Seq - Fix update ini files
       2eed96b38 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       417ff67ab DNA-Seq - Fix update bai output
       114d004a5 DNA-Seq - Fix update indexing
       fc9901c1f DNA-Seq - Fix update including input_file_revision_eh changes
       0016ba2b2 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       77e504140 MGI-Seq - First commit
       9447bed74 Renaming outputs and default genome
       6f456559f Fixing consensus renaming
       87a626484 Fixing consensus renaming
       cb955341a Renaming + flaging consensus seq
       6c7058d39 Fixing quast options
       1ec14ff9b Fixing qualimap output
       90eb6b0d5 Fixing outputs quast and snpeff
       3dbf0a1d1 Fixing mkdir
       8cd2ba342 Fixing typo
       200cb4bd4 Fixing typo
       2b01f81dc Changing to samtools bam2fq and fixing options ini
       d5dc5b77b Fixing quast
       f32dd8143 Fixing latest commit
       d58a7b237 Adding intermediate bam file on host removal step
       dbd486b3a Forcing pigz to overwritte gz
       95505ba25 fixing host removal
       a11910095 Adding pigz within bam2fq to be able to skip step
       de7633c43 quast fix
       d0383f5fb quast + snpeff + host removal fix
       f165e567b Fixing cutadapt input files
       2dfd2e621 Fixing type
       6e67b6da9 Removing indexing after name sorting
       4b0740ea2 Fixing path creation at host removal step
       7020fe0df Fixing snpeff
       5e6b7a912 Removing print files
       031d9dfa1 Fixing input choosing
       256cdc790 Fixing choosing inout files mapping
       cd7e52e9f Fixing host removal
       c11d2200a Fixing host removal
       a3d58c9c5 Fixing quast step
       8abe6c089 Fixing quast
       5a68e072f Fixing quast step
       b141fa8f4 Fixing param requirements
       a7318fdc9 Fixing typo
       aa32b62e3 Fixing pigz
       617f817bc Not using kraken anymore
       df8f32cef Switching steps order
       10c3586d6 Adding 3 steps
       eb7d6120e Fixing picard multiple metrics raw
       8d852a76a Insert Size metrics on raw bam
       67dcc5f49 sambamba flagstat fix
       08fb9264b Renaming snpeff step/job
       068d072a4 Fixing flagstat
       770c233bc Fixing flagstat
       4ff0fb6ed flagstat on all the bams
       1e3ef28e8 Zipping output of snpeff
       5c219e64f Fixing snpeff
       c0faa3758 Fixing cutadapt
       effce8a8f Flagstat on raw bam
       8c7dd5a00 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       221dea23b Adding snpEff step for mgi pipeline
       86a0acf3b Fixing sambamba flagstat
       dd8163f7e Fixing sambamba module
       609e6ca85 Adding module sambamba to ini
       71316e973 Fixing sambamba indexing
       904072204 Switching from picard to sambamba withing methylseq
       9d3bdfd17 Fix
       e5a7f40df fix
       59aeaff89 Fixing input output files
       da3b1cd0e Fixing filtering bam
       cb775da33 Fixing renaming step
       f78909cb2 Fix renaming step
       49dba7dbd fix
       50ac6ee71 Fix
       860591e55 Fixing filtering
       cc0f319cb Adding filtering step
       4d36888f5 Switch to fgbio
       6066f5265 ivar triming switch
       5790ea2eb Switching to fgbio trim primer
       45c244912 Updating metrics
       ffc8fd16f MGI init commit
       13d4f105d cutadapt bfx
       2a61965b9 DNA-Seq - Fix update
       7c8d9e343 DNA-Seq - Fix update
       aad95a269 DNA-Seq - Fix update
       db01cdb0a DNA-Seq - Fix update
       c8fd1f4a6 DNA-Seq - Fix update
       19b7a74cc DNA-Seq - Fix update
       386a60487 DNA-Seq - Fix update
       6063ddda7 DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e1f0b545 DNA-Seq - Fix update
       30121a774 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       4b9d8e1cc DNA-Seq - Fix update
       a83ee80a5 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       033dd3ec8 DNA-Seq - Fix update including input_file_revision_eh changes
       a9102d385 DNA-Seq - Fix update
       a4ee47fe8 DNA-Seq - Fix update
       18d832e51 DNA-Seq - Fix update
       6623f8c3a DNA-Seq - Fix update
       1d002d1f1 DNA-Seq - Fix update
       fe5c4057d DNA-Seq - Fix update
       dbd8034b7 DNA-Seq - Fix update
       74d2bdde8 DNA-Seq - Fixmate with samtools and sorting with sambamba
       167863986 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       5f66227b8 DNA-Seq - Fix update including input_file_revision_eh changes
       562180b90 MGI-Seq - First commit
       a16695b5a Switch to fgbio
       1d92909dd Fixing alignment
       01f5df72d ivar triming switch
       e49e632b7 Switching to fgbio trim primer
       057e2e8f8 Fixing ivar trim
       d9638b557 Filtering reads
       c3ab6f146 Updating metrics
       dc9a43eee fgbio
       579e0dd47 Fixing ivar trim
       603bb4a09 Using ivar trim instead of fgbio
       26a2da86d Adding ivar primer trimming
       f25c75d25 Default bwa parameters to include pairs
       6a8c3018d Zipping and indexing vcfs
       cf6c19f1a Adding parameter to ivar consensus
       c9231e5c0 Changing caller
       b9bfed08d ivar bfx
       70cc54a38 Addition of consensus step
       5a990e2b7 ivar module
       cd17c0d05 MGI init commit
       5bef158f3 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       f7a155d21 Adding trim_primers to fgbio bfx
       8177037ae cutadapt bfx
       a702d2779 DNA-Seq - Fix update
       f6fe0e2db DNA-Seq - Fix update
       4a4243975 DNA-Seq - Fix update
       a880e3bf3 DNA-Seq - Fix update
       2cd5d0128 DNA-Seq - Fix update
       9f7e04801 DNA-Seq - Fix update
       75b18caf9 DNA-Seq - Fix update
       9a0adea97 DNA-Seq - Fix update
       bb259fe40 DNA-Seq - Fixmate with samtools and sorting with sambamba
       0c969370d DNA-Seq - Fix update
       75e589b3e DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       dd3999c1d DNA-Seq - Fix update including input_file_revision_eh changes
       2a2a43c82 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       cb09a7b3e DNA-Seq - Fix update
       31653bd6b DNA-Seq - Fix update
       513818caf DNA-Seq - Fix update
       b202a5b68 DNA-Seq - Fix update
       c87c2cd2f DNA-Seq - Fix update
       a4926cbe5 DNA-Seq - Fix update
       c06ee52ba DNA-Seq - Fix update
       6f894d2c9 DNA-Seq - Fix update
       3d5daeb10 DNA-Seq - Fix update
       ed9e7feb1 DNA-Seq - Fixmate with samtools and sorting with sambamba
       4215929b8 DNA-Seq - Fix update
       3dca79c7e DNA-Seq - Fix update
       7bf5141ad DNA-Seq - Fix update
       ff24ef917 DNA-Seq - Fix update
       f5352026f DNA-Seq - Fix update
       a6b428c0c DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       7e351683c DNA-Seq - Fix update
       302e6c9af DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       03fc8c296 DNA-Seq - Fix update ini files
       0cd488577 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       ade827942 DNA-Seq - Fix update bai output
       8ee326719 DNA-Seq - Fix update indexing
       4f25a9848 DNA-Seq - Fix update including input_file_revision_eh changes
       314107e45 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       6d82d3f5b MGI-Seq - First commit
       2cbfa5c0a Adding snpEff step for mgi pipeline
       8f5d83948 Fixing sambamba flagstat
       4415a7aef Fixing sambamba module
       644bc8d04 Adding module sambamba to ini
       38cd45adc Fixing sambamba indexing
       08dc201cd Switching from picard to sambamba withing methylseq
       e89813ffd Fix
       9ef7ecbd1 fix
       0d6651a85 Fixing input output files
       192bb6521 Fixing filtering bam
       030b759fc Fixing renaming step
       1348f004a Fix renaming step
       bc42658a6 fix
       49a16fc79 Fix
       c69973dc0 Fixing filtering
       76219a861 Adding filtering step
       7ed54e8bf Switch to fgbio
       b6cfddfc4 Fixing alignment
       4796cffb3 ivar triming switch
       3ec5f7f3c Switching to fgbio trim primer
       3a84e9eb8 Fixing ivar trim
       b673f473e Filtering reads
       9fdafbd90 Updating metrics
       3f63d5ca1 fgbio
       56bec970b Fixing ivar trim
       30475caa0 Using ivar trim instead of fgbio
       d5b40548a Adding ivar primer trimming
       925645e1b Default bwa parameters to include pairs
       314a86953 Zipping and indexing vcfs
       191bfb3d9 Adding parameter to ivar consensus
       db7c70a3f Changing caller
       d9e302b1b ivar bfx
       09d9f6b50 Addition of consensus step
       d9c3f32f3 ivar module
       f6a8e1723 MGI init commit
       49c7f1d05 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       eae7037e3 Adding trim_primers to fgbio bfx
       bfc47be8b cutadapt bfx
       ffcfe262a DNA-Seq - Fix update
       0f1742656 DNA-Seq - Fix update
       34b096584 DNA-Seq - Fix update
       b7323de8b DNA-Seq - Fix update
       e0909f148 DNA-Seq - Fix update
       48e203ed3 DNA-Seq - Fix update
       0cd7e5e63 DNA-Seq - Fix update
       b3d740100 DNA-Seq - Fix update
       6ff06b006 DNA-Seq - Fixmate with samtools and sorting with sambamba
       b000f97af DNA-Seq - Fix update
       52deaf398 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       2acf3c0b8 DNA-Seq - Fix update including input_file_revision_eh changes
       0fe44df6d DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       3772b9bcd DNA-Seq - Fix update
       479e9ed40 DNA-Seq - Fix update
       5f2512077 DNA-Seq - Fix update
       aa1687cb1 DNA-Seq - Fix update
       8632426c9 DNA-Seq - Fix update
       c5d384b8e DNA-Seq - Fix update
       746d2f01d DNA-Seq - Fix update
       eef3b276d DNA-Seq - Fix update
       f1d1b5513 DNA-Seq - Fix update
       d2d4430ec DNA-Seq - Fixmate with samtools and sorting with sambamba
       568329727 DNA-Seq - Fix update
       a7e450cc6 DNA-Seq - Fix update
       5cb957404 DNA-Seq - Fix update
       e888e0986 DNA-Seq - Fix update
       d366f6753 DNA-Seq - Fix update
       ba6aa73c4 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       9dc3bf6ff DNA-Seq - Fix update
       c01dcf527 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ca81d44f2 DNA-Seq - Fix update ini files
       78519f294 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       0c922ccde DNA-Seq - Fix update bai output
       6117bc989 DNA-Seq - Fix update indexing
       5e9b26857 DNA-Seq - Fix update including input_file_revision_eh changes
       519beea32 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       43bb771a0 MGI-Seq - First commit
       f2cd11aa3 Switch to fgbio
       c3aa00dd6 Fixing alignment
       8d157fd4f ivar triming switch
       42a5e2d7d Switching to fgbio trim primer
       6030fa732 Fixing ivar trim
       6bb7e087e Filtering reads
       330e345c0 Updating metrics
       b1153cd53 fgbio
       4450c70e6 Fixing ivar trim
       a68648b03 Using ivar trim instead of fgbio
       8a2bdf36d Adding ivar primer trimming
       8856367bb Default bwa parameters to include pairs
       4e968647f Zipping and indexing vcfs
       6e9529f8f Adding parameter to ivar consensus
       4b826ef9e Changing caller
       b8bc79a0a ivar bfx
       dab5579f7 Addition of consensus step
       c91d06767 ivar module
       ec2a17602 MGI init commit
       260059a7c Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       7760a795b Adding trim_primers to fgbio bfx
       3180972e5 Fixing haplotype_caller
       66d8cdda3 Fixing haplotype_caller
       3bdb0aa97 cutadapt bfx
       2d274eb20 Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       f083b6680 MGI-Seq - First commit
       538907e99 Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into dnaseq_picard_fix
       b0197fe12 DNA-Seq - Fix update
       54b4741cb DNA-Seq - Fix update
       452f86923 DNA-Seq - Fix update
       d0cb631d2 DNA-Seq - Fix update
       6e3580e67 DNA-Seq - Fix update
       222b79b02 DNA-Seq - Fix update
       5cf626990 DNA-Seq - Fix update
       de8d9017b DNA-Seq - Fix update
       d87cd96c9 DNA-Seq - Fix update
       b3b9b1b6e DNA-Seq - Fixmate with samtools and sorting with sambamba
       5e55bfc9a DNA-Seq - Fix update
       9f20fe780 DNA-Seq - Fix update
       123b73033 DNA-Seq - Fix update
       5ae600acf DNA-Seq - Fix update
       8119523be DNA-Seq - Fix update
       d6007d20e DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       1ff8cacc9 DNA-Seq - Fix update
       58d662a7e DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       83ba60d0d DNA-Seq - Fix update ini files
       f7bef7df8 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       10d50aa8a DNA-Seq - Fix update bai output
       1319c26b9 DNA-Seq - Fix update indexing
       fe4e7b99e DNA-Seq - Fix update including input_file_revision_eh changes
       80def2add DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       894566ccf DNA-Seq - Fix update
       1413d4519 DNA-Seq - Fix update
       ab945dd03 DNA-Seq - Fix update
       97186dd25 DNA-Seq - Fix update
       fc158a472 DNA-Seq - Fix update
       86b2e803f DNA-Seq - Fix update
       d2ecd4b92 DNA-Seq - Fix update
       98928a551 DNA-Seq - Fix update
       bdee1f800 DNA-Seq - Fix update
       8afb5005c DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e8699ba5 DNA-Seq - Fix update
       1e47ab113 DNA-Seq - Fix update
       bea5e948c DNA-Seq - Fix update
       eae982e83 DNA-Seq - Fix update
       15d745e22 DNA-Seq - Fix update
       0553a5de6 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       ef79ced10 DNA-Seq - Fix update
       afaa484b0 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7dfa7303f DNA-Seq - Fix update ini files
       3d22e6728 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       6dd6e0c67 DNA-Seq - Fix update bai output
       54a5fc055 DNA-Seq - Fix update indexing
       20cee1f0a DNA-Seq - Fix update including input_file_revision_eh changes
       a1579f469 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       c5e99661b Merged in ihec_metrics (pull request #126)
       7bceefb1b Merge branch 'ihec_metrics' of bitbucket.org:mugqic/genpipes into ihec_metrics
       388c9bf4d ChIP-Seq - Fix IHEC metrics
       b5683f252 ChIP-Seq - Fix IHEC metrics
       7c0f68fb3 Merged in chipseq_atacseq_mode (pull request #121)
       7a101a707 Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       e83caef01 MethylSeq - Trimmomatic resources revisited for methylseq
       f9bbdb79e General - Genome installation with cvmfs grep
       a8c4d934d ChIP-Seq - Fixing GenPipes metrics
       d82ed8993 ChIP-Seq - Sambamba loading properly in mito calcul
       e6c4274b1 ChIP-Seq - Reducing ihec report ressources
       f3755e930 ChIP-Seq - Fixing typo
       315301796 ChIP-Seq - Increasing IHEC preprocess ressources
       a02fd3228 ChIP-Seq - Fix IHEC metrics
       7f1f4a964 ChIP-Seq - Fix IHEC report template md file
       805c58b61 ChIP-Seq - Fix IHEC report template md file
       19b0b6292 ChIP-Seq - Fix IHEC report template md file
       3c573f8d6 ChIP-Seq - Fixing merge metrics
       1bb96cb2d ChIP-Seq - Fixing merge metrics
       971431309 ChIP-Seq - Fixing merge metrics
       ad5e53139 ChIP-Seq - Fixing merge metrics
       f3ebb74c9 ChIP-Seq - Fixing merge metrics
       b5658b8bb ChIP-Seq - Adding IHEC report template md file
       3ad67e06e ChIP-Seq - Fixing merge metrics
       bef0f1873 ChIP-Seq - Fixing merge metrics
       ec3dea001 ChIP-Seq - Fixing merge metrics
       4883d4c3e ChIP-Seq - Fixing merge metrics
       51bd5f2cd ChIP-Seq - Adding merge IHEC metrics
       c310d1e2c ChIP-Seq - Fixing metrics
       b1ec49d44 ChIP-Seq - Fixing metrics
       508417135 ChIP-Seq - Fixing metrics
       fa035a103 ChIP-Seq - Fixing metrics
       85c8c53b2 ChIP-Seq - Fixing metrics
       41081d625 ChIP-Seq - Fixing metrics
       1acdfa5f4 ChIP-Seq - Fixing metrics
       cfe40591d ChIP-Seq - Fixing metrics & report
       d29881f12 ChIP-Seq - Fixing metrics
       d1d438525 ChIP-Seq - Fixing metrics
       5e801f947 ChIP-Seq - Fixing metrics
       a1a478424 ChIP-Seq - Fixing metrics
       0f847f18f ChIP-Seq - Fixing metrics
       2f9084afb ChIP-Seq - Fixing metrics
       e545fb483 ChIP-Seq - Fixing metrics
       56f6af0df ChIP-Seq - Fixing metrics
       878335c1c ChIP-Seq - Fixing metrics
       380e1224d ChIP-Seq - Fixing metrics
       f6f26f4bb ChIP-Seq - Fixing metrics
       343c9ffa5 ChIP-Seq - Fixing metrics
       3da997d38 ChIP-Seq - Fixing metrics
       db54da14a ChIP-Seq - Fixing metrics
       6a91f0d3f ChIP-Seq - Adding metrics
       b6a1e6085 ChIP-Seq - Fixing bwa missing import
       319655b79 ChIP-Seq - Fixing chipseq pipeline
       1a47e31b9 ChIP-Seq - Adding ATAC-Seq protocol
       5f113fce5 Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       a440d4e6a General - Genome installation with cvmfs grep
       d7cd595e8 ChIP-Seq - Fixing GenPipes metrics
       bbfb6bb48 ChIP-Seq - Sambamba loading properly in mito calcul
       466041cb7 Merge branch 'dev' into chipseq_atacseq_mode
       15611bfaf ChIP-Seq - Reducing ihec report ressources
       56c39f6f7 ChIP-Seq - Fixing typo
       f88ee5938 ChIP-Seq - Increasing IHEC preprocess ressources
       0bd494dd2 Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       1d727d08d ChIP-Seq - Fix IHEC metrics
       09d7280f1 Merged dev into chipseq_atacseq_mode
       364e3811e ChIP-Seq - Fix IHEC report template md file
       24b64ed84 ChIP-Seq - Fix IHEC report template md file
       9f8c10e76 ChIP-Seq - Fix IHEC report template md file
       0d811bd43 ChIP-Seq - Fixing merge metrics
       70a5ed0d5 ChIP-Seq - Fixing merge metrics
       b1da60049 ChIP-Seq - Fixing merge metrics
       aee5b149d ChIP-Seq - Fixing merge metrics
       74f44e239 ChIP-Seq - Fixing merge metrics
       b096b6350 ChIP-Seq - Adding IHEC report template md file
       12b2f31fc ChIP-Seq - Fixing merge metrics
       a807a5e3e ChIP-Seq - Fixing merge metrics
       5131a3231 ChIP-Seq - Fixing merge metrics
       60b0f6fc2 ChIP-Seq - Fixing merge metrics
       05c262315 ChIP-Seq - Adding merge IHEC metrics
       9664631ae ChIP-Seq - Fixing metrics
       1e2bc38a5 ChIP-Seq - Fixing metrics
       c6a08a6fb ChIP-Seq - Fixing metrics
       e9e96da2e ChIP-Seq - Fixing metrics
       809cb52b8 ChIP-Seq - Fixing metrics
       78736476f ChIP-Seq - Fixing metrics
       4b1bb18e3 ChIP-Seq - Fixing metrics
       e806abc6c ChIP-Seq - Fixing metrics & report
       ca01c9cdd ChIP-Seq - Fixing metrics
       3c72d0f93 ChIP-Seq - Fixing metrics
       1199e5f95 ChIP-Seq - Fixing metrics
       1ee3d3643 ChIP-Seq - Fixing metrics
       59a282840 ChIP-Seq - Fixing metrics
       51570690f ChIP-Seq - Fixing metrics
       4a8a555fa ChIP-Seq - Fixing metrics
       2cb7a7bea ChIP-Seq - Fixing metrics
       a48c88ee4 ChIP-Seq - Fixing metrics
       235552a25 ChIP-Seq - Fixing metrics
       5b47faf96 ChIP-Seq - Fixing metrics
       31535af01 ChIP-Seq - Fixing metrics
       2d751ea57 ChIP-Seq - Fixing metrics
       85ea7e4ed ChIP-Seq - Fixing metrics
       b62f51937 ChIP-Seq - Adding metrics
       bbe8b2c23 ChIP-Seq - Fixing bwa missing import
       5bdb8aa22 ChIP-Seq - Fixing chipseq pipeline
       2ae58031a ChIP-Seq - Adding ATAC-Seq protocol

  Paul Stretenowich <pstretenowich@CYPRUS.local>      4 commits

       20aa88fd2 Cleaning mgi.py
       cc25234ce Cleaning mgi.py
       397efae3c Cleaning mgi.py
       97bd77659 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      15 commits

       406e38b91 Merged in fix_monitor (pull request #190)
       e885338e3 Merged in monitor_bug (pull request #189)
       af56bf447 Merged in chunk_slurm_submit (pull request #185)
       bfe072adb Merged in chunk_slurm_submit (pull request #182)
       b10eac60e Merged in rsync_in_chipseq (pull request #180)
       31826dedf Merged in remove_pacbio (pull request #176)
       68d5195ad Merged in poq/fast_module_check (pull request #164)
       652e207d3 Merged poq/graham_ini into dev
       9c6da9a87 Merged poq/graham_ini into dev
       1f66c6c99 Merged update_beluga_ini into dev
       70167f13e Merged update_beluga_ini into dev
       e7c7f3493 Merged update_beluga_ini into dev
       a54d824c0 Merged in add_container (pull request #167)
       8a076c758 Merged in poq/graham_ini (pull request #168)
       f18307f49 Merged in master (pull request #161)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      57 commits

       67c0af4b2 Merge branch 'release_3.2' of bitbucket.org:mugqic/genpipes into release_3.2
       90052a65a update readme for container install
       1294b9636 Remove tumor_pair pipeline from release
       04c0c5a34 revert on indel aligner
       d75d45849 cleanup ini for beluga
       cd1255139 remove dbSNP for oxog ini block fix gatk 4 module name
       7d84e335d fix picard_collect_oxog_metrics dbSNP vcf
       97d6a47f6 remove HOME_DEV paths
       70d131090 make job2json more robust
       a05a2430e remove mugqic_dev vawk
       7baea3295 Fix README conflict
       3bda6d70e Fix out.out in monitor
       f0b99b865 rename to monitor.sh
       1ec8ba819 update control batch description
       bd24c9af9 rename deley_sbatch script
       c0f801243 fix script usage
       c579431b7 add export to sourced files
       f00a2d12c make sure already submited jobs id are sourced
       0a4444d4e if error
       285860809 remove debbug line
       5f6e768b7 add delay and chunking script
       88026d507 update for cit testing
       92edc027f replace -a by -r in rsync
       2551c7878 fix master conflic with deleted readme file
       241f706fa remove pacbio from the repo/release
       fd02e500e no mail in cit.ini
       68d8a669f add sh script for steps in pbs too
       d79b27368 Add SARS-CoV2 genome file
       b4c47402c tweek memory usage beluga denovo
       7fcd09ae5 update cluster for rnaseq star index
       473294db2 use 1.1.0 genpipes_in_container release
       0d591bbaa make module show sure it raise with older version
       1aa57cbdd chipseq cedar and graham ini
       1329a13cc Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6f7be7f23 more graham ini
       b455c4cf3 more graham ini
       826760335 copy dnaseq ini from cedar to graham
       b4485973b fix fastq symlink on graham
       153aa9e27 less log
       40ab0ff13 cleanup
       0ac31218a speedup module check
       20ff188d9 one module call
       1fb3ae757 add mem to star_* on rnaseq
       4d682b610 update again
       19a81cf19 mem-per-cpu
       e5ce42025 Add generic ini for chipseq update README
       d10376104 create graham ini file
       05ef09e56 log when all sbatch submit is dont
       c1f2c5313 feedback at submit time
       0ac30b35e feedback at submit time
       fd755d815 fix ini typo
       77b2c9b70 wrapper for slum, pbs and batch
       bc80a4d83 update wrap options.
       304c7eb5f add default wrapper
       22d9efac7 remove docker
       004a9a438 put --wrap import at the top
       3516da40d add wrapper to all pipelines

  P-O Quirion <pioliqui@gmail.com>      4 commits

       0e486370a tweek monitor and chunk
       4f50d9e8a Fix autorestart and cleanup on failes interrupts
       328a00efd make nanopore executable
       ab1de9a91 Add automatic wrapper option

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      26 commits

       4f04bf4e0 covseq mpileup command fix
       acfe7b9a4 covseq qualimap fix, and high coverage ini fix
       9e0c8d790 add covseq - dnaseq consistency
       c6e5099e1 dnaseq high coverage trimmomatic to skewer
       8409952a9 tumor_pair fixes
       0c79ba9fc import fixes
       5fb88cb6d bam.bai to .bai fix
       8894fc91c cit fixes to b38 - samtools and other b38 annotations
       87b0da868 conflict fixes to dnaseq and tumor pair after dev merge
       5a4dc7c42 fix mpileup dependency chain
       f8cda00c9 corrections to mpileup bcftools merge and tumor_pair dependencies
       b8ee48005 bam.bai to .bai fix
       4c7b3f7a4 fixes to indentation in sambamba merge
       cb807fc6e fixes to sambamba merge
       3bc8f0ed8 further gatk4 hc fixes
       0cdd2dfca update gatk4 hc arguments with gatk4 suite update
       8d7b9b828 fixes to gatk4 mark dup
       e7029c61a gatk4 mark dup walltime fix
       1d4995456 variant recal fix
       9095ef11b variant recal dependency fix
       927011a33 cit fixes to b38 - samtools and other b38 annotations
       03ee662cd minor dnaseq.py fixes
       4ffeb590f conflict fixes to dnaseq and tumor pair after dev merge
       cd935c7f2 fix mpileup dependency chain
       ef39fb7b7 corrections to mpileup bcftools merge and tumor_pair dependencies
       4d72c00d8 conforming deliverables to cit conventions

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      29 commits

       91f765514 exome specific fixes
       e6d4caf02 updates and fixes for cit
       72d072d40 remove briaree ini and update dnaseq base
       1547aa725 updates to beluga.ini and base.ini for dnaseq
       8e64f79ca ini updates
       c925a6123 gatk4 vsqr cit fix and baf plot for b38
       dd8c84473 add cram to input selection
       8b2d9c30e multiqc fix
       55d8a64ca argument fixes for picard imported functions in gatk4 and vqsr fixes
       12f55ac8b indel realignment consistency issues between dnaseq and tumor_pair
       1c6f1bcbc add mark dup metric file to multiqc
       4964f1657 cit b38 fixes
       01f6a48ee fix to bash.ln for cit
       2e1a585c8 fixes to multiqc
       d5cbe666e exome specific fixes
       e407f1d8f updates and fixes for cit
       a3bd732a8 Updates to bcftools/samtools for dnaseq_mpileup
       d96272302 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       fc430949f remove briaree ini and update dnaseq base
       72e8439f5 adding 1 job processing to specific steps for cit.ini
       f4911b3f5 fix of dev genome reference
       2447075ce issues between dnaseq.base.ini and dnaseq.beluga.ini
       cc7c806c5 updates to modules for beluga
       d5f030157 updating beluga ini
       67758960e added tumor pair beluga ini
       16bcac2f6 updates to beluga.ini and base.ini for dnaseq
       8e8f22e90 updates to beluga.ini
       f8ea2637d ini updates
       de84a173b module updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      15 commits

       369bd8a46 picard2 and high coverage fixes
       84885daeb updates to b38 variant recal files
       3db1c5aff fixes to tumor_pair on beluga
       a40514c67 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       e04c3f881 fixes to fixmate input file name
       9fe39b097 picard2 and high coverage fixes
       30cf93955 tumor_pair beluga ini fix
       c77c8210f tumor_pair qualimap part 2
       5c039c640 qualimap tumor_pair fix
       b980f2a0e fixes to vsqr gatk4
       ed9960349 gatk4 fixes to callable loci and DoC
       6162e8e4c sym link dnaseq.base into tumor pair
       950a28544 updates to b38 variant recal files
       f6c42cba6 fixes to tumor_pair on beluga
       d8eaf5823 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      24 commits

       2f926d675 fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       800ecb253 cit fixes after rebasing
       a37298c3c dependency fixes
       ba400c9b1 updates to GRCh38 annotation file, module updates, cit fixes
       01f024967 updates to cit and fixes to one job mpileup steps
       3e213bdef updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8b745c15e major fixes to deliverables and completion of beluga test
       3078a9ac3 fix modules dev to cvmfs
       2345cafc7 fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       3c284035b fixes to merge_filter_bcf
       662e007a8 gatk4 bsqr fixes and mpileup cat
       0f3aefcfc mpileup protocol fix
       aefbb370e cit fixes after rebasing
       6f47ad31e update to gatk4 mutect2 filtering procedures
       419f376d4 create cit for gatk4 due to deviation from argument usage
       04fcaacce cit fixes to gatk4 + sym links for recalibration
       802cc6583 dependency fixes
       7c3b8d3b3 updates to GRCh38 annotation file, module updates, cit fixes
       28828267c fixes to deliverable and b38 ini
       81204f691 updates to cit and fixes to one job mpileup steps
       dfa32887a updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8bc536d3c fixes to metasv annotations
       556e48b1e major fixes to deliverables and completion of beluga test
       973087a7f fix modules dev to cvmfs

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      1 commits

       b5c148684 final PR fixes

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      7 commits

       9ac20daaf code cleaning and fixes to exome interval list
       f55a73846 fixes to cedar ini
       deca50fb8 fixes to symlinks for paired indel realignment
       31525eea6 code cleaning and fixes to exome interval list
       fa7328f4b fixes to cedar ini
       2d3d36761 fixes to symlinks for paired indel realignment
       5475ca4a6 cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      5 commits

       7fbd4f9f5 fixes to metasv, adding metasv germline
       5a6a9d75d cedar fixes and GRCh38 fixes
       b42ed1669 cedar dnaseq updates and svaba germline added
       2e389e923 cedar germline sv updates
       67a2e47d4 sequence dictionary and module updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      35 commits

       7d71caf83 updates to metasv - somatic
       eada99c81 updates to SV germline and reference genome tweaks
       0aa9d42b7 somatic sv fixes: lumpy and svaba annotations
       0f2779026 fix to germline SV: breakseq2
       df8096aca merge fixes
       a96622a55 json and folder updates
       ee184e24f fixes to sCNAphase
       fa8f99fcb merging snv and sv, adding protocols
       6cae335e8 Fixes to indel realigner
       16099054e Updates and debug
       04229bad8 Add set somatic and actionable mutations
       937e5438d added multiqc and other tweaks
       2be1a9742 add metrics step for metasv
       a9697c7b6 updates to metasv - somatic
       f622b9220 Fixes and updates to reference files
       909c2efbd remove testing steps
       0cb7a4b5e updates to SV germline and reference genome tweaks
       de92d5759 somatic sv fixes: lumpy and svaba annotations
       573c45dda fix to germline SV: breakseq2
       f2fc9d66a merge fixes
       77bc11902 GATK4 fixes - bam indexing and markDupSpark
       21051bce0 bcftools fixes for tumor pair
       44f3d044f fingerprint and bug fixes
       c9c204959 dnaseq - vcftools qc addition: --missing_indv and --depth
       9428baa6b select input for variant caller and fixes to one job calling
       1c78b285e json and folder updates
       3731f2053 fixes to sCNAphase
       402bccbe2 Added json info
       036606a81 Bug fixes prior to json additions
       bd6e6cb92 merging snv and sv, adding protocols
       d441d4714 Fixes to indel realigner
       c07293c05 Add deliverables module
       f674ec9f2 Updates and debug
       3161cf6d4 Add set somatic and actionable mutations
       e8d90db03 added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      14 commits

       23f5be7bd remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       614ee8f45 fixes to metasv for tumor pair
       9f1b500bb Bug fixes and modification derived from initial PROFYLE benchmarking
       245df59d1 remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       079fd7e61 fixes to metasv for tumor pair
       48bb3c5ef Single job bug fixes
       b77cfaaa6 manta I/O fix and other bug fixes
       d14ba6d78 config updates and b38DH added
       6b5287904 dnaseq germline SV updates
       b82f0f67b gatk4 updates and bug fixes
       14310b441 gatk4 updates and bug fixes
       77126e8e1 Fix line break types
       9fe9453ed Json related bug fixes
       dcf0b8ac5 Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       04b549010 fixes to breakseq2 and metasv
       b0626ecac cit-based fix - verifyBAMid
       b897ef1d2 dnaseq qc additions: NGScheckmate and peddy
       0cbfb0df5 fixes to breakseq2 and metasv
       612c44664 cit-based fix - verifyBAMid
       52023e1cc cit-based fixes to NGScheckmate
       4dec63d8e dnaseq - multiqc fix
       29e26f720 dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      4 commits

       f3c972e60 Merged in rebasing_tp (pull request #186)
       f6e74bd27 Debugging and Guillimin specfic fixes
       4810e9a4c Debugging and Guillimin specfic fixes
       0d16e4dc6 updates to config

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       123f6c577 Merged in fix-watch-folder (pull request #165)

  Rom Grk <romgrk.cc@gmail.com>      3 commits

       1d5912c7b Merged in fix-watch-portal-folder (pull request #187)
       c22599c51 fix: use of undeclared variables
       1baeef5ab watch_portal_folder.py: fix undefined variable

  ufgauthi <ulysse.fortiergauthier@mcgill.ca>      1 commits

       0a60ef29f Bug Fix by replacing sacct delimiter | by ^

  Ulysse Fortier Gauthier <ulysse.fortiergauthier@mcgill.ca>      1 commits

       bf54dbf99 Merged in ufg_log_report_fix (pull request #157)

3.1.5        Wed Jan 15 11:58:16 2020 -0500        424 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      78 commits

       e0844c309 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       986f2a3c5 HICSeq pipeline - improving trimmomatic resources for SULRM + added graham config file for hicseq pipeline
       6714a5d41 Genome installation - added America Mink genome (Neovison_vison.NNQGG.v01) installation script
       d60d13688 Software installation - added 'demuxlet' installation script
       6822cbfb3 Software installation - updated ucsc.sh with lat4est version i.e. v387
       b350d2da4 Software installation - replaced call of lsb_release by universal commands (avoid 'lsb_release command not found' error)
       5c87fd720 Software installation - corrected regular expression within genome_digest function
       23e93153c removed some test code introduced by latest commit...
       9445ee982 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dev
       00577d4fc GenPipes RNA-Seq de novo - corrected bugs introduced in commit 4a72735
       f3b9a92c4 GenPipes BFX - added bash_cmd python wrapper to wrap basic bash commands
       4a727351e GenPipes RNA-Seq de novo - updated pipeline with better sample assignemnt to jobs for better JSON building
       d89d4fc30 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1b4821503 GenPipes JSON - updated pipeline.py and scheduler.py so that the copy of the original JSONs to the PORTAL is made before building the pipeline steps
       6360349eb GenPipes Trimmomatic - corrected assignment of samples to the job to correct analysis JSON generation per sample
       fbaa84a9f MethylSeq - corrected assignment of samples to the jobs to correct analysis JSON generation per sample
       d1cd79b3c GenPipes analysis JSON - updated jsonator.py in reagrds of the updated .base.ini files, now containing 'source' & 'version' to ease the process + minor updates to core/pipeline.py
       fc1b6c0f8 GenPipes config files - updated most of the .base.ini files with missing 'source' & 'version' to avoid issues with the 'jsonator'
       f44ccb156 GenPipes Anlalysis JSON file - corrected core/pipeline.py regarding the use of PORTAL_OUTPUT_DIR
       9751513e5 GenPipes Analysis JSON file - added the system to update the submission_date
       a9ba83209 GenPipes Anlalysis JSON file - pipelines now create JSON analysis file as a default behavior
       7b42d155a GenPipes Analysis JSON file - added the submission_date to the JSON
       7e04404f1 Analysis JSON - updated jsonator.py with project_name and submission_date. Also updated version to 1.0.1
       f422649a6 GenPipes utils - minor updates : updated some comments
       38937e611 Software install - updated install_module.sh with finer LIBDIR for patching and better patching to avoid overwritting potential pre-existing RPATH
       fb2bffc05 Software install - updated R_Bioconductor.sh with new R packages and finer LIBDIR for patching
       721d9416e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       074045aac GenPipes Sanity Check - Refining the report
       0d0e0bdd1 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       1478455a8 GenPipes Sanity Check - Adjusting the log level
       cf3ec9e32 GenPIpes Sanity Check - set log level to 'warning' instead of info when running sanity check, thus removing some useless if statements and refining display of the messages
       6cf36ee5f GenPipes Sanity Check mode : removed some useless comments in core/pipeline.py
       b46df3f81 GenPIpes Sanity Check mode - created the SanitycheckError class in config.py to use instead of Error
       f4d6323cd Merge branch 'master' of bitbucket.org:mugqic/genpipes into rnaseq_denovo_jhg
       989996a1a GenPipes MethylSeq - updated base.ini with [samtools_cram_output] section
       066c24e3d Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       433950d27 GenPipes Sanity Check - DNASeq high coverage, RNASeq denovo assembly & RNASeq light pipelines are updated regarding the sanity-check mode
       0fbecd2eb Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       6b00ce5db GenPipes Sanity Check - updated pipelines for sanity-check mode responsiveness
       afa406ab5 GenPipes Sanity Check - updated common.py to be sanity-check responsive
       4f25a0fa4 GenPipes Sanity Check - update _raise() function in config.py & refined code of sanity-check mode in pipeline.py
       6a457d530 GenPipes Sanity Check - updated readset.py design.py to be sanity-check responsive & tools.py with proper import statements
       6f11825cb MethylSeq pipeline - updated pipeline to make SINGLE-END mode actually work
       b3ac5a9e4 GenPipes Sanity Check mode : updated PacBio Assembly pipeline as the first try to test sanity check mode
       f28854499 GenPipes Sanity Check mode : updated pipeline.py with the sanity check mode fully functionnal
       b5a78ddad GenPipes Sanity Check mode : updated sample.py to reflect the updates that have been done in config.py
       3203ce901 GenPipes Sanity Check mode : updated config.py to avoid raising errors when sanity-check mode is on, logging them instead
       daabbf45d GenPipes code cleaning - cleaned some 'import' calls : stop using 'import *' and specify which modules/classes/functions to import instead
       2699bc945 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       975e0ff33 GenPipes Core - updated pipeline.py to catch , relative to changes in config.py. Added the sanity-check mode : the pipeline now does not stop at fist met error/exection : wait untill the end to print the errors
       1f628d353 GenPipes Core - updated config.py wth PO's code : catches  instead of
       0f4c8605d Software upadte - updated mugqic_tools installation script with version 2.2.3
       68f32a979 Software update - updated R_Bioconductor.sh : added binless and DoubletFinder packages to the installation, updated the installer depending on the version of R
       2bca735e3 Saoftware update : update cellranger-atac.sh with the latest versoin 1.1.0
       147832184 Software update - added th shebang to the C3G wrappers for installed binaries
       697eaec2a Genome installation - corrected typo in 'lambda_page'
       e7010eb8d Software update - added script to install gemBS in C3G softwatre stack
       9c4f29546 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       93f92cdc2 Software update - updated installation scripts for LUMPY-SV, mugqic_tools & SMRTLink
       39b691a34 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       06fdb584f Software update - updated MultiQC installation script : shebang of MultiQC scripts now uses #!/usr/bin/env python
       e66f717b9 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1e4c96fe2 Software update : added the installation of methylKit packages in R_Bioconductor.sh - also added the generation of a file listing all the packages installed for all the new installed R versions
       251f418c5 GenPipes utils - updated nanuq2mugqic_pipelines script : added support for iSeq projects
       742aa73dc GenPipes utils - updated nanuq2mugqic_pipelines script : added -rn/--run parameter, standing for Nanuq run ID, to fetch only readsets procesed in specified run(s)
       555517196 Software update - updated install_module.sh so that it does not wrap nor patch the executable binaries when installing on DEV space
       d4fb87c77 Software update - added MiXCR v3.0.5 installation script - MiXCR: a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data
       f020d9206 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3e7d563e2 Software update - updated some installation scripts with the latest version of software : bowtie v1.2.2 - bowtie2 v2.3.5 - MultiQC v1.7
       09335ff98 install_genome.sh - added the generation of the 100-bin GC content file for bisulfite genome references
       4e1331fa5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       04047dfc4 Genome installation : added Bisulfite Genome Reference creation and indexing with Bismark
       8916af044 updated LIBDIR in R_Bioconductor.sh
       0ad84563a update Bismark version to 0.21.0
       eff800250 updated python installation script with version 3.7.3 and added pip as a symling to pip3
       db27d82a6 Added fastq software installation script
       e1f627ff3 Version bump to 3.1.5-beta
       1786fb37d Version bump to 3.1.4

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      1 commits

       650797a2d GenPipes - DEBUGGING - DNASeq SamToFastq & SymLink steps corrected + working bfx/bash_cmd.py

  Édouard Henrion <henrione@gra-login2.graham.sharcnet>      4 commits

       9e9adc788 GenPipes JSON - debugged call to job2json when output_dir is different that '.'
       b7eb0ed23 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       754a7d981 GenPipes software update - updated patching in install_module.sh
       5bf7663e8 GenPipes software update - updated R_Bioconductor.sh with the latest requested libraries & updated patching

  Édouard Henrion <henrione@ip18.m>      3 commits

       140b411d4 GenPipes Scheduler - Corrected bug in job2json call
       c77a53525 Merge branch 'eh_methylseq_single_end' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       a592ddfd6 GenPipes MethylSeq - updates for single-mode

  Édouard Henrion <henrione@ip20.m>      1 commits

       41ed7a2c6 Genome Installation - updated install_genome.sh with common grep version

  ehenrion <edouard.henrion@mcgill.ca>      83 commits

       25f72f000 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1578931756922 (pull request #163)
       64b7414f9 DNASeq.py - removed the use of 'os/path.abspath' in call of 'ln()'
       3c750f004 bash_cmd.py - added import os
       7a05d36e2 DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable
       e9213a632 bash_cmd.py - remove call to "sleep" in the ln() function. Replaced it by a call to "ls" in order to flush the cache
       01154915a GenPipes - BUG correction - trailing space removed in output files of DNASeq skewer trimming step
       1fdf5799d Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578602990982 (pull request #162)
       b74f38dcc BUG correction : corrected bash_cmd.py "ln" function with correct format keys when buiding command
       07d826fd2 Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578588554527 (pull request #158)
       828d5736c Bug Correction : corrected 'ln' function in bash_cmd.py to avoid a python "TypeError: cannot concatenate 'str' and 'int' objects" exception
       420558272 Merged in ehenrion/dnaseqbelugaini-edited-online-with-bitbu-1576097030188 (pull request #147)
       fcaef48a4 GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       520210092 Merged in eh_quick_fixes (pull request #143)
       815632f02 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918024470 (pull request #141)
       92a7c1170 GenPIpes - dnaseq.py : bug correction - typo removed
       44d9eb52c Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575399438745 (pull request #132)
       152664276 Merged in ehenrion/commonpy-edited-online-with-bitbucket-1575398053532 (pull request #131)
       482a0bba9 GenPipes - bug correction in pipelines/common.py : corrected name for bam and fastq files created by the pipeline : don't depend and raw file names anymore, but built from sample & readset name given in the readset file
       2a25960a9 GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       72e7f1e5e GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       c946bb0e9 Merged in ehenrion/schedulerpy-edited-online-with-bitbucket-1575392160612 (pull request #129)
       9bf49cbe1 Merged in ehenrion/nanuq2mugqic_pipelinespy-edited-online-w-1575392552042 (pull request #130)
       3893e10c2 GenPipes : nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       2f38d5d86 GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       a876d020a GenPipes - corrected scheduler.py with missing argument line
       8b8bca82c Merged in dev (pull request #125)
       7cd7b1cd8 RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - dowgraded trinity version to 2.0.4_patch
       656db125f RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - corrected typo in trinity version...
       01c2bee3f RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - updated trinity version to 2.2.0_patch which contains a patch (from C3G developppers) to avoid  printing buffer 'hiccup'
       b9325d16d CIT - RNASeq Denovo Assembly - update default cluster_walltime to 4:00:00 in cit.ini
       8235190e2 Merged ehenrion/cit-dnaseq-updated-cit_cluster_walltim-1574180387297 into dev
       238abec9f CIT - DNASeq - updated 'cit_cluster_walltime' for gatk_callable_loci step in cit.ini
       6a891c207 DNASeq - gatk_callable_loci - adjusted memory and cpu allocation in dnaseq.beluga.ini
       dec27c8eb CIT - DNASeq - Adjusted trimmomatic resource allocation through java_other_options, threads settings and mem-per-cpu use
       3d3424333 CIT - RNASeq - corrected 'cluster_walltime' for wiggle  step in cit.ini
       a24e7b182 RNASeq - wiggle step - adjusted memory allocation using 'mem-per-cpu' instead of 'mem'
       3e3ee1d13 CIT - RNASeq - updated java threads to 5 through 'java_other_options' in rnaseq.beluga.ini
       1a1d0f1e0 CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       bae00cddc CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       cc5df685a CIT - HiCSeq - redefined hicup_align walltime and mem-per-cpu in cit.ini
       5cd0d02ff CIT - Pacbio assembly - set specific walltime for smrtanalysis_summarize_polishing in cit.ini
       1c1979135 HICSeq pipeline - trimmomatic resources revisited in hicseq.beluga.ini : removed buffer_size
       29024f2cd CIT - RNASeq denovo Assembly - set localfit to true in cit.ini for differential_expression_deseq
       622b41ad2 CIT - Pacbio Assembly - redefined walltime for pacbio_tools_split_reads in cit.ini
       09e223549 RNASeq denove Assembly - edited rnaseq_denovo_assembly.beluga.ini adjusted memory assignment for insilico_read_normalisation_readsets
       3a25fef04 RNASeq denovo Assembly : edited rnaseq_denovo_assembly.beluga.ini with beter resources assignement
       b4a3e981d Software update - trimmomatic.py - updated trimmomatic command with the use of the  java_other_options parameter provided by the ini files
       47e5e41c1 Merged in ehenrion/rnaseq-metricspy-ihec_metrics_rnaseq--1573489397505 (pull request #124)
       87c3190b5 RNASeq - metrics.py - ihec_metrics_rnaseq : added file in the input_file list to correct job dependency
       fc680b198 RNASeq cufflinks - correcting dependencies for ihec_metrics : needed rnaseqc report file to be added to the output_file list
       25f6cc3fd CIT - Pacbio assembly - set threads in beluga.ini
       aa3d185d1 CIT - RNASeq_light - set threads for trimmomatic in beluga.ini
       b19924252 CIT - RNASeq_light - redefined walltime for callisto_count_matrix in cit.ini
       24d79ff45 CIT - Pacbio Assembly - redefined walltime for preassembly in cit.ini
       3c8f736bf CIT - RNASeq denovo Assembly - redefined walltime for align_estimate_abundance in cit.ini
       3134f05c9 CIT - RNASeq denovo Assembly - updates cit_cluster_walltime to 24h
       24351cb02 CIT - RNASeq denovo Assembly - redefined walltime for transdecoder and align_estimate_abundance in cit.ini
       c17d7918b CIT - Pacbio Assembly - redefined walltime for smrtanalysis_filtering in cit.ini
       19ed125f5 CIT - DNASeq - redefined walltime for sambamba_merge_realigned in cit.ini
       d92aa276d CIT - DNASeq High Coverage - redefined trimmomatic walltime in cit.ini
       74979fa59 CIT - HiCSeq - redefined trimmomatic walltime in cit.ini
       ed4fb136c CIT - RNASeq_light - redefined walltime for trimmomatic in cit.ini
       0cf351d12 RNASeq denovo Assembly pipeline - rnaseq_denovo_assembly.py - corrected command quoting in trinity step
       8386936b4 RNASeq denovo Assembly Pipeline - rnaseq_denovo_assembly.beluga.ini - corrected max_memory parameter for trinity step : using Gb instead of Mb
       373bd60c5 RNA-Seq de novo assembly pipeline - rnaseq_denovo_assembly.beluga.ini - corrected Jellyfish memory parameter : using Gb instead of Mb unit
       d2e4c46b3 HiC-Seq pipeline - hicseq.base.ini - corrected typo in HiCUP module name
       ceaae5345 HiC-Seq pipeline - hicseq.base.ini -  update HiCUP version to v0.7.2
       7860ea8b8 dnaseq.py - corrected file names dependencies
       c198c42e1 dnaseq.py - corrected wrong file extension in metrics_vcf_stats step
       bc952b663 HiC-Seq pipeline - corrected typo in cram_output options definition - hicseq.base.ini
       d8c5f5810 rnaseq_light.cedar.ini : corrected typo in cluster_server assigned value
       19bf83c41 rnaseq_light.mp2b.ini corrected minor typo on cluster_server
       27f7080bd rnaseq_denovo_assembly.py : corrected differential expression jobs with samples
       5010b026a corrected typo in ampliconseq.beluga.ini
       1f1601724 DNASeq High Coverage README.md updated with `picard_fixmate` step documentation
       eda0617f0 GenPipes - Readme updated for `picard_fixmate` step in dnaseq_high_coverage.py
       a05f490fd Merged in eh_methylseq_single_end (pull request #115)
       6c88ebf37 Merged in sanity_check_mode (pull request #114)
       fa256237d rnaseq_denovo_assembly.beluga.ini :  updated [insilico_read_normalization_all] with missing core allocation
       305d96e93 README.md updated : removed Guillimin settings section from README
       56d668c8c README.md updated : corrected the setting of MUGQIC_INSTALL_HOME_DEV for mp2b
       8d4f1b471 rnaseq_denovo_assembly.mp2b.ini - corrected resources requirments : all the steps now run on one single node !
       33909cfbe chipseq.base.ini - updated MultiQC version to 1.6 to support Homer

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       6c76ddc37 Merged in lefebvref/rnaseq_denovo_assemblybaseini-edited-onl-1557891680322 (pull request #110)
       53bd32ae0 Previous default expression in base.ini would end up basically only retaining transmembrane proteins. Not good.

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      2 commits

       5bea83ca9 Corrected error with beluga ini for RNAseq denovo
       02b131cff Corrected minor error in the help message that said that the default protocol was cuflinks. Starting from version 3.1.4, stringtie is the default protocol

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      10 commits

       153d4c4dc Merged in Jose-Hector-Galvez/update-to-the-rnaseq_light-readmemd-file-1569265173884 (pull request #120)
       e24a06615 Update to the rnaseq_light README.md file to address issue brought up by GSoD collaborator.
       5f4086f35 Merged in Jose-Hector-Galvez/readmemd-edited-online-with-bitbucket-1565184718156 (pull request #118)
       ed95482eb README.md edited to remove warning labels.
       341cab2f0 Merged in rnaseq_denovo_jhg (pull request #103)
       05e7dfe87 Merged in rnaseq_jhg (pull request #112)
       9a24beafe Merged dev into rnaseq_jhg
       6f98c7691 Merged master into rnaseq_denovo_jhg
       550514749 Merged master into rnaseq_jhg
       160079daf Merged master into rnaseq_denovo_jhg

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      3 commits

       e618184d0 Corrected error in the beluga ini for RNAseq de novo
       e9610a9f7 Corrected error in the beluga ini for RNAseq de novo
       ceade4228 Corrected error in the beluga ini for RNAseq de novo

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       cc5eac7fa add header in the cram file (all base.ini)
       92b5c4f2c add cram to tumor_pair pipeline
       78224133f implement the launch of the entire pipeline if no step argument is given
       689732a1e add cram to dnaseq_high_coverage pipeline
       71569ce81 add cram to methylseq pipeline
       c242a05a4 add cram to rnaseq pipeline
       70ccb0370 add cram to hicseq pipeline
       0e788734d add cram to chipseq ini files
       830cd38d6 add cram to dnaseq ini files
       f9d780a60 add cram to ChipSeq
       5f1d736f4 correct bugs in dnaseq light
       0313f9cc0 add cram creation to dnaseq pipeline
       bfde21459 add generic function to create CRAM from BAM
       195651499 make samtools view's output not mandatory removable

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      2 commits

       52cd608e8 temporary fix for tbi missing output in haplotypecaller when running only 1 job
       6a98813ca Merged in cram (pull request #111)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      9 commits

       ee86893af Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0fd7a8894 Nextflow install script
       cc8694c1d Merged in rnaseq_rnaseqc (pull request #116)
       d59ea5ee0 General - Correcting genome ini path
       caa989b8e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into rnaseq_rnaseqc
       8e90fcf02 RNA-Seq Pipeline - Typo update
       6bfe800b0 RNA-Seq Pipeline - Typo update
       bb8a0cdd0 RNA-Seq Pipeline - Removing verify_bam_id step
       7eb13c4ac Chip-Seq Pipeline - Adding homer_make_ucsc_file to mp2b.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      12 commits

       a8d1aae7e Merged in add_container (pull request #122)
       9b20dde92 Merged in cit_fix_cedar_2019-04 (pull request #109)
       4eb85d86c gatk4.py added 2 spaces in multi line string
       ef09195b7 ampliconseq.py edited online with Bitbucket
       9e04a50db pacbio_assembly.cedar.ini edited online with Bitbucket
       2b0230622 Merged in log_repport_ci (pull request #108)
       04e9b6ba5 Merged in log_repport_ci (pull request #107)
       f171ba188 Merged in log_repport_ci (pull request #106)
       0da3f3ec5 Merged in log_repport_ci (pull request #105)
       f9458d8cc Merged in log_repport_ci (pull request #104)
       1b7e74f67 Merged in fix_multiqc (pull request #89)
       0e66455be Merged in fix_multiqc (pull request #77)

  Pierre-Olivier Quirion <poq@beluga3.int.ets1.calculquebec.ca>      2 commits

       feeaf91da Log when loading Job. More robust output parser.
       5181b5cce Log when loading Job. More robust output parser.

  P-O Quirion <pioliqui@gmail.com>      158 commits

       161d1d192 fixed motifMaker ini: used default memory
       ba32a4785 typo in waltime kallisto_count_matrix
       50b824bd2 trimmomatic thread ini
       20c0f5675 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       795edd94b string tie ini
       d4e255de5 sync thread in java for trimomatic
       9bd93fa93 differential_expression_and_goseq_rsem return list list
       42bf1d808  fix scheduler logging; reduce sleep to 0.1 sec
       0902f427e ignore .pyc
       bca87cc0e remove the methylkit_differential_analysis step form the pipeline
       f9cef1b3f add gitignore
       dd0267e21 missing import
       3d6b4e702 log regression scheduler
       f4ab24a87 fix beluga ini
       6d5d769cf fix beluga ini
       ef68329f0 Merge branch 'add_container' into merge
       a872d634d Merge branch 'master' into dev
       d31246f41 fix batch system account for container
       313fcb207 remove commented code
       281e3cd78 Merge branch 'add_container' of bitbucket.org:mugqic/genpipes into add_container
       ebb67ba72 Add default config for container ini
       116615f7c fix regression
       7542a0986 fix scheduler when container not in use
       8c9d6e3ce use mugqic stack java in rnaseq on cedar.ini
       bb6fc694a simple default ini for container
       b82cacedd default ini for container
       331715dec Basic container ini file
       b787ff7ed Working singularity version tested on graham
       2fe8f4756 WIP exec line prototype
       b2b4d970b add container option
       215d35277 Merge branch 'cit_fix_cedar_2019-04' of bitbucket.org:mugqic/genpipes into cit_fix_cedar_2019-04
       4092e778e amplicon R tweak
       7b389fcc6 missing quotes in pacbio
       138f3356c amplicon memory
       32627c87d amplicon memory
       141554d56 sleuth again wrong name
       8c34b2be6 renaming sleuth to differential_expression_sleuth
       963d105ad fix typo in ampliconseq ini
       ff92d21fe fix amplicon method typo
       da2cdf800 tweak default mem for rnaseq denovo and light
       75e2f696d more ampliconseq deadlock: otu_assigning
       1cfac7351 parallel_assign_taxonomy_uclust.py deadlock
       6ecfef8eb rna fix
       cecb51446 update cit
       b4997b3af packbio tmpdir again
       c39257690 chr19 option for cit chicagorun
       46d9149f8 fix log report on output id
       5cfc33cb9 fix log report on output id
       0cc3e1793 update waltime in amplicon and regex in log
       5cc0ce57e migrate form testdata ini to cit
       2c02e926e mem teak
       51bf13daa tweak to get in bynode queue
       3f05dbe89 tweak ini for cit
       078d0e0a7 dnaseq ini tweak for slurm; sh_create_baitmap \\
       cf79f461e timout amplicon
       410185e88 remove mpileup_metrics_vcf_stats step from mpileup
       9d100d8ce pbutgcns_wf need full path an tmp need to exist
       ad31e6e86 ini file tweak
       ac6bf87ae tweak mem usage
       b4f2cf801 ini tweak plus perl module load in mummer
       3ae1b09b1 tweak mem usage per-cpu for ampliconseq
       bb8fc310f Ensembl v 90 for GRCh38
       890cd4841 mummer using deprecated perl syntax
       16326cfb1 more default time for rnaseq cufflink
       3d1aa2479 mugqic_tools/2.2.4 in rnaseq_light
       ffcf62f62 samll sample for both gatk call in variant_recalibrator
       7b25d7bcb change R_Bioconductor version/cleanup code
       bc97eae00 update rnaseq_light references
       505ffeb8f tweek memory usage in dnaseq_high_coverage
       3ce4bb612 add time to ampliconseq
       c1b5a6a3a reformat errors in the core.config module
       05aa94b43 redirect small_sample_option to the right function
       3d51765c4 tweak cit for dnaseq
       720c3fa79 fix ampliconseq dada2
       390ceee05 tweak VariantRecalibrator for testing
       f74afdff8 fix slurm_time_to_datetime
       733f258a3 set cit ini file and config
       7d9c898bf add runtime to csv output
       86e84f3eb Add GENPIPES_CIT env var to run in testing mode
       4cb6e2409 fix all slumtmpdir sting renaming problem
       8f3652ffa update mugqic_tools in hicseq
       2107daaf9 second run of pipeline fixing
       458e919eb remove extra gz extention from anotation
       690acf415 flash_log is a list of one
       b829c3166 update hicseq.base
       81a4690ea add interactive option in smrtanalysis bash cmd
       785521720 fix regression
       18aa48bfa fix scheduler when container not in use
       ae7824c03 use mugqic stack java in rnaseq on cedar.ini
       b8eee7fc5 simple default ini for container
       2522d29dd default ini for container
       60de9f181 Basic container ini file
       a5b2901d2 Working singularity version tested on graham
       5229f0ee6 WIP exec line prototype
       56e0a9b03 add container option
       27a2cb37b amplicon R tweak
       738c59dae missing quotes in pacbio
       d73e16599 amplicon memory
       a770708a6 amplicon memory
       80498da4a sleuth again wrong name
       53c5205b3 renaming sleuth to differential_expression_sleuth
       1e5d20e08 fix typo in ampliconseq ini
       7bb900da5 fix amplicon method typo
       85b51a0d6 tweak default mem for rnaseq denovo and light
       7df9cae89 more ampliconseq deadlock: otu_assigning
       777bb337e parallel_assign_taxonomy_uclust.py deadlock
       24bd06ed4 rna fix
       a9e2370d1 update cit
       b120e4cc1 packbio tmpdir again
       8e1835fdb chr19 option for cit chicagorun
       eec3f9530 fix log report on output id
       5488adc09 fix log report on output id
       34960c034 update waltime in amplicon and regex in log
       ca1fe4d33 migrate form testdata ini to cit
       eca1836ad mem teak
       058f07179 tweak to get in bynode queue
       10ff44a81 tweak ini for cit
       c7568c6e5 dnaseq ini tweak for slurm; sh_create_baitmap \\
       1c50e41de timout amplicon
       ff1a9c7d1 remove mpileup_metrics_vcf_stats step from mpileup
       d8d8a3dd3 pbutgcns_wf need full path an tmp need to exist
       48ff75bbf ini file tweak
       68d6f4689 tweak mem usage
       160a12345 ini tweak plus perl module load in mummer
       ee6a7aa38 tweak mem usage per-cpu for ampliconseq
       f0acafe63 Ensembl v 90 for GRCh38
       23788c3bc mummer using deprecated perl syntax
       14d0500af more default time for rnaseq cufflink
       4e3a8b07d mugqic_tools/2.2.4 in rnaseq_light
       624781cf5 samll sample for both gatk call in variant_recalibrator
       8307e5963 change R_Bioconductor version/cleanup code
       069dd4aa1 update rnaseq_light references
       3e01c15db tweek memory usage in dnaseq_high_coverage
       c0c5905f7 add time to ampliconseq
       d29d92675 reformat errors in the core.config module
       43fdedaaf redirect small_sample_option to the right function
       cfefb02a5 tweak cit for dnaseq
       509d3330a fix ampliconseq dada2
       927a732f4 tweak VariantRecalibrator for testing
       0a2659d14 fix slurm_time_to_datetime
       0f48c2df2 set cit ini file and config
       20c852142 add runtime to csv output
       02ce535db Add GENPIPES_CIT env var to run in testing mode
       dbbb0cd30 fix all slumtmpdir sting renaming problem
       19c6c30a3 update mugqic_tools in hicseq
       27e516456 second run of pipeline fixing
       32a1a40d6 remove extra gz extention from anotation
       4edbfc572 flash_log is a list of one
       a396c5ff7 update hicseq.base
       d806fced2 add interactive option in smrtanalysis bash cmd
       d9567ac79 add space
       31ddd6cfb add check for output job id
       63431c423 get real path in log_report
       bec1796c1 file exist file missing mixup
       4120c3aba bug fixes
       3e94acda1 get pipeline info in tsv file can also mute stdout
       29f5670cc run dnaseq multiqc on all samples at once
       8d7321737 add indentation

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       43f70c631 fix to preprocess step

  Robert Syme <robsyme@beluga4.int.ets1.calculquebec.ca>      1 commits

       85814368d We don't need to track compiled python .pyc in git

  Robert Syme <rob.syme@gmail.com>      7 commits

       30ea46dca Merged in spacesfix (pull request #117)
       0c4d6ce47 Remove spaces from Rscript commands.
       1f3c0f66d More trailing space cleanup.
       534124c36 Cleanup trailing spaces
       16b7d01d2 Merged in amplicon-pathfix (pull request #113)
       63d804af5 Add support for --output-dir argument
       eb8389c37 Trailing whitespace cleanup

  Rola Dali <rola.dali@mail.mcgill.ca>      31 commits

       92202e880 dnaseq_high_coverage.base.ini edited Varscan version to 2.4.3
       3f0b73205 dnaseq.mp2b.ini edited online with Bitbucket
       07b6376bb dnaseq.beluga.ini edited online with Bitbucket
       a0bad0770 dnaseq.cedar.ini edited online with Bitbucket
       f379676f1 dnaseq.beluga.ini edited online with Bitbucket
       2805dc692 dnaseq.cedar.ini edited online with Bitbucket
       a4c341871 pacbio_assembly.beluga.ini edited online with Bitbucket
       bdd8b7244 dnaseq.beluga.ini edited online with Bitbucket
       b648fb5ae dnaseq.cedar.ini edited online with Bitbucket
       287e82b89 dnaseq.beluga.ini edited online with Bitbucket
       a760a50a0 dnaseq.cedar.ini edited online with Bitbucket
       62de3fb77 hicseq.beluga.ini edited online with Bitbucket
       1fd778781 hicseq.cedar.ini edited online with Bitbucket
       42eb77ac0 ampliconseq.cedar.ini edited online with Bitbucket
       fabb3d5c5 ampliconseq.mp2b.ini edited online with Bitbucket
       193c97bae ampliconseq.beluga.ini edited online with Bitbucket
       5c4c908a9 ampliconseq.cedar.ini edited online with Bitbucket
       59491ffe4 dnaseq.cedar.ini edited online with Bitbucket
       9c9a75fd0 dnaseq.beluga.ini edited online with Bitbucket
       a2e7e80f2 rnaseq.beluga.ini edited online with Bitbucket
       d1777dd77 hicseq.beluga.ini edited online with Bitbucket
       5a6a19872 hicseq.cedar.ini edited online with Bitbucket
       855a02c3e hicseq.cedar.ini edited online with Bitbucket
       a4cb4c16e hicseq.beluga.ini edited online with Bitbucket
       efe78e4d1 README.md edited online with Bitbucket
       c0b67c17f README.md edited online with Bitbucket
       5c7825129 README.md edited online with Bitbucket
       982f934f3 README.md edited online with Bitbucket
       4402aeffa README.md edited online with Bitbucket
       ed0fd9e44 README.md edited online with Bitbucket
       b97304ef6 README.md edited online with Bitbucket

3.1.4        Tue Mar 26 14:03:32 2019 -0400        198 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      69 commits

       cea585601 Updated all the GenPipes .base.ini files with the latest verison of mugqic_R_packages i.e. 1.0.6
       892cef59b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       cc01f4ee8 MethylSeq - added reference to Bismark in markdown config_and_references.md file for a better report document - also updated the metrics table assignement in MethylSeq.ihec_sample_metrics_report.md and bfx/report/MethylSeq.ihec_sample_metrics_targeted_report.md
       540ee0ed0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into slurm_report
       3ab4a6f08 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f8c58f3db updated mugqic_tools version to 2.2.2 in ampliconseq.base.ini
       4a100ec8b updated R_mugqic_packages.sh with version 1.0.6 of the mugqic R packages
       1d490a7a3 updated pipeline inis file with the latest version of mugqic_tools i.e. v2.2.2
       d1ad83446 MethylSeq - MethylKit DMR - to fit with C3G code standards, getting rid of bfx/methylkit.py and use bfx/tools.py instead to call methylkit.R, which is part of mugqic_tools
       d5e1184da MethylSeq - methylKit DMR - correted typo in bfx/methylkit.py
       a989966c0 MethylSeq - methylKit DMR - correted typo in methylseq.base.ini
       ba9c8b135 MethylSeq - MethylKit DMR - change the call to methylKit.R : now using R instead of Rscript
       c6c82bf39 Pac Bio Assembly pipeline - updated the bfx wrapper circlator.py : make sure th eoutput is removed before launching circlator, otherwise it will throw an error
       a763f6d88 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       9650b322b PacBio Assembly pipeline - minor update in pacbio_assembly.py for code reading purposes
       e41a72ee8 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       0ddffdd59 updated svaba.sh installation script with latest version of SvABA : 1.1.0
       43e752701 MethylSeq pipeline - DMR analysis - updated R module version in the base.ini file to make sure methylKit is available
       3576f8d00 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       90d9206ab updated mugqic_tools version to 2.2.2 in methylseq.base.ini
       34f4ccb43 update mugqic_tools.sh with the latest version : 2.2.2
       f1a34d08f updated install_module.sh to better handle LD_LIBARAY_PATH in the LIBDIR
       0fb5b2cc2 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0aa9f92ce corrected Octopus installation script
       a5a7bd9ec Pacbio Assembly pipeline - corrected basemodification and MotifMaker bfx subroutine wrappers : now sourcing /etc/setup.sh
       bfefcf727 modified R module in install_genome.sh : now using mugqic/R_Bioconductor/3.5.1_3.7
       805760f79 updated mugqic_tools.sh with version 2.2.1
       02d4ef837 Adding installation scripts for CMake and Octopus
       10a90429d Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       193fe2846 MethylSeq pipeline - DMR analysis - change file namings in prepare_methylkit step
       f6f936ffb MethylSeq pipeline - DMR analysis - adjusting walltime for filter_snp_cpg step
       aca1e34ba MethylSeq pipieline - DMR analysis - clean code in bfx/tools.py
       732a39acd updated wrapping and patching in install_module.sh
       56891fac8 updating software install scripts with newer versions
       82723b7aa bedtools.py - removed unnecessary argument to bedtools.coverage subroutine
       d126c1b77 MethylSeq pipeline - methylseq.py - fixing inputs assignment in IHEC metrics step
       4eead82c0 MethylSeq pipeline - methylseq.py - correcting sort_sam call and affecting job name to bismark methylation call job
       27f85ae63 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e23a9821d BFXDEV-429 - updated syntax of temporary diretory after sourcing of /etc/setup.sh for smrtanalysis tools
       53e51310b updated resources/modules/R_Bioconductor.sh script : refine LIBDIR definition - revised wrappers creation
       3f2091282 updated install_module.sh : refine create_c3g_wrappers subroutine to catch executable binaries - updated LIBDIR definition for centOS cases
       58286be3b updated delly.sh script : lasted version 0.8.1 & added DELLY_PATH in the modulefile
       8a1084f94 udpated lumpy-sv.sh script : added a sed command to make lumpyexpress.config compliant with our environment
       1c7b2c449 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       58bde8f8e PacBio Assembly - updated bfx/smrtanalysis.py as part of the solution of BFXDEV-429
       7ee6531af Updated syntax to tmp_dir to finish fixing BFXDEV-429 - first usefull related update was made in commit ad3caf1, when bumping to version 3.1.3
       807bee4fd minor - updated syntax
       e93e06f28 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       dafb3c08e removing occurences of DEV paths within the ini files of GenPipes
       3d256c1b1 Adjusting .ini section for some htslbfx subroutines which were pointing to [DEFAULT] section of the .ini file
       3fe0d06d2 BFXDEV-752
       25bc30d86 Adjusting .ini section for some qualimap ray sambamba and other subroutines which were pointing to [DEFAULT] section of the pipeline .ini files
       6d4dd16ad corrected some typo in jsonator.py to have fastq & bam file paths properly handled
       630de05a3 Adjusting .ini section for some htslib subroutines which were pointing to [DEFAULT] section of the .ini file
       5b03880f4 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       891780f32 comment line removed in install.module.sh
       f6720fb4c Added ctpr.sh script to install CTPR package in C3G software stack
       0a3e8ab12 Adjusting .ini section for some blast * blat subroutines which were pointing to [DEFAULT] section of the .ini file
       3ebb9cc8d BFXEDV-529 - edited RNAseq to consider star_cycle_number instead of cycle_number for star.index step - also ajusted read_size value for sjdbOverhang value completion
       54fd04013 Homo_sapiens.GRCh38.sh - Added the creation of a Ribosomal RNA interval list file, for use of Picard Tools CollectRnaSeqMetrics
       fe136f799 Adjusting .ini section for some bcftools subroutines which still were pointing to [DEFAULT] section
       709347d25 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       4344c947d changing samples assignment in trimmomatic jobs for a better analysis JSON creation
       321c65d9f Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       24e90fa0f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9cc90a8c0 added genpipes.sh
       0157c9234 updated GATK installation script
       ea084430f Version bump to 3.1.4-beta
       7e318dc18 Version bump to 3.1.3

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       9cc62b8b9 MethylSeq Pipiline - MethylKit DMR analysis - updated Cedar .ini file with relevant cluster resources, changed filter_snp_cpg subroutine int bfx/tools.py to use bedops instead of bedtools, added bedops module in the .base.ini file
       e26c2afd3 MethylSeq Pipeline - updated filter_snp_cpg call with more cluster resources and removed pipes to avoid 'Broken pipe' error...
       7d77cd8b7 MethylSeq pipeline - DMR analysis - fixing cedar ini file + cleaning code before pull request

  ehenrion <edouard.henrion@mcgill.ca>      14 commits

       3ef3a0439 pipeline.py  : edited the help so that it actually shows that SLURM is the default scheduler used by GenPipes
       d3edaa688 methylseq.py corrected report files in metrics step
       86646c43d Merged in slurm_report (pull request #73)
       b2dbb4693 job.py : corrected typo
       0961c1ef4 job.py : added missing setter class for name and multiqc_files attributes
       4f48cc688 Merged in methylseq (pull request #67)
       10bcac61f dnaseq.cedar.ini : resolved conflicts
       a18369a40 rnaseq_denovo_assembly.base.ini : corrected missing variables
       d80c17883 tumor_pair.base.exome.ini - updated COSMIC path
       1a66a54e9 tumor_pair.base.ini - updated COSMIC path
       a3c68a371 rnaseq_light.base.ini : commetn adapter_fasta because it was pointing to some old place
       0f799c37b BFXDEV-737 - updated ucsc.py for a better handling of Danio Rerio GRCz11 genome build
       b768d04f5 BFXDEV-734 - updated ucsc.py for a better handling of Mus Musculus GRCm38 genome build
       54dab24f7 AmpliconSeq - reassign silva_db to correct path for dada2 analysis - ampliconseq.base.ini

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      28 commits

       7da61f053 Corrected slurm parameters for stringtie jobs, including stringtie_abund which was causing errors
       fb2d51da7 Removed DEV from sleuth configuration
       68878b13f Added note to RNAseq README regarding switch to stringtie protocol by default
       94342a5cc Merge branch 'rnaseq_light_jhg' of bitbucket.org:mugqic/genpipes into rnaseq_light_jhg
       68bf88616 Specified genome support for RNAseq_light in the pipeline's README
       5088cbdcd Corrected INI file to inlcude latest version of Mugqic_tools
       e455bbc9d Modified INI file to include latest version of MUGQIC_tools as well as more memory for kallisto steps by default
       e974d4506 Merge branch 'rnaseq_jhg' of https://bitbucket.org/mugqic/genpipes into rnaseq_jhg
       49bbed545 Corrected minor parameters for abacus on the RNAseq ini
       330c861ca corrected merge issue with the rnaseq_ligth cedar INI
       0d7a7d7fc Corrected minor issues with rnaseq_light pipeline, and added documentation to the README
       160fd60ed Corrected some issues with the UCSC script that was causing issues with non GRCh38 genomes
       52daaeaa7 Added documentation for stringtie protocol. Corrected some problems with the base inif file. Minor corrections to the stringtie script and rnaseq.py script.
       214930b43 Fixed merge issues
       497139c6e merged to master and resolved conficts
       95be45a49 Merged to latest version of master branch
       90e40bad9 Corrected a mistake on the methylseq ini files for Cedar and Mp2b
       5b6841782 Corrected CVMFS Stringtie module in INIs
       615ad8588 Fixed merge conflict with the cedar ini
       548a27f8f Merge branch 'master' into hgalvez
       b6a9e4e89 Second testing version stringtie/ballgown pipeline. Includes ballgown now, as well as corrections for stringtie. Additionally, stringtie is now default protocol for RNA-seq.
       08bbb3cc7 Fixed a small typo with the last commit
       1abef45db Fixed bug for the order or arguments passed to kallisto
       abbbfb2f0 Updated version of Bioconductor used by default in the rnaseq_light pipeline (modified base.ini file).
       a97c715a6 Merge branch 'hector-cedar' of https://bitbucket.org/mugqic/genpipes into hector-cedar
       7a877b314 Minimally functional rnaseq_light pipeline with the added sleuth step. Still requires further testing on clusters beyond abacus. Also, some harmonization of the required genome reference files would be advisable.
       3c142b81f Merge branch 'master' into hector-cedar
       2ccccdea8 Eliminated warning message from rnaseq_light pipeline for single read samples

  jose.hector.galvez@computationalgenomics.ca <hgalvez@abacus2.ferrier.genome.mcgill.ca>      2 commits

       15b6ecbd0 Fixed issues with stringtie merge and abund. Ready to test in other servers
       4f3e5e8b4 Corrected job.py to allow for the definition of the multiqc_files parameter, should fix error with multiqc and Kallisto

  Jose Hector Galvez <Hector@rhododendron.genome.mcgill.ca>      10 commits

       4a5614834 latest modifications to the tools.py script
       9363e82a1 Fixed minor typo on rnaseq.py
       397a903e2 Fixed indent issue in the stringtie job creation
       35f290489 Fixed minor bug in the stringtie job creation
       ebe3d347d Fixed minor bug in the stringtie job creation
       ffd2a51ee Fixed stringtie module location in the rnaseq.base.ini file
       a76c8e4f3 Fixed sample name issue within stringtie.py
       c6b938f52 corrected minor typo
       18c69678c correct minor typo
       01cc53d07 First commit adding stringtie functionality, still testing

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      22 commits

       447c061a5 Merged in rnaseq_light_jhg (pull request #75)
       daabf762c Merged in rnaseq_jhg (pull request #76)
       de0d01eba Merged master into rnaseq_jhg
       7293ec07e Merged master into rnaseq_light_jhg
       e032d8fdb Merged in rnaseq_jhg (pull request #69)
       fecf7c71e Merged master into rnaseq_jhg
       d687b0fe4 Merged in rnaseq_light_jhg (pull request #68)
       e78a608f0 Merged master into rnaseq_light_jhg
       9225ebbd0 Merged master into rnaseq_jhg
       3868e89f5 Merged master into rnaseq_jhg
       573a0e0e0 Merged master into rnaseq_light_jhg
       640d57c5a Merged master into rnaseq_jhg
       77c6491f6 Merged master into rnaseq_jhg
       16c4eb7d0 Merged in missing-cedar-ini-files (pull request #64)
       c95739897 Merged master into hgalvez
       072c23ded Merged hector-multiqc into hector-cedar
       753edc2cc Merged master into hector-cedar
       456948661 Merged hector-cedar into hector-multiqc
       0585acc48 Merged master into hector-multiqc
       ff9fcfb2f Merged master into hector-cedar
       05daa1db3 Merged in hector-multiqc (pull request #38)
       eb7c1fad9 Merged master into hector-cedar

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      3 commits

       b89789779 Added cedar and mp2b ini files for four pipelines: ampliconseq, dnaseq_high_coverage, pacbio_assembly, rnaseq_light
       289f74aa0 modified cedar ini to support stringtie
       ac01fb0c9 Commit of my ini file for rnaseq_light in cedar and modifications to the pipeline to allow for bootstraps in kallisto

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      36 commits

       eb4829440 update to last version of mugqic_tools 2.2.1
       4a71fcf57 add dnaseq_high_coverage.mp2b.ini file
       af718ac75 correct select_input_file call
       ea4250e96 correct select_input_file call
       5fbb45686 correct select_input_file call
       c1fcfd2af fix haplotype list parrelization
       84001d2ca fix haplotype list parrelization
       260a33630 correct select_input_file call
       950c98088 fix haplotype list parrelization
       bd92b74e2 correct select_input_file call
       6866c66af fix haplotype list parrelization
       a4427e4d5 fix number of haplotype and bed compatibilty
       32f78b6bb fix number of haplotype and bed compatibilty
       c9bc68a0d correct select_input_file call
       9a4a376e7 remove file format error
       35b3dcb4e correct select_input_file call
       e6ef88196 include both gatk and gatk4 lib import
       0a34bc750 correct bgzip_tabix libnrary call
       068983ea3 correct flash log input file
       6a2b07011 correct Mutect2 new arguments
       4da201af1 correct realigned bam name changed in inherited function
       bccc5cb83 add sequence_dictionary_variant proprety in dnaseq
       dbcd5a434 correct typo in DNAseq base ini
       12dbf3ed1 remove typo in deilverable command
       dc11a6fba remove verifyBamId for methylSeq pipeline
       637606eb7 correct typo in the tabix_bgzip call
       266752327 correct a typo
       de9c8b3e3 add md5sum generation
       d905e598d add required=F for markDup other_option
       908ffc36b change htslib.bgzip_tabix_vcf to htslib.bgzip_tabix as thereis no _vcf function in the lib
       c98c362cb make the add UMI steps skip automatically if no UMI in  the readset file
       043ba6542 add known_mills option
       40e65362d add required=F for markDup other_option
       51428804c primers addition bug correction
       c859b1a2b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       79ce228f2 add other_option paramter to compair in order to support other reference than GRCh37

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       c084a8eca Merged in slurm_report (pull request #72)
       4258a4538 Merged in new_style_class (pull request #63)

  P-O Quirion <pioliqui@gmail.com>      6 commits

       899d587d0 parse exit status, and only 50 lines after pro/epilogue
       a50700c8d log error when output logs have the wrong job number
       6b28d34ee log report for slurm
       a6c55532e id setter for job
       01c96cfd1 All class are new style, add setter to some getter
       8bf27ae9d working on Job

  Rola Dali <rola.dali@mail.mcgill.ca>      3 commits

       b72b06b3f Merged in beluga_inis (pull request #74)
       6077a5880 adding draft beluga inis
       258876e90 README.md edited online with Bitbucket

3.1.3        Tue Dec 18 16:37:25 2018 -0500        108 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      52 commits

       e239ca0f6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       6cd0351cf Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       43c5f25be Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1b3f70c5c Software installation scripts : new ones (mainly for Tumor SV pipeline) and updated ones Genome installation scripts : updated Homo_sapiens.GRCh37.sh (uncommented some useful commands) & Homo_sapiens.hg19.sh (updated vcf annotation process)
       8d3a636b2 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0a312f316 updated version of mugqic_tools
       395c7dcda Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       1506eac46 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       84fa5170b MethylSeq - DMR - adusted 'other+options' for dmr step
       e6a306d6d new addings to AMP_Scanner installation script
       5e5094f17 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       fd32bb542 MethylSeq pipeline - updated base.ini file with MethylKit differential analysis step requirements
       26d842c69 MethylSeq pipeline - added MethylKit differential analysis steps & cleaned the code
       e25b313ba added methylkit.py as a bfx wrapper to call methylKit.R from mugqic_tools
       5869184d4 debugged 'filter_snp_cpg' & 'prepare_methylkit' subroutines used by MethylSeq pipeline
       0c798db11 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       f6047c44f updated VarScan installation script with lates VarScan version v2.4.3
       efa28c889 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       869f99cf4 updated Conpair installation script to latest Conpair version v0.2
       7e2942135 added Manta installation script
       2233a27f9 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ad29c1503 added Ruby, LAST & Picky installation scripts
       f085968be Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       0683e0469 updated AMP_Scanner.sh installation script
       498f19032 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       41a042ae8 add AMP_Scanner installation script
       2de2537ea updated regexp in install_module.sh
       934fdec2b update supernova.sh with latest version 2.1.1
       6a82f17f4 update cellranger.sh with latest version 3.0.0
       82783dbc5 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       36e6b1c3e Merge branch 'master' of bitbucket.org:mugqic/genpipes
       bad579250 CCCG-1143 - added Caenorhabditis_elegans.ce11 genome installation script
       703383837 BFXDEV-578 - updated FasTree installatino script
       c35f7b27e BFXDEV-756 - changed the '-j' pipeline parameter default value to 'slurm' instead of 'pbs'
       4847e413c added StringTie installation script
       6cc4a37a4 added RNA-Seq Light section in the README
       90ea2232f Version bump to 3.1.3-beta
       0b4dcefab Version bump to 3.1.2
       2e259ec69 MethylSeq - updating ini files after testings on guillimin
       19568c7ca Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       cbefa4177 MethylSeq - cleaned methylkit_dmr subroutine
       8be6e376d MethylSeq - added MethylKit DMR steps
       2ef5d07ee MethylSeq - adding the bfx wrappers for DMR steps
       35aac401a Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e6d9edf0e Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       685e8a519 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e7f4335d4 fixing ini for trimmomatic16S and updating mammouth ini for dada2 protocol steps
       b716290bf AmpliconSeq : resolving asva step dependencies
       57f08b852 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e8af44bdc AmpliconSeq - updates to Dada2 protocol
       523d1da66 AmpliconSeq - updated last step of dada2 protocol for better dependencies and coding standards
       8a9605ed5 BFXDEV-674 - updates : flagstats calls are now made after alignment and deduplication instead of during the metrics step

  Édouard Henrion <henrione@cedar1.cedar.computecanada.ca>      2 commits

       7670b7e1f ChIP-Seq pipeline - correcting cedar.ini file for missing 'homer_make_ucsc_file' step requirements
       b66e93544 AmpliconSeq - call 'zless' instead of 'less' to avoid issues on Graham and Cedar systems

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       fcafa26a9 Merge branch 'dada2' of bitbucket.org:mugqic/genpipes into dada2
       fa5e298b8 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       a92edf35c AmpliconSeq - adding cedar ini file

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      6 commits

       d6515ace9 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       c3180f272 update database references in ampliconseq.base.ini
       37685bd41 asva.py pipeline merged to ampliconseq as another protocol asva.R deleted because moved to mugqic_tools new trimmomatic16S function added in bfx/trimmomatic.py
       518105081 update nanuq2mugqic_pipelines.py to also fetch the primer sequences ; needed by the new AmpliconSeq protocol (dada2)
       069489991 removed cutPrimer and set sys.path correctly
       dfe15133d Creating dada2 branch content

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      4 commits

       c54fe2b8e added pool parameter to bfx/dada2.py & pipelines/ampliconseq/ampliconseq.base.ini
       94f9800f4 removing deprecated pipelines/ASVA/asva.py
       9b09c6644 first working verson of dada2 protocol for ampliconseq pipeline
       5acf52721 minor updates (line break & indentation)

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       5627b8e68 GATK4 Updates - gatk4.ini - removed mugqic_dev module
       0f32dbfe7 Updates GATK4 - dnaseq.cedar.ini - removed mugqic_dev module
       75b915317 Updatres GATK4 - dnaseq.base.ini - remove mugqic_dev module
       dfda31ebd MethylSeq UMI - methylseq.base.ini - updated mugqic_tools to mugqic/mugqic_tools/2.2.0 for fgbio tools
       7c9c703aa MethylSeq UMI - methylseq.base.ini - updated module_fgbio to cvmfs version of the module
       b29b1171a methylseq.py - commented subroutine all_sample_metrics_report as it has been remove from the pipeline (because useless)
       b8f3fb62c fgbio.py - correted typo in addumi surbroutine

  Emmanuel Gonzalez <emmanuel.gonzalez@mcgill.ca>      1 commits

       223052d2a Merged in dada2 (pull request #49)

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       e1c7ef875 README.md edited online with Bitbucket
       6d1bd5cc6 README.md edited online with Bitbucket

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a6fe363b0 remove all_sample_metrics_report as it does the same as ihec_metrics
       b446c8667 change the ini section for module in tools.methylseq_ihec_metrics_report
       598028bc5 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       cee046448 adding umi and on-target metrics
       96f4decde adding the missing count parameter -c to samtools.count
       ba6545f3b correct typos
       0d03902ba correct small bugs for UMI integration
       317fe6a3c Merge branch 'methylSeq_UMI' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       acb55363a corect typo

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      11 commits

       1d2a7d321 add on target metrics and UMI
       86cea7062 add target_cpg_profile function in metrics lib
       ad2dd57e7 add samtools mapped_count function
       b322c94e3 add samtools mapped_count function
       9d49c2a6b correst typo
       6ccf4ee41 integrate UMI step in ini file
       cf4c008f6 include bam UMI annotation step
       f74f1ae8f include UMI annotated bam to the select input file of the picard_merge_sam_files step
       5c52b3dd0 Add other_option  system to MarkDuplicate
       2efad2206 add the UMI field in the readset file
       9ef861066 add the UMI field in the readset file

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       7f369e1c6 Merged in methylSeq_UMI (pull request #52)

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      1 commits

       cf47de1d7 updates to cedar ini

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      2 commits

       5b912aa9d Cedar resource and gatk4 metric fixes
       672ad3308 Updates to cedar.ini

  Robert Eveleigh <eveleigh@ip16.m>      2 commits

       417782216 GATK4 mp2b file added and improvements to alt contig exclusions
       91fb34960 Added mp2b ini and exome improvements

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       e9c66c64d Updates GATK4 and annotations

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       3c0ac3118 Resolve conflicts
       71b749618 improvements to exome handling

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       ad3caf1ed Merged in dnaseq_gatk4 (pull request #47)

  Rola Dali <rola.dali@mail.mcgill.ca>      1 commits

       dc7996dd8 README.md edited online with Bitbucket

3.1.2        Wed Nov 21 15:05:01 2018 -0500        30 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      17 commits

       be396f04b updated mugqic_tools version 2.1.12
       f187cdf17 removed _dev modules from ini files
       2bfc30777 Updated R_Bioconductor.sh : pointing to system libraries improved
       e80f72961 PacBio Assembly pipeline - corrected incompatibility bug with Guillimin
       1be0124f9 updated install_module.sh to accomodate both CentOS7 & Ubuntu16.04 system libraries
       ff709f1fc update HiCUP install script with latest version 0.7.0
       659fbc958 corrected typo in ampliconseq.py causing pipeline crach at plot heatmap
       c8b5634ed updated README-release.txt with correct URL for GenPipes download page
       09fe36979 Changing resources requirements for [gatk_merge_and_call_individual_gvcfs] in DNA=-Seq pipeline
       30a9d33a1 corrected reference to kallisto tx2gene file in the RNA-Seq light pipeline
       0e7dae6fe corrected reference to kallisto index in the RNA-Seq light pipeline
       e8ec0ac47 Corrected bug in HicSeq pipeline when running in capture mode
       58841e5ca updated RNASeq light base ini file with kallisto version on CVMFS
       78a4d4b34 corrected core/pipeline.py to avoid pipeline erroring when using '--json' parameter...
       1e2ef94d2 corrected typo in R_Bioconductor.sh
       e4f493fe5 Version bump to 3.1.2-beta
       c6b48be69 Version bump to 3.1.1

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      4 commits

       929e37c24 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f4479be0e freebayes install script
       9fadba98d platypus install script
       7909a4e37 vcfanno fgbio & delly install scripts

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       362f10166 Merged in hicup_arima (pull request #51)
       f3593a675 tumor_pair.base.ini : changed gemini version to 0.20.1
       31d2d00e5 tumor_pair.guillimin.ini : removed some more _dev modules
       b689763a3 tumor_pair.guillimin.ini : removed mugqic_dev module
       6e9b3f995 hicseq.base.ini : updated HiCUP version to 0.7.0
       07b9a21e1 dnaseq_high_coverage.base.ini : setting ram parameter witin section igvtools_compute_tdf
       d21089324 Tumor_pair pipeline : bug fix in bfx/bcbio_variation_recall.py Correted typo in the executable call

  Rola Dali <rola.dali@mail.mcgill.ca>      2 commits

       6d9447418 editing hicseq.py for Arima compatibility
       33d597d3c adding HiC Arima digest to install_genome

3.1.1        Thu Nov 1 15:32:25 2018 -0400        161 commits

  David Bujold <david.bujold@mail.mcgill.ca>      1 commits

       b0adf94dd Merged in pipeline_stats (pull request #20)

  dbujold <david.bujold@mail.mcgill.ca>      2 commits

       86c8a724c Display JSON log statistics into tables and figures on the log VM.
       6f716cb93 Python CGI script to tranform pipelines stats log file into a JSON document.

  Edouard Henrion <edouard.henrion@mcgill.ca>      32 commits

       24957210a Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b14189dfd update install_module.sh script with integration of apt along yum as system libraries resources
       0de17ddfd update R installation script with integration of apt along yum as system libraries resources
       7a57dfe73 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       de73f87e5 update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       cd6f62afd update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       052e9fcbe Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e8af086e2 updating genome ini file generation with versioning
       e2f53d808 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b92974218 update smrtlink.sh with latest SMRTLink version
       8130d0594 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       13080e22a updated R_Bioconductor.sh with new packages and corrections
       63e12b542 new version of multiqc.sh to install latest version
       f499c1b91 aded popoolation2 installation script
       211dc6225 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       d58b1c3cb updated gatk.sh so it can handle both versions 3 & 4 installation
       4b5565f63 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       c4004ed5f updated R_Bioconductor.sh - new libraries & patches
       38e01cdef updated python.sh : make -j12 & configure command
       b06c41080 resources/modules/install_module.sh
       d1fe8bf66 update flash.sh with parallele make -j12
       6eca5f3e1 watch_portal_folder: removed sample_name in filename
       f0e62220d Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1ef165c91 corrected gatk.sh : paths in the modulefile were wrong
       749cd972c updated rnaseq.py : report steps are now included in the analysis JSON file
       8978a5925 updated ampiconseq.py : both steps merge_flash_stats & merge_uchime_stats are now inclded in the analysis JSON file
       b4d4b046f updated rnaseq_light.py : samples added to each job including report jobs, reviewed indentations and spacing
       8e75e026f updated rnaseq_denovo_assembly.py : samples added to each job (some were still remaining) including report jobs
       61af5bf26 updated common.py : added samples to the call of rmarkdown.render
       cc974a38a updated bfx rmarkdown wrapper : JSON file generation added
       1d9b8347b Version bump to 3.1.1-beta
       bd721f1d8 Version bump to 3.1.0

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      1 commits

       b3e456011 removing all the remaining MUGQIC_INSTALL_HOME_DEV in all the ini files

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      6 commits

       4a4d0e572 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ca805a0f0 updated install_module.sh : new installation procedure integration, adding the patching the C3G executatbles with patchelf : making sure all system libraries are now searched in /cvmfs/soft.mugqic/yum/centos7/1.0
       f771e9bc6 updated python.sh with some new package : umap-learn
       63a855eb2 udpated R_Bioconductor.sh with new packages. Also, now integrates the new installation procedures including pathcing the C3G executables (i.e. use of patchelf)
       41fcad8cf Merge branch 'master' of bitbucket.org:mugqic/genpipes
       4a5958c0f updated archive URL in flash.sh install script

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      7 commits

       c243ec8f7 updated R_Bioconductor.sh script : set the PAGER variable to /usr/bin/less
       f315771dc updated pipeline READMEs : all the steps are now shown independently of the available pipeline protocols
       5fa7bb431 updated R_Bioconductor.sh with some new packages in the install list
       1c468e09f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b24d7f89c updated silva.sh to install latest vesion of silva DB
       04140e3a9 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e9fcbff4b added Glycine_max.v2.sh installatino script for Glycine (Soybean) genome installation

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      16 commits

       bdf93cc7f updated hicup.sh
       795b9e5b1 updated smrtlink.sh with the version 6.0.0 of SMRTLink
       0a9fb4b7b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e141ec79a updated pipeline cedar ini files to avoid using java from CC software stack
       efbf07ab5 modified core/pipeline.py to avoid having to set  when not generating the anaylsis JSON file
       256907f30 modified DNA-Seq README with better step descriptions
       76a26cb30 Updated the pipeline workflow diagram download links : path of the full-size picture instead of the resized one
       d0fd87451 Added a link to download the pipeline workflow diagram along with the diagram picture itself
       11cbf4a83 Updated pipeline README.md files, with workflow diagram pictures embeded
       30c53cc5f deleted pipelines/tumor_pair.base.exome.ini
       bd6f0000b moving tumor_pair.base.exome.ini from 'pipelines/' to 'pipelines/tumor_pair/'
       9362d720b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       baf53baec updated ampliconseq pipeline following Emmanuel gonzalez comments
       6a85b8fb7 added source EnsemblPlants to install_genome.sh script
       690f20b3b added skewer installation script
       59996b55d BFXDEV-591 - updated ampliconseq.base.ini based on Emmanuel Gonzalez feedback

  ehenrion <edouard.henrion@computationalgenomics.ca>      2 commits

       3b3ff687f updated mugqic_pipelines.sh so that it now refers to genpipes repository on bitbucket
       97a4552c6 corrected typo within gatk.sh installation script

  ehenrion <edouard.henrion@mcgill.ca>      18 commits

       12388828f methylseq.py : corrected dependencies assignment for ihect_sample_metrics_report step
       85bfc1853 rnaseq_light.py : typo correction
       2057f1e6e rnaseq_light.py : corrected typo
       91be7b895 rnaseq_light.py : corrected typo in job parameter assignement
       5de9c4e2e dnaseq.cedar.ini removed 'module_java=java/1.8.0_121' from cedar.ini
       b4ce9ef9e methylseq.base.ini : modified [bismark_align] section : maximum insert_size now set to 1000
       146eaa735 methylseq.mammouth.ini : modified bismark align section
       a4eea15a3 methylseq.cedar.ini : modified bismark_align walltime
       a22490f35 methylseq.base.ini modified [bismark_methyl_call] section within the base.ini
       a5276a3da methylseq.base.ini - modified bismark_align parameters and resources within the base.ini file
       4fd03a107 dnaseq.py edited online with Bitbucket Correct typo at line 724 : "jobs" replaces "obs"
       425c086f1 dnaseq.cedar.ini : added variant_recalibrator section to defined resources on cedar
       d55362f02 Merged in revert-pr-41 (pull request #46)
       17aa18fd9 Revert "Can run on HPC with slurm and  containers (pull request #41)"
       60be1aa26 rnaseq.py edited online with Bitbucket added missing job (metrics.wigzip) to the json analysis file
       1c5f56dde README.md edited online with Bitbucket
       6dccf5b4b README.md edited online with Bitbucket
       3b3407be3 README.md edited online with Bitbucket

  Éloi Mercier <emercier@cedar5.cedar.computecanada.ca>      1 commits

       c6593cf5e in cedar.ini and graham ini of rnaseq, chipseq and dnaseq: change assembly_dir to MUGQIC_INSTALL; in dnaseq.graham.ini: uncomment assembly_dir variable

  emercier <eloi.mercier@mcgill.ca>      10 commits

       ec3183ff3 all pyc files removed
       f0568f47e Update R_module in ini files to mugqic/R_Bioconductor/3.5.0_3.7 (except for illumina_run_processing)
       5a14b4d81 in ampliconseq, pacbio and rnaseq.guillimin.ini: change lm queue (depreciated) to meta queue
       e13be161e in modules/weblogo.sh: remove whitespaces at the beginning of the echo blocks
       4675ad672 in dnaseq and rnaseq.base.ini: change R_Bioconductor to stable version 3.4.3_3.6
       fcad0258f in install_all_genome.sh: add Danio_rerio.GRCz11.sh
       d353ecf95 add install script for genome zebrafish Danio_rerio.GRCz11
       941615658 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0b53f6455 in install_genome.sh: small fix to create_kallisto_index and create_transcripts2genes_file functions
       f5681c4d5 in install_genome.sh: fix bug in create_transcripts2genes_file to work with recent version of Ensembl

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      4 commits

       53d11f3e4 rnaseq.base.ini edited online with Bitbucket: commented explicit adapter file parameter
       8653d0272 R_Bioconductor.sh edited online with Bitbucket: added PopSV. Currently commented out since not tested
       0b42f96c6 R_Bioconductor.sh edited online with Bitbucket: Added a few dependencies to the list
       3a268e667 R_Bioconductor.sh edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       a23ba480f Added dev install scripts for delly, lumpy, sv, vcfanno

  José Héctor Gálvez López <hgalvez@ip16.m>      2 commits

       88e492937 Further refinements to Mp2b based on feedback from mammouth admins
       39bf7606c Added ini files tailored for Mp2b based on Cedar ini files for the following pipelines: RNA-seq, RNA-seq de novo, DNA-seq, Hi-C seq, Methylseq, and ChIP-seq.

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       1f95d6806 Add cedar ini file for methylSeq
       8b071760c Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b1ed77cd1 removing pipe empty of jobs in tumor_pair
       7e70dc572 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       db2095776 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9e0178b24 correct bug in ucsc bedGraphToBigWig which raise an error for specie with MT chromosome name except for GRCh37 - BFXDEV-737

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      7 commits

       2e021f676 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       85ead2b6a resolve issue with multiple inputs in DNAseq picard markduplicates
       1e00830b0 resolve pull conflict
       b0042a0e7 update dnaseq.cedar.in file
       4074ab5f8 adjust dnaseqto add select input file to some of the steps - need to be continued
       e38920e44 modify Slurm scheduler delay (sleep) from 0.5 to 0.2
       b1d3d4b5f ChipSeq - add mutliqc param in the cedar ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       9c7c8158d Merged in add_container (pull request #41)

  P-O Quirion <pioliqui@gmail.com>      5 commits

       cfa63b208 Merge branch 'master' into add_container
       44274bdbf Basic container ini file
       4296d3174 Working singularity version tested on graham
       b14afd93f WIP exec line prototype
       8b9241612 add container option

  Rola Dali <rola.dali@mail.mcgill.ca>      8 commits

       795776b90 methylseq.cedar.ini edited online with Bitbucket: added cluster_walltime to ini to avoid errors
       47addd28b chipseq.base.ini edited online with Bitbucket: bigwig and run_spp resources are not enough; edited them to avoid failure
       50747a731 Merged in IHEC_metrics (pull request #44)
       bdb0b76e0 dnaseq.base.ini edited online with Bitbucket: change threads for GATK due to errors
       142921b06 update csvToreadset.R in utils to use column names due to changes in nanuq csv
       40cdff0ff Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       8709b9d81 README.md edited online with Bitbucket
       13f08a5f5 README.md edited online with Bitbucket: -j slurm and tutorial

  Romain Grégoire <romgrk.cc@gmail.com>      3 commits

       6f8ea86f9 Merged in add-json-version (pull request #50)
       1271bc19e Merged in dashboard-display-launching-user (pull request #43)
       73e768a4f Merged in logs-add-checksum (pull request #42)

  Rom Grk <romgrk.cc@gmail.com>      28 commits

       ff11dbc52 jsonator.py: add version number
       bb72124f4 watch_portal_folder: read sample_name from filename but stay backward-compatible
       921e54df6 job2json: add sample_name to file for genpipes dashboard
       784848f61 copy sample json files during script run
       b414b0455 scheduler.py: clean unused value
       5a9d8afc5 Sample: remove .json_dump property
       7d45ce1ed lint pipeline.py
       337a1c8f3 lint job2json
       014c19d1c job2json: guard __main__
       26b9ef986 lint job2json
       bf4f976d8 job2json: use $USER of user running script
       8980701ed common.py: remove unused import
       622b61119 common.py: fix log issues
       676f7bee9 common.py: add unique md5 checksum to logs
       56545fecc Revert "watch_portal_folder: put sample_name in filename"
       a8badcd0f Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0526eccad watch_portal_folder: put sample_name in filename
       13b56384a Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0411ae98b watch_portal_folder: fix memory errors
       e5e357955 watch_portal_folder: add logging
       46d8a35dd watch_portal_folder: fix typo
       f54a1277f watch_portal_folder: add cache option
       7f89e828d Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       2fbdca883 watch_portal_folder: implement update by diff
       45cc1899b watch_portal_folder: update script
       67c14714c watch_portal_folder: skip missing files
       8a9dd08a3 watch_portal_folder: more resilient to network errors
       32bf92a9b watch_portal_folder.py: dont watch if no interval is provided

3.1.0        Wed Mar 28 15:46:33 2018 -0400        188 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      77 commits

       ff52e7c98 MethylSeq - added 'ram' parameter to the igvtools_compute_tdf step in methylseq.base.ini
       68661e8c9 updated jb2json.py with a better locking system : now creates a folder instead of a file in order to create the lock
       55bf1b0fc Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       43c9cc1d7 updated core/scheduler.py to remove job2json call when --json parameter is omitted
       21b3bdfdd updated gatk.sh with version gatk-4.0.2.1
       99087312d updated longranger.sh with LongRanger version 2.2.2
       d03abd52b updated hicseq.py : calling of the newly built bfx libraries and minor indentations changes
       c52b58206 reviewed topdom.py wrapper for a better handling of input and putput files, so that job dependencies do not break...
       8ae89bdf4 aded locking file system to jsonator.py and job2json.py to avoid multiple and synchronous writing attempts to JSON file
       4d1eaa680 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       70ddf870d updated hicseq.base.ini with newer version of mugqic_tools (2.1.11) and with revised resource allocation for hic_align step
       3f798270f removed useless loading of mugqic_tools module when calling job2json.py
       4d7da4579 added locking file system to job2json to avoid multiple and synchronous writing attempts to the JSON file
       6aafbf04c removed useless loading of mugqic_tools module when calling job2json.py
       205d577d8 mugqic_tools.sh : swith to mugqic_tools version 2.1.11
       0ea6336dd improved the process regarding the update of the analysis JSON file when resuming a pipeline execution
       62dd71d30 removed useless import from core/scheduler.py
       ac6e9ce66 added the locking file system to avoid multiple & simultaneous writing attempts on the same file
       ffd35a6df Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       e6545082d added new bash tools wrappers to be called by hicseq.py pipeline script + some minor indentation & syntaxe updates
       843d6573f added new bfx libraries to be called by hicseq.py pipeline script
       1f6be1b7a Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e1dd68f77 added some new software installation scripts
       e2abf45c1 updated 'core/job.py' : added 'samples' parameter to 'concat_jobs' and 'pipe_jobs'
       8cf84f4e5 added the initialisation of 'portal_output_dir' to '/lb/project/mugqic/analyste_dev/portal_output_dir' within all of the pipeline .base.ini files
       9bae17c10 added the generation of the analysis JSON file to the SLURM scheduler class
       6e6363d6d Merge branch 'master' of bitbucket.org:mugqic/genpipes into cedar
       87fb0142f updated core/scheduler.py : call of job2json.py is now done without /home/ehenrion/work/portal_repo/genpipes
       1681c6f6d Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       f95e69097 corrected trimmomatic input selection to reflect the file names what might be outputed from picard_sam_to_fastq
       f22ee01ea update verifyBamID subroutine within common.py so that it now accepts .vcf.gz files (formerly was only accepting .vcf files)
       998e1bf46 updated version of cellranger to 2.1.1 within the bash installation script
       ba5099f4a updated pipeline READMEs
       b5cecbd08 updated core/pipeline.py to make the analysis JSON file creation optional : no JSON file created by default
       3fae6e86c updated bfx/macs2.py so that it follows the C3G coding standards
       3b426eaf1 added installation script for hdf5 and zlib libraries
       2e82f0795 corrected typo in install_genome.sh
       ff59b0a4f updated install_genome.sh script with newer software version as well as minor indentation fixes
       3d26798f9 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1bd0b9883 updated kallisto installation script
       cf3219ad3 updated STAR version in the installation script
       19997fafa updated dbSNP version to 150 for Homo_sapiens GRCh37 installation script
       1e6f03200 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f89e056f2 updated kallisto version to 0.44.0 in the genome installation script - corrected typo within create_transcripts2genes_file subroutine
       7bf13a920 updated dbSNP and dbNSFP in Homo_sapiens genome installation scipts
       30299e2bd added cDNA fasta to the genome installation scripts
       e17fcb1b8 added the portal_ouput_dir to in base ini file
       acea9bf90 Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       fe6162758 Merge branch 'portal-integration' of bitbucket.org:mugqic/genpipes into portal-integration
       e0010e280 Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       1f63dc5d3 Added README-GenAP_coding_standards.txt which gives the coding guidelines for who may want to participate in the C3G/GenAP developments
       ed345d284 BFXDEV-721 - updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       9a59536b1 updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       0ed0ef60a updated Prokka installation script so that 'prokka --setupdb' is launched afer the installation
       743668f16 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       22b260132 updated Supernova bash installation script to version 2.0.0
       c62cf42e5 updated the bash installation script of python : installation of QIIME has been updated so that emperor is well installed before QIIME
       96c98a443 added the initialization of the LFS variable, otherwise it is recognize on guillimin...
       b02491e05 updated the bash installation script of python : after the installation of the specified version of Python is successfull, the script now also takes care of the installation of all the needed python libraries
       7609147ad update type within python_lib.sh script
       1ab7a7ea5 RNASeq pipeline : value of parameter 'localfit' used by DESeq, is set to default i.e. empty (which also means 'false'), so Parametric Dispertion Fit is performed by default instead of Local Fit
       896c3b2bf DNASeq pipeline : updated versions of mugqic_tools (to 2.1.10) and of samtools (to 1.4.1) within the base.ini file
       74a5c9ad3 updated picard installation script : now installs version 2.17.3
       8d2dde0ce added the installation script of Prokka, a tool for prokaryotic genome annotation
       cbfe9d755 BFXDEV-673 - remove some verbose during the execution of job2json.py
       66893f89b updated mugqic_tools.sh so it now installs version 2.1.10
       170f39bf6 added the bash installation script for the Illumina InterOp parser
       dc0ebcb0a BFXDEV-674 - MethylSeq pipeline - minor updates within the metrics .md report template file, for standardization purpose
       d9e7b0b4a updated ortograph within sample metrics .md report files
       e0b41ff94 BFXDEV-673 - updated scheduler.py to generalize the use of job2json to all schedulers
       5250d807d BFXDEV-674 - added a report .md file for IHEC metrics reports for targeted anaylsis
       cb466e700 BFXDEV-674 - updated .md files for metrics reporting : revised headers & descriptions
       460f24a67 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       2833f6a95 updated gemini.sh bash installation script : added PYTHONHOME setting when loading the gemini modules
       f483e37e4 update release instructions to generate proper README for all the pipelines including HicSeq
       2d70be4d0 Version bump to 3.0.1-beta
       3cb8610ee Version bump to 3.0.0 - updated

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      2 commits

       89831808e Merge branch 'master' of bitbucket.org:mugqic/genpipes
       49c6b3204 updated kallisto bash installation script

  ehenrion <edouard.henrion@computationalgenomics.ca>      3 commits

       3f60795b1 updated R_Bioconductor.sh with good indentation and new packages installation
       b1215c322 corrected wrong environment variable name within blast.sh
       dc4fbb7b6 updated the version to v359 within ucsc.sh script

  ehenrion <edouard.henrion@mcgill.ca>      6 commits

       39fea264b rnaseq_denovo_assembly.base.ini edited online with Bitbucket
       f9278f904 README.md edited : "MUGQIC pipelines" replaced by "GenPipes"
       e6abfd90d apis_mellifera.sh deleted : was the exact replicate of Apis_mellifera.sh
       8e9ab1cdf README.md edited Updated links to the pipeline pages
       469b9b482 README.md edited online with Bitbucket updated some links to reflect the repository renaming to genpipes
       87524a9d3 smrtanalysis.py - standardized command-line format and indentation

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      3 commits

       99b1a5dcc removing pyc files from rnaseq
       1c7731a52 in rnaseq.mammouth.ini add section for bed_graph to set up ppn to 1
       6d910887d In nanuq2mugqic: change readset file name to readset_<seq_type>.tsv

  Eric Fournier <ericfournier2@yahoo.ca>      2 commits

       6fafc8d07 Merged in ericfournier2/genpipes (pull request #35)
       4138a7289 Fix list within list bug which breaks chipseq pipeline.

  Gary Leveque <gary.leveque@gmail.com>      1 commits

       685730766 Merged in gary_pacbio_assembly (pull request #29)

  gary.leveque@mail.mcgill.ca <gleveque@abacus2.ferrier.genome.mcgill.ca>      7 commits

       fd1e91cc6 additions made to smrtanalysis.py for basemodification and motifMaker steps
       b88dd9e8d Addressed issues commented by Edouard; tested on abacus and mammouth
       39f57bce1 added cluster_server= to pacbio_assembly.mammouth.ini
       b135ae0bd revised versions of .base and mammouth.ini files; I was changing between sw and lm nodes
       6bac3683e revision of pacbio_assembly.mammouth.ini, back to qwork
       356aadc0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into gary_pacbio_assembly
       a24451f22 Addition of base modification detection and generation of a motif_summary.csv steps to the pacbio HGAP assembly pipeline;  see BFXDEV-703

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       2455b6db7 Merged in hector-cedar (pull request #37)

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      1 commits

       36b892fa2 Corrected minor bug in the create_scheduler() function that was creating errors when using slurm

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       05bbcec55 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3b818dc9d externalize excluded chormosome to be specified in the genome ini
       a2d807b97 remove conflict
       3896d7729 uniformize verify_bam_id vcf path
       f4bf8cc75 add a delay (sleep) after the symlink creation to avoid issue (Invalid Job Dependency) during submission
       322348df5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f734dfb45 test new url
       37529400c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0834b95f add NovaSeq to nanuq2mugqic script
       0dd2ba3e1 DNAseq - split the pipeline into 2 different pipeline GATK best_practice or old mpileup - BFXDEV-716
       5c586fffd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9d48a589e chnage RNAseq_denovo inheritance to RNAseq
       8a250c814 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7716f961a include a test of inputs type for some library function which does not enforce the list type

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      24 commits

       e2b2a7ab4 resolve conflict with master
       82d1bee61 Merge branch 'cedar' of bitbucket.org:mugqic/genpipes into cedar
       443fcf649 restore low memory ini now the slurm bug is corrected
       78d156b38 update rnaseq for cedar
       62520dad0 update rnaseqDN on cedar
       9a2740231 update hicseq on cedar
       b9842b9da update dnaseq ini on cedar
       b695bcfea add more more RAM to compensate slurm bug
       924fff54e add 0.5s sleep to let slurm submiting the job correctly
       282fa9eb3 Increase memory request to avoid I/O buffering hit the memory limit
       5f157cd32 Changing I/O block size
       3b395c5e9 Changing I/O block size
       614564d5b remove confilct pulling master
       be8157940 Modifying I/O block size for cedar
       496e10d2e update ini file
       7f77b1953 update log with slurm scheduler
       818060988 add chipseq ini file for cedar
       0e279092c update RNAseq and link to mugqic_dev
       d26fef63a update DNAseq and link to mugqic_dev
       10c20693c remove lattency for jobsubimission to slurm
       d192a9152 remove conflict
       77e5d3596 update cedar scheduler
       67b80063d update cedar ini
       120f92d3f DNAseq -ini file for cedar

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       804ae6cf6 Merged in cedar (pull request #34)

  Mathieu Bourgey <mbourgey@cedar5.cedar.computecanada.ca>      9 commits

       9a9ab773d DNAseq - update ini file
       910be7692 remove scheduler module testing
       9e416b180 change igvtools exe to igvtools jar in order to have control of the ram usage
       9d3c7bdca adjust RNA and DNA ini files;  generate fake prologue and epilogue; add 2s delay between each job submission- BFXDEV-683
       95f20c639 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       f247a6a92 create RNAseq ini file - BFXDEV-683
       6f5cca53f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       297a2c6bc create a DNAseq ini for cedar; remove typo in base DNaseq ini  - BFXDEV-683
       1dce04b85 create a version of the scheuler class for slurm  - BFXDEV-683

  Mathieu Bourgey <mbourgey@gra-login4.graham.sharcnet>      1 commits

       fea150101 Add Graham ini for RNA and DNA

  Rola Dali <rola.dali@mail.mcgill.ca>      17 commits

       6978e141c Merged in IHEC_metrics (pull request #33)
       d66737bcb added .hic file generation to capture hic
       3a9b9a382 adding RobusTAD TAD scoring to hicseq.py
       3c7972102 adding multiqc to chipseq. should customize yaml file and check homer module
       ad0ff99a2 changing rnaseq dependencies to ensure no repeats when job is complete
       e581c3ab7 adding mugqicValidator.py to utils to validate basic structure of readset and design file
       8ee10839c removing job outputs from chipseq ihec metric method to accomodate samples running without input
       befab29e1 Merged in IHEC_metrics (pull request #31)
       2a0d785c5 fixing the tmp_dir for macs
       6a363de2f Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       a10ad5188 changing genome to fit with merged homer.maketagdir and changing ihec matrics output name
       bba10e0ba Merged in IHEC_metrics (pull request #28)
       28850c508 config.py edited online with Bitbucket
       00b9e34e0 adding query_module=spider to mammouth ini
       00a746faa allowing module spider on mammouth to reduce module loading time--committing to test on guillimin and abacus
       e3c9fabbf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fea082b7f fixing MACS bug BFXDEV-712

  romain.gregoire@mcgill.ca <rgregoir@abacus2.ferrier.genome.mcgill.ca>      2 commits

       53fe3b919 pass config files directly to job2json.py
       94d21ed8a Merge branch 'portal-integration' of https://bitbucket.org/mugqic/genpipes into portal-integration

  Romain Gregoire <romgrk.cc@gmail.com>      1 commits

       655f7c93e watch_portal_folder.py: fix --username argument

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       2bcae74bb Merged in portal-integration (pull request #32)

  Rom Grk <romgrk.cc@gmail.com>      15 commits

       c3a008098 watch_portal_folder.py: remove username option
       f3c3f34e2 send $USER to portal integration
       bdb24f8ad watch_portal_folder.py: fix file path
       4981525b2 watch_portal_folder.py: sort files sent by modification time
       e1f361a84 export CONFIG_FILES to allow loading config from pipeline
       430d0c487 export CONFIG_FILES to allow loading config from pipeline
       446122c61 job2json.py: add mugqic dir to python path
       c35b6f93f watch_portal_folder.py: safer handling of response
       13196bfde watch_portal_folder.py: fix data posting
       bbbfe8a20 watch_portal_folder.py: fix script exit
       bff9c6495 watch_portal_folder.py: fix argument passing
       47b1073dc use uuid to avoid collisions when buffering JSON files
       7f99bdc26 job2json.py: make a copy of the updated JSON files for the portal
       629cef619 add watch_portal_folder.py
       63e4de625 jsonator.py: make a copy of the JSON files to a buffer folder

3.0.0        Thu Dec 7 14:19:49 2017 -0500        444 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      247 commits

       8c5a5c61f MethylSeq pipeline - BFXDEV-674 - updated wiggle_tracks step with more comprehensive output file names
       2024838e6 updated ucsc.py with simplified if-else statement for more clarity and corrected the 'chr' prefixing behavior
       bcf3d8a6a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       859428810 updated Homo_sapiens.GRCh38.sh installation script with vcf indexes
       cb43d498a updated jellyfish installation script with version 2.2.3
       9e3c955aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7e14f0918 Updated README for all the pipelines
       a6b244e4c updated mugqic_pipelines.sh script with the most recent version of the mugqic_pipelines, now called GenAP_Pipes, version 3.0.0
       86092c6d6 Version bump to 3.0.1-beta
       b8f43102b Version bump to 3.0.0
       3844eee84 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c9ad36803 added tagInfo.txt in homer.py
       33702603a DNA-Seq pipeline - update base.ini file by removing all the paths which were still pointing to a '_dev' location
       7b15520fa BFXDEV-674 - MethylSeq pipeline - updated bedtools.py intersect function to carry header from input bed to output bed
       da1ece6b5 added tagInfo.txt
       1b84ec6d3 Version bump to 3.1.0-beta
       03f3dda38 Version bump to 3.0.0
       5d0a77490 added the README file for RNA-Seq Light Pipeline
       837b4fa93 slightly updated release instructions
       1343ab644 Version bump to 3.0.1-beta
       188037f91 Version bump to 3.1.0-beta
       03638bd37 Version bump to 3.0.0
       4c027b6f8 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8768a7788 updated version to 2.0.2 for cellranger in the bash installation script
       00eb96b9a RNA-Seq pipeline - updated .base.ini file : removed specific module loading for differential expression step + minor indentation updates within trinity.py for RNA-Seq de-novo Assembly pipeline
       a4ab69750 RNA-Seq pipeline - update .base.ini file with newer software versions
       626856f98 DNA-Seq pipeline - updated base.ini with GATK version 3.7
       a7b56ba17 updated base.ini with GATK version 3.7
       dde4f146c master VERSION changed to 3.0.0-beta
       bfdcbae60 PacBio Assembly Pipeline - updated circlator step : created a circlator.py within the bfx folder and review the coding/calling of the circlator step within the pipeline python wrapper
       132f2f249 BFXDEV-674 - updated MethylSeq pipeline to correct wrong dependencies causing some jobs to always be consider as 'NOT up to date'
       8cbbff9d4 BFXDEV-673 - corrected error in bfx/jsonator.py occuring when modifying the list of softwares from an already existing JSON file
       942804617 some more minor updates on bfx/tools.py regarding standard indentation and parameters naming
       1aaf4f94d MethylSeq - IHEC metric report jobs are now labelled with the sample name
       29e74ac58 minor updates on bfx/tools.py especially to make indentation uniform across the whole file
       573dba592 updated & added software install scripts
       c8fa5cd49 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d964b2b8 commit prior merging to avoid conflict
       2baf98265 AmpliconSeq pipeline - adding mammouth walltime for qiime_otu_assigning step
       fb8c75db4 resolving conflicts in bfx/tools.py
       fdddddb7a BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq pipeline
       c961e1580 corrected typo introduced after resolving conflicts...
       f4116b0d0 BFXDEV-673 - updated jsonator.py to handle Illumina Run Processing ini file entries when generating JSON analysis file for Illumina Run Processing pipeline
       d573722c9 BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       fb346c9da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       763aac032 BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       8b4012908 BFXDEV-673 - adding analysis JSON file generation to the Tumor Pair pipeline
       965325051 BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq De Novo Assembly pipeline
       f961ef10c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1dcae3f78 BFXDEV-673 - added the analysis JSON file generation to the RnaSeq pipeline
       f107e0490 BFXDEV-673 - added the analysis JSON file generation to the PacBio Assembly pipeline
       b806ac674 BFXDEV-673 - updated jsonator.py to handle PacBio Assembly readset files when generating JSON analysis file for PacBio Assembly pipeline
       fdd8577d9 BFXDEV-673 - added JSON analysis file generation to the DnaSeq high Coverage Pipeline
       d48d64242 BFXDEV-673 - adding JSON analysis file generation to DnaSeq pipeline
       b8824eff0 BFXDEV-673 - adding JSON analysis file generation ito ChipSeq pipeline
       d19b3044e BFXDEV-673 - updated AmpliconSeq pipeline with analysis JSON file generation
       102f66cb5 BFXDEV-673 - updating jsonator.py to generalized the way dbsnp_version and server entries are handled
       054e41093 BFXDEV-674 - MethylSeq pipeline - adding samtools to the loaded modules for methylseq_metrics_report
       0169e5096 BFXDEV-673 - minor update to the help content
       18e8e6893 BFXDEV-674 - updated wiggle tracks generation tools by spliting bedgraph and wiggle traks into 2 different jobs, also managing .bw file for GRCh37 build to make it UCSC compatible
       cc7af63b3 updated 'tmp_dir' in guillimin .ini files : now using  space which is automatically cleaned up after the job ends
       04d9fed3d BFXDEV-674 - updated methylseq_ihec_metrics_report to reflect changes in mugqic_tools IHEC_methylseq_metrics.sh
       6e854c2a3 Picard mark_duplicate minor update
       8958faef7 BFXDEV-673 - corrected error in variable assignment within jsonator.py
       487227639 BFXDEV-673 - jasonataor.py now handles the case where user has used the 'force' argument for his pipeline execution in order to entirely re-create the analysis JSON file
       22c477d9f BFXDEV-673 - review job2json.py passed parameters to handle start & end job date/time
       61225eeb9 BFXDEV-673 - JSON analysis file key 'hpc_center' changed to 'server'
       4b2da7f93 BFXDEV-674 - MethylSeq pipline - updated mugqic_tools from dev to prod within pipelines/methylseq/methylseq.base.ini
       3e8afc373 BFXDEV-673 - updated core/scheduler.py to handle both job_start_date and jobs_end_date for the json analysis file
       c777a6479 BFXDEV-673 - updated core/pipeline.py for a better generation of the json analysis file
       1d735a0e2 BFXDEV-673 - updated bfx/jsonator.py
       1b0b956f8 BFXDEV-673 - updated utils/job2json.py
       98b0d7f5e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       625ca4fc1 added optparse library installation in R_Bioconductor.sh
       f751eb507 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c41c514b4 removing MULTIQC from mugqic tools and updating version to 2.1.9
       00118d611 version 1.0.5 of R_mugqic_packages
       18ab8513b BFXDEV-674 - removed useless bedtools parameters from mammouth and guillimin ini files
       a17600f99 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       b98dcf91e adding some new bash install scripts
       89a48fdf4 add chicago module
       0c6d44a0d fix conflicts
       9919d1a7b keep locale module change on libuser version of the repo
       6a3eff86e keep locale module change on libuser version of the repo
       2b55b3fdc BFXDEV-674 - corrected input flagstat file name for metrics steps
       9ac6e8134 BFXDEV-674 - adding on_target context handling
       69eb2bb39 BFXDEV-673 - corrected log file assignment to the JSON file
       7dd7c0de7 BFXDEV-673 - BFXDEV-674 - updated syntax and added IHEC metrics report
       5fccb8b46 importing samtools from common.py
       18e87a581 BFXDEV-668 - added ihec_metrics step to rnaseq.py
       167e9c305 BFXDEV-674 - removed useless parameter settings from ini file
       f7e2c9468 BFXDEV-673 - updated cluster_hpc_center parameter within the ini files
       5ca6c5213 removed unnecessary trailing comma dnaseq.py
       62c711a8b removed unnecessary parameters from dnaseq.mammouth.ini
       e3186d4eb BFXDEV-674 - corrected call of the bash script which generate the metrisc for IHEC
       10c042ddc BFXDEV-673 - updated jsonator for a better behavior of json file updates
       bfc0042d7 updated syntax for unrequired paramaters in bedtools.py
       8e42663cd BFXDEV-674 - updated bedtools with proper syntax standards
       f9538e808 BFXDEV-668 - BFXDEV-675 - updated bedtools graph other_options parameter to fit with group syntax standards
       790234ed8 BFXDEV-668 - BFXDEV-675 - corrected typo in ucsc.py and bedtools.py
       af1696011 BFXDEV-674 - adjusted walltime for bissnp step
       ba015ef89 BFXDEV-673 - corrected error in job2json.py
       3f5c34e51 BFXDEV-674 - corrected typo in ucsc.py
       47cd3f03a BFXDEV-674 - updated guillimin.ini
       fb5d03cf7 BFXDEV674 - updated base ini
       910066625 unset batch.ini file
       4b8594c04 BFXDEV-674 - briaree ini file
       14ae3c1c8 BFXDEV-674 - mammouth ini file
       58f6b2801 BFXDEV-673 - updated some 'if' statements to avoid syntax errors...
       f036046f1 BFXDEV-674 - updated call of ucsc.bedGraphToBigWig within the pipeline wrapper
       8e187ab86 BFXDEV-674 - added walltime to bismark align step
       a2bcb0c22 BFXDEV-674 - updated bedGraphToBigWig subroutine in ucsc.py, added bedToBigBed to ucsc.py, and updated minor things in bedtools.py
       050654783 BFXDEV-668 - BFXDEV-675 - updated bedtools.py graph subroutine which nows calls ucsc.py to avoid code redundancy of bedGraphToBigWig
       4a6e9cdb1 BFXDEV-668 - BFXDEV-675 - updated ucsc.py to handle cases where bedGraph is in .gz format
       895d8202d BFXDEV-674 - type correction in job2json.py
       a32270549 BFXDEV-674 - cancel the creation of one folder per sample for the json files as on file per sample is enough
       12a085c0d BFXDEV-674 - updated json file path in scheduler.py
       c2f2827ef BFXDEV-674 - updated bismark align & dedup outputs in order to create good dependencies for all_sample_metrics_report and ihec_sample_metrics_report steps
       607ee2227 BFXDEV-674 - changed the json file location to a more simple one and got rid of the resume subroutine since not used anymore
       5b6f982fe BFXDEV-674 - corrected error in scheduler when launching job2json command
       25a4edede BFXDEV-674 - another typo correted in core/scheduler.py and removed the dev references from pipelines/methylseq/methylseq.base.ini
       d5aa74d34 BFXDEV-674 - correcting typo in core/scheduler.py during Json generation
       1dd7bd43b BFXDEV-674 - added some missig report files and updated job2json file command and tool
       c8266f3e9 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       da1de3cdb BFXDEV-674 - updated pipeline wrapper with report files, ihec metrics and updated metrics computing
       c59a3025f BFXDEV-674 - added HPC center name to the ini files
       a128e7048 BFXDEV-674 - updated bedtools when working on MethylSeq pipeline to a better handling of job piping from bamtobed to coverage
       57b022916 MethylSeq pipeline - updated pipeline python wrapper - BFXDEV-674
       807a26bd4 BFXDEV-673 - add a python script in utils to handle the appending of information to the JSON file
       31a092d77 BFXDEV-673 - add the bfx python wrapper to take care of the creation of the JSON file
       4c575e941 BFXDEV-673 - added a sample list object to the Job class to handle JSON file generation
       c9e342e0e ucsc bfx tools - bedgraph_to_bigbwig : removal of the temporary sorted bedgraph after bigWig is created
       f7f011372 MethylSeq pipeline - updated version of the pipeline : steps have been condensed, json file for each sample has been added, more metrics have been added to fit with HiHeq requirements, whole genome datasets are now handle correctly
       be4179ed9 MethylSeq pipeline - updated guillimin ini file for bismark_dedup and bissnp requested resources
       3ee60cd02 MethylSeq pipeline - very minor changes within the base ini file
       34cddb6fc DNA-Seq pipeline - added the verifyBamID step to the pipeline (BFXDEV-619) and json file generation add-on
       a638e4cb7 DNA-Seq pipeline - added verifyBamID settings to the base ini file - BFXDEV-619
       52fb55a25 MethylSeq pipeline - updated common.py so that all the common functions (i.e. sam_to_fastq, timmomatic, merge_trimmomatic_stats, verify_bam_id) now generate a json dump to be append to each sample json file
       1425d42c9 MethylSeq pipeline - updated scheduler.py to handle the json file generation while the pipeline is running, i.e. adds a json section sample json files as soon as a job successfully ends
       5dae0403e MethylSeq pipeline - updated pipeline.py to take care of the creation of the json file for each sample
       8a915230a MethylSeq pipeline - modified sample.py to handle json file creation during pipeline execution
       16c8f934d MethylSeq pipeline - minor update of the verifyBamID python wrapper
       c330a3786 MethylSeq pipeline - reviewed picard add_read_group command
       585c92baa MethylSeq pipeline - reviewed parameters passed to bismark align
       fc34adf7b MethylSeq pipeline - removed jsonator import from bfx, waiting for it to be fully implemented
       3674baf7f MethylSeq pipeline - updated the pipeline with new metrics calculations as well as merged some steps together
       b579365e7 BFXDEV-619 - added verifyBamID .Rmd template file for report
       91071fe00 MethylSeq pipeline - added a python wrapper for all the tools related to methylation profiling
       9a2defcb9 BFXDEV-619 - added verifyBamID in pipelines/common.py
       a3b94d336 MethylSeq pipeline - added GC Bias R tools to metrics
       551335834 BFXDEV-619 - added verifyBamID in bfx/tools.py
       168fb6928 MethylSeq pipeline - Added ucsc.py, with 'bedgraph_to_bigbwig' which is now called by new bedtools.graph function. Also called from methylseq.py
       c3ce44972 MethylSeq pipeline - added bedtools 'coverage' & 'bamtobed' to bfx/bedtools.py, also updated 'graph' & 'intersect'
       946bd7080 BFXDEV-661 - DNAseq - added the use module R_Bioconductor within picard_collect_multiple_metrics
       a5b24ac40 BFXDEV-642 - RNAseq - passed the adjusted p-value as a markdown parameter to the differential expression report
       6e5470ada BFXDEV-644 - RNAseq - added the loading of module_python and the use of 'bash pipefail' for step differential_expression_goseq_report
       38db21934 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d8a33d9af execute permissions updated
       1d8f3a77e Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       dee964a41 MethylSeq pipeline - create specific ini file to handle capture datasets
       3ce2005fd MethylSeq pipeline - create specific ini file to handle capture datasets
       c042cb20b corrected typo in python3.sh install script
       6406b9ce4 added index file as part of the output files in add_read_groups function, in picard.py
       8a940f317 added index file as part of the output files in add_read_groups function
       1c2f36a62 MethylSeq pipeline - updated wiggle track step with splitted forward and reverse strand reads
       040eefb5e MethylSeq pipeline - corrected typo in bfx/bedtools.py
       a1c2b6c81 new CHANGELOG coming along release v2.3.1
       0f801b87e updated Homo_sapiens.GRCh38.sh install script with Ensembl87, dbSNP_149 and dbNSFPv3.4
       ec6489940 updated rnammer_transcriptome resources within rnaseq_denovo_assembly.guillimin.ini
       ae6119c6b updated surpi.sh install script
       ae8ebb522 updated snap.sh install script
       08908b5cf updated seqtk.sh install script
       b6d06b857 updated bowtie2.sh install script
       f040ef648 updated RAPSearch2.sh install script
       ad4b019cc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       41fdf29a4 corrected bowtie.sh install script
       e9772774b Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       79bb77eb5 updated versions of some module install scripts and added some other
       3f90a979e updated versions of some module install scripts
       f6ffda188 MethylSeq pipeline - adjuted resource allocation for bismark align step
       53c3d8528 resolving conflicts and merging
       e5e007433 minor changes to resolve conflicts before merging mMaster with MethylSeq
       eb2157518 Small changes and updates before merging MethylSeq branch to Master to prevent conflicts
       93d9f01fc MethylSeq pipeline - updates done after pur pull request review
       e5eb8fa1e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       adb1cbeb6 updated version to v346 in UCSC module installation script
       faa82153e AmpliconSeq pipeline - updated version of R_Bioconductor to 3.3.3_3.4
       23d2966c7 new genome installation script : Apis_mallifera.sh
       9173ee88c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       396af2525 updated mirbase installation script : removed useless loading of mugqic/R_Bioconductor/3.2.3_3.2
       15d7243b2 Apis_mellifera.sh
       834765f90 updated module installation scripts
       996b10f95 new bash scripts for new modules
       08b52b070 MethylSeq pipeline - few minor corrections before the pull request
       dced67c00 MethylSeq pipeline - adding of the README.md file
       212451efd MethylSeq pipeline - reducing ram to 25G for bissnp step to avoid 'out of memory' issue
       e3d09b59e MethylSeq pipeline - added index creation at then end of picard_add_read_groups step to avoid error in next step (picard_merge_sam_files) when only one readset per sample
       571d0f882 MethylSeq pipeline - corrected ini file for bismark_methyl_call and bismark_coverage2cytosine section where assembly_dir was used instead of bismark_assembly_dir
       f837d91d6 MethylSeq pipeline - added the creation of the bam index after filtering is done, and updated picard.build_sam_index regarding the ini_section parameter
       dafcd65c7 MethylSeq pipeline - reducing ram for bissnp step to avoid 'out of memory' issue
       0f4507d3d MethylSeq pipeline - pipeline with the mapping quality filter step as well as with the picard_calculate_hs_metric step to get ontarget mapping rate
       67e38424c MethylSeq pipeline - updated ini file with newer module versions as well as removing all MUGQIC_INSTALL_HOME_DEV reference and setting them to MUGQIC_INSTALL_HOME
       8f2a7dc12 MethylSeq pipeline - updates after correction of the section name within the base.ini file
       e6dfe0c51 MethylSeq pipeline - corrected section name within the base.ini file
       9c5c314b4 MethylSeq pipeline - corrected typo in the bismark python wrapper
       af1f9f96b MethylSeq pipeline - added tmp_dir parameter to bismark_align step
       ccbf77503 MethylSeq pipeilne - addings of .md files for the first steps
       aaf3ba3fb MethylSeq pipeline - corrected dependency problem within bismark_align step
       57db20f57 MethylSeq pipeline - corrected metrics input and step order
       6248dbff8 MethylSeq pipeline - addings for report generation
       deb45a637 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       02259ac01 MethylSeq - corrected typo within the picard2_collect_multiple_metrics which caused to always restart the step
       83cecace6 MethylSeq - Updates and corrections, pipeline is now fully functionnal on test data
       8bd10dbcb MethylSeq - BisSNP step implemented
       77176b42c MethylSeq - debugging bed_graph step
       b58b12d64 MethyleSeq - correct output files directory for methylation_call step
       d8d7c74f8 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       1422706b9 MethylSeq - debugging methylation_profile step
       ae3a9df72 MethylSeq - updates and corrections, wiggle traks are now generated properly
       4b96e31da MethylSeq - updated pipeline up to 'methylation_profile', still need to work on BisSNP (last step)
       7a2a721e8 MethylSeq - minor bug corrections
       6b0f82010 MethylSeq - First working version of the pipeline : generates bismark alignment and metrics
       89f4a608d MethylSeq - updates of pipeline up to methylation call & also removed/merged some metrics steps
       86a9ae4c3 MethylSeq - refine metrics calling and preparation, up to pUC19 & lambda reads step
       ce3c86d57 MethylSeq - updated pipeline until metrics step
       ff5bee30f Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       78c8fd08d MethylSeq - updates and corrections up to step 10 flagstat
       5a554729f Methylseq - updates and corrections up to step 10 flagstat
       7b0a9a985 MethylSeq - setps 7 & 8 added
       88dfebdd4 update bedtools version in the /bedtools.sh module installation script
       8f9c19840 BEDTools - added intersect function to bedtools.px wrapper
       63b46281f MethylSeq - debugging wrong dependencies du to wrong output file names
       17746b595 MethyleSeq - debugging file names in step 4 & 5 inorder to make step6 picard_merge_sam_files work
       5057bc04b MethyleSeq - defining parameter of step6 picard_merge_sam_files within the .ini file
       b63b02c2e MethylSeq - preparing step 6 picard_merge_sam_file by reviewing/changing file names upstream (steps 4 & 5)
       c5896841c MethyleSeq - bug fix : quoted parameters of add_or_replace_read_groups function in picard(2).py wrappers
       206f1280a MethyleSeq - corrected wrong reference to .ini file section regarding step5 picard_add_read_groups
       d0006b1dc MethylSeq - missing output file is now passed to the picard.add_or_replace_read_groups function
       cae1cdf4d MethylSeq - corrected step5 picard_add_read_groups
       6294b650e MethylSeq - updated step5 of the pipeline : AddOrReplaceReadGroups
       0019f59b7 updated module installation scripts for bowties2 htslib python & samtools
       a8e9c4a7c PICARD - bug correction in python wrappers (picard.py & picard2.py)
       acd02572f GENAP MODULES - added bismark.sh script for Bismark module installation
       9c7010357 MethylSeq - redefined output file for bismark_align step and added AddOrReplaceReadGroups function to picard(2).py wrappers
       222ce028f MethylSeq - python wrappers and scripts updates
       edb1d9bfc MethySeq - corrected typo in __init__.py
       41f69bdf9 MethyleSeq - creation of the (empty) files as a first commit to the branch

  ehenrion <edouard.henrion@mcgill.ca>      19 commits

       4e7c266c7 BFXDEV-673 - updated jsonator.py for a better handling of module names & versions
       39aac2da0 Merged in methylseq (pull request #23)
       11c2e27da chipseq.base.ini : edited module_deeptools to remove reference to mugqic_dev
       bc8a9bed2 README.md updated RAC_ID export line
       a29213cf0 README.md edited : added the export of the $RAC_ID variable which will be used on Cedar for job sumission
       5f6623019 README.md edited : added $MUGQIC_INSTALL_HOME_DEV setting for cedar
       27a43ac96 MethylSeq pipeline - edited guillimin ini file : more walltime for gatk_depth_of_coverage & more cores for bismark_dedup
       0c5f3a0bc MethylSeq pipeline - changed flagstat output name so it is more obviously related to ontarget bam
       3c1f73086 MethylSeq pipeline - inverted order of input file for methylation_call step : "*.readset_sorted.dedup.bam" is now set before "*.sorted.dedup.bam" to avoid unnecessary sorting of bam files...
       f79719efe MethylSeq pipeline - GCbias and mapping_qual_fileter jobs have been added to metrics, their commands were created but not submitted...
       e77d93e21 MethylSeq pipeline - methylseq.guillimin.ini adjusted some ppn values for guillimin
       793836ea8 MethylSeq pipeline - methylseq.base.ini edited bismark align parameters
       4d6b5d04f MethylSeq pipeline - methylseq.base.ini added module_R
       c6a2446f9 methylseq.py : added missing variable
       4f247e670 README.md edited online with Bitbucket
       321816ee7 README.md edited online with Bitbucket
       0bb51e131 README.md edited online with Bitbucket
       a42811f9d bedtools.py edited online with Bitbucket removed unused and unfinished function genomecov...
       e76018c97 bedtools.py edited online with Bitbucket

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      43 commits

       33f242de9 removed depreciated reasignment of module_snpeff in snpsift_annotate step in dnaseq.base.ini
       bad3d894d renamed and duplicated snp_effect section in mammouth.ino to mpileup_snp_effect and haplotype_caller_snp_effect in order to correctly set ppn=1
       54251b5b5 Moved report.kallisto job after kallisto_count_matrix since the it needs the output of kallisto_count_matrix; added dependancies to copy_tx2genes_file so it waits until kallisto_count_matrix is done
       19a3beb09 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       73b4c9786 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eac410701 fixing local merge conflict resolution issue
       e4a5abe00 (hopefully) fixed conflict
       159618830 added a few more os.path.join; added an import to tools
       d4e5e0f25 removing kallisto.sh from module/dev
       5f0a9b875 remove all compiled files
       220ed5372 adding os.path.join where possible
       030d09da2 remove unused sections in mammouth and guillimin ini files
       87289ad67 resolve merge conflicts
       46cb6453d remove commented lines
       2941bfada move rmarkdown calls to appropriate step
       b25aef3a5 fix ini file to correctly set walltime
       76e4c73ba fix an issue with dash in sample names in report
       fa04e0a28 minor improvments to kallisto report
       b55eb983a change path to file in kallisto.Rmd
       6d3e20f86 change job names to fit sections in ini file; replace mention of samples by readsets in md files
       4fd75b29b Changed R code in kallisto.Rmd; add more columns and fix error in the stat table; changed kallisto method text
       1728a1c22 add step to generate transcript count matrix
       341b0fdd3 Change text RNAseq report for differential expression
       1431e2e4a update report files
       9bb9dbfbf added option for single reads, added new parameters in ini file
       f9310f0fe added trimmomatic stats to report
       f81641c55 added first version of RNAseq Light report
       7e53c0f39 added a mkdir command to the merge script
       efb008fca Change path of merged abudance file
       d8dedd1d6 fix job dependancies
       1639f593b fix exploratory function
       9d875ebb8 change path for call to rnaseq_light_kallisto.sh
       4fc4f3bcf adding new step for exploratory analysis
       3831ba609 added a step for merging individual abundance files
       d7208762b add call to module_tools
       d2996e11c fix path to abundanceTranscript2geneLevel function
       5dc4760be adding configuration ini files for mammouth and guillimin
       6c2900a75 make it clear transcrptome must end by idx
       feca8f723 fix path, disable exploratory
       aed4e03e8 new RNAseq_light pipeline with kallisto, dev version
       053dbf0ee added RSEM 1.3.0
       2c86a7f6a updated kallisto
       74e5c0cce update sailfish to 0.9.2

  eloi.mercier@mcgill.ca <eloi.mercier@mcgill.ca>      1 commits

       c82520a4f Approval granted by Rola. Merged in RNAseq_light_dev (pull request #24)

  eloi.mercier@mcgill.ca <emercier@abacus1.ferrier.genome.mcgill.ca>      3 commits

       22084d19e revert sailfish version change
       c8ecd2c7c ajout Salmon 0.8.0
       e2986c581 mise a jour kallisto 0.43.0

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      10 commits

       34124b697 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a9f69a48b add protocole compatibility to rnaseq_light
       cd1caeb38 Make other pipeline supporting several prtocols - BFXDEV-692
       ffdf4f6fd Create 2 protocols with different steps for hicseq - BFXDEV-692
       554f8ad4a Make other pipeline supporting several prtocols - BFXDEV-692
       5ac65fc28 allow several prtocols with different step list - BFXDEV-692
       e634cd2e0 Create 2 protocols with different steps for hicseq - BFXDEV-692
       ece6ab0b1 test multi protocole pipeline
       b10981419  RNAseq & ChIPseq-   Update ini file for the new release of mugqic_tools 2.1.9 - BFXDEV-668  - BFXDEV-675
       e05711475 ChipSeq - finish ihec metrics, preprocess and reformat -  BFXDEV-675

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      21 commits

       6150a8ab7 remove bad output in ihec_metrics_report
       07709bfff Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cc3d29d69 MethylSeq - add MethylSeq.picard_remove_duplicates.md report file
       f8a97db1d fixing conflict in resources/genomes/install_genome.sh
       2a9190be8 change Dedup from bismark to remove duplicat from picard
       9a9b70f61 change ihec methylseq metrics to work on single sample
       720678767 methylseq - generate ihexc report per sample and move methyl_profile lib to tools lib
       900f112bc add specific temp dir to the sort step while generating bigwig
       a9b1ed3e8 increase general recalibration walltime in  Dnaseq to 96h
       764bb1d93 Allow bedtools.graph to support not having the other_options set in the ini
       961e0e198 removing .DS_store file
       fc1defd4a ChIPseq - address reviewer coments - BFXDEV-675
       f8bf3042a Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       cb48b4514 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       1162758ef RNAseq - add concat job  - BFXDEV-668
       c8591209b  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ee2eab022  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       188462d56  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ba1f73ce3  RNAseq & ChIPseq-   Update (chip) and debug (Rna) IHEC metrics steps - BFXDEV-668  - BFXDEV-675
       9e71317e2  RNAseq -  implement IHEC RNA metrics step - BFXDEV-668
       7cfbb6d83  RNAseq - BFX - implement mugqic_tools module for the IHEC RNA metrics generation script - BFXDEV-668

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       f19e749c5 Merged in IHEC_metrics (pull request #21)

  pascale.marquis2@mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      2 commits

       852ef1c00 /#cluster_queue
       eb94727e8 python/2.7.12

  Pascale Marquis <pmarquis@lg-1r17-n03.guillimin.clumeq.ca>      1 commits

       1052ba180 update tumor_pair.base.ini

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       df577fb63 Minor bug fixes and addition of base exome ini

  Rola Dali <rdali@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       46371a980 starting the hicup_align step

  Rola Dali <rola.dali@mail.mcgill.ca>      88 commits

       9be347485 Merged in IHEC_metrics (pull request #27)
       130d0b56c changes to mammouth.ini to set all ppn=1; changed module spider back to module show since it is incompatible with abacus
       71cbe96aa Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       f93c39e5a module show changed to module spider in config.py to accelerate module checking
       a64c34fa3 fixing homer dependencies
       d27429396 editing homer tag directory output back to folder
       78dc2fc7e Merged in IHEC_metrics (pull request #26)
       6b9787115 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       ac73327c8 homer edits to generalize methods
       649af9726 resolving merge conflicts
       fd5cb3876 dependencies in rnaseq metrics
       23921c9c1 rnaseq metrics
       62da3e6d4 fix ihec_metric job names
       fb831740c adding chip_type to ini to avoid pipeline crash
       b7e20902f run_spp fuctional
       85dd69ab6 added run_spp to calculate dnsc and rsc. TO TEST on mp2
       1459b137c adding TMPDIR to bigWig creation
       a058aa11e sample name change in ihec_metrics
       31c83eb58 fixing TMPDir issues on guillimin-to test
       80e315fb4 pipeline to run without input & chnage in ihec_metrics format
       14fde80f5 moving homer code to homer.py moduke
       d33094e27 Merged in hicseq (pull request #25)
       fe6d2618f resolving juicer conflict
       3cd068585 resolve merge conflict
       eff83a22b resolving merge conflict
       574b511e9 adding C3G logo
       400b590a3 commit for pull request
       18c1112fc commit for pull request
       27b87a360 adding bedops to annotate capture file
       198416bb1 adding runchicago to capture hic
       742376e2f Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c26ea1df5 Resolved merge conflict
       cc289e6c1 commit before Mathieu's pull
       64cd09085 adding capture seq
       46b4f6bba merging capture hicseq with hicseq
       4b965527f chicseq
       b19d167c4 Merged in hicseq (pull request #22)
       a25f08bd6 added new chromosome contigs based on ensemble and NCBI genomes BFXDEV-670
       11b2bb953 changing perl to perl/env  BFXDEV-670
       a790fa186 added genome.py for genome related methods BFXDEV-670
       6a6312457 changes for pull request BFXDEV-670
       f9a255254 changes for pull request BFXDEV-670
       0d65329dd changes for pull request BFXDEV-670
       aa2e874e0 edits for merge request:wrappers. BFXDEV-670
       de38b5a36 BFXDEV-670 hicseq merge edits
       d69045fe7 creat_hic_file edits BFXDEV-670
       1649642e9 adding juicer.sh installation script BFXDEV-670
       665ea381e Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c4adfa90c testing create_hic_file BFXDEV-670
       4a67afe51 modified hicseq.briaree.ini and batch BFXDEV-670
       a37fe8a8c deleted files moved to mugqic_tools BFXDEV-670
       bcf5b1a1d moved genome_digest before ncrna method which is failing in install_genome.sh BFXDEV-670
       861cc7998 added module install files and genome digest BFXDEV-670
       15230602e changes for pull request
       ad4a14128 commit before merge changes
       98d794719 added samtools_bam_sort BFXDEV-670
       6fd3b9381 edited ini files BFXDEV-670
       18e077477 samtools_bam_sort and multiqc_report in testing
       3c9468bc7 Fixed hicup/TAD/HiCPlotter restrart bugs BFXDEV-670
       6a58aa4ba split chr interaction matrices from plotting BFXDEV-670
       dca827d2e HiC v1.0 is ready for testing on guillimin and abacus BFXDEV-670
       2b40a8c73 bam merging in testing BFXDEV-670
       65763b14c samtools_merge_bams in testing BFXDEV-670
       d2ff7116f added petagDistDistribution to HomerQc plots BFXDEV-670
       e8edda6ac added output_dir property to reduce code redundancy
       ded82f6f7 fixed path in identify_compartments
       1b31f4003 added identify_peaks BFXDEV-670
       d212eba1d compartments and TAD identification now working BFXDEV-670
       384d836eb chr and genome interaction system now working BFXDEV-670
       6714aabc6 interaction matrix plots testing BFXDEV-670
       edadc7445 homer archiving and Qc plotting now working
       8c94b93f8 first 6 steps working BFXDEV-670
       5c3eb2876 resolving git merge issues
       e1d46d6fc Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       399cc72e7 syntax changes
       c545fc629 added fastq_readName_Edit() BFXDEV-670
       2b7ee9eae mammouth.ini cpu set to 1
       4afc530b1 genome digests on mammouth
       e4e52d38b update genome digest files
       92773981b homer_tag_directory archiving in testing BFXDEV-670
       acb862f1b basic formatting changes
       53257d6b5 make_tag_directory in testing
       148eb65fc works to produce hicup bam and library Qc. ITS ALIVE :)
       958242d30 ITS ALIVE
       c446402c5 need to expand  variables
       099e5ef2a pipeline now accepts enzyme
       ffbdf83a4 hicup_align producing script; to test tomo
       884af2600 initialising hicseq analysis pipeline

  Xiaojian SHAO <xshao@lg-1r14-n04.guillimin.clumeq.ca>      2 commits

       6ddab51be add ppn to wiggle_tracks step
       06bd6a7fd add ppn to wiggle_tracks step

  Xiaojian SHAO <xshao@lg-1r17-n03.guillimin.clumeq.ca>      3 commits

       6fd6ec5f0 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cd938a3ac add walltime to bismark aligner
       8a5b82631 ppn_Changes_in_Guillimin.ini

  Xiaojian SHAO <xshao@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       db2901fa2 methylseq: edit on ppn setting. -xiaojian

2.3.0        Mon Feb 27 13:40:01 2017 -0500        82 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      13 commits

       7717c35eb Pre-release - adding the installation scripts for all the new software used by tumor_pair pipeline
       8eada6709 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8e17bca94 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       59cf94076 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       026ccb111 BFXDEV-602 - PacBio Assembly pipeline now contains a new optional step (step 10) to circularize the successful/clean assembly contigs with circlator - also comes with some minor unrelated updates
       06e095f69 updates brought to many module install scripts, and adding of new modules
       9d565ba2c RNASeq - corrected a typo inserted after correting bedtools bug...
       982a0145a RNASeq - corrected a bug bedtools.graph function : samtools_options now handles reverse strand specific parameters, avoiding an empty begGraph for reverse strand
       72d1a3f23 updating python & python libraries installation bash scripts
       6242c97c2 BUG correction within mpileup function : parameters were shift after introduction of 'ini_section' parameter for tumor_pair prpose
       1ca38fa77 DNASeq - bug correction after merging tumor_pair branch to master
       e6b482bf4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d8876d70 FastQC - updated bash install script to match our standards

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       0e875a566 README.md edited online with Bitbucket

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       b5b1d51e9 tumor_pair - add comments to dict2beds function - BFXDEV-521
       e0590a428 tumor_pair - add feature to dict2beds function - BFXDEV-521
       7c0efb4dc Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       bd7e169a3 tumor_pair - modify ini file to take into account the new version of scalpel (CVMFS) - BFXDEV-477
       8e735ac6c tumor_pair - add space charter before scalpel option - BFXDEV-478
       33a9347ae pull origin Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       5ed57839e tumor_pair - add the two pass option to scalpel - BFXDEV-478
       bf963a388 tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing - BFXDEV-476
       18110f56d tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      39 commits

       d6be75aa0 remove .gitignore
       3821384f2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8e31719ae rnaseq - increase abacus ressource for STAR
       6f5a17b7b remove conflict with tumor_pair changes
       07cab4fcd remove bfx/samtools.py conflict with tumor_pair changes
       abecc6c7c remove bfx/picard.py conflict with tumor_pair changes
       4368b7caa merge and remove conflicts
       2d5b7c568 merge and remove conflicts
       dfcbf627f tumor_pair- update germline_loh ensemble
       376ff9263 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       07b9e8b56 add encode  sample+readaset+pipelien to log system
       f4c9b602e Tumor pair - add comment on the bed system
       ba4770244 Tumor pair - change vardict input output to correct deepndency
       0e20097f5 Tumor pair - change varscan input output to correct deepndency
       4ab037ec3 sequence_dictionary.py - split_by_size - correct bug create a first job with everybody when the nb of split was too high
       b5a7b3371 tumor_pair - add output bam index file when only 1 realignment is produced
       a5b4491bd tumor_pair - remove file name bug in indel realignement mv
       4429b197e tumor_pair - rebuild gatk indel recalibration input/output scheme in tumor
       b9030e544 bfx/gatk - indel realigner add target list file as input for dependency
       440001c6d tumor_pair - remove duplicated lines in guillimin.ini file (from base.ini)
       6c1e60602 tumor_pair - correct symlink creation
       b15f11845 tumor_pair - remove few issues (paths, dependency, input/output)
       329a56c7b tumor_pair - put WGS as default (instead of WES)
       271ba292c remove issue with uncorrect interval list when no bed is attached to the project
       517e8c3a9 adjusting tumor_pair
       c049fe240 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       c2c0986df tumor_pair - debug bfx/bcbio_variation.py
       5c6034d29 tumor_pair - remove conflicts
       d2c00151d rtumor_pair - emove conflicts
       7c4bd8a1f remove conflict
       c67b83756 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       9cbd2220a tumor_pair - correct typo in scapel download adress - BFXDEV-477
       760faaa09 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       8d2a24339 tumor_pair - add scalpel install script - BFXDEV-477
       e6beef4c0 tumor-pair - mv tools.py from highcov branch to tumor_pair - BFXDEV-475
       1bc346fa4 tumor-pair - adding bedfile spliting process - BFXDEV-476
       5acad27d1 tumor_pair - extract tumor_pair code for the high coverage branch - BFXDEV-475
       6c73c2b3c tumor-pair - adding bedfile spliting process & start implementation in tumor_pair.py - BFXDEV-476
       905c1acdd tumor_pair - extract tumor_pair code for the high coverage branch

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      5 commits

       192f4d355 Bug squashes and speed improvements to ensemble processes
       40df5d932 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       93a2142f6 Modified to remove analyses of alt contigs + split Varscan2 by chromosome
       8d3f03fc9 Beta version: module files created and code tested on wes and wgs on abacus and guillimin
       755b0c063 update to Mutect2, added Vardict, re-added samtools and added ensemble approach

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      4 commits

       e94352f07 Fixes to samtools modules calling for pairs
       9729909b7 Numerous speed improvements and addition of fast variant calling
       3f1aceb16 Additional fixes to resource allocation issues
       0f04f2aff Dealt with comments, added paired indel realignment, varscan2, and seperate somatic and germline call files

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      11 commits

       67b5c5326 Debugging samtools pair bcftools
       3a4996e3f Fixes to samtools due to dnaseq changes
       a920ef82f Update README.md
       2c7ce4083 Bug squashes and preformance improvements to ensemble process
       3349a8186 Bug squashing and speed improvements to ensemble process
       15a7f87fc Dependency and other bug fixes. Guillimin ini and install scripts updated
       a21682a63 tumor_pair - fixes to tumor_pair README.md
       196643408 tumor_pair - fixes to README.md and resolving of conflicts
       c53612c47 tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       dd6df5561 tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       c708f88e5 Updates/fixes from guillimin test

2.2.1        Mon Dec 19 10:57:33 2016 -0500        212 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       cd7c8e746 Added proper error message when running script with too old Python version.

  Edouard Henrion <edouard.henrion@mcgill.ca>      138 commits

       68f036667 GenAP Pipelines - updated tmp_dir variable within all the .base.ini files for a better use of the memory in the compute nodes on abacus
       c364c3ec1 modules - updated version of VSEARCH (2.3.4) within the module installation script
       ed99e0e8b GENAP PIPELINES - updated all the .base.ini files to set tmp_dir to /lb/scratch/ehenrion instead of /lb/scratch/ on abacus
       d0faf2edc RNASeq - added __init__.py
       a8d110117 PICARD - bug correction in python wrappers (picard.py & picard2.py)
       7d2e9b2f1 install_genome.sh - corrected MUGQIC_INSTALL_HOME_DEV link
       b2ec14eeb install_genome.sh - set abacus2 variable to manage mammouth execution properly
       c3549edf2 install_genome.sh - updated for a better handling of mammouth cases
       472693245 Homo_sapiens.hg19.sh install script - updated URL to retrieve dbSNP vcf
       b0dfe8265 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       465b913af RNASeq denovo Assembly - updated mugqic_tools version to 2.1.6
       55e645641 MUGQIC MODULES - updated bedtools version to 2.26.0 within the bash install script
       797fe54f1 MUGQIC MODULES - updated samtools version to 1.3.1 within the bash install script
       f84674658 MUGQIC MODULES - updated bowtie2 version to 2.2.9 within the bash install script
       7ef8596bb MUGQIC MODULES - updated bowtie version to 1.1.2 within the bash install script
       eca9f4d3d MUGQIC MODULES - updated bismark version to 0.16.3 within the bash install script
       b92f63f27 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2b0f2e01c RNASeq_denovo_assembly - updated guillimin specific .ini file for trinity & insilico_read_normalization steps
       5c1160937 DNASeq - updated guillimin specific .ini file for snp_effect step
       639668a1a DNASeq - updated .ini file for snp_effect step
       14f4c0f5a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b3a054e53 PICARD 2.0.1 - updated picard2.py wrapper with proper syntax
       7a20e2f85 briaree pipeline .ini files
       b4b5e8c86 RNASeq_denovo_assembly - updated mammouth .ini file
       b3ffcc7f1 python.sh - updated python version (2.7.12) within the bash installation script & added all the python library installation steps : no need for python_lib.sh anymore
       731b02993 Homo_sapiens.GRCh38.sh - updated Ensembl version (85) & dbNSFP version (3.2c)
       ba150756c install_genome.sh - update cmd_or_job function to automatically handle mammouth environment cases
       5d07d3100 exonerate.sh - update version to 2.4.0
       c1ac0e432 RNASeq_denovo_assembly - back to hmmer version 2.3.2 (from mugqic modules)
       1f6475796 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       16c346402 exonerate - within the bash install script : corrected broken link to the archive
       875e3eb0d Genome reference .ini files updated for each reference installed on CVMFS
       8c4c2ed03 module install script - corrected emboss.sh : unresoved conflicts are now resolved
       21a06d7f5 module install script - corrected emboss.sh : unresoved conflicts are now resolved
       5f46cf695 Python packages - added TEToolkit package to the installation scripts
       f4a2adb51 FASTQC - updated installation script with newer version and integration of /lb/project/mugqic/analyste_dev as a possible installation path
       a396762d4 RNASeq_deNovo_assembly - updated ini file with integration of picard v2.0.1 and working version of trinity (2.0.4)
       916b8d6a2 PICARD 2.0.1 - added picard2.py to the bfx tools
       ef03d6bfe commiting minor updates, prior the new release
       aa0f55fb2 RNASeq_denovo - briaree specific config for insilico_read_normalization_all
       29fc52abe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd77e98e7 corrected version of mugqic_tools in rnaseq.base.ini
       daad645b3 corrected version of mugqic_tools in rnaseq.base.ini
       e2d140ece RNASeq de-novo assembly - bug correction when calling trinity with new default parameters
       87d68339a RNASeq de-novo assembly - updated ini file for briaree
       1605b8fe8 modules - new STAR version handled in the installation script
       8f88daae6 RNASeq de-novo assembly - updated trinity call with newer version of trinity and new default parameter
       73ca888ac RNASeq - update differenial expression to handle local fit option if needed (i.e. if parametric dispersion fit fails)
       ec71e94a3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d0d765e2 new reference genome installation scripts (new species)
       0029c2a3a metaMarkerSeq - guillimin .ini file adjustments regarding qiime
       e4c1a5b97 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b5cc05d5a added butter installation script
       8ec7a74db Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       16f2ae774 qiime db installation scripts and some new module scipts
       8fe19f83c chimera_unite.sh added
       ad09a1ee8 new module installation scripts
       a4f4b053f some updated & new module installation scripts
       8fd6a71dd updated genome installation scripts
       083ff4e51 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       846d93199 metaMarkerSeq - removed bug with amplicon_type within the .ini file
       2ba2d4c10 metaMarkerSeq - updated python & bash script as well as .ini with the newer silva db version
       4504c46b9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e68ffa31a RNASeq - fixed bug in gq_seq_utils_exploratory_analysis_rnaseq
       903ef8db6 mammouth ini files added for dnaseq_high_coverage
       def0206ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c822e0d73 mammouth ini files updated for pacbio & dnaseq
       887e51897 adding pipeline ini files specific to brairee environment
       e5e3585cd DNASeq - review some parameters within ini file during the pre-release pipeline testing
       cece203d7 RNASeq - allow sample names to have '.' within their name without crashing during gq_seq_utils_exploratory_analysis_rnaseq
       46d069f2b fixed typo
       54f260b40 RNASeq - fix mammouth setting for htseq_count step in rnaseq.mammouth.ini
       ada673813 RNASeq - restoring bfx/report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd
       35c5cad3a ChIPSeq - metaMarkerSeq - adding briaree specific ini file
       f0112e131 RNASeq - correct default genome in the ini file
       4fb234d88 pushing ampliconseq.guillimin.ini after testing on guillimin
       90ff7dca8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4ec8c8298 pushing ampliconseq.mammouth.ini after testing on mammouth
       4d8e58322 Corrected typo in dnaseq.base.ini
       6e9c3b2a4 correct mammouth master branch divergence
       8ab80b30e updating module and genome scripts from mammouth
       7e2849060 remove conflict
       2c439fe11 conflicts resolution
       4db2a9e4b some residual diffs to commit
       de4c04364 redo the call of bcftools module within samtools.py
       b69fa6e8c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ad24f0cd3 resolving conflicts after merging...
       6330a7f45 commits before merging
       5e4eac401 DNASeq/RNASeq - code correction after testing on briaree
       580589a9e metaMarkerSeq - Merging branch to master git add resources/genomes/Bos_taurus.UMD3.1.sh resources/modules/perl.sh resources/modules/star.sh resources/modules/vcftools.sh
       0fae67355 metaMarkerSeq - reformat the report by removing some redundancies in the paragraphs and ensuring correct links to the tables/figures - BFXTD-26
       8873109b7 BFXTD-26 - corrected the template and links used for report generation - MetaMarkerSeq
       a0c0188ae Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e18bb247c correct a bug in wiggle step which causes *forward.bedGraph and *reverse.bedGraph files to be identical
       850b148cd metaMarkerSeq - commits after reviewing pull request comments sent by Marc Michaud - BFXDEV-453
       23901609d threshold instead of treshold
       51848d9fb metaMarkerSeq - debug Krona input generation and updated md files for report - BFXDEV-453
       475035b1e meatMarkerSeq - get rid of remaining  variables - BFXDEV-453
       0fd827d51 AmpliconSeq - debugging report template files
       076e6f71b bug fixes and file names correction
       044141557 ampliconseq - tuning of the ini file - BFXDEV-434 BFXDEV-448
       feac3c052 ampliconseq - new/updated module install scripts related to the pipeline needs - BFXDEV-434 BFXDEV-448
       3865cfaf8 ampliconseq - code review & debug prior to relase - BFXDEV-434 BFXDEV-448
       664910aaa AmpliconSeq - modified some parameters in .ini file & other minor syntax changes within the wrapper - BFXDEV-434
       b22e6e2ff ampliconseq - conflict resolution - BFXDEV-434
       28788ceb8 ampliconseq - code review of ampliconseq.py - BFXDEV-434
       d07b3d86b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ba4bca971 BFXDEV-523 - error correction & version update in the gemini installation script
       4e7e0a604 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f22df9524 BFXDEV-516 - Macaque genome installation bash scripts + other minor updates
       90d34cd63 BFX-4532 - added prepend-path entry for dnaseq_high_coverage in mugqic_pipelines.sh
       7f41180aa RNA-Seq - remove useless parameter from rnaseq.base.ini
       7442c01e4 RNA-Seq - increase default walltime for star_align & star_index steps
       c9d36941d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a9e67a58 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       839b257d4 genomes - updated versions of genome installation scripts, essentially fixing STAR indexes installation - BFXDEV-494
       a3da825c3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       192a9717e STAR - updated STAR version (to 2.5.1b) in the installation script (star.sh) as well as in the RNA-Seq pipeline .ini files - BFXDEV-514
       646ad6ac6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       394913df1 Version bump to 2.1.0-beta
       0dac262fe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7a78685a1 BFXDEV-490 - recommit after conflict resolution
       8fc43e3cb BFXDEV-490 - updating chipseq.base.ini, dnaseq.base.ini, pacbio_assembly.base.ini, & rnaseq.base.ini
       cfd56b995 ampliconseq - module - adding vsearch.sh in the module script folder - BFXDEV-434
       fb1f4d2a8 ampliconseq - base ini file updated - BFXDEV-434
       ed7611f63 ampliconseq - updated ampliconseq ini with new module versions - BFXDEV-434
       8c047cccc ampliconseq - updated module versions - BFXDEV-434
       565ab9989 Merge branch 'ampliconseq' of bitbucket.org:mugqic/mugqic_pipelines into ampliconseq
       d2ef4f43e AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEV-434
       b15af326e AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEX-434
       ace14e4da AmpliconSeq - revised and standardized code of ampliconseq.py - BFXDEV-434
       76c924e68 AmpliconSeq - updated ampliconseq.base.ini - BFXDEV-434
       83d3b0e43 AmpliconSeq - new VSearch library added to bfx - BFXDEV-434
       96ce6075a AmpliconSeq - updated tools.py with AmpliconSeq functions in bfx - BFXDEV-434
       b886764d7 AmpliconSeq - new Qiime library added to bfx - BFXDEV-434
       bb008af5c AmpliconSeq - new Krona library added to bfx - BFXDEV-434
       d9016c867 update of resources/modules/mugqic_tools.sh - BFXDEV-490
       061158490 ampliconseq - updated syntax and some corrections

  Edouard Henrion <henrione@briaree2.rqchp.qc.ca>      2 commits

       91a47a540 RNASeq - adding .ini file for briaree
       ce7f2aedc small adjustments after cloning/before testing on briaree

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      1 commits

       5806c4a6d ampliconseq - bump README.md to the latest version

  ehenrion <edouard.henrion@mcgill.ca>      3 commits

       62133d094 README.md edited online with Bitbucket
       d988c477b README.md edited online with Bitbucket
       5edde6884 README.md edited online with Bitbucket

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      1 commits

       0c8624061 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       722734200 mini and spades modules
       be71d1c87 nxtrim and quest modules updates
       5ec026e61 dev install scripts for MAKER

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       2f9db74ed HighCoverage - add missing README.md file
       21f73de26 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5e66e13a1 Removing unmaintained pipelines (PUURE & rRNATAGGER) from master
       435f7d53f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       33667103c RNAseq synchromize mammouth ini with the base ini (missing tuxedo_hard_clip) - BFXDEV-515

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a1574375d AmpliconSeq - remove ini recursive interpolation
       84a3e9b9b remove conflict
       0e5bd2a71 High coverage - debug vt and gemini step - BFXDEV-542 BFXDEV-541
       28ce1e1e2 ampliconseq - remove dependencies issues and remove try instances - BFXDEV-531
       342169d8c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ed67908b4 DNASEQ - bwa always set LB in RG tags: when library barecode not given used sample name instead - BFXDEV-354
       e28eff7e4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       73bce2e79 ressource - python lib: add futures for multithearding
       a7c0b13ee update install module general script

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      5 commits

       8a0bb7683 dnaseq.py  - gatk DoC job name wasn't following the rule: <stepName>.<sampleName>.<OtherInfo> . corrected
       5056c01da Merged in irp_modifications (pull request #17)
       4991ce6dd change mail for bug & issue in  README.md edited online with Bitbucket
       ebe9e9896 README.md edited online with Bitbucket
       0cdc302a0 README.md edited online with Bitbucket

  mmichaud <marc.michaud@mail.mcgill.ca>      20 commits

       ed1aece18 Increase max jobs per step to 24 (cover most case without overloading the cluster) (cherry picked from commit 4b06703)
       583664980 BFXDEV-559 Error when a lane contains only double-index librairies on a single-index run (cherry picked from commit 19bd532)
       b2f9c7ee1 Increase even more STAR walltime to allow fetching the reference index from cvmfs.
       a8b9ffe3c Don't mark the blast job as failed if there is no blast results. If grep can't find any results, it exits with a code 1. After, when we test the $PIPESTATUS, we intercept that error code and mark the job as failed. By adding a " && true", we preserve exit codes but we reset the $PIPESTATUS.
       5c5d9336d Increase STAR walltime to allow fetching the reference index from cvmfs.
       61b608ba5 BFXDEV-538 Add perl for bed2interval_list
       4d7247181 BFXDEV-538 Update modules version
       e39b20274 BFX-4696 Use a newer version of python having all the needed packages.
       d141b2fd7 BFXDEV-544 VerifyBamId: Don't filter the variant file with the genomic target file.
       0237e8ac2 Merge branch 'master' into irp_modifications
       5849c790b Illumina Run Processing: Accept "." in asembly name in the genomic database. (fix regex)
       e020c6455 Illumina Run Processing: Accept "." in asembly name in the genomic database.
       8806df0ef Illumina Run Processing: Use cvmfs version of verify_bam_id.
       6735cae38 BFXDEV-512 Illumina Run Processing: Use "Parc" silva database.
       8aa74f21e Illumina Run Processing: Add "nanuq_environment" configuration variable to be able to ease testing on the QC environment.
       3d2bec95a BFXDEV-533 Illumina Run Processing: The exclude bam option doesn't exclude the sorted bam.
       776850b78 BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       7d1694dce BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       edfb2e4be BFXDEV-527 Illumina Run Processing: Merge jobs of a task when their count exceed a configurable threshold
       cd6cfde65 Tweak cluster parameters : Use parallelGCThread for all available cores. BvaTools run faster with only 4 threads.

  ptranvan <patrick.tranvan@mail.mcgill.ca>      24 commits

       46ee612bb Eulidean distance option for PCOA
       8ca5e0ec9 README link correction
       18b9818ab Module perl correction
       d63b31394 Module perl correction
       b4fc989f5 Format edition
       39c2882b9 . file deletion
       c2cd2c3e0 BFXDEV-434 README update
       e0d27cc4a BFXDEV-434 Create inis for other clusters
       f4475b93c BFXDEV-434 Deploy db files in $MUGQIC_INSTALL_HOME
       dc7f9fa1c - Other pipelines have been added to the branch. - Create inis for other clusters, leverage overlaading of parameters. - put db parameters in DEFAULT section,. - QIIME section has been exploded in the config file. - Tutorial in README file.
       27f3fa58a README and configuration file correction
       5a9fdb726 Configuration file modification.
       3419fac41 Adding alternative normalization method
       a7ade41ef Dependance correction
       d4e383d38 add a step (close ref)
       6b63f3f51 CLuster algorithm change
       a5d887614 tutorial.txt modification
       72ef558d4 Minor modifications (report and option)
       9cf47dcb6 Amplicon-Seq pipeline + report upload for test.
       c9c6125c8 Change name: metagenomic pipeline to Amplicon-Seq pipeline Adding new features (alpha, beta plots) Final step implemented
       d7db13767 Adding metagenomics pipeline (1st part)
       075ec74db Adding metagenomics (amplicon) pipeline.
       990ba4b41 Merge branch 'metagenomics' of https://bitbucket.org/mugqic/mugqic_pipelines into metagenomics
       f779d06eb test

2.2.0        Mon Feb 8 12:04:44 2016 -0500        405 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       dc03d86e0 Added link to the GenAP project in front page.

  Edouard Henrion <edouard.henrion@mcgill.ca>      63 commits

       477b33246 Version bump to 2.2.0
       caa197dde Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d87ea7e9 BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       371dd66a3 add report section in DNA-Seq High Coverage pipeline ini file
       514d5b613 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       be6d52a1f updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       61ac85c27 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       042bb0381 modules - updated vcftools VERSION in vcftools.sh -  BFXDEV-490
       f71fc0149 minor update in module file created by perl.sh
       d3c2bbfe0 add report section in DNA-Seq High Coverage pipeline ini file
       cd5d58c05 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       7d0e6f005 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       7bf0a3925 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       349ae58e6 modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       ab4b16db4 module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       9e04f6a86 modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       e93f812af minor update in module file created by perl.sh
       fe2a9cbaf rnaseq & rnaseqdn ini files updated after testings on guillimin - BFXDEV-490
       39e91f13c add report section in DNA-Seq High Coverage pipeline ini file
       e3bce2937 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       ecc918b5e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a4f79773e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       442830db3 updated ini files - BFXDEV-490
       175b162e4 BFXDEV-490 - resolving conflict on guillimin
       36b6f5db4 BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       d1b51cbd3 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       a3e2ab82d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ec80af334 modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       b6e687aec minor update in module file created by perl.sh
       97bd1559d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       682edbee2 mugqic_tools.sh VERSION=2.1.5, again...
       fd3f0e3a5 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       424768a78 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4b4c08ac5 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       f6415e121 modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       045c237f7 module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       d2072f58c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c22a59618 modules - corrected typo in perl.sh when creating module file - BFXDEV-490
       f89dbb336 modules/genomes - updated module version calls as well as database releases - BFXDEV-490
       3970c215c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       17531c78e module - debugged scalpel.sh script - BFXDEV-490
       a05c71b14 module - updated scalpel version - BFXDEV-490
       1c9d8682d module - debug  in snap.sh  - BFXDEV-490
       eee7e5674 module - updated version of Python (to 2.7.11) in gemini.sh - BFXDEV-490
       8018094b1 module - removed module load calls in shortstack.sh - BFXDEV-490
       8c699539d module - updated star.sh - BFXDEV
       071a17b25 module - updated shortstack.sh snap.sh - BFXDEV
       3fa5330c9 module - corrected typo in sphinx.sh - BFXDEV-490
       976b8565b module - vcftools updated version to 0.1.14 - BFXDEV-490
       0d5709bc2 module - debugged gemini.sh - BFXDEV-490
       add9770bc module - samtools updated to 1.3 - BFXDEV-490
       15428f890 module - debug the name of the archive - BFXDEV-490
       7ab3bbfbf module - UCSC version set to v326 instead of 'latest' - BFXDEV-490
       16571fccb update version of bwa module - BFXDEV-490
       bf6422a8b some more updated modules - BFXDEV-490
       0cf2a4d5d module & genome updates - BFXDEV-490
       39c17582a module updates for to the release - BFXDEV-490
       bb09f04c1 resources/modules/dev/epacts.sh has been removed (really)
       0815aec62 resources/modules/dev/epacts.sh has been removed (really)
       27058bcf1 BFXDEV-490 - update of resources/modules/mugqic_tools.sh
       0021aaa25 Merge branch 'gatk_variants' of bitbucket.org:mugqic/mugqic_pipelines into gatk_variants
       a41229e67 gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted - BFXDEV-481
       b3ae471e4 gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted

  Edouard Henrion <henrione@ip03.m>      2 commits

       839c46445 updated chipseq.base.ini after mugqic_tools update - BFXDEV-501
       476279f77 updated chipseq.base.ini after mugqic_tools update - BFXDEV-501

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       565ea8b27 BFXDEV-491, missing argument to re.sub for the SINGLE_END STAR commends case
       bc9c0dfeb Reverte back to older scheme. Was assuming bacteria only.

  Francois Lefebvre <lefebvrf@gmail.com>      37 commits

       ae9667fb8 Mammouth rnaseq ini required cluster_cpu=-l nodes=1:ppn=1 for new steps related to rRNA estimation
       e21b18cf4 -S flag was missing from the last samtools view command in the hard clip command
       f5734af02 Install scripts
       6eabc352d minor mod to R_Bioconductor (removed tabs in HERE DOCS)
       432b3abcb sleuth R package
       93c3b1da1 Removed old R installation scripts
       7758ea3c7 Added slash to URL to be able to retrieve latest Bioc version
       b80611a23 sspace-longread dev install script
       b626deb2e Kallisto install script (abacus only)
       0f9449a9b popoolation install scripts modifications
       dd9e03134 more notes on pacts
       1e9442f60 Added notes to EPACTS install script
       ee112de27 Fastx and EPACTS install scripts
       a1b52899b install scripts for ray, proved, sortmerna
       a08674428 Added celera and loRDEC install scripts
       c56d2979f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01d05df10 Merge branch 'knitr_test'
       62efc5fec Added PBJelly install script + networkx python module
       169100299 Fixed typos
       394fca5e9 Replaced library calls in rmarkdown.py with direct ::render call. Added missing section header to exploratory report.
       7d614146a no message
       cf13cafba Modified exploratory analysis report text
       a8fd54f45 import rmarkdown in wrapper
       abd9f895d damn dots
       336b15691 Testing out markdown::render paradigm
       95df30155 row.names=F for exploratory summary table
       0f11b5573 Working directory madness fixed
       4de25fc02 Trying out rmarkdown::render() instead
       22d5cc6f9 Trying html output
       2c3caa858 knitr would set the wd to input document location...
       c62489776 dir.create problems
       960c5e515 quotes missing
       01e5ea509 no message
       94b1e9b25 EOF is simpler
       5742f004a Had forgotten the module calls
       65e698e53 knit for exploratory + other changes
       39f86fe8f Created Platanus install script

  Gary Leveque <gary.leveque@gmail.com>      2 commits

       248c455c0 pacbio_assembly.base.ini edited online with Bitbucket
       93cfad19d pacbio_assembly.base.ini edited online with Bitbucket --changed smrtanalysis_version = 2.3.0.140936.p2 to: smrtanalysis_version = 2.3.0.140936.p4

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      1 commits

       2d582b48d changes necessary for bacterial RNAseq using STAR --see BFXDEV-449

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      52 commits

       f7e38cbad BFXDEV-405 correct single library issue in insilico read normalization
       a7457994b Merged in rnaseq_denovo_assembly_new_features (pull request #10)
       b1ecdcc0c BFXDEV-397 resolved rebase conflict resources/modules/verifyBamID.sh
       23a6e08e0 BFXDEV-397 PRJBFX-1187 added genome and software installs for verifyBamID, resolved issues from pull request #10
       296aab1a4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       91707a583 BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features removed dev modules
       f80579a07 BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features
       8ebbedb99 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d492790ad BFXDEV-397 added spaces to table captions to goseq and dge filtered outputs
       c03656d69 BFXDEV-397 added DGE and GOseq using filtereg contigs (based on trinotate annotations)
       c22ce84cf README.md edited online with Bitbucket
       706303f1a BFXDEV-439 force installation or rmarkdown and knitr
       ab8526cf3 BFXDEV-439 resolved rebase conflicts rnaseq_denovo_assembly_new_features
       716c646b5 BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       028430fc0 BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       1905e38b0 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, commit new libraries
       13cefa952 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report, regenerated rnaseq_de_novo pipeline README markdown page
       8d5737a9c BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report
       b40353bcd BFXDEV-397 correct dependency problem in exploratory analysis using filtered isoforms
       59fdc4e22 BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       107fe2a9a BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       b85f60b45 BFXDEV-396 added exploratory analysis using filtered transcripts
       17e2143d3 BFXDEV-423 BFXDEV-399 BFXDEV-397 corrected blastx, abundance estimates to matrix and generated tabular and fasta filtered assembly using python SeqIO
       de22cfa8c BFXDEV-432 define module picard rnaseq_denovo_assembly base ini
       38207383f Merge branch 'rnaseq_denovo_assembly_new_features' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       f249f9cf0 BFXDEV-397 rnaseq_denovo_assembly.py and rnaseq_denovo_assembly.base.ini rebased from master
       9639c9f50 BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       1922f9226 BFXDEV-397 tested filtered trinity output using real data (1 sample)
       afad4375a BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d223c7456 detected bugs during tests
       6d0e428a6  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       f91e42726 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       c4597fa93 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       b2afc9e00 detected bugs during tests
       bd91a5833  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       2939c0e27 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       5ef76ae63 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       641a6532c BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       cb1e022a8 BFXDEV-397 tested filtered trinity output using real data (1 sample)
       d9a9c464c BFXDEV-397 merging conflicting versions of rnaseq_denovo_assembly.base.ini and rnaseq_denovo_assembly.py files
       5c5b83e97 Merge branch 'rnaseq_denovo_assembly_new_features', remote branch 'origin' into rnaseq_denovo_assembly_new_features
       62acde9a7 BFXDEV-397 fixing differences when rebasing from master
       33d0dd0dd BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d3e9ff6be detected bugs during tests
       73ec34152  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       8a0bb017f BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       4d099e5a0 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       126102e7c detected bugs during tests
       f9d902ac1  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       0cb2d5322 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       5579bc4d0 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       ff1f6daef BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo

  lletourn <louis.letourneau@mail.mcgill.ca>      18 commits

       74503b2b6 Merge branch 'master' into highCoverageVariants
       bfdf02329 fixed io_buffer default
       8b2cfcc96 BFXDEV-392 First implementation of high depth calling
       d1249587f BFXDEV-370 Fixed merging and output naming bugs
       7381a036d Merge branch 'tumorPair' into highCoverageVariants
       da6e1fb27 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fffb12b62 BFXDEV-392 need varscan for high coverage
       4bdd55303 Merge branch 'master' into tumorPair
       8ae0dd3da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       23c3e596e set default ppn for igv to 2
       bba9c5659 BFXDEV-379 remove -gpfs from the ini files
       6cf6ac409 Version bump to 2.2.0-beta
       8e1345df7 BFXDEV-370 added merging step
       27983c2ff Merge branch 'master' into tumorPair
       61e35f357 BFXDEV-370 Added indels and COSMIC
       ffe5e6d7f BFXDEV-370 fixed scalpel script, added LD_LIB_PATH
       8c0cfc757 BFXDEV-370 Added the scalpel module
       7452be9b0 bvatools version bump

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      32 commits

       a57fd6d38 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       129d564cb RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       46958811c RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       226c59b40 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       1cbb0f175 chipseq - update module -  BFXDEV-490
       9e2448622 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       501c3ac02 DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       bc5920a45 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       849d732a3 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       db6c64a7a chipseq - update module -  BFXDEV-490
       65e38616d dnaseq - add correct dependency in metrics snv - BFXDEV-508
       e1637bf8a DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       c9e8db09b RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       0d077aa60 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       32d1d5bd9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4c43bb27c chipseq - update module -  BFXDEV-490
       7a6816025 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       7e2f61b36 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a95aa12b DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       58f9477e0 BFXDEV-465 - correct GATK DoC small bug to support when no bed fileis given in the readset file
       269b60fba DNAseq - Adding variant recalibration BFXDEV-436
       557f75ade DNAseq - gatk DoC will use the bed file as intervalls if the bed file is in the readset shett  - BFXDEV-465
       ce3014ac8 DNAseq - implement haplotype caller and mpileup annotation and filtering using old foinction as background - BFXDEV-463
       2aa000885 DNAseq - create new pipeline steps - BFXDEV-463
       f46876ac5 DNAseq - Add gvcf Combining the set of sample and genotyping - BFXDEV-440
       157b7e9db DNAseq - starting to implement GATK gvcf merging - BFXDEV-440
       3c407b263 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       74c6c2a32 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66cfdc682 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcfdd2c88 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c83d71ca8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1b1741151 bump pacbio module to patch 4 -  BFXDEV-415

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      105 commits

       db54b51a4 resolving conflict
       e7dc25160 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       3a2c22faf resolving conflict
       c554b5dc0 resolving conflict
       93a55f8e9 RNAseqDN - update module and remove trinity version check - BFXDEV-490
       48fafb907 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       2be21e182 Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       4103b9fb3 DNASEQ - update SnpEff command for versuion 4.2 - BFXDEV-490
       8e1bbc27f DNASEQ - support single end libray for metrics and dnasample metrics steps & allow filter_nstrech to use non recalibrated compress data - BFXDEV-503 BFXDEV-505
       4e0cd6cb4 resolving conflict
       228c23afb resolving conflict
       6e42b3df5 resolving conflict
       fea21f4ac RNAseqDN - update module and remove trinity version check - BFXDEV-490
       9cccf23b9 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       50f9346f3 Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       71c233cc4 DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       3a3c119ab resolving merging conflict
       a0491a4ff RNAseqDN - update module and remove trinity version check - BFXDEV-490
       ca27f3f1f RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       5b7e920f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7481dc791 Rnaseq_denovo_assembly - correct single library type bugs - BFXDEV-405
       b041650c5 DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       132e479e3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ffcee7d8c DNASEQ - support single end libray for metrics and dnasample metrics steps - BFXDEV-503
       3f22b9d99 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       10e4a0a07 updating modules (test on abacus human done) - BFXDEV-490
       5de58e204 rnaseq - support fr-stranded single end library in wiggle tracks and metrics - BFXDEV-499
       106788a2d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9b97e8853 rnaseq - make  wiggle tracks working for non UCSC genome ref - BFXDEV-498
       9bb4b278b correct the version of pandoc - BFXDEV-490
       e3f47cd7a update chipseq ini to latest version of modules - BFXDEV-490
       ae112c0fa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f7b90b570 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0e0ee575 update human genome sh files - BFXDEV490
       e3c632527 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       acc1d01e6 module - remove type in macs2 that make the module file invisible to the module system - BFXDEV-490
       4c413a261 module - add gemini and bwakit modules - BFXDEV-490
       c9a2f313d module - add star indexes for read length 75 and 150 bp - BFXDEV-490
       14464b164 module - adding a step print for the module creation - BFXDEV-490
       39f89bcdd module - add gemini and bwakit ; update bcftools htslib java picard prinseq-lite rnaseqc samtools - BFXDEV-490
       b50d111d0 BFXDEV-490 - updating modules gnuplot ; hmmer; igvtools; macs2; pandoc; python_lib
       61770a483 fixing conflict between master and highCoverageVariants branches
       a3491eb1f BFXDEV-490 - updating modules bedtools; blast; gatk; ucsc
       bd023c2a5 bfx/gatk.py - removing conflict to allow gatk_variants branch to be merged
       292d51caf Ressource - bump vcftools version to 0.1.14 - BFXDEV-489
       13aa92550 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5b8a705cd genomes - add rat Rnor_6 - BFXDEV-487
       9ef56ce60 genomes - install genome by default in dev and in prod if argument MUGQIC_INSTALL_HOME is given - BFXDEV-485
       4ac216ddc ressource - fix version of matplot lib to 1.4.3 for Qiime compatibility - BFXDEV-483
       28175da77 gatk_variants - danseq -  add mark_duplicates cleaning files - BFXDEV-471
       1c06f0c08 gatk_variant - module - bump mugqic_tools install script to 2.1.3 - BFXDEV-484
       62c067e56 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       da7245bf0 ressource - update python to 2.7.10 and add scikit-bio lib - BFXDEV-483
       7c774e6a3 remove tumor_pair from highCoverage
       5a20f1fbe DNAseq - add cleaning list to the new steps - BFXDEV-471
       9266e030f DNAseq - change destionation path in the copy of SNV metrics zip file - BFXDEV-473
       71b439e9e Pacbio assembly - modify whitelist option usage to do not generate additional ini and xml - BFXDEV-456
       1f4eb9d2d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be622168d Bump  module to the new smartanalysis patch 5
       15b99f439 gatk_variants - add baits intervals file for picard_hs_metrics - BFXDEV-467
       ad0a5c9c4 ingstats.py - Adding Array callrate check (Missing or Low) and add support for different manifest version through the -m option - BFXDEV-445 - BFXDEV-466 - BFXDEV-451
       75379ef31 DNAseq - synchronize code and ini and correct runnning issue - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       7a5290563 DNAseq - remove some smalls dev typos and bugs - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       a19894cfe PacBio - add require=False to the whitelist param - BFXDEV-456
       e8223e659 PacBio - incorpore whitelist option during intial filtering - BFXDEV-456
       53b2fbf10 PacBio - add addtionnal filtering config xml file for  whitelisted assembly - BFXDEV-456
       064ae2e6a PacBio - add specific ini file for whitelisted assembly - BFXDEV-456
       8e298ef41 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7769b771d update bwa install script to the newer version
       726137304 change star index place in resources/genomes/install_genome.sh
       0904455ef update install_genome.sh to use the coirrect version of R_packages
       a3e034012 bump bos_taurus install script to ensembl 80
       c57be7210 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       94223ff10 bump bos_taurus install script to ensembl 81
       17bea25ea make ignstats poiting on nanuq (not QC nanuq)
       382c9528f Merge remote-tracking branch 'origin/processingReport'
       051e37825 RNAseq - repair missing dependency - BFXDEV-427
       26d5d41ec debug paired_tumor scalpel vcf merging
       d2c011b8a removed merge conflicts
       29a319b54 DNAseq - regenarate updated README.md
       85017c1b2 DNAseq - change fixemate description text: picard -> bvatools
       ba184ccd2 High_coverage - needs dict2BEDs.py which is only in dev version of mugqic_tools - modifiy the ini file
       3d887984f DNAseq - change post fixmate sorting to generate and indexed bam - BFXDEV-421
       83e7e1614 modify igvtools lib to allow modifiy paramter through the ini file - BFXDEV-419
       d5d86bf80 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       71c3d26f8 DNAseq - change coverage metrics report to remove the CCDS part of the table - BFXDEV-381
       0de4e6565 rnaseq add bam_hard_clip step : Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature - BFXDEV-377 + change rRNA metrics output to avoid name conflict - BFXDEV-401 - BFXDEV-407
       8a6368104 update python and python lib to include Qiime ligth version installation - BFXDEV-412
       8bb30eea5 bump trinity install script to version 2.0.6
       1b2a34dba rnaseq - remove tophat/botwie commented section and modules
       f53ed32c4 rnaseq_denovo_assembly - correct single library issue and add in comment the correponding change in the base.ini - BFXDEV-405
       d37b2c89f nanuq2mugqic_pipelines.py - support nanuq group info for Miseq/Hiseq/Pacbio technology - BFXDEV-418
       b03537ada update install module general script
       2d1dbe489 bump smrtanalysis module to 2.3.0 patch 4
       27bad469e update install module general script
       006f59cfb update install module general script
       08acf20a7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcb95896f bump smrtanalysis module to 2.3.0 patch 4
       ba4b9f140 update bvatools module for release 1.6
       c72acd4b1 core/pipeline.py remove duplicates sections in reports - BFXDEV-388
       4cb0ba7a7 ChipSeq - remove expected input files for broad peaks when generating homer annotation report - BFXDEV-393
       797c96812 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3a54c8564 modify resources/modules/smrtanalysis.sh for patch 3
       ce104cd10 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d309b7da9 RESSOURCES - R_bioconductor - Addthe package 'circlize' to the list of R package to automatically install

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      6 commits

       5db43a8b0 Merged in irp_end_fastq_notification (pull request #14)
       de7d7f847 Merged in irp_bcl2fastq2 (pull request #13)
       d3acddd6e Merged in illumina_run_processing_sprint (pull request #11)
       82747c860 Merged in pacBio_whitelist (pull request #7)
       e16e90df1 Merged in ignstats (pull request #8)
       8e2cf4c43 rnaseq.base.ini increase star io limit to 4G by default

  mmichaud <marc.michaud@mail.mcgill.ca>      62 commits

       4d820a7cd BFXDEV-504 - Notify nanuq when fastqs are completed for a lane (fix url) - illuminaRunProcessing
       7ab569eee BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       cdafb2390 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       fd8579e95 BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       7282b301a BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       2d9a8965e Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       e0be33183 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       6e04c6294 BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       34ca19a2f BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       09ff7003b Increase resources for qc.
       3dfce03b6 BFXDEV-488 Illumina Run Processing: Use bcl2fastq2. Fix for dual-indexing.
       5e8133ee2 Use 2 core for rsync.
       48a4bd862 Merge branch 'master' into irp_bcl2fastq2
       a9560d4f6 Fine tune cluster settings according to the history of the jobs of the last three months.
       b917ee551 BFXDEV-413 Fix code according to code-review.
       a8fd8bfe3 BFXDEV-488 Illumina Run Processing: Use bcl2fastq2.
       be6ec3f13 Merge branch 'master' into illumina_run_processing_sprint
       e0c3e20b6 BFXDEV-468 VerifyBamId: Fix job name when not using bed file.
       ec19be9d9 BFXDEV-468 VerifyBamId: Run even we don't have a bed file.
       1ffafbb6c BFXDEV-468 VerifyBamId: Use a version supporting "chr" in chromosome name.
       dec09d848 BFXDEV-468 VerifyBamId: Fix output for nanuq.
       eed1cc169 BFXDEV-468 VerifyBamId: Shorter job name.
       750495f7b BFXDEV-468 VerifyBamId: Properly concat jobs
       ce01a76f3 BFXDEV-468 VerifyBamId: Nanuq friendly output.
       9d5fa0e47 BFXDEV-468 VerifyBamId: Don't use a compressed VCF + Nanuq friendly output.
       c4bd81b15 BFXDEV-468 VerifyBamId: Add missing configuration.
       72c77d237 BFXDEV-468 VerifyBamId: Fix annotation file name.
       000109977 BFXDEV-468 VerifyBamId: Fix wrong genome for sample.
       15a3a2bca BFXDEV-468 VerifyBamId: Changes to the annotation filename.
       7fc91bd9a BFXDEV-468 Don't run verifyBamId when there is no "dbnfp_af_field" on the genome.
       d747ef6f4 Illumina Run Processing: Fix barcode counter path with cvmfs
       68c008938 BFXDEV-468 IRP - VerifyBamId: Optional dbsnp_version and dbnsfp_af_field in genome ini file.
       5a3fd82b1 Use cvmfs blast and bwa.
       75af739c8 BFXDEV-468 Add verifyBamID in Illumina Run Processing
       94ccde0c3 Merge branch 'master' into illumina_run_processing_sprint
       7d6cf517c BFXDEV-447 Add rRNA estimate using silva blast database. (Fix change output file name)
       dd9f3f7a8 BFXDEV-464 Don't use an hardcoded genome list for the alignment: parse the value from the sample sheet and validate that we have the genome on disk.
       015435eda BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name (fix redo of pipeline always restarting the STAR align).
       f0ecd6394 BFXDEV-447 Add rRNA estimate using silva blast database. (Fix error when creating result file)
       2d39be9bb BFXDEV-413 Add Estimated Genome Size from the nanuq CSV.
       c73b3fc9d BFXDEV-447 Add rRNA estimate using silva blast database.
       fdcb45eed BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name.
       1819bf2cb Add sample sheets examples and add description about "Genomic Database"
       be10d09b1 BFXDEV-457 BFXDEV-459 Update documentation.
       b337a061d Don't warn when skipping alignment on a sample without "Genomic Database"
       a73ca70a3 BFXDEV-459 Merge start_copy and copy steps.
       583237b9a BFXDEV-457 Add the possibility to regroup all the md5 jobs into one (add "one_job=1" in "[md5]" config section).
       56cdf2d33 BFXDEV-458 Fix mask calculation problem on lanes with nextera/truseq samples.
       041d4e816 More useful usage message.
       3cd065cc0 Merge branch 'master' into irp_genomic_database
       d90c137d7 Merge branch 'master' into irp_genomic_database
       2f0ef65a3 Illumina Run Processing: Increase memory for BVATool DoC and Readsqc
       a864f3182 BFXDEV-369 Illumina Run Processing: Fix GRCm38 genome detection.
       7eae59d77 BFXDEV-369 Illumina Run Processing: Use the reference specified in the request form
       64fb218e2 Increase processors per node for Fastq and bwa alignment. Ignore setfacl exit code for the copy step.
       fe4457726 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0) - Fix passFail threshold according to the new value.
       cafa62c30 BFXDEV-417 IGNStats: Output Array Call Rate as ratio, not percentage (0.95 not 95.0%)
       eb810e069 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0)
       37245835d BFXDEV-417 Modifiy IGN stats parser to allow the upload of data to a nanuq server.
       ce5285496 Use 13 core for BWA job to avoid excessive memory usage.
       8505fab84 BFXDEV-385 Fix bed files handling in Illumina Run Processing.
       29d5c0f60 BFXDEV-380 Change permissions of the run files on the destination.

  noreply@clumeq.ca <libuser@lg-1r14-n01.guillimin.clumeq.ca>      15 commits

       ad55029f4 genomes - updated bash script for human genome GRCh37 and GRCh38  - BFXDEV-490
       4d653bf93 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       d2fd863a5 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       e4b6339f1 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       1c79f7e84 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       36b260f67 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       565502762 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       7e8870a0c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0a73d79f3 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eb3010723 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8c32446d4 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0b9b70c62 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0233d7024 mugqic_tools is now up to date && resources/modules/dev/epacts.sh has been removed
       64c821d52 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a65de9e00 temporary modify install_genome.sh to run the job in batch (job submission is blocked

  ptranvan <patrick.tranvan@mail.mcgill.ca>      1 commits

       2a1917b25 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      4 commits

       f18496611 Addition of htslib module and resolved pull request comments
       457c3695f addition of vcf preprocessing, snpeff, gemini and ini adjustments
       09ad7f1f9 addition of vcf preprocessing, snpeff, gemini and ini adjustments
       05f05b228 Merge branch 'highCoverageVariants' of https://bitbucket.org/mugqic/mugqic_pipelines

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      2 commits

       60fd23052 Merge branch 'highCoverageVariants' of bitbucket.org:mugqic/mugqic_pipelines into highCoverageVariants
       76802ee0d Updates to pairedTumor: addition of CombineVariants

2.1.1        Mon Apr 13 22:23:46 2015 -0400        172 commits

  Francois Lefebvre <lefebvrf@gmail.com>      4 commits

       65710dc83 QUAST and Minia dev install scripts
       5b404f38a Updated kmergenie version in install script and moved to dev
       6db368e9a Added NxTrim mugqic_dev install script
       45b9fb07b Updated a bunch of module dev install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      109 commits

       4db053f68 Added Gorilla_gorilla.gorGor3.sh in resources/genomes/old/
       deb0b822a Changed report job names with '_' instead of '.' to avoid scheduler cpu setting > 1
       4d8b38a77 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b471f6f43 Fixed bug in RNA-Seq De Novo Assembly : use RSEM plugin in Trinity instead of external one
       6c79dd4c0 README.md edited online with Bitbucket
       2e53a7659 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       455a50c6d Fixed missing pandoc module in rnaseq
       0ffb7eadc Added variables total and average read length in pacbio assembly stats report + changed locale from fr_FR to en_CA
       70150e791 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       121faa00e Added sample_cutoff_mer_size to pacbio report job names
       fcb6b897f README.md edited online with Bitbucket
       272830bee README.md edited online with Bitbucket
       89c55cbc9 README.md edited online with Bitbucket
       8f6eeb79a Fixed cpu bug for blastp_transdecoder_uniprot job in RNA-Seq De Novo assembly
       cae10e5de Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       28eb9540e Minor comment updates
       2c4a463a5 README.md edited online with Bitbucket
       c4262fa99 README.md edited online with Bitbucket
       3f664a577 README.md edited online with Bitbucket
       dc8ec12ab README.md edited online with Bitbucket
       c60d74246 README.md edited online with Bitbucket
       b38e19507 README.md edited online with Bitbucket
       6d4206c1b README.md edited online with Bitbucket
       4a6982c61 README.md edited online with Bitbucket
       39ee0b9b0 README.md edited online with Bitbucket
       1f13501a1 README.md edited online with Bitbucket
       0376786ec README.md edited online with Bitbucket
       43d0ac6dd README.md edited online with Bitbucket
       be50cd0ac README.md edited online with Bitbucket
       f85506e81 README.md edited online with Bitbucket
       5e335b491 Added GNU Lesser General Public License (LGPL) to MUGQIC Pipelines
       0b79b46bd  BFXDEV-59 Updated ChIP-Seq report redesign with sample metrics without trimming
       ec88ea54f Fixed bug cluster_cpu for blastx_trinity_uniprot
       5962a9fe3 Removed UniRef BLAST in RNA-Seq De Novo Assembly since it is too long; factorized differential expression code; renamed design variable into contrast
       611223e9e Adjusted Trinity butterfly cluster resources for RNA-Seq De Novo Assembly on abacus
       371d60c68  BFXDEV-59 Completed ChIP-Seq report redesign
       7eedd054f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       01a41aac0  BFXDEV-59 ChIP-Seq report redesign
       f8bf89fad BFXDEV-59 Completed PacBio Assembly report redesign
       479ea0b0a BFXDEV-59 Differential expression RNA-Seq De Novo Assembly report redesign commit
       76c9bc36e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       c7a9036f8 BFXDEV-59 More RNA-Seq De Novo Assembly report redesign commit
       ac2daf387 README.md edited online with Bitbucket
       e7f156a51 README.md edited online with Bitbucket
       3caf504bb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       6771f039c BFXDEV-59 First RNA-Seq De Novo Assembly report redesign commit
       5c858ad99 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5d82ab02a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       94e9930d0 Decreased rnaseq cufflinks default pmem to 2700 for guillimin
       79ea01b15 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       80cfde80d Increased rnaseq cufflinks default pmem to 3700 for guillimin
       9e655b98d Fixed BED file abspath in bvatools
       c183fadef BFXDEV-59 Completed RNA-Seq report redesign
       61676ebb2 BFXDEV-59 Added metrics steps for RNA-Seq report
       e760ac7ee BFXDEV-59 Fixed merge conflict
       fbedae5ed BFXDEV-59 More and more commit for partial HTML report
       7a95c8442 Increased dnaseq compute_effects ram
       c6c7ead0e BFXDEV-59 Even more commit for RNA-Seq report
       75b2a7ee3 BFXDEV-59 Even more commit for RNA-Seq report
       5da988599 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       0463d0097 BFXDEV-348 Fixed with readset name = <sample>.<library_barcode>.<run>.<lane>
       3072da14b BFXDEV-59 Even more commit for RNA-Seq report
       57f2aef6b BFXDEV-59 More commit for RNA-Seq report
       ad7e670e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be924071b Updated mugqic_tools to version 2.1.0
       672b68143 BFXDEV-59 First commit for RNA-Seq report
       f194ec395 BFXDEV-59 DNA-Seq report minor fix
       7f98cbfeb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       40a04d08f BFXDEV-59 DNA-Seq report redesign done + fix
       b4b408f7e BFXDEV-59 DNA-Seq report redesign done
       426e8788f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       a0b4223d1 Minor doc fix
       bc1d662e9 BFXDEV-59 Even more report redesign commit
       09de8ac43 BFXDEV-59 More report redesign commit
       91e8e95b4 BFXDEV-59 First report redesign commit
       77f66c40d Added UCSC genomes in install_all_genomes.sh
       aa6a8e0a0 Create genome rrna.fa with grep -i 'rrna'... + remove variation sequences containing '.' in Ensembl vcf.gz which make GATK crash
       525faa860 BFXDEV-295 Updated mugqic_R_packages to 1.0.3 for RNA-Seq De Novo pipeline
       da8253e47 BFXDEV-295 minor fix for RNA-Seq De Novo pipeline
       e161ad93e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4d177beb4 Updated Python pysam to version 0.8.2
       e7ca27f56 README.md edited online with Bitbucket
       48073caf6 README.md edited online with Bitbucket
       70a0304b6 README.md edited online with Bitbucket
       1e24257cc Increased cores for homer_annotate_peaks in chipseq
       3c2c24def Fix new section names for blast on uniprot in RNA-Seq De Novo pipeline
       e699c6427 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       aeedee41f Added chipseq in PATH of module mugqic_pipelines
       7ce201c60 Added ccds filepath in Rattus_norvegicus.Rnor_5.0.ini
       7b08143be BFXDEV-295 Moved trinotate Pfam + UniProt DBs in /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod/genomes/[blast|pfam]_db/ for RNA-Seq De Novo pipeline
       c47a5190e BFXDEV-295 Update RNA-Seq De Novo pipeline with modules in prod and trinotate updated
       33812be2e Minor mammouth adjustment regarding increased ram and core in snp_effect job for DNA-Seq Pipeline
       9689e6dbe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c345049c7 Increased ram and core in snp_effect job for DNA-Seq Pipeline
       8743bbb91 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f5b42293 BFXDEV-295 Updated RNA-Seq De Novo assembly pipeline with Trinity 2.0 and Trinotate 2.0
       b484ff7f9 README.md edited online with Bitbucket
       4b19c2007 README.md edited online with Bitbucket
       1b0f09cd5 Minor fix in python lib install
       918061242 Separated python install and python lib install + minor weblogo install update
       e06560dd8 Minor fix in Perl lib install
       b202849f7 Fixed missing bvatools.depth_of_coverage other_options + snpsift_annotate module_snpeff=mugqic/snpEff/4.0 + minor uppercased '.insert_size_Histogram.pdf' for picard.collect_multiple_metrics output file in dnaseq pipeline
       f8b35eb3c Fixed missing 'nodes=1' in pacbio_assembly.guillimin.ini smrtanalysis_run_ca
       e7d4a0e1c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       62a0e0059 Updated bedtools to version 2.22.1
       9328bc99b Standardized perl module install and moved CPAN libs in a different script
       0fb4c648e BFXDEV-335 Removed adapter FASTA files except adapters-fluidigm.fa
       4aba5c164 BFXDEV-335 Create adapter FASTA from readset file for Trimmomatic, if not defined in config file
       ba787480d Version bump to 2.1.1-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      14 commits

       1792f2c8f Version bump to 2.1.1
       881082206 BFXDEV-375 Fixed ram sorting issues when using star
       e9fcf5434 Merge branch 'master' into rna_metrics
       a7e3465c8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ba666200 BFXDEV-368 IGN script to extract stats
       8f6914002 Fixed code for tools that need to be downloaded manually like gatk
       f6ee04ebb Updated gatk
       4622df3b1 Updated samtools, bcftools, htslib
       6eac31587 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a44561f3f BFXDEV-351 removed md5 from markdup and added it to recal
       8abc6ccf2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1db78f982 BFXDEV-355 Removed CCDS from the options BFXDEV-356 added compression to haplotype caller output
       9a0a69848 BFXDEV-351 Removed MD5 from markdup, added it to recalibration
       b33195c0a BFXDEV-346 Split jobs in a more uniform way

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       70099f3fc rnaseq - include correlation matrix manual estimation using an utils.py new function
       fd5380189 RNAmetrics - include ini_section argument to some bfs picard function to allow to use them several time with different parameter in the same ini
       e7c5a010b rnaMetrics - remove confilct pipelines/rnaseq/rnaseq.base.ini
       9ad7fbb1b rnaMetrics - update pipelines/rnaseq/rnaseq.py pipelines/rnaseq/rnaseq.base.ini
       9d3b87da8 RNAseq - metrics update ini
       1913bc103 RNAseq - remove conflict
       deca058d5 RNAseq - metrics RNA - update base ini
       2fc2006be RNAseq - remove rnaseqc; add picard_rna_metrics ; partial add estimate_ribosomal_rna - BFXDEV-345

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       05937d609 RNA-seq - rna_metrics : test are ok; new files annoation files are created ; point to the production assembly folder -  BFXDEV-345
       b2c948b88 ressource -Genome install script add the creation of the ref_flat file format of the annotation from the gtf and correct a path in the GRCh38 file (All.vcf.gz) - BFXDEV-374
       7d8f72fc7 RNAseq- update module version not found in CVMFS (bowtie, bvatools, mugqic_tools) - BFXDEV-373
       e57ed5edc RNAseq- remove redundant step in the step initialization - BFXDEV-371
       3138867ba release new version of resources/modules/mugqic_tools.sh
       9aa58f7ae RNAseq - rna_metrics - finish the rRNA and rna metrcis modifications - BFXDEV-345
       cd5f2391f RNAseq -rna_metrics -fit the new bam2fq from bvatools_dev
       d2048ebe4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rna_metrics
       d651feeb7 RNAseq - rRNA metrics add bvatools | bwa | picard && rrnaBMAcount.py

  mmichaud <marc.michaud@mail.mcgill.ca>      27 commits

       2a821a4b8 bcl2fastq module name change.
       00b481f57 Merge branch 'master' into irp_rna_metrics
       c2c33a3f2 MPS-1740 Use updated production genomes.
       5f3c74e5c MPS-1740 Use released version 1.5 of bvaTools
       d191969e0 Increase STAR sort memory
       dd936d0af BFXDEV-363 IRP: Don't copy phasing and matrix files.
       ae9713bce BFXDEV-353 IRP: Standardize job name
       e83e25fa2 Merge branch 'master' into irp_rna_metrics
       87aaf0db4 Merging rna_metrics on irp-rna_metrics
       b0dff3b12 BFXDEV-353 Use dev version of mugqic tools
       f94c13578 Merge branch 'master' into irp_rna_metrics
       75398a389 BFXDEV-353 Use new version of rRNABAMcounter.
       5f443d4b6 BFXDEV-353 Use a nanuq friendly name for the rRNA metrics file
       ea661af49 BFXDEV-353 Add the rnaseqc '-rRNA' option to set an optional ribosomal rna interval list. Setting an empty file will skip the rRNA count.
       91c0d87d5 Simplify illumina run processing RG tag logic. As library, run and lane are mandatory there is no need to validate that they exists.
       4abaa7b14 BFXDEV-353 Fix bwa other_options
       81e7423bb BFXDEV-353 Use DEV versions of bvatools and genome. Update mugqic_tools to 2.1.0.
       79a02a723 BFXDEV-353 Update rna-seq metrics according to rna-seq pipeline
       081c6b814 Merge branch 'master' into irp_rna_metrics
       fcc8392b5 Code format
       d62f8550c Merge branch 'master' into irp_rna_metrics
       0bc4cbe70 Add missing thread parameter for bvatools_depth_of_coverage.
       221eaae0c BFXDEV-339 Use rrna file in rnaseqc
       8e95c9943 BFXDEV-339 Use rrna file as ribosomal annotation (instead of ncrna)
       5efa1042d BFXDEV-338 New Nanuq MPS run association
       7aab26f95 BFXDEV-338 New Nanuq MPS run association
       6e6de0755 BFXDEV-338 Add run_id in all commands available parameters

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       559973d83 correction other_options

2.1.0        Wed Feb 4 14:40:16 2015 -0500        123 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       eefcc474f Added qualimap installa script + updated population
       d9397b016 Changes to sailfisj , R install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      104 commits

       09344ee60 Version bump to 2.1.0
       b3e856181 Changed procs= to nodes=1:ppn= in rnaseq and rnaseq de novo guillimin config
       a63500dbc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a195837e6 Added libgd in PacBio mammouth config
       a863cdd5e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e13ac413c Added pacbio_assembly.mammouth.ini + removed unnecessary wgs module
       1cb2ea5df Added pacbio_assembly.guillimin.ini + adjusted guillimin cluster settings
       b9f126420 README.md edited online with Bitbucket
       ff58b26e4 Removed 'daemon' from scheduler options
       2d4459d78 BFXDEV-292 Removed optional trimming dependencies for report in chipseq pipeline
       6b9cdad56 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f40f39acc BFXDEV-292 Fix metrics and report dependencies in chipseq pipeline
       c1fb3245f README.md edited online with Bitbucket
       d8ade06f6 Removed tmp_dir validation since some directories are available on exec nodes but not on login nodes
       44c931e08 Removed memtime from PacBio pipeline + updated smrtanalysis to version 2.3.0 patch 2
       698d3f701 BFXDEV-292 Added cleaning in chipseq pipeline
       cac630551 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       086ee05ed Updated smrtanalysis module 2.3.0 with patch 2
       0c6d970fc README.md edited online with Bitbucket
       f808abe06 Updated ChIP-Seq README + main READMEs
       e789ffd63 Added BiSNP module install script
       26a062937 Update mugqic_tools to 2.0.3 in rnaseq config
       7af099daa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d10a05841 BFXDEV-292 Fix metrics tmp file preventing job to be up to date + module ImageMagick to use convert command on mammouth
       367f86eb7 Removed '\' before /ltmp/fillon in mammouth config files
       0c1988f05 Added picard tmp_dir type validation
       cf8206563 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0884bf53f BFXDEV-292 Added config files for all clusters in chipseq pipeline
       2a5721520 Updated mugqic_tools install to version 2.0.3
       387c225b6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       5017938f4 BFXDEV-292 Last steps implementation in chipseq pipeline
       52d1f782d BFXDEV-292 Added UCSC genomes hg19, mm9, mm10, rn5
       751e701af BFXDEV-319 Adjusted cluster settings in rnaseq.mammouth.ini + .ini absolute path for report
       edc26514d BFXDEV-35 Standardized memtime module install
       90f4aa5c6 Minor update in README-release.txt
       ab65f5da2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       f37bed0a4 BFXDEV-292 Set peak annotation with internal homer genome in chipseq pipeline (incomplete)
       9b3560281 BFXDEV-35 Moved wgs-assembler module install in dev
       77859cc24 BFXDEV-35 Standardized vcftool module install
       a07ea0b5c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       2ee431548 BFXDEV-35 Standardized ucsc module install
       284138b55 BFXDEV-35 Standardized trinotate module install
       9f31173a5 BFXDEV-35 Standardized trimmomatic module install
       1fab82a12 BFXDEV-35 Standardized tophat module install
       00b9355bb BFXDEV-35 Standardized tmhmm module install
       2917504dc BFXDEV-35 Standardized tabix module install
       6c1ec5b53 BFXDEV-35 Standardized star module install
       55547462f BFXDEV-35 Standardized snpEff module install
       25f53a1c2 BFXDEV-35 Fix smrtanalysis archive exec permission
       275894fb1 BFXDEV-35 Standardized smrtanalysis module install
       6f3b27ab8 BFXDEV-35 Standardized signalp module install
       a335bfa58 BFXDEV-35 Standardized samtools module install
       7cfda90da BFXDEV-319 Added rnaseq.mammouth.ini
       dac40de25 BFXDEV-319 Removed java option -Dsamjdk.use_async_io=true in all .ini files
       e394bf78c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       705f2c6bc BFXDEV-292 Adjusted cluster settings in chipseq pipeline
       a87203d36 README.md edited online with Bitbucket
       2524cb9d2 README.md edited online with Bitbucket
       af39fb42a README.md edited online with Bitbucket
       e180a8387 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0daae6ed5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       93c1a2883 Fixed bam sub arguments in dnaseq bwa_mem_picard_sort_sam
       934b42f67 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       a5931e42a BFXDEV-35 Standardized rsem module install
       dbaebd1ee Moved repeatmasker module install script in dev
       0c30ff422 BFXDEV-292 Minor docstring change in chipseq pipeline
       dc0b98886 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ef6c09e41 BFXDEV-292 Added GTF for homer_annotate_peaks in chipseq pipeline
       26894dccd README.md edited online with Bitbucket
       c7cdd1aef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       4f0b08581 BFXDEV-292 Fixed input/output files bugs in chipseq pipeline
       4c1cf7d27 BFXDEV-35 Standardized weblogo module install
       bf36cc0e3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       7952a5fb4 BFXDEV-292 homer_find_motifs_genome and annotation_graphs steps in chipseq pipeline
       b29cb3d91 BFXDEV-35 Standardized rnaseqc module install
       05208ed06 BFXDEV-35 Standardized rnammer module install
       e33776a58 BFXDEV-35 Standardized more prinseq-lite module install
       fcc8bc9e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d15d9b003 BFXDEV-35 Minor fix in AMOS module
       07bba1e88 BFXDEV-35 Fixed path bugs in MUMmer and AMOS module install
       99363eed9 BFXDEV-35 Standardized MUMmer module install
       423c3c7c8 BFXDEV-35 Minor fix in module install template
       c00bdfcf0 BFXDEV-35 Standardized java module install
       d185a5072 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       8bed5d97a BFXDEV-292 homer_annotate_peaks step in chipseq pipeline
       0bca1960b BFXDEV-35 Standardized igvtools module install
       1581fc97b BFXDEV-35 Standardized hmmer module install
       2acdc4926 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       cddbce73e BFXDEV-292 More macs callpeak in chipseq pipeline
       ed7cf4ec6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b2d02b10b Fixed macs2 with generic shebang
       cf61d7e0e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       bbe5703e6 BFXDEV-292 Fixed missing samtools module for MACS2
       ff1a4ec95 Standardized gnuplot module install
       dcdc21141 Standardized exonerate module install
       5c4774625 Minor change in README-release.txt
       076c02b67 Version bump to 2.1.0-beta
       953894cb4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       1090d3830 More chipseq deelopment
       0841b02b7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ace0e308c BFXDEV-292 Beginning of macs2_callpeak step in chipseq
       b404e2226 Added qc_plots_R step in chipseq
       43e70dcef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       51f34444b BFXDEV-292 First draft of chipseq pipeline

  lletourn <louis.letourneau@mail.mcgill.ca>      3 commits

       9f9f60cf9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d2d65f01c BFXDEV-327 Used only one thread for haplotypecaller because of a race condition
       639f65027 BFXDEV-327 Used only one thread for haplotypecaller because of a race condition

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       12ff20e98 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       233244b40 RNAseq - correct stdin input issue of htseq-count - BFXDEV-318
       fa335288e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       47c4cc811 RNAseq - htse-count: pipe samtools view -F 4 output in htseq-count instead of using the bam to remove error due to unmapped reads - BFXDEV-318

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       7ec034492 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1905eb60a COMMON - add fastq2 = None in sam_to_fastq pipeline step wehen the read are single - BFXDEV-321

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       953de3d8f BFXDEV-332 Fix depth_of_coverage 'other_options' by removing extra quotes
       abca18089 BFXDEV-332 Base the jobs configuration on the perl pipeline
       619c0760b BFXDEV-329 Add run and lane number to the alignment and metrics jobs name
       6f4440f31 BFXDEV-202 Increase RAM for barcode counting
       eac645f2b BFXDEV-328 Increase mem walltime to 48h (and set ram for SortSam)
       ccb73f33f BFXDEV-331 Run Processing: Use a RG id that is unique across multiple lanes
       c2218dafc BFXDEV-317 Run processing: Use the old name for the coverage and onTarget metrics
       11ec9e77c BFXDEV-316 Fix errors when using a custom Illumina sheet file

2.0.2        Mon Jan 12 16:56:16 2015 -0500        80 commits

  Joël Fillon <joel.fillon@mcgill.ca>      54 commits

       104778ab3 Version bump to 2.0.2
       bbdf438cf Version bump to 2.0.2
       9d984cecf Fixed cluster resources in rnaseq_denovo_assembly.mammouth.ini
       ee96fe3bb Minor fix in pacbio_tools_split_reads config section name
       8831a3ae0 Even more standardization of module install
       a2c582486 More standardization of module install
       2b334e015 Standardized mugqic_pipelines module install
       395d3ec60 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66230bd2d Standardized mugqic_tools module install
       d9a56ab7e Standardized trinity module install
       0549be9dc Standardized cufflinks module install
       ddbd8f2ee Standardized MACS2 module install
       386921ec4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9ba8d3907 Standardized picard module install
       55a5ccb76 Standardized cd-hit module install
       d693b9a2b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f05fe1b1 BFXDEV-310 Updated picard to version 1.123
       a118e0d9f Standardized bwa module install
       4dc7bb7ac Updated bvatools to version 1.4 for dnaseq and puure
       cca258b1a Standardized bvatools module install
       94910f99c Standardized bowtie2 module install
       3e8f20164 Standardized bowtie module install
       1b70a257e In picard_sam_to_fastq, updated skip test if FASTQ column present and BAM missing
       345909af5 BFXDEV-290 Added job_input_files, job_output_files in JSON export
       f6e623e73 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       39e9ebfa6 More JSON export for Daemon Scheduler
       729567194 Standardized blast module install
       947111a8c Standardized bedtools module install
       9325049b7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c52d27755 Standardized GATK module install
       a69bd8207 More module install generalization; first test with prinseq-lite
       d116ebd60 Minor wget output file fix in prinseq-lite.sh and module install template
       c30ad6656 Added prinseq-lite module install script
       b2059f5ff Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ca901772 Fixed pb wget --spider in install_genome.sh + various up-to-date checks in Homo sapiens assemblies install
       21e59af06 BFXDEV-284 In nanuq2mugqic_pipelines.py, use native httplib, retrieve BED files and create adapters file
       e9620ff04 BFXDEV-290 First draft of daemon scheduler
       b778e5dbf BFXDEV-74 Finished cleaning in rnaseq denovo assmbly pipeline
       a171581ed BFXDEV-74 Started cleaning in rnaseq denovo assmbly pipeline
       c1cee2951 Reorganised lib functions more compact
       1cc8e33b8 BFXDEV-74 Added cleaning in rnaseq pipeline
       13475e2a1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       84954f327 BFXDEV-74 Added cleaning for dnaseq pipeline + reformatted various bfx modules
       1948d9d61 README.md edited online with Bitbucket
       848a0abbd Fixed bug use lstat instead of stat to check job up-to-date status without following symlinks
       fdfcbdde5 Minor aesthetic updates on pipelines doctrings + README.md
       da66a410b Generated pipelines README.md from docstrings using --help
       89718b0ac Added steps docstrings in RNA-Seq De Novo Assembly Pipeline
       ab5e595db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       232bed51f BFXDEV-296 Added steps docstring for RNA-Seq pipeline
       8d6f9a087 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b665c4cc5 Minor release doc update
       bf9afb6ca Update mugqic_pipelines module to 2.0.1
       530edcf45 Version bump to 2.1.0-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       ccb812e84 Fixed picard installer, reverted back to 123

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       fcd7134d7 pull before pushing Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       05b352d8b  PacbioQC - Fix code vs ini section name - BFXDEV-315
       d08526903 RNAseq - correct discrepency in hsteq_count ini calls - link to BFXDEV-312
       6773ef753 pull before pushing
       d5d82f0de RNAseq - correct errounous section header in ini files - BFXDEV-312
       98a456abf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fdd92aa23 RNAseq -update guillimin cluster ini file - BFXDEV-307
       b4a12a5a5 RNAseq - fix report dependencies - BFXDEV-306

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       693ea1dc6 PacBIo - add input files to the first pacBio job file (cp job) - BFXDEV-305
       ce7f0b369 Revert "PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305"
       fc6eb26f3 PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305
       162af5e5b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2ff273beb RNAseq - modify docstrings - BFXDEV-303
       c58b6d3dd DNAseq - correct report  dependency (typo) - BFXDEV-304
       ab81bc347 DNAseq - correct report  dependency - BFXDEV-304
       12830d322 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b6f4e8fd6 PacBio de novo - to allows running multiples sampes in parralele in the same analysis - BFXDEV-301

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c14a28187 BFXDEV-314 Illumina Run Processing: Output nanuq friendly metrics files
       3141cccc2 BFXDEV-313 Run processing: Don't depend on the start fastq notification
       5956b717e BFXDEV-310 Use Picard 1.123 with a bug fix in MarkDuplicate
       d01887933 BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       793e91cc7 BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       000d043e7 BFXDEV-309 Use the hiseq qos by default on abacus
       91cfb2627 Run processing: Check the library type when the library source is 'library' to determine if the sample is from RNA
       94f196eb7 Run processing: Fix configuration for miseq (copy's destination folder)

2.0.1        Wed Dec 17 09:56:23 2014 -0500        33 commits

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       76a44bc1f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fe96c5809 no message
       0f6a3ec68 Modified R_Bioconductor script to accommodate older gcc's

  Joël Fillon <joel.fillon@mcgill.ca>      22 commits

       bc72fc836 Version bump to 2.0.1
       5868b8573 Updated mugqic_R_packages to 1.0.1 and mugqic_tools to 2.0.2
       3230102ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd49feddf Updated /mugqic_tools.sh to 2.0.2
       2b979a2d4 Fixed dnaseq species_vcf_format_descriptor with absolute path (no /sb/programs/analyste)
       51fef9a14 BFXDEV-299 Fixed bug GATK realign with unmapped parameter not skipping other chromosomes
       6890d45e3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9a53a1158 BFXDEV-296 Added DNA-Seq step docstring documentation
       a820bbc29 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4475a6b97 BFXDEV-296 Added detailed --help option with output in markdown format
       cb49e7b10 Added detailed --help option with output in markdown format
       bf7bb0ab0 More PacBio README.md update
       cd046510b README.md edited online with Bitbucket
       93e28c7cb README.md edited online with Bitbucket
       48619e05f README.md edited online with Bitbucket
       2e79a5519 README.md edited online with Bitbucket
       3a17272ef PacBio Assembly README.md generated by --help
       93e81247e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       69073ee65 Trinity + normalization settings for guillimin in rnaseq_denovo_assembly pipeline
       8d8571f3e BFXDEV-26 Added cluster_max_jobs param in config files + added warning in PBS scheduler if this limit is reached
       14194c6e9 Minor update in README-release.txt
       38f24d73d Version bump to 2.1.0-beta

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c661248df Run processing: Copy the rnaseqc metrics file to the bam directory to be compatible with nanuq
       e9ecd7ab1 BFXDEV-165 Add missing Perl module for bcltoFastq
       44c06a01b BFXDEV-296 Changes to run processing documentation to use the generic help
       ecb244735 Run processing: Put rnaseqc results in a specific folder by librairy
       336a26521 Run processing config: Only send email on abort and increase cluster wall-time to 48h for the fastq job
       3e8733a7f Run processing: Add symbolic link to STAR created bam to fix output file check on re-run. The original file is renammed to follow nanuq naming conventions
       7e8dddbbb Run processing: Fix rnaseqc sample list
       7e38acffb BFXDEV-297 Add RNA-SeQC in run processing pipeline

2.0.0        Thu Dec 11 18:12:54 2014 -0500        669 commits

  Eric Audemard <audemard@ip03.m>      1 commits

       dc40fe20e update mammouth .ini

  Eric Audemard <eaudemard@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       62e246479 fix bug to find lib in pipeline python

  Eric Audemard <eaudemard@lg-1r17-n02.guillimin.clumeq.ca>      1 commits

       02b7ea60a fixe bug on puure before Ray execution

  Eric Audemard <eaudemard@lg-1r17-n03.guillimin.clumeq.ca>      2 commits

       3165f9aa2 Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       bbe9d0108 Bug fixed. Pipeline tested and validated on guillimin (step 1 to 21)

  Eric Audemard <eaudemard@lg-1r17-n04.guillimin.clumeq.ca>      2 commits

       cd0dfbcfa Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       ad3c17ea6 bug fixed : 1) unmapped read in realign 2) small script error in puure

  eric audemard <eaudemar@imac6-ub.(none)>      15 commits

       3b2258dfc bug puure
       c0acd51bb bug puure
       0a33731c7 bug puure
       c4ca9549c update base.ini for guillimin
       b47c7c9f1 add base.ini for guillimin
       ab7535acb 3nd script of PUURe done + some bug
       5e7c1fe3b Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a0b5418c8 2nd script of PUURe done + some bug
       9576d8453 add file for puure pipeline
       844021a9c bug correction on 1st  script of PUURe
       c49f5583c first script of PUURe done
       89793044c add art software (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) : Simulation Tools add perl lib path for SVDetect
       315bfc11a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_resources
       2b250c7a7 add script for install: SVMerge, RetroSeq, TIGRA, cmake
       f89034f1c add biopython

  eric.audemard@mail.mcgill.ca <eaudemar@abacus2.ferrier.genome.mcgill.ca>      4 commits

       caa90df11 add ini files. bug fixed on mkdir cov directory
       477148257 add ini files
       46f547e27 debuging step: fixe some bug
       00d5fb8a4 bug fixed on puure + add all output in dnaseq metrics (gatk, picard)

  Francois Lefebvre <lefebvr3@ip03.m>      1 commits

       e84abf37b picard version 1.125 module

  Francois Lefebvre <lefebvrf@gmail.com>      54 commits

       22fd75c50 Updated R install scripts to reflect mugqic_pipelineS new name and deprecation of mugqic_resources
       ccb88f595 fix to deploy script
       3e061c3fd update R packages list in resources
       22d054545 Fixed bug when -v is specified
       cd63425da more cleanup
       a3e5bbe28 more cleanup
       2be18dfe0 no message
       8536bf30d Fixed bug in prefix mecanism
       2e83482e1 Fixed default behaviour of update mode
       2bb28121b R_Bioconductor.sh: update mode is default
       7b04a5f2f First commit of the new R install scripts
       6e908acc9 STAR install script update to 2.4.0e, accounting for new folder structure in the archive (bin/,source/, etc)
       39922d9b1 shortStack install script
       b032d1db1 tcltk is part of base install
       c80bb97cb Updated star install dev script
       46888c058 magittr and other packages added to install list
       d8773a33e Guillimin before mammouth
       b819fb602 cufflinks version change
       a62b51950 Updated dependencies list (gqMicroarrays)
       0a5aeba11 Added sailfish module install script
       dfb486d34 bowtie 2.2.2
       97999d491 minor changes to top hat and vienna install script
       00d65aa64 Newer org.MeSH.* packages are too numerous and their installation take forever. Excluding them from the org.* packages list.
       7d51b6a92 Created gmap-gsnap install script (dev)
       af8b6c08e Added HTSFilter to del list
       21cc2016b more package dependencies
       f34a34dc0 Fixes to Trinotate related install scripts
       b937e6f73 Removed G phase1 from R deploy script
       8ae1b01d7 Multiple install scripts related to trinotate
       d3a8fe241 module install script for BEERS, RNA-seq Illlumina reads simulator
       7aca566e8 no message
       ba527a2c8 GCC module call for phase1 only: before compilation and called by R module for run time.
       915b729e6 modified example R.sh calls; need sh -l for module call to work on guillimin phase 1 ....
       2402dee2c bash synthax fix
       34d5126ef entrezdirect installation script
       2de4f3040 Added -p option to R.sh: name of an env variable that specifies a prefix to -i and -m.
       ba1638a45 no message
       9903c2a2a cron modified, next clip module
       c071c1e4c no message
       6c8c26f41 Nextclip install script
       99a908445 wrong file name in cron script
       c4ed552f2 no message
       d21308ef1 temporary commits.. sourcing from Dropbox because problems with abacus
       42a0781c7 Cleaning and and chomping module files too
       21ef6d8fc Changed chmod commands
       0e7cacb08 Problem with roxygen2 on CRAN. install.packages does not find it for R<3.0.2. Msg posted to r-help. R.sh will not roxygenize until this is fixed.
       67164c5db Last commit to R install scripts. Three files added: - R.sh is the workhorse which checks installs R/Bioc if necessary and performs updates. List of dependencies is hard-coded within this script. By default it will install/update the latest R version to $MUGQIC_INSTALL_HOME_DEV, hacking build.R to insure any subsequent package installation is performed with umask 0002. It can also install specifyc versions with -v, install to a specific location with -i and create the module file at a specific location with -m. The -f option forces rm -rf on any previous installation. The latest R version number is obtained by parsing the VERISON file from a R-latest.tar.gz download.
       3c6a6c60b tar is silent, chmod done outside R
       baacdadf0 Combined installation/update in one script. package list and old install script will be removed in later commits
       d186f78bc Added R/Bioc Update script
       453fdbd08 added varscan install script
       8a2f54322 Updated R install script
       a3a16f386 list of R packages deps now linking to mugqic_ressources
       22b024283 no message

  François Lefebvre <lefebvrf@gmail.com>      3 commits

       4ae7affbb tophat and bowtie2 according to template install script
       35e4743c8 Minor changes to deploy script
       cc01370aa no message

  Joël Fillon <joel.fillon@mcgill.ca>      363 commits

       bd175defd Updated README-release.txt with more info
       21cbf6555 Version bump to 2.0.0
       c3aea54f1 Added utils/ to PATH in module mugqic_pipelines
       ae6ae813a Updated module mugqic_R_packages to version 1.0.0
       b91ac7bc0 Standardized mugqic_pipelines module install script
       493da6efd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       6eb4df27c Updated PacBio pipeline steps as docstrings
       2ebef9da5 Updated module mugqic_tools to version 2.0.1
       718961a57 Fixed bug metrics.matrix sort tmpMatrix.txt before join
       ab9e00178 Fixed bug Picard sort_sam: add BAM index as output file only if sort_order='coordinate'
       661081442 Added trimmomatic adapters FASTA files in bfx
       975beb6ea Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01f101c05 Minor update in rnaseq de novo pipeline
       4b11101df README.md edited online with Bitbucket
       d5b5827fe README.md edited online with Bitbucket
       e897da548 Updated SnpEff to 3.6 and snpeff_genome=GRCh37.75 by default in dnaseq pipeline
       078f739f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d0e8ec030 Separate bed2interval_list job from picard.calculate_hs_metrics in dnaseq pipeline
       347f88842 Updated mugqic/tools/... to mugqic/mugqic_tools/... in config files
       5b93fc371 Updated and standardized module install mugqic_tools-2.0.0
       b7626d277 README.md edited online with Bitbucket
       fd69eeb4b Renamed mugqic_pipeline into mugqic_pipelines
       36aed5473 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       66a477b5d Added multiple input files for star step in rnaseq
       35ac8b059 Added star memory/cpu settings in rnaseq.batch.ini
       9990a91bb Minor output changes in star.py
       dc44a1dc4 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f7667c98d Added step multiple candidate input file support + updated trimmomatic and dnaseq bwa_mem_sort_sam, picard_merge_sam_files accordingly
       a6dda2616 Added BAM index file as output of Picard merge_sam_files and sort_sam
       8c3c87f20 Removed unused trimmomatic skip option
       bc4bf7393 Adjusted picard_merge_sam_files ppn value in dnaseq.guillimin.ini
       11c164cd4 Added symlink to BAM index (*.bai) in nanuq2mugqic_pipelines.py
       ae0a63f35 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6cca375a2 Fixed BED Files split bug in nanuq2mugqic_pipelines.py
       5edef0900 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       98f2dd85f Updated dnaseq dbNSFP2.0.txt to dbNSFP2.4.txt.gz
       a5926c5d8 Replaced dnaseq dbsnp and known_sites parameters by known_variants + updated genome config .ini files with dbsnp_version + minor fixes
       11f624583 Added vcf.gz.tbi in genome install
       4cf7e4352 Added log_report.pl and nanuq2mugqic_pipelines.py
       86167d7ff Reverted to previous pacbio_assembly name
       fd9fd0f5e README.md edited online with Bitbucket
       3e0de926a Added RNA-Seq report step (needs update for cuffdiff section)
       165073d55 Removed differential_expression.goseq for cuffdiff and added self.gq_seq_utils_exploratory_analysis_rnaseq jobs in RNA-Seq pipeline
       3ea16d35f Merged in jfillon/readmemd-edited-online-with-bitbucket-1417115922753 (pull request #6)
       7dbb18835 README.md edited online with Bitbucket
       8c70917dc Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a3f02c373 Recovered previous gq_seq_utils_exploratory_rnaseq modifs
       11863cb81 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       62964d7c9 Update rnaseq_denovo_assembly.guillimin.ini with generic metaq, proc, pmem cluster settings
       29a808767 Upgraded smrtanalysis module 2.3.0 with patch 1
       08552cbb6 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f604b0fc6 Recovered rnaseq.base.ini with new differential_expression and goseq settings
       6fb27160a Updated all R modules with mugqic/R_Bioconductor/3.1.2_3.0 and mugqic/mugqic_R_packages/0.1 if necessary
       4bca9119f Fixed merge conflict
       8aab691bc Added exploratory rnaseq step (partial)
       fe1bcae26 Added genes file in rnaseq config
       11016523c Updated Star module to 2.4.0f1
       f64a60c3e README.md edited online with Bitbucket
       9c0818aa7 README.md edited online with Bitbucket
       c0fc855f6 README.md edited online with Bitbucket
       45c78c092 Merged in jfillon/readmemd-edited-online-with-bitbucket-1416841132369 (pull request #5)
       a7dfa9723 README.md edited online with Bitbucket
       6583d2f04 Completed implementation of cleaning feature
       cb4c0c6a8 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       54378d51c More on pipeline cleaning feature
       8113c5176 Fixed merge conflicts + first implementation of pipeline cleaning feature
       8c3c120e6 First implementation of exploratory step in rnaseq pipeline
       4ccebd1c8 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       9c33f31a4 Added Gene Ontology files in genomes + added goseq in rnaseq pipeline
       5b2a7a97e First attempt to install genome Gene Ontologies annotations
       1a34018e6 Fixed bug forgot super() call in Illumina and PacBioAssembly pipelines
       a46d77287 Updated all config base.ini with module_mugqic_tools=mugqic/tools/1.10.5
       97243c7ac Tagged mugqic_tools 1.10.5
       10f05ae02 README.md edited online with Bitbucket
       d624f8044 README.md edited online with Bitbucket
       b52b592fb Added design file description + minor changes in README.md
       c79b9c934 Updated rnaseq_denovo_assembly differential expression cluster settings
       e8d968cfa Updated rnammer cluster resources
       b882786a1 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e703a7387 Added Python interpreter version check
       9fb047286 README.md edited online with Bitbucket
       fe5d9541f Updated ini files with RAP_ID comment and rnammer settings
       6fd646d2b Removed openmpi settings on abacus
       1486f6150 Design parsing support for empty or '0' fields
       171bf2c5e Removed command dump in .done file due to bug
       f4226eb3a Renamed readset and design file 'SampleID' column into 'Sample'
       74bef6d07 Updated cluster queue settings for rnammer_transcriptome in rnaseq_denovo_transcriptome
       a8ba35e05 Fixed version _argparser bug
       f9796a56b Dumped command in .done file + moved version in MUGQICPipeline
       a3d02bb76 Added differential_expression in RnaSeq pipeline
       38731979d Fixed syntax errors in rnaseq
       b0494b97f Added readset path normalization and variable expansion in readset file parsing
       531b02b2b Fixed PacBio readset parsing bug
       f17853b3d Added guillimin config file for rnaseq_denovo_assembly pipeline
       bf1e01f9e Moved genome_configs into resources/genomes/config; updated README.md accordingly.
       7eb3adabc Imported mugqic_resources repository as 'resources' subtree
       be12079e8 In PacBio summarizePolishing, removed redundant cmph5 sort and used a copy of cmph5 file to prevent Job from being always out of date
       a3cf8443d Updated PacBio pipeline with smrtanalysis 2.3.0
       76f9f3327 Improved Job debug log
       d658c8429 Minor cosmetics changes in README.md and PacBio readset file
       7ce1851c1 Added Bitbucket URL in command line help
       aab8a99e6 README.md edited online with Bitbucket
       0c66e6f0f Added pipelines/pacbio_assembly/README.md; removed pipelines/pacbio_assembly/README.txt
       49b0f0c5e Added
       48eb8360d README.md edited online with Bitbucket
       e0ad9310b Added pipelines/rnaseq/README.md
       329a627f9 README.md edited online with Bitbucket
       744fe945f Added pipelines/dnaseq/README.md
       e76ca730e README.md edited online with Bitbucket
       3426e21bc README.md edited online with Bitbucket
       0e4226f03 README.md edited online with Bitbucket
       24a3bf37d README.md edited online with Bitbucket
       415922ff3 README.md edited online with Bitbucket
       8bd144ec6 README.md edited online with Bitbucket
       34db02a21 README.md edited online with Bitbucket
       44024ca66 README.md edited online with Bitbucket
       a8a022798 README.md edited online with Bitbucket
       0d57ec026 README.md edited online with Bitbucket
       f6bb2f47c README.md edited online with Bitbucket
       5442dec95 README.md edited online with Bitbucket
       d738c55ec README.md edited online with Bitbucket
       73a1d8c10 README.md edited online with Bitbucket
       47ae71cf9 More on README.md
       0ae6c2d24 README.md edited online with Bitbucket
       9a9751df4 More general documentation on pipelines
       cce883328 Replaced internal RAP ID with generic  variable in dnaseq.[guillimin|mammouth].ini
       49f5e04e3 Updated igv_genome config param to match generic <genome>.fa.fai index + updated all genome configs with Ensembl 77
       bd09fa8a2 Module versions update in config files
       ace8e6f63 Fixed pacbio filtering dependency bug
       c971597db Minor comment update in homer module install
       0ba8de741 Changed DnaSeq snpeff_genome to hg19
       8abed9c1a Fixed bug readset.fastq[12] attributes should be writable
       dd7828733 Reverted SMRT Analysis module install script to version 2.3.0
       9a34f4300 Added SMRT Analysis 2.2.0.133377-patch-3 in module install script
       68c77e556 Update SMRT Analysis module instal with 2.3.0 patch 0
       b06999f6c Fixed dbSNP download_path bug
       4612110b1 Updated dbSNP to build 142 for Homo sapiens genomes GRCh37 and GRCh38
       b8e26ad0e In genomes/install_genome.sh, added functions skipping if up to date + added rRNA FASTA creation and BWA indexing if present in ncRNA fasta
       49c8cf309 Moved star_index/ into genome/
       e0fe746f6 Added install_all_genomes.sh ; fix STAR module version for index creation
       7e92a39ea Added STAR index
       d67a0c4d2 Added version number + updated README.md (incomplete)
       515ab99af Added star module install script
       39d4b84d7 Updated gqSeqUtils report calls and parameters
       162ec0bc2 Fixed merged conflicts
       93a94e937 Updated RNAseq nozzle reports with Python version
       e4bae26b1 Updated DNAseq nozzle reports with Python version
       5ad7b1648 Renamed bio module into bfx
       ad6c8b17d Added config file list support in gq_seq_utils report  and pacbio_assembly
       83cd49112 Added report and compile steps in pacBioAssembly pipeline
       59fafec44 Added log.debug for job.is_up2date + set config.filepath with config.trace.ini file
       c90118a76 Added raise NotImplementedError for PUUre and RRNATagger pipelines
       470221533 Updated config genome paths with new genome files organization
       13471eb5c Fixed bug smrtanalysis.reference_uploader transforms '-' into '_' in fasta filename
       1d948bc1e Moved /sb/programs/analyste/genomes/blast_db/README.txt on bitbucket genomes/blast.sh + move oryCun2.sh in genomes/old/
       46062c589 Updated default dnaseq genome paths + readset library, run, lane and config sequencing_center are optional for dnaseq bwa mem
       45ec20772 Fixed permissions bug for Homo sapiens genomes install
       eb7ea9cec Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d052f61c4 RELEASE -> VERSION; DBSNP from NCBI for Homo sapiens; symlink DBSNP to Ensembl VCF for other species
       f103474dc Added pacbio mummer step
       f2e84d748 Added pacbio steps up to blast
       29af2de9b Added job input/output files debug in pipeline
       3ae18d9b7 Fix cluster settings from job name instead of step name
       1da5b2888 Fixed pacbio pbutgcns cmd missing &&
       412287e93 Check job input/output files as valid paths instead of regular files
       a0727600a Update genome installs with Ensembl release 77
       a38bdfd75 Fixed rnaseq cluster_submit_cmd_suffix bug
       7e0c6e842 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       799ccb437 Fixed permissions on genome files
       7b1b697d9 Use Ensembl genome instead of 1000genome
       c7dd3496f RNASeqQC now uses GTF with transcript_id only
       7aa6712c8 Removed abacus config .ini
       b6689cc14 Updated ini files with proper genome keys
       2e9c839b0 Added genome configs
       802c39db2 More on pacBio assembly steps
       49bdc5c85 Updated rnaseq genome config names
       e91a300f9 Increased abacus cores and memory for insilico_read_normalization_readsets
       59e7accba Added DBD::SQLite in perl install dependencies
       3fae17c4b Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       ca5b12590 Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       f551b415e Fixed trimmomatic headcrop length order + insilico_read_normalization_all job name
       eeb1ff623 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       37a126061 More pacbio steps + here-document for bash command to avoid escape characters in output
       4cdd7a9ec Added DBI and PDF::API2 as CPAN dependencies in Perl install
       0c4d6cd70 Fix RSEM reference dependencies
       59f01f960 Standardized perl module install with reduced CPAN dependencies
       f51ae0d21 Minor trimmomatic thread update in rnaseq_denovo_assembly.base.ini
       c3be9ab35 Added rnaseq_denovo_assembly.mammouth.ini + minor path fix
       74e88c074 module_tools -> module_mugqic_tools + pacbio filtering step (incomplete)
       1870e1ab4 Added rat,dog,chicken genomes + python fix compiled in ucs4 + matplotlib manual install
       359971604 Standardized chipSeq modules as separate weblogo, homer, macs install scripts
       07f6d3bd6 Tagged mugqic_tools 1.10.2
       f29ef0019 Various fixes in rnaseq + pacbio first draft
       ea16a1d2a Removed picard_index and sam_index subdirectories + error tweak for gunzip human_g1k_v37.fasta.gz
       2cb5f0ba8 Added dnaseq alignment readset subdirectory + mugqic_log not sent if pipeline has no jobs
       94590ba4f Cleaned genome_configs
       17a3b31ff Cleaned genome_configs
       ec62c5817 Minor documentation change
       a62ab9260 Removed old RRNATagger files
       3d5278035 Added MUGQIC parent pipeline with call_home log function
       a409b72ac Added config trace creation from merged input config files
       178d1c12a Change readset file column name Sample -> SampleID; remove -V from qsub; updated gtf tophat index path in rnaseq
       0b1bd44fe Creation of the major model genomes with new standard organisation
       8a6a093ff Fixed puure.abacus.ini filename
       dba20ba52 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       cf3c3bac1 Fixed merge_and_call_gvcf output file name + updated dnaseq with mugqic/GenomeAnalysisTK/3.2-2
       8614cd1f5 Updated rnaseq config files and parameters + minor dnaseq config update
       ab1c6d142 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a1466e78f GENOME_INSTALL_TEMPLATE with Ensembl support
       0b2880b71 Fixed mkdir metrics in illumina.py
       66b529d2c Minor fix in dnaseq.batch.ini
       e00235131 Minor fix trinotate nb of columns + scheduler date formatting
       bff030190 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       bc2b9fc20 RNASeq De Novo pipeline implemented except nozzle report
       23a4a7225 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e6212716d Fixed samtools_sort and baf_plot options for dnaseq on mammouth
       4491132db Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       7b211be37 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       c91f8351d Added fastasplit in RNASeq De Novo pipeline
       dc5044225 Minor fix in dnaseq.mammouth.ini tmp_dir and R version
       a6570ef6a Minor scheduler formatting
       2a9c29851 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       edaf6f101 Added normalization steps in RNA-Seq De Novo pipeline + fixed python lib bug with absolute path
       f88a4483a Modified -j and -c options descriptions + removed all .pyc files
       47a1c9527 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95b82a391 Restandardized mugqic_tools install template
       2ab7d2ac4 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       10f5ecce7 Added first draft of Puure pipeline
       82006acb0 Added dnaseq mammouth config file
       e4b99eea6 Added dnaseq guillimin config file
       6f1cf2424 More on GENOME_INSTALL_TEMPLATE.sh
       5dc06903a First draft of GENOME_INSTALL_TEMPLATE.sh
       4eb982341 Minor fix in Python module install
       c950b250d Fixed python module install + minor fix in MODULE_INSTALL_TEMPLATE.sh
       ddd1e5d6a Standardized python module install
       8eede7f88 Moved AMOS module install in dev and standardized it
       e2912e038 More config standardization
       e22b5c4da Minor fix in ucsc module install
       0c9abe7f9 Standardized UCSC module installation + minor template improvements
       2873b94ad Fixed merge conflicts in bio/rrna_amplicons.py
       475e92cb4 Reorganized all config parameters
       b8d0673d3 Minor fixes
       d6094312f Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       d73fe4958 Fixed minor tophat bugs + scheduler proper exit code
       e19aaae85 Added rnaseq steps up to cufflinks
       6cf0a0653 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b65411f18 Added cufflinks module
       104113f2c Added gtf in htseq + argparser factorisation in pipeline classes
       8b68b0914 Added bio/htseq.py
       dd3896724 Fixed newline and output file bugs in R mugqicPipelineReport + batch scheduler with job.done
       96fcc8e4d Fixed merge conflict in samtools.py
       631d37f61 First htseq draft in rnaseq + snv fix in dnaseq
       655c6bb0e Added multiple ini files argument feature + first draft rnaseq raw_counts
       7f92b4905 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a351008c8 Added rnaseq wiggle step
       3d3d84408 Fixed input/output files bug in concat_jobs()
       bdf6c8bbf Removed all .pyc from git
       f6995a9b4 Added rnaseqc step in rnaseq + other minor fixes
       cec67b7a5 Grouped common dnaseq and rnaseq features in a parent pipeline Illumina
       192486b20 Fixed minor bugs
       4441fb945 First rnaseq draft; dnaseq with os.path.join
       443c9c92c Update input/output files for metrics_snv_graph and deliverable
       f776268fb Modified sys.path according to dir changes
       710c0bd67 Reorganised directories
       55e7904c3 All dnaseq steps implemented (need to test though)
       7cde825a5 Added dnaseq steps: rawmpileup, rawmpileup_cat, snp_and_indel_bcf, merge_filter_bcf, filter_nstretches
       af2c35a04 Added dnaseq steps merge_and_call_gvcf, dna_sample_metrics
       b0ea4b8cf Added dnaseq haplotype_caller step
       952ebecb8 Formatted batch scheduler + calculate_hs_metrics, callable_loci, extract_common_snp_freq, baf_plot steps
       2a3b223eb More dnaseq steps + job input files validation
       dba57675f Moved python scripts in mugqic_pipeline directory
       c941e47bf Removed step loop + standardized step names
       96e6e96b7 Added indel_realigner step + fixed job is_up2date with dependencies
       6cedd3344 Added first GATK functions + extended dnaseq pipeline
       2c2937cd4 Added first gatk draft + parse Nanuq readset file
       397fb364f Added global config object + picard module
       4dd79a1e0 Added group, pipe, concat job functions
       fef73cd38 Added DnaSeq trim and mem + force option
       24057c2ac Torque scheduler test OK
       b484d113f Added torque submit header
       dc45f445b Added scheduler and log
       b6fada28f Updated mugqic_tools install with version 1.10
       16e9aa355 Added argument parsing and more
       2e6cc0ebc First Python version with basic args parsing
       3ee3e3edc Core classes coded
       97d89408b More on python prototype
       fa8dd9ad4 Added DB_File Perl module + minor permissions fix
       a24fd3a90 Added prerequisite for modules rnammer and trinotate
       e6d3ded6f Standardized trinotate and dependencies module install
       af7cea5db Standardized BWA module install
       53b960598 Moved R&D modules in dev/
       fbccba11b Added archive check in MODULE_INSTALL_TEMPLATE.sh to avoid redundant archive download
       1e72c28d3 Minor fix on template install: archive dir with /
       8898ddf85 Restored blat module in dev
       54192b2fd Removed blat module (already part of ucsc module)
       6a70836b3 Prefix java module by openjdk
       d5575c0d4 Updated Java module install with OpenJDK
       3e9546b7b Updated RSEM and Trinity module install with latest version
       f80d0d29a First Python draft
       65176a94c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       28f5ab6fe Fixed conflicts in merge master into bam2fastq branch
       26c226acf Resolved conflicts with master branch
       1729a4a3f More on object-oriented dnaSeq pipeline
       2378320c4 More on object-oriented design
       a3dabbed2 More on object-oriented redesign
       3d942d885 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       455ebedec More on object-oriented pipelines
       20bc7af9e Minor variable name changes
       b2c5d0f61 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       d017c57de Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       9b1ae9041 Incomplete version of bam2fastq dnaseq pipeline
       e87e5b921 Added Parse::Range in Perl modules
       29bd3d93d Solved merging conflicts in rnaSeqDeNovo.mammouth.ini
       a91dbdc4c First object-oriented version of RNA-Seq De Novo Pipeline
       2255e9c9d Minor RSEM install script fix
       a7c127596 Updated RSEM install script by modifyin Perl script shebangs with /usr/bin/env/perl
       a37edd04c Updated mugqic_tools to version 1.8
       063caa27e Added PerlIO::gzip as module dependency
       2fa1750e6 Added dev/prod/ environment in MODULE_INSTALL_TEMPLATE.sh + modified all shebangs with bash
       ae1028995 Merged dnaSeq.pl conflict from master to bamToFastq
       a6d0d7619 More on bamToFastq process (uncomplete)
       f911e8a1b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       89a306620 Merged master in local branch and resolved conflict
       ced737e99 First bamToFastq attempt in dnaSeq.pl
       338afbec6 Minor fix in snpEff archive name
       8154c76e8 Standardized snpEff, VCFtools install scripts; added Java install script; added various permissions
       17a4f9533 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       c351f0995 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       b64d2e761 Merged with master
       27c29de88 Even more on bam2fastq
       9c38e4011 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       04f66e44a More on bam2fastq
       ca1172ced Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       e37f19aed Some steps toward bam2fastq support
       ab015b832 Standardized tophat and cufflinks module install script
       201cc7d70 Standardized bowtie, exonerate, samtools module install scripts
       bbc946ef1 Fixed module directory named "pipeline" instead of "mugqic_pipeline"
       25141fc19 Standardized mugqic_pipeline module install script
       732d095fa Added others permissions for MUMmer and wgs-assembler
       493b28b03 Updated mugqic_tools version number to 1.6
       689141b2a Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ea2d0aeb8 Installed trinity and rsem in production
       9b4f23fe4 Standardized Trimmomatic module install
       ff81b2729 Updated RSEM version, RSEM and Trinity permissions
       e51836d6a Upgraded to mugqic_tools-1.5
       f205a7aa4 Added BLASTDB variable in blast module install script
       b3f1686d2 Minor permission fix
       9147435c1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e894200a1 Updated Trinity module install with last version
       819415a08 Updated mugqic_tools install script according to MODULE_INSTALL_TEMPLATE.sh
       5093baf08 Renamed module install template
       9473a9040 Added read/execute permissions for others in module install template
       0e6cf37b0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       99c14da5a Added ape + ctc R dependencies for Trinity differential expression analysis
       f6648d9e9 Updated Trinity edgeR PATH
       1f140b6e3 Removed comment in Trinity install script
       0348f395c Updated trinity and rsem module install scripts in dev
       5e0865f76 Moved some install script to dev/
       91784d03b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       29666f4c4 Clean up obsolete/dev install scripts
       4ec387a7a Removed Makefile draft
       90d94d3ad Minor change
       771112388 Initial import of genome and module install scripts from mugqic_pipeline repository to mugqic_resources repository

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       9f787f315 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95ac9e380 BFXDEV-133 detected perl bug on homer and weblogo scripts using our own perl module, updated perl scripts shebangs
       4b5a245a6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       72a19c1fd BFXDEV-133 Adapt the software/module installation to guillimin phase2 requirements
       90db5b548 BFXDEV-82 version 1.7 of mugqic tools, added bug correction for R-tools related to chipSEQ pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      12 commits

       0418260b9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d33f8e249 prepend popoolation path to PERL5LIBS
       c5a79a8c1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       05d653225 added Picard Mark duplicates no optical v 1.90 to the dev software repository
       5edf0ed9f pmarquis: Added genome setup for arabidopsis thaliana-TAIR10
       4bc9c17f1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a8a9a0f7c added env-variable POPOOLATION_HOME to the module file
       13c9b20e8 added an installation script + module for popoolation
       f2046f103 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       2c85b4f1b added plotrix to the list of R packages
       eaab7c17c creating R v.3.0.2 module in dev
       3a6476df3 installing R on the test environment

  jtrembla <jtremblay514@gmail.com>      24 commits

       4d370a664 rrnatagger.py continuation
       a228ba6ea completed bio/rrna_amplicons.py
       a4bbfc718 continued coding of rrna_amplicons.py.
       f82ac5228 writing rrnatagger pipeline in python. First step works.
       cf4051430 Added two functions (that actually works :-)). I need to convert the rest to python syntax.
       bc82663d9 added rrna_amplicons
       d045a6a45 renamed RRNATagger to rrnatagger
       d0a7f3418 Added rrnatagger.py and list of steps.
       9f616eed8 added lib path.
       b0d9af62c Added bamtools.
       907ef6beb Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       87e713ab9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e3ddbbb10 minor fixes.
       058fa827f Updated bwa version.
       6432b27d5 New module install scripts.
       cecd97f0d Added path (prepend-env) to root/R-tools.
       c61c904f4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       119c7add0 minor modif to module list for perl.
       2fe26541d Added comments to smrtanalysis installation (2.2.0). Tested and installed.
       9858f66a2 Added muscle aligner.
       9403ab75e Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources corrected typo, then merges..
       ac2a43e43 corrected typo.
       5c1e23a04 Added chmod at the end of the install module "scripts".
       79e2d4a65 Added list/script to install modules. Could be improved...

  Julien Tremblay <jtrembla@ip03.m>      1 commits

       72347dfc1 added module DB_File to perl module installation.

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       7d2d5a3af Added RRNATagger-tools to prepend path in the module file. BFXDEV-31
       b792e6afc added python install script.
       1be6f845b Fixed two missing paths for PERL5LIB
       0b26b2d17 Added perl to mugqic modules.
       82b9d1974 Modified install scripts according to our group's template.
       3156d84b5 Added gnuplot-4.6.4 and memtime-1.3

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      1 commits

       923d66198 Removed genelenght.r. Use mugqic_tools/R-tools or gqSeqUtils directly with HERE docs or Rscript -e

  lefebvrf <lefebvrf@gmail.com>      14 commits

       0be2b7674 Removed superfluous (already in base) dependencies from list
       4cf1272b4 violin plots package added to list of dependencees
       16882a5ec PerlIO::gzip is a dependency for Trinity clusterProfiler as R package
       0530c576b Guillimin phase 2 added to daily R deploy
       9a37c3ce2 R.sh needed a module call to a more recent gcc than system one.
       e90b9cad0 devtools::install_local() would install packages to home folder if exists...   .libPaths(.Library) solves the problem
       923b50ad9 —vanilla removed from R.sh too
       a13f58208 blurb added related to last commit
       d028e16a7 Add creation of Rprofile.site to R install script. This will force using cairo X11 backend since cairo is not always set to default when available…
       321311bc2 Added ViennaRNA and mirdeep2 install scripts
       33dfb8787 Small hack to tools::build.R when installing R allows umask 002 (!)
       74967799f Corrected exit code, first version of R_deploy
       b8db61a31 Additional setdif(reps) to avoid duplicate installations -> slow
       2082c05c9 Fixed leading tabs problem in module file by using <<-EOF here tag. Neutralized R_LIBS to insure installation proceeds correctly when R module already loaded Added vanille biocLite() call

  lletourn <louis.letourneau@mail.mcgill.ca>      45 commits

       e76f595ba Merged old perl changed into python
       0f7f0f873 Fixed haplotypeCaller output file name extension
       4098ccc3d Version bump mugqic pipeline to 1.4
       b9625abd7 Version bump to 1.5-beta
       7e8d28c9f Bumped bowtie to 1.1.1
       2528f9276 BVATools version bump
       277abf2a1 FastQC version 0.11.2
       c1d9878cd Version bumped Ray to 2.3.1, removed unneccessary fix now
       45b85da11 PIcard Version bump to 1.118
       9db0c5481 Added pysam
       fc0ea4c9c BWA version bump AND fixed the script, it was a bad merged script
       cc5aa56e1 Updated GATK and mutect
       2cd15d63c Merged diffs
       b10f91edb Added gsalib and added configure switch
       cb1777f71 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       367791619 Version bump pipeline to 1.3
       58fff5f1b Force output file name
       ef35a74e0 Version bump bvatools to 1.3
       a441cda36 Version bump of mugqic_tools and snpeff
       9b6fa74b8 Added ascat
       a66943937 BWA Version bump to 0.7.8
       610973909 Added the jellyfish tool
       610892242 Updated gatk
       09e52708c Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       44e530a1b Added Mutect home
       9a96cbfd5 Version bump
       6827d55fe Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       adda76530 Version bumped wgs-assembler, bwa, picard
       674305d44 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       9e6a4c645 Version bump bowtie2 bwa picard
       287aa2344 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c2aac9ab3 Added aspera instructions
       701df8ce1 Version bumped vcftools to 0.1.11
       6f5a6ed59 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c8f949ace Version bumped blast to 2.2.29
       d0dab7f7d Bumped version of picard
       4e6fcc484 Updated snpeff
       ff4f3974d Version Bump of the pipeline
       34c15b08f Added variable to access the repo's location
       fa87e1974 new Ray version
       a6402ac6b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ac6513f52 Updated BVATools to 1.1
       ea6aff9a1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       dca240f8c Added BVATools
       71bfcf243 Added wgs assemble

  Louis Letourneau <louis.letourneau@mail.mcgill.ca>      1 commits

       e1c5341f7 BFXDEV-246 Version bump

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       14060c198 Remove conflicts in modules/mugqic_tools.sh and modules/dev/star_dev.sh
       cb81397fb update modules/mugqic_tools.sh to 1.10.4
       b2db257fa remove decrepated python module script and add a new one
       16d306848 replace dev module in module/dev/ cufflinks_dev.sh  star_dev.sh
       fd6119e6d Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       85a452337 add genomes/oryCun2.sh  modules/cufflinks_dev.sh  modules/star_dev.sh
       b35f6decb up-date mugqic_tools.sh version to 1.7

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      47 commits

       0eff10b25 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6272cb0fb Python - RNAseq : allows more flexibity in fdr and p-value for goseq AND remove native goseq approach - BFXDEV-289
       e78717a64 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8694d0c5b Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       218dbdcac Python - RNAseq: update rnaseq.py
       a1333101a Python RNAseq - updates
       6a061db50 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       5e8df5440 Python RNAseq - updates
       61a12cf0e Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b1f9bc8c2 change RNAseq python update
       5c9e82f88 change RNAseq python cuffnorm output to include sample names
       015fab6f6 remove RNAseq python step bug exploratory v2
       582e5e3c5 remove RNAseq python step bug exploratory
       b14fbf9f6 remove RNAseq python code conflict
       62a0e4f07 PYTHON -RNASEQ: star, cufflinks, htseq updates
       40970ffb6 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e600a91d7 Python: RNA - update star + DGE anbd goseq input oputpuyt file correction
       01fccf96d resolved pipelines/rnaseq/rnaseq.py confilcts
       87f92a365 major RNA implementation: STAR, picard, cufflinks, etc...
       a39a35e8c Python RNAseq: add cufflinks > 2.2.0 new workflow - cullinks - cuffmerge step done
       d902b4d57 Python RNAseq : update STAR - add optional read sorting during alignment
       8f572bc70 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8305fe708 Python - RNAseq: add STAR alignment 2 pass && add utils folder/function for methods generic unrelated to any software or tools
       c713736cc Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       4cd30e38e update samtoiols module install to 1.0
       860eba6ba python : RNAseq - change ini
       a6ea523b8 On goinig adding star to RNAseq
       33f80d671 add dnacopy in R.sh edited online with Bitbucket
       3bdf7cb51 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       0f1910ead correct python.sh script that was not workin for the modules: BFXDEV-46
       27348a08f Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       793e633d9 add specific RNAseq packages to R package list
       1016ccd9e replace correct permissions in mugqic_tools - BFXDEV-55
       fda280248 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       40d5ef584 ngoing python module package
       03dab0981 update breakdancer module
       7f60581a2 update breakdancer module
       494455b53 pdate breakdancer module
       9deedd931 change in samtools modules and add breakdancer module
       02fde1498 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       bfdd54443 add module igv
       60e67320e Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       6615a3b5c create module install sh for the mugqic pipeline
       30983a974 update nugqic_tools to the new version 1.2
       834b58e83 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       19f75087a update mugqic_tools
       5e75b3d85 update modules/mugqic_tools to point to the new repository mugqic_tools with tag v1.0

  mmichaud <marc.michaud@mail.mcgill.ca>      54 commits

       4617cd8ed Run processing: Remove duplicate '.' in dup bam file
       f72f529a5 Run processing: Run blast on bam file, fallback on fastqs
       ae2fcbd19 Run processing: Add the suffix '.sorted' in the readset bam name
       7bb7a5066 Run processing: Merge md5 jobs (fastq+bam) and run qc graphs on BAM (or fallback to fastq)
       5603e8884 Run processing: Run DepthOfCoverage even if there is no BED file
       a81e7199f Run processing: Use bcl2fastq 'force' option instead of removing destination folder
       1890439c9 Run processing: Remove dependencies dependencies from copy job
       6f33eaae8 Run Processing: Use a kind of factory to manage the different aligners
       e006fa914 Run processing: Remove the dependency to samtools by using picard to generate the STAR BAM index.
       f46022cc7 Run processing: Move md5 step later, to optimize job allocation and to minimize condition when the md5 job is finished before the copy job is submitted
       8d0697349 Run processing: Remove unused adaptor columns parsing
       fd4705a5a Run processing: More documentation (sample sheets)
       60abc0cac Run processing: Output the index metrics file in the output folder, not the run folder
       88a090738 Run processing: Use xfer queue for downloading BED files
       ba1b4d6ba Run processing: Add samtools module version and increase number of core for STAR (a test job took 71G of ram)
       318612135 Run processing: Manually generate BAM index file when using STAR.
       1b42a02ea Run processing: Fix usage in readme
       f4e87a75d Run processing: Add documentation
       a6e183838 Run processing: Add missing '.' in bam name
       7aad80b65 Run processing: Fix copy step exclusion
       7100d7c65 Run processing: Copy output file is now in the copy destination folder
       073b8a97c Run processing: Change back to manual copy job depedencies gathering
       ff84ee02a Run processing: Initial RNA-seq support with STAR
       b5969a303 Run processing: Return empty list when there are no input for the copy job
       b5b30a6d8 Run processing: Fix copy job inputs
       58e432618 Run processing: Fix end copy notification job queue
       43c153abb Run processing: Fix qc output file
       9d9e87acd Run Processing: Fix race condition in HsMetrics interval file creation
       7dc2b9b6d Run processing: Fix configuration for blast and qc
       edb326945 Run processing: Generate sample sheet in the submit job method
       f5f23ff0d Run processing: Fix config for index step, now using 'index' category
       348972bc0 Run processing: BAM metrics are run in markdup output
       176e96f74 Run processing: Various fixes for the first test run
       63944b971 Run processing: Various fixes for the first test run
       1043bb1ce Run processing: Various fixes for the first test run
       c8ee9b0d7 Run processing: Add basic support for different aligners
       aa16b2dde Run processing: Get copy step dependencies by introspection
       c74a06b59 Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       3af310b2b Run processing: Delete existing 'Unaligned' folder when running fastqs
       583c80f4a Put Illumina configure bcl2fastq in a job
       aa93387f6 Run processing: Use . as job name separator
       ddafcfe87 Use $MUGQIC_INSTALL_HOME variable in config file to replace hardcoded paths
       79d09cce8 Fix arguments for IlluminaRunProcessing by removing readset argument from MUGQICPipeline
       cbe8ddc21 Run processing: Use run_dir instead of run_directory to follow output_dir convention.
       aa213fe57 Run processing: Improvement to the sample sheet generator to loop only one time
       7d16bcf2d Run processing: Seperate config files for MiSeq and HiSeq. The base config file still hold for a HiSeq
       df783bc2f Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       2b53d1100 Run processing: Add wget download of sample sheets and bed files. Fix Casava sample sheet generation.
       07aacb68f Add some imports from __future__ to ease the transition to python 3
       7af4b28f8 Add copy job dependencies
       82841bf45 BFXDEV-283 Fix reference for coverage calculation, add copy step.
       213e7160c Change illumina_run_processing name to follow code conventions
       8d8e1af2c Code convention changes
       f51f5e001 First draft of the illumina run processing pipeline in python.

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.ferrier.genome.mcgill.ca>      1 commits

       acde23713 adding sh tools for riboRNA

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      1 commits

       f6ddfbebd BFXDEV-112 install genomes from IGENOMES

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       2fc3f9d24 usearch.sh

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      2 commits

       481bc98a0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       7ff159c46 gene_length

1.4        Mon Nov 17 13:15:48 2014 -0500        139 commits

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       667c6f7d4 README.md edited online with Bitbucket
       ee9efd526 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvr3@ip03.m>      6 commits

       3355905ee Previous minor bug fix (dots in sample names) actually introduced major bug.
       0bc3e58d6 temp files remove
       93c2905aa rnaseqc other options possible
       c04554162 rnaseq .ini file more tophat options. BFXDEV-215
       8270318a6 Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.
       0ef6afc99 Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.

  Francois Lefebvre <lefebvrf@gmail.com>      7 commits

       242f2bc93 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c004e5163 Passing projectName to SubmitTocluster can create invalid job names. Replaced with string
       e578fb415 variable name was inappropriate
       eab031001 Updated rnaSeq mammouth template .ini file.
       202215e2f Defined a job name prefix for wig zip.  'metrics'  as a job name was not enough information
       59f152825 cuffRescolumns and dgeRescolumns now in goseq param section. Also adjusted those values for UCSC genomes hg19 and mm10  templates, original value did not work.
       524dd82c9 overwrite.sheets=TRUE to avoid common problem of updated projects

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      4 commits

       257c44d0e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4b8a23be9 patch3 applied; see BFXDEV-260
       3c6e85023 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ff5fe3dea change R module due to crashing nozzle report generation, see BFXDEV-255

  Joël Fillon <joel.fillon@mcgill.ca>      34 commits

       0860442e1 Updated pacBio .ini config file with new smrtanalysis module name: 2.2.0.133377-patch-3
       c2d500077 Updated config genome paths according to new genome organization
       030e601c0 Updated chipSeq mammouth config with new Homer 4.7
       be23877a7 Added explicit Python module in rnaseq cuffcompare step
       9d16e5a01 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c2569139d Added newline after mugqicLog command
       8f3c5f357 Updated rnaSeq.mammouth_hg19.ini with generic /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod path
       1fd8a6369 Removed explicit RAP ID in RNASeq De Novo config files
       0efe22d79 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6aac202f2 Removed moduleVersion.htseq ini ini config files
       ded1a4601 Fixed missing java module in igvtools
       2c37b3ed6 Minor fix in dnaSeq step_range help
       cf0861199 Added default RAP ID in RNA-Seq De Novo guillimin config file
       e67468afe Added --exclude <fastq> option in illuminaRunProcessing rsync for samples having BAM files
       98420e402 Uncommented phone home in dnaseq
       0e6e8a8d1 BFXDEV-203 Create one .done with checksum instead of one per output file + update config files with default adapters path + update perlpods removing -e option
       988606153 Fixed bug missing Library Source -> column not mandatory + updated pipelines/rnaseq/rnaSeq.guillimin.ini with accurate module versions
       e07ffd80f README.md edited online with Bitbucket
       d02ae67fc Added comment to update Resource Allocation Project ID
       13bf5c975 BFXDEV-221 Migrated abacus/guillimin config files from msub to qsub
       d8a411e47 MUGQIC call home is now run inside bash script after job submissions instead of bash creation.
       5c2627ae3 Updated chipseq mammouth config with python 2.7.6
       29163dd56 README.md edited online with Bitbucket
       bad364846 Added call home feature notice in README.md
       94130140f Fixed RNA-seq de novo resume jobs; added env variable WORK_DIR in pipeline; formatted bash output
       b2f3fd52b Another glob fix for prefix path check in LoadConfig::getParam
       65b992322 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       768184e84 Minor fix for config param prefix path
       74d43282f Added prefixpath check type in LoadConfig::getParam
       249a5cb16 Fixed transdecoder bug ln: creating symbolic link 'Trinity.fasta.transdecoder.pfam.dat.domtbl': File exists
       59ea4b6fb BFXDEV-32 Fixed wrong transdecoder file path for missing PFAM 'cds.' prefix
       9c6fe6f39 BFXDEV-32 Fixed pfam missingcds. ID prefix + blast clusterCPU tag for guillimin and abacus
       73ce18015 Added chipSeq pipeline cleaning
       9f6e203ff Fixed rnammer missing modules hmmer 2.3.2 and trinity

  jtrembla <jtremblay514@gmail.com>      14 commits

       6c0720baa --arg for low abundant clusters. after 99% clustering. BFXDEV-31
       d09db8983 Added low abundance cluster cutoff argument. (After 99% ID clustering step). BFXDEV-31
       9166850c1 Put more lenient parameters for itags QC to make it more 'universal' for most projects, especially those for which quality of reads 2 is low. BFXDEV-31
       05e7b99c2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f60fa07ab Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a996d3f19 Indentation correction. BFXDEV-31
       38d8da843 Fixed parsing for blastdb command. BFXDEV-30
       1b20686f7 Updated README for 454 data processing instructions. BFXDEV-31
       8bafb19a7 Added even more description. BFXDEV-31
       c003951a7 Added details to description of output. BFXDEV-31
       fcd21e5d3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9dcc9bc7 Added step to filter number of blast results. BFXDEV-30
       5d0c803cd Removed --vanilla. BFXDEV-30
       309ef5e24 replaced -num_target_seqs with -num_alignments. BFXDEV-30

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      1 commits

       a16e5e9e5 Added --sampleSheet argument to getMiSeqBarcodes.pl BFXDEV-31

  lefebvrf <lefebvrf@gmail.com>      2 commits

       b7f7d6614 Changed TopHatBowtie.pm to put an end to  the .fa.fa.fasta.fa symlinks madness when installing genomes. Parameter is now the bowtie index basename, consistent with the tool's documentation.
       d29a172f8 Fixed dots in sample names bug BFXDEV-51

  lletourn <louis.letourneau@mail.mcgill.ca>      31 commits

       3f785de99 Version bump to 1.4
       38d0a8f4f BFXDEV-39 Fixed realigner when only one job is given
       f728ac9c7 Updated bvatools
       5e168e6de Tweak guillimin parameters
       ea83ecd53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1748e1a1e Updated parameters
       69b8f2ee8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5da052afb Adjusted cluster requirements
       4cc7adfcb Removed useless param
       4dace3ec8 BFXDEV-256 Added step range to paired variants
       cdbd4a08d BFXDEV-256 Added step range to paired variants
       b53557036 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46d016cc0 BFXDEV-254 Fixed GATK 3.2 change on CatVariants
       10fb56b31 BFXDEV-252 Removed flagstat
       125945f18 Fixed BAQ from pileup and ram from fixmate
       71366af8e Changed picard to version 1.118 to fix the freeze when an exception occurs in asyncWriter
       e9112eb7d BFXDEV-216 Removed per lane metrics
       056878870 BFXDEV-245 Fixed uninitialized error when no bams are present
       682a99703 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6d2b309bd Fixed CCDS location
       85578b84b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       54fde94b2 Merged master
       407b9fee2 BFXDEV-236 optimized settings and split human builds
       419e0a0eb Added missing perl module
       976ba1f9a Fixed missing HS metric, add 2 cores to bwa
       9cf6dc5f3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c3fab1c89 Missing validation silent
       5cfb61f2e Add genome versions of ini files
       de4a383e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       dd39c54ac Fixed typo in mergeAndCallGVCF section
       33669702b Version bump to 1.4-beta

  Marc Michaud <marc.michaud@mail.mcgill.ca>      1 commits

       6b8f9f7bb Add missing use

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       3c0fb6e66 update ini to fit the new mugqic_tools tag 1.10.4 - BFXDEV-275
       e5bd47669 correct dnaseq bwa samnse dependency bug - BFXDEV-253
       6fdf80268 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f3abd8cca correct cuffdiff input double array issue when checking the job object is up to date

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      13 commits

       3bd1f63e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0e991b0b2 Update pairedVariant
       469415ff5 RNAseq replace headcrop at the good position in the trimmmomatic command; remove by default headcrop from the ini files - BFXDEV-267
       534615dc1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       402db54a9 PairedfVariant.pl: change where pindel get the insert size info - BFXDEV-266
       87a9da680 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       192600ec9 PairedVariant: allow mutec to run without cosmic file - BFXDEV-263
       62e335d3a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9696dbc6c update dnaseq and paired variant ini files
       37ad08b7a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       21e3c79e5 remove conflicys in pipelines/dnaseq/pairedVariants.abacus.ini
       a156ae54a pairedVariant.pl formAt SV and CNV to new standard - part of BFXDEV-41
       cb5b91e22 cuffdiff now should not be relaunch in a restart if it exit correctly during the previous analysis BFXDEV-212

  mmichaud <marc.michaud@mail.mcgill.ca>      14 commits

       b9ef53296 Fix genome path
       ec3104553 BFXDEV-269 Change genome files hierarchy
       0234cd15e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       5eb75aad0 Fix usage to reflect new -s argument
       f16ff095e Allow more ram for DepthOfCoverage
       403317089 Fix CPU limit usage error by using less thread for the GC
       88e50961a More RAM for QC (to avoid java heap space errors) + More threads and more walltime for mem (to avoid wall-time exceeded errors)
       d946aa77c Run processing: Add option to force download of sample sheets
       7feb92307 Fix quote escaping in filter
       88832015c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       0e60d3dab BFXDEV-171 Trim index in the generated sample sheet according to the first/last index specified as parameter
       c5a35cbaf BFXDEV-210 Use Parse::Range for steps to run, as all other pipelines
       598ff88ed Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8e8ca8ed6 BFXDEV-211 Don't send email when jobs are successfully completed

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      4 commits

       0cb571a44 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f934df3ec MiSeq ini
       1f34db5dd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31ae4509b add ini MiSeq

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      1 commits

       446b4d809 ppn=16 pour mem

  Pascale Marquis <pmarquis@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       90e7a6188 rm illuminaRunProcessingMiSeq_PM.ini

1.3        Mon Jun 2 10:02:07 2014 -0400        109 commits

  Joël Fillon <joel.fillon@mcgill.ca>      17 commits

       993125831 Partial MUGQIC remote log for RRNATagger
       da59f18d5 Remote MUGQIC Log Report in chipSeq, dnaSeq, pacBioAssembly, rnaSeq, rnaSeqDeNovoAssembly
       521765dff Updated RNA-Seq De Novo config files with latest module R/3.1.0
       f59efb988 Updated rnaseq_denovo default config ini files with trinotate steps
       019458453 More Trinotate steps in RNASeq De Novo
       9a5cd2d01 RNASeq De Novo config files conflict solved + openjdk
       1cfd28807 Added file cleaning for RNA-Seq De Novo pipeline
       4382d804a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       109247835 Solved conflicts when merging master
       42e7aca14 Updated RNASeq De Novo pipeline with Trinity version 20140413p1
       240d576e5 Beginning cleaning
       e1c902100 Cleaning of cleaning...
       ad6ed7abc Fix on samToFastq/trimming dependencies in chipSeq pipeline
       02d3cbe70 Added BAM file check in samToFastq
       5e1aeb21a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq_basic
       8ea8c98fe Basic bamToFastq support for all pieplines
       9afcd9e05 Version bump

  jtrembla <jtremblay514@gmail.com>      32 commits

       6e3d8edff Added ini file for production purposes (QC assembly among others). BFXDEV-30
       462ecaafe fixed end step for BB. BFXDEV-31
       06a442386 Added BigBrother modifs to RRNATagger pipelines. BFXDEV-31
       e0c59aca2 Added sample counting in pipeline loop. BFXDEV-31
       ebd465f66 Added a description of output files. BFXDEV-31
       e0272b5cd Fixes to itags_QC. decision if primers are present or not. BFXDEV-31
       b199a3d85 mooooooore fixes. BFXDEV-31
       d67905fe7 more fixes to ini. BFXDEV-31
       7b8a64ce8 updated ini files. BFXDEV-31
       609c2043f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       befd271ea Changed arg to adapt from percent to X hard cutoff.
       16248cda4 Replaced percent cutoff by X cov cutoff. BFXDEV-30
       d79cc3885 added / updated ini files.
       65ed0f35c ini files of RRNATagger changes. BFXDEV-31 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9e30f84e6 updated ini files for RRNATagger. BFXDEV-31
       37c1b0c59 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4cb24e16e ini file modif. BFXDEV-30
       de66a9b9c update celera specs for hgap3. BFXDEV-30
       0faaea634 Updated file check for restart at blasr ste. BFXDEV-30
       83c0cbf10 Modified way seeds.m4.filtered is handled. BFXDEV-30
       00482bd54 Updated ini file for pacbio assembly on guillimin. BFXDEV-30
       9174116c3 fixed lib path. BFXDEV-30
       f823b9370 Fixed tab indentations. BFXDEV-30
       5f8123da3 Upgrade pipeline from HGAP2 to HGAP3. BFXDEV-30
       98d20c460 Removed unused SMRTpipe .xml files. Only keep filtering xml file. BFXDEV-30
       468aacb5f Updated mugqic tools module.BFXDEV-30
       a1ae2534c Fixed output .mugqic.done file for pacbio dcmegablast. BFXDEV-30
       e7426a89a Added module load perl in subroutines. BFXDEV-31
       056cae4ab added perl in ini files. BFXDEV-30
       cf512b28f pacbio ini file for abacus. BFXDEV-31
       85ebc8aef Added blast parameters to ini files for guillimin. BFXDEV-31
       52a9f7fdc Added blast step for pacbio rRNA tags data type. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       cc1fab05f updated ini file abacus.
       9363d1fa2 Updated ini file for abacus. BFXDEV-30
       c3a118919 replaced compareSequences with pbalign for numthreads. BFXDEV-30
       7189270ee ini file for hgap3 on abacus. BFXDEV-30
       37df4a842 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1bf18fd29 Fixes to HGAP3. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      38 commits

       8f385581c Version bump to 1.3
       c051f9aaf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       42735dd4d BFXDEV-207 Fixed when interval lists are created
       102ee1ded Fixed gvcf bugs
       613a990b9 Added emtpy quotes to empty keys so they don't become ARRAYs
       b594056eb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d39b5d97d Updated BVATools version to 1.3
       9eae448d3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9cc183277 BFXDEV-204 removed varfilter
       82592ff0c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8b4f2bdc8 BFXDEV-196 Added onTarget metric
       2f57d4fd4 Completed POD documentation
       8176ab0db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1d6793d6c BFXDEV-161 Added last steps to Haplotyper caller
       240ca527b BFXDEV-198 Added simpleChrName in the ini since it was removed from the pm
       7f0d5398e Fixed params
       39a7b383d BFXDEV-196 Added onTarget metric
       8e7a04e7f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       582546148 Ram was too close to max
       db33ed9da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89c52165b BFXDEV-193 Use new fixmate from bvatools. Fixed bad step dependency
       ff9d44c7c BFXDEV-192 Added possibility to have extra flags in indel realigner
       f30f24e5a BFXDEV-182 CollectMetrics sometimes needs to GC so 2 cores are needed
       20768a455 Merged changes
       de76e02a9 BFXDEV-181 Updated java vm version
       d9468715b BFXDEV-153 Added a way to ignore readset status
       4ad3df161 BFXDEV-176 create MD5 on alignement
       52b96be80 BFXDEV-174 Test that the data is valid
       25f2a408f Merge branch 'master' into haplotypeCaller
       d4a588dd1 BFXDEV-168 added depth ratio plots BFXDEV-161 added haplotype caller
       6596bf35b BFXDEV-166 fixed no index read, but one is associated
       7b4bdc2bb BFXDEV-162 fixed hiseq recognition
       3ad331949 BFXDEV-157 Add callable region generation stats and BED in DNASeq
       9f0e145fe Added BVATools depth of coverage as a transition phase
       17c7a8fbe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c775e183b Don't use BAQ since we use recalibration, this was decided awhile ago
       1c6b7a4f6 BFXDEV-156 fixed csv encoding issues
       46573b706 Updated settings for phase1 + phase2 merge of hardware

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       35060f638 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7faafe059 Create pipeline cleaning lib; partially implemented with RNA cleaning sub only BFXDEV-74
       201124c0e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       51777364c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       27c987aea change how exit status is catched and exit the correct status in case of pipe discrepency BFXDEV-140

  mmichaud <marc.michaud@mail.mcgill.ca>      11 commits

       8c059d4ee BFXDEV-155 Add rat and mouse alignment
       4e9b1e166 BFXDEV-185 Separate config file (HiSeq, MiSeq)
       637526f92 BFXDEV-164 Fetch sample sheets when they aren't specified and they aren't on disk
       49e0b8c95 BFXDEV-183 Download each bed file only once
       f87afda6d BFXDEV-184 Don't rsync Images folder. Was used on miSeq for debuging purpose
       ae9782f8c BFXDEV-186 Don't use the BAM generated by markdup
       c5d36af92 BFXDEV-180 Use processingSheetId as dependency id instead of sample name which isn't unique
       6f3d923b8 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       ec7897689 BFXDEV-171 Add option to specify first and last index
       e3ad96737 Add option to specify first and last index
       d8eacb492 Run Processing: Add dependency on the metrics in the copy step, in case the markdup & BAM MD5 was already done in a previous pipeline, but not other metrics

1.2        Fri Mar 28 16:11:44 2014 -0400        200 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       c18dfb65c RSEM more cpus on mammouth
       9fff58f40 rnaseqdeno mammouth ini tweaks for trimming and RSEM

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.(none)>      3 commits

       00721219a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8bfdcffc5 module load python in rnaseq.guillimin.ini --for htseq-count
       bb66dccea rnaSeq.guillimin.ini changed default tmpDir  --BFXDEV-144

  Joël Fillon <joel.fillon@mcgill.ca>      46 commits

       2ca8f7374 Version bump
       0ce7d1af9 Solved conflict for merge master and bam2fastq branches
       8696c9bd6 Changed gzip to zip to compress rnaseq de novo outputs
       d86c3022c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c09c29247 Changed archive command from gzip to basic zip
       63656f179 Minor fix: removed semicolumns in chipseq default config file
       f383320ab Removed unused variables and functions in all pipelines
       3924dbe6d Added Version number in RNA-Seq De Novo Usage
       3662f2e90 Minor mugqic/tools version update in rnaseq de novo ini files
       9c9b747e4 Minor ini param adjustment + comment fix
       be76b2624 Removed hard-coded path to modules.sh: not required when invoked with Bash shell
       56881db46 Updated default project paths with new /gs/ partition
       4c24b2b6e Removed deprecated GetFastaAlias.pm
       2d92ee878 Replaced shebang #!/usr/bin/perl with #!/usr/bin/env perl in all Perl scripts
       6b7b3760e Fixed bug SequenceDictionaryParser filepath with environment variables + set param [annotateDbNSFP] dbNSFP not mandatory
       842ca5fc7 Removed redundant file existence test in SequenceDictionaryParser.pm
       851df470c Check all config modules only once when config hash is built, to reduce runtime
       b5f3c9f1a Minor fix for adapters path in RNASeq De Novo config files
       a86bd7986 Updated adapters paths in RNASeq De Novo ,ini files
       9d73b18f6 Added  [Warning/Error] in msg display
       d11cb76dd Removed old lib/LoadModules.pm
       1126cd98a More minor bug fixes for getParam validation
       670c892d9 Fixed minor getParam validation bugs
       e002a7a4c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       872d712d1 Added validation in all getParam functions
       d9621c0fa Fixed merge conflict in lib/Picard.pm
       88411ce31 Added moduleLoad validation + major style reformatting
       c130dc99e Added && between job commands to get right exit status
       080410a7a Added param definition validation and module availability in LoadConfig
       2ee8ea04e Added raw read dir and adapter file path validation in Trimmomatic lib
       a4615118f And more and more about parallel normalization
       c32e3162b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into sample_norm
       8b316d0c2 Fixed metrics parent step
       f6e131b0d Even more parallel normalization
       da10fb173 More parallel normalization
       5b5522cdb Normalization parallelized by sample
       7ef8faa95 First draft of resume-jobs
       39c49a5db Use module mugqic_dev/R until properly deployed
       813d19793 Fixed typo
       3d51ae699 Fixed normalization stats file name
       b667299f0 First stable RNA-Seq De Novo pipeline version
       de50a894b Updated blast results and DGE matrices file names
       04e1d5bc5 Merged conflicted rnaSeqDeNovoAssembly.pl
       d8717e48f Added POD documentation + fixed bug blast-longest-transcript file
       cdc15b6de Added BLAST for longest transcript only, with results header
       df03b3b4a Added metrics and deliverables steps in RNA-Seq De Novo pipeline

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      1 commits

       6a2eca819 README.md edited online with Bitbucket

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       90ac5a428 BFXDEV-133 incompatibility between /usr/bin/perl and mugqic/perl/5.18.2. Added perl HOMER_HOME/bin/ to the program execution in peak annotations and motifs usign Homer
       5e62ee014 BFXDEV-133 detected bug dependencies between trimming and alignment chipseq pipeline
       972c98b8b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0fe7a1c94 BFXDEV-133 updated software, parameters, corrected bugs in chipseq pipeline for guillimin phase2
       2401cf3a9 BFXDEV-133 adapt chipseq pipeline and configuration file to guillimin phase2

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      17 commits

       7f020642e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ba4cfb0a2 bug in ini files: the following variables are not defined : genomeSize, annotation distances, markdup, report variables
       942ef2a2e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       050cc3cc0 BFXDEV-123 add design example files to chipseq and rnaseq pipelines, update the user manual
       599e532ac BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       4006f00d0 BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       8290e24fc BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper - typo
       b79e9f7ee BFXDEV-36 Generated PODs documentation for rnaSeq.pl wrapper
       c1873cecc BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper
       5bfe14efa changed my email by johanna.sandoval@mail.mcgill.ca in standard ini files
       9d2232b2d BFXDEV-85 added flagstats/ generated a file with number of reads after filtering/mark as duplicates
       1d069aa21 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6680dcf39 BFXDEV-84 correct bug for restart when merge lanes step failed
       394b3aebb BFX-1799 wrong variable initialization for genomeSize, detected when genome is not human or mouse
       5ac34dc47 comment skip trimming step from standard ini files, added imagemagick to mammouth ini
       e2b57e1ec Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8c0e000b5 BFXDEV-77 bug in merge metrics, instructions to run Chipseq in the pipeline directory

  jtrembla <jtremblay514@gmail.com>      15 commits

       8d91c821a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline Merge latest change prior to modifs to PyNAST (unset LD_LIB).
       05f53351f Added unset LD_LIBRARY_PATH to PyNAST step. BFXDEV-31
       8b7990406 Updates to pacbio stats step. BFXDEV-30
       956d5e37a Fixed pdf report for nc1. BFXDEV-31
       c578e8380 updated ini files for RRNATagger. BFXDEV-31
       abb4b9846 Added module load openmpi to PyNAST. BFXDEV-31
       0eeea2ea4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cd179fcec Added nozzle report for RRNATagger_diversity. BFXDEV-31
       22b7fea91 Removed Iterator::Fastx from the wrapper as it is not even used here. BFXDEV-31
       884627fcb fixed file checking for restart mechanism for sub loadPulses. BFXDEV-30
       d275e03f2 fixed file to check for input in sub referenceUploader. BFXDEV-30
       cab615309 Minor fixes to restart mechanisms. Removed inputFofn for filtering step and loadPulses step. BFXDEV-30
       e9c1d8191 Corrected some parameters for celera assembly step. now on lm2 by default. BFXDEV-31
       752dd0ad1 Implement module load and getParam checks. BFXDEV-31
       8afc93fe1 added missing path for abacus ini files. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      33 commits

       3cf8824c8 fixed path to primer files. BFXDEV-31
       60e8f9234 fwd and rev primers options now optional. BFXDEV-31
       7f869723e Added modules for R and mugqic_tools to rarefactionPlots.R . BFXDEV-31
       add6ebd17 mugqic.done fixes to rarefaction subroutines. Added mugqic tools module to appropriate subroutines. BFXDEV-31
       72338c790 Updated help screen and removed appending ./scripts/ to PATH in curr dir. BFXDEV-31
       6a982c815 Removed scripts/ dir and moved it to mugqic_tools
       c99e920bd minor fix for tmpdir definition. BFXDEV-31
       8e3af423d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f806248ef minor fix to ini file (for mummer). BFXDEV-30
       dabb656a5 Modifications to README and added mapping file example. BFXDEV-30
       d3473a44f README for RRNATagger. BFXDEV-31
       82ad83466 Added RRNATagger (16S/18S/ITS rRNA amplicon) pipeline.  BFXDEV-31
       a9eeb683e fixes to ini files (pacbio pipeline). BFXDEV-30
       a6bc0cc27 put only 1 mersize (14) in the gullimin ini files. BFXDEV-30
       7bcfc6575 Fixed bug with fofns. Minor modif to main loop. BFXDEV-30
       126618dc5 Changed /dev/null in the order of commands so no empty consensus.fasta anymore. BFXDEV-30
       e9d2882d7 Forgot a && before gunzip in variantCaller  step! . BFXDEV-30
       907be7d2d Fix for uncompressed consensus.fasta.gz. BFXDEV-30
       d5e5250ee Fixes for compatibility with latest version of pacbio assembly pipeline. BFXDEV-30
       fea30cc69 Fix for .bas.h5 files. BFXDEV-30
       935aea088 Updates to PacBio assembly pipeline which now support multiple polishing rounds.
       35d5972d9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       28c7b8298 Compute assembly stats from polished assembly and not unpolised ones. BFXDEV-30
       ebdbae001 Fixed and updated pacbio ini file for abacus. BFXDEV-30
       920b75f45 changed callVariants for variantCaller . BFXDEV-30
       faf6c43de typo corr. BFXDEV-30
       0b000318b Added additional instructions in the README. and changes relative lib path. BFXDEV-30
       468184fae Forgot to update this library for pacBio pipeline. BFXDEV-30
       eda60c32d updated README. BFXDEV-30
       6217fda09 corrected for typo gnuplot . BFXDEV-30
       89f4c68fb corrected relative path of lib folder
       fdf930539 BFXDEV-30 modified parameters so they are more generic.
       27ab9137e Loop/dependency fixes to PacBio pipeline. Added compile stats step at the end of pipeline.

  lefebvrf <lefebvrf@gmail.com>      2 commits

       5e46e9028 necessary to honour cairo as X11 backend for R graphics with current module installation
       8707d7825 vanilla will hinder reading Rprofile.site, which in  our modules will not be used to force cairo as X11 backend when available

  lletourn <louis.letourneau@mail.mcgill.ca>      27 commits

       8c107ea1a BFXDEV-149 Fixed the way BVATools is called for coverage.
       b839f5aa1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c7f0a75eb Removed path test since they contain env variables
       31322fb89 BFXDEV-124 fixed center when using mem
       c2625ccbf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e420186c7 BFXDEV-116 Fixed reverse adapter
       4f5f9d210 BFXDEV-111 Fixed many module versions
       973c7b029 BFXDEV-114 Added R loading
       4be23c334 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       419f6fec4 BFXDEV-109 Configured java version
       b9989e54d Updated guillimin pipeline
       70343460f Updated depth of coverage ini
       db2b3a836 Fixed case when there are no BED files
       b79a3a0f8 Fixed module typo
       be34921f5 tmp hack so nanuq takes coverage graphs
       8cf40e60b BFXDEV-89 BFXDEV-88 Change GATK to BVATools for DepthOfCoverage and support multi bed in project
       7f45e96e4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b0a885e0c Added the missing mappability key
       5b568bab2 Added line to keep overlapping pairs
       73d736b86 Added umask to dnaSeq jobs
       7e1f584a7 Updated module versions
       67b9a7866 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       472a3cd84 BFXDEV-73 Fixed undef jobId if step is up 2 date
       a6e3c7c95 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2d9d66a85 Added details to the release guide
       70f59da57 Version bump 1.2-beta
       daabc8db2 Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       3cd0dee69 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22a87a0f5 DNAseq: correct SNVmetrics dependency BFXDEV-136; RNAseq: remove hstseq dummy/null module usage BFXDEV-137
       cbdca4194 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       355837088 DNASEQ: correct depthofCoverage missing bed file field in sample sheet BFXDEV-135
       38db1aa46 ask bash to source /etc/profile.d/modules.sh BFXDEV-125
       184fd53d6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3cbe61a62 add the -T option to cuffcompare BFXDEV-93 and validate previous commit for BFXDEV-47
       a671e06e3 README.md edited online with Bitbucket
       f9a76160d README.md edited online with Bitbucket
       8f494a855 README.md add the bioinnfo email adress
       6cb9d24b1 correct the 1st line typo in chipSeq.pl
       f1df7bab0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31f3438ee replace ; by && in pipelines except for the rm of .done
       8620232b0 replace a default ppn=20 by ppn=12 otherwise the job will never be launched

  mmichaud <marc.michaud@mail.mcgill.ca>      33 commits

       e242aa193 BFXDEV-147 Fix index mask of a single-index lane in a dual-index run
       3c527b26f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       bce066621 Use new version of BVATools with simpleChrName support
       0cf39804c RunProcessing: gentle perlcritic compliance
       6e3b82e43 Fix BED file list from SampleSheet
       4753244b3 Fix again the readsqc BVATools path
       becd68c0c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c8fc20f3d MD5 job doesn't need 5 cores, only one is ok
       7b03a8f78 Fix readsqc using current session bvatools jar file instead of the one loaded on the compute node
       230a617f4 Fix run processing when there is no target file
       cb1be91f6 Merge branch 'master' into runProcessing
       59829b2a9 Use new module loader for BVATools qc
       61c25ba0c Add cranR module version
       7e38002fd Merge master into runProcessing branch
       1261ab2ab Support spaces in bed files
       2b5b54e27 Don't align mice yet
       9bb906f79 Support spaces in bed files. Fix bvatool module version
       3e2571b0e Fix usage message (to show optional parameters) and print error message on dying
       f549388a1 Die when there is a barcode collision. Fix bwa rgCenter for mem
       172452bb4 Add a parameter for the number of mismatches
       da5e47a10 Default values for both sample sheets path
       214c3d929 Merge branch 'master' into runProcessing
       9584091d9 Update BWA module version
       1a6989d59 Fix rsync command: quote character were not escaped correctly
       a8e58fdc9 Add missing depedencies for the copy job (metrics)
       f19991cb3 Don't create qc folder when the qc job will not run
       7385c8e2d Depth of coverage: add reference parameter
       7508bd41e Change '_' to '.' as seperators in the metrics job name
       2ea610366 Enforce processingSheetId column in the sample sheet only when processing a run
       f79fd8eee BFXDEV-76 IlluminaRunProcessing Add target coverage metrics & change job ids to support multiple sample with the same name in the same lane
       8e90e6234 Merging upstream changes of the barcode metrics jobs. Less core and memory used, skip last read
       5a2b132dc Merging master into runProcessing branch
       70383b361 BFXDEV-76 Illumina Run Processing Pipeline

  Pascale Marquis <pmarquis@lg-1r17-n02.guillimin.clumeq.ca>      2 commits

       21722231d rnaSeq.guillimin_hg19.ini
       d75c72d74 rnaSeq.guillimin_mm10.ini

1.1        Mon Dec 23 14:23:34 2013 -0500        137 commits

  Joël Fillon <joel.fillon@mcgill.ca>      69 commits

       305139289 Commented code
       78e706c63 Check rawReadDir value in configuration file
       3e137a01b Minor R version fix
       2eec3e95b Added trimMetrics step, fixed blast outfmt single quote bug
       758a68b55 Added option file check, module availability check
       7fabc9bc0 Added genes/isoforms length values in edgeR results
       faefc9333 Added blast description in edgeR results
       c54b33674 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e61eb0f30 Cleanup of rnaseq_denovo files
       02ed5e69f Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       378d4d977 Updated ocnfig .ini files with trim parameters
       2da3da4a0 Fixed missing escape $R_TOOLS
       470df00c1 Fixed differential expression bug
       36711d427 Added differentialExpression step
       467076bfd Added multi-step dependency support
       eba140593 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       003c41b2c Fixed Trinity.pm merge conflicts
       6046b1b95 Added trim step + various fixes
       efc7b744e Minor fix
       5c140d0cc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e31fb3ce8 Added guillimin config file + minor modifs
       faf0d714f Updated .ini config files with cluster-dependant processing values
       19f6b6d4f Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       dabafcfc6 Added blastCPUperJob config tag
       0eb0a4a76 Added blast step
       01b495ff2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       885d5f10e Updated default rnaseq tmpdir on guillimin: /sb/scratch/jfillon
       776831664 Added first version of BLAST step
       03e5423e4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       0a05c4e8a Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       fcb3c5dec Minor fix
       dc1e136ee Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       a80a6c285 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       29f179f43 Updated config ini files
       63a2a5d2f Minor fix bowtie CPU
       b3cf0fc13 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       27fd66f70 Massive cleanup
       f2eeb87e8 Another minor fix in mammouth ini file
       60fac504d Minor fix mammouth ini file
       2822aab11 Updated mammouth ini file
       a7f477838 Added TrinityQC step + code reorganization
       340a63d03 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       9d7d75ec1 Moved normalization parameters in config file
       950b713ca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       3fd953143 Fixed bug edgeR semicolumn
       f850faf61 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       c9e35cf83 Added edgeR step
       42d400d26 Added rnaseq_denovo.mamouth.ini
       a20391f8c Fixed bug unsorted list of fastq.gz from find
       4a35ad657 Extended wall time for normalization and trinity
       73bc9e6b6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       db70b578f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       02f65b8e3 rsem prepare reference separately + job log relative path
       ef38d5f8c First version of RNA-Seq de novo pipeline with normalization, trinity, rsem
       39a95e7ea Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       74d4ac04a Further development of RNA-Seq de novo assembly
       92aa2512a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e43a3e2db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       8946fc183 Further development of rnaseq_denovo pipeline
       06fb9bcb2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       d2a0332b8 Import of old deNovoAssembly in rnaseq_denovo branch
       fc8082b0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       052b80411 Minor fix
       eddaed251 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       ab7f6212d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       9eb13833b First draft of de novo RNA-Seq normalization
       f4c998af8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       587329772 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       23f34b257 Fixed bug missing '\' before

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      12 commits

       29ff2382c Readme for chipseq pipeline
       6f007a96f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       062d6986d document mugqic pipelines setup using md - correcting bugs
       30364fc8c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       fad747c27 document mugqic pipelines setup using md - correcting bugs
       5a0a009bd document mugqic pipelines setup using md - correcting bugs
       f771f1afc document mugqic pipelines setup using md - correcting bugs
       848249a17 document mugqic pipelines setup using md - correcting bugs
       58dca89d4 document mugqic pipelines setup using md - correcting bugs
       96522902e document mugqic pipelines setup using md - correcting bugs
       4cf001dcb document mugqic pipelines setup using md
       4e9a82a71 Miscelaneous graphs chipseq pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      21 commits

       a5f8c7053 corrected bug in qcTagsjobid
       d3b16fdd2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ecd9434cb Added PODs documentation to chipseq pipeline wrapper
       3555ab408 correcting links to project's directories on README.md
       1b94392d9 link to wiki on general README file
       31db440f9 BFXDEV-63 Added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       589689e04 added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       eb5a068ca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       4e1e7df2c added variables for chipseq peaks annotation plots
       ea839b33b corrected bug loading ReadMetrics and read metrics sample/run separator
       9298aae6b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       0d69c5297 commit chipseq report branch
       95aab1a2a added annotation statistics + graphs to the chipseq pipeline
       94bcb8ba6 adapt chipseq pipeline to the new resume jobs functionality
       73ee7a7ae adapt chipseq pipeline to the new resume jobs functionality- config file
       f1dabddaf adapt chipseq pipeline to the new resume jobs functionality
       4fe99531b Adapted chipseq pipeline to resume jobs changes
       9d6f3d886 Transfer chipseq pipeline metrics libraries to the pipeline space
       df472f0a2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       2e5f6134e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       67ae37980 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report

  Johanna Sandoval <sandoval@ip03.m>      2 commits

       2a5fbe0c0 document mugqic pipelines setup using md
       6e68220f2 adding mammouth configuration file for chipseq pipeline

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      8 commits

       1971ef140 Added memtime to dcmegablast and blastdbcmd. BFXDEV-30
       9f0b66a27 Fixed unwanted mofifs to BLAST.pm. BFXDEV-30
       b430c4cea Modifications to perl packages for Pacbiopipeline. BFXDEV-30
       bb5a19ea0 Major modifications to the pacbio pipeline. Functional version tested on abacus and guillimin. BFXDEV-30
       a17cd5ac7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       97cd35a69 Added readme for pacbio assembly. BFXDEV-30
       6432aada7 Remove my info fields in ini file. BFXDEV-30
       67c11c1cb Added PacBio assembly pipeline libs/wrapper/files, etc. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      19 commits

       4a6d8af2e Version bump 1.1
       ef0d74373 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2541c188c BFXDEV-68 Added mutect to paired pipeline
       d3998b623 Removed gpfs for guillimin
       ed3de9b15 Fixed conflicts
       11d0fe512 Fixed undef on steps that are up2date
       bc778813d BFXDEV-67 use a true csv parser
       bdbb7cee4 Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline
       4d7c7fe0a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1a8cf5944 BFXDEV-45 Fixed single read handling, fixed undefs
       5d2800084 BFXDEV-15 Changed the name of lane level bam metrics
       ed46494c8 Fixed starting dependency and final job output
       4555ab025 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       92520500d Added pre dependency and final sample output
       aa452ff80 Changed BWA version for run processing and added some metrics
       a29bac1e7 Added runprocessing parameters
       ab6593218 BFXDEV-52 DNAseq realignment restart generates duplicate unmapped bam
       05272a81c BFXDEV-48 added missing close BFXDEV-49 added support for alignment+QC only RNASeqQC
       e1ecac6c5 Version bump

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       da9919d3c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fa47c1305 replace ; by && - BFXDEV-47 and update the rnaseq ini files
       e0ccc024e dnaSeq.pl correct typo (realign step) : BFXDEV-54
       371c03f48 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       799fca20f change/test replacing ; by && un command line
       0a07178ef Trimmomatic.pm remove single trimm bug

1.0        Fri Nov 8 15:03:24 2013 -0500        794 commits

  David  Morais <dmorais@ccs.usherbrooke.ca>      2 commits

       262012018 Merged in daveM_denovo (pull request #2)
       0764a2564 Merged in daveM_denovo (pull request #1)

  David Morais <moraisd2@ip03.m>      94 commits

       4c13dd19d remove entry
       f458e9897 added new entries
       d8e671b1d split butterfly command lines in chunks
       a963e53f9 creates BLAST best hits reports
       61d9bf3ba Fixed bug
       5b9c16052 Modify how the butterfly files are splitted
       0234e161b fixed bug
       5b2bbbc20 modify command line
       f29fe2e70 fixed bug
       f2bf4ab3c Added full path to groupInfo left and right samples
       e8794a0bd Added option for de Novo Assembly
       3490e4332 Implemented and tested Merge Trimming Remove duplicates de Novo assembly Blast BWA
       8094aca94 fix paths
       191bdcc8d fixed picard version and env variable
       17934d1e9 change module add trinity
       891bcf9d0 change module add blast
       60bfa9f3c change module add blast
       78f4ba6e6 change module add trinity
       59905eb9d fixed typo
       aacc3ec5b fixed typo
       e78b5b691 Modified HtseqCount::matrixMake call. Now group is the last parameter so it can be left out.
       c13eb9135 Fix $group problem. Now if group is not specified the variable $group takes $sampleName value.
       dcabd4bf1 Fixed blast db error name
       592659c35 perldoc
       a94ec1d5e adding scripts needed by de Novo Assembly
       a64a2cee2 add comments
       289ec2f7d De novo Assembly config file sample. It needs to be modified according to your each project.
       ba64039bf mammouth_deNovoAssembly.ini
       befb42449 This is the main de novo RNA assembly pipeline. This is the first commit.
       774254dc3 modified config variable
       5be92fd59 tidying up
       7d190c5f7 tidying up
       8d9704b93 tidying up
       03fc63e7e tidying up
       9b75e54f0 tidying up
       24e31b0a4 tidying up
       2cb0b70ee tidying up
       94a0d3289 tidying up
       6a8f8fcfa tidying up
       24aa7afab tidying up
       42608a83b tidying up
       2d3f8a9de tidying up
       23e5b5e18 Coded full command line. Added sub contigStats
       d0f849d99 Modified $laneDirectory
       88bb2c754 Modified $laneDirectory
       db5f50632 Modified $laneDirectory in the singleCommand
       bc774be30 Modified $laneDirectory to assembly/ and added $group variable to deal with groups on deNovoTranscriptome assembly.
       55eb46d0b Modified perldoc
       0ad8b58a7 Modified $laneDirectory
       c21cab505 Modified $laneDirectory from alignment to assembly
       2037ad69e Modified $laneDirectory
       6bb39648f Modified $laneDirectory
       229bba496 Assingned reads/ to $defaultDir
       f0b1c6199 Assigned  $laneDirectory to 'reads/'
       eb949206d Modified $laneDirectory = reads/
       921adb058 new line
       b187d45f4 Output bast results to /db directory
       7c1754a0d change from own branch
       c58b61f3e First commit. Library to create differencial expression analysis.
       4086a94d8 Modified Doc
       cc24e6b25 First Commit. Library to generate basic statistics on raw read count. Not Yet Ready for use.
       1df2ca704 Library to create generate basic statistics of each sample. First commit. Not ready for use yet.
       edf8193c3 Added the possibility of looping through more tha one DB. In This case the DB must be passed as an argument to the function.
       3dac1f13a added comment
       fccb0b113 removed quality offset from file name
       db155a074 removed quality offset from file name
       301427d03 added indexing function
       e373c99d3 add function
       34fd4a44a add function
       986c00cf9 add function
       7183abe71 modify hash groupInfo
       a77802bef fixed input values
       6bae877d7 BLAST lib, first commit
       78008bd6c modify looping by group option
       998b38a6d tiding up
       6e243f552 fist commit Trinity
       fb5182cf0 change split setting
       719e16bc0 added new entries to sampleInfo hash
       4bb98a426 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       644271f73  a simple multiformat file Spliter. First commit
       8d3ccb39a Remove duplicate reads. First commit
       cb8cfdf9a modify to work with nanuq sample sheet
       56d16c84b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       80dd2fd8f First commit
       7fbacf8d8 first attempt on STD
       736365b06 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2c33bfc19 Reads the libraries form script_dir/lib
       c3d669c1b added function that allows the script to read libs from the script dir
       d9a378dbe add header
       e32573a82 remove package folder
       754f187a8 first commit. Lib to read module list
       16372cb25 first Commit. lib to read config files
       10e24cd06 change repo name
       c56a0af6e adding space.

  eric audemard <eaudemar@imac6-ub.(none)>      5 commits

       f5360846f add software install script : igvtools virusFinder SVDetect add genome install script : virus
       2bd3239b8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f55f728b1 add install canFam3 and new bwa
       7bdcac749 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       270bac922 adding deNovoSV pipeline from the mugqic_eaudem

  Francois Lefebvre <flefebvr@abacus2.(none)>      4 commits

       6f6477c5e febvre <ddcccccZZZZlefebvr@abacus2.(none)>
       a4f60317e Integrated exploratory analysis step
       1ea865bba Changed Rscript -e template
       a3468dcbd Fixed perl formatGtfCufflinks.pl call

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       d312e26c9 fixed problem when R packages list has duplicates
       3fae3d0a7 Merged in module_install_scripts (pull request #3)

  Francois Lefebvre <lefebvrf@gmail.com>      53 commits

       2ec42e90d duplicate section names
       2a175b004 duplicate sampleOutputRoot
       e9e50a1c6 duplciate param in ini
       2abd12c54 duplicate sortSam sections  in ini file removed
       e8f87929d fixed bwa module name dnaseq.abacus.ini
       5b3c1a9ba kmergenie install script
       810d08101 updated rnaseq.abacus.ini  wiggle chrsize path
       8321d0f5c Updated R install script
       d3087f740 R install script now installing Vennerable, gqSeqUtils, gqUtils from remotes.
       9caf20221 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       b6591f3fe Exploratory analysis step, still to test
       ac2361f97 Added RSEM install script  + updated Trinity one
       789cb0f49 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       6e37c4218 no message
       b014a95f9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       7d5e0aaf3 Changed bwa to mugqic/bwa/0.7.5a, problem between RNA-seQC and 0.6.2
       165e72d6c Updates R/Bioc script, now also installing all org packages
       aa3dcc8cf Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       2e48fa55d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       92e196aa5 no message
       5d60cfa76 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8d4ef871e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       51243690a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       34000aded hg19 install script" add gtf mask for chrM, tRNA, rRNA for cufflinks
       77f6acd42 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8165a0360 Walltimes for htseq, cufflinks, cuffdiff too small  for typical data
       d99c87762 htseq walltime
       03aa1e156 Added 72h wall time align step
       f70ac2f70 Updated bwa install script
       555609163 Added blast module install script
       a4cd05678 Added few more R packages for install list
       baaacb22f Removed MANPATH prepend from mugqic/R defnition
       8eb8a256f chmod in in mugqic_tools.sh install script.
       42d4ef93b <tools,perl-tools,R-tools,java-tools> from the svn now in tool_shed/
       9b773f430 R packages lists updated
       0725e01b8 Corrected bash_profiel for mammoth and changed R install script to 3.0.0
       83cc87ff8 Added R package "outliers" to R packages dependencies master list
       5c1d4a6a0 Updated mammouth $MUGQIC_INSTALL_HOME value in guessing script.
       bea778907 Few change to abacus wall time rnaseq
       0e30bd777 Previous fix to Metrics.pm did not work
       30f0bfb91 minor fix to Metrics.pm (will be depr. anyway at some point). chrome size file for wiggle left to default in default abacus .ini
       16132edc3 Added -T sort unix sort in metrics.pm to correct problem on guillimin. Added job name suffixes in rnaseq.pl to make job names unique on Guillimin (dependencies)
       60c84ccc6 Merge remote-tracking branch 'origin/rnaseq_test'
       ca2a4ac03 ..
       59e116ebf Added hg19 installation
       c83481d18 Added mm10 installation
       c0dabe431 no message
       c99255ecb updated to top hat 2.08 (bug in 2.07) + more genome scripts
       57248b56a Added Eric's abyss install script
       b20104a65 drafted genome installation script
       e82f783e8 Corrected htseq module by adding module load python. Also started genome setup scripts
       a4c423240 Finished python/numpy script. This script will not be 100% portable, need to set locate BLAS and LAPACK
       6da084c0c Added modules/ directory to hold modules related content

  Joel Fillon <fillon@ip03.m>      4 commits

       c47405835 Added java module in BWA lib + dos2unix rnaSeq.guillimin.ini
       5bf363d1e Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       929fb10a6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       5f3e86b1e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output

  Joel Fillon <jfillon@abacus1.(none)>      4 commits

       a051fa558 Added readRestartFile  function
       a7a14042c Print out MUGQIC command exit status.
       f55bec042 Missing ";"
       a7a2f6d41 Missing ";"

  Joel Fillon <jfillon@abacus2.(none)>      1 commits

       16641e85d Minor misspelling

  Joël Fillon <joel.fillon@free.fr>      27 commits

       67995c4ea Added simple Perl script tool_shed/getLogTextReport.pl to create log reports
       f7ddba65e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       3d61e719e Simplified getLogTextReport parameters
       e12fecac7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       aac70e278 Print out full job log path + 2 exit code outputs
       d9e2d1aa5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       75de93109 Fixed lib/SubmitToCluster.pm conflicts
       253bdfbec Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       e7d552e3d Added number of jobs in log output
       020e4e942 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       8ac84841e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       4991e383b Added timestamp in job log filenames .o
       ae0b820f6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       b42133bb8 Reorganisation of job logs in an array
       2db5f5726 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       5d39c044a Added a few comments
       bc9df34b8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       f3998beee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       30cf022c3 Fixed exit status + walltime
       39079a594 Added getLogTextReport function
       3c10f407e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       964e92d2e Fixed package name
       4c6b06b24 Merged master into logReport
       5e38b4300 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       70eff6ed2 Added readRestartFile log function + Exit Status code in .o
       d7f0f1750 Fixed deliverable typo
       5fd50047c added python tools path in mugqic_tools.sh

  Joel Fillon <joel.fillon@mcgill.ca>      8 commits

       81a3ec98b Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       e19c7bd6f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2128dd21d Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2235fc5b2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       dab8f2943 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       4cf3157df Check well-defined variables
       7b0370501 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c4eaaeb63 dos2unix rnaSeq.guillimin.ini

  Joël Fillon <joel.fillon@mcgill.ca>      95 commits

       a3fddd2ba Removed Log lib (now merged in mugqic_tools/getLogReport.pl)
       32f25ea9c A bit of a cleanup
       8dbabd9f2 Synchronized rnaseq .ini files
       c8fe2ba03 Fixed bug set @INC with relative path to mugqic pipeline lib
       06ca9c5df Moved module and genome files to mugqic_resources repository
       3f1150842 Fixed bug missing '\' before
       b1635464c Standardization of rnaseq .ini file for the 3 clusters
       6e7375fe1 Removed prereq in module install script template
       d76adbf2b Removed prereq + fix R module load to compile kmergenie
       784d0c7df Updated Picard module install script now based on template
       750cacaa6 Fixed inline comment bug in module install script template + fixed permissions in ea-utils module install script
       6fc665cfe Fixed ea-utils compilation error on guillimin
       24ab3feab Minor comment changes
       dd1f73036 Added ea-utils module install script + minor README change
       a8533c2ac Renamed MODULE_INSTALL_SCRIPT_TEMPLATE.sh + deleted old archive + minor fix
       e6a1552db Added module file permissions
       7e85ddc4c Created a module install script template; modified tabix install script
       1b2f9bc1c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       40b98e7fc Deleted redundant modules/add_me_to_bashprofile.sh; added README.txt instead
       e00d04081 dos2unix all pipeline .ini files
       14673934e Updated chipSeq pipeline .ini files with latest module versions
       013305367 Updated rnaseq .ini files with latest module versions
       ce59f4f2d Merged
       e187b2db2 Minor change in R module install script
       88ebc06e9 Changed permissions in R module install script
       d5cfef179 Added vcftools module install script
       8a5b3698a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a788a15b9 Added archive storage + permissions + cleanup in snpEff module install script
       127d3277e Removed unnecessary lib BiocInstaller for R install
       f1bc1a0c4 Reorder chipseq modules in .ini files
       ecb3c6b23 Removed duplicated modules in chipSeq .ini files
       4280f36c8 Added archive storage + permissions + cleanup in Trinity module install script
       1976a5185 Added archive storage in UCSC module install script + update dnaseq ini file with bwa/0.6.2-ptx renaming
       dfc61b083 Added archive storage in Trimmomatic module install script
       64e9a661c Added archive storage in tophat module install script
       04ba4a94a Added zip archive storage in picard module install script
       d1c2d54b2 Added archive storage msg in GATK module install script
       eed094d3c Added archive local storage in igvtools module install script
       8bdeb5bc7 Removed bedGraphToBigWig module install script (now part of ucsc module)
       ff546a71d Added permissions in Tophat module install script
       9c25d167e Updated dnaseq/validation.pl and dnaseq/pairedVariants.pl to comply with new version of SubmitToCluster
       cbaf966d5 Added permissions in Picard module install script
       7a726e60a Added permissions and cleanup in IGVTools module install script + minor changes
       37095823e Added permissions and cleanup in hmmer module install script
       34df850b1 Added permissions and cleanup in gemLibrary module install script + minor aesthetic changes
       89ac0f939 Added permissions + minor fix in gatk module install script
       ed0aecac8 Added permissions and cleanup in fastx module install script
       d7dbb8186 Added permissions and cleanup in fastqc module install script
       bd890969d Added permissions and cleanup in exonerate module install script
       ff6137953 Fixed permission bug in UCSC module install script
       568ef4eab Fixed another permission bug in module install scripts
       daa7cf203 Fixed bedtools install script bug for archives with different naming system
       c3a3067a4 Added permissions and cleanup in cufflinks install script
       1cb9cc5c4 Added permissions and cleanup in bwa install script
       3ef584d44 Updated default pipeline .ini files with new mugqic/ucsc/20130924 for bedGraphToBigWig
       942bc47ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0eac0a137 Removed -j option in make; added kentUtils.Documentation.txt in bin/
       e2253bc93 New install script for UCSC genome browser 'kent' bioinformatic utilities
       7e1b1650d Reorganized install scripts for a5, bedtools, blast, blat, bowtie, bowtie2
       d097e82f4 Reorganized BEDtools install script
       99f688cff Removed last workdir parameter in chipSeq, rnaSeq; removed commented code; fixed blast+ bug in blast install module
       a985b2433 Added shebang in pipeline bash
       27c480f10 Removed leading ':' in chipSeq JOB_IDS lists
       669a97d39 Clarified bash by using job variables + header
       6e72b6784 Removed sampleOutputRoot tag in all .ini
       ce93e2576 Removed 'mkdir output_jobs' commented lines
       89d60fd83 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       68e948a63 Removed output_jobs subdirs in chipSeq pipeline + added job dependencies in log
       c47fd5d12 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       454e4aba5 Fancy install script from Johanna
       e9f61b190 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       054fbfb37 Removed final ';' in MACS2 cmd + aesthetic change
       6e049a9b3 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8471e3857 Removed ';' at the end of Homer cmds
       ae79f343f Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       66a919323 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       24459567a Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       310035f18 formatting changes + chipSeq redirect job output to specific location
       830d57b82  ->
       19f116ec2 Minor formatting changes
       890c91575 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       995a8b416 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       a05a09015 Added java module in BWA lib
       622217d96 Reorganized job_output for DNAseq pipeline
       405043d93 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8f6b4e59e Added tophat 2.0.9
       6ef6c3430 Minor formatting changes
       e01c18200 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       9e220c5c3 Minor job parsing change
       4af616a23 Added Job ID number in log
       7ab13210b Redirect job list output into a specific file
       7438d6416 Job outputs go into jobs_output directory
       967b8889a Minor coding standardization changes
       472a9b233 Minor error msg edit
       abf1067fa Moved getLogTextReport.pl in tool_shed/perl-tools/

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      49 commits

       8e63db4b4 changed config files to fit the new SubmitToCluster:initSubmit structure
       98b7290a5 changed config files to fit the new SubmitToCluster:initSubmit structure
       f2c2aa077 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b56095e13 adapters small rna from Trimmomatic-0.22
       61c25df30 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       85c99b150 added module and installation script for gemLibrary
       3998d3d1a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       15e6a92c9 added chipseq.guillimin.ini parameters for hg19
       152e3a764 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7559eafb7 changed deseq.R : added flag comment.char=" " when reading edger_results.csv
       b641ee284 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e46c9c1e6 update install genomes script in tool_shed directory
       984993255 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5bde35c0 prepend perl-tools directory to PERL5LIB
       e8194696a adapted ini file to changes in trimming parameters (headcrop)
       af929cef6 Changed readstats module, added ReadMetrics.pm library
       ef32d6610 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       db61b3988 Parsers or trimmomatic and flagstats reads statisttics, create a readstats.csv file
       a8d891759 remove unused function (old tag directory generator)
       dfc0b2e65 compute read stats using the ReadMetrics library
       7faa22527 prepend perl-tools to the PERL5LIB variable
       b4c6417dc corrected freeBayes install on guillimin
       79ee5c6a6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       e3e854946 Adapt install freebayes for Guillimin
       f46883d6c Generate QC stats using R
       3fa1222af Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       09f7ec241 tag directory generation removed tmp sam file creation - added samtools
       c3ca36b03 tag directory generation removed tmp sam file
       d045ac122 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8288285f4 Freebayes install and create module
       237ed033f filter aligned unique reads using samtools
       176e710e4 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       df89fcb00 added qc plots, corrected motif weblogo bug: seqlogo available on v 2.8.2
       51077d33b generate chIPSeq QC metrics R script
       6087d8239 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a99e13938 Increase walltime for qctags generation
       8e47a3902 Initialize BWA_JOB_IDS per sample, create design 0 directories for MACS and annotations
       143492d52 generate wiggle tracks: remove -i parameter
       9666a94ef Validate differences between design file and nanuq project file
       c606dbd9d Validate missing values in design file
       0f2715e4c adapt read statistics to changes in Metrics library
       494d67e3b pbam flag for paired end reads on MACs
       d0cfe12c7 chipseq pipeline configuration file for abacus
       ed0bf0b76 chipseq pipeline wrapper
       2665867b7 chipseq pipeline HOMER tags, annotations and motifs
       dde6603e8 chipseq pipeline MACS2 peak calls
       69e8682fe uncommited changes, will control variable names from my scripts
       aa273e996 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       70afa15f5 transform jobIdPrefix to avoid illegal bash variable names

  Julien Tremblay <jtrembla@abacus2.(none)>      2 commits

       fde3058a8 modifications to skip unnecessery bam merges.
       cbdfb711f Added dnaclust install trace

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      3 commits

       4edd2fd45 Fixed bowtietophat module, added cpu param to align, drafted template guillimin
       136c03ffa test
       d75425a54 Added Config-Simple and guillimin ini wc

  lefebvrf <lefebvrf@gmail.com>      5 commits

       c7b8ab6db Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       63e674f16 quiet option  is hard coded, makes it hard to diagnose cufflinks… -q could be enabled in previously added option parameter
       bfc8e6595 Cufflinks otherOoption parameter, to allow for instance for -M rRNA, MT, etc mask
       62b192856 abacus htseq wall time increased from 3 to 12h
       4ed3cbaa1 unused raw_count folder was created

  lletourn <louis.letourneau@mail.mcgill.ca>      174 commits

       b07568dfd Added restart implementation
       6633b73b7 Fixed pairedVariants with new structure
       50e3e85f4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       29fc49a0e Added release instructions
       58333308b Removed deprecated code.\nAdded a version
       6294737ac Fixed step order bug
       3a504533d Fixed region bug from previous commit
       490b11e03 Added a way to skip trimming. Added a mem example
       a385b07a0 Added a way to skip trimming
       c60ac775d Merged master
       a535288d2 Added lib barcode to lane data
       441b96be3 Added missing ;
       8fa3ded4d Updated old code
       399a2d960 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9e7894fc Fixed mem bugs
       1c2a330ce Fixed badly named CFG key
       df22d10a2 Fixed paths for install. Crucial for these to work since they hardcode the paths
       fbe736a7e Merged Master, some fixes and new dnaSeq report
       9bcc238d7 Updated dbsnp
       623de3ef7 Amos install
       0ffb61ebe Fixed region for torque
       abd094086 Amos install
       0370f0efc Added split for stats
       508f1910b Fixed if test
       2202926a1 Fixed dependant files
       57f85c2c6 Fixed some bugs found while testing public datasets
       e9d5b4eee Fixed dir search order
       dab8a7604 Put appending and cleaning of dones in SubmitToCluster. More central, easier to maintain
       2cd0d6b4a Fixed to support any year 2XXX
       b9aafe9b6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1eb4c1379 Merged master with outputDir fixes
       ac706266c Merged HEAD, removal of job output dir
       a4c216255 Change delimiter for regions, msub doesn't like colon
       d1c131e22 Don't generate the pileup, it takes way too much space
       9805b6400 Added hg1k install and GNU parallel
       2be9b17f2 Forced default picard buffer to the abacus optimal one
       6230215c3 Finished updated the rnaSeq.pl pipeline
       431b192f1 Merged master, added java module, dos2unix rnaSeq.ini
       b2898cb0f Merged Master
       cf0fca306 Finished implementing reset in dnaSeq
       a2fc84a58 merged
       fcd41d227 added missing calls for output dir
       708959e7e added sampleOutput param to inis
       ef0f696cd Use the reight blat tools
       5fdb7f3ab installed v4 preview
       833686518 Added beagle 3.3.2
       ba10cb48b Refactored in the new Job object and input output file tests for restart
       ea2786e6a Updated the way snp calling works to support, nb jobs instead of window.
       f28a83a14 Removed dup index generation since we added recalibration
       4a6c7013a Fixed the way output directories are initialized
       10f160f8a Updated version
       0f80b8630 Added pigz
       b9e96e540 Merged master
       984550597 Added PipelineUtil pacakge to all.\nCompleted BLASt implementation
       d218799b1 Implemented per-job indel realigner Implemented multi sample snp calling Implemented per-job/per chromosome snp calling
       3359fa300 Use ini for igvtools version Updated igvtools version for correct exit code
       a3860f4e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       608ae55e8 Fixed already processed lanes.
       66430dbc8 Version bump fixed execute bit
       4ed8abc09 Trial implementation of new job restart flags
       2e09503dd Implemented restart job handling
       9912c46bf Added link to ptx patch
       c98411736 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       daa8f5aba Added and updated modules
       5448d0908 Added tabix
       35fa5c7bf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       75c9afd39 Fixed many bugs
       ae1214776 Fixed many bugs
       fb3af3a27 Fixed hash generation
       42f88d87e Changed trimming
       e7eb18f8e Fixed index vs ref problem
       1dbc2b888 Fixed index vs ref problem
       84c1895ad Split metrics...we probably should split it further.
       7747bcd29 Added TRINITY_HOME
       70cd9b83f Added a parsable file containing trimming stats
       d032c5aed Fixed default params
       9ea121188  Changed metrics position in steps
       96ca681e9 Fixed parameter passing.
       ea8cecc80 Fixed problems in RNASeq denovo pipeline
       5cf889890 Added support for non-multiplex lanes....
       1d8a1903c Fixed trimmomatic bug Added usage to sampleSetup Fixed some ini params in rnaseq
       a4054926e Added options for skipping samples
       4d7e1fac0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d06a44fcb Fixed many issues and params for the RNASeq denovo
       162142abf Format stats from new format
       7b6196158 Version bump
       972fa53ad Added missing file warning
       42b94f210 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c347b1ff8 Made many fixes to the denovo RNA Assembly pipeline
       a99596806 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       05042bf81 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7ec950061 Modified permissions
       ceea8d292 Fixed mpileup bam location
       58152cf8e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       49b140986 Fixed parallele threads in single reads
       e311d8d21 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d7e10f5a4 Fixed bwa mem RG bug.
       8fa537845 Updated version
       01e0dbded Added gatk to module list
       8ff53e4ce Fixed special single issue
       1c6edb267 Added recalibration
       36cbabd14 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f0276798b Fixed generation for single ends
       660ceda44 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8ef1373e9 Fixed validation for new project layout Added BWA mem support
       89b8a332b Added snpEff and mutect Changed mode
       070bceef7 Fixed output redirection
       add31c710 Fixed multi runId resolution
       351be7c29 Cleaned up downloaded genomes
       40edd6e34 Sample setup scripts that works with miseq and hiseq and fetches the nanuq project by ID
       69f8ad49e Added gallus gallus 3
       a408a59e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b3ee34045 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f461b67a1 Fixed default ref.
       210079be1 Use parallele bwa by default
       1d973c603 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e83d0688b Fixed libs and dna pipeline to use official project hierarchy
       b0a33ed27 Fixed cleanup code
       2af9a96e8 Updated the ray module
       643a97b22 Added the exonerate module Fixed the ray module
       5c1c6197a Fixed append bug
       41a888f62 Change file mode
       12bdf7572 Fixed configurable module load to trimmomatic Fixed pass output dir, not hard coded in trimmomatic
       3f20c53b8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46ac74f60 Added missing rawReadsDir
       6f4161d61 Partially fixed raw_read support for dnaSeq Fixed markdup call in dnaSeq
       f861134b3 fixed conflicts
       c973658f3 Fixed whole genome coverage
       32b1e1698 Added module versioning, fixed paired pileup bug
       94c6f9728 Added SV to the paired pipeline
       2e99fd6c9 Added snv calling
       2bfad558b Rolled back changes because it broke existing pipelines
       0bc8182ab Added SV to the paired pipeline
       d67ceddea Added DNAC support
       79e909f7f Added DNAC support
       b9b11aabc Merge branch 'mugqic_lletourn'
       0cceb4906 Use job ids Added variant calling Added more flexibility when calling modules
       18a60abe4 Added thresholds to genomeCoverage
       a70bf88bf Merge branch 'master' into mugqic_lletourn
       275d08aaf Fixed usage
       cee2b1829 Added for snp calling
       4e9304313 Added metrics and steps to the dnaSeq pipeline
       8222088c6 Change job dependencies from names to ids.
       77a4c6184 Fixed 2 paths bugs
       30396f0f5 Don't use alternate contigs
       0c43b15f8 Changed executable attribute
       7b04c504b Keep dictionary ordering
       f08a3ea68 Fixed usage message
       41172652a Added paired snp variant calling
       a5d2e1605 Added more flexibility with tmp.
       67b2b0598 Completed first half of the validation pipeline Added adapters file for common paired adapters
       d0df7f456 Added metrics to validation Added IGV tootls Added targetCoverage
       d66f466ad Fixed header
       466071a3f Fixed issues with multi single+pair
       df027659e Added Validation pipeline
       b48e4b56c Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       01ce3d466 Fixed index generation
       0f46ad3c6 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       0067ebf0d Fixed timestamp checking
       08d14e015 Added trim counts statistics
       092311328 Fixed Trimmomatic MAJOR bug Added more settings to trimmomatic
       275cd84ad Use the right path for trimmomatic
       f652072f3 Fixed output job directory naming
       cc0736381 Fixed bad realigner output file name
       7c15b689c Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       992ef66ea Fixed bad argument
       c4479ab65 Added guillimin conf file
       5411b1709 Fixed job dependency code.
       53a75c0f0 Finished bam generation pipeline
       faf8ed155 Merge branch 'master' of ssh://bitbucket.org/mugqic/mugqic_pipeline
       4e3067430 Added realigner
       eda6d1257 Fixed uninitialized bugs
       2b178064d First pass for dnaSeq pipeline
       5b1fc8bff Fixed typo

  Louis Letourneau <lletourn@abacus1.(none)>      1 commits

       9286b83aa Fixed many bugs

  Mathieu Bourgey <bourgey1@ip03.m>      2 commits

       3b2af18a4 RNASEQDN: modify trinity to check for previous assembly + variable name change
       d49a53aa7 TOOLS: chenge permission of non-executable scripts

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      26 commits

       58770fe6e RNASEQ: change merge test
       73407b03c RNASEQ: change merge test
       54a772aae Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22f64a30d RNASEQ: add mkdir metrics in step 2
       be482de96 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ee774f9b1 goseq.R: check and remove results that are not reported in the GO.db R package
       220e9fbad lib GqSeqUtils report call change
       aa69e7ee8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       08a72c995 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       20ce1c962 formatDenovoCombinedGTF.py update
       3183fd0a3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f6ce3df50 RNAseq update
       c833eee7a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d1e8dccee RNASEQ update unstrand wiggle bug correct
       49d7bbc88 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e7dbfbc82 RNASEQ correct stranded wiggle array assignation
       352f1ed83 add the java module call at Metrics::rnaseqQC
       7825d84af add the java module call before at each picard function
       aeb1c37eb RNASEQ: add mamouth ini file
       2f858ec85 MODULES - add temporary download folder in several module install scripts
       e63d6498f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       44263b2bf rnaseq: resolve dependency pb
       9873fc20c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89fbe54e3 RNAseq modify some output location for the reporting
       27bbc10c8 STATS: correct metrics:readstat for using output trimmomatic
       772a5e5a1 RNASEQDN: old mamouth ini file by the new one that fit changes

  mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      96 commits

       ad595721c RNAseq update cuffdiff new dependency update
       278c9d25e RNASEQ: change variable nemae of bowtie fasta
       483e00dd3 RNAseq update metrics stats
       8b02e13a3 RNAseq update trimming metrics
       214936234 RNAseq correct DESseq wrong dependency
       398a12118 RNAseq upstae trimming merge stats; update edgeR
       b07481ea1 rna update
       bed41b69c rna update
       4c45b56a3 remove tmp test for maty in dnaseq.pl
       45c8ba261 rna update
       ab2a14991 R-tool: change saturation graph format
       7d052f149 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d0bd49987 RNASEQ: add format protability to goseq
       1cd3ae088 RNASEQ: add format protability to goseq
       2c7b901bd RNASEQ: add format protability to goseq
       be3833aaa RNASEQ: add format protability to goseq
       3289aac78 RNASEQ: add format protability to goseq
       971938cfb RNASEQ: add format protability to goseq
       6cc1155b0 RNASEQ: add format protability to goseq
       9f8420760 RNASEQ: rnaseq.pl update
       e31e700fd RNASEQ: rnaseq.pl update
       fc1faad49 RNASEQ: add format protability to goseq
       18ffd9a53 RNASEQ: add format protability to goseq
       e725e01c3 RNASEQ: edger update
       7fbbb5a6c RNASEQ: goseq update
       c6dbccbd6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       73a94bf1d RNASEQ: goseq update
       06b2e4afe RNASEQ: goseq update
       f93dba0d2 RNASEQ: allow non native GOseq analysis
       656c63c70 RNAseq: cuffdiff result fillter now include in the merge with fpkm
       02e9ef8b5 RNAseq: cuffdiff result fillter now include in the merge with fpkm
       bc787bd57 RNAseq update
       49ab91f48 RNAseq update
       e6515d3c9 change mugiqc_tools.install script to avoir facl conflict
       9676b39a3 update tools changes
       287a8e2db update tools changes
       e35a443d9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       58ad96b31 RNASEQ: resolve splited metrics dependency
       4008a3d92 Rtools edger old index name  not suported by the actual edger version now use
       483781a04 RNAseq remove fpkm stats as they are done also in rna-seqc
       c47230d6d RNAseq correct bugs - see BFXDEV-20 for details
       8b614c0fa RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence 2
       b5eee70a5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       79804b1c3 RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence
       3b0ce9887 RNASEQ: bowtie reference the new code should now be portable
       4b6bf8e1f RNASEQ: allow readstat on single library
       d079f598f RNASEQ: correct readstat output format
       a61e0050a RNASEQ: metrics changes bug correction & Cufflinks correction
       28285f9c7 Metrics add  argument
       a1032f161 RNASEQ - correct typo in the mergebam step
       68b309e69 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       318f27b50 TOOLS remove type blatbestHit.awk
       30a3793e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5c8851c61 TOOLS adding blastbestHit.awk blatbestHit.awk
       e8481fd53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       26152f77c RNASEQ: rnaSeq.abacus.ini correct rnaQC fasta variable name
       9bcdcf90a TOOLS: gtf2tmpMatrix.awk - adapt for line witout gene name
       57441b610 RNASEQ : adjust saturation plot title
       f05e6a747  RNASEQ : add saturation thread Number in the abacus .ini file
       1d072b2c2 Tools : make gtf2tmpMatrix.awk portable for various type of GTF file
       fbf1e9a6b tools : adding gtf2geneSize.awk and gtf2tmpMatrix.awk which were initially in tools but not anymore in the module
       2de9a907d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3ec125de7 RNASEQ : allow running without reference GTF .2
       a19c305ad Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a48116dfb RNASEQ : allow running without reference GTF
       4c9f8033d pairedVariant - finish controlFREEC allow either bam or mpileup input everthing is specified in the ini file -  BFXDEV-7
       4a2439baf R-tools change rpkmSaturation the count file need to be tab separeted
       93eda90e8 R-tools correct the polting bug of the rpkmSaturation - BFXDEV-3
       c67cdd061 rnaseq - replace my email by the  MAIL variable
       b1039a263 pairedVariant - add Control FREEC lib
       75351f19c update rnaseq.pl for single mode
       d1557826b update rpkmSaturation to the last version
       67148739d correct formatGTFCufflinks.pl
       2220b26da Merge branch 'master' into mugqic_mathieu
       8aef01eaa Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       7ff4088ef correct rnaseq ; remove special caharacter in subimti to cluster job ID varable
       39e7fd1d2 correct cuffdiff deNovo
       9db7cdf04 correct conflict btw mugqic_mathieu and master
       314eadc5a correct small bugs
       0307f2576 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       17dadbe17 correcting pindel in lib/pindel.pm and in pairedVariant.pl form dnaseq pipeline
       392774369 remove python module call in htseqCount::readCountPortable
       ed644b5f2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fcc0e4101 Merge branch 'master' into mugqic_mathieu
       d7d4284a2 correcting last bug in rnaseq
       4cb464e75 debug rnaSeq.pl
       7725d7816 trimmomatic add java module loading
       1c9c334e9 submitToCluster add default  value =  if undef
       70afb543c rnaSeq debug
       2a78c3752 resolve conflicts
       3c840f218 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       208d0cc83 debug
       cc57bb26d update
       8889bf16b rnaseq test debug; submitcluster modification; trimmomatic modifications
       bd4afe168 debugging RNAseq
       249ded651 debugging RNAseq

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      129 commits

       a51199f54 samtools allow pileup with nb de region = 1 => pas de region
       ea59d1fef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3e902cc5c dnaseq.pl && samtools allow pileup with nb de region = 1 => pas de region
       0ff4f4884 upadte rnaseq.pl: missing convertion to new job object for stranded wiggle printToSubmit calls
       ae5a4a805 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f79aadbfb dnaseq annotation bug correct and ini update
       d4e4586d0 rnaseq.pl: correct exploratoiy dependency if start at step 9
       1b1abe4db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       33df28ef8 dnaseq correct metrics restart bug & update dnaseq ini
       6db6f7035 remove tools_shed for the new tool repo
       c653f58f6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e1bc28574 metrics/SNVGraphs add input output test for restart
       72b31c8ef update Tools
       0e1b10e73 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b75840451 update cuffdiff stranded
       af5dba279 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       371205c49 switch ToolShed to Tools
       c0c1d65c7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       665a09863 RTOOLS update
       3887f719f snvGraphMetrics.R  update
       38b94dffd DNAseq.pl update for metrics
       473210956 DNAseq.pl update for metrics
       daceb8ab3 Merge branch 'dnaSeq_Report'
       c82147af6 DNAseq.pl update for metrics and report
       b31654c65 DNAseq.pl update for metrics and report
       75f28e5b7 DNAseq.pl update for metrics and report
       425b21e25 add tool_shed/python-tools/vcfStats.py
       ac96abbfd remove mater to branch conflict
       d133749a2 rnaSeq.pl edited online with Bitbucket
       23e60fc98 remove conflict between danSeq_report and master
       e67ddd8d9 add tool_shed/tools/splitSnpEffStat.awk
       e4cfaee23 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5ce57fe8 GENOMES: add Tetrahymena_thermophila.sh
       5f468e6d5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4bd157a90 sampleSetup.pl patch
       a2dbaccb6 samplesetup correction
       f4ba09563 Trimmomatic.pm edited online with Bitbucket
       00d0af028 merging matser in dnaseq_report
       7ff833eda rnaSeq.pl edited online with Bitbucket
       4556bdabc DNArepport update
       91d8b3cc0 goseq.R edited online with Bitbucket
       e2148e76a goseq.R edited online with Bitbucket
       5a202eba9 DNAreport update
       a7a6d1309 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cc7af8dec TOOLs gtf2geneSize.awk protability
       05fdcce60 RNAseq update
       8b603f9a9 RNAseq update
       58cc9c62c RNAseq update
       dcf225dba Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b69a50b8a RNAseq correct htseqcount for stranded RNA
       f45201001 RNAseq Strand-specificity correction
       7754be2af Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       bc5c5a516 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a5450eca6 TOOLS rpkmSautrationadd more flexibilty on the file format
       4960f3ea8 RNASEQ: add output directory creation in fpkm folder
       904595ff5 General : Get back the change lost in commit 0007ea2
       84818acda Revert "RNASEQ: correct dependency issue"
       d45c5d776 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d3b3bb1db Revert "RNASEQ: correct dependency issue"
       7aabc0093 save before revert
       0007ea2b7 RNASEQ: correct dependency issue
       a3c5a79f1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       eb6d023f0 RNASEQ : add deliverable + test exploratory
       f62067d74 Merged in rnaseq_lef (pull request #4)
       63a01633c RNASEQ - GqSeqUtils.pm update
       23099ae88 RNASEQ modify exploratory for the pull request #4
       80da05678 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       af49cae6e RNAseq update
       c258e2c66 REDO lost commit: rnaSeq.pl use absolute path of design file
       7e28edddb SubmitToCluster.pm edited online with Bitbucket
       79e8b6073 add log.pm lib
       19289aaa3 add logfile for cluster output file path
       ccf1fa8c6 RNASEQ: rnaseq.pl update
       5e8fcfca3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       28964581e RNASEQ: make design file path an absolute path to conflict if the working directory does match the relative path of the file given in argument
       363892b21 rnaSeq.pl edited online with Bitbucket
       6435ebeb3 RNASEQ: rnaseq.pl && cufflinks.pm update
       3e54a64e9 SubmitToCluster.pm edited online with Bitbucket
       9ac95ffc6 RNASEQ: rawCountMetrics update
       509b645f6 RNASEQ: rawCountMetrics update
       c484eed2e Merge branch 'logReport' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       6e725ad86 Merge branch 'master' into logReport
       a7faedc6b RNASEQ: mv dgeMetrics (step 10) as rawCountMetrics (step 8)
       04752e232 add log.pm lib
       3fcc28380 add logfile for cluster output file path
       c527ebc8c RNAseq: change wingzip call in the metrics
       cbbd97abb RNAseq update rnaseq.pl
       f4e321518 RNAseq: update rnaseq.pl
       2a8c65c6c RNAseq: update rnaseq.pl change matrixreadcount form step 10 to 7
       b43088aad RNAseq update rnaseq.pl
       f672d802c RNAseq: update rnaseq.pl
       5535e5084 TOOLS: add formatDenovoCombinedGTF.py & RNAseq: update cuffcompare
       a348a62d6 RNAseq: update rnaseq add cuffcompare
       cb485b8e5 rnaSeq.pl edited online with Bitbucket
       177aa3476 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0c0f4a12d TOOL: Adding thepython tools folder and the getRefBinedGC.py script in it
       77b49e292 Merge branch 'mugqic_mathieu'
       1ebe60bc6 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       22ad9e5d7 update RNA
       36d735db0 RNASEQ: add correct cpu request for alignment
       e17186192 correcting sv code small bugs
       553a3b00a finish breakdancer filter and  add pindel to dnaSeq
       ce1cd9179 start adding breakdancer filter and pindel to dnaSeq
       59fa224a1 resolve Htseqconfilct
       0ad635022 merge with mugqic_mathieu
       c1a278a7b debug the code
       1498b149b Merge branch 'master' into mugqic_mathieu
       c2b1cc49b resolve lib/SubmitToCluster.pm conflict
       4df938ad6 resolve HtseqCount conflict
       3268ff41a try to update local master to real master
       97b44321e everyhting except delivrable are done; debugging
       a0a55447d  goseq done; 1st test on going
       eb5785eda metrics: done; General module nomeclature conversion
       5efae8134  metrics: updates; SAMtools: add viewFilter sub function ; tophatBowtie: bug correction and simplification
       2e48636d4  metrics updates
       be331d2ed matrix done; DGE done
       4223cc79f fpkm done; DTE done ; DGE on going ; metrics on going
       08778dcce correcting wiggle; fpkm and rawCount done ; DTE on going
       1ed2d7f29 wiggle done
       cd318c11c update metrics ini.file and few correction tophat/botwite
       052889d74 merging ok ; metrics started; new functions and changes in the picard lib
       2f6709a9c merging done ; metrics started
       66ea6a475 update to the master branch
       7166f7456 Tophat-botwie done; merging started
       f9662b881 uppdating TopHat lib and rnaSeq.abacus.ini
       634a238f5 uppdating my branch
       9f81ed088 Starting bowtie/tophat lib
       43c3e7a93 Global dependency and RNAseq Triming and submit with working directory argument
       d6a6324f7 starting the rnaseq pipeline

  Mathieu Bourgey <mbourgey@abacus2.(none)>      1 commits

       f36f0e216 rnaseq.pl debug; Metrics debug; Picard debug; SaMtools debug; SubmitToCluster debug; TophatBowtie debug

  Maxime Caron <mcaron@abacus1.(none)>      2 commits

       4ad23e113 test
       0b8c7fb8b test

  Pascale Marquis <pmarquis@abacus2.(none)>      5 commits

       168ea11c6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       221a12ef4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b91d4dad8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0b27236e2 zorro.sh
       9f170dcbb Tetrahymena_thermophila.sh

