44 tags, 10124 commits

HEAD        Fri Sep 8 17:48:23 2023 +0000        0 commits

4.4.4        Fri Sep 8 18:46:58 2023 +0000        14 commits

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      11 commits

       0dc54c46 Merged in release_4.4.4 (pull request #451)
       70066e40 Merged in mutect_fix (pull request #448)
       a8545574 Merged dev into mutect_fix
       fa9ea1e6 Version bump to 4.4.4
       657dae58 GenPipes tumor_pair : interval_list job name and dependency issues
       70c286ef Merged in interval_list_fix (pull request #447)
       bd975f7a GenPipes : change job name for interval list creation in recalibration step
       6e3217b8 VERSION bump to 4.4.4.-beta
       7d0f747f Merged master into dev
       718ffd77 Merged in release_4.4.3 (pull request #445)
       7b1ac61a Version bump to 4.4.3

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      3 commits

       422e4afb Merged in kallisto_index_renaming (pull request #449)
       47db2c7c rnaseq_light - kallisto - Upgrading kallisto to latest + increasing resources for cit
       49f00cc7 rnaseq_light - kallisto - index renaming to match kallisto version

4.4.3        Wed Aug 23 18:37:11 2023 +0000        338 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      43 commits

       cf60298e Utils - added dump_ChangeLog.sh in utils
       2c745e14 DNASeq - removed `sed` call in manta_sv_calls command
       80ecce89 Tumor Pair - Fix run_pair_multiqc job name
       2f31d6bf Tumor Pair - fix multiqc step
       a20dfc84 Tumor Pair - Fixed multiqc step dependencies
       98aac8e8 RNASeq - Fixed multiqc step dependencies
       8172a2ad DNASeq - Fixed multiqc step dependencies
       58dfc13d DNASeq - Fix gcbias output deps for multiqc
       187a10cf HiCSeq - updated mugqic_tools to 2.12.1 to fix hicrep
       9fdf08ec Core - Job : fix report_files in concat_jobs
       a5f1fbb2 GenPipes - updated dates in licence
       352863a1 Tumor Pair - Updated multiqc prep steps and job
       f83b6e9a DNASeq - stop index bams with picard_mark_duplicates, now using sambamba index instead to create bam.bai files
       f265ce31 RNASeq - Fix rnaseqc2 call : now uses sambamba to index the .bam as .bam.bai instead of picard creating .bai
       b19bfb05 Tumor Pair - Fixing MultiQC jobs in cases of missing dependencies
       db6c00eb RNASeq - cancer : Updated pcrg version to 1.0.3
       51ea7fdb BFX - fixed jsonator update
       5f89f78d Utils - remove JSON file when empty
       876ad395 Tumor Pair - updated bam input selecton for manta_sv and gatk_variant_annotator steps
       c175aef3 RNASeq - varianats/cancer : update docstring for merge_vcfs step
       ea76e9da RNASeq - varinats/cancer : updated merge_hc_vcf removing contigs/super contigs when merging vcf
       7258e70b RNASeq - variants/cancer : allow gatk_indel_realigner step to be skipped
       31ddad76 DNASeq / Tumor Pair - Picard SamToFastq : improved passing of resource parameters
       e09c760c Resources - updated SMRTLink install script with latest version
       9bbaf26d Tumor PAir - fixed typo in gripss wrapper
       fe29ad4c Version bump to 4.3.2
       ca99865e Merge remote-tracking branch 'origin/dev' into release_4.3.2
       a34fda40 Core - Scheduler : corrected `memory` call in `cpu` from dev
       80afe6df Version bump to 4.3.1
       a647b011 Merge remote-tracking branch 'origin/dev' into release_4.3.1
       57667889 Version bump to 4.3.0 - for real
       937a5ca4 Version bump to 4.3.0
       69d74b05 Merge remote-tracking branch 'origin/dev' into release_4.3.0.1
       3f33a35a Version bump to 4.3.0
       76b0b626 Merge remote-tracking branch 'origin/dev' into release_4.3.0
       1e06b04a Version bump to 4.2.1
       460b6a3c Version bump to 4.2.0
       695798b1 Generating READMEs for 4.2.0
       9947d057 Version bump to 4.1.3
       d7f19424 Merge remote-tracking branch 'origin/dev' into release_4.1.3
       1ef5a9ec Version bump to 4.1.2
       7997ec59 Merge remote-tracking branch 'origin/dev' into release_4.1.2
       c83a283f Version bump to 4.1.2

  ehenrion <edouard.henrion@mcgill.ca>      77 commits

       f7b553c3 Resources - Updating Conpair and Mugqic_Tools installation scripts with latest versions
       f592c031 Merged in multiqc_cit_fix (pull request #427)
       b8244276 HiCSeq - updated mugqic_tools to 2.12.0 to use the fixed version of hicrep.R (for no down-sampling cases)
       2e15b7fa Merged in update_licence (pull request #422)
       4d8e1f59 RNASeq - fixed rnaseqc2 output dependencies
       fdfb9b43 Modules - Updated binary wrapping done by installation scripts
       369653aa Modules - Added bamixchecker installation script
       171a7445 Module - Updated mugqic_tools.sh with latest version
       a1d73e3c Utils - fixed watch_portal_folder for some cases with empty json : stop throwing error and remove the json instead
       4b8c7bde Resources - AGeNT : fixed typo in module
       dcdd1814 fix typo
       fcc76d83 GenPipes - improving job2json
       a458c390 GenPipes - strengthening the jsonator
       f869fdc4 Tumor Pair - refining resources assigned to gridss_paired_somatic jobs
       a6e9a3b1 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       a756c811 Tumor Pair - Fixed inputs for "sequenza.create_seqz" jobs
       3d9b5c00 Tumor Pair - fix linx_plot output dependency to avoid useless restart of the step
       55ceb6d9 HiCSeq - RobusTAD : fixed output name to fix restart issue
       b1d7edca Utils - made mugqicValidator.py Python3-compatible
       59f2a414 Tumor Pair - fixed gatk4.ini for cobalt, amber and purple
       5ed67ce2 Tumor Pair - fixed gripss.py arguments in bfx
       6d511de6 Tumor Pair - fixed settings for mutect2 when using gatk4
       1a0c1066 Core - Add details to some sanity check error messages
       2419d273 RNASeq - fixed stringtie-merge dependencies
       47cf66ad Merged in tumor_pair_dependency_fix (pull request #413)
       acd69760 Tumor Pair - Fixed gripss calls
       6793fdf9 Version bump to 4.4.2-beta
       f40eae05 Merged in release_4.4.1 (pull request #412)
       e9cd38bb Version bump to 4.4.1
       0d2e78fc Merged in release_4.4.1 (pull request #411)
       412c2e3d GenPipes - in prep for bug-fix release 4.4.1
       c27edf15 Resources - update mosdepth version to 0.3.3 in install script
       6a5222f2 EpiQC - correcting chromimpute_preprocess command : do not use "rm" and force symlink creation with the "-f" flag
       087fd9b3 Resources - adding / updating install scripts
       252b6ee5 RNASeq Light - Added mkdir to the kalliso command
       07f45a23 Resources - adding 'gget' install script and updating others
       ed8e359e Version bump to 4.4.1-beta
       20c470b7 Merged in release_4.4.0 (pull request #410)
       8cabdc1c Version bump to 4.4.0
       1c294768 Merged in release_4.4.0 (pull request #409)
       bc806684 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       468975cb updating README
       b0db9b96 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       67fe5bbe RNASEq Light - Fixed Kallisto step
       9e014170 Resources - added HMMcopy package in defaults R libraries
       ffd223d0 RNASeq Light - Adding kallisto wrapper in bfx
       91bca69e RNASeq Light - revamping the call to Kallisto to have per gene counts instead of per readset
       d83071de Resources - Adding install script for hmmcopy-utils and ichorCNA
       79b74d03 Resources - Adding AGeNT install script and updating script for VarDict
       d2a8fef7 RNASeq Denovo Assembly - fixed trinity job dependencies
       9629f9a3 RNASeq Light - removing kallisto output folder before runnig kallisto to avoid mixed or bad results in case of pipeline restart
       eb7df8b1 EpiQC - fixed signal_to_noise dependencies
       48f936c1 DNASeq High Coverage - increased default walltime for gemini_annotations
       43121ea9 HiCSeq - fixing "homer_tag_directory" outputs
       4f7bfea4 Merge remote-tracking branch 'origin/dev' into release_4.4.0
       cf282891 GenPipes - updating READMEs before release
       70d27f7c Resources - Adding some new install scripts
       699683ca Resources - updated version in some install scripts
       3d0c96d9 EpiQC - improved restart mechanisms
       aa742e32 Merged in release_4.3.2 (pull request #389)
       ed7c1835 Merged in release_4.3.2 (pull request #388)
       118309c8 Merged in release_4.3.1 (pull request #383)
       baeee47d Merged in release_4.3.1 (pull request #382)
       f0feb9ad Merged in release_4.3.1 (pull request #381)
       e5980e36 Merged in release_4.3.0.1 (pull request #367)
       110f8ad7 Merged in release_4.3.0.1 (pull request #366)
       2e0a7ecc Merged in release_4.3.0 (pull request #364)
       390169f3 Merged in release_4.3.0 (pull request #363)
       47523b32 Merged in release_4.2.1 (pull request #359)
       9f80a8c9 Merged in release_4.2.0 (pull request #354)
       927e40e4 Merged in release_4.2.0 (pull request #353)
       a546428e Merged in release_4.1.3 (pull request #352)
       f1e6716e Merged in release_4.1.3 (pull request #351)
       a1be2b96 Merged in release_4.1.2 (pull request #314)
       e5a64db9 Merged in release_4.1.2 (pull request #313)
       99b40bf6 Merged in release_4.1.2 (pull request #312)
       dc6f51e6 Merged in release_4.1.2 (pull request #311)

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      118 commits

       52df2e50 GenPipes tumor_pair : add md5 for recalibration output bam
       3e62f88d GenPipes rnaseq : add option to force complete restart of star_fusion
       5fd625d6 GenPipes rnaseq : adjust dependencies for run_star_fusion
       d76554e5 pipelines/rnaseq/rnaseq.py
       7b508019 bfx : fetched gatk4 updates from MOH_permissions branch
       7c49ea5d GenPipes rnaseq: fix wiggle job name error
       5055a0d2 GenPipes tumor_pair : improve documentation
       16a7cc8a GenPipes : updating documentation with Shaloos suggestions
       6e3f2c76 bfx : samtools quickcheck move output to same line as input
       40db2942 bfx : samtools quickcheck fix options
       1a7b4a28 bfx : samtools quickcheck add missing comma
       c5bf16b4 GenPipes nanopore : import samtools from bfx
       7d9e01a0 bfx : add samtools quickcheck
       f77c7c2f GenPipes nanopore : add bam check to minimap2_align to avoid false 0 exit
       b091b4cf GenPipes tumor_pair : increase cit walltime for manta
       3f27c218 Genpipes hicseq : increase cit walltime for identify_peaks
       8e4ddc05 GenPipes methylseq : increase bowtie2 threads during bismark_align for cit
       1f131de6 GenPipes dnaseq : add default walltime back in
       0bea2387 bfx : deliverables revert change, causing issue on abacus
       addefee4 GenPipes rnaseq : arriba symlink troubleshooting
       aaf7037d Genpipes rnaseq : fixed typo run_arriba
       8de3b516 GenPipes rnaseq : fix multiqc symlink in run_arriba
       ced167c8 bfx : change deliverables md5 output to fix restarts
       70b9907e GenPipes rnaseq : change multiqc symlink in run_arriba to fix restarts
       d0da2678 GenPipes : cit adjustments for walltime
       27f92861 GenPipes methylseq : reduce cit walltime for new bissnp version
       a84f76f2 GenPipes methylseq : update bissnp version
       7068c8e2 GenPipes rnaseq : increase cit walltime gemini_annotations
       5be926ab GenPipes rnaseq : revert tmp dir change for gemini, did not speed up analysis
       a588b66d GenPipes cit : small walltime and ram adjustments
       063c60d9 GenPipes hicseq : rm default cit walltime, too low for design file step
       4da926d3 GenPipes methylseq : update bissnp version
       6190d8d4 GenPipes hicseq_hic : increase cit mem for identify_peaks
       f1971b34 GenPipes rnaseq : fix error in gemini_annotations tmp dir change
       fcbbd6ec GenPipes : cit ini updates for beluga
       3e777360 GenPipes rnaseq : use tmp_dir for gemini annotations
       2999d0bc GenPipes methylseq : reduce mem for filter_snp_cpg, oom resolved
       ed53a240 GenPipes dnaseq_hcov : reduce resource request for symlink
       71d6390a Merge branch 'improve_cit_inis' of bitbucket.org:mugqic/genpipes into improve_cit_inis
       e09758ab GenPipes methylseq : use tmpdir for filter_snp_cpg to resolve OOM errors
       65d34518 Genpipes rnaseq : resolve conflict
       dbe3e99c bfx : revert removal of quiet flag in manta and strelka2
       71e949ae GenPipes hicseq : cit mem adjust
       b7e0c294 GenPipes hicseq : cit walltime adjust
       27beaad6 GenPipes hicseq : adjust cit mem for homer_tag_dir
       49206402 GenPipes hicseq : cit mem and walltime adjustment
       82a233e7 GenPipes : small cit adjustments
       33e9a345 GenPipes cit : fix errors in ini
       58a5fde5 GenPipes methylseq : increase cit mem for filter_snp_cpg
       d8fd13d6 GenPipes hicseq : adjusting cit resource requests
       a154ba68 GenPipes ampliconseq: adjust cit resource requests
       a48d26d2 GenPipes nanopore : cit resource adjustments
       be9e2cef GenPipes rnaseq_denovo : cit walltime adjustment
       b604cee3 GenPipes methylseq : adjust cit resource requests
       f8db6f7b GenPipes rnaseq_denovo : adjust cit mem
       13eab33e GenPipes rnaseq_denovo : adjust jellyfish cit mem
       789324cc GenPipes rnaseq_light : adjust cit report mem
       bc0bf9bf GenPipes rnaseq_denovo : adjust cit mem seq2fun_pathway
       944fedf2 GenPipes rnaseq_denovo : fixed typo
       1509fa4d GenPipes rnaseq_denovo : adjusting cit resources for trinity and seq2fun
       004b4328 GenPipes rnaseq : add multiqc by_sample, adjust walltime
       61cc2df9 GenPipes rnaseq : rm output dir for run_arriba to fix issue with restarts after timeout
       597c31bf GenPipes dnaseq_light : adjust cit resource requests
       c9a76714 GenPipes dnaseq : adjust cit cnvkit walltime
       9f0eaf5f GenPipes dnaseq : adjust cit cnvkit correction settings
       95f4f64d GenPipes dnaseq : adjust cit ini for dnaseq_sv
       36d09308 GenPipes dnaseq : adjust cit ini for dnaseq_light
       cb7479c1 GenPipes rnaseq : fix typo
       3e0f1776 GenPipes dnaseq : cit default cluster walltime
       c62be41e GenPipes dnaseq : adjust cit resource requests for mpileup protocol
       5b21f181 GenPipes dnaseq : reduce walltime request for sym link
       7dc2c37d GenPipes dnaseq : mugqic reduce cit resource requests
       a11e1d7f GenPipes tumor_pair : increase cit walltime vardict_paired
       040acdc6 GenPipes rnaseq : add option to generate individual multiqc reports by sample
       fe38d8d3 GenPipes tumor_pair : rm slurm_tmpdir for manta
       a07ff8bc GenPipes tumor_pair : testing slurm_tmpdir for manta
       dae1d6bd GenPipes tumor_pair : testing slurm tmpdir for manta
       33342c35 GenPipes tumor_pair : adjust cit ini for vardict
       a9c6d462 bfx : manta - remove --quiet for troubleshooting
       dba72aaa GenPipes tumor_pair : fixed typo
       376761d7 GenPipes tumor_pair : add sed cmd to remove email from manta
       474cab94 GenPipes tumor_pair : add sed to strelka to investigate timeouts
       db730355 GenPipes tumor_pair : reduce bwa_mem cpus, increase strelka time
       afe1f5b0 bfx : remove quiet flag from strelka run for troubleshooting
       3adc05a7 GenPipes rnaseq : cit ini updates
       3ced8462 GenPipes rnaseq : adjusting cit ini
       a098e83e GenPipes rnaseq : adjust cit resources in ini
       d77f50df GenPipes rnaseq : adjust cit walltimes
       4dc49323 GenPipes rnaseq : adjust star_align cit walltime
       4321b40f GenPipes tumor_pair: changed ini section skewer to skewer_trimming
       4e69f6eb Genpipes chipseq : adjustments in cit ini
       0b7c7f2e GenPipes rnaseq : reduce cit resource requests for rnaseq_cancer
       64cbe652 GenPipes chipseq : adjust cit resources
       f573d5f9 GenPipes tumor_pair : adjsut cit resources
       6e3f4de6 GenPipes rnaseq : adjust cit resource requests
       8abb0818 GenPipes rnaseq : change rnaseq-merge to rnaseq_merge job name
       6d16303e GenPipes rnaseq : adjust multiqc walltime
       7e5f91e0 GenPipes tumor_pair : adjust cit resources
       ecbf35e6 GenPipes tumor_pair : improve cit resource requests for tumor_pair_fastpass
       44dbfa29 GenPipes tumor_pair : reduce cit resource requests for tumor_pair_sv
       5c45e46f GenPipes dnaseq : reduce resource requests for symlinks, trimming
       93ba9d6c GenPipes : update cit ini resources
       1916eb07 GenPipes covseq : update resource requests in cit ini
       8f95886f GenPipes chipseq : improve resource requests
       43767b82 Genpipes : methylseq - clarifying multiqc comment
       df70dbfa GenPipes : methylseq - update mkdir commands
       9614d14f GenPipes : methylseq - add comments to explain multiqc options
       918d8ac3 bfx : dragen.py - add report_files to outputs
       1ce17811 bfx : dragen.py - extend outputs to avoid list within list
       2da1c94d GenPipes : methylseq - change automatic sample naming for trimmomatic multiqc
       60baae95 GenPipes : methylseq ini add multiqc options
       8961e31d GenPipes : methylseq - ihec table to multiqc format
       04d0de10 GenPipes : methylseq - add multiqc to bismark align
       3bdff148 bfx : make sambamba options ini section not required
       f2b9ccb0 GenPipes methylseq : add multiqc to pipeline
       37608f3d GenPipes methylseq : add symlink_dragen_metrics section to ini
       776f679f GenPipes methylseq : add multiqc to base ini
       97ba46ee bfx : add metrics files to expected dragen output and report files

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      94 commits

       1589b386 Merged in release_4.4.3 (pull request #444)
       de07ec66 version bump to 4.4.3
       08a6d65e update README
       14242b6b update pipeline READMEs prior to release
       58f611b2 Merged in cit_fixes (pull request #443)
       eac2d6a0 Merged dev into cit_fixes
       de5b7b84 GenPipes rnaseq : cit walltime fixes
       3389502c GenPipes cit : adjust mem and walltime
       ca2dd54d GenPipes cit : adjustments to fix oom and timeouts
       1dd76522 Merged in readme_improvements (pull request #442)
       87678ae9 Merged in prepare_table_fix (pull request #441)
       449261a4 Merged in star_fusion_restart (pull request #439)
       dea5d5ea Merged in tumorpair_md5 (pull request #440)
       010f817c GenPipes covseq : fix version of covseq_tools module
       75bca38c Resources : update covseq_tools install script for new version
       36799c9b GenPipes covseq : bump version of covseq tools
       9d47d35d GenPipes covseq : increase cit mem for prepare_report job
       45d1bd46 GenPipes covseq : make grep specific with -w flag, avoids clashes with closely named samples
       3eab9c34 Merged in wiggle_fix (pull request #438)
       892d9717 Merged in minimap2_align_check (pull request #437)
       007c1d41 Merged in improve_cit_inis (pull request #435)
       e25ea829 Merged in improve_cit_inis (pull request #433)
       c7b8d537 Merged dev into improve_cit_inis
       5b835c14 Merged in rnaseq_moh_multiqc (pull request #432)
       e332b3e4 Version bump to 4.4.3-beta
       765847eb Merged master into dev
       42492bf5 Merged in release_4.4.2 (pull request #431)
       b6f7e639 Version bump to 4.4.2
       485c806f Merged in release_4.4.2 (pull request #430)
       ed55718f Merge remote-tracking branch 'origin/dev' into release_4.4.2
       8888ad03 update pipeline READMEs prior to release
       b1c13898 Merged in methylseq_multiqc (pull request #429)
       15780ffb Merged in trimmomatic_multiqc (pull request #428)
       b2387dcc GenPipes : add symlink for trimmomatic log to use with multiqc
       5a068a98 GenPipes methylseq : increase mem for filter_snp_cpg step in cit
       9f728b78 Merged in rm_hard_clip (pull request #426)
       a0b21222 GenPipes tumor_pair : loop over report_files output for multiqc symlinks
       1850d23d bfx : add report_files output to picard metrics
       8726b5c3 GenPipes tumor_pair : multiqc - add sample name to qualimap symlinks
       271b96ec GenPipes : rm sequencing_center refs to wrong centre, add underscores to McGill_Genome_Centre
       7eb52d45 GenPipes rnaseq : rm bam_hard_clip step
       f40ed17b GenPipes : rm wrong sequencing centre, add underscores to McGill_Genome_Centre
       5c6fe1c3 Merged in cit_multiqc_fixes (pull request #424)
       621b0f6a Genpipes methylseq : increased mem for bissnip and filtering steps during cit
       4b739b5a GenPipes tumor_pair : cleaned up commented out and unnecessary lines
       03c3e8ac GenPipes epiqc : add mkdir to chromimpute convert step
       6e1c8747 Genpipes tumor_pair : fix path for is_gz file check
       68baa958 GenPipes tumor_pair : remove unnecessary os.mkdir call
       0eae1f99 Genpipes tumor_pair : create multiqc report dir as part of jobs
       859c1cc7 GenPipes tumor_pair : fix multiqc symlink creation for qualimap outputs
       8ee7b43d Merged in rnaseq_multiqc (pull request #423)
       21db3e9d GenPipes rnaseq : multiqc module order
       1128cb75 GenPipes rnaseq : fix relative path for input of rseqc.tin step
       477a6893 bfx : rm abspath in input for rseqc.tin
       cd04bf45 GenPipes rnaseq : fixed interactive option for multiqc in ini
       f448b924 GenPipes rnaseq : fixed symlink in run_arriba, change dependency exit status for multiqc
       916d36de GenPipes rnaseq : fixed symlink in run_arriba
       efd63a17 GenPipes rnaseq : merged and resolved conflict with dev
       13ab19dd GenPipes rnaseq : fix symlink relative path
       8fb9e184 GenPipes rnaseq : make multiqc interactive
       ba66e8b9 GenPipes rnaseq : fixed symlink names for multiqc inputs
       70db4203 GenPipes rnaseq : removed stray commas, fixed typos
       30ce0929 GenPipes rnaseq : removed stray comma
       c56ce402 GenPipes rnaseq : fix job name
       e851d87b GenPipes rnaseq : fixed another relative path command
       c02906b4 GenPipes rnaseq : added missing link directory assignment
       b1796d83 GenPipes rnaseq : fixed relative path command
       451dbf41 GenPipes rnaseq : added missing parenthesis
       d92e2603 GenPipes rnaseq : added missing comma
       eae40ef6 GenPipes rnaseq : add multiqc section to base ini
       97a255fd GenPipes rnaseq : add multiqc to pipeline
       9481206e bfx : add log files to expected output for multiqc integration
       baa0813b Merged in tp_ensemble_fix (pull request #421)
       72b3853f GenPipes tumor_pair : add remove command before bcbio_ensembl call to avoid silent error when output exists
       70b3e161 Merged in chipseq_blacklist_filter (pull request #418)
       63d8fe7a Chipseq : add blacklist path to ini
       02a2c734 Merged in readme_updates (pull request #417)
       1d18f5e3 GenPipes chipseq : update README to include blacklist filtering step
       6ac84851 GenPipes chipseq : added input selection to cram, variant calling steps
       5e3f6139 GenPipes chipseq : fix inputs for diffbind and run_spp depending on blacklist removal step
       2c510bd8 bfx differential_binding : add option to specify alignment file extension
       b6cbc84a GenPipes chipseq : fixed various typos, add multiqc inputs
       86f54efe GenPipes chipseq : add blacklist and bedtools module to base ini
       30bc6ff0 GenPipes chipseq : fixed typo
       3660c64d GenPipes chipseq : add bedtools intersect step to pipeline, add ini section
       45093bbf GenPipes chipseq : add bedtools intersect to remove blacklist reads prior to peak calling
       b6fb9469 GenPipes : improved README for dnaseq pipeline
       f47ca3af Merged in Mareike-Janiak/commonpy-edited-online-with-bitbucket-1681414777818 (pull request #416)
       bc8ecef5 common.py edited online with Bitbucket - fix for sambamba merge dependency
       420681f8 Merged in sambamba_merge_fix (pull request #415)
       be8788f4 GenPipes common : simplified rm command before sambamba_merge
       d6b84fb5 GenPipes common : add comments to rm command before sambamba merge
       05957e7a Common : sambamba_merge_sam_files fix typos
       f45b8300 Common sambamba_merge_sam_files : remove any existing file/link before merge

  Pascale Marquis <pascale.marquis2@mcgill.ca>      1 commits

       88ebed7a Merged in new_branch_pascale (pull request #420)

  Pascale Marquis <pmarquis@Weigela.local>      1 commits

       e622624b changed the cpulimit

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       f8296f6c Merged in MOH_permissions (pull request #436)
       666fa64b kallisto upgrade

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       cf3d7277 Merged in release_4.2.1 (pull request #358)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       c7b3aa6d Merge remote-tracking branch 'origin/dev' into release_4.2.1

DOvEE-4.3.3        Thu Jul 20 13:24:33 2023 -0400        91 commits

  ehenrion <edouard.henrion@mcgill.ca>      2 commits

       5dc6ac85 Merged in dovee_multiqc (pull request #414)
       24f8c8f0 DOvEE - set multiqc_inputs as a pipeline attribute and make changes accordingly

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      79 commits

       51d29451 Version bump to 4.3.3 changelog
       cb3ace90 DOvEE : add dump_ChangeLog.sh from dev
       df08b457 Version bump to 4.3.3
       6a2e1547 DOvEE : update job2json to work with dovee user
       a8d140f9 Genpipes dovee : update main README
       386226da GenPipes dovee : add README
       03850c76 GenPipes dovee : cleaning up old comments
       dc28a034 GenPipes dovee : update filepaths in ini
       1cc87e8a cleaning up
       01a9d91a bfx : add missing output file to locatit
       cf82e8de DOvEE : added parsing of locatit and creak output metrics, added to multiqc
       22a18b43 DOvEE : changed options for multiqc and bcftools_stats in ini
       6d46406d DOvEE : add parsing of creak metrics, picard gcbias metrics to pipeline
       616abcae DOvEE : add picard collect gcbias metrics to ini, add multiqc options
       7e156953 bfx picard2 : fixed typo in collect_gcbias_metrics - qcbias = gcbias
       e243f8fa bfx agent : add stats file to expected creak outputs
       67161bd9 DOvee : add creak and subsample sections to ini
       9f0a8e76 dovee : add creak dedup and samtools subsample steps
       caee1e3c bfx : add creak to agent wrapper
       c0ee6560 dovee : add sample source to symlink names for multiqc display
       0caca75e dovee : change mosdepth bed file used
       7a5287a2 DOvEE : add bam matching step with conpair to that samples are from same patient
       40d46af7 Merge branch 'DOvEE_genpipes_Sflag' into DOvEE_genpipes
       f18819e1 DOvEE : modified bwa mem settings in ini
       32f21dd3 DOvEE : add option to picard HS metrics to ignore validation errors
       732ebaef DOvEE : use sorted locatit output to side step sort error
       86218505 DOvEE : update multiqc to look for symlink inputs in one directory, set order of modules in report
       d04276f2 DOvEE : change multiqc job dependencies to run even if jobs fail
       8c6cecd6 DOvEE : multiqc inputs changed to directories to avoid error with failed samples
       3adb81cc DOvEE : fix error in trimmed read symlink creation
       829ff3ca DOvEE : add --interactive option to multiqc to prevent flat plots
       7dbc1405 DOvEE : fix bed files used for mosdepth
       7cf8ad0b DOvEE : add readset level directories to mapping step
       20a7fdd9 DOvEE : add tmp directory to samtools sort
       e9dc8019 DOvEE : make files executable
       7e2daa44 DOvEE : remove -K flag from locatit call to avoid keeping tmp files
       7f32db9d DOvEE vardict : add bcftools stats step
       210429f6 DOvEE : add bcftools stats to vardict protocol
       d4401227 cleaning up, fixing indents, copyrights
       eac8909c dovee : cleaning up and adding comments
       54cd4c4a adding copyright sections
       01913ce8 dovee : added option to distinguish between brush and saliva via modified version of design file
       841ae379 small dovee ini fixes
       adb834e6 add output file to mosdepth
       ecd7b09d created dovee sample pair parser based on tumor pair system
       b9b8f34f adding multiqc
       eed14380 dovee samtools sort : changed to avoid restart errors
       bd3c654f dovee metrics : add picard HS metrics
       a3984061 dovee : added metrics steps
       8764ebbc mosdepth : fixed typo
       4f0edd44 trimmer : edits to STATS rm command
       2a211153 trimmer : changed mv to ln to prevent job restarts
       dda4c1b5 dovee : adjust resources
       4f6a4283 dovee ini : ichor version change, other small fixes
       5e3db81e ichorCNA : fixed typo
       83338513 dovee : merge sort & index steps, update sort to look for inputs from previous jobs, add hybrid dedup for copy-number protocol
       cebe94f3 dovee ini : add/fix filepaths
       e3a2fe4a vardict : add module call for vardict_single
       74431a8f dovee small fixes
       6579f709 dovee resource and other ini updates
       95eb5438 samtools merge : change module call and add other_options
       eed26dd5 hmm : add module call
       cd9d0eac tweaked job names to be consistent with ini
       4cea397d updating resources and filepaths
       0d1da041 dovee : add bam merge to allow for top ups
       42146308 ichorCNA : removed # from module call
       dc9aa5d6 bwa mem : remove single end option
       2e60edbe update bwa options, update agent calls
       de382e3b combine trimmer and locatit into single agent script
       cc534520 trimmer and locatit updates now that modules are installed
       3935a130 DOvEE : add placeholder wig file paths
       08a0c453 DOvEE : fixing small errors
       7c697e2a DOvEE : add init
       7bb7402f fixed small errors and indentation throughout
       47894372 fixed indentation
       3b0fe908 moved comment
       c7f73451 fixes in pipeline and ini for new dovee pipeline
       f9495d5e new trimmer trimming function added
       912643ea new locatit dedup function added

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      9 commits

       3c7f0b42 adding fastp.py from run_processing branch to dovee
       6b549b30 mosdepth : add to genpipes
       9f2ae202 trimmer : remove STATS file if already exists
       84a8ff7b DOvEE : add tumor pair pairing system
       27a4d74c DOvEE : add pairs file option
       9db0a250 GenPipes DOvEE : add copy number protocol to pipeline
       921ca329 GenPipes DOvEE : add hmm and ichorCNA sections to ini
       41bb80d8 GenPipes DOvEE : create script for ichorCNA
       f8fdbb48 GenPipes DOvEE : create script for hmm readCounter

4.4.2        Wed Jun 21 15:59:13 2023 +0000        204 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      43 commits

       cf60298e Utils - added dump_ChangeLog.sh in utils
       2c745e14 DNASeq - removed `sed` call in manta_sv_calls command
       80ecce89 Tumor Pair - Fix run_pair_multiqc job name
       2f31d6bf Tumor Pair - fix multiqc step
       a20dfc84 Tumor Pair - Fixed multiqc step dependencies
       98aac8e8 RNASeq - Fixed multiqc step dependencies
       8172a2ad DNASeq - Fixed multiqc step dependencies
       58dfc13d DNASeq - Fix gcbias output deps for multiqc
       187a10cf HiCSeq - updated mugqic_tools to 2.12.1 to fix hicrep
       9fdf08ec Core - Job : fix report_files in concat_jobs
       a5f1fbb2 GenPipes - updated dates in licence
       352863a1 Tumor Pair - Updated multiqc prep steps and job
       f83b6e9a DNASeq - stop index bams with picard_mark_duplicates, now using sambamba index instead to create bam.bai files
       f265ce31 RNASeq - Fix rnaseqc2 call : now uses sambamba to index the .bam as .bam.bai instead of picard creating .bai
       b19bfb05 Tumor Pair - Fixing MultiQC jobs in cases of missing dependencies
       db6c00eb RNASeq - cancer : Updated pcrg version to 1.0.3
       51ea7fdb BFX - fixed jsonator update
       5f89f78d Utils - remove JSON file when empty
       876ad395 Tumor Pair - updated bam input selecton for manta_sv and gatk_variant_annotator steps
       c175aef3 RNASeq - varianats/cancer : update docstring for merge_vcfs step
       ea76e9da RNASeq - varinats/cancer : updated merge_hc_vcf removing contigs/super contigs when merging vcf
       7258e70b RNASeq - variants/cancer : allow gatk_indel_realigner step to be skipped
       31ddad76 DNASeq / Tumor Pair - Picard SamToFastq : improved passing of resource parameters
       e09c760c Resources - updated SMRTLink install script with latest version
       9bbaf26d Tumor PAir - fixed typo in gripss wrapper
       fe29ad4c Version bump to 4.3.2
       ca99865e Merge remote-tracking branch 'origin/dev' into release_4.3.2
       a34fda40 Core - Scheduler : corrected `memory` call in `cpu` from dev
       80afe6df Version bump to 4.3.1
       a647b011 Merge remote-tracking branch 'origin/dev' into release_4.3.1
       57667889 Version bump to 4.3.0 - for real
       937a5ca4 Version bump to 4.3.0
       69d74b05 Merge remote-tracking branch 'origin/dev' into release_4.3.0.1
       3f33a35a Version bump to 4.3.0
       76b0b626 Merge remote-tracking branch 'origin/dev' into release_4.3.0
       1e06b04a Version bump to 4.2.1
       460b6a3c Version bump to 4.2.0
       695798b1 Generating READMEs for 4.2.0
       9947d057 Version bump to 4.1.3
       d7f19424 Merge remote-tracking branch 'origin/dev' into release_4.1.3
       1ef5a9ec Version bump to 4.1.2
       7997ec59 Merge remote-tracking branch 'origin/dev' into release_4.1.2
       c83a283f Version bump to 4.1.2

  ehenrion <edouard.henrion@mcgill.ca>      77 commits

       f7b553c3 Resources - Updating Conpair and Mugqic_Tools installation scripts with latest versions
       f592c031 Merged in multiqc_cit_fix (pull request #427)
       b8244276 HiCSeq - updated mugqic_tools to 2.12.0 to use the fixed version of hicrep.R (for no down-sampling cases)
       2e15b7fa Merged in update_licence (pull request #422)
       4d8e1f59 RNASeq - fixed rnaseqc2 output dependencies
       fdfb9b43 Modules - Updated binary wrapping done by installation scripts
       369653aa Modules - Added bamixchecker installation script
       171a7445 Module - Updated mugqic_tools.sh with latest version
       a1d73e3c Utils - fixed watch_portal_folder for some cases with empty json : stop throwing error and remove the json instead
       4b8c7bde Resources - AGeNT : fixed typo in module
       dcdd1814 fix typo
       fcc76d83 GenPipes - improving job2json
       a458c390 GenPipes - strengthening the jsonator
       f869fdc4 Tumor Pair - refining resources assigned to gridss_paired_somatic jobs
       a6e9a3b1 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       a756c811 Tumor Pair - Fixed inputs for "sequenza.create_seqz" jobs
       3d9b5c00 Tumor Pair - fix linx_plot output dependency to avoid useless restart of the step
       55ceb6d9 HiCSeq - RobusTAD : fixed output name to fix restart issue
       b1d7edca Utils - made mugqicValidator.py Python3-compatible
       59f2a414 Tumor Pair - fixed gatk4.ini for cobalt, amber and purple
       5ed67ce2 Tumor Pair - fixed gripss.py arguments in bfx
       6d511de6 Tumor Pair - fixed settings for mutect2 when using gatk4
       1a0c1066 Core - Add details to some sanity check error messages
       2419d273 RNASeq - fixed stringtie-merge dependencies
       47cf66ad Merged in tumor_pair_dependency_fix (pull request #413)
       acd69760 Tumor Pair - Fixed gripss calls
       6793fdf9 Version bump to 4.4.2-beta
       f40eae05 Merged in release_4.4.1 (pull request #412)
       e9cd38bb Version bump to 4.4.1
       0d2e78fc Merged in release_4.4.1 (pull request #411)
       412c2e3d GenPipes - in prep for bug-fix release 4.4.1
       c27edf15 Resources - update mosdepth version to 0.3.3 in install script
       6a5222f2 EpiQC - correcting chromimpute_preprocess command : do not use "rm" and force symlink creation with the "-f" flag
       087fd9b3 Resources - adding / updating install scripts
       252b6ee5 RNASeq Light - Added mkdir to the kalliso command
       07f45a23 Resources - adding 'gget' install script and updating others
       ed8e359e Version bump to 4.4.1-beta
       20c470b7 Merged in release_4.4.0 (pull request #410)
       8cabdc1c Version bump to 4.4.0
       1c294768 Merged in release_4.4.0 (pull request #409)
       bc806684 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       468975cb updating README
       b0db9b96 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       67fe5bbe RNASEq Light - Fixed Kallisto step
       9e014170 Resources - added HMMcopy package in defaults R libraries
       ffd223d0 RNASeq Light - Adding kallisto wrapper in bfx
       91bca69e RNASeq Light - revamping the call to Kallisto to have per gene counts instead of per readset
       d83071de Resources - Adding install script for hmmcopy-utils and ichorCNA
       79b74d03 Resources - Adding AGeNT install script and updating script for VarDict
       d2a8fef7 RNASeq Denovo Assembly - fixed trinity job dependencies
       9629f9a3 RNASeq Light - removing kallisto output folder before runnig kallisto to avoid mixed or bad results in case of pipeline restart
       eb7df8b1 EpiQC - fixed signal_to_noise dependencies
       48f936c1 DNASeq High Coverage - increased default walltime for gemini_annotations
       43121ea9 HiCSeq - fixing "homer_tag_directory" outputs
       4f7bfea4 Merge remote-tracking branch 'origin/dev' into release_4.4.0
       cf282891 GenPipes - updating READMEs before release
       70d27f7c Resources - Adding some new install scripts
       699683ca Resources - updated version in some install scripts
       3d0c96d9 EpiQC - improved restart mechanisms
       aa742e32 Merged in release_4.3.2 (pull request #389)
       ed7c1835 Merged in release_4.3.2 (pull request #388)
       118309c8 Merged in release_4.3.1 (pull request #383)
       baeee47d Merged in release_4.3.1 (pull request #382)
       f0feb9ad Merged in release_4.3.1 (pull request #381)
       e5980e36 Merged in release_4.3.0.1 (pull request #367)
       110f8ad7 Merged in release_4.3.0.1 (pull request #366)
       2e0a7ecc Merged in release_4.3.0 (pull request #364)
       390169f3 Merged in release_4.3.0 (pull request #363)
       47523b32 Merged in release_4.2.1 (pull request #359)
       9f80a8c9 Merged in release_4.2.0 (pull request #354)
       927e40e4 Merged in release_4.2.0 (pull request #353)
       a546428e Merged in release_4.1.3 (pull request #352)
       f1e6716e Merged in release_4.1.3 (pull request #351)
       a1be2b96 Merged in release_4.1.2 (pull request #314)
       e5a64db9 Merged in release_4.1.2 (pull request #313)
       99b40bf6 Merged in release_4.1.2 (pull request #312)
       dc6f51e6 Merged in release_4.1.2 (pull request #311)

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      14 commits

       43767b82 Genpipes : methylseq - clarifying multiqc comment
       df70dbfa GenPipes : methylseq - update mkdir commands
       9614d14f GenPipes : methylseq - add comments to explain multiqc options
       918d8ac3 bfx : dragen.py - add report_files to outputs
       1ce17811 bfx : dragen.py - extend outputs to avoid list within list
       2da1c94d GenPipes : methylseq - change automatic sample naming for trimmomatic multiqc
       60baae95 GenPipes : methylseq ini add multiqc options
       8961e31d GenPipes : methylseq - ihec table to multiqc format
       04d0de10 GenPipes : methylseq - add multiqc to bismark align
       3bdff148 bfx : make sambamba options ini section not required
       f2b9ccb0 GenPipes methylseq : add multiqc to pipeline
       37608f3d GenPipes methylseq : add symlink_dragen_metrics section to ini
       776f679f GenPipes methylseq : add multiqc to base ini
       97ba46ee bfx : add metrics files to expected dragen output and report files

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      66 commits

       485c806f Merged in release_4.4.2 (pull request #430)
       ed55718f Merge remote-tracking branch 'origin/dev' into release_4.4.2
       8888ad03 update pipeline READMEs prior to release
       b1c13898 Merged in methylseq_multiqc (pull request #429)
       15780ffb Merged in trimmomatic_multiqc (pull request #428)
       b2387dcc GenPipes : add symlink for trimmomatic log to use with multiqc
       5a068a98 GenPipes methylseq : increase mem for filter_snp_cpg step in cit
       9f728b78 Merged in rm_hard_clip (pull request #426)
       a0b21222 GenPipes tumor_pair : loop over report_files output for multiqc symlinks
       1850d23d bfx : add report_files output to picard metrics
       8726b5c3 GenPipes tumor_pair : multiqc - add sample name to qualimap symlinks
       271b96ec GenPipes : rm sequencing_center refs to wrong centre, add underscores to McGill_Genome_Centre
       7eb52d45 GenPipes rnaseq : rm bam_hard_clip step
       f40ed17b GenPipes : rm wrong sequencing centre, add underscores to McGill_Genome_Centre
       5c6fe1c3 Merged in cit_multiqc_fixes (pull request #424)
       621b0f6a Genpipes methylseq : increased mem for bissnip and filtering steps during cit
       4b739b5a GenPipes tumor_pair : cleaned up commented out and unnecessary lines
       03c3e8ac GenPipes epiqc : add mkdir to chromimpute convert step
       6e1c8747 Genpipes tumor_pair : fix path for is_gz file check
       68baa958 GenPipes tumor_pair : remove unnecessary os.mkdir call
       0eae1f99 Genpipes tumor_pair : create multiqc report dir as part of jobs
       859c1cc7 GenPipes tumor_pair : fix multiqc symlink creation for qualimap outputs
       8ee7b43d Merged in rnaseq_multiqc (pull request #423)
       21db3e9d GenPipes rnaseq : multiqc module order
       1128cb75 GenPipes rnaseq : fix relative path for input of rseqc.tin step
       477a6893 bfx : rm abspath in input for rseqc.tin
       cd04bf45 GenPipes rnaseq : fixed interactive option for multiqc in ini
       f448b924 GenPipes rnaseq : fixed symlink in run_arriba, change dependency exit status for multiqc
       916d36de GenPipes rnaseq : fixed symlink in run_arriba
       efd63a17 GenPipes rnaseq : merged and resolved conflict with dev
       13ab19dd GenPipes rnaseq : fix symlink relative path
       8fb9e184 GenPipes rnaseq : make multiqc interactive
       ba66e8b9 GenPipes rnaseq : fixed symlink names for multiqc inputs
       70db4203 GenPipes rnaseq : removed stray commas, fixed typos
       30ce0929 GenPipes rnaseq : removed stray comma
       c56ce402 GenPipes rnaseq : fix job name
       e851d87b GenPipes rnaseq : fixed another relative path command
       c02906b4 GenPipes rnaseq : added missing link directory assignment
       b1796d83 GenPipes rnaseq : fixed relative path command
       451dbf41 GenPipes rnaseq : added missing parenthesis
       d92e2603 GenPipes rnaseq : added missing comma
       eae40ef6 GenPipes rnaseq : add multiqc section to base ini
       97a255fd GenPipes rnaseq : add multiqc to pipeline
       9481206e bfx : add log files to expected output for multiqc integration
       baa0813b Merged in tp_ensemble_fix (pull request #421)
       72b3853f GenPipes tumor_pair : add remove command before bcbio_ensembl call to avoid silent error when output exists
       70b3e161 Merged in chipseq_blacklist_filter (pull request #418)
       63d8fe7a Chipseq : add blacklist path to ini
       02a2c734 Merged in readme_updates (pull request #417)
       1d18f5e3 GenPipes chipseq : update README to include blacklist filtering step
       6ac84851 GenPipes chipseq : added input selection to cram, variant calling steps
       5e3f6139 GenPipes chipseq : fix inputs for diffbind and run_spp depending on blacklist removal step
       2c510bd8 bfx differential_binding : add option to specify alignment file extension
       b6cbc84a GenPipes chipseq : fixed various typos, add multiqc inputs
       86f54efe GenPipes chipseq : add blacklist and bedtools module to base ini
       30bc6ff0 GenPipes chipseq : fixed typo
       3660c64d GenPipes chipseq : add bedtools intersect step to pipeline, add ini section
       45093bbf GenPipes chipseq : add bedtools intersect to remove blacklist reads prior to peak calling
       b6fb9469 GenPipes : improved README for dnaseq pipeline
       f47ca3af Merged in Mareike-Janiak/commonpy-edited-online-with-bitbucket-1681414777818 (pull request #416)
       bc8ecef5 common.py edited online with Bitbucket - fix for sambamba merge dependency
       420681f8 Merged in sambamba_merge_fix (pull request #415)
       be8788f4 GenPipes common : simplified rm command before sambamba_merge
       d6b84fb5 GenPipes common : add comments to rm command before sambamba merge
       05957e7a Common : sambamba_merge_sam_files fix typos
       f45b8300 Common sambamba_merge_sam_files : remove any existing file/link before merge

  Pascale Marquis <pascale.marquis2@mcgill.ca>      1 commits

       88ebed7a Merged in new_branch_pascale (pull request #420)

  Pascale Marquis <pmarquis@Weigela.local>      1 commits

       e622624b changed the cpulimit

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       cf3d7277 Merged in release_4.2.1 (pull request #358)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       c7b3aa6d Merge remote-tracking branch 'origin/dev' into release_4.2.1

DOvEE-4.3.2        Thu Jun 15 18:35:19 2023 +0000        87 commits

  ehenrion <edouard.henrion@mcgill.ca>      2 commits

       5dc6ac85 Merged in dovee_multiqc (pull request #414)
       24f8c8f0 DOvEE - set multiqc_inputs as a pipeline attribute and make changes accordingly

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      75 commits

       a8d140f9 Genpipes dovee : update main README
       386226da GenPipes dovee : add README
       03850c76 GenPipes dovee : cleaning up old comments
       dc28a034 GenPipes dovee : update filepaths in ini
       1cc87e8a cleaning up
       01a9d91a bfx : add missing output file to locatit
       cf82e8de DOvEE : added parsing of locatit and creak output metrics, added to multiqc
       22a18b43 DOvEE : changed options for multiqc and bcftools_stats in ini
       6d46406d DOvEE : add parsing of creak metrics, picard gcbias metrics to pipeline
       616abcae DOvEE : add picard collect gcbias metrics to ini, add multiqc options
       7e156953 bfx picard2 : fixed typo in collect_gcbias_metrics - qcbias = gcbias
       e243f8fa bfx agent : add stats file to expected creak outputs
       67161bd9 DOvee : add creak and subsample sections to ini
       9f0a8e76 dovee : add creak dedup and samtools subsample steps
       caee1e3c bfx : add creak to agent wrapper
       c0ee6560 dovee : add sample source to symlink names for multiqc display
       0caca75e dovee : change mosdepth bed file used
       7a5287a2 DOvEE : add bam matching step with conpair to that samples are from same patient
       40d46af7 Merge branch 'DOvEE_genpipes_Sflag' into DOvEE_genpipes
       f18819e1 DOvEE : modified bwa mem settings in ini
       32f21dd3 DOvEE : add option to picard HS metrics to ignore validation errors
       732ebaef DOvEE : use sorted locatit output to side step sort error
       86218505 DOvEE : update multiqc to look for symlink inputs in one directory, set order of modules in report
       d04276f2 DOvEE : change multiqc job dependencies to run even if jobs fail
       8c6cecd6 DOvEE : multiqc inputs changed to directories to avoid error with failed samples
       3adb81cc DOvEE : fix error in trimmed read symlink creation
       829ff3ca DOvEE : add --interactive option to multiqc to prevent flat plots
       7dbc1405 DOvEE : fix bed files used for mosdepth
       7cf8ad0b DOvEE : add readset level directories to mapping step
       20a7fdd9 DOvEE : add tmp directory to samtools sort
       e9dc8019 DOvEE : make files executable
       7e2daa44 DOvEE : remove -K flag from locatit call to avoid keeping tmp files
       7f32db9d DOvEE vardict : add bcftools stats step
       210429f6 DOvEE : add bcftools stats to vardict protocol
       d4401227 cleaning up, fixing indents, copyrights
       eac8909c dovee : cleaning up and adding comments
       54cd4c4a adding copyright sections
       01913ce8 dovee : added option to distinguish between brush and saliva via modified version of design file
       841ae379 small dovee ini fixes
       adb834e6 add output file to mosdepth
       ecd7b09d created dovee sample pair parser based on tumor pair system
       b9b8f34f adding multiqc
       eed14380 dovee samtools sort : changed to avoid restart errors
       bd3c654f dovee metrics : add picard HS metrics
       a3984061 dovee : added metrics steps
       8764ebbc mosdepth : fixed typo
       4f0edd44 trimmer : edits to STATS rm command
       2a211153 trimmer : changed mv to ln to prevent job restarts
       dda4c1b5 dovee : adjust resources
       4f6a4283 dovee ini : ichor version change, other small fixes
       5e3db81e ichorCNA : fixed typo
       83338513 dovee : merge sort & index steps, update sort to look for inputs from previous jobs, add hybrid dedup for copy-number protocol
       cebe94f3 dovee ini : add/fix filepaths
       e3a2fe4a vardict : add module call for vardict_single
       74431a8f dovee small fixes
       6579f709 dovee resource and other ini updates
       95eb5438 samtools merge : change module call and add other_options
       eed26dd5 hmm : add module call
       cd9d0eac tweaked job names to be consistent with ini
       4cea397d updating resources and filepaths
       0d1da041 dovee : add bam merge to allow for top ups
       42146308 ichorCNA : removed # from module call
       dc9aa5d6 bwa mem : remove single end option
       2e60edbe update bwa options, update agent calls
       de382e3b combine trimmer and locatit into single agent script
       cc534520 trimmer and locatit updates now that modules are installed
       3935a130 DOvEE : add placeholder wig file paths
       08a0c453 DOvEE : fixing small errors
       7c697e2a DOvEE : add init
       7bb7402f fixed small errors and indentation throughout
       47894372 fixed indentation
       3b0fe908 moved comment
       c7f73451 fixes in pipeline and ini for new dovee pipeline
       f9495d5e new trimmer trimming function added
       912643ea new locatit dedup function added

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      9 commits

       3c7f0b42 adding fastp.py from run_processing branch to dovee
       6b549b30 mosdepth : add to genpipes
       9f2ae202 trimmer : remove STATS file if already exists
       84a8ff7b DOvEE : add tumor pair pairing system
       27a4d74c DOvEE : add pairs file option
       9db0a250 GenPipes DOvEE : add copy number protocol to pipeline
       921ca329 GenPipes DOvEE : add hmm and ichorCNA sections to ini
       41bb80d8 GenPipes DOvEE : create script for ichorCNA
       f8fdbb48 GenPipes DOvEE : create script for hmm readCounter

4.4.1        Tue Mar 14 17:52:02 2023 +0000        10 commits

  ehenrion <edouard.henrion@mcgill.ca>      10 commits

       0d2e78fc Merged in release_4.4.1 (pull request #411)
       412c2e3d GenPipes - in prep for bug-fix release 4.4.1
       c27edf15 Resources - update mosdepth version to 0.3.3 in install script
       6a5222f2 EpiQC - correcting chromimpute_preprocess command : do not use "rm" and force symlink creation with the "-f" flag
       087fd9b3 Resources - adding / updating install scripts
       252b6ee5 RNASeq Light - Added mkdir to the kalliso command
       07f45a23 Resources - adding 'gget' install script and updating others
       ed8e359e Version bump to 4.4.1-beta
       20c470b7 Merged in release_4.4.0 (pull request #410)
       8cabdc1c Version bump to 4.4.0

4.4.0        Thu Mar 9 18:24:15 2023 +0000        156 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      24 commits

       afb8ae79 RNASeq - Fix arriba inputs
       18fb0411 RNASeq - fixing Arriba inputs and dependencies
       f4913ea1 RNASeq - Set the default build back to GRCh37
       658004c6 EpiQC - Fixing `chromimpute_convert` dependency with `chromimpute_preprocess`
       222c9d9f GenPipes - Resources : updated java and nextflow installation scripts
       be0334b9 GenPipes - Resources : some fixes in install_genome.sh
       77071b2a GenPipes - Linting + multiqc dependency stengthening
       efeb9069 Minor update : just removing trailing comma
       039fef26 GenPipes - MultiQC : updated all the multiqc call with latest MultiQC version and removed the loading of Python module
       28943fc6 HiCSeq - corrected typo
       8095f0fe RNASeq - Adding 'cancer' and 'variants' protocols + removing 'cufflinks' protocol
       0ba1053e BFX - Linx : minor formating updates
       d33db131 Tumor Pair - SV : correcting purple.py after rebase
       7e462a31 Tumor Pair - SV : revamping SV protocol of Tumor Pair with use of amber, purple, gridss, etc...
       5efad869 BFX - Conpair : corrected rm command
       d9b76059 Tumor Pair - fixed typo in filename
       7074e100 Tumor Pair - fixes qualimap inputs to multiqc
       483a1c3b Tumor Pair - updating multiqc to 1.14
       26ac4f1c Tumor Pair - update inputs to multiqc
       aad32e3c Tumor Pair - Strengthening the depedencies for the multiqc reports
       23212900 BFX - MultiQC : removed python module from command
       ada37b0c Tumor Pair - corrected multiqc dependencies
       ec003878 Tumor PAir - Adding Compair and Purple to multiqc call
       fe29ad4c Version bump to 4.3.2

  ehenrion <edouard.henrion@mcgill.ca>      49 commits

       1c294768 Merged in release_4.4.0 (pull request #409)
       bc806684 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       468975cb updating README
       b0db9b96 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into release_4.4.0
       67fe5bbe RNASEq Light - Fixed Kallisto step
       9e014170 Resources - added HMMcopy package in defaults R libraries
       ffd223d0 RNASeq Light - Adding kallisto wrapper in bfx
       91bca69e RNASeq Light - revamping the call to Kallisto to have per gene counts instead of per readset
       d83071de Resources - Adding install script for hmmcopy-utils and ichorCNA
       79b74d03 Resources - Adding AGeNT install script and updating script for VarDict
       d2a8fef7 RNASeq Denovo Assembly - fixed trinity job dependencies
       9629f9a3 RNASeq Light - removing kallisto output folder before runnig kallisto to avoid mixed or bad results in case of pipeline restart
       eb7df8b1 EpiQC - fixed signal_to_noise dependencies
       48f936c1 DNASeq High Coverage - increased default walltime for gemini_annotations
       43121ea9 HiCSeq - fixing "homer_tag_directory" outputs
       4f7bfea4 Merge remote-tracking branch 'origin/dev' into release_4.4.0
       cf282891 GenPipes - updating READMEs before release
       70d27f7c Resources - Adding some new install scripts
       699683ca Resources - updated version in some install scripts
       3d0c96d9 EpiQC - improved restart mechanisms
       dc3822fb RNASeq - Arriba : updated path of 'blacklist', 'known_fusions' and 'protein_domains' annotation files for GRCh37 in rnaseq.base.ini
       1eef454c Merged dev into cit_fixes_mcj
       925fa891 Nanopore CoVSeq - minor linting update
       6f4b9a9b EpiQC - fixed file path : distinction between files for job dependency and files for direct access
       d69696e6 GenPipes - resources : added ceratits_capitata.EGII-3.2 install script + version update in java.sh and picard.sh
       f9ca3c1c Tumor_pair - corrected java19 module
       07f2bfde PCGR report - added 'ls' command on the html report output to confirm job success
       991a65db Fixed mark_duplicate calls in all pipelines being from picard or gatk + Fixed pcgr report argument order in tumor_pair
       b6a29235 GenPipes - RNASeq : correctied calls to gatk_mark_dupliates + fixed arriba annotations for b38
       691ba5b6 GenPipes - Tumor Pair : updating cit.ini for GRCh37 default
       582055e4 GenPiopes - RNASeq : fixed picard/gatk_mark_duplicate
       362b7243 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       c5c436a0 RNASeq - Stringtie : updating resources for picrad mark duplicates
       0cec83cd RNASeq - Fixing sym link in star_align
       7faa99ae GenPipes - resources : update of GRCh38 install scripts
       9edf7b59 Resources - modules : adding pandas to the default python libraries
       bccb7031 Merged in tpsv (pull request #401)
       598e4c3f Resources - version updates for bcftools, htslib, samtools, python3 and multiqc
       fc5d59ff Adding 'core/__pycache/' to .gitignore
       317283ef EpiQC - minor change in chromimpute
       834752e2 Resources - Genome : updated 'create_bwa_index' with sym links to .fai and .dict
       aad898f3 Merged in tp_fix_multiqc (pull request #399)
       74f7b451 Merged dev into tp_fix_multiqc
       6017946a RNASEQ - Wiggle : corrected deeptools call
       1f444bca Merged dev into Gerardo_rnaseq_bigwig
       fd2b64af Merged in hotfix_dev (pull request #391)
       f64dc2b0 DNASEQ - ensemble_metasv dependency fix
       94efcdb7 GenPipes - Resources : Adding sv-prep installation script + updating gridss, gripss, linx, and purple scripts
       aa742e32 Merged in release_4.3.2 (pull request #389)

  Gerardo Zapata Abogado <gerardo_za_94@hotmail.com>      52 commits

       7b7706d4 Merged in Gerardo_rnaseq_bigwig (pull request #390)
       ba3ce5c7 rnaseq_denovo [in silico] small fix of the base.ini
       36d1ac84 rnaseq.py [wiggle] _32 addition of strand2
       841e7af7 rnaseq.py [wiggle] _30 addition of strand
       2003ed75 rnaseq.py [wiggle] _30 bamcoverage to wiggle
       e7f2e491 rnaseq.py [wiggle] _29 added the addition of merged + forward + reverse - missing "s"
       354ffb02 rnaseq.py [wiggle] _29 added the addition of merged + forward + reverse
       99a979a4 readed deeptools.py
       eb2834dc removed Deeptool.py
       4bf5befa rnaseq.py [wiggle] _28 small changes to .ini and .py
       28750083 rnaseq.py [wiggle] _27 small changes to .ini and .py
       36a123e2 rnaseq.py [wiggle] _26 small changes to .ini and .py
       a09aabd2 rnaseq.py [wiggle] _24 small changes to .ini and .py
       24279a3b rnaseq.py [wiggle] _24 small changes to .ini and .py
       12f3d6a7 rnaseq.py [wiggle] _22 small changes to .ini and .py
       39c12757 rnaseq.py [wiggle] _22 small changes to .ini and .py
       296ed2b9 rnaseq.py [wiggle] _21 small changes to .ini and .py
       4602dab6 rnaseq.py [wiggle] _20 changes to add option for for & rev.
       8a149010 rnaseq.py [wiggle] _19 need to implement changes to add option for for & rev.
       b56bbc0c rnaseq.py [wiggle] _18 location of ifelse loop
       3834821f rnaseq.py [wiggle] _17 set up option to run separate .bw??
       ecf6046c rnaseq.py [wiggle] Deeptools _16 add forward and reverse options
       7265f011 rnaseq.py [wiggle] Deeptools & Ini changes _15 small changes
       e2f2d558 rnaseq.py [wiggle] Deeptools - _14 string errors
       b46b26b1 rnaseq.py [raw_countMatrix] -  _12 for loop removal
       56ee8742 rnaseq.py [raw_countMatrix] -  _11 bigwigs .zip
       d53d4f4b rnaseq.py [wiggle] - Deeptools _11 "missing comma" :|
       112882f7 rnaseq.py [wiggle] - Deeptools _11 "output_location error" + "archive"
       6cdd381a rnaseq.py [wiggle] - Deeptools _10 "name" "samplename" location
       4877bab6 rnaseq.py [wiggle] - Deeptools _9 "name" "samplename" location
       340669c2 rnaseq.py [wiggle] - Deeptools _8 missing comma
       35c25218 rnaseq.py [wiggle] - Deeptools _t "name" "samplename"
       218a9789 rnaseq.py [wiggle] - Deeptools _6 missing comma
       a361eb14 rnaseq.py [wiggle] - Deeptools _5 "name" "samples"
       89a3e478 rnaseq.py - [wiggle] - concate_jobs
       ca916a9d rnaseq.py - fix Deeptools.bamCoverage()
       26b58695 Deeptools.py - fix bamCoverage
       97733826 rnaseq.py [wiggle] - Deeptools
       08e56f63 rnaseq.py [wiggle] - fix "tracks_directory"_4
       f5b8293b rnaseq.py [wiggle] - fix "tracks_directory"_3
       c2b03d99 rnaseq.py [wiggle] - fix "tracks_directory"_2
       91073b16 rnaseq.py [wiggle] - fix "tracks_directory"
       ed2f1c48 ranseq.py - added new wiggle step, removed previous
       0ae71bb2 rnaseq.py Modify [wiggle] ini
       bb9b9120 Added - Deeptols.py
       fe825b14 rnaseq.py - readcounts - ".csv" to ".tsv" - ALL
       0bb883c0 rnaseq.py - rawCountMatrix - ".csv" to ".tsv" - ALL
       1b247225 rnaseq.py - readCounts - ".csv" to ".tsv" 2
       90bb9e29 rnaseq.py - readCounts - ".csv" to ".tsv"
       355a5191 rnaseq.py - rawCountMatrix - ".csv" to ".tsv"
       0ed06e09 test commit Gerardo
       2b424a13 test commit

  mareike.janiak@computationalgenomics.ca <mareike.janiak@computationalgenomics.ca>      4 commits

       2121b4f2 methylseq : do not require methylation_protocol and mapping_implementation parameters to be set for hybrid protocol
       791aeb88 GenPipes methylseq : make single-pass dragen align option compatible with bismark in hybrid protocol
       b0507cf3 GenPipes - methylseq : fix cp from dragen
       2d69b61f issue caused by qiime_catenate loading two python modules

  Mareike Janiak <mareike.janiak@computationalgenomics.ca>      23 commits

       33a2685f Merged in cit_fixes_mcj (pull request #408)
       9e6f6e80 Rnaseq cancer : update genome version used for cpsr and pcgr
       a48b0ed5 Merged dev into cit_fixes_mcj
       a9567608 Merged in cit_fixes_mcj (pull request #407)
       4cab1006 gatk_indel_realigner : set -fixMisencodedQuals flag only when QualityOffsets are 64 (solexa/illumina1.5+)
       0e459caa skewer trimming : adjust trim settings based on quality offset in readset file
       b28269b7 Merge branch 'cit_fixes_mcj' of bitbucket.org:mugqic/genpipes into cit_fixes_mcj merge cit updates with local changes
       575d3052 skewer trimming : adjust trimming settings based on quality offset provided in readset files
       88128229 Merged in EpiBrain_methylseq (pull request #406)
       2ab0e92b Merged in cit_fixes_mcj (pull request #405)
       1b9d47c7 Merged in EpiBrain_methylseq (pull request #404)
       ab58ff99 GenPipes - rnaseq : adding cd command to run_arriba
       5801352e removing outfile options from arriba
       34982275 Merged in cit_fixes_mcj (pull request #403)
       25f08c9b GenPipes arriba : added specific paths for arriba output files to run_arriba command
       33515ecd GenPipes - rnaseq : removed -fixMisencodedQuals from ini for gatk_indel_realigner
       1ffffc2e GenPipes - rnaseq : fixed file extension in ini for run_arriba
       83f0d22f GenPipes - epiqc : add mkdir to chromimpute_convert
       2c6d4c71 Merged in ampliconseq_abacus_mcj (pull request #398)
       ce9f2755 Merged in rnaseq_star_1pass (pull request #394)
       fae754ee Merged dev into rnaseq_star_1pass
       cbfc7659 Merged in tumor_pair_mcj (pull request #393)
       c9343cf6 Merged dev into rnaseq_star_1pass

  Mareike Janiak <mcj43@narval1.narval.calcul.quebec>      1 commits

       0522542f ability to skip 2-pass star alignment and only do 1 pass

  Mareike Janiak <mcj43@narval2.narval.calcul.quebec>      1 commits

       cd4a4556 gatk_indel_realigner : added two cd commands

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      1 commits

       38858436 test paul

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       638b81a5 Merged in rebase_rnaseq_variant (pull request #370)

4.3.2        Thu Dec 8 17:21:23 2022 +0000        132 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      54 commits

       ca99865e Merge remote-tracking branch 'origin/dev' into release_4.3.2
       118e25a3 GenPipes - Updating READMEs for release
       fdb289a0 Tumor Pair SV - properly setting the svaba indels vcf file in the config
       f568fe4d Tumor Pair - rearrangement of the ini for cit
       475e2b94 Tumor Pair SV - fixing metasv.ensemble call + fixed scones step for GRCh38
       eb7112e0 Tumor Pair SV - fixing step order
       16cd02b9 minor - replacing 'CentOS6' by 'root'
       53e3e05f Some more fixes induced by cit
       29e46972 RNASeq denovo asembly - fixed filter_annotated_components_exploratory
       c4d242b2 RNASeq denovo assembly - fixed differential expression steps
       77f4ca47 AmpliconSeq - enuring the use of Python2 for all the qiime commands
       143bc39a Tumor Pair Ensemble - fixed strelka dependencies to manta_indels
       ab569261 RNASeq - updated mugqic_tols to later dev version
       0aa6a1b2 DNASeq SV - updating cit.ini
       5c6358fc DNASeq - corrected typo in gatk.py
       ba6439b8 RNASeq denov assembly - updated mugqic_tools to later dev version
       73a4da52 AmpliconSeq - setting mugqic_tools to later dev version
       7227a940 AmpliconSeq - some more formating
       93f057f7 RNASeq - correction of rnaseqc and htseq_count after cit
       2fdb99b5 AmpliconSeq - correctig job names
       152103b7 DNASeq SV - fix breaseq2 call
       8c738322 RNASeq - test batch effect corecction in DESeq2
       06a8eeb7 RNASeq - fixed typo
       b367867f RNASeq - fix after cit
       74bc4202 HiCSeq - fix quality_scores afeter cit
       6ef0a219 DNASeq SV - fix svtyper after cit
       6f4cb2de pipelines/hicseq/hicseq.base.ini
       8264c71f DNASeq High Coverage - fix preprocess_vcf after cit
       7174d35f Ampliconseq - fix flash step after cit
       a008986e Tumor Pair - Fixing strlka dependencies to manta_sv
       d6094118 RNASEQ + RNASeq Denovo : fixed batch corretion for RNASeq and added it to RNASeq denovo
       40cf4033 Corrections after cit + formating
       185ed831 GenPipes - Fixes after cit on Béluga
       2bebb47e Fixes after cit + code standardization
       c9b9e2bb DNASeq - SV : updating metasv module to 0.5.5
       b7d0ce4d DNASeq - SV : re-using delly VCF in meta_sv ensemble analysis, depending on meta_sv version
       3004136d DNASeq - SV : corrected manta_sv call by removing delly_vcf
       7c044401  DNASeq - correcting typo
       efa9aa38 DNASeq - correcting typo
       341de84d DNASeq - generalize use of self.output_dirs
       6775aad4 DNASeq - specify self.output_dir when setting files and folders path
       0284d9aa DNAseq - SV : correcting ini sectino definition
       313dc002 BFX - correcting typo in `cpsr.py` and `htslib.py`
       224f3eaf DNASeq - bug fixes for SV protocol
       48c925f9 DNASeq High Coverage - update sambamba_merge_sam_files in config
       ed1b33a6 Tumor Pair - format standardization
       3c82ce6f Tumor Pair - minor fix
       b21669d5 Tumor Pair - bfx/pcgr.py updated to fit with new version and keep backward conpatibility.
       e830d582 Tumor Pair - Report : updated with new version of PCGR (1.0.3) for both 'pcgr' and 'cpsr' reports
       c024bd6b DNASeq - corrected merge_sam_files + extract unmapped steps
       8279cfa4 Core - Scheduler : corredcted PBS mem-per-cpu settings when --force_mem_per_cpu
       a34fda40 Core - Scheduler : corrected `memory` call in `cpu` from dev
       4a045d5a Core - Scheduler : corrected `memory` call in `cpu`
       80afe6df Version bump to 4.3.1

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      45 commits

       88e0ba82 CovSeq - Improving restart mechanisms
       a1ef9352 ChIP-Seq - improved restart mechanisms
       962c90d9 AmpliconSeq - Improving the restart mechanisms
       aea76ebf AmpliconSeq - Improving restart/resume mechanisms
       4e401e9f Tumor Pair SV - correcting svaba ref configuration
       d63c56d5 Tunmor Pair SV - CNV detection now done with cnvkit : removing SCoNes, not as well maintained as cnvkit.
       952693d8 Resources - corrected gridss install script
       3031bd63 DNASeq - returning to compressed filtered vcf for cnvkit
       705ace9a BFX - cnvkit : updated inputs in scatter and segment
       025165cc CORE - removing useless log.debug() in core/pipeline.py
       29e81cc3 DNASeq SV - testing with no compression during cnvkit_batch.vcf_flt
       cefc1dbb Tunmor Pair - Fix cit ini
       b411c0c9 RNASeq light - correcting R version for kalisto differential expression
       b5e1d2eb Tumor Pair - fixing cit ini
       5b440da9 Resources - put back DEV instead of wrongly modified GENFS
       ca574580 Modules - adding install script for tools needed in Tumor Pair SV
       7bc7fc39 DNASeq SV - corrected bcftools filtering before cnvkit + corrected metasv_ensemble command
       b67c4aaa DNASeq SV - improved bcftools options for vcf filtering before cnvkit_batch
       8536530c DNASeq SV - corrected filtered vcf outputed by cnvkit_batch.vcf_flt
       ae84cc42 DNASeq SV - improved filetering in cnvkit_batch.vcf_flt
       56f73fd1 DNASeq + Tumor Pair - fixed breakseq2 wrapper
       4dd168e0 Core - Pipeline : corrected log.error to log.debug
       f8eaf1dd Resources - updated mugqic_tools instalation script, adding java-tools in the module file
       be5efebd RNASeq denovo assembly - corrected outputs of edger and deseq2 jobs
       f6a7bb8d INI - updating mugqic_tools to latest production version for AmpliconSeq, RNASeq and RNASeq denovo assembly pipelines
       2b12dea5 DNASeq - correcting input dependencies for manta_sv jobs for better restart/resume of the pipeline
       c9493aff RNASeq - refined ballgown job outputs to improved pipeline resume/restart mechanisms
       6300baff RNASeq - imporoved restart/resume mechanisms
       a5df32f9 DNASeq + Tumor Pair - imporoved restart/resume mechanisms
       1c7a63f4 AmpliconSeq - fixes  on python modules and restart mechanisms
       9b21440d RNASeq denovo assembly - correcting arguments in calls of py_parseMergeCsv
       fc88f6a2 DNASeq - bug fix
       1417bd8e Tools - correcting typo in bfx/tools.py
       096206be Config - corrected pipelines/common_ini/Homo_sapiens.GRCh38.ini
       4deb294b GenPipes - Config : reorganisation of the ini files contained within the repository
       a7883c1a GenPipes - config : use common_ini for GRCh38 specific config
       e784b359 GenPipes - DNASeq & Tumor Pair : cleaning ini files
       723b201a DNASeq SV - cleaning resources/genomes/config/Homo_sapiens.GRCh38.ini a bit...
       0b1e36a3 Resources - Adding some genome installation scripts
       b316c215 RNASeq denovo assembly - Corre ting calls to py_parseMergeCsv
       cadad220 AmpliconSeq - Correcting python call in qiime_otu_table
       9f3f3e75 Resources - adding some more module installation scripts
       20c0951e Install script - updates in module install scripts
       3341cbc8 Tumor Pair - Improving the restart mecanisms
       df5aa276 GenPipes - DEV : updating CHANGELOG and setting beta version

  ehenrion <edouard.henrion@mcgill.ca>      17 commits

       ed7c1835 Merged in release_4.3.2 (pull request #388)
       23f764a4 GenPipes - Resources : updated mugqic_tools install script with latest version
       15091a00 DNASeq SV - updated pair_diretory path with non-absolute path
       b18b9d3a EPIQC - worked on file path to avoid depencencies errors and missing files
       5678259d HicSeq - stop uding mugqic_dev
       552c0fba RNASeq Light - fixing kallisto_count_matrix dependencies
       ee70b175 Tumor Pair SV - fix delly sv annotations dependencies
       71f95897 Tumor Pair SV - fixing SVannot call
       3e68d78a GenPipes - more formatting and linting
       90b95353 GenPipes - code cleaning and formating
       2dcd2e6e CoVSeQ - using report directory to output freebayes and ivar reports
       ef89a8f4 GenPipes - improved averall restart/resume mechanisims with relative path
       10727fed HicSeq - working on restart/resume mechanisms + doing some reformating
       ce5e36f2 Merged in tp_hotfix (pull request #385)
       03d0ab28 DNASeq - gnomad : corrected link to vcf
       118309c8 Merged in release_4.3.1 (pull request #383)
       baeee47d Merged in release_4.3.1 (pull request #382)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      16 commits

       0148c73f General - Linting
       dacf35a0 Merged in HotFix_dev (pull request #387)
       b9419c5b tumour_pair - artifacts_metrics - Fixing output file name for GenPipes to be able to skip step
       e799a375 covseq - sambamba_filtering - Adding SE options in ini file
       27b4f550 covseq - sambamba_filtering - Debug for SE mode
       cafd9749 rnaseq_light - base ini - Updating mugqic_R_packages
       7ee0df69 Merged in bwa_sambamba_splitting (pull request #386)
       0182fc0f dnaseq_high_coverage - bwa_mem_sambamba_sort_sam - Adding bwa_mem_sambamba_sort_sam in ini
       5c2dccb9 dnaseq/tumour_pair - sambamba_sort - Minor improvement
       8b57df09 dnaseq/tumour_pair - sambamba_sort - deleting bai if existing and chmod the output bai
       1050404c dnaseq/tumour_pair - sambamba_sort - Removing index from sorting as the bai is generated by sorting
       506605f4 Merged in bwa_sambamba_splitting (pull request #384)
       656678db dnaseq/tumour_pair - base ini - Cleaning commented out code
       dd4ca690 dnaseq/tumour_pair - sambamba_sort_index - Fixing index input
       bef6394e dnaseq/tumour_pair - alignment - Switching in protocols
       2fceb104 dnaseq/tumour_pair - alignment - First commit

4.3.1        Tue Oct 4 14:11:09 2022 +0000        74 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       a647b011 Merge remote-tracking branch 'origin/dev' into release_4.3.1
       759bbfb9 GenPipes : updating READMEs for release
       24107718 GenPipes - PBS Scheduler : stop  using -mem/-pmem parameters anymore in PBS job submission commands

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       f0feb9ad Merged in release_4.3.1 (pull request #381)
       94d7a6cb GenPipes - ChIPSeq : Updated mugqic_tools to 2.10.10 in chipseq.base.ini to prevent crash at differential_expression step
       6da854f3 GenPipes - ChIPSeq : updated mugqic_tools to 2.10.9 in chipseq.base.ini
       db7a4485 tumor_pair.extras.ini edited online with Bitbucket
       285568f9 GenPipes - README : replacing former "CentOS6" by "root" in cvmfs path and variables
       4d824c56 Merged in ehenrion/tumor_pairpy-edited-online-with-bitbucke-1655753579034 (pull request #369)
       f322ebd7 tumor_pair.py edited online with Bitbucket

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      1 commits

       c4a59f50 Changed the core pipeline script to modify the trace config ini file to include the timestamp. Additionally, the full command is now added to the ini header as a comment for future reference.

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      7 commits

       468bbfec Merged in config_timestamp (pull request #374)
       ed5738de Merged in covseq_nanopore_GPUqueue_bugfix (pull request #373)
       e1a91121 nanopore_covseq.base.ini edited to re-add the GPU_QUEUE value. Without it the pipeline crashes on Abacus.
       bd66f86d Merged in MetaSV-Delly-bug (pull request #372)
       693ee1d9 Removed Delly input file requirement from dependency list in metasv.py
       bf99d700 Removed delly referneces from the metasv step in the tumor pair pipeline.
       4183eb1f MetaSV has a delly argument that should not be there since it does not support Delly (see https://github.com/bioinform/metasv/issues/110). I removed the delly parameter.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      27 commits

       e2f0de69 Merged in HotFix_dev (pull request #378)
       56b2fcc5 rnaseq_light - base ini - Improving syntax
       511f59a2 methylseq - General - updating way to mkdir via bfx
       dec5d506 methylseq - filter_snp_cpg - Debug
       908b158f methylseq - filter_snp_cpg - Allowing to change cpg threshold for cit minimal
       1b708649 methylseq - methylkit - Changing the way Rscript is called
       5d2395f5 rnaseq - gq_seq_utils_exploratory_analysis_rnaseq - Upgrading to latest R_packages + fixes
       210a4737 rnaseq - gq_seq_utils_exploratory_analysis_rnaseq - Upgrading to latest R_packages
       2fe1ac6e chipseq - differential_binding - Adding parameter contrastnb for cit atacseq only (set it to "cit" for cit minimal)
       794a0663 hicseq - homer_tag_directory - Adding chrList from ini file
       1619b24b hicseq - quasar_qc - Allowing more than 1 "_" in chr name
       9117fe42 hicseq - General - Upgrading mugqic_tools and trimmomatic
       5269243a hicseq - identify_TADs_RobusTAD - Adding min and max window parameters
       cf641882 hicseq - multiqc - Upgrading to latest multiqc version
       b1dead1d hicseq - hicup - Moving remove before hicup conf
       0a72d80f hicseq - hicup - remove recursive
       28b16290 hicseq - hicup - Re-adding removing of results otherwise the step is not able to restart
       d837f8da hicseq - General - Improving the way genome assembly name is handled
       6e2d1aa0 hicseq - hic - Improving juicer call
       9e539f6b hicseq - hicup - Debug
       0c1b223c hicseq - hicrep - Code cleaning
       3018b34f hicseq - hicup - No longer creating folder and removing things in the step
       3ccf8078 hicseq - hicup - Debug
       42245c79 hicseq - General - Improving hicup job and removing path hardcoded in ini
       2424c6aa hicseq - General - Upgrading default genome to hg38
       c58d7953 epiqc - General - Improving code
       e6776b45 epiqc - General - Improving code

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       eb2067d4 Merged in force_mem_per_cpu (pull request #380)
       9822d934 Merged in fail_on_log_pattern (pull request #376)
       e1eeb177 Merged in sample_readset (pull request #338)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       c3e418bf support G and M only
       13c7c51f Remove explicit request for memory in pbs call when --force_mem_per_cpu option is used
       a74d9547 add a fail_on_pattern ini option
       4b796dea mv sample_tumor_pairs.py to core
       92049a77 move sample readset and design to core from bfx

  pubudumanoj <pubudumanoj@gmail.com>      2 commits

       d2990b24 modified inputinfofile and chr_sizes path
       847d4765 modified inputinfofile and chr_sizes path

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      16 commits

       c6bcd423 remove comments in tumor pair sv
       7a85f1d8 fix more bugs in tumor pair sv protocol
       2f76d8fc fixed more bugs
       ca6c2826 fixed major issues in sv protocol
       3908f5e0 modified the methylseq/cit.ini
       663b11f4 fixed a regression
       3166cdf1 chnaged resources in cit.ini
       d5aa2387 changed mugqic_tools version to 2.10.9
       ba295645 fixing issues
       a93935fc added methylkit differential binding back to the analysis
       830b3175 modified base.ini
       6e7dd92d Merge branch 'epiqc_code_improvement' of bitbucket.org:mugqic/genpipes into epiqc_code_improvement
       0c1426bc Corrected the inputinfo file paths in ini
       f787ac78 expanded env_var in inputinfo file path and chr_sizes path
       4e7e8e7a Corrected the inputinfo file paths in ini
       191ed379 expanded env_var in inputinfo file path and chr_sizes path

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       12be5b9e Merged in tumor_pair_hotfix_metsv (pull request #375)
       ceadb4e8 Merged in methylseq_cit_update (pull request #371)
       f2678bed Merged in epiqc_code_improvement (pull request #368)

4.3.0        Wed Jun 15 15:44:27 2022 +0000        129 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      25 commits

       57667889 Version bump to 4.3.0 - for real
       937a5ca4 Version bump to 4.3.0
       69d74b05 Merge remote-tracking branch 'origin/dev' into release_4.3.0.1
       ee071fd8 Updating READMEs
       3f33a35a Version bump to 4.3.0
       76b0b626 Merge remote-tracking branch 'origin/dev' into release_4.3.0
       ec90b7ab GenPipes - prep for new release
       86367757 Updating CHANGELOG from master after release
       1e06b04a Version bump to 4.2.1
       4e907416 GenPipes - DNASeq : fixing dependencies in recalibratino step for unmapped_reads
       46042e1c GenPipes - remove useless sections in ini
       40d6fe2f Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7dff261f GenPipes - move 'sambamba_merge_sam_files' into common.py + renamed 'sambamba_merge_sam_files' to 'sambamba_merge_sam_extract_unmapped' in DNASeq
       34ea28a1 Cleaning bfx/sambamba from samtools dependencies
       adf5e7ff GenPipes - DNASeq : revisited recal jobs to avoid useless restart of the step
       ed49153b GenPipes - DNASeq : Fixing dependencies for recal .bam.bai
       b7ef560e GenPipes - DNASeq : Recreating of the index of the recal bam file after merging of the unmapped reads
       373fb134 GenPipes - DNASeq : correcting the job which merges the unmapped reads to the recalibrated fastqs
       d30127dd GenPipes - DNASeq - Fix job queueing for merge unmapped to recalibrated
       1fc52437 fix typo
       afffee2b GenPipes - DNASeq : corrected unmapped_reads job queueing
       379dfdfb GenPipes - DNASeq : removing useless part of code...
       4b897f01 GenPipes - DNASeq : extracting unmapped reads from the merged bam, then re-instering them into the recalibrated bam
       460b6a3c Version bump to 4.2.0
       5cc640ce Version bump to 4.1.4-beta

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       5690df40 GenPipes - C3G Software Stack : showing only the 3 top versions of the softwares and adding '...' when necessary to indicate more older versions are also available
       34b06e68 GenPipes - C3G Software Stack : reverse-sorting the software versions when combining the metadata JSONs
       416fe4a4 GenPipes - C3G Software Stack : improved sorting of JSON output
       864465e6 GenPipes - C3G Software Stack : sort JSON output
       d3f03f6d GenPipes - C3G Software Stack : parametrized the redirection of the output to a file or to stdout
       701298b1 GenPipes - C3G Software Stack : removed the sending of the JSON
       85b70091 GenPipes - C3G Software Stack : corrected typo
       833d3caf GenPipes - C3G Software Stack : corrected message when file is sent
       3303a59f GenPipes - C3G Software Stack : Added sending of the combined JSON to url + reformating the module helpers
       0b8b7d59 GenPipes - module helper : updating the pythonSearcher class

  ehenrion <edouard.henrion@mcgill.ca>      11 commits

       e5980e36 Merged in release_4.3.0.1 (pull request #367)
       110f8ad7 Merged in release_4.3.0.1 (pull request #366)
       2e0a7ecc Merged in release_4.3.0 (pull request #364)
       390169f3 Merged in release_4.3.0 (pull request #363)
       47523b32 Merged in release_4.2.1 (pull request #359)
       5f690ffe Merged in soft_jsondb_gsoc2020_eh (pull request #181)
       ba84d8c1 dnaseq.base.ini edited online with Bitbucket
       943e103c Merged in unmapped_reads (pull request #355)
       ef817ea6 fixed typo
       ecdcb0ae dnaseq.base.ini edited online with Bitbucket
       9f80a8c9 Merged in release_4.2.0 (pull request #354)

  Moonshroom <yatharthrai16@ducic.ac.in>      4 commits

       a51bba0f Some documentation edits, README updates, fixed bad variable naming
       8d223cab Migration to argparse
       c6bde0f2 Added -h/--help flag
       810321f0 1. Pushing pre-final code 2. Added CLI args 3. Added CLI Documentation

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      23 commits

       ce997355 Merged in dragen_update1 (pull request #365)
       0c39884c rnaseq - General - Updating doc.
       30799528 General - General - Formalizing README doc.
       cd2845a8 Merged in HotFix_dev (pull request #357)
       a5caedfd General - General - Fixing UnboundLocalError: local variable 'job_name_prefix' referenced before assignment
       100567c6 Merged in HotFix_dev (pull request #356)
       a35809f5 General - General - Fixing typo
       edc01ae2 General - General - Fixing typo
       20801c45 General - General - Adding step_wrapper to slurm and batch mode (fixing regression)
       e7869220 General - general - Standardizing inis + fixing methylseq.base.ini
       343ec56d Merged in HotFix_dev (pull request #350)
       401be790 ampliconseq - General - Upgrading to latest release of mugqic_tools
       a52c6d62 covseq - General - Upgrading to latest release of covseq_tools
       92e0a5f7 ampliconseq - flash - Debug previous flash commit
       d73c35df ampliconseq - ampliconLengthParser - Fixing head with pipefail
       c74fe1ca ampliconseq - flash - Fixing None appearing in cmd
       720962b8 covseq - qualimap - Fixing qualimap parameter for threads
       92081f8e covseq - prepare_report - Fixing pipefail issue for ncovtools: we always want ncvotools error to be ignored and continue the job
       e40e7cc1 covseq - prepare_report - Fixing pipefail issue when grep doesn't find anything
       2c6b46c5 covseq - prepare_table - Adding threads param
       5387d298 covseq - prepare_table - Adding samtools as module required
       b8cfc5ea methylseq - metrics - Changing ini path of file
       f20a986a General - general - Standardizing inis + fixing methylseq.base.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       cf3d7277 Merged in release_4.2.1 (pull request #358)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       c7b3aa6d Merge remote-tracking branch 'origin/dev' into release_4.2.1
       b0920027 Update readme and version for release
       5469a903 remove -nt completely from indel realigner
       ccf85290 abacus ini typo
       720fd57f typo in tabix option, and no threads in gatk_indel_realigner
       68183a84 always force overwrite on tabix
       ef425e7e tweak and fix base ini

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      39 commits

       79ef5180 addressed comments from Paul
       2011bfb5 removed methylseqraw class
       a726b977 Fix issues in input and output dependencies
       66e55e3b resolved differences in dev branch files
       9c8dd415 Modified dragen.ini file
       c5acaceb fixed documentation typo
       5949d7bc fixed issues in output file tracking added documentation
       fea7a92c fixed issue in dragen hybrid protocol in single-pass
       8fca6d1c added dragen_bedgraph step
       454abe68 added split_dragen_methylation_report
       313ce780 fixed issues in dragen_methylation_call step working pipeline
       e2a952e2 trying to fix dependencies
       0c481543 added missing job to the dragen_aling
       15dda516 fixed a bug
       d3e98749 Add changes in common and methylseq
       c8531899 Implement jsonator disabling at job level
       4d64d734 Debug dragen protocol
       d1c7bf2a update dragen protocol
       8ed0ab89 Implement Modified dependency of dragen command to allo mv commands Added dragen function
       2ce9490b Implement Modified dependency of dragen command to allo mv commands Added dragen function
       dfa349d6 fixing issues with the done file tracking system
       86120112 Trying to fix dependency issue in dragen
       61945dc1 resolved merge conflicts
       f22f0f66 resolved merge conflicts
       6664a857 added dragen mark_duplication
       427d7a02 fixed a bug
       33caf904 update dragen protocol
       3cb1965f added missing job to the dragen_aling
       e0e21a72 Resolved incorrect tool link
       93a084b3 fixed a bug
       30ef76c4 fixed a bug
       fa534de9 Add changes in common and methylseq
       e8833f50 Removing debug outputdir log info
       8485ae7d Implement jsonator disabling at job level
       98cc9725 Debug dragen protocol
       5ec4bfc8 update dragen protocol
       c32117b8 pipelines/methylseq/methylseq.dragen.ini
       5e71366e Implement Modified dependency of dragen command to allo mv commands Added dragen function
       179cabf8 Implement Modified dependency of dragen command to allo mv commands Modify dependency of bash commands Added bfx/dragen.py

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       486907ab modified readme file
       be1755ab modified readme file with dragen update
       d82bab01 added comments and fix hybrid mode single pass directional-complement
       865a65f5 fixed ihec metric by adding a job to create an empty metric file for estimated_library_size
       cd94f556 changed dragen structure and added a new protocol

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       9684c5e2 fixed some issues in epiQC md file
       b772f84a Added multiline comments to methylseq pipeline

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       42112207 Merged in dragen_update2 (pull request #362)
       7faffc51 Merged in dragen_update_version2 (pull request #346)

4.2.0        Wed Jun 1 20:01:30 2022 +0000        230 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      16 commits

       695798b1 Generating READMEs for 4.2.0
       9947d057 Version bump to 4.1.3
       d7f19424 Merge remote-tracking branch 'origin/dev' into release_4.1.3
       e12f9ed2 in prep for a release
       e135fce3 GenPipes - Tumor Pair : corrected location of the coverage_bed in strelka jobs + prepend strelka job with 'rm' command to ensure folder is clean
       b7ed6142 GenPipes - Tumor Pair : some code cleaning before correcting strelka errors
       1d34130e GenPipes - DNASseq : updating base.ini
       20453e2b GenPipes - DNASEQ : fixing base.ini for SV protocol
       fe3ccfcd GenPipes - Scheduler : Re-using the default python, from ini, to call job2json.py
       61218acd GenPipes - DNASeq SV : fix sambamba_merge_realigned
       fa7cc656 GenPipes - DNASeq SV : using Python2 for breakseq2
       b790c913 GenPipes - BFX : CNVkit now uses the muqgic/CNVkit module, instead of a pyhton module which would
       9db4af25 GenPipes - DNASeq SV : set use of Python2 for MantaSV in base.ini
       97513bd2 GenPipes - DNASeq SV : fixing resource dependencies
       5685c7b8 GenPipes - DNASeq SV : fixed delly call
       1ef5a9ec Version bump to 4.1.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      4 commits

       da1d31b4 Resources - install scripts : update to mugqic_tools 2.10.4 + install_modules updates for improved installation in dev space
       5616bb75 GenPipes - Resources : ensuring the use of rpath when patching binaries
       b15769cb GenPipes - Resources : adding some new installation scripts
       802c02d2 GenPipes - Resources : Updating some install scripts with newer version + adeind bioinfokit and lofreq install scripts

  Édouard Henrion <henrione@narval1.narval.calcul.quebec>      3 commits

       13b4a69a GenPipes - Tumor Pair : setting Python2 for Strelka2 in the config file
       69201a36 GenPipes - Tumor Pair : set the output_dirs global variable to better handle outputs, thus dependencies
       a9bf30ce GenPipes - Resourcess : adding pycoqc install script + updating other scripts

  Édouard Henrion <henrione@narval2.narval.calcul.quebec>      4 commits

       70c9b89a GenPipes - Tumor Pair : fixed dependency issues by using relative path everywhere in the pipeline, no more absolute path. Also setting 'samples' for all the jobs.
       1156e33a GenPipes - Tumor Pair : updating mugqic_tools to latest version 2.10.2
       2a541dec GenPipes - BFX : updated GATK wrappers to remove print_reads outputs before running
       576724b0 GenPipes - Tumor Pair : fixing Strelka & Manta SV step

  Édouard Henrion <henrione@narval3.narval.calcul.quebec>      1 commits

       bab3a852 GenPipes - Tumor Pair : Fixing conpair with specific call to Python2

  ehenrion <edouard.henrion@mcgill.ca>      16 commits

       927e40e4 Merged in release_4.2.0 (pull request #353)
       a546428e Merged in release_4.1.3 (pull request #352)
       f1e6716e Merged in release_4.1.3 (pull request #351)
       43e87b72 dnaseq.base.ini edited online with Bitbucket
       cdcebcb1 Merged in tumor_pair_hotfix_dev (pull request #339)
       02050579 chipseq.base.ini edited online with Bitbucket
       f52b3e24 Merged in tumorpair_hotfix_strelka (pull request #331)
       4db4e81b Merged in dnaseqsv_fix_eh (pull request #329)
       66eca370 Merged in fix_dnaseqsv_eh (pull request #320)
       e634cddf Merged in fix_delly_call_dnaseqsv_eh (pull request #317)
       2c20036d GenPipes - Scheduler : stop using python from ini, use default mugqic python3 instead to ensure job2json is working in al cases
       2ee5d49e GenPipes - COVSEQ : corrected covseq.graham.ini for dna_sample_qualimap
       0b5b8967 GenPipes - COVSEQ : corrected covseq.cedar.ini for dna_sample_qualimap
       41b755e3 GenPipes - COVSEQ : corrected covseq.beluga.ini for dna_sample_qualimap
       2364a6c8 GenPipes - Readset : corrected "writer" call in checkDuplicateReadets
       a1be2b96 Merged in release_4.1.2 (pull request #314)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       97a0701d Merged in JoseHector-Galvez-Lopez/nanopore_covseqbaseini-edited-online-wit-1645559267177 (pull request #316)
       7c1b1f53 Final adjustment to the nanopore_covseq base ini to optimize job submission on Abacus.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      54 commits

       9c240995 General - general - Adding cpulimit to be used by default on Abacus for all jobs (discussed with PO)
       65971dc0 Merged in covseq_stretenp (pull request #341)
       a03cb250 Merged in HotFix_dev (pull request #342)
       10b13c0a covseq - prepare_report - Adding pandoc from cvmfs
       9f12d3ec General - ini - Updating rnaseq and dnaseq ini
       058e9d04 rnaseq - star_index - Auto-Calculating star genomeSAindexNbases based on genome length
       8e4d4d5a rnaseq - bedtools - Not normalizing if input reads too low
       d1261c0c Merged in HotFix_dev (pull request #340)
       129c7090 methylseq - bismark - Hardcoding the multicore value for now
       4c787bde General - ini - Standardizing ini base section
       5c024003 Merged in Paul-Stretenowich/covseqbaseini-edited-online-with-bitbuck-1649436403200 (pull request #336)
       e29820b1 covseq - ini - Fixing bam2fq number of threads hard coded
       75863a05 Merged in HotFix_dev (pull request #334)
       15069c5e chipseq - macs2_callpeak + homer_annotate - making genomre_size not required in ini
       2456dec9 Merged in HotFix_dev (pull request #330)
       3f46d5c3 chipseq - homer_annotate_peaks - Adding way to customize genome size in ini for cit
       47ab87a4 chipseq - macs2_callpeak - Adding way to customize genome size in ini for cit
       f9878828 General - scheduler - Addressing PO's comment
       835b0f56 chipseq - homer_annotate_peaks - Debug genome_size
       2942015c chipseq - Misc - Adding genome size + log to .o file in batch mode
       a57654c5 chipseq - homer_annotate_peaks + homer_find_motifs_genome - Adding way of using genome fasta and not only uscs naming
       ee79ed56 chipseq - homer_make_ucsc_file_bigWig - Adding ini_section variable + IMPORTANT: adding way to skip exit code 141 due to (z)cat | head sending SIGPIPE despite not being an error
       4c9fbc22 chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       4fbf2740 chipseq - homer_make_tag_directory - Adding way of using genome fasta and not only uscs naming
       7c972fea General - scheduler - Fixing batch mode job2json
       81ceffe2 General - scheduler - Fixing batch mode job2json
       484353ca General - scheduler - Fixing batch mode
       a9d439cc General - scheduler - Fixing batch mode
       81187562 chipseq - base.ini - Fixing cluster_cpu old way
       3d5c6a37 Merged in covseq_nanopore (pull request #325)
       5d8cadb6 nanopore_covseq - General - Code cleaning
       564fc035 nanopore_covseq - General - Code cleaning
       06c695b6 covseq - rename_consensus_header - Updating date from 2021 to 2022 by default
       20df473f nanopore_covseq - General - Adding step_wrapper to slurm and batch mode + Adding pipefail in jobs to catch errors in pipes
       5231d853 nanopore_covseq - General - Changing default resources
       9db1946f nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       cba39ac5 nanopore_covseq - prepare_report - If step_wrapper not defined in ini it'll be empty by default
       f91a244b nanopore_covseq - jsonator - Adding NanoporeCoVSeq for json Nanopore
       e818ccf6 nanopore_covseq - General - Upgrading to python3 ini by default
       27453fbe nanopore_covseq - prepare_report - Re-implementing cpulimit to be set at ini level and not automatically
       73b5b5bd Merged in HotFix_dev (pull request #328)
       2d75e438 dnaseq - General - HotFix
       01c2b49b Merged in chipseq_urgent_fix (pull request #326)
       96a5f32f chipseq - General - Upgrading to python3
       17f337f2 chipseq - differential_binding - Adding pandoc to module load
       da53a9de Merged in chipseq_urgent_fix (pull request #322)
       97dfd745 Merged in covseq_nanopore (pull request #321)
       0215ccf6 chipseq - macs2_callpeak - Fixing typo error
       32d6175b nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       c4e2ad5d nanopore_covseq - prepare_report - Adding cpulimit for Abacus debug
       a8017e37 nanopore_covseq - prepare_report - Adding cpulimit for Abacus
       20f047ed Merged in Fixing-Issue-#144 (pull request #315)
       2b3c2b55 chipseq - macs2_callpeak - Fixing other_options to be used via ini
       38d0c9a4 chipseq - macs2_callpeak - Fixing other_options to be used via ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       868b49de Merged in tp_reports (pull request #348)
       36498184 Merged in tp_reports (pull request #347)
       ea6260fe Merged in new_ini (pull request #319)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      20 commits

       98bc236c tabix ovewright existing index in tumor pair
       166b252f fix ram in depth coverage gatk
       433f7a32 fix format2pcgr config
       42307914 up mugqic tool version in tumor pair
       5a9f72c7 fix compute_effect with cancer
       e42f7782 revert R for gq_seq_utils_exploratory_analysis_rnaseq
       4eb6fd65 python 3.10.4 for filter ensemble
       edc55691 remove sed from manta_sv
       de2d4252 text file in not compress
       7fabc99b gatk_indel_realigner does not support multithread
       34d2e11f Cleanup dev files
       8e0e9fe3 typo
       cee03217 typo ini files
       d40a229d let cpulimit do its stuff on abacus
       4c6514fa let cpulimit do its stuff on abacus
       0e25a16d fix rnaseq denovo crashes
       3a884782 disable reporting to mgcill.hpc with NO_MUGQIC_REPORT variable set
       2924b41b fix nodes set in pbs
       6abbcc4e more robust log_report.py
       d93e0947 new ini file setup with common ini for clusters

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       9ea6e292 fixed a bug in input file list of haplotypecaller
       81cce4e4 Added new parameters to change the FDR and P-value of the differential binding analysis (chip-seq)
       c7ff9976 fixed a bug in input file list of haplotypecaller
       90027831 fixed a bug in input file list of haplotypecaller
       75c58e42 fixed a bug in input file list of haplotypecaller
       4d7ba2cd fixed a bug in input file list of haplotypecaller
       c3b27d8b fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       193ef5a5 Added variant calling using GATK4 to chipseq pipeline

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar1.cedar.computecanada.ca>      6 commits

       c5d8d3e9 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       7605e7ce fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       681a6b86 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       ed358a8d Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       66e6e90c fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       3e29640a Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts

  Pubudu Nawarathna Mudiyanselage <pubudu@cedar5.cedar.computecanada.ca>      5 commits

       40cf8d3f fixed a bug in input file list of haplotypecaller gatk
       5c903d45 fixed a bug in input file list of haplotypecaller gatk
       7417a637 fixed a bug in input file list of haplotypecaller gatk
       2aca8660 fixed a bug in input file list of haplotypecaller gatk
       9c2c6988 fixed a bug in input file list of haplotypecaller gatk

  Pubudu Nawarathna Mudiyanselage <pubudu@narval1.narval.calcul.quebec>      52 commits

       1510c0f7 changed job name for merge_gatk
       0442a69a removed duplicated haplotype calling functions
       08af8f18 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       ec3bd09e Modified interval_padding in dnaseq, exome-seq and chipseq
       9f0b3d37 fix padding issue
       36996431 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding resolved merge conflicts
       a40ed3af Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d111f5f7 fixed a bug in input file list of haplotypecaller gatk
       4c686e83 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       f595d719 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       eeda83bf Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       0b0f2a23 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       d84db852 Added variant calling using GATK4 to chipseq pipeline
       abdc2162 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       127c4b4d Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       7a650a7b Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e74ac076 fixed a bug in input file list of haplotypecaller gatk
       af8bb910 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       3e1af780 Added variant calling using GATK4 to chipseq pipeline
       896157b1 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       f523e62d Added variant calling using GATK4 to chipseq pipeline
       67f4788c corrected diffBind path
       fc76ea5f Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       95592a42 Modified interval_padding in dnaseq, exome-seq and chipseq
       7cee1016 fix padding issue
       1e716736 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       102b90a3 fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       4ca359b1 Added variant calling using GATK4 to chipseq pipeline
       588b3fb7 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       89e4ed42 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       56bcbc6b Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       1879f7eb fixed a bug in input file list of haplotypecaller gatk
       a33712db fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       faec78c3 Added variant calling using GATK4 to chipseq pipeline
       d5189ed7 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       a120b8f6 Modified interval_padding in dnaseq, exome-seq and chipseq
       13b93e9f fix padding issue
       6f6e56c4 Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       f2e52991 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       f5d7e5ec fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       1ba29fb5 Added variant calling using GATK4 to chipseq pipeline
       b2132fa7 fixed a bug in input file list of haplotypecaller resolved merge conflict extended the function to change inter-padding
       6d18c042 Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       d4616a36 changed mugqic_tools version to version 2.10
       9e33161b Merge branch 'chip_snp_pubudu' of bitbucket.org:mugqic/genpipes into chip_snp_pubudu
       7a7a77f7 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       960a2c9d fixed a bug in input file list of haplotypecaller gatk
       7143c0dd fixed a bug in input file list of haplotypecaller extended the function to change inter-padding
       c79776a4 Added variant calling using GATK4 to chipseq pipeline
       b31f874a Added variant calling using GATK4 to chipseq pipeline resolved merge conflicts
       13318779 Added some comments to bfx/gatk4.py to explain why the new changes have been done.
       e051b83d changed default python version to python3

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      4 commits

       6174fcb8 Merged in chipseq_snp_hotfix (pull request #335)
       9213853e Merged in chip_snp_pubudu (pull request #324)
       72ef4844 Merged dev into chip_snp_pubudu
       df9593d9 Merged in chip_seq_issue_150_fix (pull request #327)

  Robert Eveleigh <eveleigh@narval1.narval.calcul.quebec>      12 commits

       b43c90d4 snpEff cancer pair file fix
       7c4ce711 panel snpeff fix
       853d0b2a merge varscan2 fixes
       88c2b0e2 add ram options to gatk_interval_list2bed
       5757d2ca dependency fixes at various steps
       f3d82fcc bwa to indel realignment dependency fix
       85f21e45 changed mugqic tools version
       e5a12355 use purple annotated strelka2 output for somatic ensemble calling
       99b328ca update module verions - issues with bcftools
       41b57cea cit mugqic_tools version fix
       87768ac5 resource fixes: tumor_pair.exome.ini
       8124f350 resource fixes: dnaseq.base.ini, tumor_pair.extras.ini and cit.ini

  Robert Eveleigh <eveleigh@narval2.narval.calcul.quebec>      9 commits

       2c41ddc9 bcftool -i/-e fix for tumor pair fastpass varscan2 merge
       0c23107b abacus resource fixes for merging and filtering vcf steps
       8a4df12f memory fixes to merge vcf steps
       85860f89 vardict and varscan2 germline filter fixes
       00dfdf32 fixes to conpair for multiqc intergration
       ff85c189 tumor_pair.dev.ini fix
       a5c24ac5 change file name : tumor.dev.ini to tumor_pair.dev.ini
       75bc9dcd reverted back to python2 for strelka2
       35e5c27f adding cpsr/pcgr reporting with addition of manta, cnvkit to ensemble protocol.  Step clean up as well

  Robert Eveleigh <eveleigh@narval3.narval.calcul.quebec>      7 commits

       653c5c63 more germline filter fixes
       dc401001 resource fix to manta and strelka2, and germline filter improvement
       0a3aef5d add snpEff QCs and fixed conpair concordance output
       f0749eb8 resource fixes
       9632453d add python2 for varscan2 mugqic tool
       6733ac4f varscan2 one job fix - snp/indel merged vcf now in right place
       12726b38 cpsr/pcgr dependency fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       8948532d python3 fixes for support scripts
       f2ea81d3 conpair and ensemble filtering fixes

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      2 commits

       54b97eed resource fixes
       889346c3 python3 fixes when running on abacus

4.1.2        Thu Feb 17 17:40:12 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       7997ec59 Merge remote-tracking branch 'origin/dev' into release_4.1.2
       c83a283f Version bump to 4.1.2
       dd9f1853 Updating pipeline READMEs before release
       4bcda971 Version bump to 4.1.1

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       e5a64db9 Merged in release_4.1.2 (pull request #313)
       99b40bf6 Merged in release_4.1.2 (pull request #312)
       dc6f51e6 Merged in release_4.1.2 (pull request #311)
       7dd399ff GenPipes - Scheduler : Correcting cluster_cpu parsing for PBS
       7a2aa6db Merged in release_4.1.1 (pull request #306)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       c07382b9 Merged in ont_covseq_report_ini (pull request #308)
       46e25d14 Adjustment to the `prepare_report` step resources in the base ini file, to try to prevent this job getting stuck in the queue.

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       4d43f87f Merged in patch_4.1.2 (pull request #309)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       98951465 remove cluster_submit_cmd from slurm and pbs ini
       df8c9ba5 fix hic ini
       f45fdba4 remode empty quotes
       a955e743 cleanup queue
       44175526 PBS gets nodes and ppn together

4.1.1        Thu Feb 10 15:12:29 2022 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      10 commits

       dc3ff3c2 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       d2a5407b minor README uedit, again...
       3cd7d347 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       21ee999d Another minor correction in the README...
       549f6dc2 Merge remote-tracking branch 'origin/dev' into release_4.1.1
       bbc6156d Updating READMEs before release
       150f4a01 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       77e0698d Adding missing space in container help message
       be595604 GenPipes - prep for minor release
       11b8ec2f Version bump to Release 4.1.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       aeccdc34 Merged in release_4.1.1 (pull request #305)
       c756a500 GenPipes - Tumor Pair : Replacing "xrange" calls by "range" calls to resolve Python3-compatibility issues
       95f8ce82 README.md edited online with Bitbucket
       691b2fdc Merged in release_4.1 (pull request #304)
       0204b7cb README.md edited online with Bitbucket
       ea152c3e corrected typo in README.md for Nanopore_covseq
       ae1c04a9 Merged in release_4.1 (pull request #303)

4.1.0        Mon Feb 7 21:39:25 2022 +0000        77 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      11 commits

       d00f24a3 GenPipes - Release : creation of the release branch
       03e69192 Merge remote-tracking branch 'origin/dev' into release_4.1
       185575e8 Correcting typo in READMEs
       ebc18e3f GenPipes - Release : updating READMEs and VERSION
       30e05e5a GenPipes - Nanopore : correcting minimap2 step
       7fb2d7c0 GenPipes - BFX : cleaning bash_cmp.py a bit
       f7538c5a GenPipes - Nanopore : fixing minimap2 and dependencies + allowing duplicate samples (if no duplicate readset) in the readset file
       07287137 GenPipes - Utils : renaming of  to
       cf7b0c65 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1da114b5 adding .gitignore to .gitignore...
       880fbf48 Version bump to 4.0.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       e503967d Merged in release_4.1 (pull request #302)
       8fe012c7 GenPipes - AmpiconSeq : updating resources for dada2 in cedar.ini
       8a0ea558 GenPipes - BFX : fixing "ln" in bash_cmd.py
       44a9dec2 Merged in fix_nanopore_eh (pull request #298)
       44ddf6bc Merged in ehenrion/genpipes-chipseq-bashini-edited-with-u-1642626170885 (pull request #292)
       0f84b782 GenPipes - ChIP-Seq : bash.ini edited with updated versions of software, fixing Issue #127
       6c8219d2 Merged in release_4.0.0 (pull request #286)

  jgalvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       e0c9f5f7 Full squash of covseq_ont branch

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      5 commits

       2841e22b Final correction to the reference genome symlink
       323b5d20 Corrected issue with the reference genome link command
       35ad8b75 Added compatibility to ARTIC primer versions 4 and 4.1
       3f3fd3f1 Update version of ncov_tools to 1.8
       8bb76cd8 Added changes to the ini files to allow for ARTIC V4 schemes (and any potential future schemes).

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       fa5b71bc Merged in covseq_v4 (pull request #287)
       ce250554 Added force option to reference symlink creation to avoid crash when re-running report.

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      10 commits

       8915d42a Merged in covseq_nanopore (pull request #300)
       a242cd1b nanopore - Cit update - Updating cit ini
       6ab04881 nanopore - pycoqc - Fixing step to match with new bfx
       1f4b6307 covseq_nanopore - Cit update - Updating cit ini
       803b08af covseq_nanopore - Cit update - Updating cit ini + removing module load in report
       34e981c1 Merged in covseq_nanopore (pull request #299)
       eaa85d8e covseq_nanopore - Cit - Code cleaning and cit.ini addition
       97869710 covseq_nanopore - Cit update - Updating cit ini
       97267840 covseq_nanopore - General - Renaming covseq and nanopore_covseq class for consistency
       fff0c8b1 covseq_nanopore - General - Fixing argparse type

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       6eddc14a Using cit magic for walltime
       c84a3b4c Merged in cpu_scheduler_agnostic (pull request #294)
       f73f008c Full squash of covseq_ont branch
       06702781 Merged in config_formater (pull request #289)
       f0047ac3 Merged in watchdog (pull request #290)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      30 commits

       2c5f270b Add mkdir in FFPE steps
       2860a3b7 cpu_str to node_str in scheduler class
       42052473 update cit.ini
       d67a799e update nanopore cit.ini
       2c679b9c make cluster_cpu scheduler agnostic
       5694384d x on tumor_pair.py
       13bd411e get symlink back in tumor pair
       b919a3e3 get bash_cmd back to old self
       9d722c25 fix regressions
       7d068c2f ad tmp to gitignore
       5c198199 update perl in tumor pair
       19aaf711 cluster_memmory -> cluster_mem
       88a70bc7 tweak tumor pair extra ini
       3fe6b0bb missing int casting
       ce9a2e81 pmem for pbs/torque
       f3c406dd cleanup sanitycheck log
       4357a19c SanitycheckError args and kargs
       fe6b300b add .message to SanitycheckError
       ddc14010 add cluster_mem to ini option
       e0aec1b8 force int on walltime read
       87f84cd2 force int on walltime read
       220b6e0b slurm to time_to_datetime
       853ee6ff fix hicseq ini
       07f2b77c walltime method from slurm
       ad88e69d time delta for pbs walltime
       1227182e update utile time delta
       2f679b7b update utile time delta
       2fe4ab3e update utile for torque
       88abd413 add walltime format for pbs
       333c32a5 rename monitor.sh to watchdog

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      5 commits

       80e3e4d3 chipseq diffbind pca plot update
       dc0680a2 methylseq methylkit min cpg issue fix
       f87dde23 testing updates to DiffBind.R in mugqic_tools
       acaaa651 modified differential expression variable names to something meaningful
       da9bf667 fixed issue #129 added differential binding to the atacseq protocol

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       5919bcf9 Merged in chipseq_atac_diffbind (pull request #288)

4.0.0        Thu Dec 9 19:13:55 2021 +0000        232 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      14 commits

       a34c5a19 Merge remote-tracking branch 'origin/dev' into release_4.0.0
       d2f25eb3 updated .gitignore
       c02c96d7 Merge remote-tracking branch 'origin/dev' into release_4.0.0
       9485b938 GenPipes - README : update epiqc REAMDE
       18dc4713 Merge remote-tracking branch 'origin/dev' into release_4.0
       8ac68476 GenPipes - README : Updated READMEs before new release
       b1e2df04 GenPipes - config : correcting parsing of walltime from config files
       b83fe9fc Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0865f206 GenPipes - DNASeq : correction of sym_link_final_bam dependencies - fixes Issue #119 and Issue #125
       85f8054f GenPipes - Tumor Pair : fixed some enumerate loop...
       5199ee7f GenPipes - Tumor Pair : fixing sym_links steps to avoid duplicate job names
       0e63d7fa GenPipes - DNASeq : fixing .bai dependencies
       3453cf00 Version bump to 3.6.3-beta
       d2361a56 Version bump to 3.6.2

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       963a49a7 GenPipes - python3 : fixes following py2to3 recommandations
       a2274713 GenPipes - Python3 : fixes after last rebase with dev
       4b5f9bda Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       5ca19c5f GenPipes - start switching to python3
       bfa2a48f GenPipes - start switching to python3
       e3f7e928 GenPipes - Python3 : fixing file opening in bfx/readset.py
       77896893 GenPipes - Python3 : fixing conflicts with bfx/readset.py in dev
       7c41602d GenPipes - start switching to python3
       978a291c GenPipes - start switching to python3
       233c1086 GenPipes - start switching to python3

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      4 commits

       2601b7f7 GenPipes - Python3 : removing the shebang from all the bfx scripts
       a99a2e42 GenPipes - Python3 : removing the shebang from all the bfx scripts
       c43eb204 Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh
       025c3c2f Merge branch 'genpipes_python3_eh' of bitbucket.org:mugqic/genpipes into genpipes_python3_eh

  Édouard Henrion <henrione@beluga5.int.ets1.calculquebec.ca>      2 commits

       4463bace GenPipes - Config : removed useless debugging messages
       d47adcfe GenPipes - Resources : updated install scripts for fgbio and mugqic_tools

  ehenrion <edouard.henrion@mcgill.ca>      8 commits

       baf52f94 Merged in release_4.0.0 (pull request #285)
       540e030d Merged in release_4.0 (pull request #284)
       fe49a832 GenPipes - BFX : corrected typo in gatk4.py
       bdea57fe Merged in tumor_pair_sym_link_fix (pull request #282)
       f5947541 Merged in hotfix_dnaseq (pull request #280)
       34e54c90 GenPipes - RNASeq : minor typo fixes in README.md
       eca4475d Merged in genpipes_python3_eh (pull request #270)
       34608bc7 Merged in release_3.6.2 (pull request #269)

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       b28b4a8c remove recalibrated bam name error from some input jobs

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       2b234225 Merged in dnaseq_fixRecal (pull request #276)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      8 commits

       73443c7a Making pandoc report working with pandoc version 2.16.1
       60093e6b Making pandoc report working with pandoc version 2.16.1
       c6ab31bd Merged in Paul-Stretenowich/log_reportpy-edited-online-with-bitbucke-1637772592796 (pull request #277)
       e8692599 log_report.py switching back to py3 and adding Narval as remote
       0ef01e96 EpiQC - First commit after Rami Coles internship
       4e55e0fb General - Fixing genome installation grep issue
       06260f0d EpiQC - First commit after Rami Coles internship
       d458b422 Merge branch 'dev' into epiqc

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      9 commits

       4f89f915 copy contents of readset into inputinfo file(have errors)
       d51b1cea copy contents of readset into inputinfo file(have errors)
       97bc3da3 creating inputinfo file after copying original file
       45c917cf epiqc - completed chromimpute development - right before remove design file and modify readset file
       95008563 epiqc - completed chromimpute development - right before remove design file and modify readset file
       2bb45e0a [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       3a86f948 copy contents of readset into inputinfo file(have errors)
       acfcbc04 copy contents of readset into inputinfo file(have errors)
       7b5b552e creating inputinfo file after copying original file

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      7 commits

       43f45759 update GIAContainer version
       9fec3337 changing error log to debug in cit config
       9329e2e8 add cit options for epiqc
       10199ccd send telemetry downloded file to /dev/null
       e942d077 fix warning
       99ec4c92 make log_report.py more robust
       c30189c9 update R bioconductor in dnaseq_high_coverage

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      47 commits

       4bee1cab bug fix #121. missing module load mugqic_tools 2.8.2
       5e8f3dab fixed file path issue in seq2fun corrected the mugqic_tools version in epiQC
       a317b6bd corrected a typo
       7f780053 added the reason to use os.remove in a comment
       c4721bc9 addressed Hector's comments - changed seq2fun db path
       ffd0abe5 addressed Hector's comments
       6e2ce1f9 fixed Issue #118: hicseq hic error in interaction_matrices_Chr (mugqic/genpipes) by adding floor division to places which return floating points where intgers needed
       c333a3e8 modified incorrect seq2fun cvmfs db path
       5e928cfc added seq2fun cvmfs db paths
       0c9418a9 modified seq2fun db paths
       1569435d updated readme file
       22d90bf7 modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       dc4c6088 seq2fun pathway_analysis completed
       de5f87b3 started pathway analysis-seq2fun
       e1627475 improved seq2fun processing
       5d6df189 added seq2fun to the function as the final step. working well but needs improvements
       bdd622d0 fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       54bf3aa8 added seq2fun protocol completed merge_fastq partially completed seq2fun function
       0c2d5693 addressed Paul's comments
       ef8ae360 fixed a bug
       2e3c7428 changed mugqic_tools version to 2.8.3
       06e458ef fixed bugs when transferring to python3
       fb9093e1 fixed bug in chromimpute convert corrected uuid issue in job2jason.py
       0bd67730 fixed an issue for python 3
       0d1622c5 updated readme file
       9e625460 testing whether env varibales are retrievable dynamically
       ff73039b fixed input files not found error when running the pipeline outside of the project directory
       78af37b9 added cit.ini
       a4dc68c9 fixed issues in final report added comments
       ca3dc8aa fixed bugs in -o option added IHEC data path
       0bcde687 modified documentation. added extend jobs for differential expression modified ini files for cedar and graham
       affd3b06 seq2fun pathway_analysis completed
       129b56b5 started pathway analysis-seq2fun
       7ef4a6fd improved seq2fun processing
       92ac7120 added seq2fun to the function as the final step. working well but needs improvements
       97bdedbf fixed fastq concatenate issue with adding zcat if the file is gz changed wall time for merge_fastq and seq2fun
       2d0359b4 added seq2fun protocol completed merge_fastq partially completed seq2fun function
       93811b3f testing whether env varibales are retrievable dynamically
       72d5f9bc fixed input files not found error when running the pipeline outside of the project directory
       b97891f2 added cit.ini
       b29e7177 fixed issues in final report added comments
       cea9f45f fixed bugs in -o option added IHEC data path
       15c56959 completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       b86e6e56 completed fixing errors-working pipeline
       d07286d8 corrected up to epiqc final report
       9304cba6 corrected up to epigeek
       562607a9 corrected up to global dist after chipseq design changes

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      59 commits

       e8371668 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       5683bc99 Merge branch 'seq2fun_denovo_rnaseq' of bitbucket.org:mugqic/genpipes into seq2fun_denovo_rnaseq
       faab0eed Addressed Paul's comments
       b3e81850 addressed Paul's comments
       4130054c Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       0225132c space changes in epiqc.py
       ca91e09b completed all the steps, working pipeline. documentation is 90% completed. might need to do some fix on chromimpute
       dbf82549 completed fixing errors-working pipeline
       fc17d37a General - Fixing genome installation grep issue
       26837e3c EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       e3ca2e65 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       9bda3533 EpiQC - Fixed bugs
       e4faee46 EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       1dcd8dd6 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       851b26da added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       29139a98 General - Fixing genome installation grep issue
       b95fc32e resolve merge conflicts
       701f1099 resolved merge conflicts
       d03737e5 resolve merge conflicts
       8852d61c resolved merge conflicts
       54634868 resolved merge conflicts
       82ff1002 resolved merge conflicts
       d833412d resolve merge coonflicts
       e047a3e3 resolved merge conflicts
       dba6b339 space changes in epiqc.py
       2fe70467 Merge branch 'epiqc' of bitbucket.org:mugqic/genpipes into epiqc
       429caf1e General - Fixing genome installation grep issue
       06113d17 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       34d3cc94 EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       c00875a3 EpiQC - Metrics thresholds can be modified in base.ini
       35a63962 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       9e3e1f4e \Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       698c46ff EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       2086d5af EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c2af489d EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       1089984b EpiQC - Fixed bugs
       bf60df96 EpiQC - Added signal to noise and epigeec steps
       bced2bd5 EpiQC - Parallelized chromimpute
       c1b462ab EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.
       4d12401f EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       0e763c87 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7f8337e7 [epiqc] - developed up to generatetrain data - need to run through alll histone marks in the inputinfo file
       2b769bff resolve merge conflicts
       05b7cd3a resolved merge conflicts
       2e219fa5 resolve merge conflicts
       aaed941e resolved merge conflicts
       cc0bccf1 resolved merge conflicts
       231f6264 resolved merge conflicts
       a7e5a61e resolve merge coonflicts
       42a741bb resolved merge conflicts
       1e5e8af1 Merge branch 'epiqc' of https://bitbucket.org/mugqic/genpipes into epiqc
       6eef5c8d resolve merge conflicts
       d679a978 resolved merge conflicts
       eb880f5d resolve merge conflicts
       06fe31e2 resolved merge conflicts
       22965391 resolved merge conflicts
       0e2c9cdb resolved merge conflicts
       301d2eec resolve merge coonflicts
       f285a7c1 resolved merge conflicts

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      8 commits

       b3194fa8 Merged in seq2fun_fix_bug121 (pull request #281)
       ffc9c01f Merged in epiqc_mugqic_tools (pull request #279)
       a9ac4d61 Merged in epiqc_comments (pull request #278)
       4bd2b26d Merged in epiqc (pull request #271)
       c38208b8 Merged in seq2fun_denovo_rnaseq (pull request #272)
       e12b1118 Merged dev into seq2fun_denovo_rnaseq
       0df8920d Merged in hic_python3_issue_fix (pull request #275)
       1d5f2b07 Merged epiqc into epiqc_ss

  rami.coles@mail.mcgill.ca <rcoles@abacus1.ferrier.genome.mcgill.ca>      2 commits

       70f8ec03 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec
       5feada10 EpiQC - report step now creates a heatmap from the correlation matrix obtained with epigeec

  rami.coles@mail.mcgill.ca <rcoles@abacus2.ferrier.genome.mcgill.ca>      11 commits

       f20a0000 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       74e2a32a Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       e911f955 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       7a5da006 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       bc66beb5 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       d0a7b9d8 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6695197a Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       a040d0b2 added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py
       6a5e9958 EpiQC - Added bigwig column to readset.py, bigwiginfo checks if there is bigwig column in readset, if not searches for bigwig file chipseq run
       53c747f5 Removed .pyc and .DS_Store files, updated epiqc.base.ini and epiqc.py
       f56d086f added epiqc pipeline, bfx/bigwiginfo.py. bfx/chromimpute.py, bfx/signal_noise.py, bfx/epigeec.py

  rami.coles@mail.mcgill.ca <rcoles@abacus3.ferrier.genome.mcgill.ca>      21 commits

       69012379 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       88c422d4 EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       f85d37a7 EpiQC - Metrics thresholds can be modified in base.ini
       b1465486 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       51b6bdb3 Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       9e118e51 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       a68aad38 EpiQC - Added signal to noise and epigeec steps
       fe957a6b EpiQC - Parallelized chromimpute
       ab1425d2 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       1cbf6d00 Added comments to the utility tools. ihec_json_parser.py creates a readset for epiqc and a list of marks to add in the base.ini file in the marks field of chromimpute section.
       86c04611 EpiQC - Stable version on abacus. Added utility tools in epiqc directory
       ef3acdea EpiQC - Changed how wigSignalNoise.py and epiqc_report.py are called
       64250463 EpiQC - Metrics thresholds can be modified in base.ini
       e9e9f2f8 EpiQC - Modified epiqc.py and epiqc inis, removed epiqc_report and wigSignalNoise from bfx (adding to mugqic_tools)
       3ecaebb6 Added documentation and comments on epiqc.py, bigwiginfo.py, chromimpute.py, epigeec.py and epiqc_report.py
       2a773e75 EpiQC - Created bigwig_to_bedgraph step and seperated chromimpute into 2 steps : chromimpute_train_step and chromimpute_compute_metrics.
       c800a679 EpiQC - Added epiqc_report in bfx and in main file. Modified chromimpute and bigwiginfo
       d3bd3e12 EpiQC - Fixed bugs
       402b34ed EpiQC - Added signal to noise and epigeec steps
       a6643145 EpiQC - Parallelized chromimpute
       4bf48b97 EpiQC - Added jobs for chromimpute in bfx/chromimpute and epiqc.py.

  Robert Syme <rob.syme@gmail.com>      1 commits

       f8dfff54 Merged in servername-fix-dev (pull request #274)

  Rob Syme <rob.syme@gmail.com>      1 commits

       a76c958a Fix example server name in dev

  Shaloo Shalini <shaloo.shalini@gmail.com>      4 commits

       db1afc37 Merged in ss_mermaid_96 (pull request #247)
       dab640b6 Merged dev into ss_mermaid_96
       6193cdac Merged in dev_covseqdoc_ss (pull request #258)
       6dcc0cb7 Merged in epiqc_ss (pull request #253)

  shaloo <shalz@hotmail.com>      14 commits

       9e547621 Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       4ec66ffc Refs #102 removed ChIP-seq pptx link
       78376887 Fixes #102 issues in README.md for epiqc pipeline
       e08c73d8 Refs #103 Paul's feedback in workflow incorporated
       52ca4d12 Fixes #103 covseq workflow added
       f356ffe4 Merge remote-tracking branch 'refs/remotes/origin/dev_covseqdoc_ss' into dev_covseqdoc_ss
       878fcf5e Refs #103 covseq doc freebayes update
       8f5625c8 Typo fix
       81ee462d Fixes typo Refs #103
       a62e20bf Refs #103 covseq doc freebayes update
       fb32de59 Fixes #101 epiQC pipeline workflow added - both manual as well as mermaid generated
       86e5a38f Refs #102 removed ChIP-seq pptx link
       29ae522d Fixes #102 issues in README.md for epiqc pipeline
       43501a4b Fixes #96 all GenPipes workflows are now coded as mermaid flowcharts and png generated from them

3.6.2        Thu Nov 11 14:14:18 2021 +0000        28 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      4 commits

       3787dea8 Cleaning before release
       c1aef1d9 GenPipes - Readset : improved readset parser so that it creates a temptative readset file with unique IDs when the readset file provided has duplicate readert IDs
       cead55d8 GenPipes - MethylSeq : correted bed2interval_list calls + some minor code reformating
       499a0f63 Version bump to 3.6.1

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      3 commits

       55074b1d GenPipes - README : updating READMEs and install scripts before release
       448d7c85 GenPipes - updating version of mugqic_tools in all the pipelines config files
       1c2d3223 GenPipes - Resources : adding some new genome and software install scripts

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       eb019794 Merged in release_3.6.2 (pull request #267)
       2fb1dc1f GenPipes - Readset : correcting the parsing of readset files to stop allowing duplicates headsets in file - fixes Issue #113
       72170d96 GenPipes - RnaDeq denovo Assembly.py : updated merge_trimmomatic_stats outputs to fix insilico_read_normalization_all_report dependencies
       50519a28 GenPipes - RnaDeq denovo Sssembly.py : fixed insilico_read_normalization_all_report dependencies
       2e429b1a Merged in release_3.6.1 (pull request #262)

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       2eca4503 Merged in Pierre-Olivier-Quirion/ampliconseqbaseini-edited-online-with-bi-1635366894762 (pull request #266)
       3d108f93 ampliconseq.base.ini edited online with Bitbucket

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       cc3f279d update all pipelines to pandoc 2.16.1
       cd9c041f new pandoc

  Robert Syme <rob.syme@gmail.com>      3 commits

       4ab8e179 Merged in wget_warn_on_fail (pull request #265)
       d20df3ea Merged dev into wget_warn_on_fail
       b809a05c Merged in rnaseq_protocol_switch_fix (pull request #263)

  Rob Syme <rob.syme@gmail.com>      9 commits

       20154ee6 Fix server name
       31fd5fe5 Test with incorrect server name
       17aea667 switch quote style
       8c4410e9 Intentially introduce an error in the reporting server.
       beb177ee Warn when wget command fails.
       ca54581c Another tiny whitespace fix
       cc33046a Tiny whitespace fix
       b3d620c8 Switch at correct location
       2cbdcee1 Fix protocol mixup for rnaseq pipeline

3.6.1        Wed Sep 29 20:53:55 2021 +0000        17 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      6 commits

       27ccb3b1 correct log_report.py before release
       39e1b94d Merge branch 'master' of bitbucket.org:mugqic/genpipes into release_3.6.1
       0fae24c1 Updating GenPipes README before release
       df8f9f45 Merge branch 'dev' into release_3.6.1
       96ce81fb Re-creating READMEs before release
       530576e1 Version bump to 3.6.0

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       1cda65c2 Merged in release_3.6.1 (pull request #261)
       b472e728 Merged issue_105_fix_eh into dev
       7cd2f094 GenPipes - Tumo Pair: correcting job name in sym_link_fastq_pair
       4df10424 GenPipes - Tumor Pair.py : fixing typo in sym_link_fastq_pair
       d4613be2 GenPipes - Tumor Pair : fixing issue #106, no more job name overwriting at gym_link_fastq_pair step
       3079e608 GenPipes - DNASeq : fixed bed2interval_list call in gatk_haplotype_caller step
       03c03bf8 Merged in release_3.6 (pull request #259)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      2 commits

       79470425 Merged in Paul-Stretenowich/sambambapy-edited-online-with-bitbucket-1630531961392 (pull request #260)
       69e1a121 Fixing sambamba sort output files: there were a bai set as output that was blocking the "resume" for covseq pipeline.

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       a0262071 more robust log_report
       b941202d remove echo debbug in get wrapper

3.6.0        Mon Aug 30 17:55:36 2021 +0000        370 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       cdaafcd7 Merge branch 'dev' into release_3.6
       bea48bd4 GenPipes - updating CHANGELOG and VERSION files (after release)
       2ea1988f Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       4554eaff GenPipes - README : updating version to 3.5.0 in all pipeline READMEs
       1f71288e GenPipes - Release 3.5.0 : updating CHANGELOG up to tag 3.5.0

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      16 commits

       9eda417a Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       79bc0fdc GenPipes - Tumor Pair : no more attemps to write where the reference files are (in case of readonly FS, e.g. CVMFS...)
       79aa7dc2 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       86c40fdb GenPipes - Tumor Pair : fixed Strelka2 jobs with mkdir
       ed205e15 GenPipes - Tumor Pair : fixed coverage_bed in strelka steps
       b7873b93 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0eafe86b GenPipes - BFX : fixed GATK4 cluster_crosscheck_metrics command
       a9eb687f Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       ae4cc084 GenPipes - Tumor Pair ; deleted dev.ini from dev branch
       b6fd7399 GenPipes - Tumor Pair : fixed symlink_fastq_pair command
       41ba9845 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       84f350f3 GenPipes - Install script : version change
       7e07739d GenPipes - RNASeq light : fixing typo
       e1c1b504 RNASeq light : importing bash_cmd fix
       3e29063b GenPipes - Release 3.6 : updating the READMEs in dev
       16be1706 GenPipes - Pipelines : Fixed dnaseq_high_coverage inheritance with dnaseq and rnaseq_light & rnaseq_denovo_assembly inheritance with rnaseq

  ehenrion <edouard.henrion@mcgill.ca>      13 commits

       677e9cf8 Merged in release_3.6 (pull request #250)
       3591ec92 GenPipes - ChipSeq : updated mugqic_tools with latest version
       da2b4102 GenPipes - Ampliconseq: updated mugqic_tools version in base.ini
       40cb3b42 GenPipes - AmpliconSeq : updated mugqic_tools version in base.ini
       92e49194 Rnaseq_light.py: kallisto_count_matrix dependency
       18d37ad1 fixing picard2 add_read_groups call
       952ac0ae GenPipes - MethylSeq.py : Fixed pipeline inheritance with DNASeq
       de571d6b GenPipes - CoVSeq : setting the gatk4_java_options in base.ini
       11755c1f GenPipes - BFX : fixing add_read_groups in picard2.py
       fa40033e GenPipes - MethlySeq : fixing sambamba_merge_sam_files dependency
       2f631833 Merged in pipeline_inherit_fix_eh (pull request #249)
       aa3ca795 Merged dev into pipeline_inherit_fix_eh
       273fe548 Merged in release_3.5 (pull request #243)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      24 commits

       a1f3040b Merged in covseq_stretenp (pull request #248)
       8374fdb9 Threshold rename consensus change for LSPQ
       d3c3b5fe Code cleaning
       7bc51893 Merged dev into covseq_stretenp
       376368f6 Upgrading covseq_tools version
       e473abe5 Changing regex for finding negative controls for LSPQ
       d0291feb covseq - prepare_table - debug
       0ed94bb6 Merged dev into covseq_stretenp
       b768d921 Adding bfx for covseq_tools and ncov-tools
       0c8407e6 Upgrading cvoseq_tools version
       6b3f8dec Debug
       7dcbf5d2 Degrouping rename_consensus step for ivar and freebayes
       d476c118 Debug
       00fafc33 Un-grouping freebayes report and ivar report
       1fe319eb Un-grouping freebayes report and ivar report
       6b32016a Debug
       3023357a Changing resources for ncovtools
       58d04972 Debug
       29d800f5 Debug
       8a773f9c Debug
       5ede73b0 Debug
       ce60e365 Debug
       84cda2be Adding freebayes metrics and reporting
       d0592cc6 Renaming prefix parameter in ini_section

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       eeee6b14 Merged in tp_debug (pull request #257)
       fae13493 Merged in monitor_limit (pull request #256)
       6283fd90 Merged in kallisto_rnaseq_light (pull request #254)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      13 commits

       8d09333a Bump GiaC to v2.0.3
       55dd9415 patch stolen from Ed to fix recalibration step
       a89de003 RETRY ret_code corrected
       64b9f79a RETRY limit on monitor
       c638ab3d RETRY limit on monitor
       bf872974 RETRY limit on monitor
       2c1c8952 RETRY limit on monitor
       207a0694 fix ini for graham covid
       203ea908 update tp ini for cedar
       64db31b2 reducing constraint on fit for integration testing
       a17345dc minor fixes before release
       248078c2 fix kallisto dependency
       5583f3fe giac updtated to v2.0.2

  P-O Quirion <pioliqui@gmail.com>      2 commits

       fb635da8 more robust log repport
       55b928dc make log report more robust went jobs are failing

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       e56c39d8 corrected a typo
       9110f867 modified readme files for chipseq differential binding step

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       fbdbe38a Merged in chipseq_readme_fix (pull request #245)

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      40 commits

       ac5a5c59 Homo_sapiens.GRCh38.ini fixes
       675d5ddc samtools mpileup spacing fix
       f1827795 cit fixes after dev.ini removal
       15e98828 multiqc ini fix
       1cdfe282 samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       5977acdc add portal dir to dnaseq
       3f93719a hs37d5 and cit strelka2 fixes
       5a20aa9f removal of cit_gatk4.ini, contents moved to cit.ini
       59473ceb tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       d98dd0ce merge vardict fix
       a64d34c2 job.sample fixes for jsonator
       63b37d85 vardict exome fixes
       13b7d0f7 germline variant annotator fixes
       d366e782 strelka2 germline and sequenza fixes
       19b187c9 cit adjustments to tumor_pair
       b9993756 beluga.ini fixes
       46f7416a strelka2 germline fixes
       3d5693d2 removal of samtools, substitute samtools germline for strelka2
       ec213cec conflict fixes to dnaseq and tumor pair after dev merge
       8cad167d fix mpileup dependency chain
       e6d6aa5b corrections to mpileup bcftools merge and tumor_pair dependencies
       1dc88b91 samtools bug fixes to mugqic protocol, variant recalibrator argument fix for gatk3
       9e3ac20d add portal dir to dnaseq
       4a59db33 portal_dir fix
       4dfa21b9 added portal_dir for cit testing
       528a77be hs37d5 and cit strelka2 fixes
       9c66ec28 removal of cit_gatk4.ini, contents moved to cit.ini
       03a602b2 tumor pair converted to gatk4 only, mutect2 LearnReadOrientationModel added to help with FFPE samples
       30d7b075 merge vardict fix
       83b75bf7 job.sample fixes for jsonator
       ec145472 vardict exome fixes
       e5a61457 germline variant annotator fixes
       fdd20bfb strelka2 germline and sequenza fixes
       a4f8581e cit adjustments to tumor_pair
       7240e3e9 beluga.ini fixes
       849625a4 strelka2 germline fixes
       80f55487 removal of samtools, substitute samtools germline for strelka2
       e92837b1 conflict fixes to dnaseq and tumor pair after dev merge
       9aedd698 fix mpileup dependency chain
       e55e1ffd corrections to mpileup bcftools merge and tumor_pair dependencies

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      45 commits

       6ec0367c collect metrics fix - dnaseq
       fd8ee638 PR fixes to sambamba and multqc
       afc4b8f4 PR conflict fixes
       04c160e4 updated dev.ini
       53cde245 cit fixes
       359c1ae0 purple sanity check fixes
       6d124355 cit fixes to purple and vardict exome
       5da77605 fix for multiple tumors with the same control
       d7724ad5 deliverable fixes - sym link final bams
       63c8bc87 non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       39a26519 bcftools update and filter argument fixes
       67410487 strelka2 germline and ensemble fixes
       b1f73968 sequenza and cit fixes
       c9cef169 fixes to gatk4 for tumor pair
       bc784a60 gatk4 vsqr cit fix and baf plot for b38
       e148bd4f argument fixes for picard imported functions in gatk4 and vqsr fixes
       036809d4 fixes to multiqc
       7b0712be exome specific fixes
       059c613a updates and fixes for cit
       dbe0ff4f cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       3c52e666 remove briaree ini and update dnaseq base
       9343e040 updating beluga ini
       d79e994c updates to beluga.ini and base.ini for dnaseq
       97793f0a ini updates
       f716b0e9 updated dev.ini
       0f7cccc3 cit fixes
       45c9d0bc purple sanity check fixes
       882096d5 cit fixes to purple and vardict exome
       97af93a5 fix for multiple tumors with the same control
       49f4f38c deliverable fixes - sym link final bams
       17f001eb non-diploid GT fixes for mutect2 gatk4 and seqeunza fixes
       6729217a bcftools update and filter argument fixes
       afba10fe strelka2 germline and ensemble fixes
       63264495 sequenza and cit fixes
       e8a4ac9d fixes to gatk4 for tumor pair
       3348617d gatk4 vsqr cit fix and baf plot for b38
       b8ec9bd9 argument fixes for picard imported functions in gatk4 and vqsr fixes
       9b714d09 fixes to multiqc
       ca97ea73 exome specific fixes
       2e20be2f updates and fixes for cit
       67933483 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       8b46abeb remove briaree ini and update dnaseq base
       4795a96f updating beluga ini
       04d851fe updates to beluga.ini and base.ini for dnaseq
       ab85a846 ini updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      37 commits

       6cec0eff module fix
       7b697cc3 removal of dev files for CVMFS version
       46ef162c genome ini fixes
       466eb237 dev and hs37d5 fixes
       367666b7 tumor pair - sambamba markdup fixes
       d5b857fb module fix
       e7a28494 sambamba markdup fix
       f7e4fdbc strelka cit fixes
       c320d5a4 added strelka2 input dependency for purple purity
       d3c32925 cit fixes to purple and vardict exome
       89c80a9f muliple normal fix for fastqc
       00531ff3 dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       d1840a09 fixes and strelka conversion addition to purple
       09d91352 added purity estimate program PURPLE
       aaaff33b strelka bed manipulation fix
       f257d83b vardict exome fix
       f74092d9 conpair and collectHS metric fixes
       2822b5f7 tumor_pair qualimap part 2
       51e88ab5 sym link dnaseq.base into tumor pair
       03e567cd updates to b38 variant recal files
       4757879e fixes to tumor_pair on beluga
       52e48cef cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       febb8234 strelka cit fixes
       42fd97c6 added strelka2 input dependency for purple purity
       54fd0cd3 cit fixes to purple and vardict exome
       858af146 muliple normal fix for fastqc
       bd1f0346 dependency mkdir fix, purple fixes and fix to muliple pairs with same control
       1b2c03fe fixes and strelka conversion addition to purple
       0957242a added purity estimate program PURPLE
       37529e6c strelka bed manipulation fix
       cc317819 vardict exome fix
       2091afb7 conpair and collectHS metric fixes
       8689d63c tumor_pair qualimap part 2
       a9888f74 sym link dnaseq.base into tumor pair
       e340eb04 updates to b38 variant recal files
       7a7b46da fixes to tumor_pair on beluga
       a449ed5d cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      28 commits

       58152ae7 varscan version downgrade
       7fd72e18 cit fixes
       3a63b4ad fixes to bcftools mpileup for dnaseq mpileup protocol
       b92ebd8a hs37d5/GRCh37 fixes for purple
       51c9658c fastpass varscan merge fix
       2cc7f002 cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       93ed2f30 cit fix for conpair
       9672af28 sym link fixes
       c4d0581e cit fixes to tumor pair
       41ccead3 fixes to dnaseq
       b07afef9 updates to GRCh38 annotation file, module updates, cit fixes
       066fb499 fixes to deliverable and b38 ini
       034447e9 updates to cit and fixes to one job mpileup steps
       17d4c791 updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       4668a0bc major fixes to deliverables and completion of beluga test
       9ffad226 fixes to bcftools mpileup for dnaseq mpileup protocol
       180ae35c hs37d5/GRCh37 fixes for purple
       87dcc11f fastpass varscan merge fix
       f6627158 cit fixes to fastpass nb_job=1 for panel, soft include of split/scatter intervals and GenomicsDBImport for dnaseq
       233bc03b cit fix for conpair
       8084b8f7 sym link fixes
       c5840c6b cit fixes to tumor pair
       e4193cae fixes to dnaseq
       14aa4a30 updates to GRCh38 annotation file, module updates, cit fixes
       de7ed14d fixes to deliverable and b38 ini
       95100617 updates to cit and fixes to one job mpileup steps
       510efd1f updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       cad6781f major fixes to deliverables and completion of beluga test

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      43 commits

       05d65f79 dnaseq sym_link fix
       c49e4692 dnaseq ini fixes
       6b63aea5 dnaseq beluga and cedar ini fixes
       04007457 sequenza bug fix
       362539f0 revert back to old sequenza-utils
       3714cd06 fixes to sequenza for exomes, and minor SV fixes
       444e6c26 exome cit fix to sequenza
       9b97f7e8 converted -n to -c
       bf70c513 cit fixes - runtimes
       bfd86d95 config file fixes
       75fd2198 further sanity-check fixes
       30c20126 sanity check fixes
       adb0724e multiple pair same control fixes and vardict exome fixes
       f0c3f899 cit fixes to vardict
       34b052de cit fixes for tumor exome
       21e17a24 better cram compression of base recalibrated bams with variantBam
       baf4196f sequenza and germline ensemble fixes
       b736a4cf sambamba mark_dup fix
       667a4fe7 module fixes
       cee17748 sym link dnaseq cit to tumor_pair cit
       1e7e2e35 tumor pair cit updates
       8d35d7f6 dnaseq sym_link fix
       7da9cee8 dnaseq ini fixes
       3d911187 dnaseq beluga and cedar ini fixes
       e1a83a36 sequenza bug fix
       f494d73c revert back to old sequenza-utils
       ee70ffc8 fixes to sequenza for exomes, and minor SV fixes
       9ffeafa3 exome cit fix to sequenza
       cb1a59cf converted -n to -c
       00943b32 cit fixes - runtimes
       29c33947 config file fixes
       d3eedbbf tumor base fix
       bed7a1d4 further sanity-check fixes
       d36c1e8e sanity check fixes
       e901171a multiple pair same control fixes and vardict exome fixes
       165e341b cit fixes to vardict
       10aa2b10 cit fixes for tumor exome
       ceef654e better cram compression of base recalibrated bams with variantBam
       3d476d22 sequenza and germline ensemble fixes
       1edc6d7f sambamba mark_dup fix
       a3b1208f module fixes
       caabe630 sym link dnaseq cit to tumor_pair cit
       b9ed61eb tumor pair cit updates

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      8 commits

       e7bd7451 code cleaning and fixes to exome interval list
       a7d64071 fixes to cedar ini
       cca1a4c6 fixes to symlinks for paired indel realignment
       211a29e2 cedar ini and exome update
       f65b070d code cleaning and fixes to exome interval list
       febded3f fixes to cedar ini
       43d40444 fixes to symlinks for paired indel realignment
       09c93a52 cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      4 commits

       c7f3fa9e cedar fixes and GRCh38 fixes
       2baaf564 cedar germline sv updates
       7367777a cedar fixes and GRCh38 fixes
       bd82ef13 cedar germline sv updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      46 commits

       21b88961 fixes to sambamba for lumpy - query sorting instead of corrdinate
       dc9750df updates to metasv - somatic
       c19ad156 Fixes and updates to reference files
       08f00272 remove testing steps
       aad782e6 updates to SV germline and reference genome tweaks
       491914df somatic sv fixes: lumpy and svaba annotations
       d71a64de json fix
       03bcbc69 fix to germline SV: breakseq2
       97a4bc18 deleting hidden files
       b24f8f3d merge fixes
       9052826c GATK4 fixes - bam indexing and markDupSpark
       fc33ed9b bcftools fixes for tumor pair
       0deb8a94 fingerprint and bug fixes
       57f218ae dnaseq - vcftools qc addition: --missing_indv and --depth
       5c82925c GRCh38 fixes
       bb43292f select input for variant caller and fixes to one job calling
       b87e9f50 json and folder updates
       f15fd740 fixes to sCNAphase
       34c9e47a Bug fixes prior to json additions
       ae851f40 merging snv and sv, adding protocols
       c82ce685 Updates and debug
       57883607 Add set somatic and actionable mutations
       bc0ca0ac added multiqc and other tweaks
       c51db536 fixes to sambamba for lumpy - query sorting instead of corrdinate
       8da9446d updates to metasv - somatic
       be214eb3 Fixes and updates to reference files
       f717bb79 remove testing steps
       4c642237 updates to SV germline and reference genome tweaks
       541f6715 somatic sv fixes: lumpy and svaba annotations
       c795a24e json fix
       8cc8e79b fix to germline SV: breakseq2
       628ce3d2 deleting hidden files
       2a83ec4b merge fixes
       90ce1f49 GATK4 fixes - bam indexing and markDupSpark
       4e8c7501 bcftools fixes for tumor pair
       4c874cc4 fingerprint and bug fixes
       3cec8c22 dnaseq - vcftools qc addition: --missing_indv and --depth
       6f3ae37a GRCh38 fixes
       ebfc7643 select input for variant caller and fixes to one job calling
       4e0b29b5 json and folder updates
       fc717428 fixes to sCNAphase
       a8034e24 Bug fixes prior to json additions
       2b232671 merging snv and sv, adding protocols
       3c8b861c Updates and debug
       92db2299 Add set somatic and actionable mutations
       39425ef1 added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      24 commits

       3e96d818 dnaseq - fixes for VariantBam
       2cd2eb5d fixes to purple and mkdir vardict fix
       3ba332b8 sequenza WGS fixes
       86f1fa32 fixes to metasv for tumor pair
       3d446d9c fixes to sv germline calls
       6f26955c Single job bug fixes
       893626b5 manta I/O fix and other bug fixes
       9f4109b5 config updates and b38DH added
       9098c409 dnaseq germline SV updates
       8d737df1 gatk4 updates and bug fixes
       843ef92b Json related bug fixes
       54fdadc7 Bug fixes and modification derived from initial PROFYLE benchmarking
       5cd67334 dnaseq - fixes for VariantBam
       047e9613 fixes to purple and mkdir vardict fix
       776c53a9 sequenza WGS fixes
       dc8f8202 fixes to metasv for tumor pair
       81b508ba fixes to sv germline calls
       46528cd5 Single job bug fixes
       be5f8d27 manta I/O fix and other bug fixes
       fb73873f config updates and b38DH added
       9c2a8b0b dnaseq germline SV updates
       e86e74ae gatk4 updates and bug fixes
       e9b3a2f8 Json related bug fixes
       0201d386 Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       e9cc4056 seqeuenza fixes and sv fixes
       3786a036 gatk4 mutect2 updates
       4ae76b75 cit-based fixes to NGScheckmate
       58c3db32 dnaseq qc additions: NGScheckmate and peddy
       9ae42103 seqeuenza fixes and sv fixes
       dc793204 gatk4 mutect2 updates
       0f9e09a5 cit-based fixes to NGScheckmate
       45837b45 dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      6 commits

       f00e0b39 Merged in dev_reveleig (pull request #252)
       00fb9e97 Merged in rebasing_tp (pull request #246)
       08361bc0 Debugging and Guillimin specfic fixes
       f4699979 updates to config
       07d625a8 Debugging and Guillimin specfic fixes
       87f1ef41 updates to config

  Shaloo Shalini <shaloo.shalini@gmail.com>      1 commits

       78f59d52 Merged in ss_wf_chipseq_93 (pull request #244)

  shaloo <shalz@hotmail.com>      1 commits

       8c742ce2 Updates workflow for chipseq -t chipseq case with differential binding step and dependency update Fixes #93

3.5.0        Mon Jul 12 16:41:20 2021 +0000        1071 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      8 commits

       fd8c1cd5 GenPipes - Release 3.5.0 : correcting typo in chipseq README
       038dcb2c GenPipes - Release 3.5.0 : correcting typos in READMES
       6e4f4de1 Version bump to 3.5.0
       48f97ec1 GenPipes - Release 3.5.0 : correcting typos in pipeline README files
       f892b2e9 GenPipes - Release 3.5.0 : updating VERSION and README-release files
       82b28e31 GenPipes - Release 3.5.0 : corrected typo in README
       0cf907cd GenPipes - Release 3.5.0 : updating version in the pipeline README files
       ddb421a4 Merge branch 'dev' into release_3.5

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      6 commits

       08bd50a8 GenPipes - Resources : updated picard index command in genome installation script to run command on a compute node instead of on the login node
       02124fc0 correcting typo genome install script
       ae646a79 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       69818199 GenPipes - Resources : updates in software and reference install scripts
       0659d35b GenPipes - Resources : adding software installation scripts
       53dcb830 GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      47 commits

       339172cf Merged in release_3.5 (pull request #242)
       525bfe10 Merged in release_3.5 (pull request #241)
       0d0f7b40 Merged in bash_cmd_for_covseq_update_eh (pull request #236)
       5a54cf99 GenPipes - CoVSeq pipeline : updated call to "zcat" in covseq.py
       2a056bd5 GenPipes - BFX : updated bash_cmd.py
       2809e5ef GenPipes - Install script : added Nanopore pipeline to PATH, fixing Issue #74
       daebeb7a Merged in ehenrion/changelogmd-edited-online-with-bitbucket-1620161608327 (pull request #220)
       210d516a Version bump to 3.4.0
       a36b2cbb GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       7835760c GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8e2018d4 GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       b4cc95ec GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       f4534a4a GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       e2193e68 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       4b416a0d GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b34a1a56 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       05a90a65 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       5748a6be GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       eb4a375b GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       c1a8b6bb VERSION edited online with Bitbucket
       4a0537a8 GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       2fd4f5ea GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       40125ca0 GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       30917fea GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       810feb8c GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       8911666d GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       3d76f2c1 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       b2b9c458 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       929424c5 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       879f109c GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       98c77304 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ad82c94e GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       31968eea GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       8dff24af GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       dbbaa112 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       4b617ba2 GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       04f7714a GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       f1a49af8 GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       f9c2aefc GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       72cd3e9d GernPipes - DNASeq : corrected merge_filter_bcf outputs
       3b635b23 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       04b24d58 GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       55ab0170 GenPipes - DNASeq SV : fixing delly input error
       09c1bad3 GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       a756e77b dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       b514d9f0 GenPipes - Bug fix : corrected dnaseq.cedar.ini
       26b421e5 GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       3297e974 Fix Issue #86
       ff3110dc Quick correction to address misleading headcrop parameter in the RNAseq base ini.
       e3d7b04f rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      344 commits

       ee6fd308 Merged in covseq_stretenp (pull request #239)
       76bfba93 Fixing sambamba version
       4bc5b997 hic code improve
       2f1529ee Debug pdf rendering
       b76a994c Fixing ini
       7a6e9c1d Merge branch 'dev' into covseq_stretenp
       985a5c5f Increasing resources for ncovtools and adding run_name for cit
       83caadfe Upgrading bcftools
       cdbb4141 Fixing merge_hicrep_scores
       b7ecfc42 Switching to samtools 1.12
       8b29ef3a Fixing merge_hicrep_scores
       beab8662 Samtools version update in all pipelines
       00509e83 Merged in covseq_stretenp (pull request #237)
       74a12033 Merged dev into covseq_stretenp
       63cb7a2f Fixing indentation error
       22c619af Merged in covseq_stretenp (pull request #235)
       16854eb7 Switching to released version of covseqtools
       a69ddeb3 Cleaning comments
       4d85dff4 Removing guillimin and mammouth inis
       f25cf786 Merged dev into covseq_stretenp
       e955aa2d Merged in chipseq_design_change (pull request #234)
       5ca220d8 Merged dev into chipseq_design_change
       4c4678bb Removing name in report/metrics job name
       0b1c5622 Debug
       8c659793 Debug
       13f9971d Debug
       9cb76693 Merged dev into chipseq_design_change
       bb6e4210 Merged dev into covseq_stretenp
       a1650eb4 Forcing ncovtools execution
       5f908a77 Fixing metadaya ncovtools file creation
       debf790c Adding filtered bam to ncovtools
       13b80b05 Fixing ncovtools
       296dfea8 Fixing peak file naming
       0133eb08 Forcing snakemake execution
       9e3647ed Forcing snakemake execution
       0b32f4bd Forcing snakemake execution
       c85a9fe4 Debug
       a562dd2c Changing resources for ncovtools
       59fdae5b Debug
       2e4dcbef Adding dependencies for reporting job
       558cebc4 Purging modules
       2147b859 Adding more resources for ncovtools
       2366c903 Debug
       355958da Debug
       4c6955bb Debug
       dff338fe Debug
       acca4096 Debug
       e18a0b7d Debug
       c97470b5 Debug
       7f82b371 Debug
       4ee91902 Debug
       a50112bd Debug
       23233fb2 Debug
       e143ff19 Debug
       9c6db14c Debug
       98c02317 Debug
       4e8b8e6a Debug
       290e9abd Debug
       032d4258 Debug
       72cc7860 Adding GenPipes version accseible inside GenPipes
       49210d2a Degrouping ncvotools from report generation because of cd for ncvotools
       f8846315 Adding modules reading for reporting
       a3209fda Debug
       4717b023 Debug
       bb6b07c7 Debug
       8587d125 Debug
       612afbf7 Debug
       f7e17ec2 Debug
       b9de420c Debug
       2fde5e2f Debug
       7a124f18 Debug
       feb23641 Debug
       d15188b2 Debug regex prepare ncovtools
       de760753 Debug
       4b83dab5 Debug
       8b1e5917 Reporting init ncovtools yaml
       b226f89c Debug finding files for report
       bc78d4ae Debug test dependencies reporting
       d519ac78 Adding missing sections to beluga ini
       60638b96 Debug
       f1109187 Debug
       e482afa1 Debug
       d77e3166 Debug
       11ff7bce First commit on gathering metrics for report
       51df0684 Merged dev into covseq_stretenp
       58fd24f4 Merge
       cc0c7615 Fixing minor bug and typo
       2376a973 Update READMEs
       83015311 Updating inis and fixing homer output folder location
       706ccb32 Debug
       c59a207a Debug
       bf45b1d1 Debug
       7b30adfb Debug
       e79e74ac Debug
       c78d0ea5 started to edit the hicseq.py added new file for hicrep
       090de772 Fixing minor bug and typo
       01d52a21 Fixing multiqc dependencies
       8d17dce0 Debug
       f0a2d017 Debug
       9dc813e0 Debug
       a3f06f2b Debug
       f8dd6516 Debug
       d831843c Debug
       9f8d02b3 Debug
       39463cdb Debug
       daebd07f Debug
       898a5230 Debug
       3bdc004a Debug
       333524bd Debug
       c017225d Debug
       a75dc5ab Changing name of header for a report table
       1eb26a73 Debug
       9eb8c957 Typo
       cbd8f1e0 Addinx x permission to job2json.py
       6f6a3817 Fixing report naming too long
       9c91a7e3 Fixing report naming too long
       d95ab926 Fixing report naming too long
       ce26275f Fixing mugqicValidator.py for new readset format (chipseq only)
       27598c0d Switching to latest mugqic_tools release
       fe2ab509 Changing Ressources
       74010231 Changing Ressources
       7e112120 Changing Ressources
       0a52ab63 Changing module versions
       3b027e9a Debug
       e5de8a76 Fix
       d1304d01 Increasing picard collect multiple metrics RAM
       fff64384 Debug
       8a0463be Debug
       cb4a7e3e Fix ressources
       268550ba Fixing annotation peaks
       305bec91 Fixing annotation peaks
       7c56538a Fixing annotation peaks
       94090176 Adding annotation peaks merge for Narrow peaks only
       8420a20f Iterating only on Narrow peaks for annotation
       2d4e2e77 Iterating only on Narrow peaks for annotation
       366142d4 Debug report
       9253d109 Debug
       aac7c9cc Debug
       9c6dc244 Debug
       c9b7191b Debug
       4ac6ca3f Renaming IHEC metrics section
       36880b94 Debug
       c3858f71 Debug
       b1b68b91 Debug Report
       08b77c98 Debug Report
       916d53f9 Debug
       745fbda1 Debug
       0c8164f0 Fixing report
       c918cb1e Debug
       f578ce16 Debug
       ec566be3 Debug
       a3b196de Debug
       7768187c Debug
       ddf75886 Debug
       4944d36d Debug
       e7d8d28f Debug
       b9f35c3b Debug
       847dddcd Debug
       66a5ca17 Debug
       93ad22ea Debug
       9d3cdd81 Changing Output metrics
       723f9ed0 Debug
       fbf10b92 Debug
       a68f0a65 Debug
       68afdd3e Debug
       1ac59b20 Debug
       edcfbbfd Debug
       d5548d1b Debug
       a3be70d0 Debug
       3b9d879f Debug
       3b5d819b Debug
       9215f338 Debug
       1d7732f9 Debug
       096571af Debug
       1e22dc82 Debug
       28bbc7f2 Changing macs2 atacseq param & add macs2 param
       982fb12b Debug
       a86d776b Debug
       9a86cbd7 Debug
       5c6e6b4d Debug
       f778cc3b Debug
       05187de2 Debug
       2732d265 Debug
       54f47a2d Debug
       46922378 Debug
       187edbe2 Debug
       6392d86b Debug
       493624f2 Debug
       553b9baf Debug
       eca0813f Debug
       721c8d10 Debug
       619d6e4c Debug
       773569d7 Debug
       7836d837 Debug
       8ab15788 Debug
       db9c7868 Debug
       2e8f5d2f Changing R version for mugqic_tools R scripts
       c90cec99 homer_annotate_peaks
       d4586662 qc_metrics report
       a4c0fbfa Fix ihec_metrics outputs
       66af5287 Fix ihec_metrics outputs
       2bca4810 Fix ihec_metrics outputs
       d33b9e23 Fixing MultiQC report
       cab6ca63 Fixing MultiQC report
       fe600e9a Fixing MultiQC report
       5754975f Fixing MultiQC report
       f31bf711 Fixing MultiQC report
       fcdc378b Fixing MultiQC report
       d3559e6f Improving MultiQC report
       f424ffac Debug
       a1dcb44a Debug
       44b2647f Debug
       fc188191 Debug
       777c716c Debug
       1217e0b6 Debug
       cbfcd261 Debug
       a5a3ca35 Debug
       479504f5 Debug
       44a9b808 Debug
       425dc3cd Debug
       a853292f Debug
       eaec07b0 Debug
       f80b4f2e Debug
       81f990a4 Debug
       ecc34289 Debug
       89b57b8a Debug
       36e4a2af Debug
       bbf0fdc1 Debug
       d4d8d648 Debug
       33521807 Debug
       98500465 Debug
       1f291caf Major changes IHEC metrics test
       047fe5f3 Major changes IHEC metrics test
       92407469 Major changes test
       23c92564 Macs2 report changes
       acfc7822 Major change test
       76bdd83e Major change test
       392234ff Major change test
       a27e5b69 debug
       d54a5903 debug
       4ab48dbc debug
       cfaad212 debug
       5daadd72 debug
       4ebada3f debug
       5a16691c debug
       807e2645 debug
       fd604bd8 debug
       c6347e79 debug
       aa9e88a5 debug
       440c2904 debug
       96d01f10 debug
       e5253567 debug
       7785d116 debug
       dc09b154 debug
       a5c60c6c debug
       7be46d4d debug
       0faaea9a debug
       4db52480 debug
       74b73f35 debug
       2df9a1c5 debug
       6b583d5a debug
       e3961cb8 debug
       622832e3 debug
       8b2d57a1 debug
       fa973d0e debug
       64fd4ed9 debug
       6a9071f5 debug
       ab1f0240 debug
       40cb701b Fix test
       f1dce58b Major readset change test
       349a1a3d Fix
       095b13e2 Fix
       17519c7f Fix
       0baff781 Filtering after merging and changing naming to fit with other pipelines
       cab4a18b Increasing default resources and adding options for markdup
       ffae951c Fixing beluga ini
       33c1ab90 Switching to sambamba markdup for IHEC
       b689ca7f Fix
       9a59883e Fix
       c23b1938 Fix
       9d1aa0a6 Fix
       fb2ad2a9 Fix
       879f850a Fix
       8cc4c624 Options becomes optional
       7cd7a814 Fix
       804de5cc Fixing typo
       9795af9d Fix
       11cd56d4 Fixing sambamba merge
       fb4bbcfa Typo
       823693c4 Adding mkdir
       11e06a2f Fixing temp typo to tmp
       27386371 Fixing job
       6fe44469 Fixing bash import
       8fbe44ae Fixing minor bug and typo
       d9dbf29a Adding missing sections to beluga ini
       5f74c2f6 Renaming freebayes consensus and ncovtools quick_align
       4e6a1bdc Degrouping rename consensus jobs because of quast issue creating infinite dependency loop
       25f9e362 Debug
       a0ed01a7 Debug
       07ab3d11 Adding missing sections to beluga ini
       cc21ed98 Debug
       a113c663 Merged dev into mgi_stretenp
       490857bf hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       c844ceaf hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       4a793a05 started to edit the hicseq.py added new file for hicrep
       0e991a7d started to edit the hicseq.py added new file for hicrep
       72c9387c hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9458e13d hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ac314f69 hicseq completed adding basic features of the hicrep analysis.
       cd8fbb9e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       086a8cbd hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       f9893e59 hicseq completed adding basic features of the hicrep analysis.
       4c811af2 Cleaning
       e87bf289 Debug
       a360050d Debug
       d9559252 Changing freebayes calling to skip bgziping
       cdee3523 Reducing Ressources
       8856ee9b Debug
       e77bd4ba Debug
       30e04a4e Debug
       670c71e5 Debug
       292c0734 Fixing Jared's script
       d799b7a2 Reducing Ressources
       e2bb1dbd Fixing Jared's script
       e6f4df81 Fixing freebayes variant calling
       a7e2b43e Fixing Jared's script
       7a6ef999 Fixing bcftools consensus
       5a88dcd5 Fixing freebayes
       114d9bba Fixing freebayes
       1a3db5e3 Fixing freebayes
       46d28391 Fixing freebayes/bcftools
       612e5068 Adding bcftools consensus creation following Jared's script
       92df26f2 Adding bcftools consensus creation following Jared's script
       e73f67df Adding freebayes variant calling
       5ce0e180 Grouping rename jobs into 1 job
       ccae3393 Grouping rename jobs into 1 job
       f3daf102 Grouping rename jobs into 1 job
       6d850fb5 Adding freebayes variant calling
       edf5a762 Grouping rename jobs into 1 job
       76fbe6bf Adding freebayes variant calling
       7ab5da22 Adding freebayes variant calling
       00acf27b Adding freebayes variant calling
       c700579b Adding freebayes variant calling
       eefed60d Fixing flagging bug and updatinh year to 2021

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      3 commits

       1719a1f0 Merged in update_container (pull request #238)
       906d0fd3 Merged in slurm_cgroup_problem (pull request #233)
       1dd042cd Merged in temp (pull request #232)

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      387 commits

       773c8d4b hicseq completed adding basic features of the hicrep analysis.
       ffe7becc Completed developing hicrep and quasar analysis
       9e790f67 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6e6b3659 Completed hicrep analysis except the creation of the graph
       dbb237f8 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       72bde2f1 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       167673f5 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2471a46c hicseq completed adding basic features of the hicrep analysis.
       68aece4f hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       5903eca7 Added pairwise combination for samples
       e7322592 started to edit the hicseq.py added new file for hicrep
       af7bfd46 hicseq completed adding basic features of the hicrep analysis.
       ad7a803e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       23746b5f Added pairwise combination for samples
       e50669f3 started to edit the hicseq.py added new file for hicrep
       96ea65eb hicseq [hicrep.py] - Corrected typo
       a952c1f0 hicseq [hicrep.py] - corrected R_TOOLS path
       3dee5c65 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2e52c82a hicseq [hicrep.py] - Corrected typo
       fc9ccc44 hicseq [hicrep.py] - corrected R_TOOLS path
       e93c722d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       50760457 hicseq [hicrep.py] - Corrected typo
       1404bb4b hicseq [hicrep.py] - corrected R_TOOLS path
       70a9c006 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5e3e5ffc hicseq [hicseq.py] - corrected file after rebase
       b96c0e71 hicseq [hicrep.py] - Corrected typo
       dc6f85b8 hicseq [hicrep.py] - corrected R_TOOLS path
       750cc6c1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       92263db4 hicseq [hicseq.py] - corrected file after rebase
       b708ea3c hicseq [hicrep.py] - Corrected typo
       b209f6a1 hicseq [hicrep.py] - corrected R_TOOLS path
       3ba0747f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2bf8e950 hicseq [hicrep.py] - Corrected typo
       8936ce6d hicseq [hicrep.py] - corrected R_TOOLS path
       9ac5f101 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1a5bc63c hicseq [hicrep.py] - Corrected typo
       c01ee521 hicseq [hicrep.py] - corrected R_TOOLS path
       1805da83 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       19df2bab hicseq [hicrep.py] - Corrected typo
       33e3949d hicseq [hicrep.py] - corrected R_TOOLS path
       77ebe87e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       85294d65 hicseq [hicrep.py] - Corrected typo
       64840a7c hicseq [hicrep.py] - corrected R_TOOLS path
       3013d61d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       734c2793 hicseq [hicrep.py] - Corrected typo
       9a33ab8a hicseq [hicrep.py] - corrected R_TOOLS path
       83306119 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       94b562bf hicseq [hicrep.py] - Corrected typo
       005b6cc7 hicseq [hicrep.py] - corrected R_TOOLS path
       28076a5c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       858b0570 Added pairwise combination for samples
       e76217b2 started to edit the hicseq.py added new file for hicrep
       8e835eb8 hicseq [hicseq.py, readme.md] - modified readmefile
       6bf04f8d hicseq [hicseq.py] - corrected file after rebase
       589ae72f hicseq [hicrep.py] - Corrected typo
       33129961 hicseq [hicrep.py] - corrected R_TOOLS path
       2c6c2116 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       54108f49 hicseq [hicrep.py] - Corrected typo
       ca68cf40 hicseq [hicrep.py] - corrected R_TOOLS path
       c1d6c9c0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       01817b64 hicseq [hicrep.py] - Corrected typo
       a68d82b4 hicseq [hicrep.py] - corrected R_TOOLS path
       310b5eda [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       785bd9db hicseq [hicrep.py] - Corrected typo
       4db60e27 hicseq [hicrep.py] - corrected R_TOOLS path
       0baffa5c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fc107ccc hicseq [hicrep.py] - Corrected typo
       8dbe281d hicseq [hicrep.py] - corrected R_TOOLS path
       f2721650 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f57e58a hicseq [hicrep.py] - Corrected typo
       11c5af45 hicseq [hicrep.py] - corrected R_TOOLS path
       e544d4f3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0c90edc1 hicseq [hicrep.py] - Corrected typo
       ec0feb84 hicseq [hicrep.py] - corrected R_TOOLS path
       165aeeb4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d4d21a3f hicseq [hicrep.py] - Corrected typo
       f665976e hicseq [hicrep.py] - corrected R_TOOLS path
       9a1295e5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9778a4a8 hicseq [hicseq.py] - corrected file after rebase
       10365ba0 hicseq [hicrep.py] - Corrected typo
       f470412c hicseq [hicrep.py] - corrected R_TOOLS path
       a7f0ba32 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0e1ddc65 hicseq [hicrep.py] - Corrected typo
       956a72de hicseq [hicrep.py] - corrected R_TOOLS path
       1cc54983 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9b7e8b32 hicseq [hicrep.py] - Corrected typo
       06cdfcf4 hicseq [hicrep.py] - corrected R_TOOLS path
       d4c14d23 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0ce13dfe hicseq [hicrep.py] - Corrected typo
       279dc529 hicseq [hicrep.py] - corrected R_TOOLS path
       7c304c72 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       81ad33d2 hicseq [hicrep.py] - Corrected typo
       87bbe668 hicseq [hicrep.py] - corrected R_TOOLS path
       0c5127b1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b67b6b04 hicseq [hicrep.py] - Corrected typo
       95a4a976 hicseq [hicrep.py] - corrected R_TOOLS path
       4c25fc25 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6d7a744c hicseq [hicrep.py] - Corrected typo
       d96ad6fd hicseq [hicrep.py] - corrected R_TOOLS path
       a1768c57 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c52093e Added pairwise combination for samples
       5aff10cf started to edit the hicseq.py added new file for hicrep
       55987dd4 hicseq pipeline [changed the step order]
       525c7d39 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       f1811aa3 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       5fe3a493 hicseq [hicseq.py, readme.md] - modified readmefile
       df1ac4b2 hicseq [hicseq.py] - corrected file after rebase
       cfd6f524 hicseq [hicrep.py] - Corrected typo
       2b5daaca hicseq [hicrep.py] - corrected R_TOOLS path
       796fe27b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e9d2d2db hicseq [hicrep.py] - Corrected typo
       34b96297 hicseq [hicrep.py] - corrected R_TOOLS path
       beb3b076 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       af6d8aca hicseq [hicrep.py] - Corrected typo
       c97468af hicseq [hicrep.py] - corrected R_TOOLS path
       658e574d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5b2fedf9 hicseq [hicrep.py] - Corrected typo
       f3c653f2 hicseq [hicrep.py] - corrected R_TOOLS path
       1690ce43 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       89371a2d hicseq [hicrep.py] - Corrected typo
       ef565948 hicseq [hicrep.py] - corrected R_TOOLS path
       ff4b8ceb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       161296e2 hicseq [hicrep.py] - Corrected typo
       61e8b115 hicseq [hicrep.py] - corrected R_TOOLS path
       a147aa98 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       db542ae8 hicseq [hicrep.py] - Corrected typo
       229f0656 hicseq [hicrep.py] - corrected R_TOOLS path
       919463a2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       077e2ba0 hicseq [hicrep.py] - Corrected typo
       8c2d82be hicseq [hicrep.py] - corrected R_TOOLS path
       7dd995d6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       78ee1f78 hicseq [hicseq.py] - corrected file after rebase
       42120ffd hicseq [hicrep.py] - Corrected typo
       e21766cd hicseq [hicrep.py] - corrected R_TOOLS path
       2aebe64d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       51e0f2aa hicseq [hicrep.py] - Corrected typo
       3d8ad797 hicseq [hicrep.py] - corrected R_TOOLS path
       59ea5d00 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       20ea3a82 hicseq [hicrep.py] - Corrected typo
       df2b34a3 hicseq [hicrep.py] - corrected R_TOOLS path
       6bdaef1f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c4727f54 hicseq [hicrep.py] - Corrected typo
       bd80a99e hicseq [hicrep.py] - corrected R_TOOLS path
       7d004943 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2c056190 hicseq [hicrep.py] - Corrected typo
       593a82da hicseq [hicrep.py] - corrected R_TOOLS path
       60581eb3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2f0362a1 hicseq [hicrep.py] - Corrected typo
       a2dfa4f4 hicseq [hicrep.py] - corrected R_TOOLS path
       04faafc2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       85569620 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       57dba1bd hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       61fce92a hicseq [quasar_qc - corrected module loading in matrix restructuring
       6bc2c0ce hicseq [hicrep.py] - Corrected typo
       23217cac hicseq [hicrep.py] - corrected R_TOOLS path
       fc77bfe8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8190866b hicseq [hicrep.py] - Corrected typo
       70897943 hicseq [hicrep.py] - corrected R_TOOLS path
       e5be8fdd [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       70754df0 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       a83bbe3f [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5b5ac2f3 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b9a9dfd [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       29fefa5e [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9c80cb01 Completed developing hicrep and quasar analysis
       14a6698e [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d4f176ef [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       69310a4f created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       38d65fda Completed hicrep analysis except the creation of the graph
       97d73811 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       a7272ddb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       99142a9a hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3c16aa5c hicseq completed adding basic features of the hicrep analysis.
       366611dc hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       a26866b4 Added pairwise combination for samples
       54e99f84 started to edit the hicseq.py added new file for hicrep
       ad051374 hicseq [hicseq.py, readme.md] - modified readmefile
       57f13dfb hicseq [hicseq.py] - corrected file after rebase
       727d579f hicseq [hicrep.py] - Corrected typo
       cabd6b11 hicseq [hicrep.py] - corrected R_TOOLS path
       4f95216f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d9d53d59 hicseq [hicrep.py] - Corrected typo
       2d50a471 hicseq [hicrep.py] - corrected R_TOOLS path
       0b6541ef [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b56a96a8 hicseq [hicrep.py] - Corrected typo
       a41dc59d hicseq [hicrep.py] - corrected R_TOOLS path
       d6a751fc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       146c8921 hicseq [hicrep.py] - Corrected typo
       552d8f97 hicseq [hicrep.py] - corrected R_TOOLS path
       e529d26b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       743948f1 hicseq [hicrep.py] - Corrected typo
       bed01546 hicseq [hicrep.py] - corrected R_TOOLS path
       c334f10b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b286fcd5 hicseq [hicrep.py] - Corrected typo
       548cfc05 hicseq [hicrep.py] - corrected R_TOOLS path
       316c44f7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       95bb7ed9 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2cd7fd62 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       dad2a18f hicseq [quasar_qc - corrected module loading in matrix restructuring
       912cca58 hicseq [hicrep.py] - Corrected typo
       6f1f5207 hicseq [hicrep.py] - corrected R_TOOLS path
       8120389b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3816804a hicseq [hicrep.py] - Corrected typo
       8581ac0b hicseq [hicrep.py] - corrected R_TOOLS path
       5371ed0e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       226a4f75 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       be5df368 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       f5876b20 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4465ebe0 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c5c922fe [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       24c34f0f Completed developing hicrep and quasar analysis
       c7d0aa89 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       72282088 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d1c440a6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       b4bc9bfe Completed hicrep analysis except the creation of the graph
       324fbc2f hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       ff16ded0 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       3d54c9fd hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1e3be791 hicseq completed adding basic features of the hicrep analysis.
       e35ac3d8 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ac2c2082 Added pairwise combination for samples
       61bd4a7a started to edit the hicseq.py added new file for hicrep
       c950a435 hicseq [hicseq.py] - corrected file after rebase
       9d755feb hicseq [hicrep.py] - Corrected typo
       9cfb4eef hicseq [hicrep.py] - corrected R_TOOLS path
       9bb97058 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b744cede hicseq [hicrep.py] - Corrected typo
       4477d467 hicseq [hicrep.py] - corrected R_TOOLS path
       74fafa09 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2a044933 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       acb090f5 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       a7578f99 hicseq [quasar_qc - corrected module loading in matrix restructuring
       0e8398af hicseq [hicrep.py] - Corrected typo
       9d27afa6 hicseq [hicrep.py] - corrected R_TOOLS path
       7fdc12cf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8448782a hicseq [hicrep.py] - Corrected typo
       09b74208 hicseq [hicrep.py] - corrected R_TOOLS path
       e5f68499 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cfdb9d46 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5af4d805 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       1ae80800 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e48da77f [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       46ba5a85 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       8a20dbad Completed developing hicrep and quasar analysis
       3e6f7634 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b9ae8485 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       60e5df89 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3c2f20 Completed hicrep analysis except the creation of the graph
       e1dea78c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b086ca53 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       114fbd17 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       edd76c8c hicseq completed adding basic features of the hicrep analysis.
       3549fd07 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e3720267 Added pairwise combination for samples
       8d80a5b0 started to edit the hicseq.py added new file for hicrep
       18f4e153 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4fbc1fdc hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       c12f62a8 hicseq [quasar_qc - corrected module loading in matrix restructuring
       5e9ddcd4 hicseq [hicrep.py] - Corrected typo
       19479308 hicseq [hicrep.py] - corrected R_TOOLS path
       c71ed68b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       157c1b51 hicseq [hicrep.py] - Corrected typo
       5937a022 hicseq [hicrep.py] - corrected R_TOOLS path
       b99ff65f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71c1160f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       03f2537a [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c9403408 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a8e16f94 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       298874c4 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       bf24afb6 Completed developing hicrep and quasar analysis
       e5fe17ff [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e2fc3604 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f67be0b6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f1b5f5ad Completed hicrep analysis except the creation of the graph
       619bbf16 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       c74201f5 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0b72f780 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f998fcf hicseq completed adding basic features of the hicrep analysis.
       a06d7268 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       009d2868 Added pairwise combination for samples
       9bc074eb started to edit the hicseq.py added new file for hicrep
       bcc23d6b hicseq [base.ini] - updated mugqic tools version to 2.3.1
       20949c26 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       81afdcce hicseq [quasar_qc - corrected module loading in matrix restructuring
       08c67579 hicseq [hicrep.py] - Corrected typo
       51c8f664 hicseq [hicrep.py] - corrected R_TOOLS path
       9c2ee30d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       adc92e7a hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b020b249 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       623d8bbc [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       7e1ecf64 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       6e22be6d [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f2815e41 Completed developing hicrep and quasar analysis
       ce39693c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       f2e21fcf [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4895f68c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       7a13a038 Completed hicrep analysis except the creation of the graph
       0288fb9a hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       708eeec3 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e97faa54 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1ae27685 hicseq completed adding basic features of the hicrep analysis.
       646858db hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dab3e7d4 Added pairwise combination for samples
       1e9e150b started to edit the hicseq.py added new file for hicrep
       dd821443 hicseq [hicrep.py] - Corrected typo
       0c8400be hicseq [quasar_qc.py] - Corrected deletion by mistake
       18a704df hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       18747b85 hicseq [hicrep.py] - corrected R_TOOLS path
       5759e582 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cc9a3184 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       dcfaa23c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       df458a53 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       6718520d [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       ab0c23aa [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       63ae5faa Completed developing hicrep and quasar analysis
       506c95bb [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       78c57d89 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b87ff3f9 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       18258411 Completed hicrep analysis except the creation of the graph
       643905bd hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       aa8395fb hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c2851616 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f10a6b3 hicseq completed adding basic features of the hicrep analysis.
       662a52c7 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d43a3c19 Added pairwise combination for samples
       3472e7a5 started to edit the hicseq.py added new file for hicrep
       4a08538b hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       e2ede562 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97d248a1 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       bb6887e8 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       815e22a6 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       854c1c0d Completed developing hicrep and quasar analysis
       c648da3c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0f850588 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54cefb76 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       fb3912eb Completed hicrep analysis except the creation of the graph
       c1e10353 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       dc414bd2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       943a7c18 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       6be727c2 hicseq completed adding basic features of the hicrep analysis.
       1e5d4446 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       3af40a81 Added pairwise combination for samples
       3180796f started to edit the hicseq.py added new file for hicrep
       e67b0797 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6e359557 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f0dd7199 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       9556bef7 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f17594ed Completed developing hicrep and quasar analysis
       a384d9a5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       339a8017 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       aefd4bbd created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       98a1cac6 Completed hicrep analysis except the creation of the graph
       2ae35431 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       255526b8 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ced9d138 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       48a96471 hicseq completed adding basic features of the hicrep analysis.
       3b2f387e hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       77e52e21 Added pairwise combination for samples
       0894881c started to edit the hicseq.py added new file for hicrep
       ffb9b96a [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7507b5a8 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       21f3a9cb [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       2c4d8bbe [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       72686d18 Completed developing hicrep and quasar analysis
       7cb650ca [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       abd1855f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       045f8233 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       70c513a4 Completed hicrep analysis except the creation of the graph
       1cf3bb2e hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       06f58f8a hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       adef366f hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       42888b6f hicseq completed adding basic features of the hicrep analysis.
       bf251fd6 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d156fd81 Added pairwise combination for samples
       0ff82520 started to edit the hicseq.py added new file for hicrep
       dd17ce73 Completed developing hicrep and quasar analysis
       db311862 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       467e6829 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1e11460c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       e9fac872 Completed hicrep analysis except the creation of the graph
       a2afec73 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       234855f2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       db9095ff hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b8f0c553 hicseq completed adding basic features of the hicrep analysis.
       a96d6c1d hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       66d63f73 Added pairwise combination for samples
       a69557d9 started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      55 commits

       b0eb60b8 update metylseq ini
       7ae859df crosscheck_fingerprint with EXIT_CODE_WHEN_MISMATCH 0
       f2171d23 add new line on genpipes_file.write()
       74130cd0 line explicite new lines at the end of fp.write()
       ea3c24b7 fix for pbs
       f8645cc5 ajust graham and cedar  ini
       129a4273 ajust graham ini
       c5a68759 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       b4cad788 cedar and graham ini maintenance
       8c6d99c4 tweak monitoring
       b81ace0d force -c on tumor pair ini files
       d5e8b078 remove unused import from hicseq
       8a46c262 log_report formating
       b94022d2 remove the verbose stuff that I am not monitoring anymore.
       6b44921e add missing args* in super call
       06ff02ba typo
       6d6b28bc update get_wrapper to v2.0.1
       b65dd6e1 remove decomissioned hpc ini files
       8fcf649e from output_file to genpipes_file
       ba89a19e --command_file options
       99292d87 remove -o option for output file
       96ee3056 used image v2.0.0.rc in for container
       3ecb2fb3 force TMPDIR to default to /tmp if not set
       1e3cc77a make sure output file is executable
       a6e17309 added option --output for command file
       66cf0f53 update form -n to -c
       e10f0903 removing useless shebang
       77d1a0dc extra space in samtools sort option
       10a9ffcf force pace between stringtie options
       047f332d make sure gemini annotation use the ini defined tmp dir
       fbf55522 force log report to python 2 :
       80d95085 error in path for chunk_genpipes csplit call
       b7eb42a0 force space before bash line continuation
       164e6f8a missing space in varscan script
       1d3bf917 update form -n to -c
       a40dce4f removing useless shebang
       350f8988 extra space in samtools sort option
       5e67fa4c force pace between stringtie options
       ba056af7 make sure gemini annotation use the ini defined tmp dir
       944fcc7d force log report to python 2 :
       a47fe860 error in path for chunk_genpipes csplit call
       9b19e7aa force space before bash line continuation
       f50600c7 missing space in varscan script
       a087236b update form -n to -c
       5a8d5121 excluding wget call from monitoring loop
       1dd1bdbf force loading mugqic python 2.7 on cuffmerge
       26183abe force loading mugqic python 2.7 on multiqc 1.7
       29735e86 more verbose went -d/--design is needed for contrast
       0a356017 revert on indel aligner
       c47384b0 cleanup ini for beluga
       c65298ea revert on indel aligner
       b18ca362 cleanup ini for beluga
       85147169 remove dbSNP for oxog ini block fix gatk 4 module name
       d5140585 fix picard_collect_oxog_metrics dbSNP vcf
       d630a603 remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       019c2572 Corrected hicseq.py file deletions
       0c780b0d Corrected hicseq.py file deletions
       9a30832a Corrected hicseq.py file deletions
       8f1db9ac Corrected hicseq.py file deletions
       621cddab Corrected hicseq.py file deletions
       a9dce2de Corrected hicseq.py file deletions
       94bb5a84 Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      14 commits

       cabe4d02 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       772bb686 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ef45d3a7 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       9bda9db8 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       93eb03e9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a7f50bfd hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       f79cacf7 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86ddffc2 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       7e5a282a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4de276e8 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3c940bb4 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       757f80fc hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       56c88759 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       86e066cb hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudumanoj <pubudu.nawarathna@mail.mcgill.ca>      82 commits

       982fdd83 addressed ED and PaulS's feedback on PR- step 3
       10920cb2 addressed ED and PaulS's feedback on PR- step 2
       498544e2 addressed ED and PaulS's feedback on PR- step 1
       d18311ca bugfix - unsuccesful
       804bc6c2 added step information to the md file and chipseq.py
       f7365e14 bugfix - unsuccesful
       f5decb53 chipseq.py - remove space in first line
       b85ff241 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       ed99e0a8 passed more parameters to R script
       a6696377 chipseq.py - changed output directory
       06a3f01b chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       0893090c adjusted alignment in functions
       b67ba87b skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       48e42ca1 added Rmarkdown rendering added more paramters
       775730d4 chipseq.py - changed output directory
       1bd300e1 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       101e2f61 chipseq.py - remove space in first line
       0aff721b differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b84e2375 fixed results files creating at custom output folder
       5830eaee added Rmarkdown rendering added more paramters
       cd4a5f10 chipseq.py - changed output directory
       8f70abf3 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       17b3232d bugfix - unsuccesful
       1d3b4639 chipseq.py - remove space in first line
       07ef003a differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       80f6cc0d fixed results files creating at custom output folder
       0d13c4f9 added Rmarkdown rendering added more paramters
       6f50be2c chipseq.py - changed output directory
       81695c6f chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       e3132b27 add diff_bind.py in bfx/
       47d33c2f added step information to the md file and chipseq.py
       53dc9235 bugfix - unsuccesful
       c2efc7d2 chipseq.py - remove space in first line
       b6416d60 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       bb677ee6 passed more parameters to R script
       15f23138 chipseq.py - changed output directory
       cc24da1a chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       48559699 adjusted alignment in functions
       196706e5 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0843c3a7 added Rmarkdown rendering added more paramters
       d25e5d0d chipseq.py - changed output directory
       8c17263d chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bf900c73 chipseq.py - remove space in first line
       62bd5a5e differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       b88a2e3b fixed results files creating at custom output folder
       84979112 added Rmarkdown rendering added more paramters
       6c780b25 chipseq.py - changed output directory
       ed75a1ee chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       bc9d0f2c added step information to the md file and chipseq.py
       43ab77e1 bugfix - unsuccesful
       994e168a chipseq.py - remove space in first line
       7413ea5b differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       241404f3 fixed results files creating at custom output folder
       a908dc37 added Rmarkdown rendering added more paramters
       26497f52 chipseq.py - changed output directory
       8e8d8e96 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       f71777e0 add diff_bind.py in bfx/
       1d1d959a bugfix - unsuccesful
       528bd387 chipseq.py - remove space in first line
       7805b91a differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       1a7e9576 fixed results files creating at custom output folder
       3f58f296 added Rmarkdown rendering added more paramters
       39732342 chipseq.py - changed output directory
       414155ce chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       9f12db9b add diff_bind.py in bfx/
       56028a81 chipseq.py - remove space in first line
       662aefaa chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       a8d2fcf0 differential analysis - fix parallel processing issue of rmarkdown by copying the R file with different name and delete it once the processing is finished. chipseq.py - skip analsysis for comparisons without having two samples for each group
       456d1143 fixed results files creating at custom output folder
       92e6fd77 adjusted alignment in functions
       4f665327 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       e09c78fb added Rmarkdown rendering added more paramters
       b65cf7e7 chipseq.py - changed output directory
       303efcda added paramters to ini files
       f927e563 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       abd9dd76 add diff_bind.py in bfx/
       96170dff passed more parameters to R script
       28ceb0ef chipseq.py - changed output directory
       696560f7 added paramters to ini files
       f46eb293 chipseq.py - corrected the input files and dependencises (added peak file instead of xls) -working version
       b9262488 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       641d393a add diff_bind.py in bfx/

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      14 commits

       eabb5bb7 refer previous commit
       5f8957f3 refer previous commit
       c225c6d1 refer previous commit
       1a091f54 refer previous commit
       540862d7 refer previous commit
       818f474a refer previous commit
       21ce258d refer previous commit
       94281fc6 refer previous commit
       c7774f2f refer previous commit
       d005007e refer previous commit
       a1cf1c0d refer previous commit
       4adc9ca6 refer previous commit
       3d25f609 refer previous commit
       25cddb0a refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      63 commits

       e34a5f38 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       94b62b91 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       304b5347 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       580338df fixed results files creating at custom output folder
       8fbd60e4 added paramters to ini files
       06b9ddf5 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       ec6762a2 passed more parameters to R script
       58aa030e added paramters to ini files
       99101ea0 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       3cc86223 rebasing
       1ed10791 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       61870419 adjusted alignment in functions
       c6c31f2a skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       038c730d passed more parameters to R script
       5502008f rebasing and resolving conflicts
       bb6343d5 added paramters to ini files
       ee7af7d2 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       f1498d8c rebasing
       745719ba chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       30ab7142 adjusted alignment in functions
       d86864ff skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       1699baef passed more parameters to R script
       9ba71cf0 rebasing and resolving conflicts
       ebc7a88b added paramters to ini files
       3055bd9e started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       bfe5a10d Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       63735774 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       044565e5 fixed results files creating at custom output folder
       ae1abe78 added paramters to ini files
       ad27f00b started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file
       68782121 passed more parameters to R script
       41ba25a3 added paramters to ini files
       b350cf2b started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2010b544 rebasing
       993dfa20 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ffd81543 adjusted alignment in functions
       80dcb4c7 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       186b7c7c passed more parameters to R script
       ac4e38d5 rebasing and resolving conflicts
       d412a8a7 added paramters to ini files
       428296b8 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       2d681ff4 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       36aa8693 rebasing
       20ea21c8 chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       8949c9eb adjusted alignment in functions
       f9c99d13 skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       ac7dbf09 passed more parameters to R script
       3c253135 rebasing and resolving conflicts
       e4d75d2c added paramters to ini files
       4c3b9401 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       d8c2155e Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis
       0d24cc41 rebasing
       35fff1fa chipseq.py - changed diff analysis step order, added new parameter for changing the method. bfx/differential_binding.py - changed R_script path chipseq.base.ini - added new paramter for method
       ca86b46c adjusted alignment in functions
       1076487b skipped when not samples are not enough fixed issue with -o outdir fixed an issue with design.py
       0bb499ba passed more parameters to R script
       716493df rebasing and resolving conflicts
       1e38ab5c added paramters to ini files
       4e22bde8 started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       78edf353 Merge branch 'chipseq_diff_analysis' of bitbucket.org:mugqic/genpipes into chipseq_diff_analysis resolved merge conflicts after pull
       0ff1cfd6 passed more parameters to R script
       3195dece started the development. changed the bfx.design.py -Added new function for chipseq design chipseq/chipseq.py - modified the file to call the new chip seq design function added code to create the input list for passing to the R file rebase
       9cceef8d corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      2 commits

       6199e349 Merged in chipseq_diff_analysis (pull request #240)
       ce590061 Merged dev into chipseq_diff_analysis

  Shaloo Shalini <shaloo.shalini@gmail.com>      9 commits

       b8f6ff31 Merged in ss_wf_83 (pull request #229)
       63c5ba1f Merged in ss_wf_82 (pull request #228)
       e38d0142 Merged in ss_wf_81 (pull request #227)
       0e0b2a23 Merged in ss_wf_80 (pull request #226)
       afdf55ab Merged in ss_wf_79 (pull request #225)
       55262398 Merged in ss_wf_78 (pull request #224)
       cd6f047a Merged in ss_wf_77 (pull request #223)
       02a08a42 Merged in ss_wflow_76 (pull request #222)
       eff450c0 Merged in ss_wflow_75 (pull request #221)

  shaloo <shalz@hotmail.com>      27 commits

       0c5149a0 Fixes #83 rnaseq denovo workflow diagram updated for v3.4.0
       116e1018 Fixes #82 rnaseq_light workflow updated for v3.4.0
       0e00fcdb Fixes #81 rnaseq workflows are now current to 3.4.0
       94b163c4 Fixes #80 methylseq pipeline workflow is now current with v3.4.0
       2811acfe Fixes #79 hicseq pipeline update for v3.4.0
       7679be51 Fixes #78 dnaseq workflow for -t mpileup updated v3.4.0
       37c5ca7f Fixes #77 dnaseq pipeline -t mugqic workflow update for genpipes v3.4.0
       531bbf5f Fixes #76 dnaseq_highcov pipeline workflow updated v3.4.0
       94315154 incorrect update should be for #76 cleaning up
       b5dd068c Fixes #76 dnaseq_higcov pipeline workflow updated in sync gpv3.4.0
       a9998d40 Fixes #75 updated amplicon sequence qiime and dada2 workflows for v3.4.0
       3bd1caa3 Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       8cbae5ef Refs #67 dnaseq -t light option color feedback addressed
       1befe532 Fixes #67 Refs #54 dnaseq -t light workflow schema added
       da225756 Refs #65 chipseq workflow color updated, report links added as per feedback
       0f63c1ef Fixes #65 Refs #54 chipseq pipeline updated
       c3298477 Refs #68 color added as per feedback
       78895de9 Fixes #68 Refs #54 nanopore pipeline workflow schema added
       1103637f Refs #66 cleanup after color update
       742de766 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       18010fb4 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       dcbf14d8 Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       97aa943a Refs #53 covseq workflow arrow step 6 -> 9 removed
       2d4a60d3 Refs #53 Paul's review inputs incorporated
       08e8068b Fixes #53 covseq.py workflow diagram added
       8911c4cb Refs #53 added covseq pipeline schema workflow draft under review by Paul
       fb8f6f6e Fixes #51 update covseq pipeline readme to v3.3.0

3.4.0        Thu Apr 29 19:24:01 2021 +0000        784 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      1 commits

       ff50d950 GenPipes - RNASeq : corrected genome_index_folder refencing in star_align

  Édouard Henrion <henrione@beluga2.int.ets1.calculquebec.ca>      10 commits

       2b72afb3 GenPipes - Release : removing Tumor pair prior release 3.4
       9de7a779 Merge remote-tracking branch 'origin/dev' into release_3.4
       ed27f331 GenPipes - Resources : adding software installation scripts
       8c9f6139 GenPipes - Resources : adding software installation scipts
       c2d24183 GenPipes - Resources : adding reference genome installation scripts
       3bc6a801 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       67cc4994 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       823f6264 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       e42f490a GenPipes - Resources : adding Sus_scrofa genome and ivar & kent installation scripts
       d76314e0 GenPipes - Resources : updates of solftware and genome installation scripts

  ehenrion <edouard.henrion@mcgill.ca>      79 commits

       ccc89c11 Merged in release_3.4 (pull request #219)
       9667e558 Resolve Issue #31
       41c771b5 Merged dev into eh_cit_correction
       5fd600c3 GenPipes - BFX : corrected ini section for star_index --outTmpDir
       982a3050 GenPipes - BFX : update STAR calls to use --outTmpDir
       7a6549de Merged dev into chipseq_design_change
       0ca16958 Merged dev into chipseq_design_change
       23f6e6d4 Merged in eh_cit_correction (pull request #207)
       daf9eb8a GenPipes - HiCSeq : corrected typo in CHICAGO makeDesignFiles call
       2e1b75e9 GenPipes - HiCSeq : updated base.ini with explicit loading of mugqic/python/2.7.14 in chicago create_design_files step
       0db8b444 GenPipes - HiCSeq : corrected CHICAGO makeDesigFiles call with explicit load of python2 module
       9d5c8af8 Merged dev into chipseq_design_change
       09b50e33 Merged dev into eh_cit_correction
       af89867e GenPipes - RNASeq : corrected genome_index_folder referencing in  star_align
       024b2517 Merged eh_RNAseq_star_correct into dev
       81e394bf GenPipes - RNASeq : corrected genome path in star_align
       0dbd8958 Merged eh_fix_callhome_fail_exit_code into dev
       a8bfffda GenPipes - Call Home : fixed wget command in common.py to always exit 0 in order to avoid crash of GenPipes execution - Issue #63
       b652f60d Merged eh_RNAseq_star_correct into dev
       14f5ca51 VERSION edited online with Bitbucket
       4a2a45df GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       ef26b52d GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       01d00902 GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       73f22e32 GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       6ae11395 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       a3781df2 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       fde1dbc9 GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       533b3fa8 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       43337528 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       9cda20cd GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       ceb95945 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       ed1a35b4 GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       6c61ca59 GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       d250236c GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       a5dc00cd GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       15511570 GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       404c9312 GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       ba277e42 GernPipes - DNASeq : corrected merge_filter_bcf outputs
       1028ccc3 VERSION edited online with Bitbucket
       82e6690e Merged in ehenrion/version-edited-online-with-bitbucket-1617908341194 (pull request #205)
       e2249aad VERSION edited online with Bitbucket
       0f75fb63 Merged eh_samtools_cram_output_ini_fix into dev
       74d7c56c GenPipes - Config : fixed samtools_cram_output in rnaseq.graham.ini
       a45ed725 GenPipes - Config : fixed samtools_cram_output in rnaseq.cedar.ini
       0da38f2c GenPipes - Config : fixed samtools_cram_output in methylseq.graham.ini
       5ded7767 GenPipes - Config : fixed samtools_cram_output in methylseq.cedar.ini
       b6437252 GenPipes - Config : fixed samtools_cram_output in methylseq.beluga.ini
       33a4b963 GenPipes - Config : fixed samtools_cram_output in hicseq.graham.ini
       053b208b GenPipes - Config : fixed samtools_cram_output in hicseq.cedar.ini
       da787e92 GenPipes - Config : fixed samtools_cram_output in hicseq.beluga.ini
       713bd648 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.graham.ini
       716110de GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.cedar.ini
       68c12941 GenPipes - Config : fixed samtools_cram_output in dnaseq_high_coverage.beluga.ini
       26ca812c GenPipes - Config : fixed samtools_cram_output in chipseq.graham.ini
       065df8d3 GenPipes - Config : fixed samtools_cram_output in chipseq.cedar.ini
       98f53913 GenPipes - Config : fixed samtools_cram_output in dnaseq.graham.ini
       2d7892d0 GenPipes - Config : fixed samtools_cram_output in dnaseq.cedar.ini
       cbe399cd GenPipes - Config : fixed samtools_cram_output in dnaseq.beluga.ini
       b74a4763 GenPipes - Config : fixed samtools_cram_output in chipseq.beluga.ini
       938318a8 Merged eh_cit_correction into dev
       4d076784 GenPipes - DNASeq : corrected iteration on samples in cnvkit_sv_annotation
       c095c6ca GernPipes - DNASeq : corrected merge_filter_bcf outputs
       c5d36c96 GenPipes - RNASeq : corrected samtools_cram_output in beluga.ini - Issue #62
       0d9ab3ea GenPipes - DNASeq SV : fixing delly call in dnaseq.py
       138081b8 GenPipes - DNASeq SV : fixing delly input error
       79f6e949 Merged in ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004 (pull request #202)
       ce2014cc Merged eh_RNAseq_star_correct into ehenrion/genpipes-rnaseq-updated-starpy-to-test-1616421770004
       f1a6ffd1 GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       f53266bd GenPipes - RNAseq : updated star.py to test if the version of STAR is included in the genome_index_folder. If not, use the old format of path (i.e. without the version of STAR)
       8a51f4b3 Merged in ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543 (pull request #200)
       e0cf3d83 Merged eh_RNAseq_star_correct into ehenrion/genpipes-ranseq-updated-starpy-to-add--1616420528543
       66c603f1 GenPipes - RANSeq : updated star.py to add the version of STAR in the genome index folder path
       77ced18c GenPipes - RNASeq : star.align updated base.ini with the version of star in the path of index
       176669b4 Merged in eh_fix_delly_issue52 (pull request #199)
       2fde5cb4 GenPipes - BFX : corrected delly.py 'call' input handling [Issue 52](https://bitbucket.org/mugqic/genpipes/issues/52/version-330-dnaseq-t-sv)
       277904f5 dnaseq.py edited online with Bitbucket : corrected protocol assgignation
       fdbef3ad GenPipes - Bug fix : corrected dnaseq.cedar.ini
       5824422c GenPipes - Bug fix : correcting indentation in illumina_run_processing.py
       9cfee2e7 dnaseq.py edited online with Bitbucket : corrected protocol assgignation

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      3 commits

       a0a8e94f rnaseq.base.ini updated to a newer version of STAR
       22f3dfad Merged in rnaseq_star_update (pull request #195)
       06a8d473 rnaseq.base.ini updated to a newer version of STAR

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      219 commits

       d0835ca2 Merged in chipseq_design_change (pull request #210)
       f7655cbd Merged dev into chipseq_design_change
       ce5c0079 Fixing atacseq
       2e0ee42a Merged in chipseq_design_change (pull request #206)
       ef03626c Update READMEs
       376858ba Merged dev into chipseq_design_change
       182612e7 Updating inis and fixing homer output folder location
       f27a211b Merged dev into chipseq_design_change
       0ce59917 Debug
       4c0b11f1 Debug
       3f6fcc6d Debug
       9cbe87c7 Debug
       586708f0 Merge branch 'dev' into chipseq_design_change
       5cb6a403 Debug
       1980df20 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       318284c3 started to edit the hicseq.py added new file for hicrep
       03e1ebc6 Fixing minor bug and typo
       48290070 Fixing multiqc dependencies
       9acc4eb7 Debug
       9b72ff08 Debug
       03c9cb67 Debug
       a159e49d Debug
       2c9cf566 Debug
       76664e89 Debug
       7e90a866 Debug
       4a4ae351 Debug
       5718604f Debug
       4c01435c Debug
       bb6ee3f1 Debug
       aaf4809b Merged dev into chipseq_design_change
       ae98ed2e Debug
       d105f719 Debug
       9c42079d Changing name of header for a report table
       786ac814 Debug
       d98c7ba8 Typo
       8b1b565f Addinx x permission to job2json.py
       2de4ae65 Fixing report naming too long
       d53e173d Fixing report naming too long
       a3d412e0 Merged dev into chipseq_design_change
       40bcb0f2 Fixing report naming too long
       a399b3d6 Fixing mugqicValidator.py for new readset format (chipseq only)
       2249d6a0 Switching to latest mugqic_tools release
       3b77185d Changing Ressources
       34524b8c Changing Ressources
       af360e36 Changing Ressources
       4738bda0 Changing module versions
       6d5e26a1 Debug
       860b44d8 Fix
       ef12823e Merged dev into chipseq_design_change
       6512f2f1 Changing qualimap bamqc section naming to a parameter in bfx
       9b4d61ce Increasing picard collect multiple metrics RAM
       b264e7dc Debug
       030fb31a Debug
       8de52dd6 Fix ressources
       dacaec86 Fixing annotation peaks
       12e5d866 Fixing annotation peaks
       9d9a50f5 Fixing annotation peaks
       351f5c31 Adding annotation peaks merge for Narrow peaks only
       0a2e098c Iterating only on Narrow peaks for annotation
       fd3bd746 Iterating only on Narrow peaks for annotation
       58da5b0b Debug report
       687969cd Debug
       955fe4bb Debug
       d7560d90 Debug
       e9e09de8 Debug
       4b8ab24c Renaming IHEC metrics section
       8f9ee567 Debug
       b73ce498 Debug
       43bebdfb Debug Report
       cbf2b319 Debug Report
       d8943245 Debug
       5dc81e88 Debug
       98d61cfc Fixing report
       3b179a45 Debug
       74cd55db Debug
       d3e03b28 Debug
       6545cbd6 Debug
       9f90aacf Debug
       b94c8a1e Debug
       03279f65 Debug
       830e4565 Debug
       79cb6542 Debug
       19a0024e Debug
       cf9d2527 Debug
       505b1cae Debug
       bbf28520 Changing Output metrics
       e654a974 Debug
       43105ff5 Debug
       2092de56 Debug
       3c6c7b45 Debug
       0ab7349a Debug
       2ee5418c Debug
       a30cf4a7 Debug
       ec6d257a Debug
       610181d7 Debug
       a85c9742 Debug
       bfd2e2fe Debug
       f21e83ea Debug
       3fb54184 Debug
       d77307e1 Debug
       ee397ef8 Changing macs2 atacseq param & add macs2 param
       8c4274b4 Debug
       a28d2758 Debug
       a64ce532 Debug
       a689539e Debug
       9587ca81 Debug
       e370f820 Debug
       0d3741d7 Debug
       db78c589 Debug
       f5d866d1 Debug
       79656256 Debug
       fdf3d9bf Debug
       7847ac05 Debug
       47da9ef2 Debug
       6a8b1ba9 Debug
       dddcf884 Debug
       161353f8 Debug
       73790f90 Debug
       188496e0 Debug
       361b1754 Debug
       f5940eda Debug
       e502b6fb Changing R version for mugqic_tools R scripts
       45ef768f homer_annotate_peaks
       6a70faf5 qc_metrics report
       49e8bfe8 Fix ihec_metrics outputs
       cc35a5dc Fix ihec_metrics outputs
       29c13621 Fix ihec_metrics outputs
       7c680718 Fixing MultiQC report
       12212bd6 Fixing MultiQC report
       c46e42dc Fixing MultiQC report
       73ced5d5 Fixing MultiQC report
       08c1f435 Fixing MultiQC report
       b403bfbd Fixing MultiQC report
       bdc3e467 Improving MultiQC report
       cde658a1 Debug
       528f0476 Debug
       cccc2e92 Debug
       a545b4ab Debug
       2b79de34 Debug
       acd95d27 Debug
       c0d1a393 Debug
       92d584e2 Debug
       0c54b4ed Debug
       c919c900 Debug
       e9c59559 Debug
       64657e3b Debug
       69857a49 Debug
       2683f885 Debug
       e8a45bf5 Debug
       2570791b Debug
       6bf3d107 Debug
       524a5c0c Debug
       985fa272 Debug
       3b6ce4a5 Debug
       ee40fde5 Debug
       6eedc38f Debug
       34383e75 Major changes IHEC metrics test
       2385c51a Major changes IHEC metrics test
       72f30ff7 Major changes test
       9dfe33e8 Macs2 report changes
       705e66a7 Major change test
       e23426a5 Major change test
       7d787c14 Major change test
       58b992c5 debug
       b9721511 debug
       9f4c020b debug
       4a14c426 debug
       e69cf403 debug
       20a8959a debug
       d800207a debug
       d0cebc14 debug
       38616a43 debug
       3f0fea1a debug
       fcca2856 debug
       61ff7d36 debug
       d147d39a debug
       10e30302 debug
       70c0c5d5 debug
       acc9d2c8 debug
       68ea5f43 debug
       a19d64f9 debug
       76f51557 debug
       8a922595 debug
       9e6c5349 debug
       3b9bf7a4 debug
       1e35666d debug
       b2f165b7 debug
       bf557912 debug
       4677380f debug
       5cc5350d debug
       f0d9f385 debug
       b306bbd4 debug
       1d7d8716 debug
       89123424 Fix test
       166b57b5 Major readset change test
       8db202b7 Fix
       abfcb3d8 Fix
       ae7964ac Fix
       c9bc7506 Filtering after merging and changing naming to fit with other pipelines
       aae46094 Increasing default resources and adding options for markdup
       9f540c4d Fixing beluga ini
       dd795e63 Switching to sambamba markdup for IHEC
       27e1a087 Fix
       ecc127db Fix
       e01ccf2a Fix
       b299062a Fix
       f8cf1225 Fix
       264059ca Fix
       df2604b8 Options becomes optional
       1748eb6c Fix
       c09bf2f2 Fixing typo
       adf2e7fa Fix
       6c4cbd38 Fixing sambamba merge
       9cf59889 Typo
       fc6f4cb8 Adding mkdir
       f060825c Fixing temp typo to tmp
       aca35ab9 Fixing job
       95e56619 Fixing bash import
       f5508117 Fixing minor bug and typo

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      385 commits

       03dbf551 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2bd98087 hicseq completed adding basic features of the hicrep analysis.
       7c19ba7a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3b8a191 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       32ffa563 hicseq [hicrep.py] - Corrected typo
       5321e34a hicseq [hicrep.py] - corrected R_TOOLS path
       7d56a211 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fe8103ab started to edit the hicseq.py added new file for hicrep
       9f24b59a hicseq [hicrep.py] - Corrected typo
       78d9e596 hicseq [hicrep.py] - corrected R_TOOLS path
       ac8395d5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       56f9aa36 hicseq [hicrep.py] - Corrected typo
       f2f843b5 hicseq [hicrep.py] - corrected R_TOOLS path
       2e3a2e80 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b242698f hicseq [hicseq.py] - corrected file after rebase
       b9c2899f hicseq [hicrep.py] - Corrected typo
       56fb9ebb hicseq [hicrep.py] - corrected R_TOOLS path
       bd160f60 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3726aea9 Completed developing hicrep and quasar analysis
       ef6f77b4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e3426 Completed hicrep analysis except the creation of the graph
       e28dd55a started to edit the hicseq.py added new file for hicrep
       5d82058d hicseq [hicseq.py] - corrected file after rebase
       9f494ca9 hicseq [hicrep.py] - Corrected typo
       a5f7a903 hicseq [hicrep.py] - corrected R_TOOLS path
       e098a617 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3defd1e3 hicseq [hicrep.py] - Corrected typo
       3acb37df hicseq [hicrep.py] - corrected R_TOOLS path
       f18061a5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ec44a6b hicseq [hicrep.py] - Corrected typo
       ac36fca9 hicseq [hicrep.py] - corrected R_TOOLS path
       911d6693 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8b851f36 hicseq [hicrep.py] - Corrected typo
       aee4c68a hicseq [hicrep.py] - corrected R_TOOLS path
       a4ac85ce [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       766f8db6 hicseq [hicrep.py] - Corrected typo
       8f9a52d3 hicseq [hicrep.py] - corrected R_TOOLS path
       e2a85e87 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eca4b356 hicseq [hicrep.py] - Corrected typo
       cae9d46a hicseq [hicrep.py] - corrected R_TOOLS path
       c742e6c9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9bd1c6ee hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       ea5095ff hicseq [hicrep.py] - Corrected typo
       0d7d6c08 hicseq [hicrep.py] - corrected R_TOOLS path
       99e1052f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5deb5597 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       31d1b1bb hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       dc822e8d hicseq [hicseq.py, readme.md] - modified readmefile
       b828c729 hicseq [hicseq.py] - corrected file after rebase
       142dfc30 hicseq [hicrep.py] - Corrected typo
       1fd82388 hicseq [hicrep.py] - corrected R_TOOLS path
       b0d0045a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ffaba4c6 hicseq [hicrep.py] - Corrected typo
       16866fb8 hicseq [hicrep.py] - corrected R_TOOLS path
       5c52f5ab [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1ca85421 hicseq [hicrep.py] - Corrected typo
       07ce5ba9 hicseq [hicrep.py] - corrected R_TOOLS path
       9886c9e8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       83e13141 hicseq [hicrep.py] - Corrected typo
       4ee15eb6 hicseq [hicrep.py] - corrected R_TOOLS path
       8ae9494c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7de77365 hicseq [hicrep.py] - Corrected typo
       cc6564d5 hicseq [hicrep.py] - corrected R_TOOLS path
       ebcc1a11 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1eafa49 hicseq [hicrep.py] - Corrected typo
       380dfd47 hicseq [hicrep.py] - corrected R_TOOLS path
       4d9f8531 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c2f186d7 hicseq [hicrep.py] - Corrected typo
       654198f7 hicseq [hicrep.py] - corrected R_TOOLS path
       c96c4757 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bf2449c9 hicseq [hicrep.py] - Corrected typo
       342b6f58 hicseq [hicrep.py] - corrected R_TOOLS path
       5831718a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6f6e3c5d hicseq [hicseq.py] - corrected file after rebase
       bac4d8d1 hicseq [hicrep.py] - Corrected typo
       81688407 hicseq [hicrep.py] - corrected R_TOOLS path
       e4fdcd0e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7250a4f2 hicseq [hicrep.py] - Corrected typo
       3942f566 hicseq [hicrep.py] - corrected R_TOOLS path
       10b31a81 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b6805e64 hicseq [hicrep.py] - Corrected typo
       03ed962c hicseq [hicrep.py] - corrected R_TOOLS path
       fa3d60b0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       9c4734f5 hicseq [hicrep.py] - Corrected typo
       dccfad08 hicseq [hicrep.py] - corrected R_TOOLS path
       907e9325 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b3e8d9a9 hicseq [hicrep.py] - Corrected typo
       cbc0e992 hicseq [hicrep.py] - corrected R_TOOLS path
       a15ba0c9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3a6a1358 hicseq [hicrep.py] - Corrected typo
       7d17401f hicseq [hicrep.py] - corrected R_TOOLS path
       3c495bf1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e8333655 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       5dd2d25e hicseq [hicrep.py] - Corrected typo
       62b39434 hicseq [hicrep.py] - corrected R_TOOLS path
       9a8f22a7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d95fe9f2 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff71f6fd hicseq completed adding basic features of the hicrep analysis.
       35c32698 Added pairwise combination for samples
       85571a9d started to edit the hicseq.py added new file for hicrep
       72be0065 hicseq pipeline [changed the step order]
       0800b925 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       3363d075 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       e5b3070e hicseq [hicseq.py, readme.md] - modified readmefile
       d4dbe19a hicseq [hicseq.py] - corrected file after rebase
       adf7ceb8 hicseq [hicrep.py] - Corrected typo
       9f9e0712 hicseq [hicrep.py] - corrected R_TOOLS path
       84815043 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d286d5dd hicseq [hicrep.py] - Corrected typo
       8f8c6128 hicseq [hicrep.py] - corrected R_TOOLS path
       85a1cf96 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8eca7d4a hicseq [hicrep.py] - Corrected typo
       9229a936 hicseq [hicrep.py] - corrected R_TOOLS path
       5a60eff2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef77769e hicseq [hicrep.py] - Corrected typo
       467b8e94 hicseq [hicrep.py] - corrected R_TOOLS path
       c0559eb9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c8ff0042 hicseq [hicrep.py] - Corrected typo
       6e1367f1 hicseq [hicrep.py] - corrected R_TOOLS path
       2266dd99 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       999455da hicseq [hicrep.py] - Corrected typo
       240eaf33 hicseq [hicrep.py] - corrected R_TOOLS path
       f2ab622d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d0755352 hicseq [hicrep.py] - Corrected typo
       28e544c5 hicseq [hicrep.py] - corrected R_TOOLS path
       5312b560 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4161b94f hicseq [hicrep.py] - Corrected typo
       7ddd445a hicseq [hicrep.py] - corrected R_TOOLS path
       5a24278c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1dd26f40 hicseq [hicseq.py] - corrected file after rebase
       3d526b79 hicseq [hicrep.py] - Corrected typo
       10196c8f hicseq [hicrep.py] - corrected R_TOOLS path
       66cf3f7d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2b72cc9a hicseq [hicrep.py] - Corrected typo
       0a1ae4bd hicseq [hicrep.py] - corrected R_TOOLS path
       ec8ec063 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       52ceb831 hicseq [hicrep.py] - Corrected typo
       5691c948 hicseq [hicrep.py] - corrected R_TOOLS path
       33907621 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       071ed297 hicseq [hicrep.py] - Corrected typo
       6f23e5fd hicseq [hicrep.py] - corrected R_TOOLS path
       29dd0599 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bdd9b878 hicseq [hicrep.py] - Corrected typo
       ee032b05 hicseq [hicrep.py] - corrected R_TOOLS path
       0325e110 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f9cb9baf hicseq [hicrep.py] - Corrected typo
       e2f2d6bc hicseq [hicrep.py] - corrected R_TOOLS path
       5c92379e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0f4de0bb hicseq [base.ini] - updated mugqic tools version to 2.3.1
       f1361d63 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       fe828589 hicseq [quasar_qc - corrected module loading in matrix restructuring
       589286e2 hicseq [hicrep.py] - Corrected typo
       9dcbc833 hicseq [hicrep.py] - corrected R_TOOLS path
       85b4bbd5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4cc5884e hicseq [hicrep.py] - Corrected typo
       26d8fde4 hicseq [hicrep.py] - corrected R_TOOLS path
       689a76b6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5529f4d1 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       2d0af167 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c60cd157 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       a5cfa02c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c34a220a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       89a55c1d Completed developing hicrep and quasar analysis
       0b74a0aa [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       54c0ea89 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       45c545a2 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       3ea32860 Completed hicrep analysis except the creation of the graph
       f6f2a519 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       034732ad hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       777005f7 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       36ab2e9b hicseq completed adding basic features of the hicrep analysis.
       f9a4ef6b hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7fff19ee Added pairwise combination for samples
       5ae1358e started to edit the hicseq.py added new file for hicrep
       875a2afd hicseq [hicseq.py, readme.md] - modified readmefile
       afb029fe hicseq [hicseq.py] - corrected file after rebase
       ae2799dd hicseq [hicrep.py] - Corrected typo
       22c4aeaa hicseq [hicrep.py] - corrected R_TOOLS path
       a8b525d6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5a84543e hicseq [hicrep.py] - Corrected typo
       29ef5779 hicseq [hicrep.py] - corrected R_TOOLS path
       8e99a729 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0af5ff6 hicseq [hicrep.py] - Corrected typo
       2c514bf0 hicseq [hicrep.py] - corrected R_TOOLS path
       9c984a06 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65e4495a hicseq [hicrep.py] - Corrected typo
       d6b6bdab hicseq [hicrep.py] - corrected R_TOOLS path
       15c7831c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a98b3265 hicseq [hicrep.py] - Corrected typo
       bfc0ee4f hicseq [hicrep.py] - corrected R_TOOLS path
       2dd0caa5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6070378d hicseq [hicrep.py] - Corrected typo
       b7f35aee hicseq [hicrep.py] - corrected R_TOOLS path
       c739410f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       473304ea hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b536a374 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2e9b2661 hicseq [quasar_qc - corrected module loading in matrix restructuring
       e40bf89f hicseq [hicrep.py] - Corrected typo
       ec4b642d hicseq [hicrep.py] - corrected R_TOOLS path
       6ebce209 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8a7e4f11 hicseq [hicrep.py] - Corrected typo
       4ec40e07 hicseq [hicrep.py] - corrected R_TOOLS path
       d7e39f9a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1cc96d3b hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5e52fc32 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       59c65ac5 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fdfbe14 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       e5eb0973 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e95cd4db Completed developing hicrep and quasar analysis
       c5c79082 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       08ff8a78 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d6abcd4d created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a6156a17 Completed hicrep analysis except the creation of the graph
       4eea55cc hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       43501ba2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       25c34b43 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a164d8bf hicseq completed adding basic features of the hicrep analysis.
       e8aa04c3 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       491fa1ce Added pairwise combination for samples
       b42241ad started to edit the hicseq.py added new file for hicrep
       20627d69 hicseq [hicseq.py] - corrected file after rebase
       884d4e4a hicseq [hicrep.py] - Corrected typo
       622d5f3a hicseq [hicrep.py] - corrected R_TOOLS path
       c4a9c60c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       fb0ec53b hicseq [hicrep.py] - Corrected typo
       7082fd6a hicseq [hicrep.py] - corrected R_TOOLS path
       f90b5be4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       02600a93 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       456e6495 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d96bfb1c hicseq [quasar_qc - corrected module loading in matrix restructuring
       86c12350 hicseq [hicrep.py] - Corrected typo
       79fca771 hicseq [hicrep.py] - corrected R_TOOLS path
       c28c3f88 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef796f98 hicseq [hicrep.py] - Corrected typo
       30375570 hicseq [hicrep.py] - corrected R_TOOLS path
       b5b40fe1 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e538fb0c hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       28dee2f7 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d7e7c279 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       82aa89a0 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       00320364 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db025271 Completed developing hicrep and quasar analysis
       edf03c15 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7b1b7605 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2520bb90 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       314a5054 Completed hicrep analysis except the creation of the graph
       71f75b6d hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       723cecf5 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       38ec796a hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       fb02b115 hicseq completed adding basic features of the hicrep analysis.
       56648793 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       696a52bd Added pairwise combination for samples
       886ea4ed started to edit the hicseq.py added new file for hicrep
       10d4d24c hicseq [base.ini] - updated mugqic tools version to 2.3.1
       9d7b8658 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       2c3b54d6 hicseq [quasar_qc - corrected module loading in matrix restructuring
       dbdf0ca6 hicseq [hicrep.py] - Corrected typo
       e70f0ee0 hicseq [hicrep.py] - corrected R_TOOLS path
       92a5667a [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c05b038f hicseq [hicrep.py] - Corrected typo
       1d57da7e hicseq [hicrep.py] - corrected R_TOOLS path
       09b27e05 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       565b4f46 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       f9f52826 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       79da62c2 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       4518e456 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       11e583d0 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f7d773d0 Completed developing hicrep and quasar analysis
       b7f60de8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b3af6960 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       0034cea7 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2a5492d4 Completed hicrep analysis except the creation of the graph
       d6e8b9b5 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       f1f7d98c hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ba1d617d hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       76e10ba8 hicseq completed adding basic features of the hicrep analysis.
       72a3137a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c2af6c67 Added pairwise combination for samples
       6f0965b2 started to edit the hicseq.py added new file for hicrep
       e9f3fe40 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       2c206ace hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3e2a53ed hicseq [quasar_qc - corrected module loading in matrix restructuring
       f880a497 hicseq [hicrep.py] - Corrected typo
       809f95e6 hicseq [hicrep.py] - corrected R_TOOLS path
       aaddec59 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       251c0c6f hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       68889f95 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       68003e05 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2fef2208 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8782a74b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e2caf449 Completed developing hicrep and quasar analysis
       30008d63 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7a640583 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       02e39e7c created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ea8610d2 Completed hicrep analysis except the creation of the graph
       48173c51 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       419d60e0 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       b26f6236 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2f2eba58 hicseq completed adding basic features of the hicrep analysis.
       92b291a5 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c3dcb3ec Added pairwise combination for samples
       67488a04 started to edit the hicseq.py added new file for hicrep
       a22ca567 hicseq [hicrep.py] - Corrected typo
       e28aedee hicseq [quasar_qc.py] - Corrected deletion by mistake
       5a1f5025 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       5b6fff8d hicseq [hicrep.py] - corrected R_TOOLS path
       b02d2280 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4ed3b593 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       c361bc62 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       56223893 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       e4758a19 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       c5641494 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       760bfd64 Completed developing hicrep and quasar analysis
       275a2243 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       627a5016 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       33ffee0b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       df084a84 Completed hicrep analysis except the creation of the graph
       f3bebd7d hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       b7fe3828 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       805c4a96 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       dbfde764 hicseq completed adding basic features of the hicrep analysis.
       10e2dfa7 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       7c2e89b8 Added pairwise combination for samples
       07099a0a started to edit the hicseq.py added new file for hicrep
       300ba1bb hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       688a4956 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9aa6c614 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       0f393c68 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       b64480f8 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       353a55d1 Completed developing hicrep and quasar analysis
       7aa1f48b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       008904a0 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b27275b6 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       66c49d4e Completed hicrep analysis except the creation of the graph
       a09f0387 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5c370b55 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       de693aec hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       a504b171 hicseq completed adding basic features of the hicrep analysis.
       fe551f3c hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       52208162 Added pairwise combination for samples
       17e1b5c0 started to edit the hicseq.py added new file for hicrep
       8d592db3 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c37d365f [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8c152446 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       7f55aebe [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       f371fb96 Completed developing hicrep and quasar analysis
       14b332d0 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       66ded908 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c8cabe6b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       f066550e Completed hicrep analysis except the creation of the graph
       0a733f01 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       445a4684 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       ebcf685e hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8e3602c4 hicseq completed adding basic features of the hicrep analysis.
       21b08c35 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       59771c35 Added pairwise combination for samples
       98194825 started to edit the hicseq.py added new file for hicrep
       8fc473c6 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       673bad3e [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       979085be [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       5392e006 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       96b5c41e Completed developing hicrep and quasar analysis
       8d38fc18 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9b2e517d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dade053a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       14247d25 Completed hicrep analysis except the creation of the graph
       79e675ea hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       59125936 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       c19d1be7 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       cc2f34c1 hicseq completed adding basic features of the hicrep analysis.
       e7e729af hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8495df61 Added pairwise combination for samples
       d789632a started to edit the hicseq.py added new file for hicrep
       d8a1537b Completed developing hicrep and quasar analysis
       c13cb74d [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       316db2c7 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a382fb8 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       68fc66e6 Completed hicrep analysis except the creation of the graph
       f2cb6a0b hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       702fadb2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4fefe1e9 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       294e3aa3 hicseq completed adding basic features of the hicrep analysis.
       5edc80e2 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       fc95853d Added pairwise combination for samples
       15925106 started to edit the hicseq.py added new file for hicrep

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      10 commits

       4b08c93c force loading mugqic python 2.7 on cuffmerge
       ac7cd868 force loading mugqic python 2.7 on multiqc 1.7
       b0875ae9 more verbose went -d/--design is needed for contrast
       0dc71a9f revert on indel aligner
       37af5701 cleanup ini for beluga
       3c0fad0e revert on indel aligner
       54635678 cleanup ini for beluga
       2080dfda remove dbSNP for oxog ini block fix gatk 4 module name
       3a1ab223 fix picard_collect_oxog_metrics dbSNP vcf
       1c5e2ed4 remove HOME_DEV paths

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      7 commits

       7bbe4faf Corrected hicseq.py file deletions
       0ca8a123 Corrected hicseq.py file deletions
       e2e57af0 Corrected hicseq.py file deletions
       9d90aeff Corrected hicseq.py file deletions
       312d3dba Corrected hicseq.py file deletions
       a95d8b8a Corrected hicseq.py file deletions
       52c9240c Corrected hicseq.py file deletions

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      11 commits

       39db7bb7 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       49432742 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       452015d9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       61a86ece hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       e7985498 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       185c81b6 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21929bf1 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       fe03ec69 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cdca1f90 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       a146b251 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4025c716 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      11 commits

       c646e675 refer previous commit
       30b0838e refer previous commit
       0e4a6cb1 refer previous commit
       138fafc1 refer previous commit
       0fd5742d refer previous commit
       962cc3e5 refer previous commit
       bb6b1b7d refer previous commit
       a1f0ba28 refer previous commit
       2b79a1bd refer previous commit
       ae7e5c82 refer previous commit
       3b6c393a refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       30b2d60a corrected md file

  Shaloo Shalini <shaloo.shalini@gmail.com>      17 commits

       f0d10456 Merged in ss_issu_69_rnaseq_strintie_wf (pull request #214)
       6935d81c Merged in ss_issu_67_wf_dnaseq_light (pull request #216)
       2b8fe8c2 Merged in ss_issu_66_dnaseq_sv_schema (pull request #215)
       17c8fcc4 Merged dev into ss_issu_69_rnaseq_strintie_wf
       0b1fe1eb Merged dev into ss_issu_66_dnaseq_sv_schema
       0f9def7c Merged dev into ss_issu_67_wf_dnaseq_light
       37504fcd Merged in ss_issu_67_wf_dnaseq_light (pull request #212)
       b217c763 Merged in ss_issu_65_chipseq_workflow (pull request #211)
       890dc407 Merged in ss_issu_68_nano_wf (pull request #213)
       1c7a288c Merged dev into ss_issu_67_wf_dnaseq_light
       0a6fd274 Merged dev into ss_issu_69_rnaseq_strintie_wf
       c7272d59 Merged dev into ss_issu_68_nano_wf
       0a1ffe39 Merged dev into ss_issu_65_chipseq_workflow
       3f9add6e Merged in ss_issu_66_dnaseq_sv_schema (pull request #209)
       2bddc593 Merged in ss_issu_60_wf_ampllconseq (pull request #203)
       71461066 Merged in ss_covseq_wflow (pull request #197)
       192f7a30 Merged in covseq_readme (pull request #196)

  shaloo <shalz@hotmail.com>      30 commits

       14cb58cf Refs #69 Hector and Ed's feedback incorporated
       ea4ea09a Refs #66 Feedback from Ed and Rob wrt step dependency has been addressed
       de861a5d Refs #66 cleanup after color update
       1b686195 Refs #67 cleanup after color update
       26733d9d Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       75101db1 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       47c87054 Refs #66 color updated as suggested for dnaseq light same for dnaseq sv option
       dc52d268 Refs #67 dnaseq -t light option color feedback addressed
       0c95d465 Refs #68 color updated as per feedback
       796cf88c Refs #68 color added as per feedback
       1ad517e6 Refs #65 chipseq workflow color updated, report links added as per feedback
       179955f4 Fixes #69 Refs #54 rnaseq -t stringtie workflow schema added
       c1a88613 Merge branch 'ss_issu_65_chipseq_workflow' of bitbucket.org:mugqic/genpipes into ss_issu_65_chipseq_workflow
       ce49d390 Fixes #65 Refs #54 chipseq pipeline updated
       8beba66d Fixes #68 Refs #54 nanopore pipeline workflow schema added
       cc91732d Fixes #67 Refs #54 dnaseq -t light workflow schema added
       9f4b8ce6 Update chipseq pipeline schema for -t chipseq and -t atacseq options to reflect latest dev branch pipeline code
       6f884568 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       53cc1d36 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       e4fd2a28 Merge branch 'ss_issu_66_dnaseq_sv_schema' of bitbucket.org:mugqic/genpipes into ss_issu_66_dnaseq_sv_schema
       647d39a2 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       7b041d92 cleanup DS_Store file
       44db7072 Fixes #66 Dnaseq pipeline -t sv option workflow schema created
       48626df0 Fixes #60 ampliconseq -t dada2 workflow schema diagram created
       52a22866 Refs #53 covseq workflow arrow step 6 -> 9 removed
       faf94499 Refs #53 Paul's review inputs incorporated
       51870352 Fixes #53 covseq.py workflow diagram added
       bacd6ed2 Refs #53 added covseq pipeline schema workflow draft under review by Paul
       0f7b09a6 Fixes #51 update covseq pipeline readme to v3.3.0
       2fed25b4 Fixes #51 update covseq pipeline readme to v3.3.0

3.3.0        Fri Feb 19 16:37:40 2021 -0500        641 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      3 commits

       a4d9a14e GenPipes - removing tumor_apir before release
       1a28942a Version bump to 3.2.1-beta
       63d211df Version bump to 3.2.0

  ehenrion <edouard.henrion@mcgill.ca>      4 commits

       0756af3e Merged in release_3.3 (pull request #194)
       7b17e70b GenPipes - Bug fix : corrected dnaseq.cedar.ini
       d62b0c01 GenPipes - Bug fix : removed buggy line in dnaseq.cedar.ini
       30bdd0a4 GenPipes - Bug fix : correcting indentation in illumina_run_processing.py

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       ebef4a20 Release 3.2

  pnawarathna <pubudu.nawarathna@mail.mcgill.ca>      572 commits

       e6d1d8aa hicseq pipeline [changed the step order]
       165e0600 hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       498209ba hicseq [hicseq.py] - corrected output -o issue fastq_readset
       8d382dc5 hicseq [hicseq.py, readme.md] - modified readmefile
       638e16e1 hicseq [hicseq.py] - corrected file after rebase
       7fe877f4 hicseq [hicrep.py] - Corrected typo
       2c26336a hicseq [hicrep.py] - corrected R_TOOLS path
       91dd29c5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f44f83d7 hicseq [hicrep.py] - Corrected typo
       124ce0e3 hicseq [hicrep.py] - corrected R_TOOLS path
       678b9cf3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bb10a4ed hicseq [hicrep.py] - Corrected typo
       3e2b5502 hicseq [hicrep.py] - corrected R_TOOLS path
       8772abdb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       77cbfb90 hicseq [hicrep.py] - Corrected typo
       6f820f10 hicseq [hicrep.py] - corrected R_TOOLS path
       2528718f [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65191abb hicseq [hicrep.py] - Corrected typo
       28b43776 hicseq [hicrep.py] - corrected R_TOOLS path
       c82045f9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6b1b9d62 hicseq [hicrep.py] - Corrected typo
       be461ea1 hicseq [hicrep.py] - corrected R_TOOLS path
       26cf9254 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b06d2eaf hicseq [hicrep.py] - Corrected typo
       eff718e0 hicseq [hicrep.py] - corrected R_TOOLS path
       8c7b0eb2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d1f43e78 hicseq [hicrep.py] - Corrected typo
       97174e23 hicseq [hicrep.py] - corrected R_TOOLS path
       d82f6771 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       aa9e9215 hicseq [hicseq.py] - corrected file after rebase
       96a94204 hicseq [hicrep.py] - Corrected typo
       68f348bf hicseq [hicrep.py] - corrected R_TOOLS path
       284f0e1e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac094de hicseq [hicrep.py] - Corrected typo
       ec553bd9 hicseq [hicrep.py] - corrected R_TOOLS path
       e32d3bbf [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c1f326ef hicseq [hicrep.py] - Corrected typo
       4b99cbbc hicseq [hicrep.py] - corrected R_TOOLS path
       da58a44b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8985879c hicseq [hicrep.py] - Corrected typo
       0d1f31b6 hicseq [hicrep.py] - corrected R_TOOLS path
       d5aab522 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1d3a959b hicseq [hicrep.py] - Corrected typo
       2ec58eb5 hicseq [hicrep.py] - corrected R_TOOLS path
       aaf3aee3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2dd3de83 hicseq [hicrep.py] - Corrected typo
       0eb19568 hicseq [hicrep.py] - corrected R_TOOLS path
       c2627543 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       944ba8e7 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4021b7d7 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       e8567612 hicseq [quasar_qc - corrected module loading in matrix restructuring
       fc92143e hicseq [hicrep.py] - Corrected typo
       a13ed40c hicseq [hicrep.py] - corrected R_TOOLS path
       ee6bc35d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3b82c00f hicseq [hicrep.py] - Corrected typo
       9a509dcc hicseq [hicrep.py] - corrected R_TOOLS path
       edb308a6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ca9fa5df hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       d9b7dc4e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       aa8821b6 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       9419efd7 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8d61964c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       2e7870dd Completed developing hicrep and quasar analysis
       d54c390a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       3ca0a12a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       845d4734 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       42f947ad Completed hicrep analysis except the creation of the graph
       0274f322 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       1c1f819f hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e27ded70 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       53408440 hicseq completed adding basic features of the hicrep analysis.
       5b078cfe hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       94230e0e Added pairwise combination for samples
       39b34b00 started to edit the hicseq.py added new file for hicrep
       9e2040f4 hicseq [hicseq.py, readme.md] - modified readmefile
       44659521 hicseq [hicseq.py] - corrected file after rebase
       ccb344f8 hicseq [hicrep.py] - Corrected typo
       35df970d hicseq [hicrep.py] - corrected R_TOOLS path
       240f9c5d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       eecc3b29 hicseq [hicrep.py] - Corrected typo
       78e6c8ab hicseq [hicrep.py] - corrected R_TOOLS path
       616e8d37 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       84e8ec91 hicseq [hicrep.py] - Corrected typo
       b142de6e hicseq [hicrep.py] - corrected R_TOOLS path
       dd6cb48b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5bad885b hicseq [hicrep.py] - Corrected typo
       030d5584 hicseq [hicrep.py] - corrected R_TOOLS path
       e29ae96c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       97622c6d hicseq [hicrep.py] - Corrected typo
       a477d712 hicseq [hicrep.py] - corrected R_TOOLS path
       d05a6b62 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       12a4e2f4 hicseq [hicrep.py] - Corrected typo
       45695cbd hicseq [hicrep.py] - corrected R_TOOLS path
       e7041d31 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       302c4e96 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       154f9b27 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       8e55fb57 hicseq [quasar_qc - corrected module loading in matrix restructuring
       f1534e3b hicseq [hicrep.py] - Corrected typo
       45572379 hicseq [hicrep.py] - corrected R_TOOLS path
       c2063ea9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e41bb32f hicseq [hicrep.py] - Corrected typo
       ebbbf879 hicseq [hicrep.py] - corrected R_TOOLS path
       44063144 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       bd36df09 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       40765477 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       6f40129e [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       789fb24b [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       83ea31bd [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       0eafabda Completed developing hicrep and quasar analysis
       dc1cf6aa [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c0974eb5 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       82c023bb created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       0ea49fdb Completed hicrep analysis except the creation of the graph
       55f2b850 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bb223680 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       19aec506 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       7cd07859 hicseq completed adding basic features of the hicrep analysis.
       7bdf5b49 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       ae313c6d Added pairwise combination for samples
       08dae47d started to edit the hicseq.py added new file for hicrep
       eff71142 hicseq [hicseq.py] - corrected file after rebase
       b6931059 hicseq [hicrep.py] - Corrected typo
       386ca833 hicseq [hicrep.py] - corrected R_TOOLS path
       5a412f71 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       0b7fe2a2 hicseq [hicrep.py] - Corrected typo
       1765bef4 hicseq [hicrep.py] - corrected R_TOOLS path
       5851499d [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       1555fed4 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       d3199ece hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       01d239a9 hicseq [quasar_qc - corrected module loading in matrix restructuring
       46e7d3a5 hicseq [hicrep.py] - Corrected typo
       5e92f8e6 hicseq [hicrep.py] - corrected R_TOOLS path
       f2dbfd7b [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       4fa81e9b hicseq [hicrep.py] - Corrected typo
       b284455e hicseq [hicrep.py] - corrected R_TOOLS path
       b12fdbe5 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a78a623e hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       125a6e30 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       67272ac2 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       19f3edc7 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       a8b8d85c [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       9d37f638 Completed developing hicrep and quasar analysis
       acce3835 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       92bfb414 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1a3bcdf2 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       812ba6d0 Completed hicrep analysis except the creation of the graph
       27467f36 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3cc9f393 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2684371d hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       92744be2 hicseq completed adding basic features of the hicrep analysis.
       06d8e53b hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       31b01578 Added pairwise combination for samples
       6775791d started to edit the hicseq.py added new file for hicrep
       38a937ec hicseq [base.ini] - updated mugqic tools version to 2.3.1
       4e9dea90 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       6efe674b hicseq [quasar_qc - corrected module loading in matrix restructuring
       b918766e hicseq [hicrep.py] - Corrected typo
       0219397e hicseq [hicrep.py] - corrected R_TOOLS path
       187cba37 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       cac69b20 hicseq [hicrep.py] - Corrected typo
       a80a6bcb hicseq [hicrep.py] - corrected R_TOOLS path
       06e2631e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f903c4b6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ebc75da3 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       c9acbc74 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       dcbe50b4 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bd63e8e6 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4f9a3087 Completed developing hicrep and quasar analysis
       6f0d1a8a [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1052c8db [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       704d7f19 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       a7f8b2a9 Completed hicrep analysis except the creation of the graph
       f1fd9af6 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       712cb462 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       979429b8 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       3930cd3f hicseq completed adding basic features of the hicrep analysis.
       b3c1160a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       f493eae5 Added pairwise combination for samples
       6f98cb3a started to edit the hicseq.py added new file for hicrep
       f58f27f5 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       b576afd0 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       0eb37c3f hicseq [quasar_qc - corrected module loading in matrix restructuring
       b718a5b2 hicseq [hicrep.py] - Corrected typo
       5f8cf474 hicseq [hicrep.py] - corrected R_TOOLS path
       1fde1089 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6478d079 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       b262df5b [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b340c0d5 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       3b6b3578 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       d740ed8a [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c1b8eef7 Completed developing hicrep and quasar analysis
       6bc89a86 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       95c477e8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       53fc03ce created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       cf6a3af0 Completed hicrep analysis except the creation of the graph
       d8542242 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b961772 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       8438ff65 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       b18b1915 hicseq completed adding basic features of the hicrep analysis.
       c4ba0941 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bd650a8b Added pairwise combination for samples
       6f4136c3 started to edit the hicseq.py added new file for hicrep
       9364f1e4 hicseq [hicrep.py] - Corrected typo
       301afce0 hicseq [quasar_qc.py] - Corrected deletion by mistake
       bd2a8e12 hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       96761788 hicseq [hicrep.py] - corrected R_TOOLS path
       a5ca0bab [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       432f5862 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       3e68e18c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5c841f74 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       ab730be1 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       3bce8239 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       e1cee2e1 Completed developing hicrep and quasar analysis
       45b526c8 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6cf84507 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       21e2e6ce created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       34416892 Completed hicrep analysis except the creation of the graph
       6570bc97 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       46a28b49 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       645c2f89 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       0189b808 hicseq completed adding basic features of the hicrep analysis.
       26069513 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       6fd22e84 Added pairwise combination for samples
       d664e0e5 started to edit the hicseq.py added new file for hicrep
       1e51cf6d hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       4cca0610 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2f8c42d4 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       1c49e458 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cfa04f30 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c0bd3651 Completed developing hicrep and quasar analysis
       70d8bdfe [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6711a689 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e8effea3 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       5fe98a01 Completed hicrep analysis except the creation of the graph
       0bc9a5dd hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bea29b63 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       eb89bc04 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       d7bb65df hicseq completed adding basic features of the hicrep analysis.
       ccc376e4 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       93f957d3 Added pairwise combination for samples
       1c14ddf2 started to edit the hicseq.py added new file for hicrep
       0da3f4c4 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       9ae9d5ed [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       30a4e554 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       8a5ac66b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       4fb194fa Completed developing hicrep and quasar analysis
       7f2917d2 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2492118c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       9a8cfcd5 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       54c7c37f Completed hicrep analysis except the creation of the graph
       c35566cc hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fe1171a0 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       e9db27dd hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8584e29d hicseq completed adding basic features of the hicrep analysis.
       dccac155 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       56890510 Added pairwise combination for samples
       f609b2f9 started to edit the hicseq.py added new file for hicrep
       31ba309e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       55f77588 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       f794061c [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       84e93b1b [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       915d9ace Completed developing hicrep and quasar analysis
       769c97be [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       01ed5378 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       b8a79f4f created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       61632401 Completed hicrep analysis except the creation of the graph
       7efc2ba3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bf7412fd hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       14256a14 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       042082e5 hicseq completed adding basic features of the hicrep analysis.
       fa1e7d0a hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       eba9e468 Added pairwise combination for samples
       9bd17ff3 started to edit the hicseq.py added new file for hicrep
       99d62228 Completed developing hicrep and quasar analysis
       719115c3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       7c3222a9 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       966e814b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       ece06e33 Completed hicrep analysis except the creation of the graph
       4ceda720 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a650327 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       91f3b4b7 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       15e34d05 hicseq completed adding basic features of the hicrep analysis.
       92baadc3 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       d9d98a54 Added pairwise combination for samples
       c4853c36 started to edit the hicseq.py added new file for hicrep
       78ab0105 hicseq pipeline [changed the step order]
       c2ac146d hicseq pipeline [removed user comments, delete guillimin and mammouth ini files]
       8bf6a4d7 hicseq [hicseq.py] - corrected output -o issue fastq_readset
       a016c1eb hicseq [hicseq.py, readme.md] - modified readmefile
       7fbb757c hicseq [hicseq.py] - corrected file after rebase
       f8fc53cb hicseq [hicrep.py] - Corrected typo
       46c6016a hicseq [hicrep.py] - corrected R_TOOLS path
       fb83a924 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       00a5b87f hicseq [hicrep.py] - Corrected typo
       3e840af9 hicseq [hicrep.py] - corrected R_TOOLS path
       c61122bc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d3659cb4 hicseq [hicrep.py] - Corrected typo
       230b7b68 hicseq [hicrep.py] - corrected R_TOOLS path
       9a454ab8 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       857f8c9a hicseq [hicrep.py] - Corrected typo
       5c402450 hicseq [hicrep.py] - corrected R_TOOLS path
       ba545bd4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       3ff1e509 hicseq [hicrep.py] - Corrected typo
       a0098d0e hicseq [hicrep.py] - corrected R_TOOLS path
       2c820d82 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d23933c7 hicseq [hicrep.py] - Corrected typo
       b0645643 hicseq [hicrep.py] - corrected R_TOOLS path
       54bb894c [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b81d7550 hicseq [hicrep.py] - Corrected typo
       74e48a25 hicseq [hicrep.py] - corrected R_TOOLS path
       40fe0ccb [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a3680d65 hicseq [hicrep.py] - Corrected typo
       798c19f6 hicseq [hicrep.py] - corrected R_TOOLS path
       06ec7886 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       2fd36b1d hicseq [hicseq.py] - corrected file after rebase
       be817b44 hicseq [hicrep.py] - Corrected typo
       5248d325 hicseq [hicrep.py] - corrected R_TOOLS path
       39617b23 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec979601 hicseq [hicrep.py] - Corrected typo
       24e7c7c9 hicseq [hicrep.py] - corrected R_TOOLS path
       4d6d3f88 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       c17d5f12 hicseq [hicrep.py] - Corrected typo
       aeae41c1 hicseq [hicrep.py] - corrected R_TOOLS path
       77dabaf9 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       032cfe90 hicseq [hicrep.py] - Corrected typo
       07a72bd1 hicseq [hicrep.py] - corrected R_TOOLS path
       1039f672 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       590a6004 hicseq [hicrep.py] - Corrected typo
       8f9b5466 hicseq [hicrep.py] - corrected R_TOOLS path
       574336c6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       6c20f492 hicseq [hicrep.py] - Corrected typo
       ff45f2df hicseq [hicrep.py] - corrected R_TOOLS path
       85d661e3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ec67234d hicseq [base.ini] - updated mugqic tools version to 2.3.1
       463da45c hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       7701580c hicseq [quasar_qc - corrected module loading in matrix restructuring
       0c5ee2fb hicseq [hicrep.py] - Corrected typo
       1a221b48 hicseq [hicrep.py] - corrected R_TOOLS path
       8966cfbc [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b02d1a7f hicseq [hicrep.py] - Corrected typo
       4d37159f hicseq [hicrep.py] - corrected R_TOOLS path
       22fbd606 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5541fbc2 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       fdeb31d1 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       d58d40f2 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       fa1f4b49 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       cdcc7906 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7b9f6466 Completed developing hicrep and quasar analysis
       31c4df57 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c1cc394c [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       80dc0ab5 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       72cbd003 Completed hicrep analysis except the creation of the graph
       0a938267 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       85cc99e5 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2e455782 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       46f2fa31 hicseq completed adding basic features of the hicrep analysis.
       8c3e8d90 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8e5c6545 Added pairwise combination for samples
       55862cd8 started to edit the hicseq.py added new file for hicrep
       e79b132d hicseq [hicseq.py, readme.md] - modified readmefile
       ffea320a hicseq [hicseq.py] - corrected file after rebase
       fed62535 hicseq [hicrep.py] - Corrected typo
       d285ff0e hicseq [hicrep.py] - corrected R_TOOLS path
       e04cb392 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       de2d8e17 hicseq [hicrep.py] - Corrected typo
       9c6f2d21 hicseq [hicrep.py] - corrected R_TOOLS path
       a8ae3ab4 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       866b11d8 hicseq [hicrep.py] - Corrected typo
       82db1223 hicseq [hicrep.py] - corrected R_TOOLS path
       4b9e5e65 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       dda53524 hicseq [hicrep.py] - Corrected typo
       452f1db6 hicseq [hicrep.py] - corrected R_TOOLS path
       db683cc7 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       f2e36129 hicseq [hicrep.py] - Corrected typo
       2accb63c hicseq [hicrep.py] - corrected R_TOOLS path
       f8413f07 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       7938420d hicseq [hicrep.py] - Corrected typo
       e78e318b hicseq [hicrep.py] - corrected R_TOOLS path
       80eff204 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       a0b00628 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       1160fe25 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       3ba11742 hicseq [quasar_qc - corrected module loading in matrix restructuring
       b27f4021 hicseq [hicrep.py] - Corrected typo
       01ac15da hicseq [hicrep.py] - corrected R_TOOLS path
       80ad7082 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       ef6d0348 hicseq [hicrep.py] - Corrected typo
       045b540a hicseq [hicrep.py] - corrected R_TOOLS path
       f6422cd3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       8808f055 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       65f84c4e [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7adcfa1e [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       14c0a789 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       935862a9 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       fd6eef8a Completed developing hicrep and quasar analysis
       8c8248e4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       85f0d9ef [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       dc225b23 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       dfe7725f Completed hicrep analysis except the creation of the graph
       db50a5d9 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       10222fe0 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       5dcda07c hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       062551bf hicseq completed adding basic features of the hicrep analysis.
       766c1b49 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       e979eee4 Added pairwise combination for samples
       cfff9164 started to edit the hicseq.py added new file for hicrep
       4443e3d2 hicseq [hicseq.py] - corrected file after rebase
       ca8b09db hicseq [hicrep.py] - Corrected typo
       584dd184 hicseq [hicrep.py] - corrected R_TOOLS path
       afa77dc2 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5c2535f6 hicseq [hicrep.py] - Corrected typo
       8db0da7e hicseq [hicrep.py] - corrected R_TOOLS path
       6f264bd6 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       5960be4c hicseq [base.ini] - updated mugqic tools version to 2.3.1
       537be829 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       1313e2b7 hicseq [quasar_qc - corrected module loading in matrix restructuring
       9bcc430c hicseq [hicrep.py] - Corrected typo
       8cacfa8a hicseq [hicrep.py] - corrected R_TOOLS path
       c8c3cef0 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       71b02a39 hicseq [hicrep.py] - Corrected typo
       b98db26d hicseq [hicrep.py] - corrected R_TOOLS path
       3698d950 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       d461e1e8 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       5c728c64 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       7a9cfa6b [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       2b665729 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       312c08ea [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       7e6f9aa0 Completed developing hicrep and quasar analysis
       334abafc [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       87973bb4 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       03aab022 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8359e46d Completed hicrep analysis except the creation of the graph
       8a845b7b hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       eee8a212 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       9a55b5c2 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       2755cba7 hicseq completed adding basic features of the hicrep analysis.
       3c517c46 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       8852df1f Added pairwise combination for samples
       54e5f9fd started to edit the hicseq.py added new file for hicrep
       9d541009 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       95e4b22f hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       d63f7b4d hicseq [quasar_qc - corrected module loading in matrix restructuring
       56c0d8c4 hicseq [hicrep.py] - Corrected typo
       80acf528 hicseq [hicrep.py] - corrected R_TOOLS path
       7cf9e946 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       65bafda3 hicseq [hicrep.py] - Corrected typo
       990b8fae hicseq [hicrep.py] - corrected R_TOOLS path
       b86b171e [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       e5fdc3a9 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       ff495116 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       5da9af43 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       8113e902 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       845ec626 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       09a27509 Completed developing hicrep and quasar analysis
       1e1ab4c3 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       4f0cba11 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6ba3f9dd created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       be91a406 Completed hicrep analysis except the creation of the graph
       e2189f7c hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       fc54a3a2 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       09ec9842 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       c58406d4 hicseq completed adding basic features of the hicrep analysis.
       36dd6fab hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       1861910a Added pairwise combination for samples
       ef297a33 started to edit the hicseq.py added new file for hicrep
       5d7c80f9 hicseq [base.ini] - updated mugqic tools version to 2.3.1
       868a7110 hicseq [base.ini] - updated mugqic tools version, [other inis] - updated walltime for mergeing quasr stats
       9685043d hicseq [quasar_qc - corrected module loading in matrix restructuring
       f5f88d63 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu Corrected conflicts Conflicts: 	bfx/quasar_qc.py 	pipelines/hicseq/hicseq.base.ini 	pipelines/hicseq/hicseq.py
       5a838602 hicseq [hicrep.py] - Corrected typo
       b9e6c3b8 hicseq [hicrep.py] - corrected R_TOOLS path
       6d4ce4f3 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       b89f8e74 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       8beb5907 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       b78f92c1 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       18523b87 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       0ff0ff7f [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       c6f98d60 Completed developing hicrep and quasar analysis
       0b10bc0b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       8b60a2bc [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       6f50b29b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       8faa8ebf Completed hicrep analysis except the creation of the graph
       0bd1d0a3 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       3b6ab90b hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       4a7bbcaf hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       74e0bcb9 hicseq completed adding basic features of the hicrep analysis.
       9cce0060 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       183d5a95 Added pairwise combination for samples
       7f7301fd started to edit the hicseq.py added new file for hicrep
       0f33a580 hicseq [hicrep.py] - Corrected typo
       c47a6eec hicseq [quasar_qc.py] - Corrected deletion by mistake
       8b33c6fe hicseq [hicrep.py, hicseq.py, quasar_qc.py] - Added further comments to making easy to understand the code
       61be1c4d hicseq [hicrep.py] - corrected R_TOOLS path
       20c3b4a8 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       edda3a12 [hicseq.py, quasar_qc.py,hicrep.py] - modified to put all temporary files saved in to temp folder and defined them as removable files Added a new merging step to merge all quasar_qc report files in to one file Corrected a typo in hicrep output files
       18506b49 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       75cbaaa9 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       2a0e1a8a [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       15bbbb37 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       44735d2f [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       eccd7c75 Completed developing hicrep and quasar analysis
       e9078d81 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       2ac9239f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       542bad0b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       418c5156 Completed hicrep analysis except the creation of the graph
       f587c742 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       5a9980d3 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       2d870202 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       8df5abcd hicseq completed adding basic features of the hicrep analysis.
       f55e1969 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       534d63fd Added pairwise combination for samples
       756599e1 started to edit the hicseq.py added new file for hicrep
       b45e44d6 hicseq pipeline - Added parameters for reproducibility and quality scores steps of all ini files
       cf5a2f40 [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       56985bb0 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       b48157cd [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       bdbc7de8 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       d1c1ae9b Completed developing hicrep and quasar analysis
       e68d3079 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       c832c491 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1b63533a created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       2c192d1d Completed hicrep analysis except the creation of the graph
       e081f172 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       6f3560e4 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       35efcde1 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5d6fcc13 hicseq completed adding basic features of the hicrep analysis.
       ca0d1635 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       72165318 Added pairwise combination for samples
       396a06d9 started to edit the hicseq.py added new file for hicrep
       35412d6c [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       97ce4746 [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       38b339b3 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       53603e16 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       db3da4bf Completed developing hicrep and quasar analysis
       ba7bbb87 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       e6625b9b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       a7978565 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       720a6263 Completed hicrep analysis except the creation of the graph
       b5cbb3fe hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       bfbb8968 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       1099ed04 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       5f5f1685 hicseq completed adding basic features of the hicrep analysis.
       f3031630 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       bbbcf3ef Added pairwise combination for samples
       a5c29b57 started to edit the hicseq.py added new file for hicrep
       292a0b2b [hicseq.py] - changed the output directory name for reproducibilty score [hicrep.py] - changed the path for hicrep.R
       cd9d16ad [btx/quasar_qc.py]-removed sed 's/ *//g' and added additional awk to print values since there was a text passing issue. Need to be optimized further.
       83f719b4 [btx/quasar_qc.py]-added trunc() to remove other decimals than the first one
       4c23bf13 [btx/quasar_qc.py]-replaced awk with Rscript and improved code. now data frames only with decimal values are multiplied. missing values were removed and replaced with 0
       5f42b4c2 Completed developing hicrep and quasar analysis
       849b382b [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       d9ed6481 [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       1767a32b created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       068daa9d Completed hicrep analysis except the creation of the graph
       210ad291 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       411d8e10 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       0f31f813 hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       79aff176 hicseq completed adding basic features of the hicrep analysis.
       b5075cb2 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       57e8baef Added pairwise combination for samples
       4a57c1a7 started to edit the hicseq.py added new file for hicrep
       89006b6c Completed developing hicrep and quasar analysis
       87eb470f [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       477281af [hicseq.py, quasar_qc.py] completed the job for creating a fend object for the project
       10c5eb18 created new file for QUASAR-QC in directory btx and defined methods in hicseq.py for qusarqc analysis. obtained path for chromosome lengths file
       335f4c52 Completed hicrep analysis except the creation of the graph
       10daca13 hicseq [hicseq.py] - changed 23 parallel jobs for each chromosome to one job for all 23 chromosomes.
       4aeab549 hicseq [hicseq.py, hicrep.py] - modified to handle multiple paramters for hicrep and create combined csv file for each sample comparison as a new output
       509fef1c hicseq [hicseq.py, hicrep.py] - modified for additional user parameters, optimal h finding, down-sampling, save weight and corr matrices
       1f096d18 hicseq completed adding basic features of the hicrep analysis.
       b623db54 hic-seq [hic] hicseq.py Added pairwise comparison and added jobs to hicrep.py - include testing lines
       c6057427 Added pairwise combination for samples
       9ad63fde started to edit the hicseq.py added new file for hicrep

  pubudumanoj <pnawarat@abacus3.ferrier.genome.mcgill.ca>      16 commits

       a1ddafb1 Corrected hicseq.py file deletions
       e011fbe4 Corrected hicseq.py file deletions
       9c40c037 Corrected hicseq.py file deletions
       d8ded448 Corrected hicseq.py file deletions
       0a0903c5 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       1fdc6405 Corrected hicseq.py file deletions
       cbf6c692 Corrected hicseq.py file deletions
       96e5d3bd Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       41a5bc29 Corrected hicseq.py file deletions
       3ee012f3 Corrected hicseq.py file deletions
       228a9e46 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       9211f3c9 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       3f918251 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected merge conflict after rebase Conflicts: 	pipelines/hicseq/hicseq.py
       6291f898 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu corrected module loading
       296b86bf Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a757b8df Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  pubudumanoj <pubudu@gra-login1.graham.sharcnet>      20 commits

       925b082a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       3d8131d0 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       ae8f296e hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       21e76652 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       938c27c9 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       38516f39 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       08b13811 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       4f835f9b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       402a4d1a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       6b60d129 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       dd98639a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2d3a5561 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       cd42f4e5 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       26313486 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       996ba98b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       d384dc9b hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       5ce6c34d hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       bfd60bcd hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       c743b445 hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini
       2305332a hicseq hicseq.py hicrep.py hicseq.base.ini competed the hicrep basic running steps. now it calls the R script added new wall time on hicseq.base.ini changed bioconductor version on hicseq.base.ini changed python version on hicseq.base.ini

  pubudu.nawarathna@mail.mcgill.ca <pnawarat@abacus3.ferrier.genome.mcgill.ca>      1 commits

       e5fa3455 Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu

  Pubudu Nawarathna Mudiyanselage <pubudu@beluga3.int.ets1.calculquebec.ca>      20 commits

       1608e2b1 refer previous commit
       a749a350 refer previous commit
       ce31a37f refer previous commit
       4e8a71e4 refer previous commit
       2ba43453 refer previous commit
       5c06e434 refer previous commit
       bcac0989 refer previous commit
       1f855354 refer previous commit
       04967444 refer previous commit
       f0eb790e refer previous commit
       afce9450 refer previous commit
       fb33ca90 refer previous commit
       45984451 refer previous commit
       ba1f6b0f refer previous commit
       9227c2a4 refer previous commit
       3b5869f7 refer previous commit
       daccd079 refer previous commit
       34642e8d refer previous commit
       4bd83cec refer previous commit
       4ccc1064 refer previous commit

  Pubudu Nawarathna Mudiyanselage <pubudu.nawarathna@mail.mcgill.ca>      3 commits

       3c7dbddc Merge branch 'hicseq_hicrep_pubudu' of https://bitbucket.org/mugqic/genpipes into hicseq_hicrep_pubudu
       a9fe3dd4 corrected md file
       fc8ac1ad corrected md file

  Pubudu Nawarathna <pubudu.nawarathna@mail.mcgill.ca>      1 commits

       1114ea8a Merged in hicseq_hicrep_pubudu (pull request #170)

3.2.0        Mon Jan 25 12:47:42 2021 -0500        371 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      5 commits

       ed04f330 testing end-of-line character
       957c11db Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       057a8d13 GenPipes - Genomes : updated R installation script to ease installation in dev, also corrected the hicup module called in install_genome.sh
       690f56e1 Version bump to 3.1.6-beta
       369cad4a Version bump to 3.1.5

  ehenrion <edouard.henrion@mcgill.ca>      5 commits

       1ce85b32 GenPipes : illumina_run_processing.py : correcting indentation
       f2d7e724 GenPipes - Tumor Pair pipeline : removed dev vawk module in tumor_pair.base.ini
       62ac6797 GenPipes - DNASeq pipeline - removing mugqic_dev modules in gatk4.ini
       39ceaa25 GenPipes - Tumor Pair pipeline : removing mugqic_dev modules
       a84f5aa2 Merged in genpipes-3.1.5_release (pull request #166)

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      4 commits

       42981f13 Merged in rnaseq_light_docs (pull request #183)
       a190faaa rnaseq_light.py edited to adjust docstrings to address issue raised by Shaloo here : https://github.com/c3g/GenPipes/issues/63
       ceb0c8bb Merged in JoseHectorGalvezLopez/nanoporebaseini-edited-online-with-bitbu-1596652406491 (pull request #179)
       35bc91bb Edited the nanopore ini file to address the mugqic_tools error.

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       52f9015b modify cedar ini
       f38492f7 a bit ugly resolution from argparse overriding issue of the type argument

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      154 commits

       41671c60 Merged in mgi_stretenp (pull request #192)
       63a54683 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       a3cffe0a Changing qualimap bamqc section naming to a parameter in bfx
       6e33047a Changing qualimap bamqc section naming to a parameter in bfx
       ad67e66e Merged in mgi_stretenp (pull request #191)
       356bd208 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       dcbf891d igvtools ressources change
       f5026a89 Reducing ressources
       f15a4f8e Fixing awk
       c74a4c67 Changing default ressources
       fede180d igvtools ressources change
       58c298c2 Reducing ressources
       29ec475f Fixing awk
       cfca926e Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       66189928 Changing default ressources
       9f663b00 Changing default ressources
       e74c221f Merged in mgi_stretenp (pull request #188)
       21f213f0 Fix
       cec6282b Fix
       907ab2a0 Fix interval_list checking
       08c1c549 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       9ba27117 Checking if interval_list is already on genome folder
       b6f2b919 Merged dev into mgi_stretenp
       247697b6 samtools bam2fq typo
       92d206f9 Adding kraken to beluga ini
       bb375b9b Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       43d126eb Switching to ivar 1.3 and NML dehosting filters
       7b86a544 Updating ivar trim usage
       d04bd81e Adding htslib to ivar module
       3ea8a2e5 Fixing kraken output and adding warning on rnaseq htseq-count
       9f4fde94 Switching to latest kraken2 db build
       cad6ed3e Forcing pigz to overwrite
       522905b7 Switching to latest kraken2 release
       de8e702b Fixing kraken
       e2249ccd Fixing kraken bfx
       28dcdf30 Fixing kraken bfx
       63339cd2 Fixing kraken bfx
       0f47d1e7 Fixing kraken bfx
       338b71d5 Fix kraken
       3f5c373c Fix kraken
       d13c20d0 Addinf kraken analysis for metrics
       272131bb Update sambamba sort ini
       69451045 Fix rename ln test
       940b62c6 Fix
       3ac107c7 Fixing tee
       624b2c23 Test
       c469a365 Fix
       9993fca8 Adding metrics on dehosted bam
       ff29ecf6 renaming cit covseq file
       01e6e6c2 cit ini for cit test
       67b78d6a Fixing hybrid genome path and output selection
       5d8d2f7e renaming cit covseq file
       15ba0d21 cit ini for cit test
       f6f9782c Switching to ivar 1.3 and NML dehosting filters
       4e16c725 Updating ivar trim usage
       4a07fed8 Adding htslib to ivar module
       05475b4a Fixing kraken output and adding warning on rnaseq htseq-count
       82ba0847 Switching to latest kraken2 db build
       05e7d80c Forcing pigz to overwrite
       a11ec96c Switching to latest kraken2 release
       e43b7fb3 Fixing kraken
       2f5805b7 Fixing kraken bfx
       879ee559 Fixing kraken bfx
       19b050c9 Fixing kraken bfx
       2cb88182 Fixing kraken bfx
       eb614ee5 Fix kraken
       74227cca Fix kraken
       c797c2aa Addinf kraken analysis for metrics
       52b1110f Update sambamba sort ini
       a06f5886 Fix rename ln test
       398fcefc Fix
       012b267d Fixing tee
       54b5fa5b Test
       33201b80 Fix
       516add15 Adding metrics on dehosted bam
       d5b321a2 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       0b9b56a0 Reducing sambamba filtering default cpus
       77cddeb1 Fixing rename consensus
       916b0912 renaming cit covseq file
       a5509539 cit ini for cit test
       67a81742 Changing ini sections to use Illumina beds files as default
       fcfce997 Changing default seq_method
       84afdaa0 Cleaning and switching to cvmfs genome
       ed0fc196 Fix cit
       409ffcd2 Cit fix
       ccd50ffb Cit test
       c999d20c Fix
       24fd48c3 Fix
       a11172c3 Fix
       03180fbb Fix rename consensus symlink
       7533de7b Fixing rename consensus
       d33ab92b Fixing tsv
       b09b78e0 Fixing tsv
       2ab5bd4a Fixing tsv renaming
       53954899 Fixing output rename header + tsv for ncov-tools
       b0eb526d Fixing picard metrics
       13a5d310 Collecting picard metrics on raw AND filtered bam
       95dadb9b Fixing select output file
       8cfc8ea8 Fixing hybrid genome path and output selection
       35c8d492 Fixing merging step for 1 sample with multiple readsets
       210ebdd0 Update inis
       fe84c0b0 quast -> Quast
       13c60bf9 renaming cit covseq file
       9704c198 cit ini for cit test
       11fd8090 Reducing sambamba filtering default cpus
       5843d7bb Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       1fccae9a Changing ini sections to use Illumina beds files as default
       4c321cab Changing default seq_method
       ad84e075 Cleaning and switching to cvmfs genome
       5b66a1b3 Fix cit
       8ac2d785 Cit fix
       aa5609b0 Cit test
       a5a46a53 Fix
       360e7a1c Fix
       46bc3964 Fix
       93fa4525 Fix rename consensus symlink
       c48a2989 Fixing rename consensus
       816eb579 Fixing tsv
       24e97343 Fixing tsv
       679da0c2 Fixing tsv renaming
       52b8e1b8 Fixing output rename header + tsv for ncov-tools
       8727d67a Fixing picard metrics
       c6a1a2d1 Collecting picard metrics on raw AND filtered bam
       860d44bf Fixing select output file
       26be1581 Fixing hybrid genome path and output selection
       facc8cba Fixing merging step for 1 sample with multiple readsets
       1442f6a1 Update inis
       126d53dd quast -> Quast
       3a7504d5 renaming cit covseq file
       efac275b cit ini for cit test
       412da130 Changing ini sections to use Illumina beds files as default
       24ccf1f9 Changing default seq_method
       f4d97e9c Cleaning and switching to cvmfs genome
       3c8daf00 Fix cit
       e6e97385 Cit fix
       3ba059d9 Cit test
       7a29d012 Fix
       e6a86865 Fix
       e179616b Fix
       e528c792 Fix rename consensus symlink
       4de70768 Fixing rename consensus
       54301693 Fixing tsv
       2952b1ac Fixing tsv
       4ca82a10 Fixing tsv renaming
       f7df72ad Fixing output rename header + tsv for ncov-tools
       00b60522 Fixing picard metrics
       c009672a Collecting picard metrics on raw AND filtered bam
       fc6322f3 Fixing select output file
       a33b34c5 Fixing hybrid genome path and output selection
       75332a91 Fixing merging step for 1 sample with multiple readsets
       4a4b67eb Update inis
       527a801a quast -> Quast
       6bc79c0a renaming cit covseq file
       d467187e cit ini for cit test

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      5 commits

       406e38b9 Merged in fix_monitor (pull request #190)
       e885338e Merged in monitor_bug (pull request #189)
       af56bf44 Merged in chunk_slurm_submit (pull request #185)
       bfe072ad Merged in chunk_slurm_submit (pull request #182)
       b10eac60 Merged in rsync_in_chipseq (pull request #180)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      24 commits

       67c0af4b Merge branch 'release_3.2' of bitbucket.org:mugqic/genpipes into release_3.2
       90052a65 update readme for container install
       1294b963 Remove tumor_pair pipeline from release
       04c0c5a3 revert on indel aligner
       d75d4584 cleanup ini for beluga
       cd125513 remove dbSNP for oxog ini block fix gatk 4 module name
       7d84e335 fix picard_collect_oxog_metrics dbSNP vcf
       97d6a47f remove HOME_DEV paths
       70d13109 make job2json more robust
       a05a2430 remove mugqic_dev vawk
       7baea329 Fix README conflict
       3bda6d70 Fix out.out in monitor
       f0b99b86 rename to monitor.sh
       1ec8ba81 update control batch description
       bd24c9af rename deley_sbatch script
       c0f80124 fix script usage
       c579431b add export to sourced files
       f00a2d12 make sure already submited jobs id are sourced
       0a4444d4 if error
       28586080 remove debbug line
       5f6e768b add delay and chunking script
       88026d50 update for cit testing
       92edc027 replace -a by -r in rsync
       2551c787 fix master conflic with deleted readme file

  P-O Quirion <pioliqui@gmail.com>      2 commits

       0e486370 tweek monitor and chunk
       4f50d9e8 Fix autorestart and cleanup on failes interrupts

  Robert Eveleigh <eveleigh@beluga1.int.ets1.calculquebec.ca>      26 commits

       4f04bf4e covseq mpileup command fix
       acfe7b9a covseq qualimap fix, and high coverage ini fix
       9e0c8d79 add covseq - dnaseq consistency
       c6e5099e dnaseq high coverage trimmomatic to skewer
       8409952a tumor_pair fixes
       0c79ba9f import fixes
       5fb88cb6 bam.bai to .bai fix
       8894fc91 cit fixes to b38 - samtools and other b38 annotations
       87b0da86 conflict fixes to dnaseq and tumor pair after dev merge
       5a4dc7c4 fix mpileup dependency chain
       f8cda00c corrections to mpileup bcftools merge and tumor_pair dependencies
       b8ee4800 bam.bai to .bai fix
       4c7b3f7a fixes to indentation in sambamba merge
       cb807fc6 fixes to sambamba merge
       3bc8f0ed further gatk4 hc fixes
       0cdd2dfc update gatk4 hc arguments with gatk4 suite update
       8d7b9b82 fixes to gatk4 mark dup
       e7029c61 gatk4 mark dup walltime fix
       1d499545 variant recal fix
       9095ef11 variant recal dependency fix
       927011a3 cit fixes to b38 - samtools and other b38 annotations
       03ee662c minor dnaseq.py fixes
       4ffeb590 conflict fixes to dnaseq and tumor pair after dev merge
       cd935c7f fix mpileup dependency chain
       ef39fb7b corrections to mpileup bcftools merge and tumor_pair dependencies
       4d72c00d conforming deliverables to cit conventions

  Robert Eveleigh <eveleigh@beluga2.int.ets1.calculquebec.ca>      29 commits

       91f76551 exome specific fixes
       e6d4caf0 updates and fixes for cit
       72d072d4 remove briaree ini and update dnaseq base
       1547aa72 updates to beluga.ini and base.ini for dnaseq
       8e64f79c ini updates
       c925a612 gatk4 vsqr cit fix and baf plot for b38
       dd8c8447 add cram to input selection
       8b2d9c30 multiqc fix
       55d8a64c argument fixes for picard imported functions in gatk4 and vqsr fixes
       12f55ac8 indel realignment consistency issues between dnaseq and tumor_pair
       1c6f1bcb add mark dup metric file to multiqc
       4964f165 cit b38 fixes
       01f6a48e fix to bash.ln for cit
       2e1a585c fixes to multiqc
       d5cbe666 exome specific fixes
       e407f1d8 updates and fixes for cit
       a3bd732a Updates to bcftools/samtools for dnaseq_mpileup
       d9627230 cit fixes to dnaseq and test with real wes data, fixes to dependencies tumor_pair
       fc430949 remove briaree ini and update dnaseq base
       72e8439f adding 1 job processing to specific steps for cit.ini
       f4911b3f fix of dev genome reference
       2447075c issues between dnaseq.base.ini and dnaseq.beluga.ini
       cc7c806c updates to modules for beluga
       d5f03015 updating beluga ini
       67758960 added tumor pair beluga ini
       16bcac2f updates to beluga.ini and base.ini for dnaseq
       8e8f22e9 updates to beluga.ini
       f8ea2637 ini updates
       de84a173 module updates

  Robert Eveleigh <eveleigh@beluga3.int.ets1.calculquebec.ca>      15 commits

       369bd8a4 picard2 and high coverage fixes
       84885dae updates to b38 variant recal files
       3db1c5af fixes to tumor_pair on beluga
       a40514c6 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv
       e04c3f88 fixes to fixmate input file name
       9fe39b09 picard2 and high coverage fixes
       30cf9395 tumor_pair beluga ini fix
       c77c8210 tumor_pair qualimap part 2
       5c039c64 qualimap tumor_pair fix
       b980f2a0 fixes to vsqr gatk4
       ed996034 gatk4 fixes to callable loci and DoC
       6162e8e4 sym link dnaseq.base into tumor pair
       950a2854 updates to b38 variant recal files
       f6c42cba fixes to tumor_pair on beluga
       d8eaf582 cit dnaseq/tumor pair optimizations and fixes to mpileup and germline sv

  Robert Eveleigh <eveleigh@beluga4.int.ets1.calculquebec.ca>      24 commits

       2f926d67 fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       800ecb25 cit fixes after rebasing
       a37298c3 dependency fixes
       ba400c9b updates to GRCh38 annotation file, module updates, cit fixes
       01f02496 updates to cit and fixes to one job mpileup steps
       3e213bde updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8b745c15 major fixes to deliverables and completion of beluga test
       3078a9ac fix modules dev to cvmfs
       2345cafc fixes to chipseq, rnaseq_cufflinks, rnaseq_stringtie, and dnaseq_high_coverage
       3c284035 fixes to merge_filter_bcf
       662e007a gatk4 bsqr fixes and mpileup cat
       0f3aefcf mpileup protocol fix
       aefbb370 cit fixes after rebasing
       6f47ad31 update to gatk4 mutect2 filtering procedures
       419f376d create cit for gatk4 due to deviation from argument usage
       04fcaacc cit fixes to gatk4 + sym links for recalibration
       802cc658 dependency fixes
       7c3b8d3b updates to GRCh38 annotation file, module updates, cit fixes
       28828267 fixes to deliverable and b38 ini
       81204f69 updates to cit and fixes to one job mpileup steps
       dfa32887 updated wrapper bash commands to use bash_cmd and fixed indel realignment dependency bug
       8bc536d3 fixes to metasv annotations
       556e48b1 major fixes to deliverables and completion of beluga test
       973087a7 fix modules dev to cvmfs

  Robert Eveleigh <eveleigh@beluga5.int.ets1.calculquebec.ca>      1 commits

       b5c14868 final PR fixes

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      7 commits

       9ac20daa code cleaning and fixes to exome interval list
       f55a7384 fixes to cedar ini
       deca50fb fixes to symlinks for paired indel realignment
       31525eea code cleaning and fixes to exome interval list
       fa7328f4 fixes to cedar ini
       2d3d3676 fixes to symlinks for paired indel realignment
       5475ca4a cedar ini and exome update

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      5 commits

       7fbd4f9f fixes to metasv, adding metasv germline
       5a6a9d75 cedar fixes and GRCh38 fixes
       b42ed166 cedar dnaseq updates and svaba germline added
       2e389e92 cedar germline sv updates
       67a2e47d sequence dictionary and module updates

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      35 commits

       7d71caf8 updates to metasv - somatic
       eada99c8 updates to SV germline and reference genome tweaks
       0aa9d42b somatic sv fixes: lumpy and svaba annotations
       0f277902 fix to germline SV: breakseq2
       df8096ac merge fixes
       a96622a5 json and folder updates
       ee184e24 fixes to sCNAphase
       fa8f99fc merging snv and sv, adding protocols
       6cae335e Fixes to indel realigner
       16099054 Updates and debug
       04229bad Add set somatic and actionable mutations
       937e5438 added multiqc and other tweaks
       2be1a974 add metrics step for metasv
       a9697c7b updates to metasv - somatic
       f622b922 Fixes and updates to reference files
       909c2efb remove testing steps
       0cb7a4b5 updates to SV germline and reference genome tweaks
       de92d575 somatic sv fixes: lumpy and svaba annotations
       573c45dd fix to germline SV: breakseq2
       f2fc9d66 merge fixes
       77bc1190 GATK4 fixes - bam indexing and markDupSpark
       21051bce bcftools fixes for tumor pair
       44f3d044 fingerprint and bug fixes
       c9c20495 dnaseq - vcftools qc addition: --missing_indv and --depth
       9428baa6 select input for variant caller and fixes to one job calling
       1c78b285 json and folder updates
       3731f205 fixes to sCNAphase
       402bccbe Added json info
       036606a8 Bug fixes prior to json additions
       bd6e6cb9 merging snv and sv, adding protocols
       d441d471 Fixes to indel realigner
       c07293c0 Add deliverables module
       f674ec9f Updates and debug
       3161cf6d Add set somatic and actionable mutations
       e8d90db0 added multiqc and other tweaks

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      14 commits

       23f5be7b remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       614ee8f4 fixes to metasv for tumor pair
       9f1b500b Bug fixes and modification derived from initial PROFYLE benchmarking
       245df59d remove gatk2 cat_variant in favour of picard2/gatk4 mergeVcfs
       079fd7e6 fixes to metasv for tumor pair
       48bb3c5e Single job bug fixes
       b77cfaaa manta I/O fix and other bug fixes
       d14ba6d7 config updates and b38DH added
       6b528790 dnaseq germline SV updates
       b82f0f67 gatk4 updates and bug fixes
       14310b44 gatk4 updates and bug fixes
       77126e8e Fix line break types
       9fe9453e Json related bug fixes
       dcf0b8ac Bug fixes and modification derived from initial PROFYLE benchmarking

  robert.eveleigh@mcgill.ca <reveleig@abacus3.ferrier.genome.mcgill.ca>      8 commits

       04b54901 fixes to breakseq2 and metasv
       b0626eca cit-based fix - verifyBAMid
       b897ef1d dnaseq qc additions: NGScheckmate and peddy
       0cbfb0df fixes to breakseq2 and metasv
       612c4466 cit-based fix - verifyBAMid
       52023e1c cit-based fixes to NGScheckmate
       4dec63d8 dnaseq - multiqc fix
       29e26f72 dnaseq qc additions: NGScheckmate and peddy

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      4 commits

       f3c972e6 Merged in rebasing_tp (pull request #186)
       f6e74bd2 Debugging and Guillimin specfic fixes
       4810e9a4 Debugging and Guillimin specfic fixes
       0d16e4dc updates to config

  Rom Grk <romgrk.cc@gmail.com>      2 commits

       1d5912c7 Merged in fix-watch-portal-folder (pull request #187)
       c22599c5 fix: use of undeclared variables

covid_1.0        Mon Aug 3 19:56:47 2020 -0400        781 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      19 commits

       a24b8db7 GenPipes Covid Release - Renamed README.md as INFO.md and added GenPipes-like README.md for CovSeq pipeline
       73a694f3 Version bump to Covid Release 1.0
       eda88d0e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       3555e1b9 Version bump to Covid Release 1.0
       a2f2178f GenPipes - Covid Release : adding resources
       9e07fc4f GenPipes - HiC-Seq pipeline : corrected fastq_readName_edit input path
       dd582121 GenPipes - DNA-Seq pipeline : correcting symlink creation in sambamba_merge_sam_file
       7e00f1d3 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1799cea7 updated version of MUGQIC_TOOLS in installation script
       fb8ad6e0 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       7a9b872e GenPipes bug correction - corrected ln() function in bfx/bash_cmd.py
       62402add GenPipes Update - debugging use portal_output_dir variable : check for both undef and empty value
       2830d235 GenPipes Genome - added Sus_scrofa.sh (Pig genome) installation script
       81d580ad GenPipes Soft - added kent.sh installation script
       756d9a0e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       076e2008 GenPipes Update - DNASeq - updated sym_link for better handling of path reconstruction
       19cc2c06 GenPipes Update : fixing path of config files when passed to job2json script
       1d204f91 DNASeq - removed use of 'os.path.abspath' in call of 'ln()'
       399aac0b DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable

  Édouard Henrion <henrione@beluga3.int.ets1.calculquebec.ca>      4 commits

       46405691 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       8c465415 GenPipes - ChIPSeq update : a bit of code cleaning and simplifying
       b03ee2f6 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       72c75381 GENOME INSTALLATION - updated install genome.sh : added bismark genome preparation + refined genome digest command

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      2 commits

       0dc068bb Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       c7449389 Merged in eh_quick_fixes (pull request #144)

  Édouard Henrion <henrione@ip18.m>      6 commits

       b1909265 GenPIpes Update - corrected one problematic sym_link call...
       ed444a74 GenPipes Update - corrected pipeline behavior regarding PORTAL_OUTPUT_DIR environment variable : if te variable is empty or not set, then no JSON file at all will be generated
       c5407517 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6331eed2 C3G Software - added demuxlet.sh installation script
       d63b3568 Genome Update - install_genome.sh
       a899066e Genome Update - some updates to Homo_sapiens.GRCh38.sh

  ehenrion <edouard.henrion@mcgill.ca>      15 commits

       7d9f8114 BFX- Software - iVar installation script update with latest version
       3c3a80c2 common.py : genpipes replacing mugqic_pipeline....
       8658f7b5 GenPipes RNA-SEq - calling DESeq2 instead of DESeq for differential_expression
       6a0c8255 Coorected typo in README.md
       d59cb594 GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       af972e15 GenPipes - pipelines/dnaseq.py  : corrected prefix generation in SymLinkFastq step
       8ebbd14e GenPipes - pipelines/common.py corrected outputs name generation patterns for SamToFastq & Trimmomatic steps
       65fe80a3 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918129493 (pull request #142)
       97afb33b GenPipes - dnaseq.py : bug correction - typo removed
       4f9a8dd7 GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       fc67a7bc GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       9fbd085d GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       16a5635b GenPipes - nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       853c806d updated methylseq.base.ini, useless comments removed
       be965bdc updated methylseq.base.ini, useless comments removed

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      27 commits

       d4857288 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       7cb40243 final edits to the nanopore pipeline
       f1db8a94 Added support for gzipped fastq
       8da4a1c5 Added the nanopore CIT ini file
       86988dfd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       71f1658c Final corrections before merge to dev
       6994d693 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       83a6f5f7 More corrections to the INI files for other servers to allow proper running with SLURM
       caf6cb2e Corrected an error in the mp2b ini file
       e2637c0d Made corrections to the nanopore ini files for all servers that were causing the pipeline to break when running with SLURM
       a10fe05c Merge conflict resolutions
       e08cd6d7 Corrected problem with nanopore readset parsing that caused problem with the paths
       5d884c13 Added full documentation for Nanopore pipeline, as well as the Graham config file.
       e2edb704 Included commands necessary to add readgroup tag to alignments in minimap2
       278feb74 Final corrections before testing on other servers
       cc0f765a Corrected merge error related to .gitignore
       cb837885 Fixed bug caused by missing module import in nanopore.py
       649e04b2 Added minimap2 script that was missing from previous commit
       0f2a89dd First working version of the nanopore pipeline
       7acad993 Added full documentation for Nanopore pipeline, as well as the Graham config file.
       7c528bd8 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into nanopore_jhg
       56699ff4 Included commands necessary to add readgroup tag to alignments in minimap2
       3d53188d Final corrections before testing on other servers
       558f3e9b More bug corrections for nanopore pipeline after initial testing. Switched to only one protocol, with an optional first step (guppy)
       5107ed7c Fixed bug caused by missing module import in nanopore.py
       653c4e08 Added minimap2 script that was missing from previous commit
       a0fe5ebd First working version of the nanopore pipeline

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      8 commits

       50bc9f3a Merged in Jose-Hector-Galvez/rnaseq_lightbaseini-edited-online-with-b-1580485656011 (pull request #171)
       e2c91c8b Added module_perl to the rnaseq_light ini file.
       007835bd Merged in Jose-Hector-Galvez/jobpy-edited-online-with-bitbucket-1576266027164 (pull request #154)
       c500f6b4 Suggestion, add `module purge` to all jobs that load modules, to avoid conflicts between modules.
       4c35a7c8 Merged in Jose-Hector-Galvez/gatk4py-edited-online-with-bitbucket-1575403540009 (pull request #135)
       89e70fe0 gatk4.py edited to correct for inconsistencies with configuration parameters within functions.
       eefc95aa Merged in Jose-Hector-Galvez/found-a-bug-in-the-schedulerpy-script-i--1575322351581 (pull request #127)
       afae6eec Found a bug in the scheduler.py script. I am adding a line to correct it.

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      1 commits

       3ec3f1c7 Added CoVSeQ ini files for Graham and Cedar. Corrected a few errors on the Beluga ini

  José Héctor Gálvez López <jose.hector.galvez@computationalgenomics.ca>      2 commits

       10e762bb Merged in mgi_stretenp (pull request #177)
       39607170 Merged in nanopore_jhg (pull request #173)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      644 commits

       41514c00 Merged in mgi_stretenp (pull request #175)
       d6ae1290 Adding WARN for not changing R ver in rnaseq denovo
       f51f6938 Fixing rnaseq denovo R issue by creating deseq and a deseq2
       4a1a5bff Fixing back rnaseq denovo
       56fb3fbf Fixing rnaseq denovo
       1d29eb6d Fixing rnaseq denovo assembly R versions
       50fb4704 Fixing trinity report on cit
       1ae8b177 deliverables fix
       c51f7e18 Switching to 10% as minor variants threshold
       dca8feac Removing kraken module & Changing R version for rnaseq
       c49a0733 Including sambamba changes in ini
       29b5aec3 Including sambamba changes into other ini
       c4a7ef31 Fixing gatk_indel_realigner if 1 job
       8fd747f1 Including sambamba modifs in ini
       d826d819 sambamba merge realigned high coverage
       f135bafc Switching to sambamba merge sam
       c55297c3 Including bwa sambamba into high coverage
       466906d1 Inluding a with sambamba bam.bai to picard mark duplicates
       858c1f53 Merge branch 'dnaseq_picard_fix' into mgi_stretenp
       a3354a52 Typo
       dd39dedd Fixing for cit run
       1b4d0149 Merge branch 'dev' into mgi_stretenp
       d9566ff3 Using Illumina as default
       34cd4223 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       2502abac Renaming outputs and default genome
       8ccb242e Fixing consensus renaming
       075ab5a6 Fixing consensus renaming
       9dec6685 Renaming + flaging consensus seq
       621449ab Fixing quast options
       da263847 Fixing qualimap output
       cb98c5ef Fixing outputs quast and snpeff
       6d220605 Fixing mkdir
       d15425c8 Fixing typo
       921d5657 Fixing typo
       53c68fdd Changing to samtools bam2fq and fixing options ini
       4f5c5ceb Fixing quast
       f7042ec7 Fixing latest commit
       04f29165 Adding intermediate bam file on host removal step
       d931a893 Forcing pigz to overwritte gz
       b0711a8e fixing host removal
       4cfa142e Adding pigz within bam2fq to be able to skip step
       7916935b quast fix
       4788baed quast + snpeff + host removal fix
       2e3bc597 Fixing cutadapt input files
       75d63826 Fixing type
       308ace13 Removing indexing after name sorting
       f76fdadc Fixing path creation at host removal step
       d5aa6f0e Fixing snpeff
       92024780 Removing print files
       afa134ef Fixing input choosing
       979eafe2 Fixing choosing inout files mapping
       9b6ae435 Fixing host removal
       f7ba7254 Fixing host removal
       f0ac000f Fixing quast step
       8ce9bf21 Fixing quast
       c07ffae7 Fixing quast step
       5415ab72 Fixing param requirements
       9218a8e7 Fixing typo
       52999aee Fixing pigz
       fec2ceae Not using kraken anymore
       3363c61c Switching steps order
       7e51be51 Adding 3 steps
       863d61b1 Fixing picard multiple metrics raw
       4d8cc43b Insert Size metrics on raw bam
       f2e8b42b sambamba flagstat fix
       0a8c640f Renaming snpeff step/job
       d4734044 Fixing flagstat
       8163c9a3 Fixing flagstat
       8461e198 flagstat on all the bams
       a380f681 Zipping output of snpeff
       28969ba9 Fixing snpeff
       a899bdd5 Fixing cutadapt
       2f86db55 Flagstat on raw bam
       f2b8a9f0 Fix
       f4fba5a7 Fixing renaming step
       b0fcd1cb Switching to fgbio trim primer
       e9c450f8 Updating metrics
       6582bc9b Addition of consensus step
       319e50d6 MGI init commit
       6df0cca4 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       350131aa DNA-Seq - Fix update
       5cce8187 DNA-Seq - Fix update
       838cc920 DNA-Seq - Fix update
       990fc564 DNA-Seq - Fix update
       7e7b8107 DNA-Seq - Fix update
       0a8d7eee DNA-Seq - Fix update
       fe22929a DNA-Seq - Fixmate with samtools and sorting with sambamba
       3c84dd1e DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       e0fa9ff8 DNA-Seq - Fix update including input_file_revision_eh changes
       d03669f1 DNA-Seq - Fix update
       5f22222d DNA-Seq - Fix update
       18dd1069 DNA-Seq - Fix update
       08c2864b DNA-Seq - Fix update
       117181c2 DNA-Seq - Fixmate with samtools and sorting with sambamba
       813c02e6 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c193068c DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       0e8d422c DNA-Seq - Fix update including input_file_revision_eh changes
       80f732a4 Zipping and indexing vcfs
       6b1179bb Adding parameter to ivar consensus
       a463f5b3 Changing caller
       75a071ac Addition of consensus step
       56399851 MGI init commit
       98118422 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       34dac0c6 cutadapt bfx
       9ccadd39 DNA-Seq - Fix update
       8733abc8 DNA-Seq - Fix update
       e007f9a4 DNA-Seq - Fix update
       93e342ce DNA-Seq - Fix update
       0160ab53 DNA-Seq - Fix update
       567bd8a9 DNA-Seq - Fixmate with samtools and sorting with sambamba
       9e4119b0 DNA-Seq - Fix update
       7ed58f89 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       aefdd2c6 DNA-Seq - Fix update including input_file_revision_eh changes
       8a94a76e DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       36fbd806 DNA-Seq - Fix update
       cb43b486 DNA-Seq - Fix update
       069de1c2 DNA-Seq - Fix update
       2775e6f9 DNA-Seq - Fix update
       348c2c61 DNA-Seq - Fix update
       809b43dd DNA-Seq - Fix update
       08c41830 DNA-Seq - Fix update
       7314301c DNA-Seq - Fixmate with samtools and sorting with sambamba
       b4728a85 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       d22a3b11 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7a608bef DNA-Seq - Fix update including input_file_revision_eh changes
       96dcc43a MGI-Seq - First commit
       c6e80af2 Adding snpEff step for mgi pipeline
       0e998c0a Fixing sambamba flagstat
       79e81149 Fixing sambamba module
       3bd9321b Adding module sambamba to ini
       812e317b Fixing sambamba indexing
       b0cbd931 Switching from picard to sambamba withing methylseq
       3f5d5a3b Fix
       271b81d5 fix
       6b0c09fa Fixing input output files
       79e49f4b Fixing filtering bam
       c4dee0fb Fixing renaming step
       f1311248 Fix renaming step
       c765138d fix
       c04a496e Fix
       eea7c7ba Fixing filtering
       c0bc132f Adding filtering step
       7d34c4e4 Switch to fgbio
       f38cd7d9 ivar triming switch
       15b7438f Switching to fgbio trim primer
       cec53cbc Updating metrics
       b8de4e31 MGI init commit
       11b01b0d cutadapt bfx
       e3858aa5 DNA-Seq - Fix update
       1eaa3dc1 DNA-Seq - Fix update
       d23743fc DNA-Seq - Fix update
       bd59ddbb DNA-Seq - Fix update
       4eb32689 DNA-Seq - Fix update
       9d8e58ba DNA-Seq - Fix update
       c93e2339 DNA-Seq - Fix update
       2934b708 DNA-Seq - Fixmate with samtools and sorting with sambamba
       ff53bb30 DNA-Seq - Fix update
       5d3fb26d DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       bb731b11 DNA-Seq - Fix update
       1b681322 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       3a6e72dc DNA-Seq - Fix update including input_file_revision_eh changes
       01551e52 DNA-Seq - Fix update
       dcaedc85 DNA-Seq - Fix update
       a01251ad DNA-Seq - Fix update
       56a50d0c DNA-Seq - Fix update
       32e070f2 DNA-Seq - Fix update
       e4778b07 DNA-Seq - Fix update
       94560456 DNA-Seq - Fix update
       083f6a69 DNA-Seq - Fixmate with samtools and sorting with sambamba
       67603dda DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       a59d9940 DNA-Seq - Fix update including input_file_revision_eh changes
       c2039121 MGI-Seq - First commit
       c5a1baa5 Switch to fgbio
       b4ca5f77 Fixing alignment
       f3378e14 ivar triming switch
       8a0045dc Switching to fgbio trim primer
       d479b856 Fixing ivar trim
       80c578a2 Filtering reads
       30f51c7d Updating metrics
       ea1149b2 fgbio
       22daf31f Fixing ivar trim
       64629082 Using ivar trim instead of fgbio
       fc6de5f6 Adding ivar primer trimming
       7049f311 Default bwa parameters to include pairs
       0e834bc9 Zipping and indexing vcfs
       0935479d Adding parameter to ivar consensus
       d8751ccc Changing caller
       1e6dd749 ivar bfx
       3399cbb1 Addition of consensus step
       f46ccb5c ivar module
       52f8a2b1 MGI init commit
       3bd09324 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       2fbd1e34 Adding trim_primers to fgbio bfx
       327b963f cutadapt bfx
       7652e5b8 DNA-Seq - Fix update
       0ab9060c DNA-Seq - Fix update
       66b240b8 DNA-Seq - Fix update
       472c6b0e DNA-Seq - Fix update
       e1e7318a DNA-Seq - Fix update
       6954eaf7 DNA-Seq - Fix update
       db7bab26 DNA-Seq - Fix update
       2a0bb004 DNA-Seq - Fix update
       60e7d9a8 DNA-Seq - Fixmate with samtools and sorting with sambamba
       7be45c10 DNA-Seq - Fix update
       1dddc4f4 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       27ac18b6 DNA-Seq - Fix update including input_file_revision_eh changes
       748e668a DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       0b537860 DNA-Seq - Fix update
       6218cdd2 DNA-Seq - Fix update
       56547aa0 DNA-Seq - Fix update
       05baa00f DNA-Seq - Fix update
       ff5a8712 DNA-Seq - Fix update
       591339e7 DNA-Seq - Fix update
       de100ee1 DNA-Seq - Fix update
       5e3e598c DNA-Seq - Fix update
       584181f9 DNA-Seq - Fix update
       02439f06 DNA-Seq - Fixmate with samtools and sorting with sambamba
       90e9ee94 DNA-Seq - Fix update
       ac566918 DNA-Seq - Fix update
       f73bce0c DNA-Seq - Fix update
       22e030d5 DNA-Seq - Fix update
       69289267 DNA-Seq - Fix update
       0facc342 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       c6bf187f DNA-Seq - Fix update
       111cbcb5 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ef0a25e2 DNA-Seq - Fix update ini files
       2eed96b3 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       417ff67a DNA-Seq - Fix update bai output
       114d004a DNA-Seq - Fix update indexing
       fc9901c1 DNA-Seq - Fix update including input_file_revision_eh changes
       0016ba2b DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       77e50414 MGI-Seq - First commit
       9447bed7 Renaming outputs and default genome
       6f456559 Fixing consensus renaming
       87a62648 Fixing consensus renaming
       cb955341 Renaming + flaging consensus seq
       6c7058d3 Fixing quast options
       1ec14ff9 Fixing qualimap output
       90eb6b0d Fixing outputs quast and snpeff
       3dbf0a1d Fixing mkdir
       8cd2ba34 Fixing typo
       200cb4bd Fixing typo
       2b01f81d Changing to samtools bam2fq and fixing options ini
       d5dc5b77 Fixing quast
       f32dd814 Fixing latest commit
       d58a7b23 Adding intermediate bam file on host removal step
       dbd486b3 Forcing pigz to overwritte gz
       95505ba2 fixing host removal
       a1191009 Adding pigz within bam2fq to be able to skip step
       de7633c4 quast fix
       d0383f5f quast + snpeff + host removal fix
       f165e567 Fixing cutadapt input files
       2dfd2e62 Fixing type
       6e67b6da Removing indexing after name sorting
       4b0740ea Fixing path creation at host removal step
       7020fe0d Fixing snpeff
       5e6b7a91 Removing print files
       031d9dfa Fixing input choosing
       256cdc79 Fixing choosing inout files mapping
       cd7e52e9 Fixing host removal
       c11d2200 Fixing host removal
       a3d58c9c Fixing quast step
       8abe6c08 Fixing quast
       5a68e072 Fixing quast step
       b141fa8f Fixing param requirements
       a7318fdc Fixing typo
       aa32b62e Fixing pigz
       617f817b Not using kraken anymore
       df8f32ce Switching steps order
       10c3586d Adding 3 steps
       eb7d6120 Fixing picard multiple metrics raw
       8d852a76 Insert Size metrics on raw bam
       67dcc5f4 sambamba flagstat fix
       08fb9264 Renaming snpeff step/job
       068d072a Fixing flagstat
       770c233b Fixing flagstat
       4ff0fb6e flagstat on all the bams
       1e3ef28e Zipping output of snpeff
       5c219e64 Fixing snpeff
       c0faa375 Fixing cutadapt
       effce8a8 Flagstat on raw bam
       8c7dd5a0 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       221dea23 Adding snpEff step for mgi pipeline
       86a0acf3 Fixing sambamba flagstat
       dd8163f7 Fixing sambamba module
       609e6ca8 Adding module sambamba to ini
       71316e97 Fixing sambamba indexing
       90407220 Switching from picard to sambamba withing methylseq
       9d3bdfd1 Fix
       e5a7f40d fix
       59aeaff8 Fixing input output files
       da3b1cd0 Fixing filtering bam
       cb775da3 Fixing renaming step
       f78909cb Fix renaming step
       49dba7db fix
       50ac6ee7 Fix
       860591e5 Fixing filtering
       cc0f319c Adding filtering step
       4d36888f Switch to fgbio
       6066f526 ivar triming switch
       5790ea2e Switching to fgbio trim primer
       45c24491 Updating metrics
       ffc8fd16 MGI init commit
       13d4f105 cutadapt bfx
       2a61965b DNA-Seq - Fix update
       7c8d9e34 DNA-Seq - Fix update
       aad95a26 DNA-Seq - Fix update
       db01cdb0 DNA-Seq - Fix update
       c8fd1f4a DNA-Seq - Fix update
       19b7a74c DNA-Seq - Fix update
       386a6048 DNA-Seq - Fix update
       6063ddda DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e1f0b54 DNA-Seq - Fix update
       30121a77 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       4b9d8e1c DNA-Seq - Fix update
       a83ee80a DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       033dd3ec DNA-Seq - Fix update including input_file_revision_eh changes
       a9102d38 DNA-Seq - Fix update
       a4ee47fe DNA-Seq - Fix update
       18d832e5 DNA-Seq - Fix update
       6623f8c3 DNA-Seq - Fix update
       1d002d1f DNA-Seq - Fix update
       fe5c4057 DNA-Seq - Fix update
       dbd8034b DNA-Seq - Fix update
       74d2bdde DNA-Seq - Fixmate with samtools and sorting with sambamba
       16786398 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       5f66227b DNA-Seq - Fix update including input_file_revision_eh changes
       562180b9 MGI-Seq - First commit
       a16695b5 Switch to fgbio
       1d92909d Fixing alignment
       01f5df72 ivar triming switch
       e49e632b Switching to fgbio trim primer
       057e2e8f Fixing ivar trim
       d9638b55 Filtering reads
       c3ab6f14 Updating metrics
       dc9a43ee fgbio
       579e0dd4 Fixing ivar trim
       603bb4a0 Using ivar trim instead of fgbio
       26a2da86 Adding ivar primer trimming
       f25c75d2 Default bwa parameters to include pairs
       6a8c3018 Zipping and indexing vcfs
       cf6c19f1 Adding parameter to ivar consensus
       c9231e5c Changing caller
       b9bfed08 ivar bfx
       70cc54a3 Addition of consensus step
       5a990e2b ivar module
       cd17c0d0 MGI init commit
       5bef158f Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       f7a155d2 Adding trim_primers to fgbio bfx
       8177037a cutadapt bfx
       a702d277 DNA-Seq - Fix update
       f6fe0e2d DNA-Seq - Fix update
       4a424397 DNA-Seq - Fix update
       a880e3bf DNA-Seq - Fix update
       2cd5d012 DNA-Seq - Fix update
       9f7e0480 DNA-Seq - Fix update
       75b18caf DNA-Seq - Fix update
       9a0adea9 DNA-Seq - Fix update
       bb259fe4 DNA-Seq - Fixmate with samtools and sorting with sambamba
       0c969370 DNA-Seq - Fix update
       75e589b3 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       dd3999c1 DNA-Seq - Fix update including input_file_revision_eh changes
       2a2a43c8 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       cb09a7b3 DNA-Seq - Fix update
       31653bd6 DNA-Seq - Fix update
       513818ca DNA-Seq - Fix update
       b202a5b6 DNA-Seq - Fix update
       c87c2cd2 DNA-Seq - Fix update
       a4926cbe DNA-Seq - Fix update
       c06ee52b DNA-Seq - Fix update
       6f894d2c DNA-Seq - Fix update
       3d5daeb1 DNA-Seq - Fix update
       ed9e7feb DNA-Seq - Fixmate with samtools and sorting with sambamba
       4215929b DNA-Seq - Fix update
       3dca79c7 DNA-Seq - Fix update
       7bf5141a DNA-Seq - Fix update
       ff24ef91 DNA-Seq - Fix update
       f5352026 DNA-Seq - Fix update
       a6b428c0 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       7e351683 DNA-Seq - Fix update
       302e6c9a DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       03fc8c29 DNA-Seq - Fix update ini files
       0cd48857 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       ade82794 DNA-Seq - Fix update bai output
       8ee32671 DNA-Seq - Fix update indexing
       4f25a984 DNA-Seq - Fix update including input_file_revision_eh changes
       314107e4 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       6d82d3f5 MGI-Seq - First commit
       2cbfa5c0 Adding snpEff step for mgi pipeline
       8f5d8394 Fixing sambamba flagstat
       4415a7ae Fixing sambamba module
       644bc8d0 Adding module sambamba to ini
       38cd45ad Fixing sambamba indexing
       08dc201c Switching from picard to sambamba withing methylseq
       e89813ff Fix
       9ef7ecbd fix
       0d6651a8 Fixing input output files
       192bb652 Fixing filtering bam
       030b759f Fixing renaming step
       1348f004 Fix renaming step
       bc42658a fix
       49a16fc7 Fix
       c69973dc Fixing filtering
       76219a86 Adding filtering step
       7ed54e8b Switch to fgbio
       b6cfddfc Fixing alignment
       4796cffb ivar triming switch
       3ec5f7f3 Switching to fgbio trim primer
       3a84e9eb Fixing ivar trim
       b673f473 Filtering reads
       9fdafbd9 Updating metrics
       3f63d5ca fgbio
       56bec970 Fixing ivar trim
       30475caa Using ivar trim instead of fgbio
       d5b40548 Adding ivar primer trimming
       925645e1 Default bwa parameters to include pairs
       314a8695 Zipping and indexing vcfs
       191bfb3d Adding parameter to ivar consensus
       db7c70a3 Changing caller
       d9e302b1 ivar bfx
       09d9f6b5 Addition of consensus step
       d9c3f32f ivar module
       f6a8e172 MGI init commit
       49c7f1d0 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       eae7037e Adding trim_primers to fgbio bfx
       bfc47be8 cutadapt bfx
       ffcfe262 DNA-Seq - Fix update
       0f174265 DNA-Seq - Fix update
       34b09658 DNA-Seq - Fix update
       b7323de8 DNA-Seq - Fix update
       e0909f14 DNA-Seq - Fix update
       48e203ed DNA-Seq - Fix update
       0cd7e5e6 DNA-Seq - Fix update
       b3d74010 DNA-Seq - Fix update
       6ff06b00 DNA-Seq - Fixmate with samtools and sorting with sambamba
       b000f97a DNA-Seq - Fix update
       52deaf39 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       2acf3c0b DNA-Seq - Fix update including input_file_revision_eh changes
       0fe44df6 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       3772b9bc DNA-Seq - Fix update
       479e9ed4 DNA-Seq - Fix update
       5f251207 DNA-Seq - Fix update
       aa1687cb DNA-Seq - Fix update
       8632426c DNA-Seq - Fix update
       c5d384b8 DNA-Seq - Fix update
       746d2f01 DNA-Seq - Fix update
       eef3b276 DNA-Seq - Fix update
       f1d1b551 DNA-Seq - Fix update
       d2d4430e DNA-Seq - Fixmate with samtools and sorting with sambamba
       56832972 DNA-Seq - Fix update
       a7e450cc DNA-Seq - Fix update
       5cb95740 DNA-Seq - Fix update
       e888e098 DNA-Seq - Fix update
       d366f675 DNA-Seq - Fix update
       ba6aa73c DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       9dc3bf6f DNA-Seq - Fix update
       c01dcf52 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       ca81d44f DNA-Seq - Fix update ini files
       78519f29 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       0c922ccd DNA-Seq - Fix update bai output
       6117bc98 DNA-Seq - Fix update indexing
       5e9b2685 DNA-Seq - Fix update including input_file_revision_eh changes
       519beea3 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       43bb771a MGI-Seq - First commit
       f2cd11aa Switch to fgbio
       c3aa00dd Fixing alignment
       8d157fd4 ivar triming switch
       42a5e2d7 Switching to fgbio trim primer
       6030fa73 Fixing ivar trim
       6bb7e087 Filtering reads
       330e345c Updating metrics
       b1153cd5 fgbio
       4450c70e Fixing ivar trim
       a68648b0 Using ivar trim instead of fgbio
       8a2bdf36 Adding ivar primer trimming
       8856367b Default bwa parameters to include pairs
       4e968647 Zipping and indexing vcfs
       6e9529f8 Adding parameter to ivar consensus
       4b826ef9 Changing caller
       b8bc79a0 ivar bfx
       dab5579f Addition of consensus step
       c91d0676 ivar module
       ec2a1760 MGI init commit
       260059a7 Adding bed2interval_list to gatk, gatk4 and picard, picard2 bfx
       7760a795 Adding trim_primers to fgbio bfx
       3180972e Fixing haplotype_caller
       66d8cdda Fixing haplotype_caller
       3bdb0aa9 cutadapt bfx
       2d274eb2 Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into mgi_stretenp
       f083b668 MGI-Seq - First commit
       538907e9 Merge branch 'dnaseq_picard_fix' of bitbucket.org:mugqic/genpipes into dnaseq_picard_fix
       b0197fe1 DNA-Seq - Fix update
       54b4741c DNA-Seq - Fix update
       452f8692 DNA-Seq - Fix update
       d0cb631d DNA-Seq - Fix update
       6e3580e6 DNA-Seq - Fix update
       222b79b0 DNA-Seq - Fix update
       5cf62699 DNA-Seq - Fix update
       de8d9017 DNA-Seq - Fix update
       d87cd96c DNA-Seq - Fix update
       b3b9b1b6 DNA-Seq - Fixmate with samtools and sorting with sambamba
       5e55bfc9 DNA-Seq - Fix update
       9f20fe78 DNA-Seq - Fix update
       123b7303 DNA-Seq - Fix update
       5ae600ac DNA-Seq - Fix update
       8119523b DNA-Seq - Fix update
       d6007d20 DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       1ff8cacc DNA-Seq - Fix update
       58d662a7 DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       83ba60d0 DNA-Seq - Fix update ini files
       f7bef7df DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       10d50aa8 DNA-Seq - Fix update bai output
       1319c26b DNA-Seq - Fix update indexing
       fe4e7b99 DNA-Seq - Fix update including input_file_revision_eh changes
       80def2ad DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       894566cc DNA-Seq - Fix update
       1413d451 DNA-Seq - Fix update
       ab945dd0 DNA-Seq - Fix update
       97186dd2 DNA-Seq - Fix update
       fc158a47 DNA-Seq - Fix update
       86b2e803 DNA-Seq - Fix update
       d2ecd4b9 DNA-Seq - Fix update
       98928a55 DNA-Seq - Fix update
       bdee1f80 DNA-Seq - Fix update
       8afb5005 DNA-Seq - Fixmate with samtools and sorting with sambamba
       7e8699ba DNA-Seq - Fix update
       1e47ab11 DNA-Seq - Fix update
       bea5e948 DNA-Seq - Fix update
       eae982e8 DNA-Seq - Fix update
       15d745e2 DNA-Seq - Fix update
       0553a5de DNA-Seq - Fix update including input_file_revision_eh changes on samtools
       ef79ced1 DNA-Seq - Fix update
       afaa484b DNA-Seq - Fix update including input_file_revision_eh changes on gatk
       7dfa7303 DNA-Seq - Fix update ini files
       3d22e672 DNA-Seq - Fix update including input_file_revision_eh changes on job.py
       6dd6e0c6 DNA-Seq - Fix update bai output
       54a5fc05 DNA-Seq - Fix update indexing
       20cee1f0 DNA-Seq - Fix update including input_file_revision_eh changes
       a1579f46 DNA-Seq - First commit for https://bitbucket.org/mugqic/genpipes/issues/21/index-file-of-bam-bai
       c5e99661 Merged in ihec_metrics (pull request #126)
       7bceefb1 Merge branch 'ihec_metrics' of bitbucket.org:mugqic/genpipes into ihec_metrics
       388c9bf4 ChIP-Seq - Fix IHEC metrics
       b5683f25 ChIP-Seq - Fix IHEC metrics
       7c0f68fb Merged in chipseq_atacseq_mode (pull request #121)
       7a101a70 Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       e83caef0 MethylSeq - Trimmomatic resources revisited for methylseq
       f9bbdb79 General - Genome installation with cvmfs grep
       a8c4d934 ChIP-Seq - Fixing GenPipes metrics
       d82ed899 ChIP-Seq - Sambamba loading properly in mito calcul
       e6c4274b ChIP-Seq - Reducing ihec report ressources
       f3755e93 ChIP-Seq - Fixing typo
       31530179 ChIP-Seq - Increasing IHEC preprocess ressources
       a02fd322 ChIP-Seq - Fix IHEC metrics
       7f1f4a96 ChIP-Seq - Fix IHEC report template md file
       805c58b6 ChIP-Seq - Fix IHEC report template md file
       19b0b629 ChIP-Seq - Fix IHEC report template md file
       3c573f8d ChIP-Seq - Fixing merge metrics
       1bb96cb2 ChIP-Seq - Fixing merge metrics
       97143130 ChIP-Seq - Fixing merge metrics
       ad5e5313 ChIP-Seq - Fixing merge metrics
       f3ebb74c ChIP-Seq - Fixing merge metrics
       b5658b8b ChIP-Seq - Adding IHEC report template md file
       3ad67e06 ChIP-Seq - Fixing merge metrics
       bef0f187 ChIP-Seq - Fixing merge metrics
       ec3dea00 ChIP-Seq - Fixing merge metrics
       4883d4c3 ChIP-Seq - Fixing merge metrics
       51bd5f2c ChIP-Seq - Adding merge IHEC metrics
       c310d1e2 ChIP-Seq - Fixing metrics
       b1ec49d4 ChIP-Seq - Fixing metrics
       50841713 ChIP-Seq - Fixing metrics
       fa035a10 ChIP-Seq - Fixing metrics
       85c8c53b ChIP-Seq - Fixing metrics
       41081d62 ChIP-Seq - Fixing metrics
       1acdfa5f ChIP-Seq - Fixing metrics
       cfe40591 ChIP-Seq - Fixing metrics & report
       d29881f1 ChIP-Seq - Fixing metrics
       d1d43852 ChIP-Seq - Fixing metrics
       5e801f94 ChIP-Seq - Fixing metrics
       a1a47842 ChIP-Seq - Fixing metrics
       0f847f18 ChIP-Seq - Fixing metrics
       2f9084af ChIP-Seq - Fixing metrics
       e545fb48 ChIP-Seq - Fixing metrics
       56f6af0d ChIP-Seq - Fixing metrics
       878335c1 ChIP-Seq - Fixing metrics
       380e1224 ChIP-Seq - Fixing metrics
       f6f26f4b ChIP-Seq - Fixing metrics
       343c9ffa ChIP-Seq - Fixing metrics
       3da997d3 ChIP-Seq - Fixing metrics
       db54da14 ChIP-Seq - Fixing metrics
       6a91f0d3 ChIP-Seq - Adding metrics
       b6a1e608 ChIP-Seq - Fixing bwa missing import
       319655b7 ChIP-Seq - Fixing chipseq pipeline
       1a47e31b ChIP-Seq - Adding ATAC-Seq protocol
       5f113fce Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       a440d4e6 General - Genome installation with cvmfs grep
       d7cd595e ChIP-Seq - Fixing GenPipes metrics
       bbfb6bb4 ChIP-Seq - Sambamba loading properly in mito calcul
       466041cb Merge branch 'dev' into chipseq_atacseq_mode
       15611bfa ChIP-Seq - Reducing ihec report ressources
       56c39f6f ChIP-Seq - Fixing typo
       f88ee593 ChIP-Seq - Increasing IHEC preprocess ressources
       0bd494dd Merge branch 'chipseq_atacseq_mode' of bitbucket.org:mugqic/genpipes into chipseq_atacseq_mode
       1d727d08 ChIP-Seq - Fix IHEC metrics
       09d7280f Merged dev into chipseq_atacseq_mode
       364e3811 ChIP-Seq - Fix IHEC report template md file
       24b64ed8 ChIP-Seq - Fix IHEC report template md file
       9f8c10e7 ChIP-Seq - Fix IHEC report template md file
       0d811bd4 ChIP-Seq - Fixing merge metrics
       70a5ed0d ChIP-Seq - Fixing merge metrics
       b1da6004 ChIP-Seq - Fixing merge metrics
       aee5b149 ChIP-Seq - Fixing merge metrics
       74f44e23 ChIP-Seq - Fixing merge metrics
       b096b635 ChIP-Seq - Adding IHEC report template md file
       12b2f31f ChIP-Seq - Fixing merge metrics
       a807a5e3 ChIP-Seq - Fixing merge metrics
       5131a323 ChIP-Seq - Fixing merge metrics
       60b0f6fc ChIP-Seq - Fixing merge metrics
       05c26231 ChIP-Seq - Adding merge IHEC metrics
       9664631a ChIP-Seq - Fixing metrics
       1e2bc38a ChIP-Seq - Fixing metrics
       c6a08a6f ChIP-Seq - Fixing metrics
       e9e96da2 ChIP-Seq - Fixing metrics
       809cb52b ChIP-Seq - Fixing metrics
       78736476 ChIP-Seq - Fixing metrics
       4b1bb18e ChIP-Seq - Fixing metrics
       e806abc6 ChIP-Seq - Fixing metrics & report
       ca01c9cd ChIP-Seq - Fixing metrics
       3c72d0f9 ChIP-Seq - Fixing metrics
       1199e5f9 ChIP-Seq - Fixing metrics
       1ee3d364 ChIP-Seq - Fixing metrics
       59a28284 ChIP-Seq - Fixing metrics
       51570690 ChIP-Seq - Fixing metrics
       4a8a555f ChIP-Seq - Fixing metrics
       2cb7a7be ChIP-Seq - Fixing metrics
       a48c88ee ChIP-Seq - Fixing metrics
       235552a2 ChIP-Seq - Fixing metrics
       5b47faf9 ChIP-Seq - Fixing metrics
       31535af0 ChIP-Seq - Fixing metrics
       2d751ea5 ChIP-Seq - Fixing metrics
       85ea7e4e ChIP-Seq - Fixing metrics
       b62f5193 ChIP-Seq - Adding metrics
       bbe8b2c2 ChIP-Seq - Fixing bwa missing import
       5bdb8aa2 ChIP-Seq - Fixing chipseq pipeline
       2ae58031 ChIP-Seq - Adding ATAC-Seq protocol

  Paul Stretenowich <pstretenowich@CYPRUS.local>      4 commits

       20aa88fd Cleaning mgi.py
       cc25234c Cleaning mgi.py
       397efae3 Cleaning mgi.py
       97bd7765 Merge branch 'mgi_stretenp' of bitbucket.org:mugqic/genpipes into mgi_stretenp

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      10 commits

       31826ded Merged in remove_pacbio (pull request #176)
       68d5195a Merged in poq/fast_module_check (pull request #164)
       652e207d Merged poq/graham_ini into dev
       9c6da9a8 Merged poq/graham_ini into dev
       1f66c6c9 Merged update_beluga_ini into dev
       70167f13 Merged update_beluga_ini into dev
       e7c7f349 Merged update_beluga_ini into dev
       a54d824c Merged in add_container (pull request #167)
       8a076c75 Merged in poq/graham_ini (pull request #168)
       f18307f4 Merged in master (pull request #161)

  P-O Quirion <pierre-olivier.quirion@computationalgenomics.ca>      33 commits

       241f706f remove pacbio from the repo/release
       fd02e500 no mail in cit.ini
       68d8a669 add sh script for steps in pbs too
       d79b2736 Add SARS-CoV2 genome file
       b4c47402 tweek memory usage beluga denovo
       7fcd09ae update cluster for rnaseq star index
       473294db use 1.1.0 genpipes_in_container release
       0d591bba make module show sure it raise with older version
       1aa57cbd chipseq cedar and graham ini
       1329a13c Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       6f7be7f2 more graham ini
       b455c4cf more graham ini
       82676033 copy dnaseq ini from cedar to graham
       b4485973 fix fastq symlink on graham
       153aa9e2 less log
       40ab0ff1 cleanup
       0ac31218 speedup module check
       20ff188d one module call
       1fb3ae75 add mem to star_* on rnaseq
       4d682b61 update again
       19a81cf1 mem-per-cpu
       e5ce4202 Add generic ini for chipseq update README
       d1037610 create graham ini file
       05ef09e5 log when all sbatch submit is dont
       c1f2c531 feedback at submit time
       0ac30b35 feedback at submit time
       fd755d81 fix ini typo
       77b2c9b7 wrapper for slum, pbs and batch
       bc80a4d8 update wrap options.
       304c7eb5 add default wrapper
       22d9efac remove docker
       004a9a43 put --wrap import at the top
       3516da40 add wrapper to all pipelines

  P-O Quirion <pioliqui@gmail.com>      2 commits

       328a00ef make nanopore executable
       ab1de9a9 Add automatic wrapper option

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       123f6c57 Merged in fix-watch-folder (pull request #165)

  Rom Grk <romgrk.cc@gmail.com>      1 commits

       1baeef5a watch_portal_folder.py: fix undefined variable

  ufgauthi <ulysse.fortiergauthier@mcgill.ca>      1 commits

       0a60ef29 Bug Fix by replacing sacct delimiter | by ^

  Ulysse Fortier Gauthier <ulysse.fortiergauthier@mcgill.ca>      1 commits

       bf54dbf9 Merged in ufg_log_report_fix (pull request #157)

3.1.5        Wed Jan 15 11:58:16 2020 -0500        424 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      78 commits

       e0844c30 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       986f2a3c HICSeq pipeline - improving trimmomatic resources for SULRM + added graham config file for hicseq pipeline
       6714a5d4 Genome installation - added America Mink genome (Neovison_vison.NNQGG.v01) installation script
       d60d1368 Software installation - added 'demuxlet' installation script
       6822cbfb Software installation - updated ucsc.sh with lat4est version i.e. v387
       b350d2da Software installation - replaced call of lsb_release by universal commands (avoid 'lsb_release command not found' error)
       5c87fd72 Software installation - corrected regular expression within genome_digest function
       23e93153 removed some test code introduced by latest commit...
       9445ee98 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dev
       00577d4f GenPipes RNA-Seq de novo - corrected bugs introduced in commit 4a72735
       f3b9a92c GenPipes BFX - added bash_cmd python wrapper to wrap basic bash commands
       4a727351 GenPipes RNA-Seq de novo - updated pipeline with better sample assignemnt to jobs for better JSON building
       d89d4fc3 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       1b482150 GenPipes JSON - updated pipeline.py and scheduler.py so that the copy of the original JSONs to the PORTAL is made before building the pipeline steps
       6360349e GenPipes Trimmomatic - corrected assignment of samples to the job to correct analysis JSON generation per sample
       fbaa84a9 MethylSeq - corrected assignment of samples to the jobs to correct analysis JSON generation per sample
       d1cd79b3 GenPipes analysis JSON - updated jsonator.py in reagrds of the updated .base.ini files, now containing 'source' & 'version' to ease the process + minor updates to core/pipeline.py
       fc1b6c0f GenPipes config files - updated most of the .base.ini files with missing 'source' & 'version' to avoid issues with the 'jsonator'
       f44ccb15 GenPipes Anlalysis JSON file - corrected core/pipeline.py regarding the use of PORTAL_OUTPUT_DIR
       9751513e GenPipes Analysis JSON file - added the system to update the submission_date
       a9ba8320 GenPipes Anlalysis JSON file - pipelines now create JSON analysis file as a default behavior
       7b42d155 GenPipes Analysis JSON file - added the submission_date to the JSON
       7e04404f Analysis JSON - updated jsonator.py with project_name and submission_date. Also updated version to 1.0.1
       f422649a GenPipes utils - minor updates : updated some comments
       38937e61 Software install - updated install_module.sh with finer LIBDIR for patching and better patching to avoid overwritting potential pre-existing RPATH
       fb2bffc0 Software install - updated R_Bioconductor.sh with new R packages and finer LIBDIR for patching
       721d9416 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       074045aa GenPipes Sanity Check - Refining the report
       0d0e0bdd Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       1478455a GenPipes Sanity Check - Adjusting the log level
       cf3ec9e3 GenPIpes Sanity Check - set log level to 'warning' instead of info when running sanity check, thus removing some useless if statements and refining display of the messages
       6cf36ee5 GenPipes Sanity Check mode : removed some useless comments in core/pipeline.py
       b46df3f8 GenPIpes Sanity Check mode - created the SanitycheckError class in config.py to use instead of Error
       f4d6323c Merge branch 'master' of bitbucket.org:mugqic/genpipes into rnaseq_denovo_jhg
       989996a1 GenPipes MethylSeq - updated base.ini with [samtools_cram_output] section
       066c24e3 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       433950d2 GenPipes Sanity Check - DNASeq high coverage, RNASeq denovo assembly & RNASeq light pipelines are updated regarding the sanity-check mode
       0fbecd2e Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       6b00ce5d GenPipes Sanity Check - updated pipelines for sanity-check mode responsiveness
       afa406ab GenPipes Sanity Check - updated common.py to be sanity-check responsive
       4f25a0fa GenPipes Sanity Check - update _raise() function in config.py & refined code of sanity-check mode in pipeline.py
       6a457d53 GenPipes Sanity Check - updated readset.py design.py to be sanity-check responsive & tools.py with proper import statements
       6f11825c MethylSeq pipeline - updated pipeline to make SINGLE-END mode actually work
       b3ac5a9e GenPipes Sanity Check mode : updated PacBio Assembly pipeline as the first try to test sanity check mode
       f2885449 GenPipes Sanity Check mode : updated pipeline.py with the sanity check mode fully functionnal
       b5a78dda GenPipes Sanity Check mode : updated sample.py to reflect the updates that have been done in config.py
       3203ce90 GenPipes Sanity Check mode : updated config.py to avoid raising errors when sanity-check mode is on, logging them instead
       daabbf45 GenPipes code cleaning - cleaned some 'import' calls : stop using 'import *' and specify which modules/classes/functions to import instead
       2699bc94 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into sanity_check_mode
       975e0ff3 GenPipes Core - updated pipeline.py to catch , relative to changes in config.py. Added the sanity-check mode : the pipeline now does not stop at fist met error/exection : wait untill the end to print the errors
       1f628d35 GenPipes Core - updated config.py wth PO's code : catches  instead of
       0f4c8605 Software upadte - updated mugqic_tools installation script with version 2.2.3
       68f32a97 Software update - updated R_Bioconductor.sh : added binless and DoubletFinder packages to the installation, updated the installer depending on the version of R
       2bca735e Saoftware update : update cellranger-atac.sh with the latest versoin 1.1.0
       14783218 Software update - added th shebang to the C3G wrappers for installed binaries
       697eaec2 Genome installation - corrected typo in 'lambda_page'
       e7010eb8 Software update - added script to install gemBS in C3G softwatre stack
       9c4f2954 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       93f92cdc Software update - updated installation scripts for LUMPY-SV, mugqic_tools & SMRTLink
       39b691a3 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       06fdb584 Software update - updated MultiQC installation script : shebang of MultiQC scripts now uses #!/usr/bin/env python
       e66f717b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1e4c96fe Software update : added the installation of methylKit packages in R_Bioconductor.sh - also added the generation of a file listing all the packages installed for all the new installed R versions
       251f418c GenPipes utils - updated nanuq2mugqic_pipelines script : added support for iSeq projects
       742aa73d GenPipes utils - updated nanuq2mugqic_pipelines script : added -rn/--run parameter, standing for Nanuq run ID, to fetch only readsets procesed in specified run(s)
       55551719 Software update - updated install_module.sh so that it does not wrap nor patch the executable binaries when installing on DEV space
       d4fb87c7 Software update - added MiXCR v3.0.5 installation script - MiXCR: a universal tool for fast and accurate analysis of T- and B- cell receptor repertoire sequencing data
       f020d920 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3e7d563e Software update - updated some installation scripts with the latest version of software : bowtie v1.2.2 - bowtie2 v2.3.5 - MultiQC v1.7
       09335ff9 install_genome.sh - added the generation of the 100-bin GC content file for bisulfite genome references
       4e1331fa Merge branch 'master' of bitbucket.org:mugqic/genpipes
       04047dfc Genome installation : added Bisulfite Genome Reference creation and indexing with Bismark
       8916af04 updated LIBDIR in R_Bioconductor.sh
       0ad84563 update Bismark version to 0.21.0
       eff80025 updated python installation script with version 3.7.3 and added pip as a symling to pip3
       db27d82a Added fastq software installation script
       e1f627ff Version bump to 3.1.5-beta
       1786fb37 Version bump to 3.1.4

  Édouard Henrion <henrione@beluga4.int.ets1.calculquebec.ca>      1 commits

       650797a2 GenPipes - DEBUGGING - DNASeq SamToFastq & SymLink steps corrected + working bfx/bash_cmd.py

  Édouard Henrion <henrione@gra-login2.graham.sharcnet>      4 commits

       9e9adc78 GenPipes JSON - debugged call to job2json when output_dir is different that '.'
       b7eb0ed2 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       754a7d98 GenPipes software update - updated patching in install_module.sh
       5bf7663e GenPipes software update - updated R_Bioconductor.sh with the latest requested libraries & updated patching

  Édouard Henrion <henrione@ip18.m>      3 commits

       140b411d GenPipes Scheduler - Corrected bug in job2json call
       c77a5352 Merge branch 'eh_methylseq_single_end' of bitbucket.org:mugqic/genpipes into eh_methylseq_single_end
       a592ddfd GenPipes MethylSeq - updates for single-mode

  Édouard Henrion <henrione@ip20.m>      1 commits

       41ed7a2c Genome Installation - updated install_genome.sh with common grep version

  ehenrion <edouard.henrion@mcgill.ca>      83 commits

       25f72f00 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1578931756922 (pull request #163)
       64b7414f DNASeq.py - removed the use of 'os/path.abspath' in call of 'ln()'
       3c750f00 bash_cmd.py - added import os
       7a05d36e DNASeq - Skewer trimming call to ln() upadted without 'sleep' variable
       e9213a63 bash_cmd.py - remove call to "sleep" in the ln() function. Replaced it by a call to "ls" in order to flush the cache
       01154915 GenPipes - BUG correction - trailing space removed in output files of DNASeq skewer trimming step
       1fdf5799 Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578602990982 (pull request #162)
       b74f38dc BUG correction : corrected bash_cmd.py "ln" function with correct format keys when buiding command
       07d826fd Merged in ehenrion/bash_cmdpy-edited-online-with-bitbucket-1578588554527 (pull request #158)
       828d5736 Bug Correction : corrected 'ln' function in bash_cmd.py to avoid a python "TypeError: cannot concatenate 'str' and 'int' objects" exception
       42055827 Merged in ehenrion/dnaseqbelugaini-edited-online-with-bitbu-1576097030188 (pull request #147)
       fcaef48a GenPIpes - DEBUGGING - added slurm-comprehensive walltime for picard_sam_to_fastq in dnaseq.beluga.ini
       52021009 Merged in eh_quick_fixes (pull request #143)
       815632f0 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575918024470 (pull request #141)
       92a7c117 GenPIpes - dnaseq.py : bug correction - typo removed
       44d9eb52 Merged in ehenrion/dnaseqpy-edited-online-with-bitbucket-1575399438745 (pull request #132)
       15266427 Merged in ehenrion/commonpy-edited-online-with-bitbucket-1575398053532 (pull request #131)
       482a0bba GenPipes - bug correction in pipelines/common.py : corrected name for bam and fastq files created by the pipeline : don't depend and raw file names anymore, but built from sample & readset name given in the readset file
       2a25960a GenPipes - bug correction in pipeline/dnaseq.py : corrected sym_link_fastq, skewer_trimming & bwa_mem_picard_sort_sam steps, regarding the path of the fastq files when they have to be determined from the readset.bam
       72e7f1e5 GenPipes - bug correction in pipelines/common.py : corrected the path where the sorted bam files as well as the raw_reads fastq files(from sam_to_fastq) should be written, i.e. within the output directory provided at the pipeline execution
       c946bb0e Merged in ehenrion/schedulerpy-edited-online-with-bitbucket-1575392160612 (pull request #129)
       9bf49cbe Merged in ehenrion/nanuq2mugqic_pipelinespy-edited-online-w-1575392552042 (pull request #130)
       3893e10c GenPipes : nanuq2mugqic_pipelines.py : bug corrected - typo in seq_type selection
       2f38d5d8 GenPipes - corrected scheduler.py : removed unwanted sed command in --no-json context
       a876d020 GenPipes - corrected scheduler.py with missing argument line
       8b8bca82 Merged in dev (pull request #125)
       7cd7b1cd RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - dowgraded trinity version to 2.0.4_patch
       656db125 RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - corrected typo in trinity version...
       01c2bee3 RNASeq Denovo Assembly - rnaseq_denovo_assembly.base.ini - updated trinity version to 2.2.0_patch which contains a patch (from C3G developppers) to avoid  printing buffer 'hiccup'
       b9325d16 CIT - RNASeq Denovo Assembly - update default cluster_walltime to 4:00:00 in cit.ini
       8235190e Merged ehenrion/cit-dnaseq-updated-cit_cluster_walltim-1574180387297 into dev
       238abec9 CIT - DNASeq - updated 'cit_cluster_walltime' for gatk_callable_loci step in cit.ini
       6a891c20 DNASeq - gatk_callable_loci - adjusted memory and cpu allocation in dnaseq.beluga.ini
       dec27c8e CIT - DNASeq - Adjusted trimmomatic resource allocation through java_other_options, threads settings and mem-per-cpu use
       3d342433 CIT - RNASeq - corrected 'cluster_walltime' for wiggle  step in cit.ini
       a24e7b18 RNASeq - wiggle step - adjusted memory allocation using 'mem-per-cpu' instead of 'mem'
       3e3ee1d1 CIT - RNASeq - updated java threads to 5 through 'java_other_options' in rnaseq.beluga.ini
       1a1d0f1e CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       bae00cdd CIT - RNASeq_light - updated default cluster_walltime to 4:00:00 in cit.ini
       cc5df685 CIT - HiCSeq - redefined hicup_align walltime and mem-per-cpu in cit.ini
       5cd0d02f CIT - Pacbio assembly - set specific walltime for smrtanalysis_summarize_polishing in cit.ini
       1c197913 HICSeq pipeline - trimmomatic resources revisited in hicseq.beluga.ini : removed buffer_size
       29024f2c CIT - RNASeq denovo Assembly - set localfit to true in cit.ini for differential_expression_deseq
       622b41ad CIT - Pacbio Assembly - redefined walltime for pacbio_tools_split_reads in cit.ini
       09e22354 RNASeq denove Assembly - edited rnaseq_denovo_assembly.beluga.ini adjusted memory assignment for insilico_read_normalisation_readsets
       3a25fef0 RNASeq denovo Assembly : edited rnaseq_denovo_assembly.beluga.ini with beter resources assignement
       b4a3e981 Software update - trimmomatic.py - updated trimmomatic command with the use of the  java_other_options parameter provided by the ini files
       47e5e41c Merged in ehenrion/rnaseq-metricspy-ihec_metrics_rnaseq--1573489397505 (pull request #124)
       87c3190b RNASeq - metrics.py - ihec_metrics_rnaseq : added file in the input_file list to correct job dependency
       fc680b19 RNASeq cufflinks - correcting dependencies for ihec_metrics : needed rnaseqc report file to be added to the output_file list
       25f6cc3f CIT - Pacbio assembly - set threads in beluga.ini
       aa3d185d CIT - RNASeq_light - set threads for trimmomatic in beluga.ini
       b1992425 CIT - RNASeq_light - redefined walltime for callisto_count_matrix in cit.ini
       24d79ff4 CIT - Pacbio Assembly - redefined walltime for preassembly in cit.ini
       3c8f736b CIT - RNASeq denovo Assembly - redefined walltime for align_estimate_abundance in cit.ini
       3134f05c CIT - RNASeq denovo Assembly - updates cit_cluster_walltime to 24h
       24351cb0 CIT - RNASeq denovo Assembly - redefined walltime for transdecoder and align_estimate_abundance in cit.ini
       c17d7918 CIT - Pacbio Assembly - redefined walltime for smrtanalysis_filtering in cit.ini
       19ed125f CIT - DNASeq - redefined walltime for sambamba_merge_realigned in cit.ini
       d92aa276 CIT - DNASeq High Coverage - redefined trimmomatic walltime in cit.ini
       74979fa5 CIT - HiCSeq - redefined trimmomatic walltime in cit.ini
       ed4fb136 CIT - RNASeq_light - redefined walltime for trimmomatic in cit.ini
       0cf351d1 RNASeq denovo Assembly pipeline - rnaseq_denovo_assembly.py - corrected command quoting in trinity step
       8386936b RNASeq denovo Assembly Pipeline - rnaseq_denovo_assembly.beluga.ini - corrected max_memory parameter for trinity step : using Gb instead of Mb
       373bd60c RNA-Seq de novo assembly pipeline - rnaseq_denovo_assembly.beluga.ini - corrected Jellyfish memory parameter : using Gb instead of Mb unit
       d2e4c46b HiC-Seq pipeline - hicseq.base.ini - corrected typo in HiCUP module name
       ceaae534 HiC-Seq pipeline - hicseq.base.ini -  update HiCUP version to v0.7.2
       7860ea8b dnaseq.py - corrected file names dependencies
       c198c42e dnaseq.py - corrected wrong file extension in metrics_vcf_stats step
       bc952b66 HiC-Seq pipeline - corrected typo in cram_output options definition - hicseq.base.ini
       d8c5f581 rnaseq_light.cedar.ini : corrected typo in cluster_server assigned value
       19bf83c4 rnaseq_light.mp2b.ini corrected minor typo on cluster_server
       27f7080b rnaseq_denovo_assembly.py : corrected differential expression jobs with samples
       5010b026 corrected typo in ampliconseq.beluga.ini
       1f160172 DNASeq High Coverage README.md updated with `picard_fixmate` step documentation
       eda0617f GenPipes - Readme updated for `picard_fixmate` step in dnaseq_high_coverage.py
       a05f490f Merged in eh_methylseq_single_end (pull request #115)
       6c88ebf3 Merged in sanity_check_mode (pull request #114)
       fa256237 rnaseq_denovo_assembly.beluga.ini :  updated [insilico_read_normalization_all] with missing core allocation
       305d96e9 README.md updated : removed Guillimin settings section from README
       56d668c8 README.md updated : corrected the setting of MUGQIC_INSTALL_HOME_DEV for mp2b
       8d4f1b47 rnaseq_denovo_assembly.mp2b.ini - corrected resources requirments : all the steps now run on one single node !
       33909cfb chipseq.base.ini - updated MultiQC version to 1.6 to support Homer

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       6c76ddc3 Merged in lefebvref/rnaseq_denovo_assemblybaseini-edited-onl-1557891680322 (pull request #110)
       53bd32ae Previous default expression in base.ini would end up basically only retaining transmembrane proteins. Not good.

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      2 commits

       5bea83ca Corrected error with beluga ini for RNAseq denovo
       02b131cf Corrected minor error in the help message that said that the default protocol was cuflinks. Starting from version 3.1.4, stringtie is the default protocol

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      10 commits

       153d4c4d Merged in Jose-Hector-Galvez/update-to-the-rnaseq_light-readmemd-file-1569265173884 (pull request #120)
       e24a0661 Update to the rnaseq_light README.md file to address issue brought up by GSoD collaborator.
       5f4086f3 Merged in Jose-Hector-Galvez/readmemd-edited-online-with-bitbucket-1565184718156 (pull request #118)
       ed95482e README.md edited to remove warning labels.
       341cab2f Merged in rnaseq_denovo_jhg (pull request #103)
       05e7dfe8 Merged in rnaseq_jhg (pull request #112)
       9a24beaf Merged dev into rnaseq_jhg
       6f98c769 Merged master into rnaseq_denovo_jhg
       55051474 Merged master into rnaseq_jhg
       160079da Merged master into rnaseq_denovo_jhg

  José Héctor Gálvez López <hgalvez@beluga3.int.ets1.calculquebec.ca>      3 commits

       e618184d Corrected error in the beluga ini for RNAseq de novo
       e9610a9f Corrected error in the beluga ini for RNAseq de novo
       ceade422 Corrected error in the beluga ini for RNAseq de novo

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       cc5eac7f add header in the cram file (all base.ini)
       92b5c4f2 add cram to tumor_pair pipeline
       78224133 implement the launch of the entire pipeline if no step argument is given
       689732a1 add cram to dnaseq_high_coverage pipeline
       71569ce8 add cram to methylseq pipeline
       c242a05a add cram to rnaseq pipeline
       70ccb037 add cram to hicseq pipeline
       0e788734 add cram to chipseq ini files
       830cd38d add cram to dnaseq ini files
       f9d780a6 add cram to ChipSeq
       5f1d736f correct bugs in dnaseq light
       0313f9cc add cram creation to dnaseq pipeline
       bfde2145 add generic function to create CRAM from BAM
       19565149 make samtools view's output not mandatory removable

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      2 commits

       52cd608e temporary fix for tbi missing output in haplotypecaller when running only 1 job
       6a98813c Merged in cram (pull request #111)

  Paul Stretenowich <paul.stretenowich@mcgill.ca>      9 commits

       ee86893a Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       0fd7a889 Nextflow install script
       cc8694c1 Merged in rnaseq_rnaseqc (pull request #116)
       d59ea5ee General - Correcting genome ini path
       caa989b8 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into rnaseq_rnaseqc
       8e90fcf0 RNA-Seq Pipeline - Typo update
       6bfe800b RNA-Seq Pipeline - Typo update
       bb8a0cdd RNA-Seq Pipeline - Removing verify_bam_id step
       7eb13c4a Chip-Seq Pipeline - Adding homer_make_ucsc_file to mp2b.ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      12 commits

       a8d1aae7 Merged in add_container (pull request #122)
       9b20dde9 Merged in cit_fix_cedar_2019-04 (pull request #109)
       4eb85d86 gatk4.py added 2 spaces in multi line string
       ef09195b ampliconseq.py edited online with Bitbucket
       9e04a50d pacbio_assembly.cedar.ini edited online with Bitbucket
       2b023062 Merged in log_repport_ci (pull request #108)
       04e9b6ba Merged in log_repport_ci (pull request #107)
       f171ba18 Merged in log_repport_ci (pull request #106)
       0da3f3ec Merged in log_repport_ci (pull request #105)
       f9458d8c Merged in log_repport_ci (pull request #104)
       1b7e74f6 Merged in fix_multiqc (pull request #89)
       0e66455b Merged in fix_multiqc (pull request #77)

  Pierre-Olivier Quirion <poq@beluga3.int.ets1.calculquebec.ca>      2 commits

       feeaf91d Log when loading Job. More robust output parser.
       5181b5cc Log when loading Job. More robust output parser.

  P-O Quirion <pioliqui@gmail.com>      158 commits

       161d1d19 fixed motifMaker ini: used default memory
       ba32a478 typo in waltime kallisto_count_matrix
       50b824bd trimmomatic thread ini
       20c0f567 Merge branch 'dev' of bitbucket.org:mugqic/genpipes into dev
       795edd94 string tie ini
       d4e255de sync thread in java for trimomatic
       9bd93fa9 differential_expression_and_goseq_rsem return list list
       42bf1d80  fix scheduler logging; reduce sleep to 0.1 sec
       0902f427 ignore .pyc
       bca87cc0 remove the methylkit_differential_analysis step form the pipeline
       f9cef1b3 add gitignore
       dd0267e2 missing import
       3d6b4e70 log regression scheduler
       f4ab24a8 fix beluga ini
       6d5d769c fix beluga ini
       ef68329f Merge branch 'add_container' into merge
       a872d634 Merge branch 'master' into dev
       d31246f4 fix batch system account for container
       313fcb20 remove commented code
       281e3cd7 Merge branch 'add_container' of bitbucket.org:mugqic/genpipes into add_container
       ebb67ba7 Add default config for container ini
       116615f7 fix regression
       7542a098 fix scheduler when container not in use
       8c9d6e3c use mugqic stack java in rnaseq on cedar.ini
       bb6fc694 simple default ini for container
       b82caced default ini for container
       331715de Basic container ini file
       b787ff7e Working singularity version tested on graham
       2fe8f475 WIP exec line prototype
       b2b4d970 add container option
       215d3527 Merge branch 'cit_fix_cedar_2019-04' of bitbucket.org:mugqic/genpipes into cit_fix_cedar_2019-04
       4092e778 amplicon R tweak
       7b389fcc missing quotes in pacbio
       138f3356 amplicon memory
       32627c87 amplicon memory
       141554d5 sleuth again wrong name
       8c34b2be renaming sleuth to differential_expression_sleuth
       963d105a fix typo in ampliconseq ini
       ff92d21f fix amplicon method typo
       da2cdf80 tweak default mem for rnaseq denovo and light
       75e2f696 more ampliconseq deadlock: otu_assigning
       1cfac735 parallel_assign_taxonomy_uclust.py deadlock
       6ecfef8e rna fix
       cecb5144 update cit
       b4997b3a packbio tmpdir again
       c3925769 chr19 option for cit chicagorun
       46d9149f fix log report on output id
       5cfc33cb fix log report on output id
       0cc3e179 update waltime in amplicon and regex in log
       5cc0ce57 migrate form testdata ini to cit
       2c02e926 mem teak
       51bf13da tweak to get in bynode queue
       3f05dbe8 tweak ini for cit
       078d0e0a dnaseq ini tweak for slurm; sh_create_baitmap \\
       cf79f461 timout amplicon
       410185e8 remove mpileup_metrics_vcf_stats step from mpileup
       9d100d8c pbutgcns_wf need full path an tmp need to exist
       ad31e6e8 ini file tweak
       ac6bf87a tweak mem usage
       b4f2cf80 ini tweak plus perl module load in mummer
       3ae1b09b tweak mem usage per-cpu for ampliconseq
       bb8fc310 Ensembl v 90 for GRCh38
       890cd484 mummer using deprecated perl syntax
       16326cfb more default time for rnaseq cufflink
       3d1aa247 mugqic_tools/2.2.4 in rnaseq_light
       ffcf62f6 samll sample for both gatk call in variant_recalibrator
       7b25d7bc change R_Bioconductor version/cleanup code
       bc97eae0 update rnaseq_light references
       505ffeb8 tweek memory usage in dnaseq_high_coverage
       3ce4bb61 add time to ampliconseq
       c1b5a6a3 reformat errors in the core.config module
       05aa94b4 redirect small_sample_option to the right function
       3d51765c tweak cit for dnaseq
       720c3fa7 fix ampliconseq dada2
       390ceee0 tweak VariantRecalibrator for testing
       f74afdff fix slurm_time_to_datetime
       733f258a set cit ini file and config
       7d9c898b add runtime to csv output
       86e84f3e Add GENPIPES_CIT env var to run in testing mode
       4cb6e240 fix all slumtmpdir sting renaming problem
       8f3652ff update mugqic_tools in hicseq
       2107daaf second run of pipeline fixing
       458e919e remove extra gz extention from anotation
       690acf41 flash_log is a list of one
       b829c316 update hicseq.base
       81a4690e add interactive option in smrtanalysis bash cmd
       78552172 fix regression
       18aa48bf fix scheduler when container not in use
       ae7824c0 use mugqic stack java in rnaseq on cedar.ini
       b8eee7fc simple default ini for container
       2522d29d default ini for container
       60de9f18 Basic container ini file
       a5b2901d Working singularity version tested on graham
       5229f0ee WIP exec line prototype
       56e0a9b0 add container option
       27a2cb37 amplicon R tweak
       738c59da missing quotes in pacbio
       d73e1659 amplicon memory
       a770708a amplicon memory
       80498da4 sleuth again wrong name
       53c5205b renaming sleuth to differential_expression_sleuth
       1e5d20e0 fix typo in ampliconseq ini
       7bb900da fix amplicon method typo
       85b51a0d tweak default mem for rnaseq denovo and light
       7df9cae8 more ampliconseq deadlock: otu_assigning
       777bb337 parallel_assign_taxonomy_uclust.py deadlock
       24bd06ed rna fix
       a9e2370d update cit
       b120e4cc packbio tmpdir again
       8e1835fd chr19 option for cit chicagorun
       eec3f953 fix log report on output id
       5488adc0 fix log report on output id
       34960c03 update waltime in amplicon and regex in log
       ca1fe4d3 migrate form testdata ini to cit
       eca1836a mem teak
       058f0717 tweak to get in bynode queue
       10ff44a8 tweak ini for cit
       c7568c6e dnaseq ini tweak for slurm; sh_create_baitmap \\
       1c50e41d timout amplicon
       ff1a9c7d remove mpileup_metrics_vcf_stats step from mpileup
       d8d8a3dd pbutgcns_wf need full path an tmp need to exist
       48ff75bb ini file tweak
       68d6f468 tweak mem usage
       160a1234 ini tweak plus perl module load in mummer
       ee6a7aa3 tweak mem usage per-cpu for ampliconseq
       f0acafe6 Ensembl v 90 for GRCh38
       23788c3b mummer using deprecated perl syntax
       14d0500a more default time for rnaseq cufflink
       4e3a8b07 mugqic_tools/2.2.4 in rnaseq_light
       624781cf samll sample for both gatk call in variant_recalibrator
       8307e596 change R_Bioconductor version/cleanup code
       069dd4aa update rnaseq_light references
       3e01c15d tweek memory usage in dnaseq_high_coverage
       c0c5905f add time to ampliconseq
       d29d9267 reformat errors in the core.config module
       43fdedaa redirect small_sample_option to the right function
       cfefb02a tweak cit for dnaseq
       509d3330 fix ampliconseq dada2
       927a732f tweak VariantRecalibrator for testing
       0a2659d1 fix slurm_time_to_datetime
       0f48c2df set cit ini file and config
       20c85214 add runtime to csv output
       02ce535d Add GENPIPES_CIT env var to run in testing mode
       dbbb0cd3 fix all slumtmpdir sting renaming problem
       19c6c30a update mugqic_tools in hicseq
       27e51645 second run of pipeline fixing
       32a1a40d remove extra gz extention from anotation
       4edbfc57 flash_log is a list of one
       a396c5ff update hicseq.base
       d806fced add interactive option in smrtanalysis bash cmd
       d9567ac7 add space
       31ddd6cf add check for output job id
       63431c42 get real path in log_report
       bec1796c file exist file missing mixup
       4120c3ab bug fixes
       3e94acda get pipeline info in tsv file can also mute stdout
       29f5670c run dnaseq multiqc on all samples at once
       8d732173 add indentation

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       43f70c63 fix to preprocess step

  Robert Syme <robsyme@beluga4.int.ets1.calculquebec.ca>      1 commits

       85814368 We don't need to track compiled python .pyc in git

  Robert Syme <rob.syme@gmail.com>      7 commits

       30ea46dc Merged in spacesfix (pull request #117)
       0c4d6ce4 Remove spaces from Rscript commands.
       1f3c0f66 More trailing space cleanup.
       534124c3 Cleanup trailing spaces
       16b7d01d Merged in amplicon-pathfix (pull request #113)
       63d804af Add support for --output-dir argument
       eb8389c3 Trailing whitespace cleanup

  Rola Dali <rola.dali@mail.mcgill.ca>      31 commits

       92202e88 dnaseq_high_coverage.base.ini edited Varscan version to 2.4.3
       3f0b7320 dnaseq.mp2b.ini edited online with Bitbucket
       07b6376b dnaseq.beluga.ini edited online with Bitbucket
       a0bad077 dnaseq.cedar.ini edited online with Bitbucket
       f379676f dnaseq.beluga.ini edited online with Bitbucket
       2805dc69 dnaseq.cedar.ini edited online with Bitbucket
       a4c34187 pacbio_assembly.beluga.ini edited online with Bitbucket
       bdd8b724 dnaseq.beluga.ini edited online with Bitbucket
       b648fb5a dnaseq.cedar.ini edited online with Bitbucket
       287e82b8 dnaseq.beluga.ini edited online with Bitbucket
       a760a50a dnaseq.cedar.ini edited online with Bitbucket
       62de3fb7 hicseq.beluga.ini edited online with Bitbucket
       1fd77878 hicseq.cedar.ini edited online with Bitbucket
       42eb77ac ampliconseq.cedar.ini edited online with Bitbucket
       fabb3d5c ampliconseq.mp2b.ini edited online with Bitbucket
       193c97ba ampliconseq.beluga.ini edited online with Bitbucket
       5c4c908a ampliconseq.cedar.ini edited online with Bitbucket
       59491ffe dnaseq.cedar.ini edited online with Bitbucket
       9c9a75fd dnaseq.beluga.ini edited online with Bitbucket
       a2e7e80f rnaseq.beluga.ini edited online with Bitbucket
       d1777dd7 hicseq.beluga.ini edited online with Bitbucket
       5a6a1987 hicseq.cedar.ini edited online with Bitbucket
       855a02c3 hicseq.cedar.ini edited online with Bitbucket
       a4cb4c16 hicseq.beluga.ini edited online with Bitbucket
       efe78e4d README.md edited online with Bitbucket
       c0b67c17 README.md edited online with Bitbucket
       5c782512 README.md edited online with Bitbucket
       982f934f README.md edited online with Bitbucket
       4402aeff README.md edited online with Bitbucket
       ed0fd9e4 README.md edited online with Bitbucket
       b97304ef README.md edited online with Bitbucket

3.1.4        Tue Mar 26 14:03:32 2019 -0400        198 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      69 commits

       cea58560 Updated all the GenPipes .base.ini files with the latest verison of mugqic_R_packages i.e. 1.0.6
       892cef59 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       cc01f4ee MethylSeq - added reference to Bismark in markdown config_and_references.md file for a better report document - also updated the metrics table assignement in MethylSeq.ihec_sample_metrics_report.md and bfx/report/MethylSeq.ihec_sample_metrics_targeted_report.md
       540ee0ed Merge branch 'master' of bitbucket.org:mugqic/genpipes into slurm_report
       3ab4a6f0 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f8c58f3d updated mugqic_tools version to 2.2.2 in ampliconseq.base.ini
       4a100ec8 updated R_mugqic_packages.sh with version 1.0.6 of the mugqic R packages
       1d490a7a updated pipeline inis file with the latest version of mugqic_tools i.e. v2.2.2
       d1ad8344 MethylSeq - MethylKit DMR - to fit with C3G code standards, getting rid of bfx/methylkit.py and use bfx/tools.py instead to call methylkit.R, which is part of mugqic_tools
       d5e1184d MethylSeq - methylKit DMR - correted typo in bfx/methylkit.py
       a989966c MethylSeq - methylKit DMR - correted typo in methylseq.base.ini
       ba9c8b13 MethylSeq - MethylKit DMR - change the call to methylKit.R : now using R instead of Rscript
       c6c82bf3 Pac Bio Assembly pipeline - updated the bfx wrapper circlator.py : make sure th eoutput is removed before launching circlator, otherwise it will throw an error
       a763f6d8 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       9650b322 PacBio Assembly pipeline - minor update in pacbio_assembly.py for code reading purposes
       e41a72ee Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       0ddffdd5 updated svaba.sh installation script with latest version of SvABA : 1.1.0
       43e75270 MethylSeq pipeline - DMR analysis - updated R module version in the base.ini file to make sure methylKit is available
       3576f8d0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       90d9206a updated mugqic_tools version to 2.2.2 in methylseq.base.ini
       34f4ccb4 update mugqic_tools.sh with the latest version : 2.2.2
       f1a34d08 updated install_module.sh to better handle LD_LIBARAY_PATH in the LIBDIR
       0fb5b2cc Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0aa9f92c corrected Octopus installation script
       a5a7bd9e Pacbio Assembly pipeline - corrected basemodification and MotifMaker bfx subroutine wrappers : now sourcing /etc/setup.sh
       bfefcf72 modified R module in install_genome.sh : now using mugqic/R_Bioconductor/3.5.1_3.7
       805760f7 updated mugqic_tools.sh with version 2.2.1
       02d4ef83 Adding installation scripts for CMake and Octopus
       10a90429 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       193fe284 MethylSeq pipeline - DMR analysis - change file namings in prepare_methylkit step
       f6f936ff MethylSeq pipeline - DMR analysis - adjusting walltime for filter_snp_cpg step
       aca1e34b MethylSeq pipieline - DMR analysis - clean code in bfx/tools.py
       732a39ac updated wrapping and patching in install_module.sh
       56891fac updating software install scripts with newer versions
       82723b7a bedtools.py - removed unnecessary argument to bedtools.coverage subroutine
       d126c1b7 MethylSeq pipeline - methylseq.py - fixing inputs assignment in IHEC metrics step
       4eead82c MethylSeq pipeline - methylseq.py - correcting sort_sam call and affecting job name to bismark methylation call job
       27f85ae6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e23a9821 BFXDEV-429 - updated syntax of temporary diretory after sourcing of /etc/setup.sh for smrtanalysis tools
       53e51310 updated resources/modules/R_Bioconductor.sh script : refine LIBDIR definition - revised wrappers creation
       3f209128 updated install_module.sh : refine create_c3g_wrappers subroutine to catch executable binaries - updated LIBDIR definition for centOS cases
       58286be3 updated delly.sh script : lasted version 0.8.1 & added DELLY_PATH in the modulefile
       8a1084f9 udpated lumpy-sv.sh script : added a sed command to make lumpyexpress.config compliant with our environment
       1c7b2c44 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       58bde8f8 PacBio Assembly - updated bfx/smrtanalysis.py as part of the solution of BFXDEV-429
       7ee6531a Updated syntax to tmp_dir to finish fixing BFXDEV-429 - first usefull related update was made in commit ad3caf1, when bumping to version 3.1.3
       807bee4f minor - updated syntax
       e93e06f2 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       dafb3c08 removing occurences of DEV paths within the ini files of GenPipes
       3d256c1b Adjusting .ini section for some htslbfx subroutines which were pointing to [DEFAULT] section of the .ini file
       3fe0d06d BFXDEV-752
       25bc30d8 Adjusting .ini section for some qualimap ray sambamba and other subroutines which were pointing to [DEFAULT] section of the pipeline .ini files
       6d4dd16a corrected some typo in jsonator.py to have fastq & bam file paths properly handled
       630de05a Adjusting .ini section for some htslib subroutines which were pointing to [DEFAULT] section of the .ini file
       5b03880f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       891780f3 comment line removed in install.module.sh
       f6720fb4 Added ctpr.sh script to install CTPR package in C3G software stack
       0a3e8ab1 Adjusting .ini section for some blast * blat subroutines which were pointing to [DEFAULT] section of the .ini file
       3ebb9cc8 BFXEDV-529 - edited RNAseq to consider star_cycle_number instead of cycle_number for star.index step - also ajusted read_size value for sjdbOverhang value completion
       54fd0401 Homo_sapiens.GRCh38.sh - Added the creation of a Ribosomal RNA interval list file, for use of Picard Tools CollectRnaSeqMetrics
       fe136f79 Adjusting .ini section for some bcftools subroutines which still were pointing to [DEFAULT] section
       709347d2 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       4344c947 changing samples assignment in trimmomatic jobs for a better analysis JSON creation
       321c65d9 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       24e90fa0 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9cc90a8c added genpipes.sh
       0157c923 updated GATK installation script
       ea084430 Version bump to 3.1.4-beta
       7e318dc1 Version bump to 3.1.3

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       9cc62b8b MethylSeq Pipiline - MethylKit DMR analysis - updated Cedar .ini file with relevant cluster resources, changed filter_snp_cpg subroutine int bfx/tools.py to use bedops instead of bedtools, added bedops module in the .base.ini file
       e26c2afd MethylSeq Pipeline - updated filter_snp_cpg call with more cluster resources and removed pipes to avoid 'Broken pipe' error...
       7d77cd8b MethylSeq pipeline - DMR analysis - fixing cedar ini file + cleaning code before pull request

  ehenrion <edouard.henrion@mcgill.ca>      14 commits

       3ef3a043 pipeline.py  : edited the help so that it actually shows that SLURM is the default scheduler used by GenPipes
       d3edaa68 methylseq.py corrected report files in metrics step
       86646c43 Merged in slurm_report (pull request #73)
       b2dbb469 job.py : corrected typo
       0961c1ef job.py : added missing setter class for name and multiqc_files attributes
       4f48cc68 Merged in methylseq (pull request #67)
       10bcac61 dnaseq.cedar.ini : resolved conflicts
       a18369a4 rnaseq_denovo_assembly.base.ini : corrected missing variables
       d80c1788 tumor_pair.base.exome.ini - updated COSMIC path
       1a66a54e tumor_pair.base.ini - updated COSMIC path
       a3c68a37 rnaseq_light.base.ini : commetn adapter_fasta because it was pointing to some old place
       0f799c37 BFXDEV-737 - updated ucsc.py for a better handling of Danio Rerio GRCz11 genome build
       b768d04f BFXDEV-734 - updated ucsc.py for a better handling of Mus Musculus GRCm38 genome build
       54dab24f AmpliconSeq - reassign silva_db to correct path for dada2 analysis - ampliconseq.base.ini

  Hector Galvez <jose.hector.galvez@computaitonalgenomics.ca>      28 commits

       7da61f05 Corrected slurm parameters for stringtie jobs, including stringtie_abund which was causing errors
       fb2d51da Removed DEV from sleuth configuration
       68878b13 Added note to RNAseq README regarding switch to stringtie protocol by default
       94342a5c Merge branch 'rnaseq_light_jhg' of bitbucket.org:mugqic/genpipes into rnaseq_light_jhg
       68bf8861 Specified genome support for RNAseq_light in the pipeline's README
       5088cbdc Corrected INI file to inlcude latest version of Mugqic_tools
       e455bbc9 Modified INI file to include latest version of MUGQIC_tools as well as more memory for kallisto steps by default
       e974d450 Merge branch 'rnaseq_jhg' of https://bitbucket.org/mugqic/genpipes into rnaseq_jhg
       49bbed54 Corrected minor parameters for abacus on the RNAseq ini
       330c861c corrected merge issue with the rnaseq_ligth cedar INI
       0d7a7d7f Corrected minor issues with rnaseq_light pipeline, and added documentation to the README
       160fd60e Corrected some issues with the UCSC script that was causing issues with non GRCh38 genomes
       52daaeaa Added documentation for stringtie protocol. Corrected some problems with the base inif file. Minor corrections to the stringtie script and rnaseq.py script.
       214930b4 Fixed merge issues
       497139c6 merged to master and resolved conficts
       95be45a4 Merged to latest version of master branch
       90e40bad Corrected a mistake on the methylseq ini files for Cedar and Mp2b
       5b684178 Corrected CVMFS Stringtie module in INIs
       615ad858 Fixed merge conflict with the cedar ini
       548a27f8 Merge branch 'master' into hgalvez
       b6a9e4e8 Second testing version stringtie/ballgown pipeline. Includes ballgown now, as well as corrections for stringtie. Additionally, stringtie is now default protocol for RNA-seq.
       08bbb3cc Fixed a small typo with the last commit
       1abef45d Fixed bug for the order or arguments passed to kallisto
       abbbfb2f Updated version of Bioconductor used by default in the rnaseq_light pipeline (modified base.ini file).
       a97c715a Merge branch 'hector-cedar' of https://bitbucket.org/mugqic/genpipes into hector-cedar
       7a877b31 Minimally functional rnaseq_light pipeline with the added sleuth step. Still requires further testing on clusters beyond abacus. Also, some harmonization of the required genome reference files would be advisable.
       3c142b81 Merge branch 'master' into hector-cedar
       2ccccdea Eliminated warning message from rnaseq_light pipeline for single read samples

  jose.hector.galvez@computationalgenomics.ca <hgalvez@abacus2.ferrier.genome.mcgill.ca>      2 commits

       15b6ecbd Fixed issues with stringtie merge and abund. Ready to test in other servers
       4f3e5e8b Corrected job.py to allow for the definition of the multiqc_files parameter, should fix error with multiqc and Kallisto

  Jose Hector Galvez <Hector@rhododendron.genome.mcgill.ca>      10 commits

       4a561483 latest modifications to the tools.py script
       9363e82a Fixed minor typo on rnaseq.py
       397a903e Fixed indent issue in the stringtie job creation
       35f29048 Fixed minor bug in the stringtie job creation
       ebe3d347 Fixed minor bug in the stringtie job creation
       ffd2a51e Fixed stringtie module location in the rnaseq.base.ini file
       a76c8e4f Fixed sample name issue within stringtie.py
       c6b938f5 corrected minor typo
       18c69678 correct minor typo
       01cc53d0 First commit adding stringtie functionality, still testing

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      22 commits

       447c061a Merged in rnaseq_light_jhg (pull request #75)
       daabf762 Merged in rnaseq_jhg (pull request #76)
       de0d01eb Merged master into rnaseq_jhg
       7293ec07 Merged master into rnaseq_light_jhg
       e032d8fd Merged in rnaseq_jhg (pull request #69)
       fecf7c71 Merged master into rnaseq_jhg
       d687b0fe Merged in rnaseq_light_jhg (pull request #68)
       e78a608f Merged master into rnaseq_light_jhg
       9225ebbd Merged master into rnaseq_jhg
       3868e89f Merged master into rnaseq_jhg
       573a0e0e Merged master into rnaseq_light_jhg
       640d57c5 Merged master into rnaseq_jhg
       77c6491f Merged master into rnaseq_jhg
       16c4eb7d Merged in missing-cedar-ini-files (pull request #64)
       c9573989 Merged master into hgalvez
       072c23de Merged hector-multiqc into hector-cedar
       753edc2c Merged master into hector-cedar
       45694866 Merged hector-cedar into hector-multiqc
       0585acc4 Merged master into hector-multiqc
       ff9fcfb2 Merged master into hector-cedar
       05daa1db Merged in hector-multiqc (pull request #38)
       eb7c1fad Merged master into hector-cedar

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      3 commits

       b8978977 Added cedar and mp2b ini files for four pipelines: ampliconseq, dnaseq_high_coverage, pacbio_assembly, rnaseq_light
       289f74aa modified cedar ini to support stringtie
       ac01fb0c Commit of my ini file for rnaseq_light in cedar and modifications to the pipeline to allow for bootstraps in kallisto

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      36 commits

       eb482944 update to last version of mugqic_tools 2.2.1
       4a71fcf5 add dnaseq_high_coverage.mp2b.ini file
       af718ac7 correct select_input_file call
       ea4250e9 correct select_input_file call
       5fbb4568 correct select_input_file call
       c1fcfd2a fix haplotype list parrelization
       84001d2c fix haplotype list parrelization
       260a3363 correct select_input_file call
       950c9808 fix haplotype list parrelization
       bd92b74e correct select_input_file call
       6866c66a fix haplotype list parrelization
       a4427e4d fix number of haplotype and bed compatibilty
       32f78b6b fix number of haplotype and bed compatibilty
       c9bc68a0 correct select_input_file call
       9a4a376e remove file format error
       35b3dcb4 correct select_input_file call
       e6ef8819 include both gatk and gatk4 lib import
       0a34bc75 correct bgzip_tabix libnrary call
       068983ea correct flash log input file
       6a2b0701 correct Mutect2 new arguments
       4da201af correct realigned bam name changed in inherited function
       bccc5cb8 add sequence_dictionary_variant proprety in dnaseq
       dbcd5a43 correct typo in DNAseq base ini
       12dbf3ed remove typo in deilverable command
       dc11a6fb remove verifyBamId for methylSeq pipeline
       637606eb correct typo in the tabix_bgzip call
       26675232 correct a typo
       de9c8b3e add md5sum generation
       d905e598 add required=F for markDup other_option
       908ffc36 change htslib.bgzip_tabix_vcf to htslib.bgzip_tabix as thereis no _vcf function in the lib
       c98c362c make the add UMI steps skip automatically if no UMI in  the readset file
       043ba654 add known_mills option
       40e65362 add required=F for markDup other_option
       51428804 primers addition bug correction
       c859b1a2 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       79ce228f add other_option paramter to compair in order to support other reference than GRCh37

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      2 commits

       c084a8ec Merged in slurm_report (pull request #72)
       4258a453 Merged in new_style_class (pull request #63)

  P-O Quirion <pioliqui@gmail.com>      6 commits

       899d587d parse exit status, and only 50 lines after pro/epilogue
       a50700c8 log error when output logs have the wrong job number
       6b28d34e log report for slurm
       a6c55532 id setter for job
       01c96cfd All class are new style, add setter to some getter
       8bf27ae9 working on Job

  Rola Dali <rola.dali@mail.mcgill.ca>      3 commits

       b72b06b3 Merged in beluga_inis (pull request #74)
       6077a588 adding draft beluga inis
       258876e9 README.md edited online with Bitbucket

3.1.3        Tue Dec 18 16:37:25 2018 -0500        108 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      52 commits

       e239ca0f Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       6cd0351c Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       43c5f25b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1b3f70c5 Software installation scripts : new ones (mainly for Tumor SV pipeline) and updated ones Genome installation scripts : updated Homo_sapiens.GRCh37.sh (uncommented some useful commands) & Homo_sapiens.hg19.sh (updated vcf annotation process)
       8d3a636b Merge branch 'master' of bitbucket.org:mugqic/genpipes
       0a312f31 updated version of mugqic_tools
       395c7dcd Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       1506eac4 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       84fa5170 MethylSeq - DMR - adusted 'other+options' for dmr step
       e6a306d6 new addings to AMP_Scanner installation script
       5e5094f1 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       fd32bb54 MethylSeq pipeline - updated base.ini file with MethylKit differential analysis step requirements
       26d842c6 MethylSeq pipeline - added MethylKit differential analysis steps & cleaned the code
       e25b313b added methylkit.py as a bfx wrapper to call methylKit.R from mugqic_tools
       5869184d debugged 'filter_snp_cpg' & 'prepare_methylkit' subroutines used by MethylSeq pipeline
       0c798db1 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       f6047c44 updated VarScan installation script with lates VarScan version v2.4.3
       efa28c88 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       869f99cf updated Conpair installation script to latest Conpair version v0.2
       7e294213 added Manta installation script
       2233a27f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ad29c150 added Ruby, LAST & Picky installation scripts
       f085968b Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       0683e046 updated AMP_Scanner.sh installation script
       498f1903 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       41a042ae add AMP_Scanner installation script
       2de2537e updated regexp in install_module.sh
       934fdec2 update supernova.sh with latest version 2.1.1
       6a82f17f update cellranger.sh with latest version 3.0.0
       82783dbc Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       36e6b1c3 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       bad57925 CCCG-1143 - added Caenorhabditis_elegans.ce11 genome installation script
       70338383 BFXDEV-578 - updated FasTree installatino script
       c35f7b27 BFXDEV-756 - changed the '-j' pipeline parameter default value to 'slurm' instead of 'pbs'
       4847e413 added StringTie installation script
       6cc4a37a added RNA-Seq Light section in the README
       90ea2232 Version bump to 3.1.3-beta
       0b4dcefa Version bump to 3.1.2
       2e259ec6 MethylSeq - updating ini files after testings on guillimin
       19568c7c Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       cbefa417 MethylSeq - cleaned methylkit_dmr subroutine
       8be6e376 MethylSeq - added MethylKit DMR steps
       2ef5d07e MethylSeq - adding the bfx wrappers for DMR steps
       35aac401 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       e6d9edf0 Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylseq
       685e8a51 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e7f4335d fixing ini for trimmomatic16S and updating mammouth ini for dada2 protocol steps
       b716290b AmpliconSeq : resolving asva step dependencies
       57f08b85 Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       e8af44bd AmpliconSeq - updates to Dada2 protocol
       523d1da6 AmpliconSeq - updated last step of dada2 protocol for better dependencies and coding standards
       8a9605ed BFXDEV-674 - updates : flagstats calls are now made after alignment and deduplication instead of during the metrics step

  Édouard Henrion <henrione@cedar1.cedar.computecanada.ca>      2 commits

       7670b7e1 ChIP-Seq pipeline - correcting cedar.ini file for missing 'homer_make_ucsc_file' step requirements
       b66e9354 AmpliconSeq - call 'zless' instead of 'less' to avoid issues on Graham and Cedar systems

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      3 commits

       fcafa26a Merge branch 'dada2' of bitbucket.org:mugqic/genpipes into dada2
       fa5e298b Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       a92edf35 AmpliconSeq - adding cedar ini file

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      6 commits

       d6515ace Merge branch 'master' of bitbucket.org:mugqic/genpipes into dada2
       c3180f27 update database references in ampliconseq.base.ini
       37685bd4 asva.py pipeline merged to ampliconseq as another protocol asva.R deleted because moved to mugqic_tools new trimmomatic16S function added in bfx/trimmomatic.py
       51810508 update nanuq2mugqic_pipelines.py to also fetch the primer sequences ; needed by the new AmpliconSeq protocol (dada2)
       06948999 removed cutPrimer and set sys.path correctly
       dfe15133 Creating dada2 branch content

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      4 commits

       c54fe2b8 added pool parameter to bfx/dada2.py & pipelines/ampliconseq/ampliconseq.base.ini
       94f9800f removing deprecated pipelines/ASVA/asva.py
       9b09c664 first working verson of dada2 protocol for ampliconseq pipeline
       5acf5272 minor updates (line break & indentation)

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       5627b8e6 GATK4 Updates - gatk4.ini - removed mugqic_dev module
       0f32dbfe Updates GATK4 - dnaseq.cedar.ini - removed mugqic_dev module
       75b91531 Updatres GATK4 - dnaseq.base.ini - remove mugqic_dev module
       dfda31eb MethylSeq UMI - methylseq.base.ini - updated mugqic_tools to mugqic/mugqic_tools/2.2.0 for fgbio tools
       7c9c703a MethylSeq UMI - methylseq.base.ini - updated module_fgbio to cvmfs version of the module
       b29b1171 methylseq.py - commented subroutine all_sample_metrics_report as it has been remove from the pipeline (because useless)
       b8f3fb62 fgbio.py - correted typo in addumi surbroutine

  Emmanuel Gonzalez <emmanuel.gonzalez@mcgill.ca>      1 commits

       223052d2 Merged in dada2 (pull request #49)

  Francois Lefebvre <francois.lefebvre@mcgill.ca>      2 commits

       e1c7ef87 README.md edited online with Bitbucket
       6d1bd5cc README.md edited online with Bitbucket

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a6fe363b remove all_sample_metrics_report as it does the same as ihec_metrics
       b446c866 change the ini section for module in tools.methylseq_ihec_metrics_report
       598028bc Merge branch 'master' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       cee04644 adding umi and on-target metrics
       96f4decd adding the missing count parameter -c to samtools.count
       ba6545f3 correct typos
       0d03902b correct small bugs for UMI integration
       317fe6a3 Merge branch 'methylSeq_UMI' of bitbucket.org:mugqic/genpipes into methylSeq_UMI
       acb55363 corect typo

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      11 commits

       1d2a7d32 add on target metrics and UMI
       86cea706 add target_cpg_profile function in metrics lib
       ad2dd57e add samtools mapped_count function
       b322c94e add samtools mapped_count function
       9d49c2a6 correst typo
       6ccf4ee4 integrate UMI step in ini file
       cf4c008f include bam UMI annotation step
       f74f1ae8 include UMI annotated bam to the select input file of the picard_merge_sam_files step
       5c52b3dd Add other_option  system to MarkDuplicate
       2efad220 add the UMI field in the readset file
       9ef86106 add the UMI field in the readset file

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       7f369e1c Merged in methylSeq_UMI (pull request #52)

  Robert Eveleigh <eveleigh@cedar1.cedar.computecanada.ca>      1 commits

       cf47de1d updates to cedar ini

  Robert Eveleigh <eveleigh@cedar5.cedar.computecanada.ca>      2 commits

       5b912aa9 Cedar resource and gatk4 metric fixes
       672ad330 Updates to cedar.ini

  Robert Eveleigh <eveleigh@ip16.m>      2 commits

       41778221 GATK4 mp2b file added and improvements to alt contig exclusions
       91fb3496 Added mp2b ini and exome improvements

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       e9c66c64 Updates GATK4 and annotations

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      2 commits

       3c0ac311 Resolve conflicts
       71b74961 improvements to exome handling

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      1 commits

       ad3caf1e Merged in dnaseq_gatk4 (pull request #47)

  Rola Dali <rola.dali@mail.mcgill.ca>      1 commits

       dc7996dd README.md edited online with Bitbucket

3.1.2        Wed Nov 21 15:05:01 2018 -0500        30 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      17 commits

       be396f04 updated mugqic_tools version 2.1.12
       f187cdf1 removed _dev modules from ini files
       2bfc3077 Updated R_Bioconductor.sh : pointing to system libraries improved
       e80f7296 PacBio Assembly pipeline - corrected incompatibility bug with Guillimin
       1be0124f updated install_module.sh to accomodate both CentOS7 & Ubuntu16.04 system libraries
       ff709f1f update HiCUP install script with latest version 0.7.0
       659fbc95 corrected typo in ampliconseq.py causing pipeline crach at plot heatmap
       c8b5634e updated README-release.txt with correct URL for GenPipes download page
       09fe3697 Changing resources requirements for [gatk_merge_and_call_individual_gvcfs] in DNA=-Seq pipeline
       30a9d33a corrected reference to kallisto tx2gene file in the RNA-Seq light pipeline
       0e7dae6f corrected reference to kallisto index in the RNA-Seq light pipeline
       e8ec0ac4 Corrected bug in HicSeq pipeline when running in capture mode
       58841e5c updated RNASeq light base ini file with kallisto version on CVMFS
       78a4d4b3 corrected core/pipeline.py to avoid pipeline erroring when using '--json' parameter...
       1e2ef94d corrected typo in R_Bioconductor.sh
       e4f493fe Version bump to 3.1.2-beta
       c6b48be6 Version bump to 3.1.1

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      4 commits

       929e37c2 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f4479be0 freebayes install script
       9fadba98 platypus install script
       7909a4e3 vcfanno fgbio & delly install scripts

  ehenrion <edouard.henrion@mcgill.ca>      7 commits

       362f1016 Merged in hicup_arima (pull request #51)
       f3593a67 tumor_pair.base.ini : changed gemini version to 0.20.1
       31d2d00e tumor_pair.guillimin.ini : removed some more _dev modules
       b689763a tumor_pair.guillimin.ini : removed mugqic_dev module
       6e9b3f99 hicseq.base.ini : updated HiCUP version to 0.7.0
       07b9a21e dnaseq_high_coverage.base.ini : setting ram parameter witin section igvtools_compute_tdf
       d2108932 Tumor_pair pipeline : bug fix in bfx/bcbio_variation_recall.py Correted typo in the executable call

  Rola Dali <rola.dali@mail.mcgill.ca>      2 commits

       6d944741 editing hicseq.py for Arima compatibility
       33d597d3 adding HiC Arima digest to install_genome

3.1.1        Thu Nov 1 15:32:25 2018 -0400        161 commits

  David Bujold <david.bujold@mail.mcgill.ca>      1 commits

       b0adf94d Merged in pipeline_stats (pull request #20)

  dbujold <david.bujold@mail.mcgill.ca>      2 commits

       86c8a724 Display JSON log statistics into tables and figures on the log VM.
       6f716cb9 Python CGI script to tranform pipelines stats log file into a JSON document.

  Edouard Henrion <edouard.henrion@mcgill.ca>      32 commits

       24957210 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b14189df update install_module.sh script with integration of apt along yum as system libraries resources
       0de17ddf update R installation script with integration of apt along yum as system libraries resources
       7a57dfe7 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       de73f87e update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       cd6f62af update dnaseq.py, dnaseq_high_coverage.py & methylseq.py so that interval_list file is created locally instead of where the bed file is (which leads to an error non read-only systems)
       052e9fcb Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e8af086e updating genome ini file generation with versioning
       e2f53d80 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b9297421 update smrtlink.sh with latest SMRTLink version
       8130d059 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       13080e22 updated R_Bioconductor.sh with new packages and corrections
       63e12b54 new version of multiqc.sh to install latest version
       f499c1b9 aded popoolation2 installation script
       211dc622 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       d58b1c3c updated gatk.sh so it can handle both versions 3 & 4 installation
       4b5565f6 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       c4004ed5 updated R_Bioconductor.sh - new libraries & patches
       38e01cde updated python.sh : make -j12 & configure command
       b06c4108 resources/modules/install_module.sh
       d1fe8bf6 update flash.sh with parallele make -j12
       6eca5f3e watch_portal_folder: removed sample_name in filename
       f0e62220 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1ef165c9 corrected gatk.sh : paths in the modulefile were wrong
       749cd972 updated rnaseq.py : report steps are now included in the analysis JSON file
       8978a592 updated ampiconseq.py : both steps merge_flash_stats & merge_uchime_stats are now inclded in the analysis JSON file
       b4d4b046 updated rnaseq_light.py : samples added to each job including report jobs, reviewed indentations and spacing
       8e75e026 updated rnaseq_denovo_assembly.py : samples added to each job (some were still remaining) including report jobs
       61af5bf2 updated common.py : added samples to the call of rmarkdown.render
       cc974a38 updated bfx rmarkdown wrapper : JSON file generation added
       1d9b8347 Version bump to 3.1.1-beta
       bd721f1d Version bump to 3.1.0

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      1 commits

       b3e45601 removing all the remaining MUGQIC_INSTALL_HOME_DEV in all the ini files

  Édouard Henrion <henrione@gra-login1.graham.sharcnet>      6 commits

       4a4d0e57 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       ca805a0f updated install_module.sh : new installation procedure integration, adding the patching the C3G executatbles with patchelf : making sure all system libraries are now searched in /cvmfs/soft.mugqic/yum/centos7/1.0
       f771e9bc updated python.sh with some new package : umap-learn
       63a855eb udpated R_Bioconductor.sh with new packages. Also, now integrates the new installation procedures including pathcing the C3G executables (i.e. use of patchelf)
       41fcad8c Merge branch 'master' of bitbucket.org:mugqic/genpipes
       4a5958c0 updated archive URL in flash.sh install script

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      7 commits

       c243ec8f updated R_Bioconductor.sh script : set the PAGER variable to /usr/bin/less
       f315771d updated pipeline READMEs : all the steps are now shown independently of the available pipeline protocols
       5fa7bb43 updated R_Bioconductor.sh with some new packages in the install list
       1c468e09 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b24d7f89 updated silva.sh to install latest vesion of silva DB
       04140e3a Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e9fcbff4 added Glycine_max.v2.sh installatino script for Glycine (Soybean) genome installation

  edouard.henrion@mcgill.ca <ehenrion@abacus2.ferrier.genome.mcgill.ca>      16 commits

       bdf93cc7 updated hicup.sh
       795b9e5b updated smrtlink.sh with the version 6.0.0 of SMRTLink
       0a9fb4b7 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e141ec79 updated pipeline cedar ini files to avoid using java from CC software stack
       efbf07ab modified core/pipeline.py to avoid having to set  when not generating the anaylsis JSON file
       256907f3 modified DNA-Seq README with better step descriptions
       76a26cb3 Updated the pipeline workflow diagram download links : path of the full-size picture instead of the resized one
       d0fd8745 Added a link to download the pipeline workflow diagram along with the diagram picture itself
       11cbf4a8 Updated pipeline README.md files, with workflow diagram pictures embeded
       30c53cc5 deleted pipelines/tumor_pair.base.exome.ini
       bd6f0000 moving tumor_pair.base.exome.ini from 'pipelines/' to 'pipelines/tumor_pair/'
       9362d720 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       baf53bae updated ampliconseq pipeline following Emmanuel gonzalez comments
       6a85b8fb added source EnsemblPlants to install_genome.sh script
       690f20b3 added skewer installation script
       59996b55 BFXDEV-591 - updated ampliconseq.base.ini based on Emmanuel Gonzalez feedback

  ehenrion <edouard.henrion@computationalgenomics.ca>      2 commits

       3b3ff687 updated mugqic_pipelines.sh so that it now refers to genpipes repository on bitbucket
       97a4552c corrected typo within gatk.sh installation script

  ehenrion <edouard.henrion@mcgill.ca>      18 commits

       12388828 methylseq.py : corrected dependencies assignment for ihect_sample_metrics_report step
       85bfc185 rnaseq_light.py : typo correction
       2057f1e6 rnaseq_light.py : corrected typo
       91be7b89 rnaseq_light.py : corrected typo in job parameter assignement
       5de9c4e2 dnaseq.cedar.ini removed 'module_java=java/1.8.0_121' from cedar.ini
       b4ce9ef9 methylseq.base.ini : modified [bismark_align] section : maximum insert_size now set to 1000
       146eaa73 methylseq.mammouth.ini : modified bismark align section
       a4eea15a methylseq.cedar.ini : modified bismark_align walltime
       a22490f3 methylseq.base.ini modified [bismark_methyl_call] section within the base.ini
       a5276a3d methylseq.base.ini - modified bismark_align parameters and resources within the base.ini file
       4fd03a10 dnaseq.py edited online with Bitbucket Correct typo at line 724 : "jobs" replaces "obs"
       425c086f dnaseq.cedar.ini : added variant_recalibrator section to defined resources on cedar
       d55362f0 Merged in revert-pr-41 (pull request #46)
       17aa18fd Revert "Can run on HPC with slurm and  containers (pull request #41)"
       60be1aa2 rnaseq.py edited online with Bitbucket added missing job (metrics.wigzip) to the json analysis file
       1c5f56dd README.md edited online with Bitbucket
       6dccf5b4 README.md edited online with Bitbucket
       3b3407be README.md edited online with Bitbucket

  Éloi Mercier <emercier@cedar5.cedar.computecanada.ca>      1 commits

       c6593cf5 in cedar.ini and graham ini of rnaseq, chipseq and dnaseq: change assembly_dir to MUGQIC_INSTALL; in dnaseq.graham.ini: uncomment assembly_dir variable

  emercier <eloi.mercier@mcgill.ca>      10 commits

       ec3183ff all pyc files removed
       f0568f47 Update R_module in ini files to mugqic/R_Bioconductor/3.5.0_3.7 (except for illumina_run_processing)
       5a14b4d8 in ampliconseq, pacbio and rnaseq.guillimin.ini: change lm queue (depreciated) to meta queue
       e13be161 in modules/weblogo.sh: remove whitespaces at the beginning of the echo blocks
       4675ad67 in dnaseq and rnaseq.base.ini: change R_Bioconductor to stable version 3.4.3_3.6
       fcad0258 in install_all_genome.sh: add Danio_rerio.GRCz11.sh
       d353ecf9 add install script for genome zebrafish Danio_rerio.GRCz11
       94161565 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0b53f645 in install_genome.sh: small fix to create_kallisto_index and create_transcripts2genes_file functions
       f5681c4d in install_genome.sh: fix bug in create_transcripts2genes_file to work with recent version of Ensembl

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      4 commits

       53d11f3e rnaseq.base.ini edited online with Bitbucket: commented explicit adapter file parameter
       8653d027 R_Bioconductor.sh edited online with Bitbucket: added PopSV. Currently commented out since not tested
       0b42f96c R_Bioconductor.sh edited online with Bitbucket: Added a few dependencies to the list
       3a268e66 R_Bioconductor.sh edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      1 commits

       a23ba480 Added dev install scripts for delly, lumpy, sv, vcfanno

  José Héctor Gálvez López <hgalvez@ip16.m>      2 commits

       88e49293 Further refinements to Mp2b based on feedback from mammouth admins
       39bf7606 Added ini files tailored for Mp2b based on Cedar ini files for the following pipelines: RNA-seq, RNA-seq de novo, DNA-seq, Hi-C seq, Methylseq, and ChIP-seq.

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       1f95d680 Add cedar ini file for methylSeq
       8b071760 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       b1ed77cd removing pipe empty of jobs in tumor_pair
       7e70dc57 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       db209577 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       9e0178b2 correct bug in ucsc bedGraphToBigWig which raise an error for specie with MT chromosome name except for GRCh37 - BFXDEV-737

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      7 commits

       2e021f67 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       85ead2b6 resolve issue with multiple inputs in DNAseq picard markduplicates
       1e00830b resolve pull conflict
       b0042a0e update dnaseq.cedar.in file
       4074ab5f adjust dnaseqto add select input file to some of the steps - need to be continued
       e38920e4 modify Slurm scheduler delay (sleep) from 0.5 to 0.2
       b1d3d4b5 ChipSeq - add mutliqc param in the cedar ini

  Pierre-Olivier Quirion <pierre-olivier.quirion@computationalgenomics.ca>      1 commits

       9c7c8158 Merged in add_container (pull request #41)

  P-O Quirion <pioliqui@gmail.com>      5 commits

       cfa63b20 Merge branch 'master' into add_container
       44274bdb Basic container ini file
       4296d317 Working singularity version tested on graham
       b14afd93 WIP exec line prototype
       8b924161 add container option

  Rola Dali <rola.dali@mail.mcgill.ca>      8 commits

       795776b9 methylseq.cedar.ini edited online with Bitbucket: added cluster_walltime to ini to avoid errors
       47addd28 chipseq.base.ini edited online with Bitbucket: bigwig and run_spp resources are not enough; edited them to avoid failure
       50747a73 Merged in IHEC_metrics (pull request #44)
       bdb0b76e dnaseq.base.ini edited online with Bitbucket: change threads for GATK due to errors
       142921b0 update csvToreadset.R in utils to use column names due to changes in nanuq csv
       40cdff0f Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       8709b9d8 README.md edited online with Bitbucket
       13f08a5f README.md edited online with Bitbucket: -j slurm and tutorial

  Romain Grégoire <romgrk.cc@gmail.com>      3 commits

       6f8ea86f Merged in add-json-version (pull request #50)
       1271bc19 Merged in dashboard-display-launching-user (pull request #43)
       73e768a4 Merged in logs-add-checksum (pull request #42)

  Rom Grk <romgrk.cc@gmail.com>      28 commits

       ff11dbc5 jsonator.py: add version number
       bb72124f watch_portal_folder: read sample_name from filename but stay backward-compatible
       921e54df job2json: add sample_name to file for genpipes dashboard
       784848f6 copy sample json files during script run
       b414b045 scheduler.py: clean unused value
       5a9d8afc Sample: remove .json_dump property
       7d45ce1e lint pipeline.py
       337a1c8f lint job2json
       014c19d1 job2json: guard __main__
       26b9ef98 lint job2json
       bf4f976d job2json: use $USER of user running script
       8980701e common.py: remove unused import
       622b6111 common.py: fix log issues
       676f7bee common.py: add unique md5 checksum to logs
       56545fec Revert "watch_portal_folder: put sample_name in filename"
       a8badcd0 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0526ecca watch_portal_folder: put sample_name in filename
       13b56384 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       0411ae98 watch_portal_folder: fix memory errors
       e5e35795 watch_portal_folder: add logging
       46d8a35d watch_portal_folder: fix typo
       f54a1277 watch_portal_folder: add cache option
       7f89e828 Merge branch 'master' of https://bitbucket.org/mugqic/genpipes
       2fbdca88 watch_portal_folder: implement update by diff
       45cc1899 watch_portal_folder: update script
       67c14714 watch_portal_folder: skip missing files
       8a9dd08a watch_portal_folder: more resilient to network errors
       32bf92a9 watch_portal_folder.py: dont watch if no interval is provided

3.1.0        Wed Mar 28 15:46:33 2018 -0400        188 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      77 commits

       ff52e7c9 MethylSeq - added 'ram' parameter to the igvtools_compute_tdf step in methylseq.base.ini
       68661e8c updated jb2json.py with a better locking system : now creates a folder instead of a file in order to create the lock
       55bf1b0f Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       43c9cc1d updated core/scheduler.py to remove job2json call when --json parameter is omitted
       21b3bdfd updated gatk.sh with version gatk-4.0.2.1
       99087312 updated longranger.sh with LongRanger version 2.2.2
       d03abd52 updated hicseq.py : calling of the newly built bfx libraries and minor indentations changes
       c52b5820 reviewed topdom.py wrapper for a better handling of input and putput files, so that job dependencies do not break...
       8ae89bdf aded locking file system to jsonator.py and job2json.py to avoid multiple and synchronous writing attempts to JSON file
       4d1eaa68 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       70ddf870 updated hicseq.base.ini with newer version of mugqic_tools (2.1.11) and with revised resource allocation for hic_align step
       3f798270 removed useless loading of mugqic_tools module when calling job2json.py
       4d7da457 added locking file system to job2json to avoid multiple and synchronous writing attempts to the JSON file
       6aafbf04 removed useless loading of mugqic_tools module when calling job2json.py
       205d577d mugqic_tools.sh : swith to mugqic_tools version 2.1.11
       0ea6336d improved the process regarding the update of the analysis JSON file when resuming a pipeline execution
       62dd71d3 removed useless import from core/scheduler.py
       ac6e9ce6 added the locking file system to avoid multiple & simultaneous writing attempts on the same file
       ffd35a6d Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       e6545082 added new bash tools wrappers to be called by hicseq.py pipeline script + some minor indentation & syntaxe updates
       843d6573 added new bfx libraries to be called by hicseq.py pipeline script
       1f6be1b7 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       e1dd68f7 added some new software installation scripts
       e2abf45c updated 'core/job.py' : added 'samples' parameter to 'concat_jobs' and 'pipe_jobs'
       8cf84f4e added the initialisation of 'portal_output_dir' to '/lb/project/mugqic/analyste_dev/portal_output_dir' within all of the pipeline .base.ini files
       9bae17c1 added the generation of the analysis JSON file to the SLURM scheduler class
       6e6363d6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into cedar
       87fb0142 updated core/scheduler.py : call of job2json.py is now done without /home/ehenrion/work/portal_repo/genpipes
       1681c6f6 Merge branch 'master' of bitbucket.org:mugqic/genpipes into IHEC_metrics
       f95e6909 corrected trimmomatic input selection to reflect the file names what might be outputed from picard_sam_to_fastq
       f22ee01e update verifyBamID subroutine within common.py so that it now accepts .vcf.gz files (formerly was only accepting .vcf files)
       998e1bf4 updated version of cellranger to 2.1.1 within the bash installation script
       ba5099f4 updated pipeline READMEs
       b5cecbd0 updated core/pipeline.py to make the analysis JSON file creation optional : no JSON file created by default
       3fae6e86 updated bfx/macs2.py so that it follows the C3G coding standards
       3b426eaf added installation script for hdf5 and zlib libraries
       2e82f079 corrected typo in install_genome.sh
       ff59b0a4 updated install_genome.sh script with newer software version as well as minor indentation fixes
       3d26798f Merge branch 'master' of bitbucket.org:mugqic/genpipes
       1bd0b988 updated kallisto installation script
       cf3219ad updated STAR version in the installation script
       19997faf updated dbSNP version to 150 for Homo_sapiens GRCh37 installation script
       1e6f0320 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f89e056f updated kallisto version to 0.44.0 in the genome installation script - corrected typo within create_transcripts2genes_file subroutine
       7bf13a92 updated dbSNP and dbNSFP in Homo_sapiens genome installation scipts
       30299e2b added cDNA fasta to the genome installation scripts
       e17fcb1b added the portal_ouput_dir to in base ini file
       acea9bf9 Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       fe616275 Merge branch 'portal-integration' of bitbucket.org:mugqic/genpipes into portal-integration
       e0010e28 Merge branch 'master' of bitbucket.org:mugqic/genpipes into portal-integration
       1f63dc5d Added README-GenAP_coding_standards.txt which gives the coding guidelines for who may want to participate in the C3G/GenAP developments
       ed345d28 BFXDEV-721 - updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       9a59536b updated install_genome.sh script : corrected create_transcripts2genes_file() subroutine with missing “then” and “fi” within the 'if/else' statement
       0ed0ef60 updated Prokka installation script so that 'prokka --setupdb' is launched afer the installation
       743668f1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       22b26013 updated Supernova bash installation script to version 2.0.0
       c62cf42e updated the bash installation script of python : installation of QIIME has been updated so that emperor is well installed before QIIME
       96c98a44 added the initialization of the LFS variable, otherwise it is recognize on guillimin...
       b02491e0 updated the bash installation script of python : after the installation of the specified version of Python is successfull, the script now also takes care of the installation of all the needed python libraries
       7609147a update type within python_lib.sh script
       1ab7a7ea RNASeq pipeline : value of parameter 'localfit' used by DESeq, is set to default i.e. empty (which also means 'false'), so Parametric Dispertion Fit is performed by default instead of Local Fit
       896c3b2b DNASeq pipeline : updated versions of mugqic_tools (to 2.1.10) and of samtools (to 1.4.1) within the base.ini file
       74a5c9ad updated picard installation script : now installs version 2.17.3
       8d2dde0c added the installation script of Prokka, a tool for prokaryotic genome annotation
       cbfe9d75 BFXDEV-673 - remove some verbose during the execution of job2json.py
       66893f89 updated mugqic_tools.sh so it now installs version 2.1.10
       170f39bf added the bash installation script for the Illumina InterOp parser
       dc0ebcb0 BFXDEV-674 - MethylSeq pipeline - minor updates within the metrics .md report template file, for standardization purpose
       d9e7b0b4 updated ortograph within sample metrics .md report files
       e0b41ff9 BFXDEV-673 - updated scheduler.py to generalize the use of job2json to all schedulers
       5250d807 BFXDEV-674 - added a report .md file for IHEC metrics reports for targeted anaylsis
       cb466e70 BFXDEV-674 - updated .md files for metrics reporting : revised headers & descriptions
       460f24a6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       2833f6a9 updated gemini.sh bash installation script : added PYTHONHOME setting when loading the gemini modules
       f483e37e update release instructions to generate proper README for all the pipelines including HicSeq
       2d70be4d Version bump to 3.0.1-beta
       3cb8610e Version bump to 3.0.0 - updated

  Édouard Henrion <henrione@cedar5.cedar.computecanada.ca>      2 commits

       89831808 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       49c6b320 updated kallisto bash installation script

  ehenrion <edouard.henrion@computationalgenomics.ca>      3 commits

       3f60795b updated R_Bioconductor.sh with good indentation and new packages installation
       b1215c32 corrected wrong environment variable name within blast.sh
       dc4fbb7b updated the version to v359 within ucsc.sh script

  ehenrion <edouard.henrion@mcgill.ca>      6 commits

       39fea264 rnaseq_denovo_assembly.base.ini edited online with Bitbucket
       f9278f90 README.md edited : "MUGQIC pipelines" replaced by "GenPipes"
       e6abfd90 apis_mellifera.sh deleted : was the exact replicate of Apis_mellifera.sh
       8e9ab1cd README.md edited Updated links to the pipeline pages
       469b9b48 README.md edited online with Bitbucket updated some links to reflect the repository renaming to genpipes
       87524a9d smrtanalysis.py - standardized command-line format and indentation

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      3 commits

       99b1a5dc removing pyc files from rnaseq
       1c7731a5 in rnaseq.mammouth.ini add section for bed_graph to set up ppn to 1
       6d910887 In nanuq2mugqic: change readset file name to readset_<seq_type>.tsv

  Eric Fournier <ericfournier2@yahoo.ca>      2 commits

       6fafc8d0 Merged in ericfournier2/genpipes (pull request #35)
       4138a728 Fix list within list bug which breaks chipseq pipeline.

  Gary Leveque <gary.leveque@gmail.com>      1 commits

       68573076 Merged in gary_pacbio_assembly (pull request #29)

  gary.leveque@mail.mcgill.ca <gleveque@abacus2.ferrier.genome.mcgill.ca>      7 commits

       fd1e91cc additions made to smrtanalysis.py for basemodification and motifMaker steps
       b88dd9e8 Addressed issues commented by Edouard; tested on abacus and mammouth
       39f57bce added cluster_server= to pacbio_assembly.mammouth.ini
       b135ae0b revised versions of .base and mammouth.ini files; I was changing between sw and lm nodes
       6bac3683 revision of pacbio_assembly.mammouth.ini, back to qwork
       356aadc0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into gary_pacbio_assembly
       a24451f2 Addition of base modification detection and generation of a motif_summary.csv steps to the pacbio HGAP assembly pipeline;  see BFXDEV-703

  Jose Hector Galvez <jose.hector.galvez@computationalgenomics.ca>      1 commits

       2455b6db Merged in hector-cedar (pull request #37)

  José Héctor Gálvez López <hgalvez@cedar5.cedar.computecanada.ca>      1 commits

       36b892fa Corrected minor bug in the create_scheduler() function that was creating errors when using slurm

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       05bbcec5 Merge branch 'master' of bitbucket.org:mugqic/genpipes
       3b818dc9 externalize excluded chormosome to be specified in the genome ini
       a2d807b9 remove conflict
       3896d772 uniformize verify_bam_id vcf path
       f4bf8cc7 add a delay (sleep) after the symlink creation to avoid issue (Invalid Job Dependency) during submission
       322348df Merge branch 'master' of bitbucket.org:mugqic/genpipes
       f734dfb4 test new url
       37529400 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0834b95 add NovaSeq to nanuq2mugqic script
       0dd2ba3e DNAseq - split the pipeline into 2 different pipeline GATK best_practice or old mpileup - BFXDEV-716
       5c586fff Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9d48a589 chnage RNAseq_denovo inheritance to RNAseq
       8a250c81 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7716f961 include a test of inputs type for some library function which does not enforce the list type

  Mathieu Bourgey <mathieu.bourgey@mcgill.ca>      24 commits

       e2b2a7ab resolve conflict with master
       82d1bee6 Merge branch 'cedar' of bitbucket.org:mugqic/genpipes into cedar
       443fcf64 restore low memory ini now the slurm bug is corrected
       78d156b3 update rnaseq for cedar
       62520dad update rnaseqDN on cedar
       9a274023 update hicseq on cedar
       b9842b9d update dnaseq ini on cedar
       b695bcfe add more more RAM to compensate slurm bug
       924fff54 add 0.5s sleep to let slurm submiting the job correctly
       282fa9eb Increase memory request to avoid I/O buffering hit the memory limit
       5f157cd3 Changing I/O block size
       3b395c5e Changing I/O block size
       614564d5 remove confilct pulling master
       be815794 Modifying I/O block size for cedar
       496e10d2 update ini file
       7f77b195 update log with slurm scheduler
       81806098 add chipseq ini file for cedar
       0e279092 update RNAseq and link to mugqic_dev
       d26fef63 update DNAseq and link to mugqic_dev
       10c20693 remove lattency for jobsubimission to slurm
       d192a915 remove conflict
       77e5d359 update cedar scheduler
       67b80063 update cedar ini
       120f92d3 DNAseq -ini file for cedar

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       804ae6cf Merged in cedar (pull request #34)

  Mathieu Bourgey <mbourgey@cedar5.cedar.computecanada.ca>      9 commits

       9a9ab773 DNAseq - update ini file
       910be769 remove scheduler module testing
       9e416b18 change igvtools exe to igvtools jar in order to have control of the ram usage
       9d3c7bdc adjust RNA and DNA ini files;  generate fake prologue and epilogue; add 2s delay between each job submission- BFXDEV-683
       95f20c63 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       f247a6a9 create RNAseq ini file - BFXDEV-683
       6f5cca53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into cedar
       297a2c6b create a DNAseq ini for cedar; remove typo in base DNaseq ini  - BFXDEV-683
       1dce04b8 create a version of the scheuler class for slurm  - BFXDEV-683

  Mathieu Bourgey <mbourgey@gra-login4.graham.sharcnet>      1 commits

       fea15010 Add Graham ini for RNA and DNA

  Rola Dali <rola.dali@mail.mcgill.ca>      17 commits

       6978e141 Merged in IHEC_metrics (pull request #33)
       d66737bc added .hic file generation to capture hic
       3a9b9a38 adding RobusTAD TAD scoring to hicseq.py
       3c797210 adding multiqc to chipseq. should customize yaml file and check homer module
       ad0ff99a changing rnaseq dependencies to ensure no repeats when job is complete
       e581c3ab adding mugqicValidator.py to utils to validate basic structure of readset and design file
       8ee10839 removing job outputs from chipseq ihec metric method to accomodate samples running without input
       befab29e Merged in IHEC_metrics (pull request #31)
       2a0d785c fixing the tmp_dir for macs
       6a363de2 Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       a10ad518 changing genome to fit with merged homer.maketagdir and changing ihec matrics output name
       bba10e0b Merged in IHEC_metrics (pull request #28)
       28850c50 config.py edited online with Bitbucket
       00b9e34e adding query_module=spider to mammouth ini
       00a746fa allowing module spider on mammouth to reduce module loading time--committing to test on guillimin and abacus
       e3c9fabb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fea082b7 fixing MACS bug BFXDEV-712

  romain.gregoire@mcgill.ca <rgregoir@abacus2.ferrier.genome.mcgill.ca>      2 commits

       53fe3b91 pass config files directly to job2json.py
       94d21ed8 Merge branch 'portal-integration' of https://bitbucket.org/mugqic/genpipes into portal-integration

  Romain Gregoire <romgrk.cc@gmail.com>      1 commits

       655f7c93 watch_portal_folder.py: fix --username argument

  Romain Grégoire <romgrk.cc@gmail.com>      1 commits

       2bcae74b Merged in portal-integration (pull request #32)

  Rom Grk <romgrk.cc@gmail.com>      15 commits

       c3a00809 watch_portal_folder.py: remove username option
       f3c3f34e send $USER to portal integration
       bdb24f8a watch_portal_folder.py: fix file path
       4981525b watch_portal_folder.py: sort files sent by modification time
       e1f361a8 export CONFIG_FILES to allow loading config from pipeline
       430d0c48 export CONFIG_FILES to allow loading config from pipeline
       446122c6 job2json.py: add mugqic dir to python path
       c35b6f93 watch_portal_folder.py: safer handling of response
       13196bfd watch_portal_folder.py: fix data posting
       bbbfe8a2 watch_portal_folder.py: fix script exit
       bff9c649 watch_portal_folder.py: fix argument passing
       47b1073d use uuid to avoid collisions when buffering JSON files
       7f99bdc2 job2json.py: make a copy of the updated JSON files for the portal
       629cef61 add watch_portal_folder.py
       63e4de62 jsonator.py: make a copy of the JSON files to a buffer folder

3.0.0        Thu Dec 7 14:19:49 2017 -0500        444 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      247 commits

       8c5a5c61 MethylSeq pipeline - BFXDEV-674 - updated wiggle_tracks step with more comprehensive output file names
       2024838e updated ucsc.py with simplified if-else statement for more clarity and corrected the 'chr' prefixing behavior
       bcf3d8a6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       85942881 updated Homo_sapiens.GRCh38.sh installation script with vcf indexes
       cb43d498 updated jellyfish installation script with version 2.2.3
       9e3c955a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7e14f091 Updated README for all the pipelines
       a6b244e4 updated mugqic_pipelines.sh script with the most recent version of the mugqic_pipelines, now called GenAP_Pipes, version 3.0.0
       86092c6d Version bump to 3.0.1-beta
       b8f43102 Version bump to 3.0.0
       3844eee8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c9ad3680 added tagInfo.txt in homer.py
       33702603 DNA-Seq pipeline - update base.ini file by removing all the paths which were still pointing to a '_dev' location
       7b15520f BFXDEV-674 - MethylSeq pipeline - updated bedtools.py intersect function to carry header from input bed to output bed
       da1ece6b added tagInfo.txt
       1b84ec6d Version bump to 3.1.0-beta
       03f3dda3 Version bump to 3.0.0
       5d0a7749 added the README file for RNA-Seq Light Pipeline
       837b4fa9 slightly updated release instructions
       1343ab64 Version bump to 3.0.1-beta
       188037f9 Version bump to 3.1.0-beta
       03638bd3 Version bump to 3.0.0
       4c027b6f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8768a778 updated version to 2.0.2 for cellranger in the bash installation script
       00eb96b9 RNA-Seq pipeline - updated .base.ini file : removed specific module loading for differential expression step + minor indentation updates within trinity.py for RNA-Seq de-novo Assembly pipeline
       a4ab6975 RNA-Seq pipeline - update .base.ini file with newer software versions
       626856f9 DNA-Seq pipeline - updated base.ini with GATK version 3.7
       a7b56ba1 updated base.ini with GATK version 3.7
       dde4f146 master VERSION changed to 3.0.0-beta
       bfdcbae6 PacBio Assembly Pipeline - updated circlator step : created a circlator.py within the bfx folder and review the coding/calling of the circlator step within the pipeline python wrapper
       132f2f24 BFXDEV-674 - updated MethylSeq pipeline to correct wrong dependencies causing some jobs to always be consider as 'NOT up to date'
       8cbbff9d BFXDEV-673 - corrected error in bfx/jsonator.py occuring when modifying the list of softwares from an already existing JSON file
       94280461 some more minor updates on bfx/tools.py regarding standard indentation and parameters naming
       1aaf4f94 MethylSeq - IHEC metric report jobs are now labelled with the sample name
       29e74ac5 minor updates on bfx/tools.py especially to make indentation uniform across the whole file
       573dba59 updated & added software install scripts
       c8fa5cd4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d964b2b commit prior merging to avoid conflict
       2baf9826 AmpliconSeq pipeline - adding mammouth walltime for qiime_otu_assigning step
       fb8c75db resolving conflicts in bfx/tools.py
       fdddddb7 BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq pipeline
       c961e158 corrected typo introduced after resolving conflicts...
       f4116b0d BFXDEV-673 - updated jsonator.py to handle Illumina Run Processing ini file entries when generating JSON analysis file for Illumina Run Processing pipeline
       d573722c BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       fb346c9d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       763aac03 BFXDEV-673 - adding analysis JSON file generation to the Illumina Run Processing pipeline
       8b401290 BFXDEV-673 - adding analysis JSON file generation to the Tumor Pair pipeline
       96532505 BFXDEV-673 - adding analysis JSON file generation to the RNA-Seq De Novo Assembly pipeline
       f961ef10 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1dcae3f7 BFXDEV-673 - added the analysis JSON file generation to the RnaSeq pipeline
       f107e049 BFXDEV-673 - added the analysis JSON file generation to the PacBio Assembly pipeline
       b806ac67 BFXDEV-673 - updated jsonator.py to handle PacBio Assembly readset files when generating JSON analysis file for PacBio Assembly pipeline
       fdd8577d BFXDEV-673 - added JSON analysis file generation to the DnaSeq high Coverage Pipeline
       d48d6424 BFXDEV-673 - adding JSON analysis file generation to DnaSeq pipeline
       b8824eff BFXDEV-673 - adding JSON analysis file generation ito ChipSeq pipeline
       d19b3044 BFXDEV-673 - updated AmpliconSeq pipeline with analysis JSON file generation
       102f66cb BFXDEV-673 - updating jsonator.py to generalized the way dbsnp_version and server entries are handled
       054e4109 BFXDEV-674 - MethylSeq pipeline - adding samtools to the loaded modules for methylseq_metrics_report
       0169e509 BFXDEV-673 - minor update to the help content
       18e8e689 BFXDEV-674 - updated wiggle tracks generation tools by spliting bedgraph and wiggle traks into 2 different jobs, also managing .bw file for GRCh37 build to make it UCSC compatible
       cc7af63b updated 'tmp_dir' in guillimin .ini files : now using  space which is automatically cleaned up after the job ends
       04d9fed3 BFXDEV-674 - updated methylseq_ihec_metrics_report to reflect changes in mugqic_tools IHEC_methylseq_metrics.sh
       6e854c2a Picard mark_duplicate minor update
       8958faef BFXDEV-673 - corrected error in variable assignment within jsonator.py
       48722763 BFXDEV-673 - jasonataor.py now handles the case where user has used the 'force' argument for his pipeline execution in order to entirely re-create the analysis JSON file
       22c477d9 BFXDEV-673 - review job2json.py passed parameters to handle start & end job date/time
       61225eeb BFXDEV-673 - JSON analysis file key 'hpc_center' changed to 'server'
       4b2da7f9 BFXDEV-674 - MethylSeq pipline - updated mugqic_tools from dev to prod within pipelines/methylseq/methylseq.base.ini
       3e8afc37 BFXDEV-673 - updated core/scheduler.py to handle both job_start_date and jobs_end_date for the json analysis file
       c777a647 BFXDEV-673 - updated core/pipeline.py for a better generation of the json analysis file
       1d735a0e BFXDEV-673 - updated bfx/jsonator.py
       1b0b956f BFXDEV-673 - updated utils/job2json.py
       98b0d7f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       625ca4fc added optparse library installation in R_Bioconductor.sh
       f751eb50 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c41c514b removing MULTIQC from mugqic tools and updating version to 2.1.9
       00118d61 version 1.0.5 of R_mugqic_packages
       18ab8513 BFXDEV-674 - removed useless bedtools parameters from mammouth and guillimin ini files
       a17600f9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       b98dcf91 adding some new bash install scripts
       89a48fdf add chicago module
       0c6d44a0 fix conflicts
       9919d1a7 keep locale module change on libuser version of the repo
       6a3eff86 keep locale module change on libuser version of the repo
       2b55b3fd BFXDEV-674 - corrected input flagstat file name for metrics steps
       9ac6e813 BFXDEV-674 - adding on_target context handling
       69eb2bb3 BFXDEV-673 - corrected log file assignment to the JSON file
       7dd7c0de BFXDEV-673 - BFXDEV-674 - updated syntax and added IHEC metrics report
       5fccb8b4 importing samtools from common.py
       18e87a58 BFXDEV-668 - added ihec_metrics step to rnaseq.py
       167e9c30 BFXDEV-674 - removed useless parameter settings from ini file
       f7e2c946 BFXDEV-673 - updated cluster_hpc_center parameter within the ini files
       5ca6c521 removed unnecessary trailing comma dnaseq.py
       62c711a8 removed unnecessary parameters from dnaseq.mammouth.ini
       e3186d4e BFXDEV-674 - corrected call of the bash script which generate the metrisc for IHEC
       10c042dd BFXDEV-673 - updated jsonator for a better behavior of json file updates
       bfc0042d updated syntax for unrequired paramaters in bedtools.py
       8e42663c BFXDEV-674 - updated bedtools with proper syntax standards
       f9538e80 BFXDEV-668 - BFXDEV-675 - updated bedtools graph other_options parameter to fit with group syntax standards
       790234ed BFXDEV-668 - BFXDEV-675 - corrected typo in ucsc.py and bedtools.py
       af169601 BFXDEV-674 - adjusted walltime for bissnp step
       ba015ef8 BFXDEV-673 - corrected error in job2json.py
       3f5c34e5 BFXDEV-674 - corrected typo in ucsc.py
       47cd3f03 BFXDEV-674 - updated guillimin.ini
       fb5d03cf BFXDEV674 - updated base ini
       91006662 unset batch.ini file
       4b8594c0 BFXDEV-674 - briaree ini file
       14ae3c1c BFXDEV-674 - mammouth ini file
       58f6b280 BFXDEV-673 - updated some 'if' statements to avoid syntax errors...
       f036046f BFXDEV-674 - updated call of ucsc.bedGraphToBigWig within the pipeline wrapper
       8e187ab8 BFXDEV-674 - added walltime to bismark align step
       a2bcb0c2 BFXDEV-674 - updated bedGraphToBigWig subroutine in ucsc.py, added bedToBigBed to ucsc.py, and updated minor things in bedtools.py
       05065478 BFXDEV-668 - BFXDEV-675 - updated bedtools.py graph subroutine which nows calls ucsc.py to avoid code redundancy of bedGraphToBigWig
       4a6e9cdb BFXDEV-668 - BFXDEV-675 - updated ucsc.py to handle cases where bedGraph is in .gz format
       895d8202 BFXDEV-674 - type correction in job2json.py
       a3227054 BFXDEV-674 - cancel the creation of one folder per sample for the json files as on file per sample is enough
       12a085c0 BFXDEV-674 - updated json file path in scheduler.py
       c2f2827e BFXDEV-674 - updated bismark align & dedup outputs in order to create good dependencies for all_sample_metrics_report and ihec_sample_metrics_report steps
       607ee222 BFXDEV-674 - changed the json file location to a more simple one and got rid of the resume subroutine since not used anymore
       5b6f982f BFXDEV-674 - corrected error in scheduler when launching job2json command
       25a4eded BFXDEV-674 - another typo correted in core/scheduler.py and removed the dev references from pipelines/methylseq/methylseq.base.ini
       d5aa74d3 BFXDEV-674 - correcting typo in core/scheduler.py during Json generation
       1dd7bd43 BFXDEV-674 - added some missig report files and updated job2json file command and tool
       c8266f3e Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       da1de3cd BFXDEV-674 - updated pipeline wrapper with report files, ihec metrics and updated metrics computing
       c59a3025 BFXDEV-674 - added HPC center name to the ini files
       a128e704 BFXDEV-674 - updated bedtools when working on MethylSeq pipeline to a better handling of job piping from bamtobed to coverage
       57b02291 MethylSeq pipeline - updated pipeline python wrapper - BFXDEV-674
       807a26bd BFXDEV-673 - add a python script in utils to handle the appending of information to the JSON file
       31a092d7 BFXDEV-673 - add the bfx python wrapper to take care of the creation of the JSON file
       4c575e94 BFXDEV-673 - added a sample list object to the Job class to handle JSON file generation
       c9e342e0 ucsc bfx tools - bedgraph_to_bigbwig : removal of the temporary sorted bedgraph after bigWig is created
       f7f01137 MethylSeq pipeline - updated version of the pipeline : steps have been condensed, json file for each sample has been added, more metrics have been added to fit with HiHeq requirements, whole genome datasets are now handle correctly
       be4179ed MethylSeq pipeline - updated guillimin ini file for bismark_dedup and bissnp requested resources
       3ee60cd0 MethylSeq pipeline - very minor changes within the base ini file
       34cddb6f DNA-Seq pipeline - added the verifyBamID step to the pipeline (BFXDEV-619) and json file generation add-on
       a638e4cb DNA-Seq pipeline - added verifyBamID settings to the base ini file - BFXDEV-619
       52fb55a2 MethylSeq pipeline - updated common.py so that all the common functions (i.e. sam_to_fastq, timmomatic, merge_trimmomatic_stats, verify_bam_id) now generate a json dump to be append to each sample json file
       1425d42c MethylSeq pipeline - updated scheduler.py to handle the json file generation while the pipeline is running, i.e. adds a json section sample json files as soon as a job successfully ends
       5dae0403 MethylSeq pipeline - updated pipeline.py to take care of the creation of the json file for each sample
       8a915230 MethylSeq pipeline - modified sample.py to handle json file creation during pipeline execution
       16c8f934 MethylSeq pipeline - minor update of the verifyBamID python wrapper
       c330a378 MethylSeq pipeline - reviewed picard add_read_group command
       585c92ba MethylSeq pipeline - reviewed parameters passed to bismark align
       fc34adf7 MethylSeq pipeline - removed jsonator import from bfx, waiting for it to be fully implemented
       3674baf7 MethylSeq pipeline - updated the pipeline with new metrics calculations as well as merged some steps together
       b579365e BFXDEV-619 - added verifyBamID .Rmd template file for report
       91071fe0 MethylSeq pipeline - added a python wrapper for all the tools related to methylation profiling
       9a2defcb BFXDEV-619 - added verifyBamID in pipelines/common.py
       a3b94d33 MethylSeq pipeline - added GC Bias R tools to metrics
       55133583 BFXDEV-619 - added verifyBamID in bfx/tools.py
       168fb692 MethylSeq pipeline - Added ucsc.py, with 'bedgraph_to_bigbwig' which is now called by new bedtools.graph function. Also called from methylseq.py
       c3ce4497 MethylSeq pipeline - added bedtools 'coverage' & 'bamtobed' to bfx/bedtools.py, also updated 'graph' & 'intersect'
       946bd708 BFXDEV-661 - DNAseq - added the use module R_Bioconductor within picard_collect_multiple_metrics
       a5b24ac4 BFXDEV-642 - RNAseq - passed the adjusted p-value as a markdown parameter to the differential expression report
       6e5470ad BFXDEV-644 - RNAseq - added the loading of module_python and the use of 'bash pipefail' for step differential_expression_goseq_report
       38db2193 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d8a33d9a execute permissions updated
       1d8f3a77 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       dee964a4 MethylSeq pipeline - create specific ini file to handle capture datasets
       3ce2005f MethylSeq pipeline - create specific ini file to handle capture datasets
       c042cb20 corrected typo in python3.sh install script
       6406b9ce added index file as part of the output files in add_read_groups function, in picard.py
       8a940f31 added index file as part of the output files in add_read_groups function
       1c2f36a6 MethylSeq pipeline - updated wiggle track step with splitted forward and reverse strand reads
       040eefb5 MethylSeq pipeline - corrected typo in bfx/bedtools.py
       a1c2b6c8 new CHANGELOG coming along release v2.3.1
       0f801b87 updated Homo_sapiens.GRCh38.sh install script with Ensembl87, dbSNP_149 and dbNSFPv3.4
       ec648994 updated rnammer_transcriptome resources within rnaseq_denovo_assembly.guillimin.ini
       ae6119c6 updated surpi.sh install script
       ae8ebb52 updated snap.sh install script
       08908b5c updated seqtk.sh install script
       b6d06b85 updated bowtie2.sh install script
       f040ef64 updated RAPSearch2.sh install script
       ad4b019c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       41fdf29a corrected bowtie.sh install script
       e9772774 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       79bb77eb updated versions of some module install scripts and added some other
       3f90a979 updated versions of some module install scripts
       f6ffda18 MethylSeq pipeline - adjuted resource allocation for bismark align step
       53c3d852 resolving conflicts and merging
       e5e00743 minor changes to resolve conflicts before merging mMaster with MethylSeq
       eb215751 Small changes and updates before merging MethylSeq branch to Master to prevent conflicts
       93d9f01f MethylSeq pipeline - updates done after pur pull request review
       e5eb8fa1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       adb1cbeb updated version to v346 in UCSC module installation script
       faa82153 AmpliconSeq pipeline - updated version of R_Bioconductor to 3.3.3_3.4
       23d2966c new genome installation script : Apis_mallifera.sh
       9173ee88 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       396af252 updated mirbase installation script : removed useless loading of mugqic/R_Bioconductor/3.2.3_3.2
       15d7243b Apis_mellifera.sh
       834765f9 updated module installation scripts
       996b10f9 new bash scripts for new modules
       08b52b07 MethylSeq pipeline - few minor corrections before the pull request
       dced67c0 MethylSeq pipeline - adding of the README.md file
       212451ef MethylSeq pipeline - reducing ram to 25G for bissnp step to avoid 'out of memory' issue
       e3d09b59 MethylSeq pipeline - added index creation at then end of picard_add_read_groups step to avoid error in next step (picard_merge_sam_files) when only one readset per sample
       571d0f88 MethylSeq pipeline - corrected ini file for bismark_methyl_call and bismark_coverage2cytosine section where assembly_dir was used instead of bismark_assembly_dir
       f837d91d MethylSeq pipeline - added the creation of the bam index after filtering is done, and updated picard.build_sam_index regarding the ini_section parameter
       dafcd65c MethylSeq pipeline - reducing ram for bissnp step to avoid 'out of memory' issue
       0f4507d3 MethylSeq pipeline - pipeline with the mapping quality filter step as well as with the picard_calculate_hs_metric step to get ontarget mapping rate
       67e38424 MethylSeq pipeline - updated ini file with newer module versions as well as removing all MUGQIC_INSTALL_HOME_DEV reference and setting them to MUGQIC_INSTALL_HOME
       8f2a7dc1 MethylSeq pipeline - updates after correction of the section name within the base.ini file
       e6dfe0c5 MethylSeq pipeline - corrected section name within the base.ini file
       9c5c314b MethylSeq pipeline - corrected typo in the bismark python wrapper
       af1f9f96 MethylSeq pipeline - added tmp_dir parameter to bismark_align step
       ccbf7750 MethylSeq pipeilne - addings of .md files for the first steps
       aaf3ba3f MethylSeq pipeline - corrected dependency problem within bismark_align step
       57db20f5 MethylSeq pipeline - corrected metrics input and step order
       6248dbff MethylSeq pipeline - addings for report generation
       deb45a63 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       02259ac0 MethylSeq - corrected typo within the picard2_collect_multiple_metrics which caused to always restart the step
       83cecace MethylSeq - Updates and corrections, pipeline is now fully functionnal on test data
       8bd10dbc MethylSeq - BisSNP step implemented
       77176b42 MethylSeq - debugging bed_graph step
       b58b12d6 MethyleSeq - correct output files directory for methylation_call step
       d8d7c74f Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       1422706b MethylSeq - debugging methylation_profile step
       ae3a9df7 MethylSeq - updates and corrections, wiggle traks are now generated properly
       4b96e31d MethylSeq - updated pipeline up to 'methylation_profile', still need to work on BisSNP (last step)
       7a2a721e MethylSeq - minor bug corrections
       6b0f8201 MethylSeq - First working version of the pipeline : generates bismark alignment and metrics
       89f4a608 MethylSeq - updates of pipeline up to methylation call & also removed/merged some metrics steps
       86a9ae4c MethylSeq - refine metrics calling and preparation, up to pUC19 & lambda reads step
       ce3c86d5 MethylSeq - updated pipeline until metrics step
       ff5bee30 Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       78c8fd08 MethylSeq - updates and corrections up to step 10 flagstat
       5a554729 Methylseq - updates and corrections up to step 10 flagstat
       7b0a9a98 MethylSeq - setps 7 & 8 added
       88dfebdd update bedtools version in the /bedtools.sh module installation script
       8f9c1984 BEDTools - added intersect function to bedtools.px wrapper
       63b46281 MethylSeq - debugging wrong dependencies du to wrong output file names
       17746b59 MethyleSeq - debugging file names in step 4 & 5 inorder to make step6 picard_merge_sam_files work
       5057bc04 MethyleSeq - defining parameter of step6 picard_merge_sam_files within the .ini file
       b63b02c2 MethylSeq - preparing step 6 picard_merge_sam_file by reviewing/changing file names upstream (steps 4 & 5)
       c5896841 MethyleSeq - bug fix : quoted parameters of add_or_replace_read_groups function in picard(2).py wrappers
       206f1280 MethyleSeq - corrected wrong reference to .ini file section regarding step5 picard_add_read_groups
       d0006b1d MethylSeq - missing output file is now passed to the picard.add_or_replace_read_groups function
       cae1cdf4 MethylSeq - corrected step5 picard_add_read_groups
       6294b650 MethylSeq - updated step5 of the pipeline : AddOrReplaceReadGroups
       0019f59b updated module installation scripts for bowties2 htslib python & samtools
       a8e9c4a7 PICARD - bug correction in python wrappers (picard.py & picard2.py)
       acd02572 GENAP MODULES - added bismark.sh script for Bismark module installation
       9c701035 MethylSeq - redefined output file for bismark_align step and added AddOrReplaceReadGroups function to picard(2).py wrappers
       222ce028 MethylSeq - python wrappers and scripts updates
       edb1d9bf MethySeq - corrected typo in __init__.py
       41f69bdf MethyleSeq - creation of the (empty) files as a first commit to the branch

  ehenrion <edouard.henrion@mcgill.ca>      19 commits

       4e7c266c BFXDEV-673 - updated jsonator.py for a better handling of module names & versions
       39aac2da Merged in methylseq (pull request #23)
       11c2e27d chipseq.base.ini : edited module_deeptools to remove reference to mugqic_dev
       bc8a9bed README.md updated RAC_ID export line
       a29213cf README.md edited : added the export of the $RAC_ID variable which will be used on Cedar for job sumission
       5f662301 README.md edited : added $MUGQIC_INSTALL_HOME_DEV setting for cedar
       27a43ac9 MethylSeq pipeline - edited guillimin ini file : more walltime for gatk_depth_of_coverage & more cores for bismark_dedup
       0c5f3a0b MethylSeq pipeline - changed flagstat output name so it is more obviously related to ontarget bam
       3c1f7308 MethylSeq pipeline - inverted order of input file for methylation_call step : "*.readset_sorted.dedup.bam" is now set before "*.sorted.dedup.bam" to avoid unnecessary sorting of bam files...
       f79719ef MethylSeq pipeline - GCbias and mapping_qual_fileter jobs have been added to metrics, their commands were created but not submitted...
       e77d93e2 MethylSeq pipeline - methylseq.guillimin.ini adjusted some ppn values for guillimin
       793836ea MethylSeq pipeline - methylseq.base.ini edited bismark align parameters
       4d6b5d04 MethylSeq pipeline - methylseq.base.ini added module_R
       c6a2446f methylseq.py : added missing variable
       4f247e67 README.md edited online with Bitbucket
       321816ee README.md edited online with Bitbucket
       0bb51e13 README.md edited online with Bitbucket
       a42811f9 bedtools.py edited online with Bitbucket removed unused and unfinished function genomecov...
       e76018c9 bedtools.py edited online with Bitbucket

  Eloi Mercier <emercier@jonquille.genome.mcgill.ca>      43 commits

       33f242de removed depreciated reasignment of module_snpeff in snpsift_annotate step in dnaseq.base.ini
       bad3d894 renamed and duplicated snp_effect section in mammouth.ino to mpileup_snp_effect and haplotype_caller_snp_effect in order to correctly set ppn=1
       54251b5b Moved report.kallisto job after kallisto_count_matrix since the it needs the output of kallisto_count_matrix; added dependancies to copy_tx2genes_file so it waits until kallisto_count_matrix is done
       19a3beb0 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       73b4c978 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eac41070 fixing local merge conflict resolution issue
       e4a5abe0 (hopefully) fixed conflict
       15961883 added a few more os.path.join; added an import to tools
       d4e5e0f2 removing kallisto.sh from module/dev
       5f0a9b87 remove all compiled files
       220ed537 adding os.path.join where possible
       030d09da remove unused sections in mammouth and guillimin ini files
       87289ad6 resolve merge conflicts
       46cb6453 remove commented lines
       2941bfad move rmarkdown calls to appropriate step
       b25aef3a fix ini file to correctly set walltime
       76e4c73b fix an issue with dash in sample names in report
       fa04e0a2 minor improvments to kallisto report
       b55eb983 change path to file in kallisto.Rmd
       6d3e20f8 change job names to fit sections in ini file; replace mention of samples by readsets in md files
       4fd75b29 Changed R code in kallisto.Rmd; add more columns and fix error in the stat table; changed kallisto method text
       1728a1c2 add step to generate transcript count matrix
       341b0fdd Change text RNAseq report for differential expression
       1431e2e4 update report files
       9bb9dbfb added option for single reads, added new parameters in ini file
       f9310f0f added trimmomatic stats to report
       f81641c5 added first version of RNAseq Light report
       7e53c0f3 added a mkdir command to the merge script
       efb008fc Change path of merged abudance file
       d8dedd1d fix job dependancies
       1639f593 fix exploratory function
       9d875ebb change path for call to rnaseq_light_kallisto.sh
       4fc4f3bc adding new step for exploratory analysis
       3831ba60 added a step for merging individual abundance files
       d7208762 add call to module_tools
       d2996e11 fix path to abundanceTranscript2geneLevel function
       5dc4760b adding configuration ini files for mammouth and guillimin
       6c2900a7 make it clear transcrptome must end by idx
       feca8f72 fix path, disable exploratory
       aed4e03e new RNAseq_light pipeline with kallisto, dev version
       053dbf0e added RSEM 1.3.0
       2c86a7f6 updated kallisto
       74e5c0cc update sailfish to 0.9.2

  eloi.mercier@mcgill.ca <eloi.mercier@mcgill.ca>      1 commits

       c82520a4 Approval granted by Rola. Merged in RNAseq_light_dev (pull request #24)

  eloi.mercier@mcgill.ca <emercier@abacus1.ferrier.genome.mcgill.ca>      3 commits

       22084d19 revert sailfish version change
       c8ecd2c7 ajout Salmon 0.8.0
       e2986c58 mise a jour kallisto 0.43.0

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      10 commits

       34124b69 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a9f69a48 add protocole compatibility to rnaseq_light
       cd1caeb3 Make other pipeline supporting several prtocols - BFXDEV-692
       ffdf4f6f Create 2 protocols with different steps for hicseq - BFXDEV-692
       554f8ad4 Make other pipeline supporting several prtocols - BFXDEV-692
       5ac65fc2 allow several prtocols with different step list - BFXDEV-692
       e634cd2e Create 2 protocols with different steps for hicseq - BFXDEV-692
       ece6ab0b test multi protocole pipeline
       b1098141  RNAseq & ChIPseq-   Update ini file for the new release of mugqic_tools 2.1.9 - BFXDEV-668  - BFXDEV-675
       e0571147 ChipSeq - finish ihec metrics, preprocess and reformat -  BFXDEV-675

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      21 commits

       6150a8ab remove bad output in ihec_metrics_report
       07709bff Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cc3d29d6 MethylSeq - add MethylSeq.picard_remove_duplicates.md report file
       f8a97db1 fixing conflict in resources/genomes/install_genome.sh
       2a9190be change Dedup from bismark to remove duplicat from picard
       9a9b70f6 change ihec methylseq metrics to work on single sample
       72067876 methylseq - generate ihexc report per sample and move methyl_profile lib to tools lib
       900f112b add specific temp dir to the sort step while generating bigwig
       a9b1ed3e increase general recalibration walltime in  Dnaseq to 96h
       764bb1d9 Allow bedtools.graph to support not having the other_options set in the ini
       961e0e19 removing .DS_store file
       fc1defd4 ChIPseq - address reviewer coments - BFXDEV-675
       f8bf3042 Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       cb48b451 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       1162758e RNAseq - add concat job  - BFXDEV-668
       c8591209  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ee2eab02  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       188462d5  ChIPseq-   Update (chip) IHEC metrics steps -  BFXDEV-675
       ba1f73ce  RNAseq & ChIPseq-   Update (chip) and debug (Rna) IHEC metrics steps - BFXDEV-668  - BFXDEV-675
       9e71317e  RNAseq -  implement IHEC RNA metrics step - BFXDEV-668
       7cfbb6d8  RNAseq - BFX - implement mugqic_tools module for the IHEC RNA metrics generation script - BFXDEV-668

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      1 commits

       f19e749c Merged in IHEC_metrics (pull request #21)

  pascale.marquis2@mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      2 commits

       852ef1c0 /#cluster_queue
       eb94727e python/2.7.12

  Pascale Marquis <pmarquis@lg-1r17-n03.guillimin.clumeq.ca>      1 commits

       1052ba18 update tumor_pair.base.ini

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      1 commits

       df577fb6 Minor bug fixes and addition of base exome ini

  Rola Dali <rdali@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       46371a98 starting the hicup_align step

  Rola Dali <rola.dali@mail.mcgill.ca>      88 commits

       9be34748 Merged in IHEC_metrics (pull request #27)
       130d0b56 changes to mammouth.ini to set all ppn=1; changed module spider back to module show since it is incompatible with abacus
       71cbe96a Merge branch 'IHEC_metrics' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       f93c39e5 module show changed to module spider in config.py to accelerate module checking
       a64c34fa fixing homer dependencies
       d2742939 editing homer tag directory output back to folder
       78dc2fc7 Merged in IHEC_metrics (pull request #26)
       6b978711 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into IHEC_metrics
       ac73327c homer edits to generalize methods
       649af972 resolving merge conflicts
       fd5cb387 dependencies in rnaseq metrics
       23921c9c rnaseq metrics
       62da3e6d fix ihec_metric job names
       fb831740 adding chip_type to ini to avoid pipeline crash
       b7e20902 run_spp fuctional
       85dd69ab added run_spp to calculate dnsc and rsc. TO TEST on mp2
       1459b137 adding TMPDIR to bigWig creation
       a058aa11 sample name change in ihec_metrics
       31c83eb5 fixing TMPDir issues on guillimin-to test
       80e315fb pipeline to run without input & chnage in ihec_metrics format
       14fde80f moving homer code to homer.py moduke
       d33094e2 Merged in hicseq (pull request #25)
       fe6d2618 resolving juicer conflict
       3cd06858 resolve merge conflict
       eff83a22 resolving merge conflict
       574b511e adding C3G logo
       400b590a commit for pull request
       18c1112f commit for pull request
       27b87a36 adding bedops to annotate capture file
       198416bb adding runchicago to capture hic
       742376e2 Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c26ea1df Resolved merge conflict
       cc289e6c commit before Mathieu's pull
       64cd0908 adding capture seq
       46b4f6bb merging capture hicseq with hicseq
       4b965527 chicseq
       b19d167c Merged in hicseq (pull request #22)
       a25f08bd added new chromosome contigs based on ensemble and NCBI genomes BFXDEV-670
       11b2bb95 changing perl to perl/env  BFXDEV-670
       a790fa18 added genome.py for genome related methods BFXDEV-670
       6a631245 changes for pull request BFXDEV-670
       f9a25525 changes for pull request BFXDEV-670
       0d65329d changes for pull request BFXDEV-670
       aa2e874e edits for merge request:wrappers. BFXDEV-670
       de38b5a3 BFXDEV-670 hicseq merge edits
       d69045fe creat_hic_file edits BFXDEV-670
       1649642e adding juicer.sh installation script BFXDEV-670
       665ea381 Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       c4adfa90 testing create_hic_file BFXDEV-670
       4a67afe5 modified hicseq.briaree.ini and batch BFXDEV-670
       a37fe8a8 deleted files moved to mugqic_tools BFXDEV-670
       bcf5b1a1 moved genome_digest before ncrna method which is failing in install_genome.sh BFXDEV-670
       861cc799 added module install files and genome digest BFXDEV-670
       15230602 changes for pull request
       ad4a1412 commit before merge changes
       98d79471 added samtools_bam_sort BFXDEV-670
       6fd3b938 edited ini files BFXDEV-670
       18e07747 samtools_bam_sort and multiqc_report in testing
       3c9468bc Fixed hicup/TAD/HiCPlotter restrart bugs BFXDEV-670
       6a58aa4b split chr interaction matrices from plotting BFXDEV-670
       dca827d2 HiC v1.0 is ready for testing on guillimin and abacus BFXDEV-670
       2b40a8c7 bam merging in testing BFXDEV-670
       65763b14 samtools_merge_bams in testing BFXDEV-670
       d2ff7116 added petagDistDistribution to HomerQc plots BFXDEV-670
       e8edda6a added output_dir property to reduce code redundancy
       ded82f6f fixed path in identify_compartments
       1b31f400 added identify_peaks BFXDEV-670
       d212eba1 compartments and TAD identification now working BFXDEV-670
       384d836e chr and genome interaction system now working BFXDEV-670
       6714aabc interaction matrix plots testing BFXDEV-670
       edadc744 homer archiving and Qc plotting now working
       8c94b93f first 6 steps working BFXDEV-670
       5c3eb287 resolving git merge issues
       e1d46d6f Merge branch 'hicseq' of bitbucket.org:mugqic/mugqic_pipelines into hicseq
       399cc72e syntax changes
       c545fc62 added fastq_readName_Edit() BFXDEV-670
       2b7ee9ea mammouth.ini cpu set to 1
       4afc530b genome digests on mammouth
       e4e52d38 update genome digest files
       92773981 homer_tag_directory archiving in testing BFXDEV-670
       acb862f1 basic formatting changes
       53257d6b make_tag_directory in testing
       148eb65f works to produce hicup bam and library Qc. ITS ALIVE :)
       958242d3 ITS ALIVE
       c446402c need to expand  variables
       099e5ef2 pipeline now accepts enzyme
       ffbdf83a hicup_align producing script; to test tomo
       884af260 initialising hicseq analysis pipeline

  Xiaojian SHAO <xshao@lg-1r14-n04.guillimin.clumeq.ca>      2 commits

       6ddab51b add ppn to wiggle_tracks step
       06bd6a7f add ppn to wiggle_tracks step

  Xiaojian SHAO <xshao@lg-1r17-n03.guillimin.clumeq.ca>      3 commits

       6fd6ec5f Merge branch 'methylseq' of bitbucket.org:mugqic/mugqic_pipelines into methylseq
       cd938a3a add walltime to bismark aligner
       8a5b8263 ppn_Changes_in_Guillimin.ini

  Xiaojian SHAO <xshao@lg-1r17-n04.guillimin.clumeq.ca>      1 commits

       db2901fa methylseq: edit on ppn setting. -xiaojian

2.3.0        Mon Feb 27 13:40:01 2017 -0500        82 commits

  Edouard Henrion <edouard.henrion@mcgill.ca>      13 commits

       7717c35e Pre-release - adding the installation scripts for all the new software used by tumor_pair pipeline
       8eada670 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8e17bca9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       59cf9407 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       026ccb11 BFXDEV-602 - PacBio Assembly pipeline now contains a new optional step (step 10) to circularize the successful/clean assembly contigs with circlator - also comes with some minor unrelated updates
       06e095f6 updates brought to many module install scripts, and adding of new modules
       9d565ba2 RNASeq - corrected a typo inserted after correting bedtools bug...
       982a0145 RNASeq - corrected a bug bedtools.graph function : samtools_options now handles reverse strand specific parameters, avoiding an empty begGraph for reverse strand
       72d1a3f2 updating python & python libraries installation bash scripts
       6242c97c BUG correction within mpileup function : parameters were shift after introduction of 'ini_section' parameter for tumor_pair prpose
       1ca38fa7 DNASeq - bug correction after merging tumor_pair branch to master
       e6b482bf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d8876d7 FastQC - updated bash install script to match our standards

  ehenrion <edouard.henrion@mcgill.ca>      1 commits

       0e875a56 README.md edited online with Bitbucket

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       b5b1d51e tumor_pair - add comments to dict2beds function - BFXDEV-521
       e0590a42 tumor_pair - add feature to dict2beds function - BFXDEV-521
       7c0efb4d Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       bd7e169a tumor_pair - modify ini file to take into account the new version of scalpel (CVMFS) - BFXDEV-477
       8e735ac6 tumor_pair - add space charter before scalpel option - BFXDEV-478
       33a9347a pull origin Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       5ed57839 tumor_pair - add the two pass option to scalpel - BFXDEV-478
       bf963a38 tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing - BFXDEV-476
       18110f56 tumor_pair - corect bugs and typo in bed file integration of scalpel & rewrite to speed up the bed parsing

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      39 commits

       d6be75aa remove .gitignore
       3821384f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8e31719a rnaseq - increase abacus ressource for STAR
       6f5a17b7 remove conflict with tumor_pair changes
       07cab4fc remove bfx/samtools.py conflict with tumor_pair changes
       abecc6c7 remove bfx/picard.py conflict with tumor_pair changes
       4368b7ca merge and remove conflicts
       2d5b7c56 merge and remove conflicts
       dfcbf627 tumor_pair- update germline_loh ensemble
       376ff926 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       07b9e8b5 add encode  sample+readaset+pipelien to log system
       f4c9b602 Tumor pair - add comment on the bed system
       ba477024 Tumor pair - change vardict input output to correct deepndency
       0e20097f Tumor pair - change varscan input output to correct deepndency
       4ab037ec sequence_dictionary.py - split_by_size - correct bug create a first job with everybody when the nb of split was too high
       b5a7b337 tumor_pair - add output bam index file when only 1 realignment is produced
       a5b4491b tumor_pair - remove file name bug in indel realignement mv
       4429b197 tumor_pair - rebuild gatk indel recalibration input/output scheme in tumor
       b9030e54 bfx/gatk - indel realigner add target list file as input for dependency
       440001c6 tumor_pair - remove duplicated lines in guillimin.ini file (from base.ini)
       6c1e6060 tumor_pair - correct symlink creation
       b15f1184 tumor_pair - remove few issues (paths, dependency, input/output)
       329a56c7 tumor_pair - put WGS as default (instead of WES)
       271ba292 remove issue with uncorrect interval list when no bed is attached to the project
       517e8c3a adjusting tumor_pair
       c049fe24 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       c2c0986d tumor_pair - debug bfx/bcbio_variation.py
       5c6034d2 tumor_pair - remove conflicts
       d2c00151 rtumor_pair - emove conflicts
       7c4bd8a1 remove conflict
       c67b8375 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       9cbd2220 tumor_pair - correct typo in scapel download adress - BFXDEV-477
       760faaa0 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       8d2a2433 tumor_pair - add scalpel install script - BFXDEV-477
       e6beef4c tumor-pair - mv tools.py from highcov branch to tumor_pair - BFXDEV-475
       1bc346fa tumor-pair - adding bedfile spliting process - BFXDEV-476
       5acad27d tumor_pair - extract tumor_pair code for the high coverage branch - BFXDEV-475
       6c73c2b3 tumor-pair - adding bedfile spliting process & start implementation in tumor_pair.py - BFXDEV-476
       905c1acd tumor_pair - extract tumor_pair code for the high coverage branch

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      5 commits

       192f4d35 Bug squashes and speed improvements to ensemble processes
       40df5d93 Merge branch 'tumor_pair' of bitbucket.org:mugqic/mugqic_pipelines into tumor_pair
       93a2142f Modified to remove analyses of alt contigs + split Varscan2 by chromosome
       8d3f03fc Beta version: module files created and code tested on wes and wgs on abacus and guillimin
       755b0c06 update to Mutect2, added Vardict, re-added samtools and added ensemble approach

  robert.eveleigh@mcgill.ca <reveleig@abacus2.ferrier.genome.mcgill.ca>      4 commits

       e94352f0 Fixes to samtools modules calling for pairs
       9729909b Numerous speed improvements and addition of fast variant calling
       3f1aceb1 Additional fixes to resource allocation issues
       0f04f2af Dealt with comments, added paired indel realignment, varscan2, and seperate somatic and germline call files

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      11 commits

       67b5c532 Debugging samtools pair bcftools
       3a4996e3 Fixes to samtools due to dnaseq changes
       a920ef82 Update README.md
       2c7ce408 Bug squashes and preformance improvements to ensemble process
       3349a818 Bug squashing and speed improvements to ensemble process
       15a7f87f Dependency and other bug fixes. Guillimin ini and install scripts updated
       a21682a6 tumor_pair - fixes to tumor_pair README.md
       19664340 tumor_pair - fixes to README.md and resolving of conflicts
       c53612c4 tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       dd6df556 tumor_pair - corrected BQSR for WGS, added README.md, added removable_files
       c708f88e Updates/fixes from guillimin test

2.2.1        Mon Dec 19 10:57:33 2016 -0500        212 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       cd7c8e74 Added proper error message when running script with too old Python version.

  Edouard Henrion <edouard.henrion@mcgill.ca>      138 commits

       68f03666 GenAP Pipelines - updated tmp_dir variable within all the .base.ini files for a better use of the memory in the compute nodes on abacus
       c364c3ec modules - updated version of VSEARCH (2.3.4) within the module installation script
       ed99e0e8 GENAP PIPELINES - updated all the .base.ini files to set tmp_dir to /lb/scratch/ehenrion instead of /lb/scratch/ on abacus
       d0faf2ed RNASeq - added __init__.py
       a8d11011 PICARD - bug correction in python wrappers (picard.py & picard2.py)
       7d2e9b2f install_genome.sh - corrected MUGQIC_INSTALL_HOME_DEV link
       b2ec14ee install_genome.sh - set abacus2 variable to manage mammouth execution properly
       c3549edf install_genome.sh - updated for a better handling of mammouth cases
       47269324 Homo_sapiens.hg19.sh install script - updated URL to retrieve dbSNP vcf
       b0dfe826 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       465b913a RNASeq denovo Assembly - updated mugqic_tools version to 2.1.6
       55e64564 MUGQIC MODULES - updated bedtools version to 2.26.0 within the bash install script
       797fe54f MUGQIC MODULES - updated samtools version to 1.3.1 within the bash install script
       f8467465 MUGQIC MODULES - updated bowtie2 version to 2.2.9 within the bash install script
       7ef8596b MUGQIC MODULES - updated bowtie version to 1.1.2 within the bash install script
       eca9f4d3 MUGQIC MODULES - updated bismark version to 0.16.3 within the bash install script
       b92f63f2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2b0f2e01 RNASeq_denovo_assembly - updated guillimin specific .ini file for trinity & insilico_read_normalization steps
       5c116093 DNASeq - updated guillimin specific .ini file for snp_effect step
       639668a1 DNASeq - updated .ini file for snp_effect step
       14f4c0f5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b3a054e5 PICARD 2.0.1 - updated picard2.py wrapper with proper syntax
       7a20e2f8 briaree pipeline .ini files
       b4b5e8c8 RNASeq_denovo_assembly - updated mammouth .ini file
       b3ffcc7f python.sh - updated python version (2.7.12) within the bash installation script & added all the python library installation steps : no need for python_lib.sh anymore
       731b0299 Homo_sapiens.GRCh38.sh - updated Ensembl version (85) & dbNSFP version (3.2c)
       ba150756 install_genome.sh - update cmd_or_job function to automatically handle mammouth environment cases
       5d07d310 exonerate.sh - update version to 2.4.0
       c1ac0e43 RNASeq_denovo_assembly - back to hmmer version 2.3.2 (from mugqic modules)
       1f647579 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       16c34640 exonerate - within the bash install script : corrected broken link to the archive
       875e3eb0 Genome reference .ini files updated for each reference installed on CVMFS
       8c4c2ed0 module install script - corrected emboss.sh : unresoved conflicts are now resolved
       21a06d7f module install script - corrected emboss.sh : unresoved conflicts are now resolved
       5f46cf69 Python packages - added TEToolkit package to the installation scripts
       f4a2adb5 FASTQC - updated installation script with newer version and integration of /lb/project/mugqic/analyste_dev as a possible installation path
       a396762d RNASeq_deNovo_assembly - updated ini file with integration of picard v2.0.1 and working version of trinity (2.0.4)
       916b8d6a PICARD 2.0.1 - added picard2.py to the bfx tools
       ef03d6bf commiting minor updates, prior the new release
       aa0f55fb RNASeq_denovo - briaree specific config for insilico_read_normalization_all
       29fc52ab Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd77e98e corrected version of mugqic_tools in rnaseq.base.ini
       daad645b corrected version of mugqic_tools in rnaseq.base.ini
       e2d140ec RNASeq de-novo assembly - bug correction when calling trinity with new default parameters
       87d68339 RNASeq de-novo assembly - updated ini file for briaree
       1605b8fe modules - new STAR version handled in the installation script
       8f88daae RNASeq de-novo assembly - updated trinity call with newer version of trinity and new default parameter
       73ca888a RNASeq - update differenial expression to handle local fit option if needed (i.e. if parametric dispersion fit fails)
       ec71e94a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3d0d765e new reference genome installation scripts (new species)
       0029c2a3 metaMarkerSeq - guillimin .ini file adjustments regarding qiime
       e4c1a5b9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b5cc05d5 added butter installation script
       8ec7a74d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       16f2ae77 qiime db installation scripts and some new module scipts
       8fe19f83 chimera_unite.sh added
       ad09a1ee new module installation scripts
       a4f4b053 some updated & new module installation scripts
       8fd6a71d updated genome installation scripts
       083ff4e5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       846d9319 metaMarkerSeq - removed bug with amplicon_type within the .ini file
       2ba2d4c1 metaMarkerSeq - updated python & bash script as well as .ini with the newer silva db version
       4504c46b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e68ffa31 RNASeq - fixed bug in gq_seq_utils_exploratory_analysis_rnaseq
       903ef8db mammouth ini files added for dnaseq_high_coverage
       def0206a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c822e0d7 mammouth ini files updated for pacbio & dnaseq
       887e5189 adding pipeline ini files specific to brairee environment
       e5e3585c DNASeq - review some parameters within ini file during the pre-release pipeline testing
       cece203d RNASeq - allow sample names to have '.' within their name without crashing during gq_seq_utils_exploratory_analysis_rnaseq
       46d069f2 fixed typo
       54f260b4 RNASeq - fix mammouth setting for htseq_count step in rnaseq.mammouth.ini
       ada67381 RNASeq - restoring bfx/report/RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq.Rmd
       35c5cad3 ChIPSeq - metaMarkerSeq - adding briaree specific ini file
       f0112e13 RNASeq - correct default genome in the ini file
       4fb234d8 pushing ampliconseq.guillimin.ini after testing on guillimin
       90ff7dca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4ec8c829 pushing ampliconseq.mammouth.ini after testing on mammouth
       4d8e5832 Corrected typo in dnaseq.base.ini
       6e9c3b2a correct mammouth master branch divergence
       8ab80b30 updating module and genome scripts from mammouth
       7e284906 remove conflict
       2c439fe1 conflicts resolution
       4db2a9e4 some residual diffs to commit
       de4c0436 redo the call of bcftools module within samtools.py
       b69fa6e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ad24f0cd resolving conflicts after merging...
       6330a7f4 commits before merging
       5e4eac40 DNASeq/RNASeq - code correction after testing on briaree
       580589a9 metaMarkerSeq - Merging branch to master git add resources/genomes/Bos_taurus.UMD3.1.sh resources/modules/perl.sh resources/modules/star.sh resources/modules/vcftools.sh
       0fae6735 metaMarkerSeq - reformat the report by removing some redundancies in the paragraphs and ensuring correct links to the tables/figures - BFXTD-26
       8873109b BFXTD-26 - corrected the template and links used for report generation - MetaMarkerSeq
       a0c0188a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e18bb247 correct a bug in wiggle step which causes *forward.bedGraph and *reverse.bedGraph files to be identical
       850b148c metaMarkerSeq - commits after reviewing pull request comments sent by Marc Michaud - BFXDEV-453
       23901609 threshold instead of treshold
       51848d9f metaMarkerSeq - debug Krona input generation and updated md files for report - BFXDEV-453
       475035b1 meatMarkerSeq - get rid of remaining  variables - BFXDEV-453
       0fd827d5 AmpliconSeq - debugging report template files
       076e6f71 bug fixes and file names correction
       04414155 ampliconseq - tuning of the ini file - BFXDEV-434 BFXDEV-448
       feac3c05 ampliconseq - new/updated module install scripts related to the pipeline needs - BFXDEV-434 BFXDEV-448
       3865cfaf ampliconseq - code review & debug prior to relase - BFXDEV-434 BFXDEV-448
       664910aa AmpliconSeq - modified some parameters in .ini file & other minor syntax changes within the wrapper - BFXDEV-434
       b22e6e2f ampliconseq - conflict resolution - BFXDEV-434
       28788ceb ampliconseq - code review of ampliconseq.py - BFXDEV-434
       d07b3d86 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ba4bca97 BFXDEV-523 - error correction & version update in the gemini installation script
       4e7e0a60 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f22df952 BFXDEV-516 - Macaque genome installation bash scripts + other minor updates
       90d34cd6 BFX-4532 - added prepend-path entry for dnaseq_high_coverage in mugqic_pipelines.sh
       7f41180a RNA-Seq - remove useless parameter from rnaseq.base.ini
       7442c01e RNA-Seq - increase default walltime for star_align & star_index steps
       c9d36941 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a9e67a5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       839b257d genomes - updated versions of genome installation scripts, essentially fixing STAR indexes installation - BFXDEV-494
       a3da825c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       192a9717 STAR - updated STAR version (to 2.5.1b) in the installation script (star.sh) as well as in the RNA-Seq pipeline .ini files - BFXDEV-514
       646ad6ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       394913df Version bump to 2.1.0-beta
       0dac262f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7a78685a BFXDEV-490 - recommit after conflict resolution
       8fc43e3c BFXDEV-490 - updating chipseq.base.ini, dnaseq.base.ini, pacbio_assembly.base.ini, & rnaseq.base.ini
       cfd56b99 ampliconseq - module - adding vsearch.sh in the module script folder - BFXDEV-434
       fb1f4d2a ampliconseq - base ini file updated - BFXDEV-434
       ed7611f6 ampliconseq - updated ampliconseq ini with new module versions - BFXDEV-434
       8c047ccc ampliconseq - updated module versions - BFXDEV-434
       565ab998 Merge branch 'ampliconseq' of bitbucket.org:mugqic/mugqic_pipelines into ampliconseq
       d2ef4f43 AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEV-434
       b15af326 AmpliconSeq - parameter change in the .ini file for otu_picking step - BFXDEX-434
       ace14e4d AmpliconSeq - revised and standardized code of ampliconseq.py - BFXDEV-434
       76c924e6 AmpliconSeq - updated ampliconseq.base.ini - BFXDEV-434
       83d3b0e4 AmpliconSeq - new VSearch library added to bfx - BFXDEV-434
       96ce6075 AmpliconSeq - updated tools.py with AmpliconSeq functions in bfx - BFXDEV-434
       b886764d AmpliconSeq - new Qiime library added to bfx - BFXDEV-434
       bb008af5 AmpliconSeq - new Krona library added to bfx - BFXDEV-434
       d9016c86 update of resources/modules/mugqic_tools.sh - BFXDEV-490
       06115849 ampliconseq - updated syntax and some corrections

  Edouard Henrion <henrione@briaree2.rqchp.qc.ca>      2 commits

       91a47a54 RNASeq - adding .ini file for briaree
       ce7f2aed small adjustments after cloning/before testing on briaree

  edouard.henrion@mcgill.ca <ehenrion@abacus1.ferrier.genome.mcgill.ca>      1 commits

       5806c4a6 ampliconseq - bump README.md to the latest version

  ehenrion <edouard.henrion@mcgill.ca>      3 commits

       62133d09 README.md edited online with Bitbucket
       d988c477 README.md edited online with Bitbucket
       5edde688 README.md edited online with Bitbucket

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      1 commits

       0c862406 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       72273420 mini and spades modules
       be71d1c8 nxtrim and quest modules updates
       5ec026e6 dev install scripts for MAKER

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       2f9db74e HighCoverage - add missing README.md file
       21f73de2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5e66e13a Removing unmaintained pipelines (PUURE & rRNATAGGER) from master
       435f7d53 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       33667103 RNAseq synchromize mammouth ini with the base ini (missing tuxedo_hard_clip) - BFXDEV-515

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       a1574375 AmpliconSeq - remove ini recursive interpolation
       84a3e9b9 remove conflict
       0e5bd2a7 High coverage - debug vt and gemini step - BFXDEV-542 BFXDEV-541
       28ce1e1e ampliconseq - remove dependencies issues and remove try instances - BFXDEV-531
       342169d8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ed67908b DNASEQ - bwa always set LB in RG tags: when library barecode not given used sample name instead - BFXDEV-354
       e28eff7e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       73bce2e7 ressource - python lib: add futures for multithearding
       a7c0b13e update install module general script

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      5 commits

       8a0bb768 dnaseq.py  - gatk DoC job name wasn't following the rule: <stepName>.<sampleName>.<OtherInfo> . corrected
       5056c01d Merged in irp_modifications (pull request #17)
       4991ce6d change mail for bug & issue in  README.md edited online with Bitbucket
       ebe9e989 README.md edited online with Bitbucket
       0cdc302a README.md edited online with Bitbucket

  mmichaud <marc.michaud@mail.mcgill.ca>      20 commits

       ed1aece1 Increase max jobs per step to 24 (cover most case without overloading the cluster) (cherry picked from commit 4b06703)
       58366498 BFXDEV-559 Error when a lane contains only double-index librairies on a single-index run (cherry picked from commit 19bd532)
       b2f9c7ee Increase even more STAR walltime to allow fetching the reference index from cvmfs.
       a8b9ffe3 Don't mark the blast job as failed if there is no blast results. If grep can't find any results, it exits with a code 1. After, when we test the $PIPESTATUS, we intercept that error code and mark the job as failed. By adding a " && true", we preserve exit codes but we reset the $PIPESTATUS.
       5c5d9336 Increase STAR walltime to allow fetching the reference index from cvmfs.
       61b608ba BFXDEV-538 Add perl for bed2interval_list
       4d724718 BFXDEV-538 Update modules version
       e39b2027 BFX-4696 Use a newer version of python having all the needed packages.
       d141b2fd BFXDEV-544 VerifyBamId: Don't filter the variant file with the genomic target file.
       0237e8ac Merge branch 'master' into irp_modifications
       5849c790 Illumina Run Processing: Accept "." in asembly name in the genomic database. (fix regex)
       e020c645 Illumina Run Processing: Accept "." in asembly name in the genomic database.
       8806df0e Illumina Run Processing: Use cvmfs version of verify_bam_id.
       6735cae3 BFXDEV-512 Illumina Run Processing: Use "Parc" silva database.
       8aa74f21 Illumina Run Processing: Add "nanuq_environment" configuration variable to be able to ease testing on the QC environment.
       3d2bec95 BFXDEV-533 Illumina Run Processing: The exclude bam option doesn't exclude the sorted bam.
       776850b7 BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       7d1694dc BFXDEV-528 PacBio Assembly: Change settings for assembly up to 15Mb.
       edfb2e4b BFXDEV-527 Illumina Run Processing: Merge jobs of a task when their count exceed a configurable threshold
       cd6cfde6 Tweak cluster parameters : Use parallelGCThread for all available cores. BvaTools run faster with only 4 threads.

  ptranvan <patrick.tranvan@mail.mcgill.ca>      24 commits

       46ee612b Eulidean distance option for PCOA
       8ca5e0ec README link correction
       18b9818a Module perl correction
       d63b3139 Module perl correction
       b4fc989f Format edition
       39c2882b . file deletion
       c2cd2c3e BFXDEV-434 README update
       e0d27cc4 BFXDEV-434 Create inis for other clusters
       f4475b93 BFXDEV-434 Deploy db files in $MUGQIC_INSTALL_HOME
       dc7f9fa1 - Other pipelines have been added to the branch. - Create inis for other clusters, leverage overlaading of parameters. - put db parameters in DEFAULT section,. - QIIME section has been exploded in the config file. - Tutorial in README file.
       27f3fa58 README and configuration file correction
       5a9fdb72 Configuration file modification.
       3419fac4 Adding alternative normalization method
       a7ade41e Dependance correction
       d4e383d3 add a step (close ref)
       6b63f3f5 CLuster algorithm change
       a5d88761 tutorial.txt modification
       72ef558d Minor modifications (report and option)
       9cf47dcb Amplicon-Seq pipeline + report upload for test.
       c9c6125c Change name: metagenomic pipeline to Amplicon-Seq pipeline Adding new features (alpha, beta plots) Final step implemented
       d7db1376 Adding metagenomics pipeline (1st part)
       075ec74d Adding metagenomics (amplicon) pipeline.
       990ba4b4 Merge branch 'metagenomics' of https://bitbucket.org/mugqic/mugqic_pipelines into metagenomics
       f779d06e test

2.2.0        Mon Feb 8 12:04:44 2016 -0500        405 commits

  dbujold <david.bujold@mail.mcgill.ca>      1 commits

       dc03d86e Added link to the GenAP project in front page.

  Edouard Henrion <edouard.henrion@mcgill.ca>      63 commits

       477b3324 Version bump to 2.2.0
       caa197dd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1d87ea7e BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       371dd66a add report section in DNA-Seq High Coverage pipeline ini file
       514d5b61 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       be6d52a1 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       61ac85c2 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       042bb038 modules - updated vcftools VERSION in vcftools.sh -  BFXDEV-490
       f71fc014 minor update in module file created by perl.sh
       d3c2bbfe add report section in DNA-Seq High Coverage pipeline ini file
       cd5d58c0 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       7d0e6f00 updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       7bf0a392 ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       349ae58e modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       ab4b16db module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       9e04f6a8 modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       e93f812a minor update in module file created by perl.sh
       fe2a9cba rnaseq & rnaseqdn ini files updated after testings on guillimin - BFXDEV-490
       39e91f13 add report section in DNA-Seq High Coverage pipeline ini file
       e3bce293 added guillimin specific ini file of DnaSeq_High_Coverage Pipeline - BFXDEV-490
       ecc918b5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a4f79773 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       442830db updated ini files - BFXDEV-490
       175b162e BFXDEV-490 - resolving conflict on guillimin
       36b6f5db BFXDEV-490 - updated base.ini files for chipseq dnaseq & rnaseq pipelines
       d1b51cbd Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       a3e2ab82 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ec80af33 modules - updated vcftools VERSION in vcftools.sh - BFXDEV-490
       b6e687ae minor update in module file created by perl.sh
       97bd1559 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       682edbee mugqic_tools.sh VERSION=2.1.5, again...
       fd3f0e3a updated version of mugqic_tools.sh - BFXDEV-490 - BFXDEV-501
       424768a7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4b4c08ac ini file for PacBio Assembly pipeline is updated with new module versions - BFXDEV-490
       f6415e12 modules - updated shebang for all the perl scripts installed by trinity.sh - BFXDEV-490
       045c237f module - corrected SQLite archive url in trinotate.sh - BFXDEV-490
       d2072f58 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c22a5961 modules - corrected typo in perl.sh when creating module file - BFXDEV-490
       f89dbb33 modules/genomes - updated module version calls as well as database releases - BFXDEV-490
       3970c215 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       17531c78 module - debugged scalpel.sh script - BFXDEV-490
       a05c71b1 module - updated scalpel version - BFXDEV-490
       1c9d8682 module - debug  in snap.sh  - BFXDEV-490
       eee7e567 module - updated version of Python (to 2.7.11) in gemini.sh - BFXDEV-490
       8018094b module - removed module load calls in shortstack.sh - BFXDEV-490
       8c699539 module - updated star.sh - BFXDEV
       071a17b2 module - updated shortstack.sh snap.sh - BFXDEV
       3fa5330c module - corrected typo in sphinx.sh - BFXDEV-490
       976b8565 module - vcftools updated version to 0.1.14 - BFXDEV-490
       0d5709bc module - debugged gemini.sh - BFXDEV-490
       add9770b module - samtools updated to 1.3 - BFXDEV-490
       15428f89 module - debug the name of the archive - BFXDEV-490
       7ab3bbfb module - UCSC version set to v326 instead of 'latest' - BFXDEV-490
       16571fcc update version of bwa module - BFXDEV-490
       bf6422a8 some more updated modules - BFXDEV-490
       0cf2a4d5 module & genome updates - BFXDEV-490
       39c17582 module updates for to the release - BFXDEV-490
       bb09f04c resources/modules/dev/epacts.sh has been removed (really)
       0815aec6 resources/modules/dev/epacts.sh has been removed (really)
       27058bcf BFXDEV-490 - update of resources/modules/mugqic_tools.sh
       0021aaa2 Merge branch 'gatk_variants' of bitbucket.org:mugqic/mugqic_pipelines into gatk_variants
       a41229e6 gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted - BFXDEV-481
       b3ae471e gatk_variants - correct getDups() in ignstats.py so it ignores '?' as a dupplication rate when library ID is omitted

  Edouard Henrion <henrione@ip03.m>      2 commits

       839c4644 updated chipseq.base.ini after mugqic_tools update - BFXDEV-501
       476279f7 updated chipseq.base.ini after mugqic_tools update - BFXDEV-501

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       565ea8b2 BFXDEV-491, missing argument to re.sub for the SINGLE_END STAR commends case
       bc9c0dfe Reverte back to older scheme. Was assuming bacteria only.

  Francois Lefebvre <lefebvrf@gmail.com>      37 commits

       ae9667fb Mammouth rnaseq ini required cluster_cpu=-l nodes=1:ppn=1 for new steps related to rRNA estimation
       e21b18cf -S flag was missing from the last samtools view command in the hard clip command
       f5734af0 Install scripts
       6eabc352 minor mod to R_Bioconductor (removed tabs in HERE DOCS)
       432b3abc sleuth R package
       93c3b1da Removed old R installation scripts
       7758ea3c Added slash to URL to be able to retrieve latest Bioc version
       b80611a2 sspace-longread dev install script
       b626deb2 Kallisto install script (abacus only)
       0f9449a9 popoolation install scripts modifications
       dd9e0313 more notes on pacts
       1e9442f6 Added notes to EPACTS install script
       ee112de2 Fastx and EPACTS install scripts
       a1b52899 install scripts for ray, proved, sortmerna
       a0867442 Added celera and loRDEC install scripts
       c56d2979 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01d05df1 Merge branch 'knitr_test'
       62efc5fe Added PBJelly install script + networkx python module
       16910029 Fixed typos
       394fca5e Replaced library calls in rmarkdown.py with direct ::render call. Added missing section header to exploratory report.
       7d614146 no message
       cf13cafb Modified exploratory analysis report text
       a8fd54f4 import rmarkdown in wrapper
       abd9f895 damn dots
       336b1569 Testing out markdown::render paradigm
       95df3015 row.names=F for exploratory summary table
       0f11b557 Working directory madness fixed
       4de25fc0 Trying out rmarkdown::render() instead
       22d5cc6f Trying html output
       2c3caa85 knitr would set the wd to input document location...
       c6248977 dir.create problems
       960c5e51 quotes missing
       01e5ea50 no message
       94b1e9b2 EOF is simpler
       5742f004 Had forgotten the module calls
       65e698e5 knit for exploratory + other changes
       39f86fe8 Created Platanus install script

  Gary Leveque <gary.leveque@gmail.com>      2 commits

       248c455c pacbio_assembly.base.ini edited online with Bitbucket
       93cfad19 pacbio_assembly.base.ini edited online with Bitbucket --changed smrtanalysis_version = 2.3.0.140936.p2 to: smrtanalysis_version = 2.3.0.140936.p4

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      1 commits

       2d582b48 changes necessary for bacterial RNAseq using STAR --see BFXDEV-449

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      52 commits

       f7e38cba BFXDEV-405 correct single library issue in insilico read normalization
       a7457994 Merged in rnaseq_denovo_assembly_new_features (pull request #10)
       b1ecdcc0 BFXDEV-397 resolved rebase conflict resources/modules/verifyBamID.sh
       23a6e08e BFXDEV-397 PRJBFX-1187 added genome and software installs for verifyBamID, resolved issues from pull request #10
       296aab1a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       91707a58 BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features removed dev modules
       f80579a0 BFXDEV-397 resolved issues from pull requests # 10 https://bitbucket.org/mugqic/mugqic_pipelines/pull-requests/10/rnaseq_denovo_assembly_new_features
       8ebbedb9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d492790a BFXDEV-397 added spaces to table captions to goseq and dge filtered outputs
       c03656d6 BFXDEV-397 added DGE and GOseq using filtereg contigs (based on trinotate annotations)
       c22ce84c README.md edited online with Bitbucket
       706303f1 BFXDEV-439 force installation or rmarkdown and knitr
       ab8526cf BFXDEV-439 resolved rebase conflicts rnaseq_denovo_assembly_new_features
       716c646b BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       028430fc BFXDEV-450 Change the mugqic_pipelines templates to add C3G logos and supporters
       1905e38b BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, commit new libraries
       13cefa95 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report, regenerated rnaseq_de_novo pipeline README markdown page
       8d5737a9 BFXDEV-439 Reorganize rnaseq_de novo assembly to send all trinotate/rsem code to bfx libraries, added trinotate annotations report
       b40353bc BFXDEV-397 correct dependency problem in exploratory analysis using filtered isoforms
       59fdc4e2 BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       107fe2a9 BFXDEV-396 added exploratory analysis using filtered transcripts. changed report markdown file
       b85f60b4 BFXDEV-396 added exploratory analysis using filtered transcripts
       17e2143d BFXDEV-423 BFXDEV-399 BFXDEV-397 corrected blastx, abundance estimates to matrix and generated tabular and fasta filtered assembly using python SeqIO
       de22cfa8 BFXDEV-432 define module picard rnaseq_denovo_assembly base ini
       38207383 Merge branch 'rnaseq_denovo_assembly_new_features' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       f249f9cf BFXDEV-397 rnaseq_denovo_assembly.py and rnaseq_denovo_assembly.base.ini rebased from master
       9639c9f5 BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       1922f922 BFXDEV-397 tested filtered trinity output using real data (1 sample)
       afad4375 BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d223c745 detected bugs during tests
       6d0e428a  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       f91e4272 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       c4597fa9 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       b2afc9e0 detected bugs during tests
       bd91a583  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       2939c0e2 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       5ef76ae6 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       641a6532 BFXDEV-397 corrected report exploratory analysis using knitr, tested using real data
       cb1e022a BFXDEV-397 tested filtered trinity output using real data (1 sample)
       d9a9c464 BFXDEV-397 merging conflicting versions of rnaseq_denovo_assembly.base.ini and rnaseq_denovo_assembly.py files
       5c5b83e9 Merge branch 'rnaseq_denovo_assembly_new_features', remote branch 'origin' into rnaseq_denovo_assembly_new_features
       62acde9a BFXDEV-397 fixing differences when rebasing from master
       33d0dd0d BFXDEV-397 component/contig filtering step based on annotations and predictions made by trinotate, added exploratory analysis based on raw counts generated by RSEM, modified differential expression analysis to use deseq and edger python libraries and to merge with annotations using mugqic tools parseMergeCsv
       d3e9ff6b detected bugs during tests
       73ec3415  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       8a0bb017 BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       4d099e5a BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       126102e7 detected bugs during tests
       f9d902ac  BFXDEV-396 add parse trinotate output to extract blast and go annotations for genes and transcripts
       0cb2d532 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rnaseq_denovo_assembly_new_features
       5579bc4d BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo
       ff1f6dae BFXDEV-396 add GOseq analysis of annotated components in RNAseq de novo

  lletourn <louis.letourneau@mail.mcgill.ca>      18 commits

       74503b2b Merge branch 'master' into highCoverageVariants
       bfdf0232 fixed io_buffer default
       8b2cfcc9 BFXDEV-392 First implementation of high depth calling
       d1249587 BFXDEV-370 Fixed merging and output naming bugs
       7381a036 Merge branch 'tumorPair' into highCoverageVariants
       da6e1fb2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fffb12b6 BFXDEV-392 need varscan for high coverage
       4bdd5530 Merge branch 'master' into tumorPair
       8ae0dd3d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       23c3e596 set default ppn for igv to 2
       bba9c565 BFXDEV-379 remove -gpfs from the ini files
       6cf6ac40 Version bump to 2.2.0-beta
       8e1345df BFXDEV-370 added merging step
       27983c2f Merge branch 'master' into tumorPair
       61e35f35 BFXDEV-370 Added indels and COSMIC
       ffe5e6d7 BFXDEV-370 fixed scalpel script, added LD_LIB_PATH
       8c0cfc75 BFXDEV-370 Added the scalpel module
       7452be9b bvatools version bump

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      32 commits

       a57fd6d3 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       129d564c RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       46958811 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       226c59b4 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       1cbb0f17 chipseq - update module -  BFXDEV-490
       9e244862 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       501c3ac0 DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       bc5920a4 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       849d732a RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       db6c64a7 chipseq - update module -  BFXDEV-490
       65e38616 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       e1637bf8 DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       c9e8db09 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       0d077aa6 RNASEQ - update ini downgrade R version for differential expression and java for rnaseqc - BFXDEV-490
       32d1d5bd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4c43bb27 chipseq - update module -  BFXDEV-490
       7a681602 dnaseq - add correct dependency in metrics snv - BFXDEV-508
       7e2f61b3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4a95aa12 DNASEQ - remove bad ppn settings for mammouth ini file - BFXDEV-490
       58f9477e BFXDEV-465 - correct GATK DoC small bug to support when no bed fileis given in the readset file
       269b60fb DNAseq - Adding variant recalibration BFXDEV-436
       557f75ad DNAseq - gatk DoC will use the bed file as intervalls if the bed file is in the readset shett  - BFXDEV-465
       ce3014ac DNAseq - implement haplotype caller and mpileup annotation and filtering using old foinction as background - BFXDEV-463
       2aa00088 DNAseq - create new pipeline steps - BFXDEV-463
       f46876ac DNAseq - Add gvcf Combining the set of sample and genotyping - BFXDEV-440
       157b7e9d DNAseq - starting to implement GATK gvcf merging - BFXDEV-440
       3c407b26 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       74c6c2a3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66cfdc68 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcfdd2c8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c83d71ca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1b174115 bump pacbio module to patch 4 -  BFXDEV-415

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      105 commits

       db54b51a resolving conflict
       e7dc2516 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       3a2c22fa resolving conflict
       c554b5dc resolving conflict
       93a55f8e RNAseqDN - update module and remove trinity version check - BFXDEV-490
       48fafb90 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       2be21e18 Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       4103b9fb DNASEQ - update SnpEff command for versuion 4.2 - BFXDEV-490
       8e1bbc27 DNASEQ - support single end libray for metrics and dnasample metrics steps & allow filter_nstrech to use non recalibrated compress data - BFXDEV-503 BFXDEV-505
       4e0cd6cb resolving conflict
       228c23af resolving conflict
       6e42b3df resolving conflict
       fea21f4a RNAseqDN - update module and remove trinity version check - BFXDEV-490
       9cccf23b RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       50f9346f Rnaseq_denovo_assembly - correct single library type bugs  - BFXDEV-405
       71c233cc DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       3a3c119a resolving merging conflict
       a0491a4f RNAseqDN - update module and remove trinity version check - BFXDEV-490
       ca27f3f1 RNAseq - add specific older version of java for rnaseqc (support only 1.7)
       5b7e920f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7481dc79 Rnaseq_denovo_assembly - correct single library type bugs - BFXDEV-405
       b041650c DNASEQ - allow filter_nstrech to use non recalibrated compress data && update SnpEff command for versuion 4.2 - BFXDEV-505 ; BFXDEV-490
       132e479e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       ffcee7d8 DNASEQ - support single end libray for metrics and dnasample metrics steps - BFXDEV-503
       3f22b9d9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       10e4a0a0 updating modules (test on abacus human done) - BFXDEV-490
       5de58e20 rnaseq - support fr-stranded single end library in wiggle tracks and metrics - BFXDEV-499
       106788a2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9b97e885 rnaseq - make  wiggle tracks working for non UCSC genome ref - BFXDEV-498
       9bb4b278 correct the version of pandoc - BFXDEV-490
       e3f47cd7 update chipseq ini to latest version of modules - BFXDEV-490
       ae112c0f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f7b90b57 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f0e0ee57 update human genome sh files - BFXDEV490
       e3c63252 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       acc1d01e module - remove type in macs2 that make the module file invisible to the module system - BFXDEV-490
       4c413a26 module - add gemini and bwakit modules - BFXDEV-490
       c9a2f313 module - add star indexes for read length 75 and 150 bp - BFXDEV-490
       14464b16 module - adding a step print for the module creation - BFXDEV-490
       39f89bcd module - add gemini and bwakit ; update bcftools htslib java picard prinseq-lite rnaseqc samtools - BFXDEV-490
       b50d111d BFXDEV-490 - updating modules gnuplot ; hmmer; igvtools; macs2; pandoc; python_lib
       61770a48 fixing conflict between master and highCoverageVariants branches
       a3491eb1 BFXDEV-490 - updating modules bedtools; blast; gatk; ucsc
       bd023c2a bfx/gatk.py - removing conflict to allow gatk_variants branch to be merged
       292d51ca Ressource - bump vcftools version to 0.1.14 - BFXDEV-489
       13aa9255 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5b8a705c genomes - add rat Rnor_6 - BFXDEV-487
       9ef56ce6 genomes - install genome by default in dev and in prod if argument MUGQIC_INSTALL_HOME is given - BFXDEV-485
       4ac216dd ressource - fix version of matplot lib to 1.4.3 for Qiime compatibility - BFXDEV-483
       28175da7 gatk_variants - danseq -  add mark_duplicates cleaning files - BFXDEV-471
       1c06f0c0 gatk_variant - module - bump mugqic_tools install script to 2.1.3 - BFXDEV-484
       62c067e5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       da7245bf ressource - update python to 2.7.10 and add scikit-bio lib - BFXDEV-483
       7c774e6a remove tumor_pair from highCoverage
       5a20f1fb DNAseq - add cleaning list to the new steps - BFXDEV-471
       9266e030 DNAseq - change destionation path in the copy of SNV metrics zip file - BFXDEV-473
       71b439e9 Pacbio assembly - modify whitelist option usage to do not generate additional ini and xml - BFXDEV-456
       1f4eb9d2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be622168 Bump  module to the new smartanalysis patch 5
       15b99f43 gatk_variants - add baits intervals file for picard_hs_metrics - BFXDEV-467
       ad0a5c9c ingstats.py - Adding Array callrate check (Missing or Low) and add support for different manifest version through the -m option - BFXDEV-445 - BFXDEV-466 - BFXDEV-451
       75379ef3 DNAseq - synchronize code and ini and correct runnning issue - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       7a529056 DNAseq - remove some smalls dev typos and bugs - BFXDEV-436 - BFXDEV-440 - BFXDEV-463 - BFXDEV-465
       a19894cf PacBio - add require=False to the whitelist param - BFXDEV-456
       e8223e65 PacBio - incorpore whitelist option during intial filtering - BFXDEV-456
       53b2fbf1 PacBio - add addtionnal filtering config xml file for  whitelisted assembly - BFXDEV-456
       064ae2e6 PacBio - add specific ini file for whitelisted assembly - BFXDEV-456
       8e298ef4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       7769b771 update bwa install script to the newer version
       72613730 change star index place in resources/genomes/install_genome.sh
       0904455e update install_genome.sh to use the coirrect version of R_packages
       a3e03401 bump bos_taurus install script to ensembl 80
       c57be721 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       94223ff1 bump bos_taurus install script to ensembl 81
       17bea25e make ignstats poiting on nanuq (not QC nanuq)
       382c9528 Merge remote-tracking branch 'origin/processingReport'
       051e3782 RNAseq - repair missing dependency - BFXDEV-427
       26d5d41e debug paired_tumor scalpel vcf merging
       d2c011b8 removed merge conflicts
       29a319b5 DNAseq - regenarate updated README.md
       85017c1b DNAseq - change fixemate description text: picard -> bvatools
       ba184ccd High_coverage - needs dict2BEDs.py which is only in dev version of mugqic_tools - modifiy the ini file
       3d887984 DNAseq - change post fixmate sorting to generate and indexed bam - BFXDEV-421
       83e7e161 modify igvtools lib to allow modifiy paramter through the ini file - BFXDEV-419
       d5d86bf8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       71c3d26f DNAseq - change coverage metrics report to remove the CCDS part of the table - BFXDEV-381
       0de4e656 rnaseq add bam_hard_clip step : Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature - BFXDEV-377 + change rRNA metrics output to avoid name conflict - BFXDEV-401 - BFXDEV-407
       8a636810 update python and python lib to include Qiime ligth version installation - BFXDEV-412
       8bb30eea bump trinity install script to version 2.0.6
       1b2a34db rnaseq - remove tophat/botwie commented section and modules
       f53ed32c rnaseq_denovo_assembly - correct single library issue and add in comment the correponding change in the base.ini - BFXDEV-405
       d37b2c89 nanuq2mugqic_pipelines.py - support nanuq group info for Miseq/Hiseq/Pacbio technology - BFXDEV-418
       b03537ad update install module general script
       2d1dbe48 bump smrtanalysis module to 2.3.0 patch 4
       27bad469 update install module general script
       006f59cf update install module general script
       08acf20a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       bcb95896 bump smrtanalysis module to 2.3.0 patch 4
       ba4b9f14 update bvatools module for release 1.6
       c72acd4b core/pipeline.py remove duplicates sections in reports - BFXDEV-388
       4cb0ba7a ChipSeq - remove expected input files for broad peaks when generating homer annotation report - BFXDEV-393
       797c9681 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       3a54c856 modify resources/modules/smrtanalysis.sh for patch 3
       ce104cd1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d309b7da RESSOURCES - R_bioconductor - Addthe package 'circlize' to the list of R package to automatically install

  mathieu bourgey <mathieu.bourgey@mcgill.ca>      6 commits

       5db43a8b Merged in irp_end_fastq_notification (pull request #14)
       de7d7f84 Merged in irp_bcl2fastq2 (pull request #13)
       d3acddd6 Merged in illumina_run_processing_sprint (pull request #11)
       82747c86 Merged in pacBio_whitelist (pull request #7)
       e16e90df Merged in ignstats (pull request #8)
       8e2cf4c4 rnaseq.base.ini increase star io limit to 4G by default

  mmichaud <marc.michaud@mail.mcgill.ca>      62 commits

       4d820a7c BFXDEV-504 - Notify nanuq when fastqs are completed for a lane (fix url) - illuminaRunProcessing
       7ab569ee BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       cdafb239 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       fd8579e9 BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       7282b301 BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       2d9a8965 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       e0be3318 Fix alignment for lanes with a mix of samples with genomic database and no genomic database.
       6e04c629 BFXDEV-504 Notify nanuq when fastqs are completed for a lane (fix url)
       34ca19a2 BFXDEV-504 Notify nanuq when fastqs are completed for a lane
       09ff7003 Increase resources for qc.
       3dfce03b BFXDEV-488 Illumina Run Processing: Use bcl2fastq2. Fix for dual-indexing.
       5e8133ee Use 2 core for rsync.
       48a4bd86 Merge branch 'master' into irp_bcl2fastq2
       a9560d4f Fine tune cluster settings according to the history of the jobs of the last three months.
       b917ee55 BFXDEV-413 Fix code according to code-review.
       a8fd8bfe BFXDEV-488 Illumina Run Processing: Use bcl2fastq2.
       be6ec3f1 Merge branch 'master' into illumina_run_processing_sprint
       e0c3e20b BFXDEV-468 VerifyBamId: Fix job name when not using bed file.
       ec19be9d BFXDEV-468 VerifyBamId: Run even we don't have a bed file.
       1ffafbb6 BFXDEV-468 VerifyBamId: Use a version supporting "chr" in chromosome name.
       dec09d84 BFXDEV-468 VerifyBamId: Fix output for nanuq.
       eed1cc16 BFXDEV-468 VerifyBamId: Shorter job name.
       750495f7 BFXDEV-468 VerifyBamId: Properly concat jobs
       ce01a76f BFXDEV-468 VerifyBamId: Nanuq friendly output.
       9d5fa0e4 BFXDEV-468 VerifyBamId: Don't use a compressed VCF + Nanuq friendly output.
       c4bd81b1 BFXDEV-468 VerifyBamId: Add missing configuration.
       72c77d23 BFXDEV-468 VerifyBamId: Fix annotation file name.
       00010997 BFXDEV-468 VerifyBamId: Fix wrong genome for sample.
       15a3a2bc BFXDEV-468 VerifyBamId: Changes to the annotation filename.
       7fc91bd9 BFXDEV-468 Don't run verifyBamId when there is no "dbnfp_af_field" on the genome.
       d747ef6f Illumina Run Processing: Fix barcode counter path with cvmfs
       68c00893 BFXDEV-468 IRP - VerifyBamId: Optional dbsnp_version and dbnsfp_af_field in genome ini file.
       5a3fd82b Use cvmfs blast and bwa.
       75af739c BFXDEV-468 Add verifyBamID in Illumina Run Processing
       94ccde0c Merge branch 'master' into illumina_run_processing_sprint
       7d6cf517 BFXDEV-447 Add rRNA estimate using silva blast database. (Fix change output file name)
       dd9f3f7a BFXDEV-464 Don't use an hardcoded genome list for the alignment: parse the value from the sample sheet and validate that we have the genome on disk.
       015435ed BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name (fix redo of pipeline always restarting the STAR align).
       f0ecd639 BFXDEV-447 Add rRNA estimate using silva blast database. (Fix error when creating result file)
       2d39be9b BFXDEV-413 Add Estimated Genome Size from the nanuq CSV.
       c73b3fc9 BFXDEV-447 Add rRNA estimate using silva blast database.
       fdcb45ee BFXDEV-462 Output STAR bam in a unique folder to support lanes with multiple samples with the same name.
       1819bf2c Add sample sheets examples and add description about "Genomic Database"
       be10d09b BFXDEV-457 BFXDEV-459 Update documentation.
       b337a061 Don't warn when skipping alignment on a sample without "Genomic Database"
       a73ca70a BFXDEV-459 Merge start_copy and copy steps.
       583237b9 BFXDEV-457 Add the possibility to regroup all the md5 jobs into one (add "one_job=1" in "[md5]" config section).
       56cdf2d3 BFXDEV-458 Fix mask calculation problem on lanes with nextera/truseq samples.
       041d4e81 More useful usage message.
       3cd065cc Merge branch 'master' into irp_genomic_database
       d90c137d Merge branch 'master' into irp_genomic_database
       2f0ef65a Illumina Run Processing: Increase memory for BVATool DoC and Readsqc
       a864f318 BFXDEV-369 Illumina Run Processing: Fix GRCm38 genome detection.
       7eae59d7 BFXDEV-369 Illumina Run Processing: Use the reference specified in the request form
       64fb218e Increase processors per node for Fastq and bwa alignment. Ignore setfacl exit code for the copy step.
       fe445772 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0) - Fix passFail threshold according to the new value.
       cafa62c3 BFXDEV-417 IGNStats: Output Array Call Rate as ratio, not percentage (0.95 not 95.0%)
       eb810e06 BFXDEV-417 IGNStats: Output identity value as ratio, not percentage (0.95 not 95.0)
       37245835 BFXDEV-417 Modifiy IGN stats parser to allow the upload of data to a nanuq server.
       ce528549 Use 13 core for BWA job to avoid excessive memory usage.
       8505fab8 BFXDEV-385 Fix bed files handling in Illumina Run Processing.
       29d5c0f6 BFXDEV-380 Change permissions of the run files on the destination.

  noreply@clumeq.ca <libuser@lg-1r14-n01.guillimin.clumeq.ca>      15 commits

       ad55029f genomes - updated bash script for human genome GRCh37 and GRCh38  - BFXDEV-490
       4d653bf9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       d2fd863a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       e4b6339f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       1c79f7e8 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       36b260f6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       56550276 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       7e8870a0 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0a73d79f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       eb301072 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       8c32446d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0b9b70c6 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines
       0233d702 mugqic_tools is now up to date && resources/modules/dev/epacts.sh has been removed
       64c821d5 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a65de9e0 temporary modify install_genome.sh to run the job in batch (job submission is blocked

  ptranvan <patrick.tranvan@mail.mcgill.ca>      1 commits

       2a1917b2 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipelines

  robert.eveleigh@mcgill.ca <reveleig@abacus1.ferrier.genome.mcgill.ca>      4 commits

       f1849661 Addition of htslib module and resolved pull request comments
       457c3695 addition of vcf preprocessing, snpeff, gemini and ini adjustments
       09ad7f1f addition of vcf preprocessing, snpeff, gemini and ini adjustments
       05f05b22 Merge branch 'highCoverageVariants' of https://bitbucket.org/mugqic/mugqic_pipelines

  Robert Eveleigh <robert.eveleigh@mcgill.ca>      2 commits

       60fd2305 Merge branch 'highCoverageVariants' of bitbucket.org:mugqic/mugqic_pipelines into highCoverageVariants
       76802ee0 Updates to pairedTumor: addition of CombineVariants

2.1.1        Mon Apr 13 22:23:46 2015 -0400        172 commits

  Francois Lefebvre <lefebvrf@gmail.com>      4 commits

       65710dc8 QUAST and Minia dev install scripts
       5b404f38 Updated kmergenie version in install script and moved to dev
       6db368e9 Added NxTrim mugqic_dev install script
       45b9fb07 Updated a bunch of module dev install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      109 commits

       4db053f6 Added Gorilla_gorilla.gorGor3.sh in resources/genomes/old/
       deb0b822 Changed report job names with '_' instead of '.' to avoid scheduler cpu setting > 1
       4d8b38a7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b471f6f4 Fixed bug in RNA-Seq De Novo Assembly : use RSEM plugin in Trinity instead of external one
       6c79dd4c README.md edited online with Bitbucket
       2e53a765 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       455a50c6 Fixed missing pandoc module in rnaseq
       0ffb7ead Added variables total and average read length in pacbio assembly stats report + changed locale from fr_FR to en_CA
       70150e79 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       121faa00 Added sample_cutoff_mer_size to pacbio report job names
       fcb6b897 README.md edited online with Bitbucket
       272830be README.md edited online with Bitbucket
       89c55cbc README.md edited online with Bitbucket
       8f6eeb79 Fixed cpu bug for blastp_transdecoder_uniprot job in RNA-Seq De Novo assembly
       cae10e5d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       28eb9540 Minor comment updates
       2c4a463a README.md edited online with Bitbucket
       c4262fa9 README.md edited online with Bitbucket
       3f664a57 README.md edited online with Bitbucket
       dc8ec12a README.md edited online with Bitbucket
       c60d7424 README.md edited online with Bitbucket
       b38e1950 README.md edited online with Bitbucket
       6d4206c1 README.md edited online with Bitbucket
       4a6982c6 README.md edited online with Bitbucket
       39ee0b9b README.md edited online with Bitbucket
       1f13501a README.md edited online with Bitbucket
       0376786e README.md edited online with Bitbucket
       43d0ac6d README.md edited online with Bitbucket
       be50cd0a README.md edited online with Bitbucket
       f85506e8 README.md edited online with Bitbucket
       5e335b49 Added GNU Lesser General Public License (LGPL) to MUGQIC Pipelines
       0b79b46b  BFXDEV-59 Updated ChIP-Seq report redesign with sample metrics without trimming
       ec88ea54 Fixed bug cluster_cpu for blastx_trinity_uniprot
       5962a9fe Removed UniRef BLAST in RNA-Seq De Novo Assembly since it is too long; factorized differential expression code; renamed design variable into contrast
       611223e9 Adjusted Trinity butterfly cluster resources for RNA-Seq De Novo Assembly on abacus
       371d60c6  BFXDEV-59 Completed ChIP-Seq report redesign
       7eedd054 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       01a41aac  BFXDEV-59 ChIP-Seq report redesign
       f8bf89fa BFXDEV-59 Completed PacBio Assembly report redesign
       479ea0b0 BFXDEV-59 Differential expression RNA-Seq De Novo Assembly report redesign commit
       76c9bc36 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       c7a9036f BFXDEV-59 More RNA-Seq De Novo Assembly report redesign commit
       ac2daf38 README.md edited online with Bitbucket
       e7f156a5 README.md edited online with Bitbucket
       3caf504b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       6771f039 BFXDEV-59 First RNA-Seq De Novo Assembly report redesign commit
       5c858ad9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       5d82ab02 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       94e9930d Decreased rnaseq cufflinks default pmem to 2700 for guillimin
       79ea01b1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       80cfde80 Increased rnaseq cufflinks default pmem to 3700 for guillimin
       9e655b98 Fixed BED file abspath in bvatools
       c183fade BFXDEV-59 Completed RNA-Seq report redesign
       61676ebb BFXDEV-59 Added metrics steps for RNA-Seq report
       e760ac7e BFXDEV-59 Fixed merge conflict
       fbedae5e BFXDEV-59 More and more commit for partial HTML report
       7a95c844 Increased dnaseq compute_effects ram
       c6c7ead0 BFXDEV-59 Even more commit for RNA-Seq report
       75b2a7ee BFXDEV-59 Even more commit for RNA-Seq report
       5da98859 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       0463d009 BFXDEV-348 Fixed with readset name = <sample>.<library_barcode>.<run>.<lane>
       3072da14 BFXDEV-59 Even more commit for RNA-Seq report
       57f2aef6 BFXDEV-59 More commit for RNA-Seq report
       ad7e670e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       be924071 Updated mugqic_tools to version 2.1.0
       672b6814 BFXDEV-59 First commit for RNA-Seq report
       f194ec39 BFXDEV-59 DNA-Seq report minor fix
       7f98cbfe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       40a04d08 BFXDEV-59 DNA-Seq report redesign done + fix
       b4b408f7 BFXDEV-59 DNA-Seq report redesign done
       426e8788 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into report
       a0b4223d Minor doc fix
       bc1d662e BFXDEV-59 Even more report redesign commit
       09de8ac4 BFXDEV-59 More report redesign commit
       91e8e95b BFXDEV-59 First report redesign commit
       77f66c40 Added UCSC genomes in install_all_genomes.sh
       aa6a8e0a Create genome rrna.fa with grep -i 'rrna'... + remove variation sequences containing '.' in Ensembl vcf.gz which make GATK crash
       525faa86 BFXDEV-295 Updated mugqic_R_packages to 1.0.3 for RNA-Seq De Novo pipeline
       da8253e4 BFXDEV-295 minor fix for RNA-Seq De Novo pipeline
       e161ad93 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4d177beb Updated Python pysam to version 0.8.2
       e7ca27f5 README.md edited online with Bitbucket
       48073caf README.md edited online with Bitbucket
       70a0304b README.md edited online with Bitbucket
       1e24257c Increased cores for homer_annotate_peaks in chipseq
       3c2c24de Fix new section names for blast on uniprot in RNA-Seq De Novo pipeline
       e699c642 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       aeedee41 Added chipseq in PATH of module mugqic_pipelines
       7ce201c6 Added ccds filepath in Rattus_norvegicus.Rnor_5.0.ini
       7b08143b BFXDEV-295 Moved trinotate Pfam + UniProt DBs in /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod/genomes/[blast|pfam]_db/ for RNA-Seq De Novo pipeline
       c47a5190 BFXDEV-295 Update RNA-Seq De Novo pipeline with modules in prod and trinotate updated
       33812be2 Minor mammouth adjustment regarding increased ram and core in snp_effect job for DNA-Seq Pipeline
       9689e6db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c345049c Increased ram and core in snp_effect job for DNA-Seq Pipeline
       8743bbb9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f5b4229 BFXDEV-295 Updated RNA-Seq De Novo assembly pipeline with Trinity 2.0 and Trinotate 2.0
       b484ff7f README.md edited online with Bitbucket
       4b19c200 README.md edited online with Bitbucket
       1b0f09cd Minor fix in python lib install
       91806124 Separated python install and python lib install + minor weblogo install update
       e06560dd Minor fix in Perl lib install
       b202849f Fixed missing bvatools.depth_of_coverage other_options + snpsift_annotate module_snpeff=mugqic/snpEff/4.0 + minor uppercased '.insert_size_Histogram.pdf' for picard.collect_multiple_metrics output file in dnaseq pipeline
       f8b35eb3 Fixed missing 'nodes=1' in pacbio_assembly.guillimin.ini smrtanalysis_run_ca
       e7d4a0e1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       62a0e005 Updated bedtools to version 2.22.1
       9328bc99 Standardized perl module install and moved CPAN libs in a different script
       0fb4c648 BFXDEV-335 Removed adapter FASTA files except adapters-fluidigm.fa
       4aba5c16 BFXDEV-335 Create adapter FASTA from readset file for Trimmomatic, if not defined in config file
       ba787480 Version bump to 2.1.1-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      14 commits

       1792f2c8 Version bump to 2.1.1
       88108220 BFXDEV-375 Fixed ram sorting issues when using star
       e9fcf543 Merge branch 'master' into rna_metrics
       a7e3465c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ba66620 BFXDEV-368 IGN script to extract stats
       8f691400 Fixed code for tools that need to be downloaded manually like gatk
       f6ee04eb Updated gatk
       4622df3b Updated samtools, bcftools, htslib
       6eac3158 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a44561f3 BFXDEV-351 removed md5 from markdup and added it to recal
       8abc6ccf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1db78f98 BFXDEV-355 Removed CCDS from the options BFXDEV-356 added compression to haplotype caller output
       9a0a6984 BFXDEV-351 Removed MD5 from markdup, added it to recalibration
       b33195c0 BFXDEV-346 Split jobs in a more uniform way

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       70099f3f rnaseq - include correlation matrix manual estimation using an utils.py new function
       fd538018 RNAmetrics - include ini_section argument to some bfs picard function to allow to use them several time with different parameter in the same ini
       e7c5a010 rnaMetrics - remove confilct pipelines/rnaseq/rnaseq.base.ini
       9ad7fbb1 rnaMetrics - update pipelines/rnaseq/rnaseq.py pipelines/rnaseq/rnaseq.base.ini
       9d3b87da RNAseq - metrics update ini
       1913bc10 RNAseq - remove conflict
       deca058d RNAseq - metrics RNA - update base ini
       2fc2006b RNAseq - remove rnaseqc; add picard_rna_metrics ; partial add estimate_ribosomal_rna - BFXDEV-345

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       05937d60 RNA-seq - rna_metrics : test are ok; new files annoation files are created ; point to the production assembly folder -  BFXDEV-345
       b2c948b8 ressource -Genome install script add the creation of the ref_flat file format of the annotation from the gtf and correct a path in the GRCh38 file (All.vcf.gz) - BFXDEV-374
       7d8f72fc RNAseq- update module version not found in CVMFS (bowtie, bvatools, mugqic_tools) - BFXDEV-373
       e57ed5ed RNAseq- remove redundant step in the step initialization - BFXDEV-371
       3138867b release new version of resources/modules/mugqic_tools.sh
       9aa58f7a RNAseq - rna_metrics - finish the rRNA and rna metrcis modifications - BFXDEV-345
       cd5f2391 RNAseq -rna_metrics -fit the new bam2fq from bvatools_dev
       d2048ebe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into rna_metrics
       d651feeb RNAseq - rRNA metrics add bvatools | bwa | picard && rrnaBMAcount.py

  mmichaud <marc.michaud@mail.mcgill.ca>      27 commits

       2a821a4b bcl2fastq module name change.
       00b481f5 Merge branch 'master' into irp_rna_metrics
       c2c33a3f MPS-1740 Use updated production genomes.
       5f3c74e5 MPS-1740 Use released version 1.5 of bvaTools
       d191969e Increase STAR sort memory
       dd936d0a BFXDEV-363 IRP: Don't copy phasing and matrix files.
       ae9713bc BFXDEV-353 IRP: Standardize job name
       e83e25fa Merge branch 'master' into irp_rna_metrics
       87aaf0db Merging rna_metrics on irp-rna_metrics
       b0dff3b1 BFXDEV-353 Use dev version of mugqic tools
       f94c1357 Merge branch 'master' into irp_rna_metrics
       75398a38 BFXDEV-353 Use new version of rRNABAMcounter.
       5f443d4b BFXDEV-353 Use a nanuq friendly name for the rRNA metrics file
       ea661af4 BFXDEV-353 Add the rnaseqc '-rRNA' option to set an optional ribosomal rna interval list. Setting an empty file will skip the rRNA count.
       91c0d87d Simplify illumina run processing RG tag logic. As library, run and lane are mandatory there is no need to validate that they exists.
       4abaa7b1 BFXDEV-353 Fix bwa other_options
       81e7423b BFXDEV-353 Use DEV versions of bvatools and genome. Update mugqic_tools to 2.1.0.
       79a02a72 BFXDEV-353 Update rna-seq metrics according to rna-seq pipeline
       081c6b81 Merge branch 'master' into irp_rna_metrics
       fcc8392b Code format
       d62f8550 Merge branch 'master' into irp_rna_metrics
       0bc4cbe7 Add missing thread parameter for bvatools_depth_of_coverage.
       221eaae0 BFXDEV-339 Use rrna file in rnaseqc
       8e95c994 BFXDEV-339 Use rrna file as ribosomal annotation (instead of ncrna)
       5efa1042 BFXDEV-338 New Nanuq MPS run association
       7aab26f9 BFXDEV-338 New Nanuq MPS run association
       6e6de075 BFXDEV-338 Add run_id in all commands available parameters

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       559973d8 correction other_options

2.1.0        Wed Feb 4 14:40:16 2015 -0500        123 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       eefcc474 Added qualimap installa script + updated population
       d9397b01 Changes to sailfisj , R install scripts

  Joël Fillon <joel.fillon@mcgill.ca>      104 commits

       09344ee6 Version bump to 2.1.0
       b3e85618 Changed procs= to nodes=1:ppn= in rnaseq and rnaseq de novo guillimin config
       a63500db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       a195837e Added libgd in PacBio mammouth config
       a863cdd5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       e13ac413 Added pacbio_assembly.mammouth.ini + removed unnecessary wgs module
       1cb2ea5d Added pacbio_assembly.guillimin.ini + adjusted guillimin cluster settings
       b9f12642 README.md edited online with Bitbucket
       ff58b26e Removed 'daemon' from scheduler options
       2d4459d7 BFXDEV-292 Removed optional trimming dependencies for report in chipseq pipeline
       6b9cdad5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       f40f39ac BFXDEV-292 Fix metrics and report dependencies in chipseq pipeline
       c1fb3245 README.md edited online with Bitbucket
       d8ade06f Removed tmp_dir validation since some directories are available on exec nodes but not on login nodes
       44c931e0 Removed memtime from PacBio pipeline + updated smrtanalysis to version 2.3.0 patch 2
       698d3f70 BFXDEV-292 Added cleaning in chipseq pipeline
       cac63055 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       086ee05e Updated smrtanalysis module 2.3.0 with patch 2
       0c6d970f README.md edited online with Bitbucket
       f808abe0 Updated ChIP-Seq README + main READMEs
       e789ffd6 Added BiSNP module install script
       26a06293 Update mugqic_tools to 2.0.3 in rnaseq config
       7af099da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d10a0584 BFXDEV-292 Fix metrics tmp file preventing job to be up to date + module ImageMagick to use convert command on mammouth
       367f86eb Removed '\' before /ltmp/fillon in mammouth config files
       0c1988f0 Added picard tmp_dir type validation
       cf820656 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0884bf53 BFXDEV-292 Added config files for all clusters in chipseq pipeline
       2a572152 Updated mugqic_tools install to version 2.0.3
       387c225b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       5017938f BFXDEV-292 Last steps implementation in chipseq pipeline
       52d1f782 BFXDEV-292 Added UCSC genomes hg19, mm9, mm10, rn5
       751e701a BFXDEV-319 Adjusted cluster settings in rnaseq.mammouth.ini + .ini absolute path for report
       edc26514 BFXDEV-35 Standardized memtime module install
       90f4aa5c Minor update in README-release.txt
       ab65f5da Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       f37bed0a BFXDEV-292 Set peak annotation with internal homer genome in chipseq pipeline (incomplete)
       9b356028 BFXDEV-35 Moved wgs-assembler module install in dev
       77859cc2 BFXDEV-35 Standardized vcftool module install
       a07ea0b5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       2ee43154 BFXDEV-35 Standardized ucsc module install
       284138b5 BFXDEV-35 Standardized trinotate module install
       9f31173a BFXDEV-35 Standardized trimmomatic module install
       1fab82a1 BFXDEV-35 Standardized tophat module install
       00b9355b BFXDEV-35 Standardized tmhmm module install
       2917504d BFXDEV-35 Standardized tabix module install
       6c1ec5b5 BFXDEV-35 Standardized star module install
       55547462 BFXDEV-35 Standardized snpEff module install
       25f53a1c BFXDEV-35 Fix smrtanalysis archive exec permission
       275894fb BFXDEV-35 Standardized smrtanalysis module install
       6f3b27ab BFXDEV-35 Standardized signalp module install
       a335bfa5 BFXDEV-35 Standardized samtools module install
       7cfda90d BFXDEV-319 Added rnaseq.mammouth.ini
       dac40de2 BFXDEV-319 Removed java option -Dsamjdk.use_async_io=true in all .ini files
       e394bf78 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       705f2c6b BFXDEV-292 Adjusted cluster settings in chipseq pipeline
       a87203d3 README.md edited online with Bitbucket
       2524cb9d README.md edited online with Bitbucket
       af39fb42 README.md edited online with Bitbucket
       e180a838 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       0daae6ed Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       93c1a288 Fixed bam sub arguments in dnaseq bwa_mem_picard_sort_sam
       934b42f6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       a5931e42 BFXDEV-35 Standardized rsem module install
       dbaebd1e Moved repeatmasker module install script in dev
       0c30ff42 BFXDEV-292 Minor docstring change in chipseq pipeline
       dc0b9888 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ef6c09e4 BFXDEV-292 Added GTF for homer_annotate_peaks in chipseq pipeline
       26894dcc README.md edited online with Bitbucket
       c7cdd1ae Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       4f0b0858 BFXDEV-292 Fixed input/output files bugs in chipseq pipeline
       4c1cf7d2 BFXDEV-35 Standardized weblogo module install
       bf36cc0e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       7952a5fb BFXDEV-292 homer_find_motifs_genome and annotation_graphs steps in chipseq pipeline
       b29cb3d9 BFXDEV-35 Standardized rnaseqc module install
       05208ed0 BFXDEV-35 Standardized rnammer module install
       e33776a5 BFXDEV-35 Standardized more prinseq-lite module install
       fcc8bc9e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d15d9b00 BFXDEV-35 Minor fix in AMOS module
       07bba1e8 BFXDEV-35 Fixed path bugs in MUMmer and AMOS module install
       99363eed BFXDEV-35 Standardized MUMmer module install
       423c3c7c BFXDEV-35 Minor fix in module install template
       c00bdfcf BFXDEV-35 Standardized java module install
       d185a507 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       8bed5d97 BFXDEV-292 homer_annotate_peaks step in chipseq pipeline
       0bca1960 BFXDEV-35 Standardized igvtools module install
       1581fc97 BFXDEV-35 Standardized hmmer module install
       2acdc492 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       cddbce73 BFXDEV-292 More macs callpeak in chipseq pipeline
       ed7cf4ec Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b2d02b10 Fixed macs2 with generic shebang
       cf61d7e0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       bbe5703e BFXDEV-292 Fixed missing samtools module for MACS2
       ff1a4ec9 Standardized gnuplot module install
       dcdc2114 Standardized exonerate module install
       5c477462 Minor change in README-release.txt
       076c02b6 Version bump to 2.1.0-beta
       953894cb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       1090d383 More chipseq deelopment
       0841b02b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       ace0e308 BFXDEV-292 Beginning of macs2_callpeak step in chipseq
       b404e222 Added qc_plots_R step in chipseq
       43e70dce Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines into chipseq
       51f34444 BFXDEV-292 First draft of chipseq pipeline

  lletourn <louis.letourneau@mail.mcgill.ca>      3 commits

       9f9f60cf Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d2d65f01 BFXDEV-327 Used only one thread for haplotypecaller because of a race condition
       639f6502 BFXDEV-327 Used only one thread for haplotypecaller because of a race condition

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       12ff20e9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       233244b4 RNAseq - correct stdin input issue of htseq-count - BFXDEV-318
       fa335288 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       47c4cc81 RNAseq - htse-count: pipe samtools view -F 4 output in htseq-count instead of using the bam to remove error due to unmapped reads - BFXDEV-318

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       7ec03449 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1905eb60 COMMON - add fastq2 = None in sam_to_fastq pipeline step wehen the read are single - BFXDEV-321

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       953de3d8 BFXDEV-332 Fix depth_of_coverage 'other_options' by removing extra quotes
       abca1808 BFXDEV-332 Base the jobs configuration on the perl pipeline
       619c0760 BFXDEV-329 Add run and lane number to the alignment and metrics jobs name
       6f4440f3 BFXDEV-202 Increase RAM for barcode counting
       eac645f2 BFXDEV-328 Increase mem walltime to 48h (and set ram for SortSam)
       ccb73f33 BFXDEV-331 Run Processing: Use a RG id that is unique across multiple lanes
       c2218daf BFXDEV-317 Run processing: Use the old name for the coverage and onTarget metrics
       11ec9e77 BFXDEV-316 Fix errors when using a custom Illumina sheet file

2.0.2        Mon Jan 12 16:56:16 2015 -0500        80 commits

  Joël Fillon <joel.fillon@mcgill.ca>      54 commits

       104778ab Version bump to 2.0.2
       bbdf438c Version bump to 2.0.2
       9d984cec Fixed cluster resources in rnaseq_denovo_assembly.mammouth.ini
       ee96fe3b Minor fix in pacbio_tools_split_reads config section name
       8831a3ae Even more standardization of module install
       a2c58248 More standardization of module install
       2b334e01 Standardized mugqic_pipelines module install
       395d3ec6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       66230bd2 Standardized mugqic_tools module install
       d9a56ab7 Standardized trinity module install
       0549be9d Standardized cufflinks module install
       ddbd8f2e Standardized MACS2 module install
       386921ec Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9ba8d390 Standardized picard module install
       55a5ccb7 Standardized cd-hit module install
       d693b9a2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9f05fe1b BFXDEV-310 Updated picard to version 1.123
       a118e0d9 Standardized bwa module install
       4dc7bb7a Updated bvatools to version 1.4 for dnaseq and puure
       cca258b1 Standardized bvatools module install
       94910f99 Standardized bowtie2 module install
       3e8f2016 Standardized bowtie module install
       1b70a257 In picard_sam_to_fastq, updated skip test if FASTQ column present and BAM missing
       345909af BFXDEV-290 Added job_input_files, job_output_files in JSON export
       f6e623e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       39e9ebfa More JSON export for Daemon Scheduler
       72956719 Standardized blast module install
       947111a8 Standardized bedtools module install
       9325049b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       c52d2775 Standardized GATK module install
       a69bd820 More module install generalization; first test with prinseq-lite
       d116ebd6 Minor wget output file fix in prinseq-lite.sh and module install template
       c30ad665 Added prinseq-lite module install script
       b2059f5f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       8ca90177 Fixed pb wget --spider in install_genome.sh + various up-to-date checks in Homo sapiens assemblies install
       21e59af0 BFXDEV-284 In nanuq2mugqic_pipelines.py, use native httplib, retrieve BED files and create adapters file
       e9620ff0 BFXDEV-290 First draft of daemon scheduler
       b778e5db BFXDEV-74 Finished cleaning in rnaseq denovo assmbly pipeline
       a171581e BFXDEV-74 Started cleaning in rnaseq denovo assmbly pipeline
       c1cee295 Reorganised lib functions more compact
       1cc8e33b BFXDEV-74 Added cleaning in rnaseq pipeline
       13475e2a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       84954f32 BFXDEV-74 Added cleaning for dnaseq pipeline + reformatted various bfx modules
       1948d9d6 README.md edited online with Bitbucket
       848a0abb Fixed bug use lstat instead of stat to check job up-to-date status without following symlinks
       fdfcbdde Minor aesthetic updates on pipelines doctrings + README.md
       da66a410 Generated pipelines README.md from docstrings using --help
       89718b0a Added steps docstrings in RNA-Seq De Novo Assembly Pipeline
       ab5e595d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       232bed51 BFXDEV-296 Added steps docstring for RNA-Seq pipeline
       8d6f9a08 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b665c4cc Minor release doc update
       bf9afb6c Update mugqic_pipelines module to 2.0.1
       530edcf4 Version bump to 2.1.0-beta

  lletourn <louis.letourneau@mail.mcgill.ca>      1 commits

       ccb812e8 Fixed picard installer, reverted back to 123

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       fcd7134d pull before pushing Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       05b352d8  PacbioQC - Fix code vs ini section name - BFXDEV-315
       d0852690 RNAseq - correct discrepency in hsteq_count ini calls - link to BFXDEV-312
       6773ef75 pull before pushing
       d5d82f0d RNAseq - correct errounous section header in ini files - BFXDEV-312
       98a456ab Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fdd92aa2 RNAseq -update guillimin cluster ini file - BFXDEV-307
       b4a12a5a RNAseq - fix report dependencies - BFXDEV-306

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      9 commits

       693ea1dc PacBIo - add input files to the first pacBio job file (cp job) - BFXDEV-305
       ce7f0b36 Revert "PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305"
       fc6eb26f PacBIo - use the readset file as fisrt input file (cp job) - BFXDEV-305
       162af5e5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       2ff273be RNAseq - modify docstrings - BFXDEV-303
       c58b6d3d DNAseq - correct report  dependency (typo) - BFXDEV-304
       ab81bc34 DNAseq - correct report  dependency - BFXDEV-304
       12830d32 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       b6f4e8fd PacBio de novo - to allows running multiples sampes in parralele in the same analysis - BFXDEV-301

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c14a2818 BFXDEV-314 Illumina Run Processing: Output nanuq friendly metrics files
       3141cccc BFXDEV-313 Run processing: Don't depend on the start fastq notification
       5956b717 BFXDEV-310 Use Picard 1.123 with a bug fix in MarkDuplicate
       d0188793 BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       793e91cc BFXDEV-308 Fix wrong index in file name when the index is truncated in the sample sheet
       000d043e BFXDEV-309 Use the hiseq qos by default on abacus
       91cfb262 Run processing: Check the library type when the library source is 'library' to determine if the sample is from RNA
       94f196eb Run processing: Fix configuration for miseq (copy's destination folder)

2.0.1        Wed Dec 17 09:56:23 2014 -0500        33 commits

  Francois Lefebvre <lefebvrf@gmail.com>      3 commits

       76a44bc1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fe96c580 no message
       0f6a3ec6 Modified R_Bioconductor script to accommodate older gcc's

  Joël Fillon <joel.fillon@mcgill.ca>      22 commits

       bc72fc83 Version bump to 2.0.1
       5868b857 Updated mugqic_R_packages to 1.0.1 and mugqic_tools to 2.0.2
       3230102a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fd49fedd Updated /mugqic_tools.sh to 2.0.2
       2b979a2d Fixed dnaseq species_vcf_format_descriptor with absolute path (no /sb/programs/analyste)
       51fef9a1 BFXDEV-299 Fixed bug GATK realign with unmapped parameter not skipping other chromosomes
       6890d45e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       9a53a115 BFXDEV-296 Added DNA-Seq step docstring documentation
       a820bbc2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       4475a6b9 BFXDEV-296 Added detailed --help option with output in markdown format
       cb49e7b1 Added detailed --help option with output in markdown format
       bf7bb0ab More PacBio README.md update
       cd046510 README.md edited online with Bitbucket
       93e28c7c README.md edited online with Bitbucket
       48619e05 README.md edited online with Bitbucket
       2e79a551 README.md edited online with Bitbucket
       3a17272e PacBio Assembly README.md generated by --help
       93e81247 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       69073ee6 Trinity + normalization settings for guillimin in rnaseq_denovo_assembly pipeline
       8d8571f3 BFXDEV-26 Added cluster_max_jobs param in config files + added warning in PBS scheduler if this limit is reached
       14194c6e Minor update in README-release.txt
       38f24d73 Version bump to 2.1.0-beta

  mmichaud <marc.michaud@mail.mcgill.ca>      8 commits

       c661248d Run processing: Copy the rnaseqc metrics file to the bam directory to be compatible with nanuq
       e9ecd7ab BFXDEV-165 Add missing Perl module for bcltoFastq
       44c06a01 BFXDEV-296 Changes to run processing documentation to use the generic help
       ecb24473 Run processing: Put rnaseqc results in a specific folder by librairy
       336a2652 Run processing config: Only send email on abort and increase cluster wall-time to 48h for the fastq job
       3e8733a7 Run processing: Add symbolic link to STAR created bam to fix output file check on re-run. The original file is renammed to follow nanuq naming conventions
       7e8dddbb Run processing: Fix rnaseqc sample list
       7e38acff BFXDEV-297 Add RNA-SeQC in run processing pipeline

2.0.0        Thu Dec 11 18:12:54 2014 -0500        669 commits

  Eric Audemard <audemard@ip03.m>      1 commits

       dc40fe20 update mammouth .ini

  Eric Audemard <eaudemard@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       62e24647 fix bug to find lib in pipeline python

  Eric Audemard <eaudemard@lg-1r17-n02.guillimin.clumeq.ca>      1 commits

       02b7ea60 fixe bug on puure before Ray execution

  Eric Audemard <eaudemard@lg-1r17-n03.guillimin.clumeq.ca>      2 commits

       3165f9aa Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       bbe9d010 Bug fixed. Pipeline tested and validated on guillimin (step 1 to 21)

  Eric Audemard <eaudemard@lg-1r17-n04.guillimin.clumeq.ca>      2 commits

       cd0dfbcf Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       ad3c17ea bug fixed : 1) unmapped read in realign 2) small script error in puure

  eric audemard <eaudemar@imac6-ub.(none)>      15 commits

       3b2258df bug puure
       c0acd51b bug puure
       0a33731c bug puure
       c4ca9549 update base.ini for guillimin
       b47c7c9f add base.ini for guillimin
       ab7535ac 3nd script of PUURe done + some bug
       5e7c1fe3 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a0b5418c 2nd script of PUURe done + some bug
       9576d845 add file for puure pipeline
       844021a9 bug correction on 1st  script of PUURe
       c49f5583 first script of PUURe done
       89793044 add art software (http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) : Simulation Tools add perl lib path for SVDetect
       315bfc11 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_resources
       2b250c7a add script for install: SVMerge, RetroSeq, TIGRA, cmake
       f89034f1 add biopython

  eric.audemard@mail.mcgill.ca <eaudemar@abacus2.ferrier.genome.mcgill.ca>      4 commits

       caa90df1 add ini files. bug fixed on mkdir cov directory
       47714825 add ini files
       46f547e2 debuging step: fixe some bug
       00d5fb8a bug fixed on puure + add all output in dnaseq metrics (gatk, picard)

  Francois Lefebvre <lefebvr3@ip03.m>      1 commits

       e84abf37 picard version 1.125 module

  Francois Lefebvre <lefebvrf@gmail.com>      54 commits

       22fd75c5 Updated R install scripts to reflect mugqic_pipelineS new name and deprecation of mugqic_resources
       ccb88f59 fix to deploy script
       3e061c3f update R packages list in resources
       22d05454 Fixed bug when -v is specified
       cd63425d more cleanup
       a3e5bbe2 more cleanup
       2be18dfe no message
       8536bf30 Fixed bug in prefix mecanism
       2e83482e Fixed default behaviour of update mode
       2bb28121 R_Bioconductor.sh: update mode is default
       7b04a5f2 First commit of the new R install scripts
       6e908acc STAR install script update to 2.4.0e, accounting for new folder structure in the archive (bin/,source/, etc)
       39922d9b shortStack install script
       b032d1db tcltk is part of base install
       c80bb97c Updated star install dev script
       46888c05 magittr and other packages added to install list
       d8773a33 Guillimin before mammouth
       b819fb60 cufflinks version change
       a62b5195 Updated dependencies list (gqMicroarrays)
       0a5aeba1 Added sailfish module install script
       dfb486d3 bowtie 2.2.2
       97999d49 minor changes to top hat and vienna install script
       00d65aa6 Newer org.MeSH.* packages are too numerous and their installation take forever. Excluding them from the org.* packages list.
       7d51b6a9 Created gmap-gsnap install script (dev)
       af8b6c08 Added HTSFilter to del list
       21cc2016 more package dependencies
       f34a34dc Fixes to Trinotate related install scripts
       b937e6f7 Removed G phase1 from R deploy script
       8ae1b01d Multiple install scripts related to trinotate
       d3a8fe24 module install script for BEERS, RNA-seq Illlumina reads simulator
       7aca566e no message
       ba527a2c GCC module call for phase1 only: before compilation and called by R module for run time.
       915b729e modified example R.sh calls; need sh -l for module call to work on guillimin phase 1 ....
       2402dee2 bash synthax fix
       34d5126e entrezdirect installation script
       2de4f304 Added -p option to R.sh: name of an env variable that specifies a prefix to -i and -m.
       ba1638a4 no message
       9903c2a2 cron modified, next clip module
       c071c1e4 no message
       6c8c26f4 Nextclip install script
       99a90844 wrong file name in cron script
       c4ed552f no message
       d21308ef temporary commits.. sourcing from Dropbox because problems with abacus
       42a0781c Cleaning and and chomping module files too
       21ef6d8f Changed chmod commands
       0e7cacb0 Problem with roxygen2 on CRAN. install.packages does not find it for R<3.0.2. Msg posted to r-help. R.sh will not roxygenize until this is fixed.
       67164c5d Last commit to R install scripts. Three files added: - R.sh is the workhorse which checks installs R/Bioc if necessary and performs updates. List of dependencies is hard-coded within this script. By default it will install/update the latest R version to $MUGQIC_INSTALL_HOME_DEV, hacking build.R to insure any subsequent package installation is performed with umask 0002. It can also install specifyc versions with -v, install to a specific location with -i and create the module file at a specific location with -m. The -f option forces rm -rf on any previous installation. The latest R version number is obtained by parsing the VERISON file from a R-latest.tar.gz download.
       3c6a6c60 tar is silent, chmod done outside R
       baacdadf Combined installation/update in one script. package list and old install script will be removed in later commits
       d186f78b Added R/Bioc Update script
       453fdbd0 added varscan install script
       8a2f5432 Updated R install script
       a3a16f38 list of R packages deps now linking to mugqic_ressources
       22b02428 no message

  François Lefebvre <lefebvrf@gmail.com>      3 commits

       4ae7affb tophat and bowtie2 according to template install script
       35e4743c Minor changes to deploy script
       cc01370a no message

  Joël Fillon <joel.fillon@mcgill.ca>      363 commits

       bd175def Updated README-release.txt with more info
       21cbf655 Version bump to 2.0.0
       c3aea54f Added utils/ to PATH in module mugqic_pipelines
       ae6ae813 Updated module mugqic_R_packages to version 1.0.0
       b91ac7bc Standardized mugqic_pipelines module install script
       493da6ef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       6eb4df27 Updated PacBio pipeline steps as docstrings
       2ebef9da Updated module mugqic_tools to version 2.0.1
       718961a5 Fixed bug metrics.matrix sort tmpMatrix.txt before join
       ab9e0017 Fixed bug Picard sort_sam: add BAM index as output file only if sort_order='coordinate'
       66108144 Added trimmomatic adapters FASTA files in bfx
       975beb6e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       01f101c0 Minor update in rnaseq de novo pipeline
       4b11101d README.md edited online with Bitbucket
       d5b5827f README.md edited online with Bitbucket
       e897da54 Updated SnpEff to 3.6 and snpeff_genome=GRCh37.75 by default in dnaseq pipeline
       078f739f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       d0e8ec03 Separate bed2interval_list job from picard.calculate_hs_metrics in dnaseq pipeline
       347f8884 Updated mugqic/tools/... to mugqic/mugqic_tools/... in config files
       5b93fc37 Updated and standardized module install mugqic_tools-2.0.0
       b7626d27 README.md edited online with Bitbucket
       fd69eeb4 Renamed mugqic_pipeline into mugqic_pipelines
       36aed547 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       66a477b5 Added multiple input files for star step in rnaseq
       35ac8b05 Added star memory/cpu settings in rnaseq.batch.ini
       9990a91b Minor output changes in star.py
       dc44a1dc Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f7667c98 Added step multiple candidate input file support + updated trimmomatic and dnaseq bwa_mem_sort_sam, picard_merge_sam_files accordingly
       a6dda261 Added BAM index file as output of Picard merge_sam_files and sort_sam
       8c3c87f2 Removed unused trimmomatic skip option
       bc4bf739 Adjusted picard_merge_sam_files ppn value in dnaseq.guillimin.ini
       11c164cd Added symlink to BAM index (*.bai) in nanuq2mugqic_pipelines.py
       ae0a63f3 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6cca375a Fixed BED Files split bug in nanuq2mugqic_pipelines.py
       5edef090 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       98f2dd85 Updated dnaseq dbNSFP2.0.txt to dbNSFP2.4.txt.gz
       a5926c5d Replaced dnaseq dbsnp and known_sites parameters by known_variants + updated genome config .ini files with dbsnp_version + minor fixes
       11f62458 Added vcf.gz.tbi in genome install
       4cf7e435 Added log_report.pl and nanuq2mugqic_pipelines.py
       86167d7f Reverted to previous pacbio_assembly name
       fd9fd0f5 README.md edited online with Bitbucket
       3e0de926 Added RNA-Seq report step (needs update for cuffdiff section)
       165073d5 Removed differential_expression.goseq for cuffdiff and added self.gq_seq_utils_exploratory_analysis_rnaseq jobs in RNA-Seq pipeline
       3ea16d35 Merged in jfillon/readmemd-edited-online-with-bitbucket-1417115922753 (pull request #6)
       7dbb1883 README.md edited online with Bitbucket
       8c70917d Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a3f02c37 Recovered previous gq_seq_utils_exploratory_rnaseq modifs
       11863cb8 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       62964d7c Update rnaseq_denovo_assembly.guillimin.ini with generic metaq, proc, pmem cluster settings
       29a80876 Upgraded smrtanalysis module 2.3.0 with patch 1
       08552cbb Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       f604b0fc Recovered rnaseq.base.ini with new differential_expression and goseq settings
       6fb27160 Updated all R modules with mugqic/R_Bioconductor/3.1.2_3.0 and mugqic/mugqic_R_packages/0.1 if necessary
       4bca9119 Fixed merge conflict
       8aab691b Added exploratory rnaseq step (partial)
       fe1bcae2 Added genes file in rnaseq config
       11016523 Updated Star module to 2.4.0f1
       f64a60c3 README.md edited online with Bitbucket
       9c0818aa README.md edited online with Bitbucket
       c0fc855f README.md edited online with Bitbucket
       45c78c09 Merged in jfillon/readmemd-edited-online-with-bitbucket-1416841132369 (pull request #5)
       a7dfa972 README.md edited online with Bitbucket
       6583d2f0 Completed implementation of cleaning feature
       cb4c0c6a Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       54378d51 More on pipeline cleaning feature
       8113c517 Fixed merge conflicts + first implementation of pipeline cleaning feature
       8c3c120e First implementation of exploratory step in rnaseq pipeline
       4ccebd1c Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       9c33f31a Added Gene Ontology files in genomes + added goseq in rnaseq pipeline
       5b2a7a97 First attempt to install genome Gene Ontologies annotations
       1a34018e Fixed bug forgot super() call in Illumina and PacBioAssembly pipelines
       a46d7728 Updated all config base.ini with module_mugqic_tools=mugqic/tools/1.10.5
       97243c7a Tagged mugqic_tools 1.10.5
       10f05ae0 README.md edited online with Bitbucket
       d624f804 README.md edited online with Bitbucket
       b52b592f Added design file description + minor changes in README.md
       c79b9c93 Updated rnaseq_denovo_assembly differential expression cluster settings
       e8d968cf Updated rnammer cluster resources
       b882786a Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e703a738 Added Python interpreter version check
       9fb04728 README.md edited online with Bitbucket
       fe5d9541 Updated ini files with RAP_ID comment and rnammer settings
       6fd646d2 Removed openmpi settings on abacus
       1486f615 Design parsing support for empty or '0' fields
       171bf2c5 Removed command dump in .done file due to bug
       f4226eb3 Renamed readset and design file 'SampleID' column into 'Sample'
       74bef6d0 Updated cluster queue settings for rnammer_transcriptome in rnaseq_denovo_transcriptome
       a8ba35e0 Fixed version _argparser bug
       f9796a56 Dumped command in .done file + moved version in MUGQICPipeline
       a3d02bb7 Added differential_expression in RnaSeq pipeline
       38731979 Fixed syntax errors in rnaseq
       b0494b97 Added readset path normalization and variable expansion in readset file parsing
       531b02b2 Fixed PacBio readset parsing bug
       f17853b3 Added guillimin config file for rnaseq_denovo_assembly pipeline
       bf1e01f9 Moved genome_configs into resources/genomes/config; updated README.md accordingly.
       7eb3adab Imported mugqic_resources repository as 'resources' subtree
       be12079e In PacBio summarizePolishing, removed redundant cmph5 sort and used a copy of cmph5 file to prevent Job from being always out of date
       a3cf8443 Updated PacBio pipeline with smrtanalysis 2.3.0
       76f9f332 Improved Job debug log
       d658c842 Minor cosmetics changes in README.md and PacBio readset file
       7ce1851c Added Bitbucket URL in command line help
       aab8a99e README.md edited online with Bitbucket
       0c66e6f0 Added pipelines/pacbio_assembly/README.md; removed pipelines/pacbio_assembly/README.txt
       49b0f0c5 Added
       48eb8360 README.md edited online with Bitbucket
       e0ad9310 Added pipelines/rnaseq/README.md
       329a627f README.md edited online with Bitbucket
       744fe945 Added pipelines/dnaseq/README.md
       e76ca730 README.md edited online with Bitbucket
       3426e21b README.md edited online with Bitbucket
       0e4226f0 README.md edited online with Bitbucket
       24a3bf37 README.md edited online with Bitbucket
       415922ff README.md edited online with Bitbucket
       8bd144ec README.md edited online with Bitbucket
       34db02a2 README.md edited online with Bitbucket
       44024ca6 README.md edited online with Bitbucket
       a8a02279 README.md edited online with Bitbucket
       0d57ec02 README.md edited online with Bitbucket
       f6bb2f47 README.md edited online with Bitbucket
       5442dec9 README.md edited online with Bitbucket
       d738c55e README.md edited online with Bitbucket
       73a1d8c1 README.md edited online with Bitbucket
       47ae71cf More on README.md
       0ae6c2d2 README.md edited online with Bitbucket
       9a9751df More general documentation on pipelines
       cce88332 Replaced internal RAP ID with generic  variable in dnaseq.[guillimin|mammouth].ini
       49f5e04e Updated igv_genome config param to match generic <genome>.fa.fai index + updated all genome configs with Ensembl 77
       bd09fa8a Module versions update in config files
       ace8e6f6 Fixed pacbio filtering dependency bug
       c971597d Minor comment update in homer module install
       0ba8de74 Changed DnaSeq snpeff_genome to hg19
       8abed9c1 Fixed bug readset.fastq[12] attributes should be writable
       dd782873 Reverted SMRT Analysis module install script to version 2.3.0
       9a34f430 Added SMRT Analysis 2.2.0.133377-patch-3 in module install script
       68c77e55 Update SMRT Analysis module instal with 2.3.0 patch 0
       b06999f6 Fixed dbSNP download_path bug
       4612110b Updated dbSNP to build 142 for Homo sapiens genomes GRCh37 and GRCh38
       b8e26ad0 In genomes/install_genome.sh, added functions skipping if up to date + added rRNA FASTA creation and BWA indexing if present in ncRNA fasta
       49c8cf30 Moved star_index/ into genome/
       e0fe746f Added install_all_genomes.sh ; fix STAR module version for index creation
       7e92a39e Added STAR index
       d67a0c4d Added version number + updated README.md (incomplete)
       515ab99a Added star module install script
       39d4b84d Updated gqSeqUtils report calls and parameters
       162ec0bc Fixed merged conflicts
       93a94e93 Updated RNAseq nozzle reports with Python version
       e4bae26b Updated DNAseq nozzle reports with Python version
       5ad7b164 Renamed bio module into bfx
       ad6c8b17 Added config file list support in gq_seq_utils report  and pacbio_assembly
       83cd4911 Added report and compile steps in pacBioAssembly pipeline
       59fafec4 Added log.debug for job.is_up2date + set config.filepath with config.trace.ini file
       c90118a7 Added raise NotImplementedError for PUUre and RRNATagger pipelines
       47022153 Updated config genome paths with new genome files organization
       13471eb5 Fixed bug smrtanalysis.reference_uploader transforms '-' into '_' in fasta filename
       1d948bc1 Moved /sb/programs/analyste/genomes/blast_db/README.txt on bitbucket genomes/blast.sh + move oryCun2.sh in genomes/old/
       46062c58 Updated default dnaseq genome paths + readset library, run, lane and config sequencing_center are optional for dnaseq bwa mem
       45ec2077 Fixed permissions bug for Homo sapiens genomes install
       eb7ea9ce Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d052f61c RELEASE -> VERSION; DBSNP from NCBI for Homo sapiens; symlink DBSNP to Ensembl VCF for other species
       f103474d Added pacbio mummer step
       f2e84d74 Added pacbio steps up to blast
       29af2de9 Added job input/output files debug in pipeline
       3ae18d9b Fix cluster settings from job name instead of step name
       1da5b288 Fixed pacbio pbutgcns cmd missing &&
       412287e9 Check job input/output files as valid paths instead of regular files
       a0727600 Update genome installs with Ensembl release 77
       a38bdfd7 Fixed rnaseq cluster_submit_cmd_suffix bug
       7e0c6e84 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       799ccb43 Fixed permissions on genome files
       7b1b697d Use Ensembl genome instead of 1000genome
       c7dd3496 RNASeqQC now uses GTF with transcript_id only
       7aa6712c Removed abacus config .ini
       b6689cc1 Updated ini files with proper genome keys
       2e9c839b Added genome configs
       802c39db More on pacBio assembly steps
       49bdc5c8 Updated rnaseq genome config names
       e91a300f Increased abacus cores and memory for insilico_read_normalization_readsets
       59e7accb Added DBD::SQLite in perl install dependencies
       3fae17c4 Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       ca5b1259 Added trimmomatic ram in config + added picard.sam_to_fastq VALIDATION_STRINGENCY=LENIENT
       f551b415 Fixed trimmomatic headcrop length order + insilico_read_normalization_all job name
       eeb1ff62 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       37a12606 More pacbio steps + here-document for bash command to avoid escape characters in output
       4cdd7a9e Added DBI and PDF::API2 as CPAN dependencies in Perl install
       0c4d6cd7 Fix RSEM reference dependencies
       59f01f96 Standardized perl module install with reduced CPAN dependencies
       f51ae0d2 Minor trimmomatic thread update in rnaseq_denovo_assembly.base.ini
       c3be9ab3 Added rnaseq_denovo_assembly.mammouth.ini + minor path fix
       74e88c07 module_tools -> module_mugqic_tools + pacbio filtering step (incomplete)
       1870e1ab Added rat,dog,chicken genomes + python fix compiled in ucs4 + matplotlib manual install
       35997160 Standardized chipSeq modules as separate weblogo, homer, macs install scripts
       07f6d3bd Tagged mugqic_tools 1.10.2
       f29ef001 Various fixes in rnaseq + pacbio first draft
       ea16a1d2 Removed picard_index and sam_index subdirectories + error tweak for gunzip human_g1k_v37.fasta.gz
       2cb5f0ba Added dnaseq alignment readset subdirectory + mugqic_log not sent if pipeline has no jobs
       94590ba4 Cleaned genome_configs
       17a3b31f Cleaned genome_configs
       ec62c581 Minor documentation change
       a62ab926 Removed old RRNATagger files
       3d527803 Added MUGQIC parent pipeline with call_home log function
       a409b72a Added config trace creation from merged input config files
       178d1c12 Change readset file column name Sample -> SampleID; remove -V from qsub; updated gtf tophat index path in rnaseq
       0b1bd44f Creation of the major model genomes with new standard organisation
       8a6a093f Fixed puure.abacus.ini filename
       dba20ba5 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       cf3c3bac Fixed merge_and_call_gvcf output file name + updated dnaseq with mugqic/GenomeAnalysisTK/3.2-2
       8614cd1f Updated rnaseq config files and parameters + minor dnaseq config update
       ab1c6d14 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a1466e78 GENOME_INSTALL_TEMPLATE with Ensembl support
       0b2880b7 Fixed mkdir metrics in illumina.py
       66b529d2 Minor fix in dnaseq.batch.ini
       e0023513 Minor fix trinotate nb of columns + scheduler date formatting
       bff03019 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       bc2b9fc2 RNASeq De Novo pipeline implemented except nozzle report
       23a4a722 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e6212716 Fixed samtools_sort and baf_plot options for dnaseq on mammouth
       4491132d Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       7b211be3 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       c91f8351 Added fastasplit in RNASeq De Novo pipeline
       dc504422 Minor fix in dnaseq.mammouth.ini tmp_dir and R version
       a6570ef6 Minor scheduler formatting
       2a9c2985 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       edaf6f10 Added normalization steps in RNA-Seq De Novo pipeline + fixed python lib bug with absolute path
       f88a4483 Modified -j and -c options descriptions + removed all .pyc files
       47a1c952 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95b82a39 Restandardized mugqic_tools install template
       2ab7d2ac Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       10f5ecce Added first draft of Puure pipeline
       82006acb Added dnaseq mammouth config file
       e4b99eea Added dnaseq guillimin config file
       6f1cf242 More on GENOME_INSTALL_TEMPLATE.sh
       5dc06903 First draft of GENOME_INSTALL_TEMPLATE.sh
       4eb98234 Minor fix in Python module install
       c950b250 Fixed python module install + minor fix in MODULE_INSTALL_TEMPLATE.sh
       ddd1e5d6 Standardized python module install
       8eede7f8 Moved AMOS module install in dev and standardized it
       e2912e03 More config standardization
       e22b5c4d Minor fix in ucsc module install
       0c9abe7f Standardized UCSC module installation + minor template improvements
       2873b94a Fixed merge conflicts in bio/rrna_amplicons.py
       475e92cb Reorganized all config parameters
       b8d0673d Minor fixes
       d6094312 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       d73fe495 Fixed minor tophat bugs + scheduler proper exit code
       e19aaae8 Added rnaseq steps up to cufflinks
       6cf0a065 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b65411f1 Added cufflinks module
       104113f2 Added gtf in htseq + argparser factorisation in pipeline classes
       8b68b091 Added bio/htseq.py
       dd389672 Fixed newline and output file bugs in R mugqicPipelineReport + batch scheduler with job.done
       96fcc8e4 Fixed merge conflict in samtools.py
       631d37f6 First htseq draft in rnaseq + snv fix in dnaseq
       655c6bb0 Added multiple ini files argument feature + first draft rnaseq raw_counts
       7f92b490 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       a351008c Added rnaseq wiggle step
       3d3d8440 Fixed input/output files bug in concat_jobs()
       bdf6c8bb Removed all .pyc from git
       f6995a9b Added rnaseqc step in rnaseq + other minor fixes
       cec67b7a Grouped common dnaseq and rnaseq features in a parent pipeline Illumina
       192486b2 Fixed minor bugs
       4441fb94 First rnaseq draft; dnaseq with os.path.join
       443c9c92 Update input/output files for metrics_snv_graph and deliverable
       f776268f Modified sys.path according to dir changes
       710c0bd6 Reorganised directories
       55e7904c All dnaseq steps implemented (need to test though)
       7cde825a Added dnaseq steps: rawmpileup, rawmpileup_cat, snp_and_indel_bcf, merge_filter_bcf, filter_nstretches
       af2c35a0 Added dnaseq steps merge_and_call_gvcf, dna_sample_metrics
       b0ea4b8c Added dnaseq haplotype_caller step
       952ebecb Formatted batch scheduler + calculate_hs_metrics, callable_loci, extract_common_snp_freq, baf_plot steps
       2a3b223e More dnaseq steps + job input files validation
       dba57675 Moved python scripts in mugqic_pipeline directory
       c941e47b Removed step loop + standardized step names
       96e6e96b Added indel_realigner step + fixed job is_up2date with dependencies
       6cedd334 Added first GATK functions + extended dnaseq pipeline
       2c2937cd Added first gatk draft + parse Nanuq readset file
       397fb364 Added global config object + picard module
       4dd79a1e Added group, pipe, concat job functions
       fef73cd3 Added DnaSeq trim and mem + force option
       24057c2a Torque scheduler test OK
       b484d113 Added torque submit header
       dc45f445 Added scheduler and log
       b6fada28 Updated mugqic_tools install with version 1.10
       16e9aa35 Added argument parsing and more
       2e6cc0eb First Python version with basic args parsing
       3ee3e3ed Core classes coded
       97d89408 More on python prototype
       fa8dd9ad Added DB_File Perl module + minor permissions fix
       a24fd3a9 Added prerequisite for modules rnammer and trinotate
       e6d3ded6 Standardized trinotate and dependencies module install
       af7cea5d Standardized BWA module install
       53b96059 Moved R&D modules in dev/
       fbccba11 Added archive check in MODULE_INSTALL_TEMPLATE.sh to avoid redundant archive download
       1e72c28d Minor fix on template install: archive dir with /
       8898ddf8 Restored blat module in dev
       54192b2f Removed blat module (already part of ucsc module)
       6a70836b Prefix java module by openjdk
       d5575c0d Updated Java module install with OpenJDK
       3e9546b7 Updated RSEM and Trinity module install with latest version
       f80d0d29 First Python draft
       65176a94 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       28f5ab6f Fixed conflicts in merge master into bam2fastq branch
       26c226ac Resolved conflicts with master branch
       1729a4a3 More on object-oriented dnaSeq pipeline
       2378320c More on object-oriented design
       a3dabbed More on object-oriented redesign
       3d942d88 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       455ebede More on object-oriented pipelines
       20bc7af9 Minor variable name changes
       b2c5d0f6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       d017c57d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq
       9b1ae904 Incomplete version of bam2fastq dnaseq pipeline
       e87e5b92 Added Parse::Range in Perl modules
       29bd3d93 Solved merging conflicts in rnaSeqDeNovo.mammouth.ini
       a91dbdc4 First object-oriented version of RNA-Seq De Novo Pipeline
       2255e9c9 Minor RSEM install script fix
       a7c12759 Updated RSEM install script by modifyin Perl script shebangs with /usr/bin/env/perl
       a37edd04 Updated mugqic_tools to version 1.8
       063caa27 Added PerlIO::gzip as module dependency
       2fa1750e Added dev/prod/ environment in MODULE_INSTALL_TEMPLATE.sh + modified all shebangs with bash
       ae102899 Merged dnaSeq.pl conflict from master to bamToFastq
       a6d0d761 More on bamToFastq process (uncomplete)
       f911e8a1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       89a30662 Merged master in local branch and resolved conflict
       ced737e9 First bamToFastq attempt in dnaSeq.pl
       338afbec Minor fix in snpEff archive name
       8154c76e Standardized snpEff, VCFtools install scripts; added Java install script; added various permissions
       17a4f953 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       c351f099 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       b64d2e76 Merged with master
       27c29de8 Even more on bam2fastq
       9c38e401 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       04f66e44 More on bam2fastq
       ca1172ce Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       e37f19ae Some steps toward bam2fastq support
       ab015b83 Standardized tophat and cufflinks module install script
       201cc7d7 Standardized bowtie, exonerate, samtools module install scripts
       bbc946ef Fixed module directory named "pipeline" instead of "mugqic_pipeline"
       25141fc1 Standardized mugqic_pipeline module install script
       732d095f Added others permissions for MUMmer and wgs-assembler
       493b28b0 Updated mugqic_tools version number to 1.6
       689141b2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ea2d0aeb Installed trinity and rsem in production
       9b4f23fe Standardized Trimmomatic module install
       ff81b272 Updated RSEM version, RSEM and Trinity permissions
       e51836d6 Upgraded to mugqic_tools-1.5
       f205a7aa Added BLASTDB variable in blast module install script
       b3f1686d Minor permission fix
       9147435c Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e894200a Updated Trinity module install with last version
       819415a0 Updated mugqic_tools install script according to MODULE_INSTALL_TEMPLATE.sh
       5093baf0 Renamed module install template
       9473a904 Added read/execute permissions for others in module install template
       0e6cf37b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       99c14da5 Added ape + ctc R dependencies for Trinity differential expression analysis
       f6648d9e Updated Trinity edgeR PATH
       1f140b6e Removed comment in Trinity install script
       0348f395 Updated trinity and rsem module install scripts in dev
       5e0865f7 Moved some install script to dev/
       91784d03 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       29666f4c Clean up obsolete/dev install scripts
       4ec387a7 Removed Makefile draft
       90d94d3a Minor change
       77111238 Initial import of genome and module install scripts from mugqic_pipeline repository to mugqic_resources repository

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       9f787f31 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       95ac9e38 BFXDEV-133 detected perl bug on homer and weblogo scripts using our own perl module, updated perl scripts shebangs
       4b5a245a Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       72a19c1f BFXDEV-133 Adapt the software/module installation to guillimin phase2 requirements
       90db5b54 BFXDEV-82 version 1.7 of mugqic tools, added bug correction for R-tools related to chipSEQ pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      12 commits

       0418260b Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       d33f8e24 prepend popoolation path to PERL5LIBS
       c5a79a8c Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       05d65322 added Picard Mark duplicates no optical v 1.90 to the dev software repository
       5edf0ed9 pmarquis: Added genome setup for arabidopsis thaliana-TAIR10
       4bc9c17f Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       a8a9a0f7 added env-variable POPOOLATION_HOME to the module file
       13c9b20e added an installation script + module for popoolation
       f2046f10 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       2c85b4f1 added plotrix to the list of R packages
       eaab7c17 creating R v.3.0.2 module in dev
       3a6476df installing R on the test environment

  jtrembla <jtremblay514@gmail.com>      24 commits

       4d370a66 rrnatagger.py continuation
       a228ba6e completed bio/rrna_amplicons.py
       a4bbfc71 continued coding of rrna_amplicons.py.
       f82ac522 writing rrnatagger pipeline in python. First step works.
       cf405143 Added two functions (that actually works :-)). I need to convert the rest to python syntax.
       bc82663d added rrna_amplicons
       d045a6a4 renamed RRNATagger to rrnatagger
       d0a7f341 Added rrnatagger.py and list of steps.
       9f616eed added lib path.
       b0d9af62 Added bamtools.
       907ef6be Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       87e713ab Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       e3ddbbb1 minor fixes.
       058fa827 Updated bwa version.
       6432b27d New module install scripts.
       cecd97f0 Added path (prepend-env) to root/R-tools.
       c61c904f Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       119c7add minor modif to module list for perl.
       2fe26541 Added comments to smrtanalysis installation (2.2.0). Tested and installed.
       9858f66a Added muscle aligner.
       9403ab75 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources corrected typo, then merges..
       ac2a43e4 corrected typo.
       5c1e23a0 Added chmod at the end of the install module "scripts".
       79e2d4a6 Added list/script to install modules. Could be improved...

  Julien Tremblay <jtrembla@ip03.m>      1 commits

       72347dfc added module DB_File to perl module installation.

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       7d2d5a3a Added RRNATagger-tools to prepend path in the module file. BFXDEV-31
       b792e6af added python install script.
       1be6f845 Fixed two missing paths for PERL5LIB
       0b26b2d1 Added perl to mugqic modules.
       82b9d197 Modified install scripts according to our group's template.
       3156d84b Added gnuplot-4.6.4 and memtime-1.3

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      1 commits

       923d6619 Removed genelenght.r. Use mugqic_tools/R-tools or gqSeqUtils directly with HERE docs or Rscript -e

  lefebvrf <lefebvrf@gmail.com>      14 commits

       0be2b767 Removed superfluous (already in base) dependencies from list
       4cf1272b violin plots package added to list of dependencees
       16882a5e PerlIO::gzip is a dependency for Trinity clusterProfiler as R package
       0530c576 Guillimin phase 2 added to daily R deploy
       9a37c3ce R.sh needed a module call to a more recent gcc than system one.
       e90b9cad devtools::install_local() would install packages to home folder if exists...   .libPaths(.Library) solves the problem
       923b50ad —vanilla removed from R.sh too
       a13f5820 blurb added related to last commit
       d028e16a Add creation of Rprofile.site to R install script. This will force using cairo X11 backend since cairo is not always set to default when available…
       321311bc Added ViennaRNA and mirdeep2 install scripts
       33dfb878 Small hack to tools::build.R when installing R allows umask 002 (!)
       74967799 Corrected exit code, first version of R_deploy
       b8db61a3 Additional setdif(reps) to avoid duplicate installations -> slow
       2082c05c Fixed leading tabs problem in module file by using <<-EOF here tag. Neutralized R_LIBS to insure installation proceeds correctly when R module already loaded Added vanille biocLite() call

  lletourn <louis.letourneau@mail.mcgill.ca>      45 commits

       e76f595b Merged old perl changed into python
       0f7f0f87 Fixed haplotypeCaller output file name extension
       4098ccc3 Version bump mugqic pipeline to 1.4
       b9625abd Version bump to 1.5-beta
       7e8d28c9 Bumped bowtie to 1.1.1
       2528f927 BVATools version bump
       277abf2a FastQC version 0.11.2
       c1d9878c Version bumped Ray to 2.3.1, removed unneccessary fix now
       45b85da1 PIcard Version bump to 1.118
       9db0c548 Added pysam
       fc0ea4c9 BWA version bump AND fixed the script, it was a bad merged script
       cc5aa56e Updated GATK and mutect
       2cd15d63 Merged diffs
       b10f91ed Added gsalib and added configure switch
       cb1777f7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       36779161 Version bump pipeline to 1.3
       58fff5f1 Force output file name
       ef35a74e Version bump bvatools to 1.3
       a441cda3 Version bump of mugqic_tools and snpeff
       9b6fa74b Added ascat
       a6694393 BWA Version bump to 0.7.8
       61097390 Added the jellyfish tool
       61089224 Updated gatk
       09e52708 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       44e530a1 Added Mutect home
       9a96cbfd Version bump
       6827d55f Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       adda7653 Version bumped wgs-assembler, bwa, picard
       674305d4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       9e6a4c64 Version bump bowtie2 bwa picard
       287aa234 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c2aac9ab Added aspera instructions
       701df8ce Version bumped vcftools to 0.1.11
       6f5a6ed5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       c8f949ac Version bumped blast to 2.2.29
       d0dab7f7 Bumped version of picard
       4e6fcc48 Updated snpeff
       ff4f3974 Version Bump of the pipeline
       34c15b08 Added variable to access the repo's location
       fa87e197 new Ray version
       a6402ac6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       ac6513f5 Updated BVATools to 1.1
       ea6aff9a Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       dca240f8 Added BVATools
       71bfcf24 Added wgs assemble

  Louis Letourneau <louis.letourneau@mail.mcgill.ca>      1 commits

       e1c5341f BFXDEV-246 Version bump

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       14060c19 Remove conflicts in modules/mugqic_tools.sh and modules/dev/star_dev.sh
       cb81397f update modules/mugqic_tools.sh to 1.10.4
       b2db257f remove decrepated python module script and add a new one
       16d30684 replace dev module in module/dev/ cufflinks_dev.sh  star_dev.sh
       fd6119e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       85a45233 add genomes/oryCun2.sh  modules/cufflinks_dev.sh  modules/star_dev.sh
       b35f6dec up-date mugqic_tools.sh version to 1.7

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      47 commits

       0eff10b2 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       6272cb0f Python - RNAseq : allows more flexibity in fdr and p-value for goseq AND remove native goseq approach - BFXDEV-289
       e78717a6 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8694d0c5 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       218dbdca Python - RNAseq: update rnaseq.py
       a1333101 Python RNAseq - updates
       6a061db5 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       5e8df544 Python RNAseq - updates
       61a12cf0 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       b1f9bc8c change RNAseq python update
       5c9e82f8 change RNAseq python cuffnorm output to include sample names
       015fab6f remove RNAseq python step bug exploratory v2
       582e5e3c remove RNAseq python step bug exploratory
       b14fbf9f remove RNAseq python code conflict
       62a0e4f0 PYTHON -RNASEQ: star, cufflinks, htseq updates
       40970ffb Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       e600a91d Python: RNA - update star + DGE anbd goseq input oputpuyt file correction
       01fccf96 resolved pipelines/rnaseq/rnaseq.py confilcts
       87f92a36 major RNA implementation: STAR, picard, cufflinks, etc...
       a39a35e8 Python RNAseq: add cufflinks > 2.2.0 new workflow - cullinks - cuffmerge step done
       d902b4d5 Python RNAseq : update STAR - add optional read sorting during alignment
       8f572bc7 Merge branch 'python' of bitbucket.org:mugqic/mugqic_pipeline into python
       8305fe70 Python - RNAseq: add STAR alignment 2 pass && add utils folder/function for methods generic unrelated to any software or tools
       c713736c Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       4cd30e38 update samtoiols module install to 1.0
       860eba6b python : RNAseq - change ini
       a6ea523b On goinig adding star to RNAseq
       33f80d67 add dnacopy in R.sh edited online with Bitbucket
       3bdf7cb5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       0f1910ea correct python.sh script that was not workin for the modules: BFXDEV-46
       27348a08 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       793e633d add specific RNAseq packages to R package list
       1016ccd9 replace correct permissions in mugqic_tools - BFXDEV-55
       fda28024 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       40d5ef58 ngoing python module package
       03dab098 update breakdancer module
       7f60581a update breakdancer module
       494455b5 pdate breakdancer module
       9deedd93 change in samtools modules and add breakdancer module
       02fde149 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       bfdd5444 add module igv
       60e67320 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       6615a3b5 create module install sh for the mugqic pipeline
       30983a97 update nugqic_tools to the new version 1.2
       834b58e8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       19f75087 update mugqic_tools
       5e75b3d8 update modules/mugqic_tools to point to the new repository mugqic_tools with tag v1.0

  mmichaud <marc.michaud@mail.mcgill.ca>      54 commits

       4617cd8e Run processing: Remove duplicate '.' in dup bam file
       f72f529a Run processing: Run blast on bam file, fallback on fastqs
       ae2fcbd1 Run processing: Add the suffix '.sorted' in the readset bam name
       7bb7a506 Run processing: Merge md5 jobs (fastq+bam) and run qc graphs on BAM (or fallback to fastq)
       5603e888 Run processing: Run DepthOfCoverage even if there is no BED file
       a81e7199 Run processing: Use bcl2fastq 'force' option instead of removing destination folder
       1890439c Run processing: Remove dependencies dependencies from copy job
       6f33eaae Run Processing: Use a kind of factory to manage the different aligners
       e006fa91 Run processing: Remove the dependency to samtools by using picard to generate the STAR BAM index.
       f46022cc Run processing: Move md5 step later, to optimize job allocation and to minimize condition when the md5 job is finished before the copy job is submitted
       8d069734 Run processing: Remove unused adaptor columns parsing
       fd4705a5 Run processing: More documentation (sample sheets)
       60abc0ca Run processing: Output the index metrics file in the output folder, not the run folder
       88a09073 Run processing: Use xfer queue for downloading BED files
       ba1b4d6b Run processing: Add samtools module version and increase number of core for STAR (a test job took 71G of ram)
       31861213 Run processing: Manually generate BAM index file when using STAR.
       1b42a02e Run processing: Fix usage in readme
       f4e87a75 Run processing: Add documentation
       a6e18383 Run processing: Add missing '.' in bam name
       7aad80b6 Run processing: Fix copy step exclusion
       7100d7c6 Run processing: Copy output file is now in the copy destination folder
       073b8a97 Run processing: Change back to manual copy job depedencies gathering
       ff84ee02 Run processing: Initial RNA-seq support with STAR
       b5969a30 Run processing: Return empty list when there are no input for the copy job
       b5b30a6d Run processing: Fix copy job inputs
       58e43261 Run processing: Fix end copy notification job queue
       43c153ab Run processing: Fix qc output file
       9d9e87ac Run Processing: Fix race condition in HsMetrics interval file creation
       7dc2b9b6 Run processing: Fix configuration for blast and qc
       edb32694 Run processing: Generate sample sheet in the submit job method
       f5f23ff0 Run processing: Fix config for index step, now using 'index' category
       348972bc Run processing: BAM metrics are run in markdup output
       176e96f7 Run processing: Various fixes for the first test run
       63944b97 Run processing: Various fixes for the first test run
       1043bb1c Run processing: Various fixes for the first test run
       c8ee9b0d Run processing: Add basic support for different aligners
       aa16b2dd Run processing: Get copy step dependencies by introspection
       c74a06b5 Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       3af310b2 Run processing: Delete existing 'Unaligned' folder when running fastqs
       583c80f4 Put Illumina configure bcl2fastq in a job
       aa93387f Run processing: Use . as job name separator
       ddafcfe8 Use $MUGQIC_INSTALL_HOME variable in config file to replace hardcoded paths
       79d09cce Fix arguments for IlluminaRunProcessing by removing readset argument from MUGQICPipeline
       cbe8ddc2 Run processing: Use run_dir instead of run_directory to follow output_dir convention.
       aa213fe5 Run processing: Improvement to the sample sheet generator to loop only one time
       7d16bcf2 Run processing: Seperate config files for MiSeq and HiSeq. The base config file still hold for a HiSeq
       df783bc2 Merge branch 'python' of https://bitbucket.org/mugqic/mugqic_pipeline into python
       2b53d110 Run processing: Add wget download of sample sheets and bed files. Fix Casava sample sheet generation.
       07aacb68 Add some imports from __future__ to ease the transition to python 3
       7af4b28f Add copy job dependencies
       82841bf4 BFXDEV-283 Fix reference for coverage calculation, add copy step.
       213e7160 Change illumina_run_processing name to follow code conventions
       8d8e1af2 Code convention changes
       f51f5e00 First draft of the illumina run processing pipeline in python.

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.ferrier.genome.mcgill.ca>      1 commits

       acde2371 adding sh tools for riboRNA

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      1 commits

       f6ddfbeb BFXDEV-112 install genomes from IGENOMES

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.ferrier.genome.mcgill.ca>      1 commits

       2fc3f9d2 usearch.sh

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      2 commits

       481bc98a Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       7ff159c4 gene_length

1.4        Mon Nov 17 13:15:48 2014 -0500        139 commits

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       667c6f7d README.md edited online with Bitbucket
       ee9efd52 README.md edited online with Bitbucket

  Francois Lefebvre <lefebvr3@ip03.m>      6 commits

       3355905e Previous minor bug fix (dots in sample names) actually introduced major bug.
       0bc3e58d temp files remove
       93c2905a rnaseqc other options possible
       c0455416 rnaseq .ini file more tophat options. BFXDEV-215
       8270318a Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.
       0ef6afc9 Added --transcriptome-index support to tophatbowtie.pm, as well as possibility for other options.

  Francois Lefebvre <lefebvrf@gmail.com>      7 commits

       242f2bc9 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c004e516 Passing projectName to SubmitTocluster can create invalid job names. Replaced with string
       e578fb41 variable name was inappropriate
       eab03100 Updated rnaSeq mammouth template .ini file.
       202215e2 Defined a job name prefix for wig zip.  'metrics'  as a job name was not enough information
       59f15282 cuffRescolumns and dgeRescolumns now in goseq param section. Also adjusted those values for UCSC genomes hg19 and mm10  templates, original value did not work.
       524dd82c overwrite.sheets=TRUE to avoid common problem of updated projects

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.ferrier.genome.mcgill.ca>      4 commits

       257c44d0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4b8a23be patch3 applied; see BFXDEV-260
       3c6e8502 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ff5fe3de change R module due to crashing nozzle report generation, see BFXDEV-255

  Joël Fillon <joel.fillon@mcgill.ca>      34 commits

       0860442e Updated pacBio .ini config file with new smrtanalysis module name: 2.2.0.133377-patch-3
       c2d50007 Updated config genome paths according to new genome organization
       030e601c Updated chipSeq mammouth config with new Homer 4.7
       be23877a Added explicit Python module in rnaseq cuffcompare step
       9d16e5a0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c2569139 Added newline after mugqicLog command
       8f3c5f35 Updated rnaSeq.mammouth_hg19.ini with generic /nfs3_ib/bourque-mp2.nfs/tank/nfs/bourque/nobackup/share/mugqic_prod path
       1fd8a636 Removed explicit RAP ID in RNASeq De Novo config files
       0efe22d7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6aac202f Removed moduleVersion.htseq ini ini config files
       ded1a460 Fixed missing java module in igvtools
       2c37b3ed Minor fix in dnaSeq step_range help
       cf086119 Added default RAP ID in RNA-Seq De Novo guillimin config file
       e67468af Added --exclude <fastq> option in illuminaRunProcessing rsync for samples having BAM files
       98420e40 Uncommented phone home in dnaseq
       0e6e8a8d BFXDEV-203 Create one .done with checksum instead of one per output file + update config files with default adapters path + update perlpods removing -e option
       98860615 Fixed bug missing Library Source -> column not mandatory + updated pipelines/rnaseq/rnaSeq.guillimin.ini with accurate module versions
       e07ffd80 README.md edited online with Bitbucket
       d02ae67f Added comment to update Resource Allocation Project ID
       13bf5c97 BFXDEV-221 Migrated abacus/guillimin config files from msub to qsub
       d8a411e4 MUGQIC call home is now run inside bash script after job submissions instead of bash creation.
       5c2627ae Updated chipseq mammouth config with python 2.7.6
       29163dd5 README.md edited online with Bitbucket
       bad36484 Added call home feature notice in README.md
       94130140 Fixed RNA-seq de novo resume jobs; added env variable WORK_DIR in pipeline; formatted bash output
       b2f3fd52 Another glob fix for prefix path check in LoadConfig::getParam
       65b99232 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       768184e8 Minor fix for config param prefix path
       74d43282 Added prefixpath check type in LoadConfig::getParam
       249a5cb1 Fixed transdecoder bug ln: creating symbolic link 'Trinity.fasta.transdecoder.pfam.dat.domtbl': File exists
       59ea4b6f BFXDEV-32 Fixed wrong transdecoder file path for missing PFAM 'cds.' prefix
       9c6fe6f3 BFXDEV-32 Fixed pfam missingcds. ID prefix + blast clusterCPU tag for guillimin and abacus
       73ce1801 Added chipSeq pipeline cleaning
       9f6e203f Fixed rnammer missing modules hmmer 2.3.2 and trinity

  jtrembla <jtremblay514@gmail.com>      14 commits

       6c0720ba --arg for low abundant clusters. after 99% clustering. BFXDEV-31
       d09db898 Added low abundance cluster cutoff argument. (After 99% ID clustering step). BFXDEV-31
       9166850c Put more lenient parameters for itags QC to make it more 'universal' for most projects, especially those for which quality of reads 2 is low. BFXDEV-31
       05e7b99c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f60fa07a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a996d3f1 Indentation correction. BFXDEV-31
       38d8da84 Fixed parsing for blastdb command. BFXDEV-30
       1b20686f Updated README for 454 data processing instructions. BFXDEV-31
       8bafb19a Added even more description. BFXDEV-31
       c003951a Added details to description of output. BFXDEV-31
       fcd21e5d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9dcc9bc Added step to filter number of blast results. BFXDEV-30
       5d0c803c Removed --vanilla. BFXDEV-30
       309ef5e2 replaced -num_target_seqs with -num_alignments. BFXDEV-30

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      1 commits

       a16e5e9e Added --sampleSheet argument to getMiSeqBarcodes.pl BFXDEV-31

  lefebvrf <lefebvrf@gmail.com>      2 commits

       b7f7d661 Changed TopHatBowtie.pm to put an end to  the .fa.fa.fasta.fa symlinks madness when installing genomes. Parameter is now the bowtie index basename, consistent with the tool's documentation.
       d29a172f Fixed dots in sample names bug BFXDEV-51

  lletourn <louis.letourneau@mail.mcgill.ca>      31 commits

       3f785de9 Version bump to 1.4
       38d0a8f4 BFXDEV-39 Fixed realigner when only one job is given
       f728ac9c Updated bvatools
       5e168e6d Tweak guillimin parameters
       ea83ecd5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1748e1a1 Updated parameters
       69b8f2ee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5da052af Adjusted cluster requirements
       4cc7adfc Removed useless param
       4dace3ec BFXDEV-256 Added step range to paired variants
       cdbd4a08 BFXDEV-256 Added step range to paired variants
       b5355703 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46d016cc BFXDEV-254 Fixed GATK 3.2 change on CatVariants
       10fb56b3 BFXDEV-252 Removed flagstat
       125945f1 Fixed BAQ from pileup and ram from fixmate
       71366af8 Changed picard to version 1.118 to fix the freeze when an exception occurs in asyncWriter
       e9112eb7 BFXDEV-216 Removed per lane metrics
       05687887 BFXDEV-245 Fixed uninitialized error when no bams are present
       682a9970 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6d2b309b Fixed CCDS location
       85578b84 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       54fde94b Merged master
       407b9fee BFXDEV-236 optimized settings and split human builds
       419e0a0e Added missing perl module
       976ba1f9 Fixed missing HS metric, add 2 cores to bwa
       9cf6dc5f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c3fab1c8 Missing validation silent
       5cfb61f2 Add genome versions of ini files
       de4a383e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       dd39c54a Fixed typo in mergeAndCallGVCF section
       33669702 Version bump to 1.4-beta

  Marc Michaud <marc.michaud@mail.mcgill.ca>      1 commits

       6b8f9f7b Add missing use

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       3c0fb6e6 update ini to fit the new mugqic_tools tag 1.10.4 - BFXDEV-275
       e5bd4766 correct dnaseq bwa samnse dependency bug - BFXDEV-253
       6fdf8026 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f3abd8cc correct cuffdiff input double array issue when checking the job object is up to date

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      13 commits

       3bd1f63e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0e991b0b Update pairedVariant
       469415ff RNAseq replace headcrop at the good position in the trimmmomatic command; remove by default headcrop from the ini files - BFXDEV-267
       534615dc Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       402db54a PairedfVariant.pl: change where pindel get the insert size info - BFXDEV-266
       87a9da68 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       192600ec PairedVariant: allow mutec to run without cosmic file - BFXDEV-263
       62e335d3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9696dbc6 update dnaseq and paired variant ini files
       37ad08b7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       21e3c79e remove conflicys in pipelines/dnaseq/pairedVariants.abacus.ini
       a156ae54 pairedVariant.pl formAt SV and CNV to new standard - part of BFXDEV-41
       cb5b91e2 cuffdiff now should not be relaunch in a restart if it exit correctly during the previous analysis BFXDEV-212

  mmichaud <marc.michaud@mail.mcgill.ca>      14 commits

       b9ef5329 Fix genome path
       ec310455 BFXDEV-269 Change genome files hierarchy
       0234cd15 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       5eb75aad Fix usage to reflect new -s argument
       f16ff095 Allow more ram for DepthOfCoverage
       40331708 Fix CPU limit usage error by using less thread for the GC
       88e50961 More RAM for QC (to avoid java heap space errors) + More threads and more walltime for mem (to avoid wall-time exceeded errors)
       d946aa77 Run processing: Add option to force download of sample sheets
       7feb9230 Fix quote escaping in filter
       88832015 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       0e60d3da BFXDEV-171 Trim index in the generated sample sheet according to the first/last index specified as parameter
       c5a35cba BFXDEV-210 Use Parse::Range for steps to run, as all other pipelines
       598ff88e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8e8ca8ed BFXDEV-211 Don't send email when jobs are successfully completed

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus1.(none)>      4 commits

       0cb571a4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f934df3e MiSeq ini
       1f34db5d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31ae4509 add ini MiSeq

  pascale.marquis@mail.mcgill.ca <pmarquis@abacus2.(none)>      1 commits

       446b4d80 ppn=16 pour mem

  Pascale Marquis <pmarquis@lg-1r17-n01.guillimin.clumeq.ca>      1 commits

       90e7a618 rm illuminaRunProcessingMiSeq_PM.ini

1.3        Mon Jun 2 10:02:07 2014 -0400        109 commits

  Joël Fillon <joel.fillon@mcgill.ca>      17 commits

       99312583 Partial MUGQIC remote log for RRNATagger
       da59f18d Remote MUGQIC Log Report in chipSeq, dnaSeq, pacBioAssembly, rnaSeq, rnaSeqDeNovoAssembly
       521765df Updated RNA-Seq De Novo config files with latest module R/3.1.0
       f59efb98 Updated rnaseq_denovo default config ini files with trinotate steps
       01945845 More Trinotate steps in RNASeq De Novo
       9a5cd2d0 RNASeq De Novo config files conflict solved + openjdk
       1cfd2880 Added file cleaning for RNA-Seq De Novo pipeline
       4382d804 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       10924783 Solved conflicts when merging master
       42e7aca1 Updated RNASeq De Novo pipeline with Trinity version 20140413p1
       240d576e Beginning cleaning
       e1c90210 Cleaning of cleaning...
       ad6ed7ab Fix on samToFastq/trimming dependencies in chipSeq pipeline
       02d3cbe7 Added BAM file check in samToFastq
       5e1aeb21 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into bam2fastq_basic
       8ea8c98f Basic bamToFastq support for all pieplines
       9afcd9e0 Version bump

  jtrembla <jtremblay514@gmail.com>      32 commits

       6e3d8edf Added ini file for production purposes (QC assembly among others). BFXDEV-30
       462ecaaf fixed end step for BB. BFXDEV-31
       06a44238 Added BigBrother modifs to RRNATagger pipelines. BFXDEV-31
       e0c59aca Added sample counting in pipeline loop. BFXDEV-31
       ebd465f6 Added a description of output files. BFXDEV-31
       e0272b5c Fixes to itags_QC. decision if primers are present or not. BFXDEV-31
       b199a3d8 mooooooore fixes. BFXDEV-31
       d67905fe more fixes to ini. BFXDEV-31
       7b8a64ce updated ini files. BFXDEV-31
       609c2043 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       befd271e Changed arg to adapt from percent to X hard cutoff.
       16248cda Replaced percent cutoff by X cov cutoff. BFXDEV-30
       d79cc388 added / updated ini files.
       65ed0f35 ini files of RRNATagger changes. BFXDEV-31 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9e30f84e updated ini files for RRNATagger. BFXDEV-31
       37c1b0c5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4cb24e16 ini file modif. BFXDEV-30
       de66a9b9 update celera specs for hgap3. BFXDEV-30
       0faaea63 Updated file check for restart at blasr ste. BFXDEV-30
       83c0cbf1 Modified way seeds.m4.filtered is handled. BFXDEV-30
       00482bd5 Updated ini file for pacbio assembly on guillimin. BFXDEV-30
       9174116c fixed lib path. BFXDEV-30
       f823b937 Fixed tab indentations. BFXDEV-30
       5f8123da Upgrade pipeline from HGAP2 to HGAP3. BFXDEV-30
       98d20c46 Removed unused SMRTpipe .xml files. Only keep filtering xml file. BFXDEV-30
       468aacb5 Updated mugqic tools module.BFXDEV-30
       a1ae2534 Fixed output .mugqic.done file for pacbio dcmegablast. BFXDEV-30
       e7426a89 Added module load perl in subroutines. BFXDEV-31
       056cae4a added perl in ini files. BFXDEV-30
       cf512b28 pacbio ini file for abacus. BFXDEV-31
       85ebc8ae Added blast parameters to ini files for guillimin. BFXDEV-31
       52a9f7fd Added blast step for pacbio rRNA tags data type. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      6 commits

       cc1fab05 updated ini file abacus.
       9363d1fa Updated ini file for abacus. BFXDEV-30
       c3a11891 replaced compareSequences with pbalign for numthreads. BFXDEV-30
       7189270e ini file for hgap3 on abacus. BFXDEV-30
       37df4a84 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1bf18fd2 Fixes to HGAP3. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      38 commits

       8f385581 Version bump to 1.3
       c051f9aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       42735dd4 BFXDEV-207 Fixed when interval lists are created
       102ee1de Fixed gvcf bugs
       613a990b Added emtpy quotes to empty keys so they don't become ARRAYs
       b594056e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d39b5d97 Updated BVATools version to 1.3
       9eae448d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       9cc18327 BFXDEV-204 removed varfilter
       82592ff0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8b4f2bdc BFXDEV-196 Added onTarget metric
       2f57d4fd Completed POD documentation
       8176ab0d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1d6793d6 BFXDEV-161 Added last steps to Haplotyper caller
       240ca527 BFXDEV-198 Added simpleChrName in the ini since it was removed from the pm
       7f0d5398 Fixed params
       39a7b383 BFXDEV-196 Added onTarget metric
       8e7a04e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       58254614 Ram was too close to max
       db33ed9d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89c52165 BFXDEV-193 Use new fixmate from bvatools. Fixed bad step dependency
       ff9d44c7 BFXDEV-192 Added possibility to have extra flags in indel realigner
       f30f24e5 BFXDEV-182 CollectMetrics sometimes needs to GC so 2 cores are needed
       20768a45 Merged changes
       de76e02a BFXDEV-181 Updated java vm version
       d9468715 BFXDEV-153 Added a way to ignore readset status
       4ad3df16 BFXDEV-176 create MD5 on alignement
       52b96be8 BFXDEV-174 Test that the data is valid
       25f2a408 Merge branch 'master' into haplotypeCaller
       d4a588dd BFXDEV-168 added depth ratio plots BFXDEV-161 added haplotype caller
       6596bf35 BFXDEV-166 fixed no index read, but one is associated
       7b4bdc2b BFXDEV-162 fixed hiseq recognition
       3ad33194 BFXDEV-157 Add callable region generation stats and BED in DNASeq
       9f0e145f Added BVATools depth of coverage as a transition phase
       17c7a8fb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c775e183 Don't use BAQ since we use recalibration, this was decided awhile ago
       1c6b7a4f BFXDEV-156 fixed csv encoding issues
       46573b70 Updated settings for phase1 + phase2 merge of hardware

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      5 commits

       35060f63 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7faafe05 Create pipeline cleaning lib; partially implemented with RNA cleaning sub only BFXDEV-74
       201124c0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       51777364 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       27c987ae change how exit status is catched and exit the correct status in case of pipe discrepency BFXDEV-140

  mmichaud <marc.michaud@mail.mcgill.ca>      11 commits

       8c059d4e BFXDEV-155 Add rat and mouse alignment
       4e9b1e16 BFXDEV-185 Separate config file (HiSeq, MiSeq)
       637526f9 BFXDEV-164 Fetch sample sheets when they aren't specified and they aren't on disk
       49e0b8c9 BFXDEV-183 Download each bed file only once
       f87afda6 BFXDEV-184 Don't rsync Images folder. Was used on miSeq for debuging purpose
       ae9782f8 BFXDEV-186 Don't use the BAM generated by markdup
       c5d36af9 BFXDEV-180 Use processingSheetId as dependency id instead of sample name which isn't unique
       6f3d923b Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       ec789768 BFXDEV-171 Add option to specify first and last index
       e3ad9673 Add option to specify first and last index
       d8eacb49 Run Processing: Add dependency on the metrics in the copy step, in case the markdup & BAM MD5 was already done in a previous pipeline, but not other metrics

1.2        Fri Mar 28 16:11:44 2014 -0400        200 commits

  Francois Lefebvre <lefebvrf@gmail.com>      2 commits

       c18dfb65 RSEM more cpus on mammouth
       9fff58f4 rnaseqdeno mammouth ini tweaks for trimming and RSEM

  gary.leveque@mail.mcgill.ca <gleveque@abacus1.(none)>      3 commits

       00721219 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8bfdcffc module load python in rnaseq.guillimin.ini --for htseq-count
       bb66dcce rnaSeq.guillimin.ini changed default tmpDir  --BFXDEV-144

  Joël Fillon <joel.fillon@mcgill.ca>      46 commits

       2ca8f737 Version bump
       0ce7d1af Solved conflict for merge master and bam2fastq branches
       8696c9bd Changed gzip to zip to compress rnaseq de novo outputs
       d86c3022 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c09c2924 Changed archive command from gzip to basic zip
       63656f17 Minor fix: removed semicolumns in chipseq default config file
       f383320a Removed unused variables and functions in all pipelines
       3924dbe6 Added Version number in RNA-Seq De Novo Usage
       3662f2e9 Minor mugqic/tools version update in rnaseq de novo ini files
       9c9b747e Minor ini param adjustment + comment fix
       be76b262 Removed hard-coded path to modules.sh: not required when invoked with Bash shell
       56881db4 Updated default project paths with new /gs/ partition
       4c24b2b6 Removed deprecated GetFastaAlias.pm
       2d92ee87 Replaced shebang #!/usr/bin/perl with #!/usr/bin/env perl in all Perl scripts
       6b7b3760 Fixed bug SequenceDictionaryParser filepath with environment variables + set param [annotateDbNSFP] dbNSFP not mandatory
       842ca5fc Removed redundant file existence test in SequenceDictionaryParser.pm
       851df470 Check all config modules only once when config hash is built, to reduce runtime
       b5f3c9f1 Minor fix for adapters path in RNASeq De Novo config files
       a86bd798 Updated adapters paths in RNASeq De Novo ,ini files
       9d73b18f Added  [Warning/Error] in msg display
       d11cb76d Removed old lib/LoadModules.pm
       1126cd98 More minor bug fixes for getParam validation
       670c892d Fixed minor getParam validation bugs
       e002a7a4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into valid_param
       872d712d Added validation in all getParam functions
       d9621c0f Fixed merge conflict in lib/Picard.pm
       88411ce3 Added moduleLoad validation + major style reformatting
       c130dc99 Added && between job commands to get right exit status
       080410a7 Added param definition validation and module availability in LoadConfig
       2ee8ea04 Added raw read dir and adapter file path validation in Trimmomatic lib
       a4615118 And more and more about parallel normalization
       c32e3162 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into sample_norm
       8b316d0c Fixed metrics parent step
       f6e131b0 Even more parallel normalization
       da10fb17 More parallel normalization
       5b5522cd Normalization parallelized by sample
       7ef8faa9 First draft of resume-jobs
       39c49a5d Use module mugqic_dev/R until properly deployed
       813d1979 Fixed typo
       3d51ae69 Fixed normalization stats file name
       b667299f First stable RNA-Seq De Novo pipeline version
       de50a894 Updated blast results and DGE matrices file names
       04e1d5bc Merged conflicted rnaSeqDeNovoAssembly.pl
       d8717e48 Added POD documentation + fixed bug blast-longest-transcript file
       cdc15b6d Added BLAST for longest transcript only, with results header
       df03b3b4 Added metrics and deliverables steps in RNA-Seq De Novo pipeline

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      1 commits

       6a2eca81 README.md edited online with Bitbucket

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       90ac5a42 BFXDEV-133 incompatibility between /usr/bin/perl and mugqic/perl/5.18.2. Added perl HOMER_HOME/bin/ to the program execution in peak annotations and motifs usign Homer
       5e62ee01 BFXDEV-133 detected bug dependencies between trimming and alignment chipseq pipeline
       972c98b8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0fe7a1c9 BFXDEV-133 updated software, parameters, corrected bugs in chipseq pipeline for guillimin phase2
       2401cf3a BFXDEV-133 adapt chipseq pipeline and configuration file to guillimin phase2

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      17 commits

       7f020642 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ba4cfb0a bug in ini files: the following variables are not defined : genomeSize, annotation distances, markdup, report variables
       942ef2a2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       050cc3cc BFXDEV-123 add design example files to chipseq and rnaseq pipelines, update the user manual
       599e532a BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       4006f00d BFXDEV-123 add design example files to chipseq and rnaseq pipelines
       8290e24f BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper - typo
       b79e9f7e BFXDEV-36 Generated PODs documentation for rnaSeq.pl wrapper
       c1873cec BFXDEV-28 added PODS documentation to the dnaseq pipeline wrapper
       5bfe14ef changed my email by johanna.sandoval@mail.mcgill.ca in standard ini files
       9d2232b2 BFXDEV-85 added flagstats/ generated a file with number of reads after filtering/mark as duplicates
       1d069aa2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       6680dcf3 BFXDEV-84 correct bug for restart when merge lanes step failed
       394b3aeb BFX-1799 wrong variable initialization for genomeSize, detected when genome is not human or mouse
       5ac34dc4 comment skip trimming step from standard ini files, added imagemagick to mammouth ini
       e2b57e1e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8c0e000b BFXDEV-77 bug in merge metrics, instructions to run Chipseq in the pipeline directory

  jtrembla <jtremblay514@gmail.com>      15 commits

       8d91c821 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline Merge latest change prior to modifs to PyNAST (unset LD_LIB).
       05f53351 Added unset LD_LIBRARY_PATH to PyNAST step. BFXDEV-31
       8b799040 Updates to pacbio stats step. BFXDEV-30
       956d5e37 Fixed pdf report for nc1. BFXDEV-31
       c578e838 updated ini files for RRNATagger. BFXDEV-31
       abb4b984 Added module load openmpi to PyNAST. BFXDEV-31
       0eeea2ea Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cd179fce Added nozzle report for RRNATagger_diversity. BFXDEV-31
       22b7fea9 Removed Iterator::Fastx from the wrapper as it is not even used here. BFXDEV-31
       884627fc fixed file checking for restart mechanism for sub loadPulses. BFXDEV-30
       d275e03f fixed file to check for input in sub referenceUploader. BFXDEV-30
       cab61530 Minor fixes to restart mechanisms. Removed inputFofn for filtering step and loadPulses step. BFXDEV-30
       e9c1d819 Corrected some parameters for celera assembly step. now on lm2 by default. BFXDEV-31
       752dd0ad Implement module load and getParam checks. BFXDEV-31
       8afc93fe added missing path for abacus ini files. BFXDEV-31

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      33 commits

       3cf8824c fixed path to primer files. BFXDEV-31
       60e8f923 fwd and rev primers options now optional. BFXDEV-31
       7f869723 Added modules for R and mugqic_tools to rarefactionPlots.R . BFXDEV-31
       add6ebd1 mugqic.done fixes to rarefaction subroutines. Added mugqic tools module to appropriate subroutines. BFXDEV-31
       72338c79 Updated help screen and removed appending ./scripts/ to PATH in curr dir. BFXDEV-31
       6a982c81 Removed scripts/ dir and moved it to mugqic_tools
       c99e920b minor fix for tmpdir definition. BFXDEV-31
       8e3af423 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f806248e minor fix to ini file (for mummer). BFXDEV-30
       dabb656a Modifications to README and added mapping file example. BFXDEV-30
       d3473a44 README for RRNATagger. BFXDEV-31
       82ad8346 Added RRNATagger (16S/18S/ITS rRNA amplicon) pipeline.  BFXDEV-31
       a9eeb683 fixes to ini files (pacbio pipeline). BFXDEV-30
       a6bc0cc2 put only 1 mersize (14) in the gullimin ini files. BFXDEV-30
       7bcfc657 Fixed bug with fofns. Minor modif to main loop. BFXDEV-30
       126618dc Changed /dev/null in the order of commands so no empty consensus.fasta anymore. BFXDEV-30
       e9d2882d Forgot a && before gunzip in variantCaller  step! . BFXDEV-30
       907be7d2 Fix for uncompressed consensus.fasta.gz. BFXDEV-30
       d5e5250e Fixes for compatibility with latest version of pacbio assembly pipeline. BFXDEV-30
       fea30cc6 Fix for .bas.h5 files. BFXDEV-30
       935aea08 Updates to PacBio assembly pipeline which now support multiple polishing rounds.
       35d5972d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       28c7b829 Compute assembly stats from polished assembly and not unpolised ones. BFXDEV-30
       ebdbae00 Fixed and updated pacbio ini file for abacus. BFXDEV-30
       920b75f4 changed callVariants for variantCaller . BFXDEV-30
       faf6c43d typo corr. BFXDEV-30
       0b000318 Added additional instructions in the README. and changes relative lib path. BFXDEV-30
       468184fa Forgot to update this library for pacBio pipeline. BFXDEV-30
       eda60c32 updated README. BFXDEV-30
       6217fda0 corrected for typo gnuplot . BFXDEV-30
       89f4c68f corrected relative path of lib folder
       fdf93053 BFXDEV-30 modified parameters so they are more generic.
       27ab9137 Loop/dependency fixes to PacBio pipeline. Added compile stats step at the end of pipeline.

  lefebvrf <lefebvrf@gmail.com>      2 commits

       5e46e902 necessary to honour cairo as X11 backend for R graphics with current module installation
       8707d782 vanilla will hinder reading Rprofile.site, which in  our modules will not be used to force cairo as X11 backend when available

  lletourn <louis.letourneau@mail.mcgill.ca>      27 commits

       8c107ea1 BFXDEV-149 Fixed the way BVATools is called for coverage.
       b839f5aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c7f0a75e Removed path test since they contain env variables
       31322fb8 BFXDEV-124 fixed center when using mem
       c2625ccb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e420186c BFXDEV-116 Fixed reverse adapter
       4f5f9d21 BFXDEV-111 Fixed many module versions
       973c7b02 BFXDEV-114 Added R loading
       4be23c33 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       419f6fec BFXDEV-109 Configured java version
       b9989e54 Updated guillimin pipeline
       70343460 Updated depth of coverage ini
       db2b3a83 Fixed case when there are no BED files
       b79a3a0f Fixed module typo
       be34921f tmp hack so nanuq takes coverage graphs
       8cf40e60 BFXDEV-89 BFXDEV-88 Change GATK to BVATools for DepthOfCoverage and support multi bed in project
       7f45e96e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b0a885e0 Added the missing mappability key
       5b568bab Added line to keep overlapping pairs
       73d736b8 Added umask to dnaSeq jobs
       7e1f584a Updated module versions
       67b9a786 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       472a3cd8 BFXDEV-73 Fixed undef jobId if step is up 2 date
       a6e3c7c9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2d9d66a8 Added details to the release guide
       70f59da5 Version bump 1.2-beta
       daabc8db Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      14 commits

       3cd0dee6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22a87a0f DNAseq: correct SNVmetrics dependency BFXDEV-136; RNAseq: remove hstseq dummy/null module usage BFXDEV-137
       cbdca419 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       35583708 DNASEQ: correct depthofCoverage missing bed file field in sample sheet BFXDEV-135
       38db1aa4 ask bash to source /etc/profile.d/modules.sh BFXDEV-125
       184fd53d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3cbe61a6 add the -T option to cuffcompare BFXDEV-93 and validate previous commit for BFXDEV-47
       a671e06e README.md edited online with Bitbucket
       f9a76160 README.md edited online with Bitbucket
       8f494a85 README.md add the bioinnfo email adress
       6cb9d24b correct the 1st line typo in chipSeq.pl
       f1df7bab Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       31f3438e replace ; by && in pipelines except for the rm of .done
       8620232b replace a default ppn=20 by ppn=12 otherwise the job will never be launched

  mmichaud <marc.michaud@mail.mcgill.ca>      33 commits

       e242aa19 BFXDEV-147 Fix index mask of a single-index lane in a dual-index run
       3c527b26 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       bce06662 Use new version of BVATools with simpleChrName support
       0cf39804 RunProcessing: gentle perlcritic compliance
       6e3b82e4 Fix BED file list from SampleSheet
       4753244b Fix again the readsqc BVATools path
       becd68c0 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       c8fc20f3 MD5 job doesn't need 5 cores, only one is ok
       7b03a8f7 Fix readsqc using current session bvatools jar file instead of the one loaded on the compute node
       230a617f Fix run processing when there is no target file
       cb1be91f Merge branch 'master' into runProcessing
       59829b2a Use new module loader for BVATools qc
       61c25ba0 Add cranR module version
       7e38002f Merge master into runProcessing branch
       1261ab2a Support spaces in bed files
       2b5b54e2 Don't align mice yet
       9bb906f7 Support spaces in bed files. Fix bvatool module version
       3e2571b0 Fix usage message (to show optional parameters) and print error message on dying
       f549388a Die when there is a barcode collision. Fix bwa rgCenter for mem
       172452bb Add a parameter for the number of mismatches
       da5e47a1 Default values for both sample sheets path
       214c3d92 Merge branch 'master' into runProcessing
       9584091d Update BWA module version
       1a6989d5 Fix rsync command: quote character were not escaped correctly
       a8e58fdc Add missing depedencies for the copy job (metrics)
       f19991cb Don't create qc folder when the qc job will not run
       7385c8e2 Depth of coverage: add reference parameter
       7508bd41 Change '_' to '.' as seperators in the metrics job name
       2ea61036 Enforce processingSheetId column in the sample sheet only when processing a run
       f79fd8ee BFXDEV-76 IlluminaRunProcessing Add target coverage metrics & change job ids to support multiple sample with the same name in the same lane
       8e90e623 Merging upstream changes of the barcode metrics jobs. Less core and memory used, skip last read
       5a2b132d Merging master into runProcessing branch
       70383b36 BFXDEV-76 Illumina Run Processing Pipeline

  Pascale Marquis <pmarquis@lg-1r17-n02.guillimin.clumeq.ca>      2 commits

       21722231 rnaSeq.guillimin_hg19.ini
       d75c72d7 rnaSeq.guillimin_mm10.ini

1.1        Mon Dec 23 14:23:34 2013 -0500        137 commits

  Joël Fillon <joel.fillon@mcgill.ca>      69 commits

       30513928 Commented code
       78e706c6 Check rawReadDir value in configuration file
       3e137a01 Minor R version fix
       2eec3e95 Added trimMetrics step, fixed blast outfmt single quote bug
       758a68b5 Added option file check, module availability check
       7fabc9bc Added genes/isoforms length values in edgeR results
       faefc933 Added blast description in edgeR results
       c54b3367 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e61eb0f3 Cleanup of rnaseq_denovo files
       02ed5e69 Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       378d4d97 Updated ocnfig .ini files with trim parameters
       2da3da4a Fixed missing escape $R_TOOLS
       470df00c Fixed differential expression bug
       36711d42 Added differentialExpression step
       467076bf Added multi-step dependency support
       eba14059 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       003c41b2 Fixed Trinity.pm merge conflicts
       6046b1b9 Added trim step + various fixes
       efc7b744 Minor fix
       5c140d0c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e31fb3ce Added guillimin config file + minor modifs
       faf0d714 Updated .ini config files with cluster-dependant processing values
       19f6b6d4 Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       dabafcfc Added blastCPUperJob config tag
       0eb0a4a7 Added blast step
       01b495ff Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       885d5f10 Updated default rnaseq tmpdir on guillimin: /sb/scratch/jfillon
       77683166 Added first version of BLAST step
       03e5423e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       0a05c4e8 Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       fcb3c5de Minor fix
       dc1e136e Merge branch 'rnaseq_denovo' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       a80a6c28 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       29f179f4 Updated config ini files
       63a2a5d2 Minor fix bowtie CPU
       b3cf0fc1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       27fd66f7 Massive cleanup
       f2eeb87e Another minor fix in mammouth ini file
       60fac504 Minor fix mammouth ini file
       2822aab1 Updated mammouth ini file
       a7f47783 Added TrinityQC step + code reorganization
       340a63d0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       9d7d75ec Moved normalization parameters in config file
       950b713c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       3fd95314 Fixed bug edgeR semicolumn
       f850faf6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       c9e35cf8 Added edgeR step
       42d400d2 Added rnaseq_denovo.mamouth.ini
       a20391f8 Fixed bug unsorted list of fastq.gz from find
       4a35ad65 Extended wall time for normalization and trinity
       73bc9e6b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       db70b578 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       02f65b8e rsem prepare reference separately + job log relative path
       ef38d5f8 First version of RNA-Seq de novo pipeline with normalization, trinity, rsem
       39a95e7e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       74d4ac04 Further development of RNA-Seq de novo assembly
       92aa2512 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       e43a3e2d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       8946fc18 Further development of rnaseq_denovo pipeline
       06fb9bcb Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into rnaseq_denovo
       d2a0332b Import of old deNovoAssembly in rnaseq_denovo branch
       fc8082b0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       052b8041 Minor fix
       eddaed25 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       ab7f6212 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       9eb13833 First draft of de novo RNA-Seq normalization
       f4c998af Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       58732977 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into deNovoRnaseq
       23f34b25 Fixed bug missing '\' before

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      12 commits

       29ff2382 Readme for chipseq pipeline
       6f007a96 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       062d6986 document mugqic pipelines setup using md - correcting bugs
       30364fc8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       fad747c2 document mugqic pipelines setup using md - correcting bugs
       5a0a009b document mugqic pipelines setup using md - correcting bugs
       f771f1af document mugqic pipelines setup using md - correcting bugs
       848249a1 document mugqic pipelines setup using md - correcting bugs
       58dca89d document mugqic pipelines setup using md - correcting bugs
       96522902 document mugqic pipelines setup using md - correcting bugs
       4cf001dc document mugqic pipelines setup using md
       4e9a82a7 Miscelaneous graphs chipseq pipeline

  johanna.sandoval@mail.mcgill.ca <johanna.sandoval@mail.mcgill.ca>      21 commits

       a5f8c705 corrected bug in qcTagsjobid
       d3b16fdd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ecd9434c Added PODs documentation to chipseq pipeline wrapper
       3555ab40 correcting links to project's directories on README.md
       1b94392d link to wiki on general README file
       31db440f BFXDEV-63 Added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       589689e0 added -singleEnd flag when calling RNASEQC if parameter libraryType indicates single end reads
       eb5a068c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       4e1e7df2 added variables for chipseq peaks annotation plots
       ea839b33 corrected bug loading ReadMetrics and read metrics sample/run separator
       9298aae6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       0d69c529 commit chipseq report branch
       95aab1a2 added annotation statistics + graphs to the chipseq pipeline
       94bcb8ba adapt chipseq pipeline to the new resume jobs functionality
       73ee7a7a adapt chipseq pipeline to the new resume jobs functionality- config file
       f1dabdda adapt chipseq pipeline to the new resume jobs functionality
       4fe99531 Adapted chipseq pipeline to resume jobs changes
       9d6f3d88 Transfer chipseq pipeline metrics libraries to the pipeline space
       df472f0a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       2e5f6134 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report
       67ae3798 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into chipseq_report

  Johanna Sandoval <sandoval@ip03.m>      2 commits

       2a5fbe0c document mugqic pipelines setup using md
       6e68220f adding mammouth configuration file for chipseq pipeline

  Julien Tremblay <julien.tremblay@mail.mcgill.ca>      8 commits

       1971ef14 Added memtime to dcmegablast and blastdbcmd. BFXDEV-30
       9f0b66a2 Fixed unwanted mofifs to BLAST.pm. BFXDEV-30
       b430c4ce Modifications to perl packages for Pacbiopipeline. BFXDEV-30
       bb5a19ea Major modifications to the pacbio pipeline. Functional version tested on abacus and guillimin. BFXDEV-30
       a17cd5ac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       97cd35a6 Added readme for pacbio assembly. BFXDEV-30
       6432aada Remove my info fields in ini file. BFXDEV-30
       67c11c1c Added PacBio assembly pipeline libs/wrapper/files, etc. BFXDEV-30

  lletourn <louis.letourneau@mail.mcgill.ca>      19 commits

       4a6d8af2 Version bump 1.1
       ef0d7437 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2541c188 BFXDEV-68 Added mutect to paired pipeline
       d3998b62 Removed gpfs for guillimin
       ed3de9b1 Fixed conflicts
       11d0fe51 Fixed undef on steps that are up2date
       bc778813 BFXDEV-67 use a true csv parser
       bdbb7cee Merge branch 'chipseq_report' of bitbucket.org:mugqic/mugqic_pipeline
       4d7c7fe0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1a8cf594 BFXDEV-45 Fixed single read handling, fixed undefs
       5d280008 BFXDEV-15 Changed the name of lane level bam metrics
       ed46494c Fixed starting dependency and final job output
       4555ab02 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       92520500 Added pre dependency and final sample output
       aa452ff8 Changed BWA version for run processing and added some metrics
       a29bac1e Added runprocessing parameters
       ab659321 BFXDEV-52 DNAseq realignment restart generates duplicate unmapped bam
       05272a81 BFXDEV-48 added missing close BFXDEV-49 added support for alignment+QC only RNASeqQC
       e1ecac6c Version bump

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      6 commits

       da9919d3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fa47c130 replace ; by && - BFXDEV-47 and update the rnaseq ini files
       e0ccc024 dnaSeq.pl correct typo (realign step) : BFXDEV-54
       371c03f4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       799fca20 change/test replacing ; by && un command line
       0a07178e Trimmomatic.pm remove single trimm bug

1.0        Fri Nov 8 15:03:24 2013 -0500        794 commits

  David  Morais <dmorais@ccs.usherbrooke.ca>      2 commits

       26201201 Merged in daveM_denovo (pull request #2)
       0764a256 Merged in daveM_denovo (pull request #1)

  David Morais <moraisd2@ip03.m>      94 commits

       4c13dd19 remove entry
       f458e989 added new entries
       d8e671b1 split butterfly command lines in chunks
       a963e53f creates BLAST best hits reports
       61d9bf3b Fixed bug
       5b9c1605 Modify how the butterfly files are splitted
       0234e161 fixed bug
       5b2bbbc2 modify command line
       f29fe2e7 fixed bug
       f2bf4ab3 Added full path to groupInfo left and right samples
       e8794a0b Added option for de Novo Assembly
       3490e433 Implemented and tested Merge Trimming Remove duplicates de Novo assembly Blast BWA
       8094aca9 fix paths
       191bdcc8 fixed picard version and env variable
       17934d1e change module add trinity
       891bcf9d change module add blast
       60bfa9f3 change module add blast
       78f4ba6e change module add trinity
       59905eb9 fixed typo
       aacc3ec5 fixed typo
       e78b5b69 Modified HtseqCount::matrixMake call. Now group is the last parameter so it can be left out.
       c13eb913 Fix $group problem. Now if group is not specified the variable $group takes $sampleName value.
       dcabd4bf Fixed blast db error name
       592659c3 perldoc
       a94ec1d5 adding scripts needed by de Novo Assembly
       a64a2cee add comments
       289ec2f7 De novo Assembly config file sample. It needs to be modified according to your each project.
       ba64039b mammouth_deNovoAssembly.ini
       befb4244 This is the main de novo RNA assembly pipeline. This is the first commit.
       774254dc modified config variable
       5be92fd5 tidying up
       7d190c5f tidying up
       8d9704b9 tidying up
       03fc63e7 tidying up
       9b75e54f tidying up
       24e31b0a tidying up
       2cb0b70e tidying up
       94a0d328 tidying up
       6a8f8fcf tidying up
       24aa7afa tidying up
       42608a83 tidying up
       2d3f8a9d tidying up
       23e5b5e1 Coded full command line. Added sub contigStats
       d0f849d9 Modified $laneDirectory
       88bb2c75 Modified $laneDirectory
       db5f5063 Modified $laneDirectory in the singleCommand
       bc774be3 Modified $laneDirectory to assembly/ and added $group variable to deal with groups on deNovoTranscriptome assembly.
       55eb46d0 Modified perldoc
       0ad8b58a Modified $laneDirectory
       c21cab50 Modified $laneDirectory from alignment to assembly
       2037ad69 Modified $laneDirectory
       6bb39648 Modified $laneDirectory
       229bba49 Assingned reads/ to $defaultDir
       f0b1c619 Assigned  $laneDirectory to 'reads/'
       eb949206 Modified $laneDirectory = reads/
       921adb05 new line
       b187d45f Output bast results to /db directory
       7c1754a0 change from own branch
       c58b61f3 First commit. Library to create differencial expression analysis.
       4086a94d Modified Doc
       cc24e6b2 First Commit. Library to generate basic statistics on raw read count. Not Yet Ready for use.
       1df2ca70 Library to create generate basic statistics of each sample. First commit. Not ready for use yet.
       edf8193c Added the possibility of looping through more tha one DB. In This case the DB must be passed as an argument to the function.
       3dac1f13 added comment
       fccb0b11 removed quality offset from file name
       db155a07 removed quality offset from file name
       301427d0 added indexing function
       e373c99d add function
       34fd4a44 add function
       986c00cf add function
       7183abe7 modify hash groupInfo
       a77802be fixed input values
       6bae877d BLAST lib, first commit
       78008bd6 modify looping by group option
       998b38a6 tiding up
       6e243f55 fist commit Trinity
       fb5182cf change split setting
       719e16bc added new entries to sampleInfo hash
       4bb98a42 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       644271f7  a simple multiformat file Spliter. First commit
       8d3ccb39 Remove duplicate reads. First commit
       cb8cfdf9 modify to work with nanuq sample sheet
       56d16c84 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       80dd2fd8 First commit
       7fbacf8d first attempt on STD
       736365b0 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       2c33bfc1 Reads the libraries form script_dir/lib
       c3d669c1 added function that allows the script to read libs from the script dir
       d9a378db add header
       e32573a8 remove package folder
       754f187a first commit. Lib to read module list
       16372cb2 first Commit. lib to read config files
       10e24cd0 change repo name
       c56a0af6 adding space.

  eric audemard <eaudemar@imac6-ub.(none)>      5 commits

       f5360846 add software install script : igvtools virusFinder SVDetect add genome install script : virus
       2bd3239b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f55f728b add install canFam3 and new bwa
       7bdcac74 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       270bac92 adding deNovoSV pipeline from the mugqic_eaudem

  Francois Lefebvre <flefebvr@abacus2.(none)>      4 commits

       6f6477c5 febvre <ddcccccZZZZlefebvr@abacus2.(none)>
       a4f60317 Integrated exploratory analysis step
       1ea865bb Changed Rscript -e template
       a3468dcb Fixed perl formatGtfCufflinks.pl call

  Francois Lefebvre <francois.lefebvre3@mail.mcgill.ca>      2 commits

       d312e26c fixed problem when R packages list has duplicates
       3fae3d0a Merged in module_install_scripts (pull request #3)

  Francois Lefebvre <lefebvrf@gmail.com>      53 commits

       2ec42e90 duplicate section names
       2a175b00 duplicate sampleOutputRoot
       e9e50a1c duplciate param in ini
       2abd12c5 duplicate sortSam sections  in ini file removed
       e8f87929 fixed bwa module name dnaseq.abacus.ini
       5b3c1a9b kmergenie install script
       810d0810 updated rnaseq.abacus.ini  wiggle chrsize path
       8321d0f5 Updated R install script
       d3087f74 R install script now installing Vennerable, gqSeqUtils, gqUtils from remotes.
       9caf2022 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       b6591f3f Exploratory analysis step, still to test
       ac2361f9 Added RSEM install script  + updated Trinity one
       789cb0f4 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       6e37c421 no message
       b014a95f Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       7d5e0aaf Changed bwa to mugqic/bwa/0.7.5a, problem between RNA-seQC and 0.6.2
       165e72d6 Updates R/Bioc script, now also installing all org packages
       aa3dcc8c Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       2e48fa55 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       92e196aa no message
       5d60cfa7 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8d4ef871 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       51243690 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       34000ade hg19 install script" add gtf mask for chrM, tRNA, rRNA for cufflinks
       77f6acd4 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       8165a036 Walltimes for htseq, cufflinks, cuffdiff too small  for typical data
       d99c8776 htseq walltime
       03aa1e15 Added 72h wall time align step
       f70ac2f7 Updated bwa install script
       55560916 Added blast module install script
       a4cd0567 Added few more R packages for install list
       baaacb22 Removed MANPATH prepend from mugqic/R defnition
       8eb8a256 chmod in in mugqic_tools.sh install script.
       42d4ef93 <tools,perl-tools,R-tools,java-tools> from the svn now in tool_shed/
       9b773f43 R packages lists updated
       0725e01b Corrected bash_profiel for mammoth and changed R install script to 3.0.0
       83cc87ff Added R package "outliers" to R packages dependencies master list
       5c1d4a6a Updated mammouth $MUGQIC_INSTALL_HOME value in guessing script.
       bea77890 Few change to abacus wall time rnaseq
       0e30bd77 Previous fix to Metrics.pm did not work
       30f0bfb9 minor fix to Metrics.pm (will be depr. anyway at some point). chrome size file for wiggle left to default in default abacus .ini
       16132edc Added -T sort unix sort in metrics.pm to correct problem on guillimin. Added job name suffixes in rnaseq.pl to make job names unique on Guillimin (dependencies)
       60c84ccc Merge remote-tracking branch 'origin/rnaseq_test'
       ca2a4ac0 ..
       59e116eb Added hg19 installation
       c83481d1 Added mm10 installation
       c0dabe43 no message
       c99255ec updated to top hat 2.08 (bug in 2.07) + more genome scripts
       57248b56 Added Eric's abyss install script
       b20104a6 drafted genome installation script
       e82f783e Corrected htseq module by adding module load python. Also started genome setup scripts
       a4c42324 Finished python/numpy script. This script will not be 100% portable, need to set locate BLAS and LAPACK
       6da084c0 Added modules/ directory to hold modules related content

  Joel Fillon <fillon@ip03.m>      4 commits

       c4740583 Added java module in BWA lib + dos2unix rnaSeq.guillimin.ini
       5bf363d1 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       929fb10a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       5f3e86b1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output

  Joel Fillon <jfillon@abacus1.(none)>      4 commits

       a051fa55 Added readRestartFile  function
       a7a14042 Print out MUGQIC command exit status.
       f55bec04 Missing ";"
       a7a2f6d4 Missing ";"

  Joel Fillon <jfillon@abacus2.(none)>      1 commits

       16641e85 Minor misspelling

  Joël Fillon <joel.fillon@free.fr>      27 commits

       67995c4e Added simple Perl script tool_shed/getLogTextReport.pl to create log reports
       f7ddba65 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       3d61e719 Simplified getLogTextReport parameters
       e12fecac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       aac70e27 Print out full job log path + 2 exit code outputs
       d9e2d1aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       75de9310 Fixed lib/SubmitToCluster.pm conflicts
       253bdfbe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       e7d552e3 Added number of jobs in log output
       020e4e94 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       8ac84841 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       4991e383 Added timestamp in job log filenames .o
       ae0b820f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       b42133bb Reorganisation of job logs in an array
       2db5f572 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       5d39c044 Added a few comments
       bc9df34b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       f3998bee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       30cf022c Fixed exit status + walltime
       39079a59 Added getLogTextReport function
       3c10f407 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       964e92d2 Fixed package name
       4c6b06b2 Merged master into logReport
       5e38b430 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       70eff6ed Added readRestartFile log function + Exit Status code in .o
       d7f0f175 Fixed deliverable typo
       5fd50047 added python tools path in mugqic_tools.sh

  Joel Fillon <joel.fillon@mcgill.ca>      8 commits

       81a3ec98 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       e19c7bd6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2128dd21 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       2235fc5b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       dab8f294 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       4cf3157d Check well-defined variables
       7b037050 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c4eaaeb6 dos2unix rnaSeq.guillimin.ini

  Joël Fillon <joel.fillon@mcgill.ca>      95 commits

       a3fddd2b Removed Log lib (now merged in mugqic_tools/getLogReport.pl)
       32f25ea9 A bit of a cleanup
       8dbabd9f Synchronized rnaseq .ini files
       c8fe2ba0 Fixed bug set @INC with relative path to mugqic pipeline lib
       06ca9c5d Moved module and genome files to mugqic_resources repository
       3f115084 Fixed bug missing '\' before
       b1635464 Standardization of rnaseq .ini file for the 3 clusters
       6e7375fe Removed prereq in module install script template
       d76adbf2 Removed prereq + fix R module load to compile kmergenie
       784d0c7d Updated Picard module install script now based on template
       750cacaa Fixed inline comment bug in module install script template + fixed permissions in ea-utils module install script
       6fc665cf Fixed ea-utils compilation error on guillimin
       24ab3fea Minor comment changes
       dd1f7303 Added ea-utils module install script + minor README change
       a8533c2a Renamed MODULE_INSTALL_SCRIPT_TEMPLATE.sh + deleted old archive + minor fix
       e6a1552d Added module file permissions
       7e85ddc4 Created a module install script template; modified tabix install script
       1b2f9bc1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       40b98e7f Deleted redundant modules/add_me_to_bashprofile.sh; added README.txt instead
       e00d0408 dos2unix all pipeline .ini files
       14673934 Updated chipSeq pipeline .ini files with latest module versions
       01330536 Updated rnaseq .ini files with latest module versions
       ce59f4f2 Merged
       e187b2db Minor change in R module install script
       88ebc06e Changed permissions in R module install script
       d5cfef17 Added vcftools module install script
       8a5b3698 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a788a15b Added archive storage + permissions + cleanup in snpEff module install script
       127d3277 Removed unnecessary lib BiocInstaller for R install
       f1bc1a0c Reorder chipseq modules in .ini files
       ecb3c6b2 Removed duplicated modules in chipSeq .ini files
       4280f36c Added archive storage + permissions + cleanup in Trinity module install script
       1976a518 Added archive storage in UCSC module install script + update dnaseq ini file with bwa/0.6.2-ptx renaming
       dfc61b08 Added archive storage in Trimmomatic module install script
       64e9a661 Added archive storage in tophat module install script
       04ba4a94 Added zip archive storage in picard module install script
       d1c2d54b Added archive storage msg in GATK module install script
       eed094d3 Added archive local storage in igvtools module install script
       8bdeb5bc Removed bedGraphToBigWig module install script (now part of ucsc module)
       ff546a71 Added permissions in Tophat module install script
       9c25d167 Updated dnaseq/validation.pl and dnaseq/pairedVariants.pl to comply with new version of SubmitToCluster
       cbaf966d Added permissions in Picard module install script
       7a726e60 Added permissions and cleanup in IGVTools module install script + minor changes
       37095823 Added permissions and cleanup in hmmer module install script
       34df850b Added permissions and cleanup in gemLibrary module install script + minor aesthetic changes
       89ac0f93 Added permissions + minor fix in gatk module install script
       ed0aecac Added permissions and cleanup in fastx module install script
       d7dbb818 Added permissions and cleanup in fastqc module install script
       bd890969 Added permissions and cleanup in exonerate module install script
       ff613795 Fixed permission bug in UCSC module install script
       568ef4ea Fixed another permission bug in module install scripts
       daa7cf20 Fixed bedtools install script bug for archives with different naming system
       c3a3067a Added permissions and cleanup in cufflinks install script
       1cb9cc5c Added permissions and cleanup in bwa install script
       3ef584d4 Updated default pipeline .ini files with new mugqic/ucsc/20130924 for bedGraphToBigWig
       942bc47a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0eac0a13 Removed -j option in make; added kentUtils.Documentation.txt in bin/
       e2253bc9 New install script for UCSC genome browser 'kent' bioinformatic utilities
       7e1b1650 Reorganized install scripts for a5, bedtools, blast, blat, bowtie, bowtie2
       d097e82f Reorganized BEDtools install script
       99f688cf Removed last workdir parameter in chipSeq, rnaSeq; removed commented code; fixed blast+ bug in blast install module
       a985b243 Added shebang in pipeline bash
       27c480f1 Removed leading ':' in chipSeq JOB_IDS lists
       669a97d3 Clarified bash by using job variables + header
       6e72b678 Removed sampleOutputRoot tag in all .ini
       ce93e257 Removed 'mkdir output_jobs' commented lines
       89d60fd8 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       68e948a6 Removed output_jobs subdirs in chipSeq pipeline + added job dependencies in log
       c47fd5d1 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       454e4aba Fancy install script from Johanna
       e9f61b19 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       054fbfb3 Removed final ';' in MACS2 cmd + aesthetic change
       6e049a9b Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8471e385 Removed ';' at the end of Homer cmds
       ae79f343 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       66a91932 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       24459567 Merge branch 'jobs_output' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       310035f1 formatting changes + chipSeq redirect job output to specific location
       830d57b8  ->
       19f116ec Minor formatting changes
       890c9157 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       995a8b41 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       a05a0901 Added java module in BWA lib
       622217d9 Reorganized job_output for DNAseq pipeline
       405043d9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       8f6b4e59 Added tophat 2.0.9
       6ef6c343 Minor formatting changes
       e01c1820 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline into jobs_output
       9e220c5c Minor job parsing change
       4af616a2 Added Job ID number in log
       7ab13210 Redirect job list output into a specific file
       7438d641 Job outputs go into jobs_output directory
       967b8889 Minor coding standardization changes
       472a9b23 Minor error msg edit
       abf1067f Moved getLogTextReport.pl in tool_shed/perl-tools/

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      49 commits

       8e63db4b changed config files to fit the new SubmitToCluster:initSubmit structure
       98b7290a changed config files to fit the new SubmitToCluster:initSubmit structure
       f2c2aa07 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b56095e1 adapters small rna from Trimmomatic-0.22
       61c25df3 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       85c99b15 added module and installation script for gemLibrary
       3998d3d1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       15e6a92c added chipseq.guillimin.ini parameters for hg19
       152e3a76 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7559eafb changed deseq.R : added flag comment.char=" " when reading edger_results.csv
       b641ee28 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e46c9c1e update install genomes script in tool_shed directory
       98499325 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5bde35c prepend perl-tools directory to PERL5LIB
       e8194696 adapted ini file to changes in trimming parameters (headcrop)
       af929cef Changed readstats module, added ReadMetrics.pm library
       ef32d661 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       db61b398 Parsers or trimmomatic and flagstats reads statisttics, create a readstats.csv file
       a8d89175 remove unused function (old tag directory generator)
       dfc0b2e6 compute read stats using the ReadMetrics library
       7faa2252 prepend perl-tools to the PERL5LIB variable
       b4c6417d corrected freeBayes install on guillimin
       79ee5c6a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       e3e85494 Adapt install freebayes for Guillimin
       f46883d6 Generate QC stats using R
       3fa1222a Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       09f7ec24 tag directory generation removed tmp sam file creation - added samtools
       c3ca36b0 tag directory generation removed tmp sam file
       d045ac12 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       8288285f Freebayes install and create module
       237ed033 filter aligned unique reads using samtools
       176e710e Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       df89fcb0 added qc plots, corrected motif weblogo bug: seqlogo available on v 2.8.2
       51077d33 generate chIPSeq QC metrics R script
       6087d823 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       a99e1393 Increase walltime for qctags generation
       8e47a390 Initialize BWA_JOB_IDS per sample, create design 0 directories for MACS and annotations
       143492d5 generate wiggle tracks: remove -i parameter
       9666a94e Validate differences between design file and nanuq project file
       c606dbd9 Validate missing values in design file
       0f2715e4 adapt read statistics to changes in Metrics library
       494d67e3 pbam flag for paired end reads on MACs
       d0cfe12c chipseq pipeline configuration file for abacus
       ed0bf0b7 chipseq pipeline wrapper
       2665867b chipseq pipeline HOMER tags, annotations and motifs
       dde6603e chipseq pipeline MACS2 peak calls
       69e8682f uncommited changes, will control variable names from my scripts
       aa273e99 Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline
       70afa15f transform jobIdPrefix to avoid illegal bash variable names

  Julien Tremblay <jtrembla@abacus2.(none)>      2 commits

       fde3058a modifications to skip unnecessery bam merges.
       cbdfb711 Added dnaclust install trace

  lefebvrf <francois.lefebvre3@mail.mcgill.ca>      3 commits

       4edd2fd4 Fixed bowtietophat module, added cpu param to align, drafted template guillimin
       136c03ff test
       d75425a5 Added Config-Simple and guillimin ini wc

  lefebvrf <lefebvrf@gmail.com>      5 commits

       c7b8ab6d Merge branch 'master' of https://bitbucket.org/mugqic/mugqic_pipeline into rnaseq_lef
       63e674f1 quiet option  is hard coded, makes it hard to diagnose cufflinks… -q could be enabled in previously added option parameter
       bfc8e659 Cufflinks otherOoption parameter, to allow for instance for -M rRNA, MT, etc mask
       62b19285 abacus htseq wall time increased from 3 to 12h
       4ed3cbaa unused raw_count folder was created

  lletourn <louis.letourneau@mail.mcgill.ca>      174 commits

       b07568df Added restart implementation
       6633b73b Fixed pairedVariants with new structure
       50e3e85f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       29fc49a0 Added release instructions
       58333308 Removed deprecated code.\nAdded a version
       6294737a Fixed step order bug
       3a504533 Fixed region bug from previous commit
       490b11e0 Added a way to skip trimming. Added a mem example
       a385b07a Added a way to skip trimming
       c60ac775 Merged master
       a535288d Added lib barcode to lane data
       441b96be Added missing ;
       8fa3ded4 Updated old code
       399a2d96 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d9e7894f Fixed mem bugs
       1c2a330c Fixed badly named CFG key
       df22d10a Fixed paths for install. Crucial for these to work since they hardcode the paths
       fbe736a7 Merged Master, some fixes and new dnaSeq report
       9bcc238d Updated dbsnp
       623de3ef Amos install
       0ffb61eb Fixed region for torque
       abd09408 Amos install
       0370f0ef Added split for stats
       508f1910 Fixed if test
       2202926a Fixed dependant files
       57f85c2c Fixed some bugs found while testing public datasets
       e9d5b4ee Fixed dir search order
       dab8a760 Put appending and cleaning of dones in SubmitToCluster. More central, easier to maintain
       2cd0d6b4 Fixed to support any year 2XXX
       b9aafe9b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       1eb4c137 Merged master with outputDir fixes
       ac706266 Merged HEAD, removal of job output dir
       a4c21625 Change delimiter for regions, msub doesn't like colon
       d1c131e2 Don't generate the pileup, it takes way too much space
       9805b640 Added hg1k install and GNU parallel
       2be9b17f Forced default picard buffer to the abacus optimal one
       6230215c Finished updated the rnaSeq.pl pipeline
       431b192f Merged master, added java module, dos2unix rnaSeq.ini
       b2898cb0 Merged Master
       cf0fca30 Finished implementing reset in dnaSeq
       a2fc84a5 merged
       fcd41d22 added missing calls for output dir
       708959e7 added sampleOutput param to inis
       ef0f696c Use the reight blat tools
       5fdb7f3a installed v4 preview
       83368651 Added beagle 3.3.2
       ba10cb48 Refactored in the new Job object and input output file tests for restart
       ea2786e6 Updated the way snp calling works to support, nb jobs instead of window.
       f28a83a1 Removed dup index generation since we added recalibration
       4a6c7013 Fixed the way output directories are initialized
       10f160f8 Updated version
       0f80b863 Added pigz
       b9e96e54 Merged master
       98455059 Added PipelineUtil pacakge to all.\nCompleted BLASt implementation
       d218799b Implemented per-job indel realigner Implemented multi sample snp calling Implemented per-job/per chromosome snp calling
       3359fa30 Use ini for igvtools version Updated igvtools version for correct exit code
       a3860f4e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       608ae55e Fixed already processed lanes.
       66430dbc Version bump fixed execute bit
       4ed8abc0 Trial implementation of new job restart flags
       2e09503d Implemented restart job handling
       9912c46b Added link to ptx patch
       c9841173 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       daa8f5ab Added and updated modules
       5448d090 Added tabix
       35fa5c7b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       75c9afd3 Fixed many bugs
       ae121477 Fixed many bugs
       fb3af3a2 Fixed hash generation
       42f88d87 Changed trimming
       e7eb18f8 Fixed index vs ref problem
       1dbc2b88 Fixed index vs ref problem
       84c1895a Split metrics...we probably should split it further.
       7747bcd2 Added TRINITY_HOME
       70cd9b83 Added a parsable file containing trimming stats
       d032c5ae Fixed default params
       9ea12118  Changed metrics position in steps
       96ca681e Fixed parameter passing.
       ea8cecc8 Fixed problems in RNASeq denovo pipeline
       5cf88989 Added support for non-multiplex lanes....
       1d8a1903 Fixed trimmomatic bug Added usage to sampleSetup Fixed some ini params in rnaseq
       a4054926 Added options for skipping samples
       4d7e1fac Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d06a44fc Fixed many issues and params for the RNASeq denovo
       162142ab Format stats from new format
       7b619615 Version bump
       972fa53a Added missing file warning
       42b94f21 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       c347b1ff Made many fixes to the denovo RNA Assembly pipeline
       a9959680 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       05042bf8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       7ec95006 Modified permissions
       ceea8d29 Fixed mpileup bam location
       58152cf8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       49b14098 Fixed parallele threads in single reads
       e311d8d2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d7e10f5a Fixed bwa mem RG bug.
       8fa53784 Updated version
       01e0dbde Added gatk to module list
       8ff53e4c Fixed special single issue
       1c6edb26 Added recalibration
       36cbabd1 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f0276798 Fixed generation for single ends
       660ceda4 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       8ef1373e Fixed validation for new project layout Added BWA mem support
       89b8a332 Added snpEff and mutect Changed mode
       070bceef Fixed output redirection
       add31c71 Fixed multi runId resolution
       351be7c2 Cleaned up downloaded genomes
       40edd6e3 Sample setup scripts that works with miseq and hiseq and fetches the nanuq project by ID
       69f8ad49 Added gallus gallus 3
       a408a59e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b3ee3404 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f461b67a Fixed default ref.
       210079be Use parallele bwa by default
       1d973c60 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e83d0688 Fixed libs and dna pipeline to use official project hierarchy
       b0a33ed2 Fixed cleanup code
       2af9a96e Updated the ray module
       643a97b2 Added the exonerate module Fixed the ray module
       5c1c6197 Fixed append bug
       41a888f6 Change file mode
       12bdf757 Fixed configurable module load to trimmomatic Fixed pass output dir, not hard coded in trimmomatic
       3f20c53b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       46ac74f6 Added missing rawReadsDir
       6f4161d6 Partially fixed raw_read support for dnaSeq Fixed markdup call in dnaSeq
       f861134b fixed conflicts
       c973658f Fixed whole genome coverage
       32b1e169 Added module versioning, fixed paired pileup bug
       94c6f972 Added SV to the paired pipeline
       2e99fd6c Added snv calling
       2bfad558 Rolled back changes because it broke existing pipelines
       0bc8182a Added SV to the paired pipeline
       d67cedde Added DNAC support
       79e909f7 Added DNAC support
       b9b11aab Merge branch 'mugqic_lletourn'
       0cceb490 Use job ids Added variant calling Added more flexibility when calling modules
       18a60abe Added thresholds to genomeCoverage
       a70bf88b Merge branch 'master' into mugqic_lletourn
       275d08aa Fixed usage
       cee2b182 Added for snp calling
       4e930431 Added metrics and steps to the dnaSeq pipeline
       8222088c Change job dependencies from names to ids.
       77a4c618 Fixed 2 paths bugs
       30396f0f Don't use alternate contigs
       0c43b15f Changed executable attribute
       7b04c504 Keep dictionary ordering
       f08a3ea6 Fixed usage message
       41172652 Added paired snp variant calling
       a5d2e160 Added more flexibility with tmp.
       67b2b059 Completed first half of the validation pipeline Added adapters file for common paired adapters
       d0df7f45 Added metrics to validation Added IGV tootls Added targetCoverage
       d66f466a Fixed header
       466071a3 Fixed issues with multi single+pair
       df027659 Added Validation pipeline
       b48e4b56 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       01ce3d46 Fixed index generation
       0f46ad3c Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       0067ebf0 Fixed timestamp checking
       08d14e01 Added trim counts statistics
       09231132 Fixed Trimmomatic MAJOR bug Added more settings to trimmomatic
       275cd84a Use the right path for trimmomatic
       f652072f Fixed output job directory naming
       cc073638 Fixed bad realigner output file name
       7c15b689 Merge branch 'master' of ssh://127.0.0.1:2222/mugqic/mugqic_pipeline
       992ef66e Fixed bad argument
       c4479ab6 Added guillimin conf file
       5411b170 Fixed job dependency code.
       53a75c0f Finished bam generation pipeline
       faf8ed15 Merge branch 'master' of ssh://bitbucket.org/mugqic/mugqic_pipeline
       4e306743 Added realigner
       eda6d125 Fixed uninitialized bugs
       2b178064 First pass for dnaSeq pipeline
       5b1fc8bf Fixed typo

  Louis Letourneau <lletourn@abacus1.(none)>      1 commits

       9286b83a Fixed many bugs

  Mathieu Bourgey <bourgey1@ip03.m>      2 commits

       3b2af18a RNASEQDN: modify trinity to check for previous assembly + variable name change
       d49a53aa TOOLS: chenge permission of non-executable scripts

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      26 commits

       58770fe6 RNASEQ: change merge test
       73407b03 RNASEQ: change merge test
       54a772aa Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       22f64a30 RNASEQ: add mkdir metrics in step 2
       be482de9 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       ee774f9b goseq.R: check and remove results that are not reported in the GO.db R package
       220e9fba lib GqSeqUtils report call change
       aa69e7ee Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       08a72c99 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       20ce1c96 formatDenovoCombinedGTF.py update
       3183fd0a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f6ce3df5 RNAseq update
       c833eee7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d1e8dcce RNASEQ update unstrand wiggle bug correct
       49d7bbc8 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e7dbfbc8 RNASEQ correct stranded wiggle array assignation
       352f1ed8 add the java module call at Metrics::rnaseqQC
       7825d84a add the java module call before at each picard function
       aeb1c37e RNASEQ: add mamouth ini file
       2f858ec8 MODULES - add temporary download folder in several module install scripts
       e63d6498 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       44263b2b rnaseq: resolve dependency pb
       9873fc20 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       89fbe54e RNAseq modify some output location for the reporting
       27bbc10c STATS: correct metrics:readstat for using output trimmomatic
       772a5e5a RNASEQDN: old mamouth ini file by the new one that fit changes

  mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      96 commits

       ad595721 RNAseq update cuffdiff new dependency update
       278c9d25 RNASEQ: change variable nemae of bowtie fasta
       483e00dd RNAseq update metrics stats
       8b02e13a RNAseq update trimming metrics
       21493623 RNAseq correct DESseq wrong dependency
       398a1211 RNAseq upstae trimming merge stats; update edgeR
       b07481ea rna update
       bed41b69 rna update
       4c45b56a remove tmp test for maty in dnaseq.pl
       45c8ba26 rna update
       ab2a1499 R-tool: change saturation graph format
       7d052f14 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d0bd4998 RNASEQ: add format protability to goseq
       1cd3ae08 RNASEQ: add format protability to goseq
       2c7b901b RNASEQ: add format protability to goseq
       be3833aa RNASEQ: add format protability to goseq
       3289aac7 RNASEQ: add format protability to goseq
       971938cf RNASEQ: add format protability to goseq
       6cc1155b RNASEQ: add format protability to goseq
       9f842076 RNASEQ: rnaseq.pl update
       e31e700f RNASEQ: rnaseq.pl update
       fc1faad4 RNASEQ: add format protability to goseq
       18ffd9a5 RNASEQ: add format protability to goseq
       e725e01c RNASEQ: edger update
       7fbbb5a6 RNASEQ: goseq update
       c6dbccbd Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       73a94bf1 RNASEQ: goseq update
       06b2e4af RNASEQ: goseq update
       f93dba0d RNASEQ: allow non native GOseq analysis
       656c63c7 RNAseq: cuffdiff result fillter now include in the merge with fpkm
       02e9ef8b RNAseq: cuffdiff result fillter now include in the merge with fpkm
       bc787bd5 RNAseq update
       49ab91f4 RNAseq update
       e6515d3c change mugiqc_tools.install script to avoir facl conflict
       9676b39a update tools changes
       287a8e2d update tools changes
       e35a443d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       58ad96b3 RNASEQ: resolve splited metrics dependency
       4008a3d9 Rtools edger old index name  not suported by the actual edger version now use
       483781a0 RNAseq remove fpkm stats as they are done also in rna-seqc
       c47230d6 RNAseq correct bugs - see BFXDEV-20 for details
       8b614c0f RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence 2
       b5eee70a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       79804b1c RTOOL: mergeCuffdifRes bug correct && adapt the Rnaseq script in consequence
       3b0ce988 RNASEQ: bowtie reference the new code should now be portable
       4b6bf8e1 RNASEQ: allow readstat on single library
       d079f598 RNASEQ: correct readstat output format
       a61e0050 RNASEQ: metrics changes bug correction & Cufflinks correction
       28285f9c Metrics add  argument
       a1032f16 RNASEQ - correct typo in the mergebam step
       68b309e6 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       318f27b5 TOOLS remove type blatbestHit.awk
       30a3793e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       5c8851c6 TOOLS adding blastbestHit.awk blatbestHit.awk
       e8481fd5 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       26152f77 RNASEQ: rnaSeq.abacus.ini correct rnaQC fasta variable name
       9bcdcf90 TOOLS: gtf2tmpMatrix.awk - adapt for line witout gene name
       57441b61 RNASEQ : adjust saturation plot title
       f05e6a74  RNASEQ : add saturation thread Number in the abacus .ini file
       1d072b2c Tools : make gtf2tmpMatrix.awk portable for various type of GTF file
       fbf1e9a6 tools : adding gtf2geneSize.awk and gtf2tmpMatrix.awk which were initially in tools but not anymore in the module
       2de9a907 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3ec125de RNASEQ : allow running without reference GTF .2
       a19c305a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a48116df RNASEQ : allow running without reference GTF
       4c9f8033 pairedVariant - finish controlFREEC allow either bam or mpileup input everthing is specified in the ini file -  BFXDEV-7
       4a2439ba R-tools change rpkmSaturation the count file need to be tab separeted
       93eda90e R-tools correct the polting bug of the rpkmSaturation - BFXDEV-3
       c67cdd06 rnaseq - replace my email by the  MAIL variable
       b1039a26 pairedVariant - add Control FREEC lib
       75351f19 update rnaseq.pl for single mode
       d1557826 update rpkmSaturation to the last version
       67148739 correct formatGTFCufflinks.pl
       2220b26d Merge branch 'master' into mugqic_mathieu
       8aef01ea Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       7ff4088e correct rnaseq ; remove special caharacter in subimti to cluster job ID varable
       39e7fd1d correct cuffdiff deNovo
       9db7cdf0 correct conflict btw mugqic_mathieu and master
       314eadc5 correct small bugs
       0307f257 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       17dadbe1 correcting pindel in lib/pindel.pm and in pairedVariant.pl form dnaseq pipeline
       39277436 remove python module call in htseqCount::readCountPortable
       ed644b5f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       fcc0e410 Merge branch 'master' into mugqic_mathieu
       d7d4284a correcting last bug in rnaseq
       4cb464e7 debug rnaSeq.pl
       7725d781 trimmomatic add java module loading
       1c9c334e submitToCluster add default  value =  if undef
       70afb543 rnaSeq debug
       2a78c375 resolve conflicts
       3c840f21 Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       208d0cc8 debug
       cc57bb26 update
       8889bf16 rnaseq test debug; submitcluster modification; trimmomatic modifications
       bd4afe16 debugging RNAseq
       249ded65 debugging RNAseq

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      129 commits

       a51199f5 samtools allow pileup with nb de region = 1 => pas de region
       ea59d1fe Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       3e902cc5 dnaseq.pl && samtools allow pileup with nb de region = 1 => pas de region
       0ff4f488 upadte rnaseq.pl: missing convertion to new job object for stranded wiggle printToSubmit calls
       ae5a4a80 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f79aadbf dnaseq annotation bug correct and ini update
       d4e4586d rnaseq.pl: correct exploratoiy dependency if start at step 9
       1b1abe4d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       33df28ef dnaseq correct metrics restart bug & update dnaseq ini
       6db6f703 remove tools_shed for the new tool repo
       c653f58f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       e1bc2857 metrics/SNVGraphs add input output test for restart
       72b31c8e update Tools
       0e1b10e7 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b7584045 update cuffdiff stranded
       af5dba27 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       371205c4 switch ToolShed to Tools
       c0c1d65c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       665a0986 RTOOLS update
       3887f719 snvGraphMetrics.R  update
       38b94dff DNAseq.pl update for metrics
       47321095 DNAseq.pl update for metrics
       daceb8ab Merge branch 'dnaSeq_Report'
       c82147af DNAseq.pl update for metrics and report
       b31654c6 DNAseq.pl update for metrics and report
       75f28e5b DNAseq.pl update for metrics and report
       425b21e2 add tool_shed/python-tools/vcfStats.py
       ac96abbf remove mater to branch conflict
       d133749a rnaSeq.pl edited online with Bitbucket
       23e60fc9 remove conflict between danSeq_report and master
       e67ddd8d add tool_shed/tools/splitSnpEffStat.awk
       e4cfaee2 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f5ce57fe GENOMES: add Tetrahymena_thermophila.sh
       5f468e6d Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       4bd157a9 sampleSetup.pl patch
       a2dbaccb samplesetup correction
       f4ba0956 Trimmomatic.pm edited online with Bitbucket
       00d0af02 merging matser in dnaseq_report
       7ff833ed rnaSeq.pl edited online with Bitbucket
       4556bdab DNArepport update
       91d8b3cc goseq.R edited online with Bitbucket
       e2148e76 goseq.R edited online with Bitbucket
       5a202eba DNAreport update
       a7a6d130 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       cc7af8de TOOLs gtf2geneSize.awk protability
       05fdcce6 RNAseq update
       8b603f9a RNAseq update
       58cc9c62 RNAseq update
       dcf225db Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b69a50b8 RNAseq correct htseqcount for stranded RNA
       f4520100 RNAseq Strand-specificity correction
       7754be2a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       bc5c5a51 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       a5450eca TOOLS rpkmSautrationadd more flexibilty on the file format
       4960f3ea RNASEQ: add output directory creation in fpkm folder
       904595ff General : Get back the change lost in commit 0007ea2
       84818acd Revert "RNASEQ: correct dependency issue"
       d45c5d77 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       d3b3bb1d Revert "RNASEQ: correct dependency issue"
       7aabc009 save before revert
       0007ea2b RNASEQ: correct dependency issue
       a3c5a79f Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       eb6d023f RNASEQ : add deliverable + test exploratory
       f62067d7 Merged in rnaseq_lef (pull request #4)
       63a01633 RNASEQ - GqSeqUtils.pm update
       23099ae8 RNASEQ modify exploratory for the pull request #4
       80da0567 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       af49cae6 RNAseq update
       c258e2c6 REDO lost commit: rnaSeq.pl use absolute path of design file
       7e28eddd SubmitToCluster.pm edited online with Bitbucket
       79e8b607 add log.pm lib
       19289aaa add logfile for cluster output file path
       ccf1fa8c RNASEQ: rnaseq.pl update
       5e8fcfca Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       28964581 RNASEQ: make design file path an absolute path to conflict if the working directory does match the relative path of the file given in argument
       363892b2 rnaSeq.pl edited online with Bitbucket
       6435ebeb RNASEQ: rnaseq.pl && cufflinks.pm update
       3e54a64e SubmitToCluster.pm edited online with Bitbucket
       9ac95ffc RNASEQ: rawCountMetrics update
       509b645f RNASEQ: rawCountMetrics update
       c484eed2 Merge branch 'logReport' of bitbucket.org:mugqic/mugqic_pipeline into logReport
       6e725ad8 Merge branch 'master' into logReport
       a7faedc6 RNASEQ: mv dgeMetrics (step 10) as rawCountMetrics (step 8)
       04752e23 add log.pm lib
       3fcc2838 add logfile for cluster output file path
       c527ebc8 RNAseq: change wingzip call in the metrics
       cbbd97ab RNAseq update rnaseq.pl
       f4e32151 RNAseq: update rnaseq.pl
       2a8c65c6 RNAseq: update rnaseq.pl change matrixreadcount form step 10 to 7
       b43088aa RNAseq update rnaseq.pl
       f672d802 RNAseq: update rnaseq.pl
       5535e508 TOOLS: add formatDenovoCombinedGTF.py & RNAseq: update cuffcompare
       a348a62d RNAseq: update rnaseq add cuffcompare
       cb485b8e rnaSeq.pl edited online with Bitbucket
       177aa347 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0c0f4a12 TOOL: Adding thepython tools folder and the getRefBinedGC.py script in it
       77b49e29 Merge branch 'mugqic_mathieu'
       1ebe60bc Merge branch 'mugqic_mathieu' of bitbucket.org:mugqic/mugqic_pipeline into mugqic_mathieu
       22ad9e5d update RNA
       36d735db RNASEQ: add correct cpu request for alignment
       e1718619 correcting sv code small bugs
       553a3b00 finish breakdancer filter and  add pindel to dnaSeq
       ce1cd917 start adding breakdancer filter and pindel to dnaSeq
       59fa224a resolve Htseqconfilct
       0ad63502 merge with mugqic_mathieu
       c1a278a7 debug the code
       1498b149 Merge branch 'master' into mugqic_mathieu
       c2b1cc49 resolve lib/SubmitToCluster.pm conflict
       4df938ad resolve HtseqCount conflict
       3268ff41 try to update local master to real master
       97b44321 everyhting except delivrable are done; debugging
       a0a55447  goseq done; 1st test on going
       eb5785ed metrics: done; General module nomeclature conversion
       5efae813  metrics: updates; SAMtools: add viewFilter sub function ; tophatBowtie: bug correction and simplification
       2e48636d  metrics updates
       be331d2e matrix done; DGE done
       4223cc79 fpkm done; DTE done ; DGE on going ; metrics on going
       08778dcc correcting wiggle; fpkm and rawCount done ; DTE on going
       1ed2d7f2 wiggle done
       cd318c11 update metrics ini.file and few correction tophat/botwite
       052889d7 merging ok ; metrics started; new functions and changes in the picard lib
       2f6709a9 merging done ; metrics started
       66ea6a47 update to the master branch
       7166f745 Tophat-botwie done; merging started
       f9662b88 uppdating TopHat lib and rnaSeq.abacus.ini
       634a238f uppdating my branch
       9f81ed08 Starting bowtie/tophat lib
       43c3e7a9 Global dependency and RNAseq Triming and submit with working directory argument
       d6a6324f starting the rnaseq pipeline

  Mathieu Bourgey <mbourgey@abacus2.(none)>      1 commits

       f36f0e21 rnaseq.pl debug; Metrics debug; Picard debug; SaMtools debug; SubmitToCluster debug; TophatBowtie debug

  Maxime Caron <mcaron@abacus1.(none)>      2 commits

       4ad23e11 test
       0b8c7fb8 test

  Pascale Marquis <pmarquis@abacus2.(none)>      5 commits

       168ea11c Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       221a12ef Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       b91d4dad Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0b27236e zorro.sh
       9f170dcb Tetrahymena_thermophila.sh

