12 tags, 2862 commits

HEAD        Mon Feb 8 12:02:54 2016 -0500        0 commits

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

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       70099f3 rnaseq - include correlation matrix manual estimation using an utils.py new function
       fd53801 RNAmetrics - include ini_section argument to some bfs picard function to allow to use them several time with different parameter in the same ini
       e7c5a01 rnaMetrics - remove confilct pipelines/rnaseq/rnaseq.base.ini
       9ad7fbb rnaMetrics - update pipelines/rnaseq/rnaseq.py pipelines/rnaseq/rnaseq.base.ini
       9d3b87d RNAseq - metrics update ini
       1913bc1 RNAseq - remove conflict
       deca058 RNAseq - metrics RNA - update base ini
       2fc2006 RNAseq - remove rnaseqc; add picard_rna_metrics ; partial add estimate_ribosomal_rna - BFXDEV-345

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

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       12ff20e Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       233244b RNAseq - correct stdin input issue of htseq-count - BFXDEV-318
       fa33528 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       47c4cc8 RNAseq - htse-count: pipe samtools view -F 4 output in htseq-count instead of using the bam to remove error due to unmapped reads - BFXDEV-318

  mathieu bourgey <mathieu.bourgey@mail.mcgill.ca>      2 commits

       7ec0344 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       1905eb6 COMMON - add fastq2 = None in sam_to_fastq pipeline step wehen the read are single - BFXDEV-321

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

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      8 commits

       fcd7134 pull before pushing Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       05b352d  PacbioQC - Fix code vs ini section name - BFXDEV-315
       d085269 RNAseq - correct discrepency in hsteq_count ini calls - link to BFXDEV-312
       6773ef7 pull before pushing
       d5d82f0 RNAseq - correct errounous section header in ini files - BFXDEV-312
       98a456a Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipelines
       fdd92aa RNAseq -update guillimin cluster ini file - BFXDEV-307
       b4a12a5 RNAseq - fix report dependencies - BFXDEV-306

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
       98f2dd8 Updated dnaseq dbNSFP2.0.txt to dbNSFP2.4.txt.gz
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

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      7 commits

       14060c1 Remove conflicts in modules/mugqic_tools.sh and modules/dev/star_dev.sh
       cb81397 update modules/mugqic_tools.sh to 1.10.4
       b2db257 remove decrepated python module script and add a new one
       16d3068 replace dev module in module/dev/ cufflinks_dev.sh  star_dev.sh
       fd6119e Merge branch 'master' of bitbucket.org:mugqic/mugqic_resources
       85a4523 add genomes/oryCun2.sh  modules/cufflinks_dev.sh  modules/star_dev.sh
       b35f6de up-date mugqic_tools.sh version to 1.7

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

  Mathieu Bourgey <mathieu.bourgey@mail.mcgill.ca>      4 commits

       3c0fb6e update ini to fit the new mugqic_tools tag 1.10.4 - BFXDEV-275
       e5bd476 correct dnaseq bwa samnse dependency bug - BFXDEV-253
       6fdf802 Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       f3abd8c correct cuffdiff input double array issue when checking the job object is up to date

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

  Johanna Sandoval <johanna.sandoval@mail.mcgill.ca>      1 commits

       6a2eca8 README.md edited online with Bitbucket

  johanna_sandoval <johanna.sandoval@mail.mcgill.ca>      5 commits

       90ac5a4 BFXDEV-133 incompatibility between /usr/bin/perl and mugqic/perl/5.18.2. Added perl HOMER_HOME/bin/ to the program execution in peak annotations and motifs usign Homer
       5e62ee0 BFXDEV-133 detected bug dependencies between trimming and alignment chipseq pipeline
       972c98b Merge branch 'master' of bitbucket.org:mugqic/mugqic_pipeline
       0fe7a1c BFXDEV-133 updated software, parameters, corrected bugs in chipseq pipeline for guillimin phase2
       2401cf3 BFXDEV-133 adapt chipseq pipeline and configuration file to guillimin phase2

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

