.. _docs_gp_ampliconseq:

.. spelling:: 

      Amplicon
      Amplicons
      amplicon
      amplicons
      seqs
      PCoA
      heatmap
      metagenomic
      oligonucleotide
      rRNA
      germline
      img
      biom
      otus
      subsamples
      config
      param
      chao
      shannon      
 
Amplicon Sequencing Pipeline
============================

Amplicon sequencing (ribosomal RNA gene amplification analysis) is a metagenomic pipeline. It is based on the established `Quantitative Insights into Microbial Ecology <https://www.ncbi.nlm.nih.gov/pubmed/22161565>`_ (QIIME) procedure for amplicon-based metagenomics. It assembles read pairs using `Fast Length Adjustment of Short Reads <https://www.ncbi.nlm.nih.gov/pubmed/21903629>`_ (FLASH), detects chimeras with `UCHIME <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150044/>`_, and picks operational taxonomic units using `VSEARCH <https://www.ncbi.nlm.nih.gov/pubmed/27781170>`_. Operational taxonomic units are then aligned using `PyNAST <https://www.ncbi.nlm.nih.gov/pubmed/19914921>`_ and clustered with `FastTree <https://www.ncbi.nlm.nih.gov/pubmed/19377059>`_. Standard diversity indices, taxonomical assignments, and ordinations are then calculated and reported graphically.

.. contents:: :local:

----

Introduction
------------

Amplicon sequencing is a highly targeted gene sequencing approach used to analyze genetic variation in specific genomic regions. Amplicons are Polymerase Chain Reaction (PCR) products and the ultra-deep sequencing allows for efficient variant identification and characterization. Amplicon sequencing uses oligonucleotide probes that target and capture genomic regions of interest and then uses next-generation sequencing techniques. 

**Uses of Amplicon sequencing**

#. Diagnostic microbiology utilizes amplicon-based profiling that allows to sequence selected amplicons such as regions encoding 16S rRNA that are used for species identification. 

#. Discovery of rare somatic mutations in complex samples such as tumors mixed with germline DNA.

See More Information section below for details. 

----

Version
-------
::

  3.1.4

For the latest implementation and usage details refer to Amplicon Sequencing implementation `README file <https://bitbucket.org/mugqic/genpipes/src/master/pipelines/ampliconseq/README.md>`_ file.

----

Usage
-----

::

  ampliconseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                      [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                      [--json] [--report] [--clean]
                      [-l {debug,info,warning,error,critical}]
                      [-t {qiime,dada2}] [-d DESIGN] [-r READSETS] [-v]

**Optional Arguments**

::

  -h                        show this help message and exit
  --help                    show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                            config INI-style list of files; config parameters are
                            overwritten based on files order
  -s STEPS, --steps STEPS
                            step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                            output directory (default: current)
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                            job scheduler type (default: slurm)
  -f, --force               force creation of jobs even if up to date (default:
                            false)
  --json                    create a JSON file per analysed sample to track the
                            analysis status (default: false)
  --report                  create 'pandoc' command to merge all job markdown
                            report files in the given step range into HTML, if
                            they exist; if --report is set, --job-scheduler,
                            --force, --clean options and job up-to-date status are
                            ignored (default: false)
  --clean                   create 'rm' commands for all job removable files in
                            the given step range, if they exist; if --clean is
                            set, --job-scheduler, --force options and job up-to-
                            date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                            log level (default: info)
  -t {qiime,dada2}, --type {qiime,dada2}
                            AmpliconSeq analysis type
  -d DESIGN, --design DESIGN
                            design file
  -r READSETS, --readsets READSETS
                            readset file
  -v, --version             show the version information and exit

----
 
Example Run
-----------

TBD - Link below needs to be updated - not listed in current README.md file

You can download `amplicon sequencing test dataset <http://www.computationalgenomics.ca>`_ and run the following commands:

::

  TBD

---- 

Pipeline Schema
---------------
Figure below shows the schema of the amplicon sequencing protocol - QIIME type. You can refer to the latest `pipeline implementation <https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/GenPipes_ampliconseq.png>`_ to download a high resolution image of the same.

.. figure:: /img/pipelines/ampseq.png
   :align: center
   :alt: amplicon schema

   Figure: Schema of QIIME Amplicon Sequencing protocol

The following figure shows the pipeline schema for `DADA2 Pipeline`_ type of amplicon sequencing protocol. 

TBD

 - Refer to the high resolution image in GenPipes sources
 - currently unavailable.

.. figure:: /img/pipelines/ampseq2.png
   :align: center
   :alt: dada2 ampseq 

   Figure: Schema of DADA2 Amplicon Sequencing protocol

----

Pipeline Steps
--------------

The table below shows various steps that constitute the Amplicon Sequencing QIIME and DADA2 type genomic analysis pipelines.

+----+--------------------------------+---------------------------------+
|    |  *QIIME sequencing Steps*      |   *DADA2 sequencing Steps*      |
+====+================================+=================================+
| 1. | |trimmomatic16S|               | |trimmomatic16S|                |
+----+--------------------------------+---------------------------------+
| 2. | |merge_trimmomatic_stats16S|   | |merge_trimmomatic_stats16S|    |
+----+--------------------------------+---------------------------------+
| 3. | |flash_pass1|                  | |flash_pass1|                   |
+----+--------------------------------+---------------------------------+
| 4. | |ampliconLengthParser|         | |ampliconLengthParser|          |
+----+--------------------------------+---------------------------------+
| 5. | |flash_pass2|                  | |flash_pass2|                   |
+----+--------------------------------+---------------------------------+
| 6. | |merge_flash_stats|            | |merge_flash_stats|             |
+----+--------------------------------+---------------------------------+
| 7. | |catenate|                     | |asva|                          |
+----+--------------------------------+---------------------------------+
| 8. | |uchime|                       |                                 |
+----+--------------------------------+                                 |
| 9. | |merge_uchime_stats|           |                                 |
+----+--------------------------------+                                 |
| 10.| |otu_picking|                  |                                 |
+----+--------------------------------+                                 |
| 11.| |otu_rep_picking|              |                                 |
+----+--------------------------------+                                 |
| 12.| |otu_assigning|                |                                 |
+----+--------------------------------+                                 |
| 13.| |otu_table|                    |                                 |
+----+--------------------------------+                                 |
| 14 | |otu_alignment|                |                                 |
+----+--------------------------------+                                 |
| 15.| |filter_alignment|             |                                 |
+----+--------------------------------+                                 |
| 16.| |phylogeny|                    |                                 |
+----+--------------------------------+                                 |
| 17.| |qiime_report|                 |                                 |
+----+--------------------------------+                                 |
| 18.| |multiple_rarefaction|         |                                 |
+----+--------------------------------+                                 |
| 19.| |alpha_diversity|              |                                 |
+----+--------------------------------+                                 |
| 20.| |collate_alpha|                |                                 |
+----+--------------------------------+                                 |
| 21.| |sample_rarefaction_plot|      |                                 |
+----+--------------------------------+                                 |
| 22.| |qiime_report2|                |                                 |
+----+--------------------------------+                                 |
| 23.| |single_rarefaction|           |                                 |
+----+--------------------------------+                                 |
| 24.| |css_normalization|            |                                 |
+----+--------------------------------+                                 |
| 25.| |rarefaction_plot|             |                                 |
+----+--------------------------------+                                 |
| 26.| |summarize_taxa|               |                                 |
+----+--------------------------------+                                 |
| 27.| |plot_taxa|                    |                                 |
+----+--------------------------------+                                 |
| 28.| |plot_heatmap|                 |                                 |
+----+--------------------------------+                                 |
| 29.| |krona|                        |                                 |
+----+--------------------------------+                                 |
| 30.| |plot_to_alpha|                |                                 |
+----+--------------------------------+                                 |
| 31.| |beta_diversity|               |                                 |
+----+--------------------------------+                                 |
| 32.| |pcoa|			      |                                 |
+----+--------------------------------+                                 |
| 33.| |pcoa_plot|                    |                                 |
+----+--------------------------------+                                 |
| 34.| |plot_to_beta|                 |                                 |
+----+--------------------------------+---------------------------------+

----

.. include:: steps_ampseq.inc

----

.. _More Information Ampliconseq:

More information
-----------------
For the latest implementation and usage details refer to Amplicon Sequencing implementation `README.md <https://bitbucket.org/mugqic/genpipes/src/master/pipelines/ampliconseq/README.md>`_ file.

* `Amplicon sequencing techniques <https://sapac.illumina.com/techniques/sequencing/dna-sequencing/targeted-resequencing/amplicon-sequencing.html>`_

* `Amplicon Sequencing Primer <http://apc.ucc.ie/pdf_old/Amplicon%20Sequencing.pdf>`_

* `High-throughput amplicon sequencing <https://www.biorxiv.org/content/10.1101/392332v2>`_.

* `Trimmomatic - flexible trimming <https://academic.oup.com/bioinformatics/article/30/15/2114/2390096>`_.

----

.. The following are html links used in this text

.. _DADA2 Pipeline: https://benjjneb.github.io/dada2/tutorial.html

.. The following are replacement texts used in this file

.. |trimmomatic16S| replace:: `Trimmomatic16S Step`_
.. |merge_trimmomatic_stats16S| replace:: `Merge Trimmomatic Stats`_
.. |flash_pass1| replace:: `Flash Pass 1`_
.. |ampliconLengthParser| replace:: `Amplicon Length Parser`_
.. |flash_pass2| replace:: `Flash Pass 2`_
.. |merge_flash_stats| replace:: `Merge Flash Stats`_
.. |catenate| replace:: `Catenate`_
.. |uchime| replace:: `UCHIME Step`_
.. |merge_uchime_stats| replace:: `Merge UCHIME Stats`_
.. |otu_picking| replace:: `OTU Picking`_
.. |otu_rep_picking| replace:: `OTU Rep Picking`_
.. |otu_assigning| replace:: `OTU Assigning`_
.. |otu_table| replace:: `OTU Table`_
.. |otu_alignment| replace:: `OTU Alignment`_
.. |filter_alignment| replace:: `Filter Alignment`_
.. |phylogeny| replace:: `Phylogeny`_
.. |qiime_report| replace:: `QIIME Report`_
.. |multiple_rarefaction| replace:: `Multiple Rarefaction`_
.. |alpha_diversity| replace:: `Alpha Diversity`_
.. |collate_alpha| replace:: `Collate Alpha`_
.. |sample_rarefaction_plot| replace:: `Sample Rarefaction Plot`_
.. |qiime_report2| replace:: `QIIME Report 2`_
.. |single_rarefaction| replace:: `Single Rarefaction`_
.. |css_normalization| replace:: `CSS Normalization`_
.. |rarefaction_plot| replace:: `Rarefaction Plot`_
.. |summarize_taxa| replace:: `Summarize Taxonomy`_
.. |plot_taxa| replace:: `Plot Taxonomy`_
.. |plot_heatmap| replace:: `Plot Heatmap`_
.. |krona| replace:: `Krona`_
.. |plot_to_alpha| replace:: `Plot to Alpha`_
.. |beta_diversity| replace:: `Beta Diversity`_
.. |pcoa| replace:: `Principal Coordinate Analysis`_
.. |pcoa_plot| replace:: `PCoA Plot`_
.. |plot_to_beta| replace:: `Plot to Beta`_
.. |asva| replace:: `ASVA`_