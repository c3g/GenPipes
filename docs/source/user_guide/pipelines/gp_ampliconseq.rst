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
      subsample
      config
      param
      chao
      shannon      

Amplicon Sequencing Pipeline
============================

:bdg-primary:`Version` |genpipes_version|

.. tab-set:: 

      .. tab-item:: About

         .. card::

            Amplicon sequencing (ribosomal RNA gene amplification analysis) is a metagenomic pipeline. It is based on the established `Quantitative Insights into Microbial Ecology <https://www.ncbi.nlm.nih.gov/pubmed/22161565>`_ (QIIME) procedure for amplicon-based metagenomics. It assembles read pairs using `Fast Length Adjustment of Short Reads <https://www.ncbi.nlm.nih.gov/pubmed/21903629>`_ (FLASH), detects chimeras with `UCHIME <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3150044/>`_, and picks operational taxonomic units using `VSEARCH <https://www.ncbi.nlm.nih.gov/pubmed/27781170>`_. Operational taxonomic units are then aligned using `PyNAST <https://www.ncbi.nlm.nih.gov/pubmed/19914921>`_ and clustered with `FastTree <https://www.ncbi.nlm.nih.gov/pubmed/19377059>`_. Standard diversity indices, taxonomical assignments, and ordinations are then calculated and reported graphically.

            Besides `QIIME Pipeline`_ GenPipes also supports the `DADA2`_ protocol, a relatively new processing workflow for recovering single-nucleotide resolved Amplicon Sequence Variants (ASVs) from amplicon data.

            See :ref:`ampschema` tab for the pipeline workflow.

      .. tab-item:: Usage

         .. dropdown:: Command

            .. code::

               ampliconseq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                              [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                              [--no-json] [--report] [--clean]
                              [-l {debug,info,warning,error,critical}]
                              [--sanity-check]
                              [--container {wrapper, singularity} <IMAGE PATH>]
                              [--genpipes_file GENPIPES_FILE]
                              [-t {qiime, dada2 (default) }] [-d DESIGN] [-r READSETS] [-v]

         .. dropdown:: Example Run

            Use the following commands to execute Amplicon sequence pipeline:

            .. include::  /user_guide/pipelines/example_runs/ampliconseq.inc

            .. include:: /user_guide/pipelines/notes/scriptfile_deprecation.inc

            .. card:: Test Dataset
               :link: docs_testdatasets
               :link-type: ref

               You can download the test dataset for this pipeline :ref:`here<docs_testdatasets>`.

         .. dropdown:: Options

            .. include:: opt_ampliconseq.inc
            .. include:: /common/gp_design_opt.inc 
            .. include:: /common/gp_readset_opt.inc 
            .. include:: /common/gp_common_opt.inc 

      .. tab-item:: Schema
         :name: ampschema        

         .. dropdown:: QIIME  

            Figure below shows the schema of the amplicon sequencing protocol - `QIIME Pipeline`_.

            .. figure:: /img/pipelines/mmd/ampliconseq.qiime.mmd.png
               :align: center
               :alt: QIIME ampseq
               :width: 100%
               :figwidth: 95% 

               Figure: Schema of QIIME Amplicon Sequencing protocol

            .. figure:: /img/pipelines/mmd/legend.mmd.png
               :align: center
               :alt: dada2 ampseq 
               :width: 100%
               :figwidth: 75%

            `Click for a high resolution image of Amplicon Sequencing Pipeline (qiime) <https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/mmd/png/ampliconseq.qiime.mmd.png>`_.

         .. dropdown:: DADA2

            .. figure:: /img/pipelines/mmd/ampliconseq.dada2.mmd.png
               :align: center
               :alt: dada2 ampseq 
               :width: 70%
               :figwidth: 95%

               Figure: Schema of DADA2 Amplicon Sequencing protocol

            .. figure:: /img/pipelines/mmd/legend.mmd.png
               :align: center
               :alt: dada2 ampseq 
               :width: 100%
               :figwidth: 75%

            `Click for a high resolution image of Amplicon Sequencing Pipeline (dada2) <https://bitbucket.org/mugqic/genpipes/src/master/resources/workflows/mmd/png/ampliconseq.dada2.mmd.png>`_.

      .. tab-item:: Steps

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
         | 14.| |otu_alignment|                |                                 |
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
         | 32.| |pcoa|                         |                                 |
         +----+--------------------------------+                                 |
         | 33.| |pcoa_plot|                    |                                 |
         +----+--------------------------------+                                 |
         | 34.| |plot_to_beta|                 |                                 |
         +----+--------------------------------+---------------------------------+

         .. card::

            .. include:: steps_ampseq.inc

      .. tab-item:: Details

         .. card::

            Amplicon sequencing is a highly targeted gene sequencing approach used to analyze genetic variation in specific genomic regions. Amplicons are Polymerase Chain Reaction (PCR) products and the ultra-deep sequencing allows for efficient variant identification and characterization. Amplicon sequencing uses oligonucleotide probes that target and capture genomic regions of interest and then uses next-generation sequencing techniques. 

            **Uses of Amplicon sequencing**

            #. Diagnostic microbiology utilizes amplicon-based profiling that allows to sequence selected amplicons such as regions encoding 16S rRNA that are used for species identification. 

            #. Discovery of rare somatic mutations in complex samples such as tumors mixed with germline DNA.

            See :ref:`More Information <More Information Ampliconseq>` section below for details. 

.. _More Information Ampliconseq:

More information
-----------------

For the latest implementation and usage details refer to Amplicon Sequencing implementation `README.md <https://bitbucket.org/mugqic/genpipes/src/master/pipelines/ampliconseq/README.md>`_ file.

* `Amplicon sequencing techniques <https://sapac.illumina.com/techniques/sequencing/dna-sequencing/targeted-resequencing/amplicon-sequencing.html>`_

* `Amplicon Sequencing Primer <http://apc.ucc.ie/pdf_old/Amplicon%20Sequencing.pdf>`_

* `High-throughput amplicon sequencing <https://www.biorxiv.org/content/10.1101/392332v2>`_.

* `Trimmomatic - flexible trimming <https://academic.oup.com/bioinformatics/article/30/15/2114/2390096>`_.

.. The following are html links used in this text

.. _QIIME Pipeline: http://qiime.org
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
