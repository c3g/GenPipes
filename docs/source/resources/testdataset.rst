.. _docs_testdatasets:

GenPipes Test Datasets
======================

You can execute various GenPipes Pipelines using the following types of data:

* The real data which is generated from your genomic analysis instruments and then measured, sampled and read into various specified bioinformatics data formats.  

* Test datasets that are available in the absence of real genomic analysis data.

**Test Dataset** in the context of GenPipes refers to the dataset that needs to be analyzed by one of the GenPipes Pipelines. It can either be real data or sample data that is used to run the pipeline. Test dataset refers to some dataset that user can use to have hands-on on the pipeline. It is typically a smaller datasets (for e.g., one chromosome only and few samples for instance) so that the test runs of the pipelines using sample data get completed quickly say for demonstration purposes.

Test dataset is different from readset file which is input to the pipeline.  For other kinds of inputs required for GenPipes pipelines, see :ref:`here<gp_terminology>`.

In contrast to the test dataset, a :ref:`Readset File<docs_readset_file>` in the context of GenPipes actually describes the dataset (test dataset or real dataset) so that the pipeline can understand the type of data and process it.  Readset file is provided as input to almost all the GenPipes pipelines. Readset file contains information about the data to analyze; the path of the raw files, the type of sequencing, the name of the samples, etc.

.. note::  

    Please remember to use the correct dataset for the respective GenPipes pipelines.  The table below lists the test dataset download link for each of the GenPipes pipeline. Do not use the test dataset specified for a different pipeline.

+-----------------------------------------------+------------------------------------------+
|  *GenPipes Pipeline*                          |    *Test Dataset*                        |
+===============================================+==========================================+
| :ref:`HiC Pipeline<docs_gp_hicseq>`           |   `Download HiC Pipeline Dataset`_       |
+-----------------------------------------------+------------------------------------------+
| :ref:`Amplicon Seq<docs_gp_ampliconseq>`      |   `Download Amplicon Seq Dataset`_       |
+-----------------------------------------------+------------------------------------------+
| :ref:`ChIP Seq<docs_gp_chipseq>`              |   `Download ChIP Seq Dataset`_           |
+-----------------------------------------------+------------------------------------------+
| :ref:`DNA Seq<docs_gp_dnaseq>`                |   `Download DNA Seq Dataset`_            |
+-----------------------------------------------+------------------------------------------+
| :ref:`epiQC<docs_gp_epiqc>`                   |   `Download epiQC Dataset`_              |
+-----------------------------------------------+------------------------------------------+
| :ref:`Nanopore<docs_gp_nanopore>`             |   `Download Nanopore Dataset`_           |
+-----------------------------------------------+------------------------------------------+
| :ref:`Nanopore Covseq<docs_gp_nanopore_cov>`  |   `Download Nanopore Covseq Dataset`_    |
+-----------------------------------------------+------------------------------------------+
| :ref:`RNA Seq<docs_gp_rnaseq>`                |   `Download RNA Seq Dataset`_            |
+-----------------------------------------------+------------------------------------------+
| :ref:`Methyl Seq<docs_methylation>`           |   `Download Methyl Seq Dataset`_         |
+-----------------------------------------------+------------------------------------------+
| :ref:`PacBio Seq<docs_gp_dnaseq>`             |   `Download PacBio Seq Dataset`_         |
+-----------------------------------------------+------------------------------------------+
| :ref:`TumorPair Seq<docs_gp_tumorpair>`       |   `Download TumorPair Seq Dataset`_      |
+-----------------------------------------------+------------------------------------------+

.. include::  /resources/cov-seq-testdataset-note.inc

.. warning:: 

     PacBio Sequencing Pipeline is no longer available in GenPipes Release 3.2.0 and beyond.

----------------------------
Test Dataset Usage Examples
----------------------------

For various GenPipes pipelines, you can refer to usage examples and commands for issuing pipeline jobs using various options in the individual pipeline reference guide listed above or a short summary :ref:`here <doc_pipeline_usage_examples>`.

------------------------
Bioinformatic resources
------------------------

If you are looking for Bioinformatic resources such as available genomes with FASTA sequence, aligner indices and annotation files listed on `Bioinformatics resources <https://computationalgenomics.ca/cvmfs-genome/>`_ C3G website page, you can download those from the public repositories using scripts provided in `GenPipes Repository <https://bitbucket.org/mugqic/genpipes/src/master/resources/genomes/>`_.

You can also download the latest test datasets from Computational Genomics website `download page <https://datahub-90-cw3.p.genap.ca>`_.

.. Test dataset archive reference

.. _Download HiC Pipeline Dataset: https://datahub-90-cw3.p.genap.ca/hicseq.chr19.tar.gz
.. _Download Amplicon Seq Dataset: https://datahub-90-cw3.p.genap.ca/ampliconseq.tar.gz
.. _Download ChIP Seq Dataset:  https://datahub-90-cw3.p.genap.ca/chipseq.chr19.new.tar.gz
.. oldchipseq file commented out https://datahub-90-cw3.p.genap.ca/chipseq.chr19.tar.gz
.. _Download CoV Seq Dataset: https://datahub-90-cw3.p.genap.ca/covseq.chr19.tar.gz
.. _Download DNA Seq Dataset: https://datahub-90-cw3.p.genap.ca/dnaseq.chr19.tar.gz
.. _Download epiQC Dataset: https://datahub-90-cw3.p.genap.ca/epiQC.tar.gz
.. _Download Nanopore Dataset: https://datahub-90-cw3.p.genap.ca/nanopore.tar.gz
.. _Download Nanopore Covseq Dataset: https://datahub-90-cw3.p.genap.ca/nanopore_covseq.tar.gz
.. _Download RNA Seq Dataset: https://datahub-90-cw3.p.genap.ca/rnaseq.chr19.tar.gz
.. _Download RNA Seq Light Dataset: https://datahub-90-cw3.p.genap.ca/rnaseq_light.chr19.tar.gz
.. _Download Methyl Seq Dataset: https://datahub-90-cw3.p.genap.ca/methylseq.chr19.tar.gz
.. _Download PacBio Seq Dataset: https://datahub-90-cw3.p.genap.ca/pacbio.tar.gz
.. _Download TumorPair Seq Dataset: https://datahub-90-cw3.p.genap.ca/tumorPair.chr19.tar.gz
