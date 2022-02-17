[TOC]


EpiQC Pipeline
==============

EpiQC is a quality control pipeline for signal files (bigwig) generated from ChIP-Seq. The pipeline does a series of calculations on
these files to assess the quality of ChIp-Seq data. Four different metrics are computed from a single bigwig file.
First, BigWigInfo will be used to initial quality check on signal tracks
Second, ChromImpute imputes signal tracks for the given chromosome (currently only chr1 is supported) and using these imputed files, EpiQC computes 2 other metrics.
Third, SignalToNoise measurement by calculating the proportion of signal in top bins
And Fourth, the pipeline creates a heatmap from the correlation matrix obtained from EpiGeEC comparing only user samples.(Please note that currently comparing
large reference database with user samples is not supported).
Finally, the pipeline executes four consecutive report steps to create the final report of the pipeline with quality control labels.

This new pipeline can be used for pre-validation in order to assess the usability of a dataset in any given study, even
in the absence of the original raw reads files. This is an advantage, for instance, in the case of
human epigenomic datasets within IHEC, as signal tracks are made publicly available, while raw
data files are stored in controlled access repositories.

You can test this pipeline with ChIP-Seq samples from the IHEC portal :
https://epigenomesportal.ca/ihec/grid.html?assembly=4&build=2018-10


Usage
-----
```
#!text

usage: epiqc.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f] [--no-json]
                [--report] [--clean] [-l {debug,info,warning,error,critical}]
                [--sanity-check]
                [--container {wrapper, singularity} <IMAGE PATH>]
                [--genpipes_file GENPIPES_FILE] [-r READSETS] [-v]

Version: 4.1.2

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

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
  -j {pbs,batch,daemon,slurm}, --job-scheduler {pbs,batch,daemon,slurm}
                        job scheduler type (default: slurm)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
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
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a valid singularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![epiqc workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_epiqc.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_epiqc.png)
```
------
1- bigwiginfo
2- chromimpute
3- signal_to_noise
4- epigeec
5- bigwig_info_report
6- chromimpute_report
7- signal_to_noise_report
8- epigeec_report
9- epiqc_final_report

```
bigwiginfo
----------
Runs the tool bigWigInfo on bigwig files (

Inspecting signal tracks to identify some obvious problems that
could have an impact on the quality of the ChIP-Seq data is performed by UCSC-bigwiginfo
https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc-bigwiginfo/README.html)
bigWigInfo is capable of identifying obvious issues such as
missing chromosomes and insufficient track coverage, which are usually symptoms of
improperly generated tracks.

If the user has specified bigwig files in the readset file under BIGIWIG column, they will be utilized by
the tool. Otherwise, the user is required to process files using ChIp-Seq pipeline to generate
bigwig files. Then paths for bigwig files are reconstructed based on the ChIP-Seq readset file
and will be used subsequently.
(Note that: the readset file should be in the same folder as the ChIp-Seq output)

chromimpute
-----------
Runs the steps (bigwig_to_bedgraph, chromimpute_preprocess, Convert, ComputeGlobalDist,
GenerateTrainData, Train apply and eval) of the ChromImpute on the bigwig files given to the pipeline.
    All the output are stored in the imputation directory.


signal_to_noise
---------------
Binned signal resolution tracks from chromipute convert step is used to determined the
percentage of the whole file signal that was located in the top 5% and 10% of the bins
The resulted output is a tsv file

epigeec
-------
Runs the epigeec pipeline (https://bitbucket.org/labjacquespe/epigeec/src/master/)
on the bigwig files given to the pipeline
Epigeec pipeline is consisted of 3 sub-steps
1. Bigwig files are first converted to the hdf5 format
2. Filter or select provided regions as a bed file (include or exclude) [optional]. The user can specify
the options and the brf file path in the ini file. Otherwise, this step will be skipped
3. Finally the correlation matrix is computed


bigwig_info_report
------------------
This step is performed to generate report on bigwiginfo result

chromimpute_report
------------------
        This step is performed to generate a report comparing ChromImpute imputed signal
track and input signal track (in bedgraph format).

signal_to_noise_report
----------------------
This step is performed to generate report on signal_to_noise result

epigeec_report
--------------
This step is performed to generate a heatmap from EpiGeEC results

epiqc_final_report
------------------
Creates a report file for each bigwig file.

Alert levels :
    High Level Alert:
        Chromosome count is under 23
    Medium Level Alert:
        Whole genome bases covered under 25,000,000
        GeEC average correlation score under 50% within the same consortium tracks
        Signal in top 10% bins below 30%
        ChromImpute OBSERVED_1.0_IMPUTE_5.0 below 30%
    Low Level Alert:
        Whole genome bases covered under 75,000,000
        Signal in top 5% bins below 20%
        ChromImpute BOTH_1.0 below 20%


