[TOC]


epiQC Pipeline
=================

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

You can test this pipeline with ChIP-Seq samples from the [IHEC portal](https://epigenomesportal.ca/ihec/grid.html?assembly=4&build=2018-10).

**Special Note: you can use the same readset file used in ChIP-seq pipeline without any modification
but it should be in the same folder as the ChIp-Seq output. Because input files for
epiQC pipeline are located based on the readset file path).**


[Here](https://bitbucket.org/mugqic/genpipes/downloads//MUGQIC_Bioinfo_ChIP-Seq.pptx)
is more information about developing epiQC pipeline that you may find interesting.


Usage
-----
```
#!text

usage: epiqc.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]                        
                  [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]                   
                  [--no-json] [--report] [--clean]
                  [-l {debug,info,warning,error,critical}] [--sanity-check]
                  [--container {wrapper, singularity} <IMAGE PATH>] [-r READSETS] [-v]

Version: 3.6.0                                                                    

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
                        Type of pipeline (default chipseq)
  -r READSETS, --readsets READSETS                                                         
                        readset file                                                       
  -v, --version         show the version information and exit

```
![chipseq workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_epiqc.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_epiqc.png)

```
epiqc:                                                                                   
1 - bigwiginfo                                                                     
2 - chromimpute                                                                             
3 - signal_to_noise                                                                 
4 - epigeec                                                                
5 - bigwig_info_report                                                                
6 - chromimpute_report                                                                
7 - signal_to_noise_report 
8 - epigeec_report
9 - epiqc_final_report                                                                   
```


bigwiginfo
-------------------
Runs the tool bigWigInfo on bigwig files (
            Inspecting signal tracks to identify some obvious problems that
            could have an impact on the quality of the ChIP-Seq data is performed by [UCSC-bigWigInfo](
            https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc-bigwiginfo/README.html).
            bigWigInfo is capable of identifying obvious issues such as
            missing chromosomes and insufficient track coverage, which are usually symptoms of
            improperly generated tracks.

If the user has specified bigwig files in the readset file under BIGIWIG column, they will be utilized by
            the tool. Otherwise, the user is required to process files using ChIp-Seq pipeline to generate
            bigwig files. Then paths for bigwig files are reconstructed based on the ChIP-Seq readset file
            and will be used subsequently.
            (Note that: the readset file should be in the same folder as the ChIp-Seq output).

    Input Files : bigwig files

chromimpute
-----------
[ChromImpute](https://ernstlab.biolchem.ucla.edu/ChromImpute/) is a Java software for large-scale systematic epigenome imputation.
ChromImpute takes an existing compendium of signal tracks (bedgraph files) 
uses it to predict signal tracks for mark-sample combinations not experimentally 
mapped or to generate a potentially more robust version of data sets that 
have been mapped experimentally. 
ChromImpute bases its predictions on features from signal tracks 
of other marks that have been mapped in the target sample and the target mark in 
other samples with these features combined using an ensemble of regression trees.
Currently, only signal tracks mapped to GRCh38 are supported with chromimpute 
analysis. GRCh37 will be added in the future. For better results, usage of 
multiple histone marks from one sample is recommended.

Chromimpute step is composed of several sub steps that are invoked internally. 
The sub steps are as follows.

### bigwig_to_bedgraph
[ucsc-bigwigtobedgraph](https://bioconda-recipes-demo.readthedocs.io/en/docs/recipes/ucsc-bigwigtobedgraph/README.html)
is used to convert bigwig files to bedgraph files (Used in ChromImpute subsequently).

    Input Files : bigwig files

### chromimpute_preprocess

This step is performed to create chromimpute directories, chromosome sizes file, inputinfo file with IHEC
            and user samples and finally link the converted IHEC bedgraph files to user directory. In order to run 
            the chromimpute, inputinfo and chromsizes file are required to be in the imputation directory.
            Please note: chromsizes and inputinfo files are created dynamically when the user runs the pipeline. It is 
            not necessary to submit jobs to create them.
Currently, only GRCh38 chr1 is used for chromimpute analysis but, it will be sufficient
for accurately evaluate the signal quality. Chromosome name can be specified in
the chromimpute_preprocess section of the ini file.

### chromimpute_convert
This step is performed to convert each unique mark and sample combination signal tracks (bedgraph files)
            in the user dataset into binned (25 bp) signal resolution tracks. Since 
the stored converted files are binned using 25bps, changing the resolution value specified in the ini
file is not recommended for chromimpute analysis.

NOTE: If you got index out of bound exception, check whether your reference genome version of the bedgraph file
            is similar to chromosome_sizes_file inside the imputation folder.

     This step takes as input files:

     1. bigwig files from bigwig_to_bedgraph
     2. Chromosome sizes file
     3. inputinfo file

### chromimpute_compute_global_dist

This steps is used to Compute the global distance based on correlation for each mark in each sample with
            the same mark in all other samples in inputinfo file. Creates a file for each mark in each sample containing
            a ranked list of the globally nearest samples.

     This step takes as input files:

     1. Converted signal tracks from chromimpute_convert
     2. Converted and linked signal tracks from IHEC data portal
     3. Chromosome sizes file
     4. inputinfo file

### chromimpute_generate_train_data

This step is performed to generate a set of training data instances taking directory of converted data 
            and global distances.
     
     This step takes as input files:

     1. Converted signal tracks from chromimpute_convert
     2. Converted and linked signal tracks from IHEC data portal
     3. Chromosome sizes file
     4. inputinfo file
     5. global distance files from chromimpute_compute_global_dist

### chromimpute_train

This step is used to train regression trees based on the feature data produced by GenerateTrainData.

     This step takes as input files:

     1. Predictor data produced by GenerateTrainData
     2. Chromosome sizes file
     3. inputinfo file

### chromimpute_apply

This step is used to apply the predictors generated in the Train command to generate the imputed data.
            A job is created for each sample and mark given in the dataset.

     This step takes as input files:

     1. Converted signal tracks from chromimpute_convert
     2. Converted and linked signal tracks from IHEC data portal
     3. Chromosome sizes file
     4. inputinfo file
     5. global distance files from chromimpute_compute_global_dist
     6. Predictor data produced by GenerateTrainData

### chromimpute_eval

The final step of the chromimpute pipeline is used to compare the agreement between an observed and imputed data set.
            A job is created for every sample-mark given in the dataset. percent1 percent2 Gives lower and upper 
percentages to use in evaluation. Default values for percent1 is 1% and percent2 is 5%

     This step takes as input files:

     1. Converted signal tracks from chromimpute_convert
     2. imputed signal tracks from chromimpute_apply

signal_to_noise
--------------------
Binned signal resolution tracks from chromipute_convert step is used to determine the
            percentage of the whole file signal that was located in the top percent1 and percent2 of the bins. 
The default percentages are 5% and 10% and the resulted output is a tsv file.


     This step takes as input files:

     1. Converted signal tracks from chromimpute_convert


epigeec
--------------------

This step is performed to run the [EpiGeEC pipeline](https://bitbucket.org/labjacquespe/epigeec/src/master/).
Epigeec pipeline is consist of 3 sub-steps
            1. Bigwig files are first converted to the hdf5 format
            2. Filter or select provided regions as a bed file (include or exclude) [optional]. 
The user can specify the options and the bed file path in the ini file. Otherwise, this step will be skipped
            3. Finally, the correlation matrix is computed.


epiqc_report
--------------------

epiQC report is consist of five sub-steps including four independent steps to generate individual
report files for each metric computed and one step to combine all of them to generate
a final report file. A user can independently run each sub step but in order to
generate final report, reports from BigWigInfo, ChromImpute and EpiGeEC are required.
This final report is a TSV file with the decision of each signal track
[i.e whether the sample has passed from all the
metrics or there are alerts that user need to be concerned]

1. ### bigwiginfo_report

This step is performed to generate report on bigwiginfo result


2. ### chromimpute_report

This step is performed to generate a report comparing ChromImpute imputed signal
track and input signal track (in bedgraph format).


3. ### signal_to_noise_report

This step is performed to generate report on signal_to_noise result

4. ### epigeec_report

This step is performed to generate a heatmap from EpiGeEC results

5. ### epiqc_final_report

Once all metrics have been obtained, they are gathered in a TSV formatted
report and compared to predetermined thresholds. A column holding EpicsQC’s verdict for a signal track is
also included, using these thresholds evaluations. The following metrics threshold were used:
Note: an Alert will be generated even if there is an issue with one metric.

* **High Level Alert:**
    * Chromosome count is under 23
* **Medium Level Alert:**
  * Whole genome bases covered under 25,000,000
  * Signal in top 10% bins below 30% 
  * ChromImpute OBSERVED_1.0_IMPUTE_5.0 below 30% 
    
* **Low Level Alert:**
  * Whole genome bases covered under 75,000,000 
  * Signal in top 5% bins below 20% 
  * ChromImpute BOTH_1.0 below 20%
    
Note: Currently EpiGeEC metric is not used for the final decision.

