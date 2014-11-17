MUGQIC Pipelines
================
This repository holds several bioinformatics pipelines developed at [McGill University and Génome Québec Innovation Centre](http://gqinnovationcenter.com) (MUGQIC).

MUGQIC pipelines consist of Python scripts which create a list of jobs running Bash commands. Those scripts support dependencies between jobs and smart restart mechanism if some jobs fail during pipeline execution. Jobs can be submitted in different ways: by being sent to a PBS scheduler like Torque or by being run as a series of commands in batch through a Bash script. Job commands and parameters can be modified through several configuration files.

On this page:

[TOC]


Software requirement
--------------------
MUGQIC pipelines have been tested with Python 2.7.


Quick setup for abacus, guillimin and mammouth users
----------------------------------------------------
Genomes and modules used by the pipelines are already installed on those clusters.
To access them, add the following lines to your *$HOME/.bash_profile*:
```
#!bash
umask 0002
    
## MUGQIC genomes and modules
    
HOST=`hostname`;
    
DNSDOMAIN=`dnsdomainname`;
    
if [[ $HOST == abacus* || $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then

  export MUGQIC_INSTALL_HOME=/sb/programs/analyste
  export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev

elif [[ $HOST == lg-* || $DNSDOMAIN == guillimin.clumeq.ca ]]; then

  export MUGQIC_INSTALL_HOME=/software/areas/genomics/phase2  
  export MUGQIC_INSTALL_HOME_DEV=/gs/project/mugqic/analyste_dev/phase2

elif [[ $BQMAMMOUTH == "mp2" ]]; then

  export MUGQIC_INSTALL_HOME=$(share_nobackup bourque)/mugqic_prod
  export MUGQIC_INSTALL_HOME_DEV=$(share_nobackup bourque)/mugqic_dev

fi
    
module use $MUGQIC_INSTALL_HOME/modulefiles $MUGQIC_INSTALL_HOME_DEV/modulefiles
```    

Also, set `JOB_MAIL` in your *$HOME/.bash_profile* to receive PBS job logs:
```
#!bash
export JOB_MAIL=<my.name@my.email.ca>
```

MUGQIC pipelines and compatible Python version are already installed as modules on those clusters.
To use them by default, add in your *$HOME/.bash_profile*:
```
#!bash
module load mugqic/python/2.7.8
module load mugqic/pipeline/<latest_version>
```
(find out the latest version with: "`module avail 2>&1 | grep mugqic/pipeline`").


### For guillimin and mammouth users
Set your `RAP_ID` (Resource Allocation Project ID from Compute Canada) in your *$HOME/.bash_profile*:
```
#!bash
export RAP_ID=<my-rap-id>
```


Usage
-----

For each pipeline, get help about usage, arguments and steps with:
```
#!bash
mugqic_pipeline/pipelines/<pipeline_name>/<pipeline_name>.py --help
```

Pipelines require as input one Readset File, one or more Configuration File(s) and possibly one design file, all described below.

For more information about a specific pipeline, visit:

* [DNA-Seq](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/dnaseq/)
* [RNA-Seq](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/rnaseq/)
* [RNA-Seq De Novo Assembly](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/rnaseq_denovo_assembly/)
* [PacBio Assembly](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/pacbio_assembly/)


Readset File
------------

The Readset File is a TAB-separated values plain text file with one line per readset and the following columns in any order:


### DNA-Seq, RNA-Seq, RNA-Seq De Novo Assembly

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Library: optional;
* RunType: "`PAIRED_END`" or "`SINGLE_END`"; mandatory;
* Run: optional;
* Lane: optional;
* QualityOffset: quality score offset integer used for trimming; optional;
* BED: relative or absolute path to BED file; optional;
* FASTQ1: relative or absolute path to first FASTQ file for paired-end readset or single FASTQ file for single-end readset; mandatory if BAM value is missing;
* FASTQ2: relative or absolute path to second FASTQ file for paired-end readset; mandatory if RunType value is "`PAIRED_END`";
* BAM: relative or absolute path to BAM file which will be converted into FASTQ files if they are not available; mandatory if FASTQ1 value is missing, ignored otherwise.

Example:

    Sample	Readset	Library	RunType	Run	Lane	QualityOffset	BED	FASTQ1	FASTQ2	BAM
    sampleA	readset1	lib0001	PAIRED_END	run100	1	33	path/to/file.bed	path/to/readset1.paired1.fastq.gz	path/to/readset1.paired2.fastq.gz	path/to/readset1.bam
    sampleA	readset2	lib0001	PAIRED_END	run100	2	33	path/to/file.bed	path/to/readset2.paired1.fastq.gz	path/to/readset2.paired2.fastq.gz	path/to/readset2.bam
    sampleB	readset3	lib0002	PAIRED_END	run200	5	33	path/to/file.bed	path/to/readset3.paired1.fastq.gz	path/to/readset3.paired2.fastq.gz	path/to/readset3.bam
    sampleB	readset4	lib0002	PAIRED_END	run200	6	33	path/to/file.bed	path/to/readset4.paired1.fastq.gz	path/to/readset4.paired2.fastq.gz	path/to/readset4.bam


### PacBio Assembly

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Smartcell: mandatory;
* NbBasePairs: total number of base pairs for this readset; mandatory;
* EstimatedGenomeSize: estimated genome size in number of base pairs used to compute seeding reads length cutoff; mandatory;
* BAS: comma-separated list of relative or absolute paths to BAS files (old PacBio format); mandatory if BAX value is missing, ignored otherwise;
* BAX: comma-separated list of relative or absolute paths to BAX files; BAX file list is used first if both BAX/BAS lists are present; mandatory if BAS value is missing;

Example:

    Sample	Readset	Smartcell	NbBasePairs	EstimatedGenomeSize	BAS	BAX
    sampleA	readset1	F_01_1	122169744	150000	path/to/readset1.bas.h5	path/to/readset1.1.bax.h5,path/to/readset1.2.bax.h5,path/to/readset1.3.bax.h5
    sampleA	readset2	F_01_2	105503472	150000	path/to/readset2.bas.h5	path/to/readset2.1.bax.h5,path/to/readset2.2.bax.h5,path/to/readset2.3.bax.h5
    sampleB	readset3	G_01_1	118603200	150000	path/to/readset3.bas.h5	path/to/readset3.1.bax.h5,path/to/readset3.2.bax.h5,path/to/readset3.3.bax.h5
    sampleB	readset4	G_01_2	104239488	150000	path/to/readset4.bas.h5	path/to/readset4.1.bax.h5,path/to/readset4.2.bax.h5,path/to/readset4.3.bax.h5


### For abacus users with Nanuq readsets
If your readsets belong to a [Nanuq](http://gqinnovationcenter.com/services/nanuq.aspx) project, use `nanuq2mugqic_pipeline.py` script in module `mugqic/tools` to automatically create a Readset File and symlinks to your readsets on abacus.


Configuration Files
-------------------
Pipeline cluster settings and command parameters can be customized using Configuration Files (`.ini` extension).
Those files have a structure similar to Microsoft Windows INI files e.g.:
```
#!ini
[my_section]
param1=my_param
dir1=my_dir/my_subdir
```

A parameter value is first searched in its specific section, then, if not found, in the special `DEFAULT` section.

Configuration files support interpolation. For example:
```
#!ini
[my_section]
dir1=my_dir
dir2=%(dir1)s/my_subdir
```
would resolve `dir2` value to `my_dir/my_subdir`.

Each pipeline has a default configuration file (`.base.ini` extension) set for running on abacus cluster using *Homo sapiens* reference genome:
```
#!bash
mugqic_pipeline/pipelines/<pipeline_name>/<pipeline_name>.base.ini
```

You can also pass a list of several configuration files to a pipeline command.
Files are read in the list order and each parameter value is overwritten if redefined in the next file.

This is useful to customize settings for a specific cluster or genome.
Each pipeline has a special configuration file for guillimin and mammouth clusters (`.guillimin.ini` and `.mammouth.ini` extensions respectively).
And various genome settings are available in `mugqic_pipeline/resources/genomes/config/`.

For example, to run the DNA-Seq pipeline on guillimin cluster with Mus musculus reference genome:
```
#!bash
mugqic_pipeline/pipelines/dnaseq/dnaseq.py -c mugqic_pipeline/pipelines/dnaseq/dnaseq.base.ini mugqic_pipeline/pipelines/dnaseq/dnaseq.guillimin.ini mugqic_pipeline/resources/genomes/config/Mus_musculus.GRCm38.ini ...
```


Design File
-----------
RNA-Seq and RNA-Seq De Novo Assembly pipelines can perform differential expression analysis if they are provided with an input design file.
The Design File is a TAB-separated values plain text file with one line per sample and the following columns:
* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* <contrast>: each of the following columns defines an experimental design contrast; the column name defines the contrast name, and the following values represent the sample group membership for this contrast:
    * '0' or '': the sample does not belong to any group;
    * '1': the sample belongs to the control group;
    * '2': the sample belongs to the test case group.

Example:

    Sample	Contrast1	Contrast2	Contrast3
    sampleA	1	1	1
    sampleB	2	0	1
    sampleC	0	2	0
    sampleD	0	0	2


Download and setup for external users
-------------------------------------


### Download

Visit our [Download page](https://bitbucket.org/mugqic/mugqic_pipeline/downloads) to get the latest stable release.

If you want to use the most recent development version:
```
#!bash
git clone git@bitbucket.org:mugqic/mugqic_pipeline.git
```


### Setup

MUGQIC Pipelines require genomes and modules resources to run properly.
First, set `MUGQIC_INSTALL_HOME` to the directory where you want to install those resources, in your *$HOME/.bash_profile*:
```
#!bash
## MUGQIC genomes and modules
    
export MUGQIC_INSTALL_HOME=/path/to/your/local/mugqic_resources
    
module use $MUGQIC_INSTALL_HOME/modulefiles
```


#### Genomes
Reference genomes and annotations must be installed in `$MUGQIC_INSTALL_HOME/genomes/`.
Default genome installation scripts are already available in `mugqic_pipeline/resources/genomes/`.
To install all of them, use the script `mugqic_pipeline/resources/genomes/install_all_genomes.sh`.


#### Modules
Software tools and associated modules must be installed in `$MUGQIC_INSTALL_HOME/software/` and `$MUGQIC_INSTALL_HOME/modulefiles/`.
Default software/module installation scripts are already available in `mugqic_pipeline/resources/modules/`.


Call home
---------
When pipeline jobs are submitted, a call home feature is invoked to collect some usage data. Those data are used to compute statistics and justify grant applications for pipeline funding support.

Data collected:

* Date and time
* Host and IP address
* Pipeline name
* Number of samples
* Pipeline steps


Contact us
----------
Please, ask questions or report bugs by sending us an email to [bioinformatics.service@mail.mcgill.ca](mailto:bioinformatics.service@mail.mcgill.ca).
