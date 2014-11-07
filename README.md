MUGQIC Pipelines
================
This repository holds several bioinformatics pipelines developed at [McGill University and Génome Québec Innovation Centre](http://gqinnovationcenter.com) (MUGQIC).

Visit our [wiki](https://biowiki.atlassian.net/wiki/display/PS/Pipeline+Space+Home) for an overview of the various pipelines.

MUGQIC pipelines consist of Python scripts which create a list of jobs running Bash commands. Those scripts support dependencies between jobs and smart restart mechanism if some jobs fail during pipeline execution. Jobs can be submitted in different ways: by creating a Bash script running a series of commands in batch or by sending those jobs to a PBS scheduler like Torque. Job commands and parameters can be modified through several configuration files.

Software requirement
--------------------
MUGQIC pipelines have been tested with Python 2.7.

Quick start for abacus, guillimin and mammouth users
----------------------------------------------------
Genomes and modules used by the pipelines are already installed on those clusters.
To access them, add the following lines to your *$HOME/.bash_profile*:

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
    

Also, set `JOB_MAIL` in your *$HOME/.bash_profile* to receive PBS job logs:

    export JOB_MAIL=my.name@email.ca


MUGQIC pipelines and compatible Python version are already installed as modules on those clusters.
To use them by default, add in your *$HOME/.bash_profile*:

    module load mugqic/python/2.7.8
    module load mugqic/pipeline/<latest_version>

(find out the latest version with: "`module avail 2>&1 | grep mugqic/pipeline`").


### For abacus users
To use parallel computing with some modules, add the following lines to your *$HOME/.bash_profile*:

    ## MPI
    export PATH=/sb/programs/mpi/mpi_pbs/openmpi-1.6/bin:$PATH
    export LD_LIBRARY_PATH=/sb/programs/mpi/mpi_pbs/openmpi-1.6/lib:$LD_LIBRARY_PATH

### For guillimin and mammouth users
Set your `RAP_ID` (Resource Allocation Project ID from Compute Canada) in your *$HOME/.bash_profile*:

    export RAP_ID=my-rap-id


Usage
-----

For each pipeline, get help about usage, arguments and steps with:

    mugqic_pipeline/pipelines/<pipeline_name>/<pipeline_name>.py --help

Pipelines require as input one Readset File and one or more Configuration File(s) described below.

For more information about a specific pipeline, visit:

* [DNA-Seq](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/dnaseq/)
* [RNA-Seq](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/rnaseq/)
* [RNA-Seq De Novo Assembly](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/rnaseq_denovo_assembly/)
* [PacBio Assembly](https://bitbucket.org/mugqic/mugqic_pipeline/src/python/pipelines/pacbio_assembly/)


Readset File
------------

The Readset File is a TAB-separated values plain text file with one line per readset and the following columns in any order:

### DNA-Seq, RNA-Seq, RNA-Seq De Novo Assembly

* SampleID: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
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

    SampleID	Readset	Library	RunType	Run	Lane	QualityOffset	BED	FASTQ1	FASTQ2	BAM
    sampleA	readset1	lib0001	PAIRED_END	run100	1	33	path/to/file.bed	path/to/readset1.paired1.fastq.gz	path/to/readset1.paired2.fastq.gz	path/to/readset1.bam
    sampleA	readset2	lib0001	PAIRED_END	run100	2	33	path/to/file.bed	path/to/readset2.paired1.fastq.gz	path/to/readset2.paired2.fastq.gz	path/to/readset2.bam
    sampleB	readset3	lib0002	PAIRED_END	run200	5	33	path/to/file.bed	path/to/readset3.paired1.fastq.gz	path/to/readset3.paired2.fastq.gz	path/to/readset3.bam
    sampleB	readset4	lib0002	PAIRED_END	run200	6	33	path/to/file.bed	path/to/readset4.paired1.fastq.gz	path/to/readset4.paired2.fastq.gz	path/to/readset4.bam


### PacBio Assembly

* SampleID: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Smartcell: mandatory;
* NbBasepairs: total number of base pairs for this readset; mandatory;
* EstimatedGenomeSize: estimated genome size in number of base pairs used to compute seeding reads length cutoff; mandatory;
* BAS: comma-separated list of relative or absolute paths to BAS files (old PacBio format); mandatory if BAX value is missing, ignored otherwise;
* BAX: comma-separated list of relative or absolute paths to BAX files; BAX file list is used first if both BAX/BAS lists are present; mandatory if BAS value is missing;

Example:

    SampleID	Readset	Smartcell	NbBasepairs	EstimatedGenomeSize	BAS	BAX
    sampleA	readset1	F_01_1	122169744	150000	path/to/readset1.bas.h5	path/to/readset1.1.bax.h5,path/to/readset1.2.bax.h5,path/to/readset1.3.bax.h5
    sampleA	readset2	F_01_2	105503472	150000	path/to/readset2.bas.h5	path/to/readset2.1.bax.h5,path/to/readset2.2.bax.h5,path/to/readset2.3.bax.h5
    sampleB	readset3	G_01_1	118603200	150000	path/to/readset3.bas.h5	path/to/readset3.1.bax.h5,path/to/readset3.2.bax.h5,path/to/readset3.3.bax.h5
    sampleB	readset4	G_01_2	104239488	150000	path/to/readset4.bas.h5	path/to/readset4.1.bax.h5,path/to/readset4.2.bax.h5,path/to/readset4.3.bax.h5


### For abacus users
If your readsets belong to a [NANUQ](http://gqinnovationcenter.com/services/nanuq.aspx) project, use `nanuq2mugqic_pipeline.py` script in module `mugqic/tools` to automatically create a Readset File and symlinks to your readsets on abacus.


Configuration (INI) Files
-------------------------
Is the standard configuration file for the pipeline. It's a plain text file composed of sections, keys and values. Section name appears on a line in square brackets ([default]). Every key has a name and a value, separated by an equals sign (=), i.e. clusterSubmitCmd=msub. Semicolons (;) or number signs (#) at the beginning of the line indicate a comment. Comment lines are ignored. Generally sections are associated to specific steps of the pipeline. If a property is associated to a specific step, the program will search for keys in the respective section. If the key is not found, the values in the default section will be used. Templates for Compute Canada's Guillimin and Mammouth clusters may already be available in the respective pipeline directory (\(root/pipeline/Pipeline_name\)).


Download
--------

To download the pipeline use git to obtain the most recent development version. Mugqic pipelines are hosted on github, and can be obtained via:

    git clone https://bitbucket.org/mugqic/mugqic_pipeline.git


If you encounter issues with the development version (master), or simply wish to obtain the most recent stable revision (1.1) then use:

    git fetch && git checkout 1.1

Whenever the latest snapshot from github is needed, use the command pull

    git pull origin

Setup
-----

In order to make pipelines that work anywhere, regardless of user path settings, and to control for versions of third party software and data associated to the analysis, the pipelines depend on environment modules. The modules created by the bioinformatics team in Compute Canada's Guillimin and Mammouth clusters are located in the $MUGCIC_INSTALL_HOME/modulefiles path and use the notation: mugqic/<software>/<version> . Click [here](https://biowiki.atlassian.net/wiki/display/CS/Software+and+Data+Dependencies) to setup modules on Guillimin and Mammouth clusters.



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
