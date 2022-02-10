GenPipes
================
This repository holds several bioinformatics pipelines developed at [McGill University and Génome Québec Innovation Centre](http://gqinnovationcenter.com) (MUGQIC), as part of the [GenAP project](https://genap.ca).

GenPipes consists of Python scripts which create a list of jobs running Bash commands. Those scripts support dependencies between jobs and smart restart mechanism if some jobs fail during pipeline execution. Jobs can be submitted in different ways: by being sent to a PBS scheduler like Torque or by being run as a series of commands in batch through a Bash script. Job commands and parameters can be modified through several configuration files.

**For a more detailed tutorial on how to use GenPipes, please visit our [tutorial page](http://www.computationalgenomics.ca/tutorials/).**  


On this page:

[TOC]


Software requirement
--------------------
GenPipes have been tested with Python 2.7.


Quick setup for abacus, Beluga, graham, cedar and mammouth users
----------------------------------------------------
Genomes and modules used by the pipelines are already installed on a CVMFS partition mounted on all those clusters in `/cvmfs/soft.mugqic/CentOS6`.
To access them, add the following lines to your *$HOME/.bash_profile*:

```
#!bash
umask 0002

## MUGQIC genomes and modules

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6

module use $MUGQIC_INSTALL_HOME/modulefiles
```

For MUGQIC analysts, add the following lines to your *$HOME/.bash_profile*:

```
#!bash
umask 0002

## MUGQIC genomes and modules for MUGQIC analysts

HOST=`hostname`;

DNSDOMAIN=`dnsdomainname`;

export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6

if [[ $HOST == abacus* || $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then

  export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev

elif [[ $HOST == ip* || $DNSDOMAIN == m  ]]; then

  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev

elif [[ $HOST == cedar* || $DNSDOMAIN == cedar.computecanada.ca ]]; then

  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev


elif [[ $HOST == beluga* || $DNSDOMAIN == beluga.computecanada.ca ]]; then

  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev

fi

module use $MUGQIC_INSTALL_HOME/modulefiles $MUGQIC_INSTALL_HOME_DEV/modulefiles
export RAP_ID=<my-rap-id>
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
module load mugqic/python/3.9.1
module load mugqic/genpipes/<latest_version>
```
(find out the latest version with: "`module avail 2>&1 | grep mugqic/genpipes`").


### For Compute Canada users
Set your `RAP_ID` (Resource Allocation Project ID from Compute Canada) in your *$HOME/.bash_profile*:
```
#!bash
export RAP_ID=<my-rap-id>
```


Download and setup for external users
-------------------------------------


### Download

Visit our [Download page](https://bitbucket.org/mugqic/genpipes/downloads) to get the latest stable release.

If you want to use the most recent development version:
```
#!bash
git clone git@bitbucket.org:mugqic/genpipes.git
```

#### GenPipes' Container:

A new installation with a better taste:

Singularity needs to be [installed on your system](https://github.com/hpcng/singularity/blob/master/INSTALL.md)

Then, make sure that you have fuse installed on your system,  if `ls /dev/fuse` returns no error, you are all set.

Once the genpipes repo has been cloned, run the following command to install the container and wrapper code for the fuse libraries.

```
#!bash
./genpipes/resources/container/get_wrapper.sh
```

You can access the Genpipes container by typing:

```
#!bash
./genpipes/resources/container/bin/container_wrapper.sh`

```

You can get more information to run [Genpipes with containers here](https://github.com/c3g/genpipes_in_a_container)

### Setup

Set `MUGQIC_PIPELINES_HOME` to your local copy path, in your *$HOME/.bash_profile*:
```
#!bash
export MUGQIC_PIPELINES_HOME=/path/to/your/local/genpipes
```

GenPipes (formerly called MUGQIC Pipelines) require genomes and modules resources to run properly.
First, set `MUGQIC_INSTALL_HOME` to the directory where you want to install those resources, in your *$HOME/.bash_profile*:
```
#!bash
## MUGQIC genomes and modules

export MUGQIC_INSTALL_HOME=/path/to/your/local/mugqic_resources

module use $MUGQIC_INSTALL_HOME/modulefiles
```


#### Genomes
Reference genomes and annotations must be installed in `$MUGQIC_INSTALL_HOME/genomes/`.
Default genome installation scripts are already available in `$MUGQIC_PIPELINES_HOME/resources/genomes/`.
To install all of them at once, use the script `$MUGQIC_PIPELINES_HOME/resources/genomes/install_all_genomes.sh`.

All species-related files are in:
`$MUGQIC_INSTALL_HOME/genomes/species/<species_scientific_name>.<assembly>/`
e.g. for *Homo sapiens* assembly *GRCh37*, the directory has the following (incomplete) hierarchy:
```
#!text
$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/
├── annotations/
│   ├── gtf_tophat_index/
│   ├── Homo_sapiens.GRCh37.dbSNP142.vcf.gz
│   ├── Homo_sapiens.GRCh37.dbSNP142.vcf.gz.tbi
│   ├── Homo_sapiens.GRCh37.Ensembl75.geneid2Symbol.tsv
│   ├── Homo_sapiens.GRCh37.Ensembl75.genes.length.tsv
│   ├── Homo_sapiens.GRCh37.Ensembl75.genes.tsv
│   ├── Homo_sapiens.GRCh37.Ensembl75.GO.tsv
│   ├── Homo_sapiens.GRCh37.Ensembl75.gtf
│   ├── Homo_sapiens.GRCh37.Ensembl75.ncrna.fa
│   ├── Homo_sapiens.GRCh37.Ensembl75.rrna.fa
│   ├── Homo_sapiens.GRCh37.Ensembl75.transcript_id.gtf
│   ├── Homo_sapiens.GRCh37.Ensembl75.vcf.gz
│   ├── ncrna_bwa_index/
│   └── rrna_bwa_index/
├── downloads/
│   ├── ftp.1000genomes.ebi.ac.uk/
│   ├── ftp.ensembl.org/
│   └── ftp.ncbi.nih.gov/
├── genome/
│   ├── bowtie2_index/
│   ├── bwa_index/
│   ├── Homo_sapiens.GRCh37.dict
│   ├── Homo_sapiens.GRCh37.fa
│   ├── Homo_sapiens.GRCh37.fa.fai
│   └── star_index/
├── Homo_sapiens.GRCh37.ini
└── log/
```
The assembly name is the one used by the download source e.g. "*GRCh37*" for [Ensembl](http://www.ensembl.org/).
Each species directory contains a `<scientific_name>.<assembly>.ini` file
which lists among other things, the assembly synonyms e.g. "*hg19*":

`Homo_sapiens.GRCh37.ini`
```
#!ini
[DEFAULT]
scientific_name=Homo_sapiens
common_name=Human
assembly=GRCh37
assembly_synonyms=hg19
source=Ensembl
version=75
dbsnp_version=142
```

##### Install a new Genome

New genomes and annotations can be installed semi-automatically from [Ensembl](http://www.ensembl.org/) (vertebrate species),
[EnsemblGenomes](http://ensemblgenomes.org/) (other species) or [UCSC](http://genome.ucsc.edu/) (genome and indexes only; no annotations).

Example for Chimpanzee:

* Retrieve the species scientific name on [Ensembl](http://useast.ensembl.org/Pan_troglodytes/Info/Index?redirect=no) or [UCSC](http://genome.ucsc.edu/cgi-bin/hgGateway): "*Pan troglodytes*"

* Retrieve the assembly name:
    - Ensembl: "*CHIMP2.1.4*"
    - UCSC: "*panTro4*"

* Retrieve the source version:
    - Ensembl: "78"
    - UCSC: unfortunately, UCSC does not have version numbers. Use [panTro4.2bit](http://hgdownload.soe.ucsc.edu/goldenPath/panTro4/bigZips/) date formatted as "YYYY-MM-DD": "2012-01-09"

* `cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh` e.g.:

    - Ensembl:

            cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/Pan_troglodytes.CHIMP2.1.4.sh

    - UCSC:

            cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/Pan_troglodytes.panTro4.sh

* Modify `$MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh` (`ASSEMBLY_SYNONYMS` can be left empty but if you know that 2 assemblies
are identical apart from `chr` sequence prefixes, document it):

    - Ensembl:

            SPECIES=Pan_troglodytes   # With "_"; no space!
            COMMON_NAME=Chimpanzee
            ASSEMBLY=CHIMP2.1.4
            ASSEMBLY_SYNONYMS=panTro4
            SOURCE=Ensembl
            VERSION=78

    - UCSC:

            SPECIES=Pan_troglodytes   # With "_"; no space!
            COMMON_NAME=Chimpanzee
            ASSEMBLY=panTro4
            ASSEMBLY_SYNONYMS=CHIMP2.1.4
            SOURCE=UCSC
            VERSION=2012-01-09

* Running `bash $MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh` will install the genome in $MUGQIC_INSTALL_HOME_DEV (by default). This will download and install genomes, indexes and, for Ensembl only, annotations (GTF, VCF, etc.).
[ADMINS ONLY] To install it in $MUGQIC_INSTALL_HOME run `bash $MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh MUGQIC_INSTALL_HOME`.

    If the genome is big, separate batch jobs will be submitted to the cluster for bwa, bowtie/tophat, star indexing.
    Check that jobs are completed OK.

* [ADMINS ONLY] If the new genome has been installed in `$MUGQIC_INSTALL_HOME_DEV`, to deploy in `$MUGQIC_INSTALL_HOME`:

        rsync -vca --no-o --no-g --no-p --size-only -I -O --ignore-times $MUGQIC_INSTALL_HOME_DEV/genomes/species/<scientific_name>.<assembly> $MUGQIC_INSTALL_HOME/genomes/species/

* Add the newly created INI file to the genome config files for further usage in pipeline command:

        cp $MUGQIC_INSTALL_HOME/genomes/species/<scientific_name>.<assembly>/<scientific_name>.<assembly>.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/


#### Modules
Software tools and associated modules must be installed in `$MUGQIC_INSTALL_HOME/software/` and `$MUGQIC_INSTALL_HOME/modulefiles/`.
Default software/module installation scripts are already available in `$MUGQIC_PIPELINES_HOME/resources/modules/`.

##### Install a new Module

New software tools and associated modules can be installed semi-automatically:

* `cp $MUGQIC_PIPELINES_HOME/resources/modules/MODULE_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh`

* Modify `$MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh` following the instructions inside.

* Run `$MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh` with no arguments. By default, it will download and extract the remote software archive, build the software and create the associated module, all in `$MUGQIC_INSTALL_HOME_DEV` if it is set.

* If everything is OK, to install it in production, run:

        $MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh MUGQIC_INSTALL_HOME
    (no `$` before `MUGQIC_INSTALL_HOME`!).

* Check if the module is available with: `module avail 2>&1 | grep mugqic/<my_software>/<version>`

Usage
-----

For each pipeline, get help about usage, arguments and steps with:

* if you use a `mugqic/genpipes/<version>` module on our clusters (or `mugqic/mugqic_pipelines/<version>`), simply:
```
#!bash
<pipeline_name>.py --help
```
* if you use your own local install:
```
#!bash
$MUGQIC_PIPELINES_HOME/pipelines/<pipeline_name>/<pipeline_name>.py --help
```

Pipelines require as input one Readset File, one or more Configuration File(s) and possibly one Design File, all described below.

For more information about a specific pipeline, visit:

### [DNA-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/dnaseq/)
### [DNA-Seq high Coverage Pipeline Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/dnaseq_high_coverage/)
### [RNA-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/rnaseq/)
### [RNA-Seq De Novo Assembly Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/rnaseq_denovo_assembly/)
### [RNA-Seq Light Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/rnaseq_light/)
### [PacBio Assembly Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/pacbio_assembly/)
### [ChIP-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/chipseq/)
### [Amplicon-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/ampliconseq/)
### [Methyl-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/methylseq/)
### [Nanopore Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/nanopore/)
### [Tumor Pair Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/tumor_pair/)
### [CoV-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/covseq/)
### [Nanopore CoV-Seq Pipeline](https://bitbucket.org/mugqic/genpipes/src/master/pipelines/nanopore_covseq/)

Readset File
------------

The Readset File is a TAB-separated values plain text file with one line per readset and the following columns in any order:


### DNA-Seq, DNA-Seq high Coverage, RNA-Seq, RNA-Seq De Novo Assembly, Amplicon-Seq, Tumor Pair, Methyl-Seq, CoV-Seq

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Library: optional;
* RunType: `PAIRED_END` or `SINGLE_END`; mandatory;
* Run: optional;
* Lane: optional;
* Adapter1 : sequence of the forward trimming adapter
* Adapter2 : sequence of the reverse trimming adapter
* QualityOffset: quality score offset integer used for trimming; optional;
* BED: relative or absolute path to BED file; optional;
* FASTQ1: relative or absolute path to first FASTQ file for paired-end readset or single FASTQ file for single-end readset; mandatory if BAM value is missing;
* FASTQ2: relative or absolute path to second FASTQ file for paired-end readset; mandatory if RunType value is "`PAIRED_END`";
* BAM: relative or absolute path to BAM file which will be converted into FASTQ files if they are not available; mandatory if FASTQ1 value is missing, ignored otherwise.

Example:

    Sample	Readset	Library	RunType	Run	Lane	Adapter1	Adapter2	QualityOffset	BED	FASTQ1	FASTQ2	BAM
    sampleA	readset1	lib0001	PAIRED_END	run100	1	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	33	path/to/file.bed	path/to/readset1.paired1.fastq.gz	path/to/readset1.paired2.fastq.gz	path/to/readset1.bam
    sampleA	readset2	lib0001	PAIRED_END	run100	2	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	33	path/to/file.bed	path/to/readset2.paired1.fastq.gz	path/to/readset2.paired2.fastq.gz	path/to/readset2.bam
    sampleB	readset3	lib0002	PAIRED_END	run200	5	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	33	path/to/file.bed	path/to/readset3.paired1.fastq.gz	path/to/readset3.paired2.fastq.gz	path/to/readset3.bam
    sampleB	readset4	lib0002	PAIRED_END	run200	6	AGATCGGAAGAGCACACGTCTGAACTCCAGTCA	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT	33	path/to/file.bed	path/to/readset4.paired1.fastq.gz	path/to/readset4.paired2.fastq.gz	path/to/readset4.bam

### ChIP-Seq

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* MarkName: name of the histone mark; mandatory
* MarkType: type of mark for MACS2 calling must be either B (Broad), N (Narrow) or I (Input); mandatory
* Library: optional;
* RunType: `PAIRED_END` or `SINGLE_END`; mandatory;
* Run: optional;
* Lane: optional;
* Adapter1 : sequence of the forward trimming adapter
* Adapter2 : sequence of the reverse trimming adapter
* QualityOffset: quality score offset integer used for trimming; optional;
* BED: relative or absolute path to BED file; optional;
* FASTQ1: relative or absolute path to first FASTQ file for paired-end readset or single FASTQ file for single-end readset; mandatory if BAM value is missing;
* FASTQ2: relative or absolute path to second FASTQ file for paired-end readset; mandatory if RunType value is "`PAIRED_END`";
* BAM: relative or absolute path to BAM file which will be converted into FASTQ files if they are not available; mandatory if FASTQ1 value is missing, ignored otherwise.

Example:

    Sample  Readset  MarkName MarkType Library RunType Run Lane    Adapter1    Adapter2    QualityOffset   BED FASTQ1  FASTQ2  BAM
    sampleA readset1 H3K27ac  N        lib0001 PAIRED_END  run100  1   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   33  path/to/file.bed    path/to/readset1.paired1.fastq.gz   path/to/readset1.paired2.fastq.gz   path/to/readset1.bam
    sampleA readset2 H3K27ac  N        lib0001 PAIRED_END  run100  2   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   33  path/to/file.bed    path/to/readset2.paired1.fastq.gz   path/to/readset2.paired2.fastq.gz   path/to/readset2.bam
    sampleB readset3 Input    I        lib0002 PAIRED_END  run200  5   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   33  path/to/file.bed    path/to/readset3.paired1.fastq.gz   path/to/readset3.paired2.fastq.gz   path/to/readset3.bam
    sampleB readset4 Input    I        lib0002 PAIRED_END  run200  6   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   33  path/to/file.bed    path/to/readset4.paired1.fastq.gz   path/to/readset4.paired2.fastq.gz   path/to/readset4.bam

### PacBio Assembly

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; mandatory;
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Smartcell: mandatory;
* NbBasePairs: total number of base pairs for this readset; mandatory;
* EstimatedGenomeSize: estimated genome size in number of base pairs used to compute seeding read length cutoff; mandatory;
* BAS: comma-separated list of relative or absolute paths to BAS files (old PacBio format); mandatory if BAX value is missing, ignored otherwise;
* BAX: comma-separated list of relative or absolute paths to BAX files; BAX file list is used first if both BAX/BAS lists are present; mandatory if BAS value is missing.

Example:

    Sample	Readset	Smartcell	NbBasePairs	EstimatedGenomeSize	BAS	BAX
    sampleA	readset1	F_01_1	122169744	150000	path/to/readset1.bas.h5	path/to/readset1.1.bax.h5,path/to/readset1.2.bax.h5,path/to/readset1.3.bax.h5
    sampleA	readset2	F_01_2	105503472	150000	path/to/readset2.bas.h5	path/to/readset2.1.bax.h5,path/to/readset2.2.bax.h5,path/to/readset2.3.bax.h5
    sampleB	readset3	G_01_1	118603200	150000	path/to/readset3.bas.h5	path/to/readset3.1.bax.h5,path/to/readset3.2.bax.h5,path/to/readset3.3.bax.h5
    sampleB	readset4	G_01_2	104239488	150000	path/to/readset4.bas.h5	path/to/readset4.1.bax.h5,path/to/readset4.2.bax.h5,path/to/readset4.3.bax.h5

### Nanopore, Nanopore CoV-Seq

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; mandatory; 
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Run: a unique ONT run name, usually has a structure similar to `PAE0000_a1b2c3d`; 
* Flowcell: code of the type of flowcell used, for example: the code for PromethION Flow Cell (R9.4) is `FLO-PRO002`;
* Library: code of the type of library preparation kit used, for example: the code for the Ligation Sequencing Kit is `SQK-LSK109`;
* Summary: path to the `sequencing_summary.txt` file outputted by the ONT basecaller; mandatory;
* FASTQ: path to the `fastq_pass` **directory**, that is usually created by the basecaller; mandatory;
* FAST5: path to the **directory** containing the raw fast5 files, before basecalling; 

Example: 

    Sample  Readset Run Flowcell    Library Summary FASTQ   FAST5
    sampleA readset1    PAE00001_abcd123    FLO-PRO002  SQK-LSK109 path/to/readset1_sequencing_summary.txt path/to/readset1/fastq_pass   path/to/readset1/fast5_pass 
    sampleA readset2    PAE00002_abcd456    FLO-PRO002  SQK-LSK109 path/to/readset2_sequencing_summary.txt path/to/readset2/fastq_pass   path/to/readset2/fast5_pass 
    sampleA readset3    PAE00003_abcd789    FLO-PRO002  SQK-LSK109 path/to/readset3_sequencing_summary.txt path/to/readset3/fastq_pass   path/to/readset3/fast5_pass 
    sampleA readset4    PAE00004_abcd246    FLO-PRO002  SQK-LSK109 path/to/readset4_sequencing_summary.txt path/to/readset4/fastq_pass   path/to/readset4/fast5_pass 


### For abacus users with Nanuq readsets
If your readsets belong to a [Nanuq](http://gqinnovationcenter.com/services/nanuq.aspx) project, use `$MUGQIC_PIPELINES_HOME/utils/nanuq2mugqic_pipelines.py` script to automatically create a Readset File and symlinks to your readsets on abacus.


Configuration Files
-------------------
Pipeline command parameters and cluster settings can be customized using Configuration Files (`.ini` extension).
Those files have a structure similar to Microsoft Windows INI files e.g.:
```
#!ini
[DEFAULT]
module_trimmomatic=mugqic/trimmomatic/0.36

[trimmomatic]
min_length=50
```

A parameter value is first searched in its specific section, then, if not found, in the special `DEFAULT` section.
The example above would resolve parameter `module_trimmomatic` value from section `trimmomatic` to `mugqic/trimmomatic/0.36`.

Configuration files support interpolation. For example:
```
#!ini
scientific_name=Homo_sapiens
assembly=GRCh37
assembly_dir=$MUGQIC_INSTALL_HOME/genomes/species/%(scientific_name)s.%(assembly)s
genome_fasta=%(assembly_dir)s/genome/%(scientific_name)s.%(assembly)s.fa
```
would resolve `genome_fasta` value to `$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa`.

Each pipeline has several configuration files in:
```
#!bash
$MUGQIC_PIPELINES_HOME/pipelines/<pipeline_name>/<pipeline_name>.*.ini
```
A default configuration file (`.base.ini` extension) is set for running on abacus cluster using *Homo sapiens* reference genome
and must always be passed first to the `--config` option.

You can also add a list of other configuration files to `--config`.
Files are read in the list order and each parameter value is overwritten if redefined in the next file.

This is useful to customize settings for a specific cluster or genome.
Each pipeline has a special configuration file for guillimin and mammouth clusters (`.guillimin.ini` and `.mammouth.ini` extensions respectively) in the same directory.
And various genome settings are available in `$MUGQIC_PIPELINES_HOME/resources/genomes/config/`.

For example, to run the DNA-Seq pipeline on guillimin cluster with *Mus musculus* reference genome:
```
#!bash
$MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.py --config $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.base.ini $MUGQIC_PIPELINES_HOME/pipelines/dnaseq/dnaseq.guillimin.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/Mus_musculus.GRCm38.ini ...
```


Design File
-----------
RNA-Seq, RNA-Seq De Novo Assembly and ChIP-Seq pipelines can perform differential expression analysis if they are provided with an input Design File.

The Design File is a TAB-separated values plain text file with one line per sample and the following columns:

### RNA-Seq and RNA-Seq De Novo Assembly
* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* Contrast: each of the following columns defines an experimental design contrast; the column name defines the contrast name, and the following values represent the sample group membership for this contrast:
    * '__0__' or '': the sample does not belong to any group;
    * '__1__': the sample belongs to the control group;
    * '__2__': the sample belongs to the treatment test case group.

Example:

    Sample	Contrast1	Contrast2	Contrast3
    sampleA	1	1	1
    sampleB	2	0	1
    sampleC	0	2	0
    sampleD	0	0	2

### Chip-Seq

* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* MarkName: Second Column; name of the histone mark; mandatory
* Contrast: each of the following columns defines an experimental design contrast; the column name defines the contrast name, and the following values represent the sample group membership for this contrast:
    * '__0__' or '': the sample does not belong to any group;
    * '__1__': the sample belongs to the control group;
    * '__2__': the sample belongs to the treatment test case group.

Example:

    Sample	MarkName Contrast1 Contrast2
    sampleA	H3K27ac 1   0
    sampleB	H3K27ac 1   0
    sampleC	H3K27ac 2   0
    sampleD	H3K27ac 2   0
    sampleA H3K4me3 0   1
    sampleB H3K4me3 0   1
    sampleC H3K4me3 0   2
    sampleD H3K4me3 0   2

HTML Analysis Report
--------------------
While pipelines are run, some jobs create a partial analysis report in [Markdown](http://daringfireball.net/projects/markdown/) format in
`<output_dir>/report/<pipeline_name>.<step_name>.md` e.g. `<output_dir>/report/DnaSeq.bwa_mem_picard_sort_sam.md`.

At any time during the pipeline processing, you can run the same pipeline command and add the option `--report`.
This will create a bash script calling the [Pandoc](http://pandoc.org/) converter to aggregate all partial Markdown reports already created into one single HTML document, which you can view in `<output_dir>/report/index.html`.

Thus, if the last pipeline steps fail, you will still get an HTML report containing sections for the first steps only.

The report title value can be overwritten in your copy of `$MUGQIC_PIPELINES_HOME/pipelines/<pipeline_name>/<pipeline_name>.base.ini` in section `[report]`.
You can also edit the partial Markdown reports before running the pandoc script, to add custom comments in your HTML report.

For developers: if you want to modify the Markdown report templates, they are all located in `$MUGQIC_PIPELINES_HOME/bfx/report/`.


PBS Job Logs
------------
When pipelines are run in PBS (Portable Batch System) job scheduler mode (default), a job list file is created in `<output_dir>/job_output/<PipelineName>_job_list_<timestamp>` and subsequent job log files are placed in `<output_dir>/job_output/<step_name>/<job_name>_<timestamp>.o` e.g.:
```
#!text
my_output_dir/job_output/
├── RnaSeqDeNovoAssembly_job_list_2014-09-30T19.52.29
├── trimmomatic
│   ├── trimmomatic.readset1_2014-09-30T19.52.29.o
│   └── trimmomatic.readset2_2014-09-30T19.52.29.o
├── trinity
│   └── trinity_2014-10-01T14.17.02.o
└── trinotate
    └── trinotate_2014-10-22T14.05.58.o
```

To view a TAB-separated values log report, use `$MUGQIC_PIPELINES_HOME/utils/log_report.pl` script by typing:
```
#!bash
$MUGQIC_PIPELINES_HOME/utils/log_report.pl <output_dir>/job_output/<PipelineName>_job_list_<timestamp>
```

which will output e.g.:
```
#!text
# Number of jobs: 41
#
# Number of successful jobs: 4
# Number of active jobs: 0
# Number of inactive jobs: 36
# Number of failed jobs: 1
#
# Execution time: 2014-09-30T19:52:58 - 2014-09-30T22:38:04 (2 h 45 min 6 s)
#
# Shortest job: merge_trimmomatic_stats (1 s)
# Longest job: insilico_read_normalization_readsets.readset2 (1 h 33 min 53 s)
#
# Lowest memory job: merge_trimmomatic_stats (0.00 GiB)
# Highest memory job: insilico_read_normalization_readsets.readset2 (31.32 GiB)
#
#JOB_ID JOB_FULL_ID    JOB_NAME    JOB_DEPENDENCIES    STATUS    JOB_EXIT_CODE    CMD_EXIT_CODE    REAL_TIME    START_DATE    END_DATE    CPU_TIME    CPU_REAL_TIME_RATIO    PHYSICAL_MEM    VIRTUAL_MEM    EXTRA_VIRTUAL_MEM_PCT    LIMITS    QUEUE    USERNAME    GROUP    SESSION    ACCOUNT    NODES    PATH
2100213.abacus2.ferrier.genome.mcgill.ca    2100213.abacus2.ferrier.genome.mcgill.ca    trimmomatic.readset1    SUCCESS    N/A    0    01:08:45 (1 h 8 min 45 s)    2014-09-30T19:52:58    2014-09-30T21:01:48    02:39:34 (2 h 39 min 34 s)    2.32    1.71 GiB    3.73 GiB    118.2 %    neednodes=1:ppn=6,nodes=1:ppn=6,walltime=24:00:00    sw    jfillon analyste    2465764    N/A    f3c10    /path/to/output_dir/job_output/trimmomatic/trimmomatic.readset1_2014-09-30T19.52.29.o
2100214.abacus2.ferrier.genome.mcgill.ca    2100214.abacus2.ferrier.genome.mcgill.ca    trimmomatic.readset2    SUCCESS    N/A    0    01:08:59 (1 h 8 min 59 s)    2014-09-30T19:52:58    2014-09-30T21:02:01    02:40:05 (2 h 40 min 5 s)    2.32    1.41 GiB    3.73 GiB    164.0 %    neednodes=1:ppn=6,nodes=1:ppn=6,walltime=24:00:00    sw    jfillon analyste    2465669    N/A    f3c10    /path/to/output_dir/job_output/trimmomatic/trimmomatic.readset2_2014-09-30T19.52.29.o
2100215.abacus2.ferrier.genome.mcgill.ca    2100215.abacus2.ferrier.genome.mcgill.ca    merge_trimmomatic_stats    2100213.abacus2.ferrier.genome.mcgill.ca:2100214.abacus2.ferrier.genome.mcgill.ca    SUCCESS    N/A    0    00:00:01 (1 s)    2014-09-30T21:04:06    2014-09-30T21:04:12    00:00:00 (0 s)    0.00    0.00 GiB    0.00 GiB    N/A    neednodes=1:ppn=1,nodes=1:ppn=1,walltime=120:00:00    sw    jfillon    analyste    3343994    N/A    f3c11    /path/to/output_dir/job_output/merge_trimmomatic_stats/merge_trimmomatic_stats_2014-09-30T19.52.29.o
2100216.abacus2.ferrier.genome.mcgill.ca    2100216.abacus2.ferrier.genome.mcgill.ca    insilico_read_normalization_readsets.readset1    2100213.abacus2.ferrier.genome.mcgill.ca    FAILED    N/A    N/A    00:38:16 (38 min 16 s)    2014-09-30T21:02:02    2014-09-30T21:40:23    04:50:10 (4 h 50 min 10 s)    7.58    30.71 GiB    32.32 GiB    5.3 %    neednodes=1:ppn=6,nodes=1:ppn=6,walltime=120:00:00    sw    jfillon    analyste    3343745    N/A    f3c11    /path/to/output_dir/job_output/insilico_read_normalization_readsets/insilico_read_normalization_readsets.readset1_2014-09-30T19.52.29.o
...
```

A Note about Cedar and Graham
------------
The default scheduler in GenPipes is the PBS scheduler. Cedar and Graham use the SLURM scheduler. To use GenPipes on Cedar or Graham, don't forget to add the "-j slurm" option.

Call home
---------
When pipeline jobs are submitted, a call home feature is invoked to collect some usage data. Those data are used to compute statistics and justify grant applications for funding support.

Data collected:

* Date and time
* Host and IP address
* Pipeline name
* Number of samples
* Pipeline steps


Contact us
----------
Please visit our [mailing list](https://groups.google.com/forum/#!forum/genpipes) to find questions and answers about GenPipes.

To subscribe to the mailing list and receive other people's messages, send an e-mail at [genpipes+subscribe@googlegroups.com](mailto:genpipes+subscribe@googlegroups.com).
You will receive an invitation which you must accept.

To use it, send us an e-mail at [genpipes@googlegroups.com](mailto:genpipes@googlegroups.com).

You can also report bugs at [pipelines@computationalgenomics.ca](mailto:pipelines@computationalgenomics.ca).

* Messages should not be sent directly to our team members. The generic e-mail addresses above are viewable by all of us and facilitate the follow-up of your request.
* Choose a meaningful subject for your message.
* Include the pipeline version number in your message (and the commit number if applicable).
* Provide the following information relevant to the problem encountered: the python command, the bash submission script, the output (job_outputs/*/*.o) file,
* An error message or code snippet illustrating your request is normally very useful.
