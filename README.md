GenPipes
================
This repository holds several bioinformatics pipelines developed at the [Canadian Centre for Computational Genomics](https://computationalgenomics.ca/) (C3G).

GenPipes consists of Python scripts which create a list of jobs running Bash commands. Those scripts support dependencies between jobs and a smart restart mechanism if some jobs fail during pipeline execution. Jobs can be submitted in different ways: by being sent to a scheduler like SLURM or PBS/Torque or by being run as a series of commands in batch through a Bash script. Job commands and parameters can be modified through several configuration files called "ini".

**For a more detailed tutorial on how to use GenPipes, please visit our [documentation page](https://genpipes.readthedocs.io/en/latest/).**  


On this page:

* [Software requirement](#software-requirement)
* [Quick setup for Abacus, Beluga, Narval, Graham and Cedar users](#quick-setup-for-abacus-beluga-narval-graham-and-cedar-users)
* [Download and setup for external users](#download-and-setup-for-external-users)
    * [Download](#download)
    * [Installation](#installation)
    * [GenPipes in a Container](#genpipes-in-a-container)
    * [Setup](#setup)
        * [Genomes](#genomes)
            * [Install a new Genome](#install-a-new-genome)
        * [Modules](#modules)
            * [Install a new Module](#install-a-new-module)
    * [Have autocompletion for GenPipes](#have-autocompletion-for-genpipes)
        * [Source the provided file](#source-the-provided-file)
        * [Generate autocomplete file yourself](#generate-autocomplete-file-yourself)
* [Usage](#usage)
* [Readset File](#readset-file)
    * [DNA-Seq, RNA-Seq, RNA-Seq De Novo Assembly, Amplicon-Seq, Methyl-Seq, CoV-Seq](#dna-seq-rna-seq-rna-seq-de-novo-assembly-amplicon-seq-methyl-seq-cov-seq)
    * [ChIP-Seq](#chip-seq)
    * [LongRead DNA-Seq, Nanopore CoV-Seq](#nanopore-nanopore-cov-seq)
    * [For abacus users with Nanuq readsets](#for-abacus-users-with-nanuq-readsets)
* [Configuration Files](#configuration-files)
* [Design File](#design-file)
    * [RNA-Seq and RNA-Seq De Novo Assembly](#rna-seq-and-rna-seq-de-novo-assembly)
    * [Chip-Seq](#chip-seq)
* [Batch File](#batch-file)
* [HTML Analysis Report](#html-analysis-report)
* [PBS/Slurm Job Logs](#pbsslurm-job-logs)
* [Contact](#contact)


Software requirement
--------------------
GenPipes has been tested with Python 3.11.1 and 3.12.2
It may work with other versions of python, but this has not been extensively tested. 


Quick setup for Abacus, Beluga, Narval, Graham and Cedar users
--------------------------------------------------------------
Genomes and modules used by the pipelines are already installed on a CVMFS partition mounted on all those clusters in `/cvmfs/soft.mugqic/root`.
To access them, add the following lines to your *$HOME/.bash_profile*:

```bash
umask 0006
# MUGQIC genomes and modules
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/root
module use $MUGQIC_INSTALL_HOME/modulefiles
```

For MUGQIC analysts, add the following lines to your *$HOME/.bash_profile*:

```bash
umask 0006
# MUGQIC genomes and modules for MUGQIC analysts
HOST=`hostname`;
DNSDOMAIN=`dnsdomainname`;
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/root
if [[ $HOST == abacus* || $DNSDOMAIN == ferrier.genome.mcgill.ca ]]; then
  export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev
elif [[ $HOST == ip* || $DNSDOMAIN == m  ]]; then
  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
elif [[ $HOST == cedar* || $DNSDOMAIN == cedar.computecanada.ca ]]; then
  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
elif [[ $HOST == beluga* || $DNSDOMAIN == beluga.computecanada.ca ]]; then
  export MUGQIC_INSTALL_HOME_DEV=/project/6007512/C3G/analyste_dev
elif [[ $HOST == narval* || $DNSDOMAIN == narval.computecanada.ca ]]; then
  export MUGQIC_INSTALL_HOME_DEV=/project/rrg-bourqueg-ad/C3G/analyste_dev
fi
module use $MUGQIC_INSTALL_HOME/modulefiles $MUGQIC_INSTALL_HOME_DEV/modulefiles
export RAP_ID=<my-rap-id>
```    

Also, set `JOB_MAIL` in your *$HOME/.bash_profile* to receive PBS job logs:
```bash
export JOB_MAIL=<my.name@my.email.ca>
```

GenPipes pipelines and compatible Python version are already installed as modules on those clusters.
To use them by default, add in your *$HOME/.bash_profile*:
```bash
module load mugqic/genpipes/<latest_version>
```
(find out the latest version with: "`module avail 2>&1 | grep mugqic/genpipes`").

For Alliance users you have to set your `RAP_ID` (Resource Allocation Project ID from Alliance) in your *$HOME/.bash_profile*:
```bash
export RAP_ID=<my-rap-id>
```


Download and setup for external users
-------------------------------------


### Download

Visit our [Download page](https://github.com/c3g/GenPipes/releases) to get the latest stable release.

If you want to use the most recent development version:
```bash
git clone https://github.com/c3g/GenPipes.git
```

### Installation

Optional, you can first create and use a virtual environment:
```bash
# For Alliance users you can load our python module by running the following. For others make sure you have python 3.11.1 or later installed
module load mugqic/python

# Create a virtual environment
python3 -m venv .genpipes_venv
# Activate the virtual environment
source .genpipes_venv/bin/activate
```

GenPipes can be installed via pip from PyPi repo:
```bash
pip install c3g-genpipes
```
Or from a GitHub clone, after downloading the repository as mentionned above:
```bash
cd genpipes
pip install .
```

The installation location may have to be added to your PATH, if it is not already on PATH. (See [Setup](#setup)).

> [!NOTE]  
> If you are using GenPipes in a PBS cluster not being Abacus, you have to ask your IT to add a config for epilogue of jobs. GenPipes uses `#PBS -T GenPipes` to generate epilogue of jobs using [this script](https://github.com/c3g/GenPipes/tree/6.0.0/genpipes/core/epilogue.py) and so your IT have to add it to the scheduler config accordingly.

### GenPipes in a Container:

You can use one of the following software to run the container: 'singularity', 'apptainer', 'docker' or 'podman'. You can get more information to run Genpipes in a Container [here](https://github.com/c3g/genpipes_in_a_container).

Once the GenPipes repo has been cloned and GenPipes has been installed, run the following command to download the container image and the settings file.

```bash
genpipes tools get_wrapper
```

Then run the usual GenPipes command but add the argument `--wrap` in order to have the genpipes being run inside the container (for sanity check to have access to all modules required). The GenPipes file will also use the container for each of his jobs when submitted to the scheduler.

To edit the config you have to edit the file `$MUGQIC_PIPELINES_HOME/resources/container/etc/wrapper.conf`, see [here](https://github.com/c3g/genpipes_in_a_container?tab=readme-ov-file#setup-a-giac-environment) for more details. For instance the container can be changed via the argument `GENPIPES_CONTAINERTYPE` and the GenPipes version to be used can be changed via the argument `GENPIPES_VERSION` (if other than latest, for more details see [here](https://github.com/c3g/genpipes_in_a_container?tab=readme-ov-file#genpipes-4-in-a-container)). You can also test a local version of GenPipes by using the option `GENPIPES_DIR`, for more details see [here](https://github.com/c3g/genpipes_in_a_container?tab=readme-ov-file#using-a-local-genpipes-version).

You can access inside the Genpipes container by running:

```bash
$MUGQIC_PIPELINES_HOME/resources/container/bin/container_wrapper.sh
```
Once inside you have access to all our modules from cvmfs and the latest version of GenPipes being loaded by default.

### Setup

Set `MUGQIC_PIPELINES_HOME` and `GENPIPES_INIS` to your local copy path, in your *$HOME/.bash_profile*:
```bash
export MUGQIC_PIPELINES_HOME=/path/to/your/local/genpipes
export GENPIPES_INIS=$MUGQIC_PIPELINES_HOME/genpipes/pipelines
```

Add the installation location to your path, if it is not already on path, in your *$HOME/.bash_profile*:
```bash
# for example:
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH
```

GenPipes (formerly called MUGQIC Pipelines) requires genomes and modules resources to run properly.
First, set `MUGQIC_INSTALL_HOME` to the directory where you want to install those resources, in your *$HOME/.bash_profile*:
```bash
# MUGQIC genomes and modules
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
```text
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
```bash
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
    * Ensembl: "*CHIMP2.1.4*"
    * UCSC: "*panTro4*"

* Retrieve the source version:
    * Ensembl: "78"
    * UCSC: unfortunately, UCSC does not have version numbers. Use [panTro4.2bit](http://hgdownload.soe.ucsc.edu/goldenPath/panTro4/bigZips/) date formatted as "YYYY-MM-DD": "2012-01-09"

* Run
```bash
cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh
```
e.g.:
* Ensembl:
```bash
cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/Pan_troglodytes.CHIMP2.1.4.sh
```
* UCSC:
```bash
cp $MUGQIC_PIPELINES_HOME/resources/genomes/GENOME_INSTALL_TEMPLATE.sh $MUGQIC_PIPELINES_HOME/resources/genomes/Pan_troglodytes.panTro4.sh
```

* Modify `$MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh` (`ASSEMBLY_SYNONYMS` can be left empty but if you know that 2 assemblies are identical apart from `chr` sequence prefixes, document it):
    * Ensembl:
    ```bash
    SPECIES=Pan_troglodytes   # With "_"; no space!
    COMMON_NAME=Chimpanzee
    ASSEMBLY=CHIMP2.1.4
    ASSEMBLY_SYNONYMS=panTro4
    SOURCE=Ensembl
    VERSION=78
    ```
    * UCSC:
    ```bash
    SPECIES=Pan_troglodytes   # With "_"; no space!
    COMMON_NAME=Chimpanzee
    ASSEMBLY=panTro4
    ASSEMBLY_SYNONYMS=CHIMP2.1.4
    SOURCE=UCSC
    VERSION=2012-01-09
    ```

* Running `bash $MUGQIC_PIPELINES_HOME/resources/genomes/<scientific_name>.<assembly>.sh` will install the genome in `$MUGQIC_INSTALL_HOME_DEV` (by default). This will download and install genomes, indexes and, for Ensembl only, annotations (GTF, VCF, etc.).
If the genome is big, separate batch jobs will be submitted to the cluster for bwa, bowtie/tophat, star indexing.
Check that jobs are completed.

* Add the newly created INI file to the genome config files for further usage in pipeline command:
```bash
cp $MUGQIC_INSTALL_HOME_DEV/genomes/species/<scientific_name>.<assembly>/<scientific_name>.<assembly>.ini $MUGQIC_PIPELINES_HOME/resources/genomes/config/
```


#### Modules
Software tools and associated modules must be installed in `$MUGQIC_INSTALL_HOME/software/` and `$MUGQIC_INSTALL_HOME/modulefiles/`.
Default software/module installation scripts are already available in `$MUGQIC_PIPELINES_HOME/resources/modules/`.

##### Install a new Module

New software tools and associated modules can be installed semi-automatically in a dev environment:

* Create a new `$MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh` file. You can have a look at existing module files in `$MUGQIC_PIPELINES_HOME/resources/modules/`.

* Modify `$MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh` following the instructions inside.

* Run `$MUGQIC_PIPELINES_HOME/resources/modules/<my_software>.sh` with no arguments. By default, it will download and extract the remote software archive, build the software and create the associated module, all in `$MUGQIC_INSTALL_HOME_DEV` if it is set.

* Check if the module is available with: `module avail 2>&1 | grep mugqic_dev/<my_software>/<version>`


Have autocompletion for GenPipes
--------------------------------

You can have tab completion for GenPipes by sourcing the right file provided by GenPipes or generate it yourself.

### Source the provided file
This will only last for the current session and you will have to source it again if you open a new terminal.
#### For bash users
```bash
source $MUGQIC_PIPELINES_HOME/resources/autocomplete/genpipes.bash
```
#### For zsh users
```bash
source $MUGQIC_PIPELINES_HOME/resources/autocomplete/genpipes.zsh
```
#### For tcsh users
```bash
source $MUGQIC_PIPELINES_HOME/resources/autocomplete/genpipes.tcsh
```

### Generate autocomplete file yourself
This will last even if you open a new terminal.
#### For bash users
```bash
genpipes -s bash > /usr/share/bash-completion/completions/genpipes
# OR if you are on a shared filesystem like The Alliance or Abacus
mkdir -p ~/.local/share/bash-completion/completions/
genpipes -s bash > ~/.local/share/bash-completion/completions/genpipes
```
#### For zsh users
```bash
genpipes -s zsh > /usr/local/share/zsh/site-functions/_genpipes
# OR if you are on a shared filesystem like The Alliance or Abacus
mkdir -p ~/.local/share/zsh/site-functions/_genpipes
genpipes -s zsh > ~/.local/share/zsh/site-functions/_genpipes
```


Usage
-----

For each pipeline, get help about usage, arguments and steps with:

* if you use a `mugqic/genpipes/<version>` module on our clusters or a local pip install, simply:
```bash
genpipes <pipeline_name> --help
```

Pipelines require as input one Readset File, one or more Configuration File(s) and possibly one Design File, all described below.

For documentation on how to use each of the pipelines, visit:
### [Pipelines Reference Guide](https://genpipes.readthedocs.io/en/latest/user_guide/pipeline_ref.html)

For more information about and source code for a specific pipeline, visit:

### [DNA-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/dnaseq)
### [RNA-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/rnaseq)
### [RNA-Seq De Novo Assembly Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/rnaseq_denovo_assembly)
### [RNA-Seq Light Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/rnaseq_light)
### [ChIP-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/chipseq)
### [Amplicon-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/ampliconseq)
### [Methyl-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/methylseq)
### [LongRead DNA-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/longread_dnaseq)
### [CoV-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/covseq)
### [Nanopore CoV-Seq Pipeline](https://github.com/c3g/GenPipes/tree/master/genpipes/pipelines/nanopore_covseq)

Readset File
------------

The Readset File is a TAB-separated values plain text file with one line per readset and the following columns in any order:


### DNA-Seq, RNA-Seq, RNA-Seq De Novo Assembly, Amplicon-Seq, Methyl-Seq, CoV-Seq

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
```text
Sample  Readset  Library RunType    Run    Lane Adapter1                          Adapter2                          QualityOffset BED              FASTQ1                            FASTQ2                            BAM
sampleA readset1 lib0001 PAIRED_END run100 1    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset1.paired1.fastq.gz path/to/readset1.paired2.fastq.gz path/to/readset1.bam
sampleA readset2 lib0001 PAIRED_END run100 2    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset2.paired1.fastq.gz path/to/readset2.paired2.fastq.gz path/to/readset2.bam
sampleB readset3 lib0002 PAIRED_END run200 5    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset3.paired1.fastq.gz path/to/readset3.paired2.fastq.gz path/to/readset3.bam
sampleB readset4 lib0002 PAIRED_END run200 6    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset4.paired1.fastq.gz path/to/readset4.paired2.fastq.gz path/to/readset4.bam
```

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
```text
Sample  Readset  MarkName MarkType Library RunType    Run    Lane Adapter1                          Adapter2                          QualityOffset BED              FASTQ1                            FASTQ2                            BAM
sampleA readset1 H3K27ac  N        lib0001 PAIRED_END run100 1    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset1.paired1.fastq.gz path/to/readset1.paired2.fastq.gz path/to/readset1.bam
sampleA readset2 H3K27ac  N        lib0001 PAIRED_END run100 2    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset2.paired1.fastq.gz path/to/readset2.paired2.fastq.gz path/to/readset2.bam
sampleB readset3 Input    I        lib0002 PAIRED_END run200 5    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset3.paired1.fastq.gz path/to/readset3.paired2.fastq.gz path/to/readset3.bam
sampleB readset4 Input    I        lib0002 PAIRED_END run200 6    AGATCGGAAGAGCACACGTCTGAACTCCAGTCA AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT 33            path/to/file.bed path/to/readset4.paired1.fastq.gz path/to/readset4.paired2.fastq.gz path/to/readset4.bam
```

### LongRead DNA-Seq, Nanopore CoV-Seq

* Sample: must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; mandatory; 
* Readset: a unique readset name with the same allowed characters as above; mandatory;
* Run: a unique ONT run name, usually has a structure similar to `PAE0000_a1b2c3d`; 
* Flowcell: code of the type of flowcell used, for example: the code for PromethION Flow Cell (R9.4) is `FLO-PRO002`;
* Library: code of the type of library preparation kit used, for example: the code for the Ligation Sequencing Kit is `SQK-LSK109`;
* Summary: path to the `sequencing_summary.txt` file outputted by the ONT basecaller; mandatory;
* FASTQ: path to the `fastq_pass` **directory**, that is usually created by the basecaller; mandatory;
* FAST5: path to the **directory** containing the raw fast5 files, before basecalling; 
* BAM: relative or absolute path to BAM file (for revio protocol), mandatory if FASTQ or FAST5 value is missing, ignored otherwise.

Example:
```text
Sample  Readset  Run              Flowcell   Library    Summary                                 FASTQ                       FAST5
sampleA readset1 PAE00001_abcd123 FLO-PRO002 SQK-LSK109 path/to/readset1_sequencing_summary.txt path/to/readset1/fastq_pass path/to/readset1/fast5_pass 
sampleA readset2 PAE00002_abcd456 FLO-PRO002 SQK-LSK109 path/to/readset2_sequencing_summary.txt path/to/readset2/fastq_pass path/to/readset2/fast5_pass 
sampleA readset3 PAE00003_abcd789 FLO-PRO002 SQK-LSK109 path/to/readset3_sequencing_summary.txt path/to/readset3/fastq_pass path/to/readset3/fast5_pass 
sampleA readset4 PAE00004_abcd246 FLO-PRO002 SQK-LSK109 path/to/readset4_sequencing_summary.txt path/to/readset4/fastq_pass path/to/readset4/fast5_pass 
```

### For abacus users with Nanuq readsets
If your readsets belong to a [Nanuq](http://gqinnovationcenter.com/services/nanuq.aspx) project, use `$MUGQIC_PIPELINES_HOME/genpipes/tools/nanuq2mugqic_pipelines.py` script to automatically create a Readset File and symlinks to your readsets on abacus.


Configuration Files
-------------------
Pipeline command parameters and cluster settings can be customized using Configuration Files (`.ini` extension).
Those files have a structure similar to Microsoft Windows INI files e.g.:
```bash
[DEFAULT]
module_trimmomatic=mugqic/trimmomatic/0.36

[trimmomatic]
min_length=50
```

A parameter value is first searched in its specific section, then, if not found, in the special `DEFAULT` section.
The example above would resolve parameter `module_trimmomatic` value from section `trimmomatic` to `mugqic/trimmomatic/0.36`.

Configuration files support interpolation. For example:
```bash
scientific_name=Homo_sapiens
assembly=GRCh37
assembly_dir=$MUGQIC_INSTALL_HOME/genomes/species/%(scientific_name)s.%(assembly)s
genome_fasta=%(assembly_dir)s/genome/%(scientific_name)s.%(assembly)s.fa
```
would resolve `genome_fasta` value to `$MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa`.

Each pipeline has several configuration files in:
```bash
$GENPIPES_INIS/<pipeline_name>/<pipeline_name>.*.ini
```
A default configuration file (`.base.ini` extension) is set for running on abacus cluster using *Homo sapiens* reference genome
and must always be passed first to the `--config` option.

You can also add a list of other configuration files to `--config`.
Files are read in the list order and each parameter value is overwritten if redefined in the next file.

This is useful to customize settings for a specific cluster or genome.
Each cluster has a special configuration file available (for example, `beluga.ini` and `narval.ini`) in the common_ini directory.
And various genome settings are available in `$MUGQIC_INSTALL_HOME/genomes/species/`.

For example, to run the DNA-Seq pipeline on the beluga cluster with *Mus musculus* reference genome:
```bash
genpipes dnaseq --config $GENPIPES_INIS/dnaseq/dnaseq.base.ini $GENPIPES_INIS/common_ini/beluga.ini $MUGQIC_INSTALL_HOME/genomes/species/Mus_musculus.GRCm38/Mus_musculus.GRCm38.ini ...
```


Design File
-----------
RNA-Seq, RNA-Seq De Novo Assembly, Methyl-Seq and ChIP-Seq pipelines can perform differential expression analysis if they are provided with an input Design File.

The Design File is a TAB-separated values plain text file with one line per sample and the following columns:

### RNA-Seq and RNA-Seq De Novo Assembly
* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* Contrast: each of the following columns defines an experimental design contrast; the column name defines the contrast name, and the following values represent the sample group membership for this contrast:
    * '__0__' or '': the sample does not belong to any group;
    * '__1__': the sample belongs to the control group;
    * '__2__': the sample belongs to the treatment test case group.

Example:
```text
Sample  Contrast1 Contrast2 Contrast3
sampleA 1         1         1
sampleB 2         0         1
sampleC 0         2         0
sampleD 0         0         2
```

### Chip-Seq

* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* MarkName: Second Column; name of the histone mark; mandatory
* Contrast: each of the following columns defines an experimental design contrast; the column name defines the contrast name, and the following values represent the sample group membership for this contrast:
    * '__0__' or '': the sample does not belong to any group;
    * '__1__': the sample belongs to the control group;
    * '__2__': the sample belongs to the treatment test case group.

Example:
```text
Sample  MarkName Contrast1 Contrast2
sampleA H3K27ac  1         0
sampleB H3K27ac  1         0
sampleC H3K27ac  2         0
sampleD H3K27ac  2         0
sampleA H3K4me3  0         1
sampleB H3K4me3  0         1
sampleC H3K4me3  0         2
sampleD H3K4me3  0         2
```

Batch File
-----------
RNA-Seq, RNA-Seq De Novo Assembly pipelines can perform batch effect correction if they are provided with an input Batch File.

The Batch File is a TAB-separated values plain text file with one line per sample and the following columns:

* Sample: first column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; the sample name must match a sample name in the readset file; mandatory;
* Batch: second (and last) column; must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only;
Example:
```text
Sample  Batch
sampleA 1
sampleB 1
sampleC 2
sampleD 2
sampleA 3
sampleB 3
sampleC 3
sampleD 3
```

HTML Analysis Report
--------------------
For most pipelines, metrics and logs are automatically collected to be displayed in a MultiQC report in HTML format.

This report will be saved under the `report` folder of the GenPipes output directory. 

If users wish to include additional log files that are not already captured, but are supported by MultiQC, we encourage you to reach out to our deverlopers by emailing [pipelines@computationalgenomics.ca](mailto:pipelines@computationalgenomics.ca). To view all tools currently supported by MultiQC, please visit the [MultiQC](https://multiqc.info/modules/) website. 


PBS/Slurm Job Logs
------------
When pipelines are run in PBS (Portable Batch System) or SLURM job scheduler mode (default), a job list file is created in `<output_dir>/job_output/<PipelineName>_job_list_<timestamp>` and subsequent job log files are placed in `<output_dir>/job_output/<step_name>/<job_name>_<timestamp>.o` e.g.:
```text
my_output_dir/job_output/
├── RnaSeqDeNovoAssembly_job_list_2014-09-30T19.52.29
├── trimmomatic
│   ├── trimmomatic.readset1_2014-09-30T19.52.29.o
│   └── trimmomatic.readset2_2014-09-30T19.52.29.o
├── trinity
│   └── trinity_2014-10-01T14.17.02.o
└── trinotate
    └── trinotate_2014-10-22T14.05.58.o
```

To view a TAB-separated values log report, use the `log_report` script by typing:
```bash
genpipes tools log_report <output_dir>/job_output/<PipelineName>_job_list_<timestamp>
```

which will output e.g.:
```text
id  name    status  user    node    priority    submit_time eligible_time   start_time  queue_time  end_time    total_time  time_limit  time_efficiency n_cpu   cpu_time    cpu_efficiency  mem ave_mem max_mem mem_efficiency  ave_diskr   max_diskr   ave_diskw   max_diskw   output_file_path
16462219 gatk_sam_to_fastq.tumorPair_COLO829T    COMPLETED   pstretenowich   f3u21c04    0   2025-01-09T10:03:30 2025-01-09T10:03:30 2025-01-09T10:04:10 00:00:40    2025-01-09T10:08:42 00:04:32    00:10:00    45.3%   5   00:14:25    318.0%  25.00 GB    Unknown 1.61 GB 6.5%    Unknown Unknown Unknown Unknown path_to/job_output/gatk_sam_to_fastq/gatk_sam_to_fastq.tumorPair_COLO829T_2025-01-09T10.03.29.o
16462220 trim_fastp.tumorPair_COLO829T   COMPLETED   pstretenowich   f3u21c01    0   2025-01-09T10:03:30 2025-01-09T10:08:42 2025-01-09T10:10:24 00:06:54    2025-01-09T10:12:12 00:01:48    00:20:00    9.0%    5   00:12:21    686.1%  25.00 GB    Unknown 1.84 GB 7.4%    Unknown Unknown Unknown Unknown path_to/job_output/trim_fastp/trim_fastp.tumorPair_COLO829T_2025-01-09T10.03.29.o
16462221 bwa_mem2_samtools_sort.tumorPair_COLO829T   COMPLETED   pstretenowich   f3u31c04    0   2025-01-09T10:03:31 2025-01-09T10:12:12 2025-01-09T10:12:27 00:08:56    2025-01-09T10:21:47 00:09:20    00:20:00    46.7%   9   00:52:34    563.2%  45.00 GB    Unknown 15.84 GB    35.2%   Unknown Unknown Unknown Unknown path_to/job_output/bwa_mem2_samtools_sort/bwa_mem2_samtools_sort.tumorPair_COLO829T_2025-01-09T10.03.29.o
...
```

A Note about non-Alliance Clusters
------------
The default scheduler in GenPipes is the SLURM scheduler. Beluga, Narval, Cedar and Graham use the SLURM scheduler. To use GenPipes on abacus, don't forget to add the "-j pbs" option.

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
