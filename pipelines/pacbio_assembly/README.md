[TOC]


PacBio Assembly Pipeline
========================

Contigs assembly with PacBio reads is done using what is refer as the HGAP workflow.
Briefly, raw subreads generated from raw .ba(s|x).h5 PacBio data files are filtered for quality.
A subread length cutoff value is extracted from subreads, depending on subreads distribution,
and used into the preassembly (aka correcting step) (BLASR) step which consists of aligning
short subreads on long subreads.
Since errors in PacBio reads is random, the alignment of multiple short reads on longer reads
allows to correct sequencing error on long reads.
These long corrected reads are then used as seeds into assembly (Celera assembler) which gives contigs.
These contigs are then *polished* by aligning raw reads on contigs (BLASR) that are then processed
through a variant calling algorithm (Quiver) that generates high quality consensus sequences
using local realignments and PacBio quality scores.

Prepare your readset file as described [here](https://bitbucket.org/mugqic/mugqic_pipelines/src#markdown-header-pacbio-assembly)
(if you use `nanuq2mugqic_pipelines.py`, you need to add and fill manually
the `EstimatedGenomeSize` column in your readset file).


Usage
-----
```
#!text

usage: pacbio_assembly.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                          [-o OUTPUT_DIR] [-j {pbs,batch,daemon,slurm}] [-f]
                          [--no-json] [--report] [--clean]
                          [-l {debug,info,warning,error,critical}]
                          [--sanity-check]
                          [--container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}]
                          [-r READSETS] [-v]

Version: 3.1.5

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

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
  --container {docker, singularity} {<CONTAINER PATH>, <CONTAINER NAME>}
                        run pipeline inside a container providing a container
                        image path or accessible docker/singularity hub path
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
```
![pacbio_assembly workflow diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_pacbio_assembly.resized.png)
[download full-size diagram](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_pacbio_assembly.png)
```
------
1- smrtanalysis_filtering
2- pacbio_tools_get_cutoff
3- preassembly
4- assembly
5- polishing
6- pacbio_tools_assembly_stats
7- blast
8- mummer
9- compile
10- circlator
11- basemodification
12- motifMaker

```
smrtanalysis_filtering
----------------------
Filter reads and subreads based on their length and QVs, using smrtpipe.py (from the SmrtAnalysis package).

1. fofnToSmrtpipeInput.py
2. modify RS_Filtering.xml files according to reads filtering values entered in .ini file
3. smrtpipe.py with filtering protocol
4. prinseq-lite.pl: write fasta file based on fastq file

Informative run metrics such as loading efficiency, read lengths and base quality are generated in this step as well.

pacbio_tools_get_cutoff
-----------------------
Cutoff value for splitting long reads from short reads is done here using
estimated coverage and estimated genome size.

You should estimate the overall coverage and length distribution for putting in
the correct options in the configuration file. You will need to decide a
length cutoff for the seeding reads. The optimum cutoff length will depend on
the distribution of the sequencing read lengths, the genome size and the
overall yield. Here, you provide a percentage value that corresponds to the
fraction of coverage you want to use as seeding reads.

First, loop through fasta sequences, put the length of each sequence in an array,
sort it, loop through it again and compute the cummulative length covered by each
sequence as we loop through the array. Once that length is > (coverage * genome
size) * $percentageCutoff (e.g. 0.10), we have our threshold. The idea is to
consider all reads above that threshold to be seeding reads to which will be
aligned lower shorter subreads.

preassembly
-----------
Having in hand a cutoff value, filtered reads are splitted between short and long reads. Short reads
are aligned against long reads and consensus (e.g. corrected reads) are generated from these alignments.

1. split reads between long and short
2. blasr: aligner for PacBio reads
3. m4topre: convert .m4 blasr output in .pre format
4. pbdagcon (aka HGAP2): generate corrected reads from alignments

assembly
--------
Corrected reads are assembled to generates contigs. Please see the
[Celera documentation](http://wgs-assembler.sourceforge.net/wiki/index.php?title=RunCA).
Quality of assembly seems to be highly sensitive to parameters you give Celera.

1. generate celera config files using parameters provided in the .ini file
2. fastqToCA: generate input file compatible with the Celera assembler
3. runCA: run the Celera assembler

polishing
---------
Align raw reads on the Celera assembly with BLASR. Load pulse information from bax or bas files into aligned file. Sort that file and run quiver (variantCaller.py).

1. generate fofn
2. upload Celera assembly with smrtpipe refUploader
3. compare sequences
4. load pulses
5. sort .cmp.h5 file
6. variantCaller.py

pacbio_tools_assembly_stats
---------------------------
blast
-----
Blast polished assembly against nr using dc-megablast.

mummer
------
Using MUMmer, align polished assembly against best hit from blast job. Also align polished assembly against itself to detect structure variation such as repeats, etc.

compile
-------
Compile assembly stats of all conditions used in the pipeline (useful when multiple assemblies are performed).

circlator
---------
Circularize the assembly contigs if possible.
User should launch this step after making sure the quality of the assembly is acceptable.

basemodification
----------------
Run ipdSummary.py for in silico detection of modified bases

motifMaker
----------
Run motifMaker to generate motif_summary.csv


