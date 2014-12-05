PacBio Assembly Pipeline
========================


Overview
--------
Contigs assembly with PacBio reads is done using what is refer as the HGAP workflow. Briefly, raw subreads generated from raw .ba(s|x).h5 PacBio data files are filtered for quality. A subread length cutoff value is extracted from subreads, depending on subreads distribution, and used into the preassembly (aka correcting step) (BLASR) step which consists of aligning short subreads on long subreads. Since errors in PacBio reads is random, the alignment of multiple short reads on longer reads allows to correct sequencing error on long reads. These long corrected reads are then used as seeds into assembly (Celera assembler) which gives contigs. These contigs are then *polished* by aligning raw reads on contigs (BLASR) that are then processed through a variant calling algorithm (Quiver) that generates high quality consensus sequences using local realignments and PacBio quality scores.


On this page:

[TOC]


Usage
-----
Prepare your readset file as described [here](https://bitbucket.org/mugqic/mugqic_pipelines/src/#markdown-header-pacbio-assembly) (if you use `nanuq2mugqic_pipelines.py`, you need to add and fill manually the `EstimatedGenomeSize` column in your readset file).

```
#!bash
mugqic_pipelines/pipelines/pacbio_assembly/pacbio_assembly.py --help
```


Steps
-----
### filtering
Filtering. This step will filter reads and subreads based on their length and QVs. Informative run metrics such as loading efficiency, readlengths, and base quality are generated in this step as well.

### getStats
Cutoff value for splitting long reads from short reads is done here using estimated coverage and estimated genome size.

You should estimate the overall coverage and length distribution for putting in the correct options in the configuration file. You will need to decide a length cutoff for the seeding reads. The optimum cutoff length will
depend on the distribution of the sequencing read lengths, the genome size and the overall yield. Here, you provide a percentage value that corresponds to the fraction of coverage you want to use as seeding reads.

First, loop through fasta sequences. put the length of each sequence in an array, sort it, loop through it again and compute the cummulative length coveredby each sequence as we loop though the array. Once that length is >
(coverage * genome size) * $percentageCutoff (e.g. 0.10), we have our threshold. The idea is to consider all reads above that threshold to be seeding reads to which will be align lower shorter subreads.

### preAssembly (aka correction)
Having in hand a cutoff value, filtered reads are splitted between short and long reads. Short reads are aligned against long reads and consensus (e.g. corrected reads) are generated from these alignments:

1. split reads between long and short;
2. blasr (aligner for PacBio reads);
3. m4topre (converts .m4 blasr output in .pre format);
4. pbdagcon (aka HGAP2) (generates corrected reads from alignments).

### assembly
Corrected reads are assembled to generates contigs. Please see the [Celera documentation](http://wgs-assembler.sourceforge.net/wiki/index.php/RunCA). Quality of assembly seems to be highly sensitive to parameters you give to Celera:

1. Generate celera config files using paramters provided in the .ini file;
2. fastqToCA: generates input file compatible with the Celera assembler;
3. runCA: run the Celera assembler.

### polishing
Align raw reads on the Celera assembly with BLASR. Load pulse information from bax or bas files into aligned file. Sort that file and run quiver (variantCaller.py):

1. Generate fofn;
2. Upload Celera assembly with smrtpipe refUploader;
3. Compare sequences;
4. Load pulses;
5. Sort .cmp.h5 file;
6. variantCaller.py.

### blast
Blast polished assembly against nr using dc-megablast.

### mummer

Using MUMmer, align polished assembly against best hit from blast job. Also align polished assembly against itself to detect structure variation such as repeats, etc.

### report
Generates summary tables and generates MUGQIC style nozzle report.

### compile (useful when multiple assemblies are performed)
Compile assembly stats of all conditions used in the pipeline.
