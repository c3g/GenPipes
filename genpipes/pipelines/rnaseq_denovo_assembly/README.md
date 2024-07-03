usage: genpipes rnaseq_denovo_assembly [-h] -c CONFIG [CONFIG ...] [-s STEPS]
                                       [-o OUTPUT_DIR]
                                       [-j {pbs,batch,daemon,slurm}] [-f]
                                       [--force_mem_per_cpu FORCE_MEM_PER_CPU]
                                       [--no-json] [--json-pt] [--report]
                                       [--clean]
                                       [--container {wrapper, singularity} <IMAGE PATH>]
                                       [--genpipes_file GENPIPES_FILE]
                                       [-l {debug,info,warning,error,critical}]
                                       [--sanity-check] [--wrap [WRAP]] -r
                                       READSETS_FILE [-d DESIGN_FILE] [-v]
                                       [-t {trinity,seq2fun}] [-b BATCH]

Version: 5.0.0

For more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/

options:
  -h, --help            show this help message and exit
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
  --force_mem_per_cpu FORCE_MEM_PER_CPU
                        Take the mem input in the ini file and force to have a
                        minimum of mem_per_cpu by correcting the number of cpu
                        (default: None)
  --no-json             do not create JSON file per analysed sample to track
                        the analysis status (default: false i.e. JSON file
                        will be created)
  --json-pt             create JSON file for project_tracking database
                        ingestion (default: false i.e. JSON file will NOT be
                        created)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  --container {wrapper, singularity} <IMAGE PATH>
                        Run inside a container providing a validsingularity
                        image path
  --genpipes_file GENPIPES_FILE, -g GENPIPES_FILE
                        Command file output path. This is the command used to
                        process the data, or said otherwise, this command will
                        "run the Genpipes pipeline". Will be redirected to
                        stdout if the option is not provided.
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  --sanity-check        run the pipeline in `sanity check mode` to verify that
                        all the input files needed for the pipeline to run are
                        available on the system (default: false)
  --wrap [WRAP]         Path to the genpipe cvmfs wrapper script. Default is g
                        enpipes/ressources/container/bin/container_wrapper.sh.
                        This is a convenience options for using genpipes in a
                        container
  -r READSETS_FILE, --readsets READSETS_FILE
                        readset file
  -d DESIGN_FILE, --design DESIGN_FILE
                        design file
  -v, --version         show the version information and exit
  -t {trinity,seq2fun}, --type {trinity,seq2fun}
                        RNAseq analysis type
  -b BATCH, --batch BATCH
                        batch file (to peform batch effect correction

Protocol trinity
0 picard_sam_to_fastq
1 trimmomatic
2 merge_trimmomatic_stats
3 insilico_read_normalization_readsets
4 insilico_read_normalization_all
5 trinity
6 exonerate_fastasplit
7 blastx_trinity_uniprot
8 blastx_trinity_uniprot_merge
9 transdecoder
10 hmmer
11 infernal_transcriptome
12 blastp_transdecoder_uniprot
13 signalp
14 tmhmm
15 trinotate
16 align_and_estimate_abundance_prep_reference
17 align_and_estimate_abundance
18 gq_seq_utils_exploratory_analysis_rnaseq_denovo
19 differential_expression
20 filter_annotated_components
21 gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered
22 differential_expression_filtered
23 multiqc
Protocol seq2fun
0 picard_sam_to_fastq
1 merge_fastq
2 seq2fun
3 differential_expression_seq2fun
4 pathway_enrichment_seq2funpicard_sam_to_fastq 
-------------------
 
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

trimmomatic 
-----------
 
        Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', get 'adapter_fasta'),
        it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
        reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
        only Adapter1 is used and left unchanged.

        This step takes as input files:

utput_dir
        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

merge_trimmomatic_stats 
-----------------------
 
The trim statistics per readset are merged at this step.

insilico_read_normalization_readsets 
------------------------------------
 
Normalize each readset, using the Trinity normalization utility.

insilico_read_normalization_all 
-------------------------------
 
Merge all normalized readsets together and normalize the result, using the Trinity normalization utility.

trinity 
-------
 
Create a de novo assembly from normalized readsets using [Trinity](http://trinityrnaseq.sourceforge.net/).

exonerate_fastasplit 
--------------------
 
Split the Trinity assembly FASTA into chunks for further parallel BLAST annotations.

blastx_trinity_uniprot 
----------------------
 
Annotate Trinity FASTA chunks with Swiss-Prot and UniRef databases using [blastx](http://blast.ncbi.nlm.nih.gov/).

blastx_trinity_uniprot_merge 
----------------------------
 
Merge blastx Swiss-Prot and UniRef chunks results.

transdecoder 
------------
 
Identifies candidate coding regions within transcript sequences using [Transdecoder](http://transdecoder.github.io/).

hmmer 
-----
 
Identifies protein domains using [HMMR](http://hmmer.janelia.org/).

infernal_transcriptome 
----------------------
 
Identify structural RNAs using cmscan function from [infernal](http://eddylab.org/infernal)
Run in parallel, using chunks created by exonerate. 

blastp_transdecoder_uniprot 
---------------------------
 
Search Transdecoder-predicted coding regions for sequence homologies on UniProt using [blastp](http://blast.ncbi.nlm.nih.gov/).

signalp 
-------
 
Predict signal peptides using [SignalP](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp).

tmhmm 
-----
 
Predict transmembrane regions using [TMHMM](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm).

trinotate 
---------
 
Perform transcriptome functional annotation and analysis using [Trinotate](http://trinotate.sourceforge.net/).
All functional annotation data is integrated into a SQLite database and a whole annotation report is created.

align_and_estimate_abundance_prep_reference 
-------------------------------------------
 
Index Trinity FASTA file for further abundance estimation using [Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).

align_and_estimate_abundance 
----------------------------
 
Estimate transcript abundance using [RSEM](http://deweylab.biostat.wisc.edu/rsem/) via
[Trinity align_and_estimate_abundance.pl utility](http://trinityrnaseq.sourceforge.net/analysis/abundance_estimation.html).

gq_seq_utils_exploratory_analysis_rnaseq_denovo 
-----------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package.

differential_expression 
-----------------------
 
Performs differential gene expression analysis using [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file. Also, performs Gene Ontology analysis for RNA-Seq denovo Assembly using the Bioconductor's R package [goseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html).
Generates GO annotations for differential genes and isoforms expression analysis, based on associated GOTERMS generated by trinotate.

filter_annotated_components 
---------------------------
 
Filter high quality contigs based on values in trinotate annotations. Recreate a high quality contigs fasta file and run Assembly statistics using the gqSeqUtils R package.

gq_seq_utils_exploratory_analysis_rnaseq_denovo_filtered 
--------------------------------------------------------
 
Exploratory analysis using the gqSeqUtils R package using a subset of filtered transcripts.

differential_expression_filtered 
--------------------------------
 
Differential Expression and GOSEQ analysis based on filtered transcripts and genes.

multiqc 
-------
 
Aggregate results from bioinformatics analyses across many samples into a single report.
MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).

picard_sam_to_fastq 
-------------------
 
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.

merge_fastq 
-----------
 
This step is performed to merge fastq files if multiple readset files for one sample is present

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

seq2fun 
-------
 
seq2fun

This step takes as input files:

1. FASTQ files from the readset file if available
2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

This step perform seq2fun analysis and generates output files including KO abundance table and KO mapped fastq files
(https://www.seq2fun.ca/manual.xhtml#sect4) and (https://www.seq2fun.ca/manual.xhtml#sect20)

For each contrast different folders and all the files for that particular contrast are
generated. Therefore, only pairwise comparisons are possible
(treatment and controls will be added according to the 1 and 2 in the design file).

differential_expression_seq2fun 
-------------------------------
 
Performs differential gene expression analysis using [DESEQ2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [EDGER](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html).
Merge the results of the analysis in a single csv file.

pathway_enrichment_seq2fun 
--------------------------
 

seq2fun pathway analysis using fgsea (https://bioconductor.org/packages/release/bioc/html/fgsea.html)
 and user provide universal pathway list as KEGG map ID. The differential KO expression results obtained
 from edgeR will be using as the input for the pathway enrichment analysis
