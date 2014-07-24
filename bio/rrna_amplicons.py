#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def merge_barcodes(reads1, reads2, outdir):
    outfile1 = outdir + "/reads1.fastq"
    outfile2 = outdir + "/reads2.fastq"
    outfile_paired = outdir + "/paired.fastq"

    job = Job(
        reads1 + reads2, 
        [outfile1, outfile2],
        [['tools', 'moduleVersion.tools'], ['memtime', 'moduleVersion.memtime']]
    )

    job.command = """\
cat \\\n  {reads1} \\\n | gunzip -c > {outfile1} && \\
cat \\\n  {reads2} \\\n | gunzip -c > {outfile2} && \\
memtime joinPairedFastqs.pl \\
  --reads_1 {outfile1} \\
  --reads_2 {outfile2} \\
  --outfile {outfile_paired}""".format(
        reads1 = " \\\n  ".join(reads1),
        reads2 = " \\\n  ".join(reads2),
        outfile1 = outfile1,
        outfile2 = outfile2,
        outfile_paired = outfile_paired
    )

    return job

def merge_barcodes_single_end_reads(reads, outdir, pair):
    outfile = outdir + "/reads" + pair + ".fastq"

    job = Job(reads, [outfile])

    job.command = """\
cat \\\n  {reads} | gunzip -c > {outfile}""".format(
        reads1 = " \\\n  ".join(reads1),
        reads2 = " \\\n  ".join(reads2),
        outfile1 = outfile1,
        outfile2 = outfile2,
        outfile_paired = outfile_paired
    )

    return job

#### TO BE CONVERTED TO PYTHON
def duk_wrapper(infile_fastq, contam, ncontam, log, db):
    
    job = Job( 
        [infile_fastq], 
        [contam, ncontam],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl'],
            ['duk', 'moduleVersion.duk']
          ]
    )

    job.command="""
memtime contamWrapper.pl \\
--infile {infile_fastq} \\
--outfile_matched {contam} \\
--outfile_unmatched {ncontam} \\
--log {log} \\
--db {db} \\
--num_threads {num_threads}""".format(
    infile_fastq = infile_fastq,
    outfile_matched = outfile_matched,
    outfile_unmatched = outfile_unmatched,
    log = log,
    db = db,
    num_threads = config.param('duk_wrapper', 'num_threads', type='int')
    )
    return job


def duk(log, ncontam, contam, db, infile):
                
    job = Job(
        [infile], 
        [contam, ncontam, log],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl'],
            ['duk', 'moduleVersion.duk']
        ]
    )
        
    job.command=""" 
memtime gunzip -c {infile} | duk \\
-o log {log} \\
-n ncontam {ncontam}\\
-m contam {contam}\\
-k {k} \\
-s {d} \\
-c {c} \\
db""".format(
    log = log,
    ncontam = ncontam,
    contam = contam,
    k = config.param('duk', 'k', 'int'),
    s = config.param('duk', 's', 'int'),
    c = config.param('duk', 'c', 'int'),
    ) 
    return job

def split_barcodes(infile, barcodes, outfile, log):
    
    job = Job(
        [infile],
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )

    job.command="""
memtime barcodes.pl \\
--infile {infile} \\
--barcodes {barcodes} \\
--outfile {outfile} \\
--num_threads {num_threads}  \\
--log {log};""".format(
    infile = infile,
    outfile = outfile,
    num_threads = config.param( 'barcodes', 'num_threads', 1, 'int'),
    log = {log}
    )
    
    return job

def removeUnpairedReads(infile, outfilePaired, unpairedR1, unpairedR2):
    job = Job(
        [infile], 
        [outfilePaired],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )
    job.command="""
memtime removeUnpaired.pl \\
--infile {infile} \\
--outfile_paired {outfilePaired} \\
--outfile_1 {unpairedR1} \\
--outfile_2 {unpairedR2} \\
--num_threads {num_threads}""".format(
    infile = infile,
    outfile_paired = outfile_paired,
    outfile_1 = outfile_1,
    outfile_2 = outfile_2,
    num_threads =  config.param( 'remove_unpaired', 'num_threads', 1, 'int');
    )
                
    return job

def splitPairs(infile, outfileR1, outfileR2):
        
    job = Job(
        [infile], 
        [outfileR1, outfileR2]
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )

    job.command="""
memtime splitPairs.pl \\
--infile {infile} \\
--outfile_1 {outfileR1} \\
--outfile_2 {outfileR2} \\
--num_threads {num_threads}""".format(
    infile = infile,
    outfile_1 = outfile_1,
    outfile_2 = outfile_2,
    num_threads = config.param( 'split_pairs', 'num_threads', 'int')
    )
                
                    
        return job

def generateQscoreSheet(infile, prefix, log, outfile, barcodes):
                
    job = Job(
        [infile],
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['fastx', 'moduleVersion.fastx'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )

    job.command="""
memtime qscoreSheets.pl \\ 
 --fastq {infile} \\
 --tmp {tmp} \\
 --prefix {prefix} \\
 --suffix {suffix} \\
 --log {log} \\
 --outfile {outfile} \\
 --phred {phred} \\
 --barcodes {barcodes} \\
 --num_threads {num_threads}""".format(
    infile = infile,
    tmp = tmp,
    suffix = suffix,
    log = log,
    outfile = outfile,
    tmp = config.param('DEFAULT', 'tmpDir', 'dirpath'),
    phred = config.param('DEFAULT', 'qual', 'int'),
    config.param( 'qscore_sheet', 'num_threads', 'int')

    )
                    
    return job

def generateQscoreGraphSingle(infile, prefix, outfile):
        
    job = Job(
        [infile] , 
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['R', 'moduleVersion.R'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )
    job.command="""
memtime qscorePlots.pl \\
--infile_1 {infile} \\
--name {prefix} \\
--pdf {outfile} \\
--display 1 \\
--single""".format(
    infile = infile,
    prefix = prefix,
    outfile = outfile
    )

    return job

def generateQscoreGraphPaired(infileR1, infileR2, outfile):
                
    job = Job(
        [infileR1, infileR2], 
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['R', 'moduleVersion.R'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )

    job.command="""
memtime qscorePlots.pl
--infile_1 {infileR1} \\
--infile_2 {infileR2} \\
--name qual_stats \\
--pdf {outfile} \\
--display 1 \\
--paired""".format(
    infileR1 = infileR1,
    infileR2 = infileR2,
    outfile = outfile
    )

    return job

def cutReads(infile, begin, end, outfile):
        
    job = Job(
        [infile],
        [outfile]
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]

    job.command="""
memtime cutFastqSeq.pl \\
--infile {infile} \\
--begin {begin} \\
--end {end} \\
--outfile {outfile}""".format(
    infile = infile,
    begin = begin,
    end = end,
    outfile = outfile
    )

    return job

def flash(infileR1, infileR2, prefix, outdir):
                
        job = Job(
            [infileR1, infileR2], 
            [outdir.'/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq'],
            [
                ['memtime', 'moduleVersion.memtime'],
                ['flash', 'moduleVersion.flash'],
                ['tools', 'moduleVersion.tools'],
                ['perl', 'moduleVersion.perl']
            ]
    job.command="""
memtime flash.pl \\
--infile_1 {infileR1} \\
--infile_2 {infileR2} \\
--prefix {prefix} \\
--outdir {outdir} \\
--n {n} \\
--m {m} \\
--M {M} \\
--x {x} \\
--p {p} \\
--num_threads {num_threads}""".format(
    infileR1 = infileR1,
    infileR2 = infileR2,
    prefix = prefix,
    outdir = outdir,
    n = config.getparam('flash', 'sampling', 'int'),
    m = config.getparam('flash', 'minOverlap', 'int'),
    M = config.getparam('flash', 'maxOverlap', 'int'),
    x = config.getparam('flash', 'percentMismatch', 'float'),
    p = config.getparam('flash', 'phred', 'int'),
    num_threads = config.getparam( 'flash', 'num_threads', 'int')
    )
    return job

def removePrimers(infile, revPrimer, fwdPrimer, outfile,  outfileFailed):
        
    job = Job(
        [infile], 
        [outfile],
        [
            ['memtime', 'moduleVersion.memtime'],
            ['tools', 'moduleVersion.tools'],
            ['perl', 'moduleVersion.perl']
        ]
    )

    job.command="""
memtime itagsQC.pl \\
 --infile {infile} \\
if(revPrimer ne 'null'){
                 --primer_3_prime revPrimer;
                 --length_3_prime ' . LoadConfig::getParam( 'itags_QC', 'length3Prime', 1, 'int'); 
            }if(fwdPrimer ne 'null'){
                             --primer_5_prime fwdPrimer;
                             --length_5_prime ' . LoadConfig::getParam( 'itags_QC', 'length5Prime', 1, 'int') ;    
                        }
# --qscore_1 ' . LoadConfig::getParam( 'itags_QC', 'qscore1');
# --qscore_2 ' . LoadConfig::getParam( 'itags_QC', 'qscore2');
 --outfile outfile;
 --outfile_failed outfileFailed;
 --num_threads ' . LoadConfig::getParam( 'itags_QC', 'num_threads', 1, 'int');
 --qual ' . LoadConfig::getParam( 'default', 'qual', 1, 'int');
# --lq_threshold ' . LoadConfig::getParam( 'itags_QC', 'lq_threshold');
 --primer_mismatch ' . LoadConfig::getParam( 'itags_QC', 'primerMismatch', 1, 'int');
# --min_length ' . LoadConfig::getParam( 'itags_QC', 'minlength');
# --N ' . LoadConfig::getParam( 'itags_QC', 'N');""".format(

         )
        return job

def itagsQC(infile, revPrimer, fwdPrimer, outfile, outfileFailed):
        
        job = Job([infile] , [outfile]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                    itagsQC.pl
                     --infile infile;
                    if(revPrimer ne 'null'){
                                     --primer_3_prime revPrimer;
                                     --length_3_prime ' . LoadConfig::getParam( 'itags_QC', 'length3Prime', 1, 'int'); 
                                }if(fwdPrimer ne 'null'){
                                                 --primer_5_prime fwdPrimer;
                                                 --length_5_prime ' . LoadConfig::getParam( 'itags_QC', 'length5Prime', 1, 'int') ;    
                                            }
                     --qscore_1 ' . LoadConfig::getParam( 'itags_QC', 'qscore1', 1, 'int');
                     --qscore_2 ' . LoadConfig::getParam( 'itags_QC', 'qscore2', 1, 'int');
                     --outfile outfile;
                     --outfile_failed outfileFailed;
                     --num_threads ' . LoadConfig::getParam( 'itags_QC', 'num_threads', 1, 'int');
                     --qual ' . LoadConfig::getParam( 'default', 'qual', 1, 'int');
                     --lq_threshold ' . LoadConfig::getParam( 'itags_QC', 'lq_threshold', 1, 'int');
                     --primer_mismatch ' . LoadConfig::getParam( 'itags_QC', 'primerMismatch', 1, 'float');
                     --min_length ' . LoadConfig::getParam( 'itags_QC', 'minlength', 1, 'int');
                     --N ' . LoadConfig::getParam( 'itags_QC', 'N', 1, 'int');""".format(

         )
                
                    
        return job

def countReport(rA_files, rA_names, analysisType, barcodesDist, OTUtable, obsTable, outfile):
        
    job = Job([OTUtable], [outfile]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                    countReport.pl
                    foreach(@rA_files){
                                     --file ' ._;
                                }
                    foreach(@rA_names){
                                     --name ' ._;
                                }
                     --analysisType analysisType;
                     --barcodesDist '. barcodesDist;
                     --OTUtable OTUtable;
                     --obsTable obsTable;
                     > ' .outfile;""".format(

         )
                
                    
        return job

def txtToPdf(infile, outfile):
                
        job = Job([infile], [outfile]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                    txtToPdf.pl
                     --infile infile;
                     --outfile outfile;""".format(

         )
        return job

def mergePdf(command):
        
        my dummyOutfile       = "mergepdf.mugqic.done";
        
        job = Job([""] , [dummyOutfile]);
        
        
                    
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['ghostscript', 'moduleVersion.ghostscript'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                     command;
                     && touch dummyOutfile; """.format(

         )

                    
        return job

def clustering1(infile, barcodes, outdir):
        
        job = Job([infile], ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['usearch', 'moduleVersion.usearch'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                     clustering1.pl
                     --infile_fastq infile;
                     --ref_db ' . LoadConfig::getParam( 'DB', 'chimeras', 1, 'filepath'); 
                     --barcodes barcodes;
                     --outdir outdir;
                     --num_threads ' . LoadConfig::getParam( 'clustering', 'num_threads', 1, 'int');
                    # --start_at 4""".format(

         )
        return job

def clustering2(infile, barcodes, outdir):
        
        job = Job([infile], ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['usearch', 'moduleVersion.usearch'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                     clustering2.pl
                     --infile_fastq infile;
                     --ref_db ' . LoadConfig::getParam( 'DB', 'chimeras', 1, 'path'); 
                     --barcodes barcodes;
                     --outdir outdir;
                     --num_threads ' . LoadConfig::getParam( 'clustering', 'num_threads', 1, 'int');
                    # --start_at 4""".format(

         )
        return job

def clustering3(infile, barcodes, outdir):
        
        job = Job([infile], ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"]);
                [
                          ['memtime', 'moduleVersion.memtime'],
                          ['usearch', 'moduleVersion.usearch'],
                          ['dnaclust', 'moduleVersion.dnaclust'],
                          ['tools', 'moduleVersion.tools'],
                          ['perl', 'moduleVersion.perl']
                        ]
    job.command="""
                     memtime 
                     clustering3.pl
                     --infile_fastq infile;
                     --ref_db ' . LoadConfig::getParam( 'DB', 'chimeras', 1, 'filepath'); 
                     --barcodes barcodes;
                     --outdir outdir;
                 --lowAbunCutOff ' . LoadConfig::getParam( 'clustering', 'lowAbunCutOff', 1, 'int');
                     --num_threads ' . LoadConfig::getParam( 'clustering', 'num_threads', 1, 'int');
                    # --start_at 4""".format(

         )
        return job

def clientReport(iniFilePath, projectPath, pipelineType, reportPath):
        
      my pipeline = 'pipeline=\"' .pipelineType .'\",
        my titleTMP = LoadConfig::getParam( 'report','projectName');
        my title = 'report.title=\"' .titleTMP .'\",
        my authorTMP = LoadConfig::getParam( 'report','report.author');
        my author = 'report.author=\"' .authorTMP .'\",
        my contactTMP = LoadConfig::getParam( 'report','report.contact');
        my contact = 'report.contact=\"' .contactTMP .'\",

        
        #job = Job([iniFilePath],[projectPath]]);
        ro_job->setUp2Date(0);

        
                    
                [
                          ['R', 'moduleVersion.R']
                        ]
    job.command="""
                     R --no-save -e \'library(gqSeqUtils) ;
                     mugqicPipelineReport(
                     pipeline=\"pipelineType . '\",
                     report.path=\"reportPath . '\",
                     ini.file.path=\"iniFilePath . '\",' ;
                     report.title=\"' . LoadConfig::getParam( 'report','projectName') . '\",' ;
                     report.author=\"' . LoadConfig::getParam( 'report','report.author') . '\",' ;
                     report.contact=\"' . LoadConfig::getParam( 'report','report.contact') . '\",' ;
                     project.path=\"projectPath . '\")\'' ;""".format(

         )

                    
        return job

def cleanup(tmpdir):
        
        job = Job([""] , [""]);
                [
                          ['memtime', 'moduleVersion.memtime']
                        ]
    job.command="""
                     memtime 
                    rm '.tmpdir.' -rf""".format(

         )
        return job

def templateSub(outdir):
        
        job = Job(undef , undef);
                [
                          ['memtime', 'moduleVersion.memtime']
                        ]
    job.command="""
                     memtime """.format(

         )
        return job


