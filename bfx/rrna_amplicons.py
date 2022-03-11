################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

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
        [['tools', 'module_mugqic_tools'], ['memtime', 'module_memtime']]
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
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl'],
            ['duk', 'module_duk']
        ]
    )

    job.command="""\
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
    num_threads = config.param('duk_wrapper', 'num_threads', param_type='int')
    )
    return job


def duk(log, ncontam, contam, db, infile):
                
    job = Job(
        [infile], 
        [contam, ncontam, log],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl'],
            ['duk', 'module_duk']
        ]
    )
        
    job.command="""\
memtime gunzip -c {infile} | duk \\
  -o log {log} \\
  -n ncontam {ncontam} \\
  -m contam {contam} \\
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
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
memtime barcodes.pl \\
  --infile {infile} \\
  --barcodes {barcodes} \\
  --outfile {outfile} \\
  --num_threads {num_threads}  \\
  --log {log};""".format(
    infile = infile,
    outfile = outfile,
    num_threads = config.param( 'barcodes', 'num_threads', 1, 'int'),
    log = log
    )
    
    return job

def removeUnpairedReads(infile, outfilePaired, unpairedR1, unpairedR2):
    job = Job(
        [infile], 
        [outfilePaired],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
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
    num_threads =  config.param( 'remove_unpaired', 'num_threads', 'int')
    )
                
    return job

def splitPairs(infile, outfileR1, outfileR2):
        
    job = Job(
        [infile], 
        [outfileR1, outfileR2]
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
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
            ['memtime', 'module_memtime'],
            ['fastx', 'module_fastx'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
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
        suffix = suffix,
        log = log,
        outfile = outfile,
        tmp = config.param('default', 'tmpDir', 'dirpath'),
        phred = config.param('default', 'qual', 'int'),
        num_threads = config.param('qscore_sheet', 'num_threads', 'int')
    )
                    
    return job

def generateQscoreGraphSingle(infile, prefix, outfile):
        
    job = Job(
        [infile] , 
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['R', 'module_R'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
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
            ['memtime', 'module_memtime'],
            ['R', 'module_R'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
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
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
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
        [outdir + "/assembly_complete/ncontam_nphix_trimmed.extendedFrags.fastq"],
        [
            ['memtime', 'module_memtime'],
            ['flash', 'module_flash'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
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
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
memtime itagsQC.pl \\
  --infile {infile} \\""".format(
    infile = infile
    )

    if(revPrimer != "null"):
        job.command+="""\
  --primer_3_prime {revPrimer} \\
  --length_3_prime {length_3_prime} \\""".format(
        revPrimer = revPrimer,
        length_3_prime = config.param('itags_QC', 'length3Prime', 'int')
        )
    
    if(fwdPrimer != "null"):
        job.command+="""\
  --primer_5_prime {fwdPrimer} \\
  --length_5_prime {length_5_prime} \\""".format(
        fwdPrimer = fwdPrimer,
        length_5_prime = config.param('itags_QC', 'length5Prime', 'int')
        )
 
    job.command+="""\
  --outfile {outfile} \\
  --outfile_failed {outfileFailed} \\
  --num_threads {num_threads} \\
  --qual {qual} \\
  --primer_mismatch {primer_mismatch}""".format(
        num_threads = config.param('itags_QC', 'num_threads', 'int'),
        qual = config.param('default', 'qual', 'int'),
        primer_mismatch = config.param('itags_QC', 'primerMismatch', 'int')
    )

    return job

def itagsQC(infile, revPrimer, fwdPrimer, outfile, outfileFailed):
        
    job = Job(
        [infile],
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
memtime  itagsQC.pl
  --infile {infile} \\""".format(
     infile = infile
    )

    if(revPrimer != "null"):
        job.command+="""\
  --primer_3_prime {revPrimer} \\
  --length_3_prime {length_3_prime} \\""".format(
            revPrimer = revPrimer,
            length_3_prime = config.param('itags_QC', 'length3Prime', 'int') 
        )

    if(fwdPrimer != "null"):
        job.command+="""\
  --primer_5_prime {fwdPrimer} \\
  --length_5_prime {length_5_prime} \\""".format(
            fwdPrimer = fwdPrimer,
            length_5_prime = config.param( 'itags_QC', 'length5Prime', 'int')
        )
    job.command+="""\
  --qscore_1 {qscore1} \\
  --qscore_2 {qscore2} \\
  --outfile {outfile}
  --outfile_failed {outfileFailed}
  --num_threads {num_threads} \\
  --qual {qual} \\
  --lq_threshold {lq_threshold} \\
  --primer_mismatch {primer_mismatch} \\
  --min_length {min_length} \\
  --N {N}""".format(
        qscore_1 = config.param('itags_QC', 'qscore1', 'int'),
        qscore_2 = config.param('itags_QC', 'qscore2', 'int'),
        outfile = outfile,
        outfile_failed = outfile_failed,
        num_threads =  config.param('itags_QC', 'num_threads', 'int'),
        qual = config.param('default', 'qual', 'int'),
        lq_threshold = config.param('itags_QC', 'lq_threshold', 'int'),
        primer_mismatch = config.param( 'itags_QC', 'primerMismatch', 'float'),
        min_length = config.param('itags_QC', 'minlength', 'int'),
        N = config.param('itags_QC', 'N', 'int')
    )
                
    return job

def countReport(rA_files, rA_names, analysisType, barcodesDist, OTUtable, obsTable, outfile):
        
    job = Job(
        [OTUtable],
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    cmd = "memtime countReport.pl \\"

    for file in rA_files:
        cmd += " --file " + file

    for name in rA_names:
        cmd += " --name " + name
        cmd += " --analysisType " + analysisType
        cmd += " --barcodesDist " + barcodesDist
        cmd += " --OTUtable " + OTUtable
        cmd += " --obsTable " + obsTable
        cmd += " > " + outfile
     
    job.command = cmd
                    
    return job

def txtToPdf(infile, outfile):
                
    job = Job(
        [infile],
        [outfile],
        [
            ['memtime', 'module_memtime'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
memtime txtToPdf.pl
  --infile {infile}
  --outfile {outfile}""".format(
        infile = infile,
        outfile = outfile
    )
    return job

def mergePdf(command):
        
    dummyOutfile = "mergepdf.mugqic.done";

    job = Job(
        [""],
        [dummyOutfile],
        [
            ['memtime', 'module_memtime'],
            ['ghostscript', 'module_ghostscript'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    
    job.command="""\
memtime {command} \\
&& touch {dummyOutfile}; """.format(
        command = command,
        dummyOutfile = dummyOutfile
    )

    return job

def clustering1(infile, barcodes, outdir):
        
    job = Job(
        [infile],
        ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"],
        [
            ['memtime', 'module_memtime'],
            ['usearch', 'module_usearch'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
memtime clustering1.pl \\
  --infile_fastq {infile} \\
  --ref_db {ref_db} \\
  --barcodes {barcodes} \\
  --outdir {outdir} \\
  --num_threads {num_threads}.""".format(
    infile = infile,
    barcodes = barcodes,
    ref_db =  config.param( 'DB', 'chimeras', 'filepath'),
    outdir = outdir,
    num_threads = config.param( 'clustering', 'num_threads', 1, 'int')
    )

    return job

def clustering2(infile, barcodes, outdir):
        
    job = Job(
        [infile],
        ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"],
        [
            ['memtime', 'module_memtime'],
            ['usearch', 'module_usearch'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )
    job.command="""\
memtime clustering2.pl \\
  --infile_fastq {infile} \\
  --ref_db {ref_db} \\ 
  --barcodes {barcodes} \\
  --outdir {outdir} \\
  --num_threads {num_threads}""".format(
        infile = infile,
        ref_db = config.param( 'DB', 'chimeras', 1, 'path'),
        barcodes = barcodes,
        outdir = outdir,
        num_threads = config.param( 'clustering', 'num_threads', 'int')
    )
    return job

def clustering3(infile, barcodes, outdir):

    job = Job(
        [infile],
        ["outdir/obs_filtered.fasta", "outdir/obs_filtered.tsv"],
        [
            ['memtime', 'module_memtime'],
            ['usearch', 'module_usearch'],
            ['dnaclust', 'module_dnaclust'],
            ['tools', 'module_mugqic_tools'],
            ['perl', 'module_perl']
        ]
    )

    job.command="""\
memtime clustering3.pl \\
  --infile_fastq {infile} \\
  --ref_db {ref_db} \\ 
  --barcodes {barcodes} \\
  --outdir {outdir} \\
  --lowAbunCutOff {lowAbunCutoff} \\
  --num_threads {num_threads}""".format(
    infile = infile,
    ref_db =  config.param('DB', 'chimeras', 'filepath'),
    barcodes = barcodes,
    outdir = outdir,
    lowAbunCutoff = config.param('clustering', 'lowAbunCutOff', 'int'),
    num_threads =  config.param('clustering', 'num_threads', 'int')
    )
    
    return job

def clientReport(iniFilePath, projectPath, pipelineType, reportPath):
        
    pipeline = "pipeline=\"" + pipelineType + "\""
    titleTMP = config.param('report', 'projectName')
    title = "report.title=\"" + titleTMP + "\""
    authorTMP = config.param('report', 'report.author')
    author = "report.author=\"" + authorTMP + "\""
    contactTMP = config.param('report', 'report.contact')
    contact = "report.contact=\"" + contactTMP + "\""

    job = Job(
        [iniFilePath],
        [projectPath],
        [
            ['R', 'module_R'],
            ['R', 'module_mugqic_R_packages']
        ]
    )

    job.command="""\
R --no-save -e 'library(gqSeqUtils) \\
mugqicPipelineReport( \\
  pipeline="{pipelineType}", \\
  report.path="{reportPath}", \\
  ini.file.path="{iniFilePath}", \\
  report.title="{project_name}", \\
  report.author="{project_author}", \\
  report.contact="{report_contact}", \\
  project.path="{projectPath}" """.format(
        pipelineType = pipelineType,
        reportPath = reportPath,
        iniFilePath = iniFilePath,
        projectName = titleTMP,
        project_author = authorTMP,
        report_contact = contactTMP,
        projectPath = projectPath
    )
 
    return job

def cleanup(tmpdir):
        
    job = Job(
        [""], 
        [""],
        [
            ['memtime', 'module_memtime']
        ]
    )
    
    job.command="""\
memtime \\
rm  {tmpdir} -rf""".format(
    tmpdir = tmpdir
    )
    return job

def templateSub(outdir):
        
    job = Job(
        ["undef"],
        ["undef"],
        [
            ['memtime', 'module_memtime']
        ]
    )
    job.command="memtime"

    return job
