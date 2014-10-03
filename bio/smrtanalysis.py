#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def blasr(
    infile,
    infile_long,
    outfile,
    outfile_fofn,
    sam=False
    ):

    outfile_filtered = outfile + ".filtered"

    return  Job(
        [infile, infile_long],
        [outfile_filtered, outfile_fofn],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ], 
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime blasr \\
  {infile} \\
  {infile_long} \\
  -out {outfile} \\
  -m {m} \\
  -nproc {threads} \\
  -bestn {bestn} \\
  -nCandidates {n_candidates} \\
  -noSplitSubreads \\
  -minReadLength {min_read_length} \\
  -maxScore {max_score} \\
  -maxLCPLength {max_lcp_length}{sam} && \\
echo {outfile} > {outfile_fofn} &&
filterm4.py {outfile} > {outfile_filtered} 2> {outfile_filtered}.log""".format(
        infile=infile,
        infile_long=infile_long,
        outfile=outfile,
        m=config.param('smrtanalysis_blasr', 'm', type='int'),
        threads=config.param('smrtanalysis_blasr', 'threads', type='posint'),
        bestn=config.param('smrtanalysis_blasr', 'bestn', type='int'),
        n_candidates=config.param('smrtanalysis_blasr', 'n_candidates', type='int'),
        min_read_length=config.param('smrtanalysis_blasr', 'min_read_length', type='int'),
        max_score=config.param('smrtanalysis_blasr', 'max_score', type='int'),
        max_lcp_length=config.param('smrtanalysis_blasr', 'max_lcp_length', type='int'),
        sam=" \\\n  -sam" if sam else "",
        outfile_fofn=outfile_fofn,
        outfile_filtered=outfile_filtered
    ))

def fastq_to_ca(
    libraryname,
    reads,
    outfile
    ):

    return  Job(
        [reads],
        [outfile],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ],
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime fastqToCA \\
  -technology sanger \\
  -type sanger \\
  -libraryname {libraryname} \\
  -reads {reads} \\
  > {outfile}""".format(
        libraryname=libraryname,
        reads=reads,
        outfile=outfile
    ))

def filtering(
    fofn,
    input_xml,
    params_xml,
    output_dir,
    log
    ):

    ref_params_xml=config.param('smrtanalysis_filtering', 'filtering_settings')
    output_prefix = os.path.join(output_dir, "data", "filtered_subreads.")
    output = output_prefix + "fastq"

    return Job(
        [fofn, ref_params_xml],
        [output, output_prefix + "fasta"],
        [
            ['smrtanalysis_filtering', 'module_memtime'],
            ['smrtanalysis_filtering', 'module_prinseq'],
            ['smrtanalysis_filtering', 'module_smrtanalysis']
        ],
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime fofnToSmrtpipeInput.py {fofn} > {input_xml} && \\
cp {fofn} {output_dir}/input.fofn && \\
memtime sed -e 's/MINSUBREADLENGTH/{min_subread_length}/g' -e 's/MINREADLENGTH/{min_read_length}/g' -e 's/MINQUAL/{min_qual}/g' \\
  < {ref_params_xml} > {params_xml} && \\
memtime smrtpipe.py \\
  -D NPROC={threads} \\
  -D TMP={tmp_dir} \\
  --params={params_xml} \\
  --output={output_dir} \\
  --debug \\
  xml:{input_xml} \\
  > {log} && \\
memtime prinseq-lite.pl \\
  -verbose \\
  -fastq {output} \\
  -out_format 1 \\
  -out_good {output_dir}/data/filtered_subreads""".format(
        fofn=fofn,
        input_xml=input_xml,
        min_subread_length=config.param('smrtanalysis_filtering', 'min_subread_length'),
        min_read_length=config.param('smrtanalysis_filtering', 'min_read_length'),
        min_qual=config.param('smrtanalysis_filtering', 'min_qual'),
        ref_params_xml=ref_params_xml,
        params_xml=params_xml,
        threads=config.param('smrtanalysis_filtering', 'threads'),
        tmp_dir=config.param('smrtanalysis_filtering', 'tmp_dir', type='dirpath'),
        output_dir=output_dir,
        log=log,
        output=output
    ))

def m4topre(
    infile,
    allm4,
    subreads,
    outfile
    ):

    outfile_filtered = outfile + ".filtered"

    return  Job(
        [infile],
        [outfile],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ], 
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime m4topre.py \\
  {infile} \\
  {allm4} \\
  {subreads} \\
  {bestn} \\
  > {outfile}""".format(
        infile=infile,
        allm4=allm4,
        subreads=subreads,
        bestn=config.param('smrtanalysis_m4topre', 'bestn', type='int'),
        outfile=outfile
    ))

def pbdagcon(
    infile,
    outfile,
    outfile_fastq
    ):

    outfile_filtered = outfile + ".filtered"

    return  Job(
        [infile],
        [outfile, outfile_fastq],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ], 
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime pbdagcon -a -j {threads} \\
  {infile} \\
  > {outfile} && \\
awk '{{if ($0~/>/) {{sub(/>/,"@",$0);print;}} else {{l=length($0);q=""; while (l--) {{q=q "9"}}printf("%s\\n+\\n%s\\n",$0,q)}}}}' {outfile} \\
  > {outfile_fastq}""".format(
        threads=config.param('smrtanalysis_pbdagcon', 'threads', type='posint'),
        infile=infile,
        outfile=outfile,
        outfile_fastq=outfile_fastq
    ))

def pbutgcns(
    gpk_store,
    tig_store,
    unitigs_list,
    prefix,
    outdir,
    outfile,
    tmp_dir
    ):

    return  Job(
        [gpk_store, tig_store],
        [outfile],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ],
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime tigStore \\
  -g {gpk_store} \\
  -t {tig_store} 1 \\
  -d properties -U | \\
awk 'BEGIN{{t=0}}$1=="numFrags"{{if ($2 > 1) {{print t, $2}} t++}}' | sort -nrk2,2 \\
  > {unitigs_list} && \\
mkdir -p {outdir} {tmp_dir} && \\
tmp={tmp_dir} && \\
cap={prefix} && \\
utg={unitigs_list} && \\
nprocs={threads} && \\
cns={outfile} && \\
pbutgcns_wf.sh""".format(
        gpk_store=gpk_store,
        tig_store=tig_store,
        unitigs_list=unitigs_list,
        outdir=outdir,
        tmp_dir=tmp_dir,
        prefix=prefix,
        threads=config.param('smrtanalysis_pbutgcns', 'threads', type='posint'),
        outfile=outfile
    ))

def run_ca(
    infile,
    ini,
    prefix,
    outdir
    ):

    return  Job(
        [infile, ini],
        [
            os.path.join(outdir, prefix + ".ovlStore.list"),
            os.path.join(outdir, prefix + ".tigStore"),
            os.path.join(outdir, prefix + ".gkpStore"),
        ],
        [
            ['smrtanalysis_blasr', 'module_memtime'],
            ['smrtanalysis_blasr', 'module_smrtanalysis']
        ],
        command = """\
set +u && source $SEYMOUR_HOME/etc/setup.sh && set -u && \\
memtime runCA \\
  -s {ini} \\
  -p {prefix} \\
  -d {outdir} \\
  {infile}""".format(
        infile=infile,
        ini=ini,
        prefix=prefix,
        outdir=outdir
    ))
