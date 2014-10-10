#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def celera_config(
    mer_size,
    infile,
    outfile
    ):

    return Job(
        [infile],
        [outfile],
        [
            ['pacbio_tools_celera_config', 'module_memtime'],
            ['pacbio_tools_celera_config', 'module_perl'],
            ['pacbio_tools_celera_config', 'module_mugqic_tools']
        ],
        command = """\
memtime pacBioAssemblyCeleraConfig.pl \\
  --infile {infile} \\
  --merylThreads {meryl_threads} \\
  --ovlThreads {ovl_threads} \\
  --overlapper {overlapper} \\
  --merCompression {mer_compression} \\
  --merSize {mer_size} \\
  --merylMemory {meryl_memory} \\
  --ovlErrorRate {ovl_error_rate} \\
  --ovlMinLen {ovl_min_len} \\
  --frgMinLen {frg_min_len} \\
  --ovlStoreMemory {ovl_store_memory} \\
  --ovlConcurrency {ovl_concurrency} \\
  --ovlCorrConcurrency {ovl_corr_concurrency} \\
  --cnsConcurrency {cns_concurrency} \\
  --frgCorrThreads {frg_corr_threads} \\
  --stopAfter {stop_after} \\
  --unitigger {unitigger} \\
  > {outfile}""".format(
        infile=infile,
        meryl_threads=config.param('pacbio_tools_celera_config', 'meryl_threads', type='int'),
        ovl_threads=config.param('pacbio_tools_celera_config', 'ovl_threads', type='int'),
        overlapper=config.param('pacbio_tools_celera_config', 'overlapper'),
        mer_compression=config.param('pacbio_tools_celera_config', 'mer_compression'),
        mer_size=mer_size,
        meryl_memory=config.param('pacbio_tools_celera_config', 'meryl_memory'),
        ovl_error_rate=config.param('pacbio_tools_celera_config', 'ovl_error_rate', type='float'),
        ovl_min_len=config.param('pacbio_tools_celera_config', 'ovl_min_len', type='int'),
        frg_min_len=config.param('pacbio_tools_celera_config', 'frg_min_len', type='int'),
        ovl_store_memory=config.param('pacbio_tools_celera_config', 'ovl_store_memory', type='int'),
        ovl_concurrency=config.param('pacbio_tools_celera_config', 'ovl_concurrency'),
        ovl_corr_concurrency=config.param('pacbio_tools_celera_config', 'ovl_corr_concurrency', type='int'),
        cns_concurrency=config.param('pacbio_tools_celera_config', 'cns_concurrency', type='int'),
        frg_corr_threads=config.param('pacbio_tools_celera_config', 'frg_corr_threads', type='int'),
        stop_after=config.param('pacbio_tools_celera_config', 'stop_after'),
        unitigger=config.param('pacbio_tools_celera_config', 'unitigger'),
        outfile=outfile
    ))

def get_cutoff(
    infile,
    coverage,
    genome_size,
    coverage_cutoff,
    outfile
    ):

    return Job(
        [infile],
        [outfile],
        [
            ['pacbio_tools_get_cutoff', 'module_memtime'],
            ['pacbio_tools_get_cutoff', 'module_perl'],
            ['pacbio_tools_get_cutoff', 'module_mugqic_tools']
        ],
        command = """\
memtime pacBioGetCutoff.pl \\
  --infile {infile} \\
  --coverage {coverage} \\
  --genomeSize {genome_size} \\
  --coverageCutoff {coverage_cutoff} \\
  > {outfile}""".format(
        infile=infile,
        coverage=coverage,
        genome_size=genome_size,
        coverage_cutoff=coverage_cutoff,
        outfile=outfile
    ))

def split_reads(
    infile,
    cutoff, # a file containing the cutoff number is actually passed in arg here.
    short_reads,
    long_reads
    ):

    return Job(
        [infile, cutoff],
        [short_reads, long_reads],
        [
            ['pacbio_tools_get_cutoff', 'module_memtime'],
            ['pacbio_tools_get_cutoff', 'module_perl'],
            ['pacbio_tools_get_cutoff', 'module_mugqic_tools']
        ],
        command = """\
memtime pacBioSplitReads.pl \\
  --infile {infile} \\
  --cutoff `cat {cutoff}` \\
  --outfileShort {short_reads} \\
  --outfileLong {long_reads}""".format(
        infile=infile,
        cutoff=cutoff,
        short_reads=short_reads,
        long_reads=long_reads
    ))
