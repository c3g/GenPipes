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

def assembly_stats(
    short_reads,
    long_reads,
    corrected_reads,
    filtered_summary,
    contigs,
    sample_name,
    suffix,
    estimated_genome_size,
    smrt_cells,
    outdir
    ):

    return Job(
        [short_reads, long_reads, corrected_reads, filtered_summary, contigs],
        [
            os.path.join(outdir, "pacBioGraph_readLengthScore.pdf"),
            os.path.join(outdir, "pacBioGraph_readLengthScore.jpeg"),
            os.path.join(outdir, "pacBioGraph_histoReadLength.pdf"),
            os.path.join(outdir, "pacBioGraph_histoReadLength.jpeg"),
            os.path.join(outdir, "summaryTableAssembly.tsv"),
            os.path.join(outdir, "summaryTableReads.tsv"),
            os.path.join(outdir, "summaryTableReads2.tsv")
        ],
        [
            ['pacbio_tools_assembly_stats', 'module_perl'],
            ['pacbio_tools_assembly_stats', 'module_R'],
            ['pacbio_tools_assembly_stats', 'module_mugqic_tools']
        ],
        command="""\
pacBioAssemblyStats.pl \\
  --shortReads {short_reads} \\
  --longReads {long_reads} \\
  --correctedReads {corrected_reads} \\
  --filteredSummary {filtered_summary} \\
  --contigs {contigs} \\
  --sampleName {sample_name} \\
  --suffix {suffix} \\
  --estimatedGenomeSize {estimated_genome_size} \\
  --smrtCells {smrt_cells} \\
  --outdir {outdir}""".format(
        short_reads=short_reads,
        long_reads=long_reads,
        corrected_reads=corrected_reads,
        filtered_summary=filtered_summary,
        contigs=contigs,
        sample_name=sample_name,
        suffix=suffix,
        estimated_genome_size=estimated_genome_size,
        smrt_cells=smrt_cells,
        outdir=outdir
    ))

def celera_config(
    mer_size,
    infile,
    outfile
    ):

    return Job(
        [infile],
        [outfile],
        [
            ['pacbio_tools_celera_config', 'module_perl'],
            ['pacbio_tools_celera_config', 'module_mugqic_tools']
        ],
        command="""\
pacBioAssemblyCeleraConfig.pl \\
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
        meryl_threads=config.param('pacbio_tools_celera_config', 'meryl_threads', param_type='int'),
        ovl_threads=config.param('pacbio_tools_celera_config', 'ovl_threads', param_type='int'),
        overlapper=config.param('pacbio_tools_celera_config', 'overlapper'),
        mer_compression=config.param('pacbio_tools_celera_config', 'mer_compression'),
        mer_size=mer_size,
        meryl_memory=config.param('pacbio_tools_celera_config', 'meryl_memory'),
        ovl_error_rate=config.param('pacbio_tools_celera_config', 'ovl_error_rate', param_type='float'),
        ovl_min_len=config.param('pacbio_tools_celera_config', 'ovl_min_len', param_type='int'),
        frg_min_len=config.param('pacbio_tools_celera_config', 'frg_min_len', param_type='int'),
        ovl_store_memory=config.param('pacbio_tools_celera_config', 'ovl_store_memory', param_type='int'),
        ovl_concurrency=config.param('pacbio_tools_celera_config', 'ovl_concurrency'),
        ovl_corr_concurrency=config.param('pacbio_tools_celera_config', 'ovl_corr_concurrency', param_type='int'),
        cns_concurrency=config.param('pacbio_tools_celera_config', 'cns_concurrency', param_type='int'),
        frg_corr_threads=config.param('pacbio_tools_celera_config', 'frg_corr_threads', param_type='int'),
        stop_after=config.param('pacbio_tools_celera_config', 'stop_after'),
        unitigger=config.param('pacbio_tools_celera_config', 'unitigger'),
        outfile=outfile
    ))

def compile(
    indir,
    sample_name,
    estimated_genome_size,
    outfile
    ):

    return Job(
        # Input files need to be specified in the wrapper since they depend on cutoffs, mer sizes and polishing rounds
        [None],
        [outfile],
        [
            ['pacbio_tools_get_cutoff', 'module_perl'],
            ['pacbio_tools_get_cutoff', 'module_mugqic_tools']
        ],
        command="""\
pacBioCompileStats.pl \\
  --indir {indir} \\
  --estimatedGenomeSize {estimated_genome_size} \\
  --sampleName {sample_name} \\
  > {outfile}""".format(
        indir=indir,
        estimated_genome_size=estimated_genome_size,
        sample_name=sample_name,
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
            ['pacbio_tools_get_cutoff', 'module_perl'],
            ['pacbio_tools_get_cutoff', 'module_mugqic_tools']
        ],
        command="""\
pacBioGetCutoff.pl \\
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
            ['pacbio_tools_get_cutoff', 'module_perl'],
            ['pacbio_tools_get_cutoff', 'module_mugqic_tools']
        ],
        command="""\
pacBioSplitReads.pl \\
  --infile {infile} \\
  --cutoff `cat {cutoff}` \\
  --outfileShort {short_reads} \\
  --outfileLong {long_reads}""".format(
        infile=infile,
        cutoff=cutoff,
        short_reads=short_reads,
        long_reads=long_reads
    ))
