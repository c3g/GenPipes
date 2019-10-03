#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import os
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
from bfx.readset import parse_pacbio_readset_file

from bfx import blast
from bfx import circlator
from bfx import gq_seq_utils
from bfx import mummer
from bfx import pacbio_tools
from bfx import smrtanalysis
from pipelines import common

log = logging.getLogger(__name__)

class PacBioAssembly(common.MUGQICPipeline):
    """
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
    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
        super(PacBioAssembly, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_pacbio_readset_file(self.args.readsets.name)
            else:
                self.argparser.error("argument -r/--readsets is required!")

        return self._readsets

    def smrtanalysis_filtering(self):
        """
        Filter reads and subreads based on their length and QVs, using smrtpipe.py (from the SmrtAnalysis package).

        1. fofnToSmrtpipeInput.py
        2. modify RS_Filtering.xml files according to reads filtering values entered in .ini file
        3. smrtpipe.py with filtering protocol
        4. prinseq-lite.pl: write fasta file based on fastq file

        Informative run metrics such as loading efficiency, read lengths and base quality are generated in this step as well.
        """

        jobs = []
        jobs.append(
            Job(
                [
                    os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), config.param('smrtanalysis_filtering', 'celera_settings')), os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), config.param('smrtanalysis_filtering', 'filtering_settings'))
                ],
                [config.param('smrtanalysis_filtering', 'celera_settings'), config.param('smrtanalysis_filtering', 'filtering_settings')],
                command="cp -a -f " + os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "protocols") + " .; touch protocols/*",
                name="smrtanalysis_filtering.config"
            )
        )

        for sample in self.samples:
            fofn = os.path.join("fofns", sample.name + ".fofn")
            input_files = [config.param('smrtanalysis_filtering', 'celera_settings'),
                           config.param('smrtanalysis_filtering', 'filtering_settings')]
            for readset in sample.readsets:
                if readset.bax_files:
                    # New PacBio format is BAX
                    input_files.extend(readset.bax_files)
                else:
                    # But old PacBio format BAS should still be supported
                    input_files.extend(readset.bas_files)
            filtering_directory = os.path.join(sample.name, "filtering")

            jobs.append(concat_jobs([
                Job(command="mkdir -p fofns", samples=[sample]),
                Job(input_files, [fofn], command="""\
`cat > {fofn} << END
{input_files}
END
`""".format(input_files="\n".join(input_files), fofn=fofn)),
                Job(command="mkdir -p " + filtering_directory),
                smrtanalysis.filtering(
                    fofn,
                    os.path.join(filtering_directory, "input.xml"),
                    os.path.join(sample.name, "filtering.xml"),
                    filtering_directory,
                    os.path.join(filtering_directory, "smrtpipe.log")
                )
            ], name="smrtanalysis_filtering." + sample.name))

        return jobs

    def pacbio_tools_get_cutoff(self):
        """
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
        """

        jobs = []

        for sample in self.samples:
            log.info("Sample: " + sample.name)
            sample_nb_base_pairs = sum([readset.nb_base_pairs for readset in sample.readsets])
            log.info("nb_base_pairs: " + str(sample_nb_base_pairs))
            estimated_genome_size = sample.readsets[0].estimated_genome_size
            log.info("estimated_genome_size: " + str(estimated_genome_size))
            estimated_coverage = sample_nb_base_pairs / estimated_genome_size
            log.info("estimated_coverage: " + str(estimated_coverage))

            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)

                log.info("COVERAGE_CUTOFF: " + coverage_cutoff + "_X_coverage")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + os.path.join(coverage_directory, "preassembly"), samples=[sample]),
                    pacbio_tools.get_cutoff(
                        os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                        estimated_coverage,
                        estimated_genome_size,
                        coverage_cutoff,
                        os.path.join(coverage_directory, "preassemblyMinReadSize.txt")
                    )
                ], name="pacbio_tools_get_cutoff." + sample.name + ".coverage_cutoff" + cutoff_x))

        return jobs

    def preassembly(self):
        """
        Having in hand a cutoff value, filtered reads are splitted between short and long reads. Short reads
        are aligned against long reads and consensus (e.g. corrected reads) are generated from these alignments.

        1. split reads between long and short
        2. blasr: aligner for PacBio reads
        3. m4topre: convert .m4 blasr output in .pre format
        4. pbdagcon (aka HGAP2): generate corrected reads from alignments
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")
                job_name_suffix = sample.name + ".coverage_cutoff" + cutoff_x

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + preassembly_directory, samples=[sample]),
                    pacbio_tools.split_reads(
                        os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                        os.path.join(coverage_directory, "preassemblyMinReadSize.txt"),
                        os.path.join(preassembly_directory, "filtered_shortreads.fa"),
                        os.path.join(preassembly_directory, "filtered_longreads.fa")
                    )
                ], name="pacbio_tools_split_reads." + job_name_suffix))

                job = smrtanalysis.blasr(
                    os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                    os.path.join(preassembly_directory, "filtered_longreads.fa"),
                    os.path.join(preassembly_directory, "seeds.m4"),
                    os.path.join(preassembly_directory, "seeds.m4.fofn")
                )
                job.name = "smrtanalysis_blasr." + job_name_suffix
                job.samples = [sample]
                jobs.append(job)

                job = smrtanalysis.m4topre(
                    os.path.join(preassembly_directory, "seeds.m4.filtered"),
                    os.path.join(preassembly_directory, "seeds.m4.fofn"),
                    os.path.join(sample.name, "filtering", "data", "filtered_subreads.fasta"),
                    os.path.join(preassembly_directory, "aln.pre")
                )
                job.name = "smrtanalysis_m4topre." + job_name_suffix
                job.samples = [sample]
                jobs.append(job)

                job = smrtanalysis.pbdagcon(
                    os.path.join(preassembly_directory, "aln.pre"),
                    os.path.join(preassembly_directory, "corrected.fasta"),
                    os.path.join(preassembly_directory, "corrected.fastq")
                )
                job.name = "smrtanalysis_pbdagcon." + job_name_suffix
                job.samples = [sample]
                jobs.append(job)

        return jobs

    def assembly(self):
        """
        Corrected reads are assembled to generates contigs. Please see the
        [Celera documentation](http://wgs-assembler.sourceforge.net/wiki/index.php?title=RunCA).
        Quality of assembly seems to be highly sensitive to parameters you give Celera.

        1. generate celera config files using parameters provided in the .ini file
        2. fastqToCA: generate input file compatible with the Celera assembler
        3. runCA: run the Celera assembler
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    sample_cutoff_mer_size = "_".join([sample.name, cutoff_x, mer_size_text])
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    assembly_directory = os.path.join(mer_size_directory, "assembly")

                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + assembly_directory, samples=[sample]),
                        pacbio_tools.celera_config(
                            mer_size,
                            config.param('DEFAULT', 'celera_settings'),
                            os.path.join(mer_size_directory, "celera_assembly.ini")
                        )
                    ], name="pacbio_tools_celera_config." + sample_cutoff_mer_size))

                    job = smrtanalysis.fastq_to_ca(
                        sample_cutoff_mer_size,
                        os.path.join(preassembly_directory, "corrected.fastq"),
                        os.path.join(preassembly_directory, "corrected.frg")
                    )
                    job.name = "smrtanalysis_fastq_to_ca." + sample_cutoff_mer_size
                    job.samples = [sample]
                    jobs.append(job)

                    jobs.append(concat_jobs([
                        Job(command="rm -rf " + assembly_directory, samples=[sample]),
                        smrtanalysis.run_ca(
                            os.path.join(preassembly_directory, "corrected.frg"),
                            os.path.join(mer_size_directory, "celera_assembly.ini"),
                            sample_cutoff_mer_size,
                            assembly_directory
                        )
                    ], name="smrtanalysis_run_ca." + sample_cutoff_mer_size))

                    job = smrtanalysis.pbutgcns(
                        assembly_directory, 
                        sample_cutoff_mer_size,
                        mer_size_directory
                    )
                    job.name = "smrtanalysis_pbutgcns." + sample_cutoff_mer_size
                    job.samples = [sample]
                    jobs.append(job)

        return jobs

    def polishing(self):
        """
        Align raw reads on the Celera assembly with BLASR. Load pulse information from bax or bas files into aligned file. Sort that file and run quiver (variantCaller.py).

        1. generate fofn
        2. upload Celera assembly with smrtpipe refUploader
        3. compare sequences
        4. load pulses
        5. sort .cmp.h5 file
        6. variantCaller.py
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    assembly_directory = os.path.join(mer_size_directory, "assembly")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))

                    for polishing_round in range(1, polishing_rounds + 1):
                        polishing_round_directory = os.path.join(mer_size_directory, "polishing" + str(polishing_round))
                        # smrtanalysis.reference_uploader transforms "-" into "_" in fasta filename
                        sample_cutoff_mer_size_polishing_round = "_".join([sample.name.replace("-", "_"), cutoff_x, mer_size_text, "polishingRound" + str(polishing_round)])
                        job_name_suffix = "_".join([sample.name, cutoff_x, mer_size_text, "polishingRound" + str(polishing_round)])

                        if polishing_round == 1:
                            fasta_file = os.path.join(assembly_directory, "9-terminator", "_".join([sample.name, cutoff_x, mer_size_text]) + ".ctg.fasta")
                        else:
                            fasta_file = os.path.join(mer_size_directory, "polishing" + str(polishing_round - 1), "data", "consensus.fasta")

                        jobs.append(concat_jobs([
                            Job(command="mkdir -p " + os.path.join(polishing_round_directory, "data"), samples=[sample]),
                            smrtanalysis.reference_uploader(
                                polishing_round_directory,
                                sample_cutoff_mer_size_polishing_round,
                                fasta_file
                            )
                        ], name="smrtanalysis_reference_uploader." + job_name_suffix))

                        job = smrtanalysis.pbalign(
                            polishing_round_directory,
                            sample.name,
                            sample_cutoff_mer_size_polishing_round
                        )
                        job.name = "smrtanalysis_pbalign." + job_name_suffix
                        job.samples = [sample]
                        jobs.append(job)

                        jobs.append(concat_jobs([
                            Job(samples=[sample]),
                            smrtanalysis.load_chemistry(
                                os.path.join(polishing_round_directory, "data", "aligned_reads.cmp.h5"),
                                os.path.join(sample.name, "filtering", "input.fofn"),
                                os.path.join(polishing_round_directory, "data", "aligned_reads.loadChemistry.cmp.h5")
                            ),
                            smrtanalysis.load_pulses(
                                os.path.join(polishing_round_directory, "data", "aligned_reads.loadChemistry.cmp.h5"),
                                os.path.join(sample.name, "filtering", "input.fofn"),
                                os.path.join(polishing_round_directory, "data", "aligned_reads.loadPulses.cmp.h5")
                            )
                        ], name = "smrtanalysis_load_chemistry_load_pulses." + job_name_suffix))

                        job = smrtanalysis.cmph5tools_sort(
                            os.path.join(polishing_round_directory, "data", "aligned_reads.loadPulses.cmp.h5"),
                            os.path.join(polishing_round_directory, "data", "aligned_reads.sorted.cmp.h5")
                        )
                        job.name = "smrtanalysis_cmph5tools_sort." + job_name_suffix
                        job.samples=[sample]
                        jobs.append(job)

                        job = smrtanalysis.variant_caller(
                            os.path.join(polishing_round_directory, "data", "aligned_reads.sorted.cmp.h5"),
                            os.path.join(polishing_round_directory, sample_cutoff_mer_size_polishing_round, "sequence", sample_cutoff_mer_size_polishing_round + ".fasta"),
                            os.path.join(polishing_round_directory, "data", "variants.gff"),
                            os.path.join(polishing_round_directory, "data", "consensus.fasta.gz"),
                            os.path.join(polishing_round_directory, "data", "consensus.fastq.gz")
                        )
                        job.name = "smrtanalysis_variant_caller." + job_name_suffix
                        job.samples=[sample]
                        jobs.append(job)

                        job = smrtanalysis.summarize_polishing(
                            "_".join([sample.name, cutoff_x, mer_size_text]),
                            os.path.join(polishing_round_directory, sample_cutoff_mer_size_polishing_round),
                            os.path.join(polishing_round_directory, "data", "aligned_reads.sorted.cmp.h5"),
                            os.path.join(polishing_round_directory, "data", "alignment_summary.gff"),
                            os.path.join(polishing_round_directory, "data", "coverage.bed"),
                            os.path.join(sample.name, "filtering", "input.fofn"),
                            os.path.join(polishing_round_directory, "data", "aligned_reads.sam"),
                            os.path.join(polishing_round_directory, "data", "variants.gff"),
                            os.path.join(polishing_round_directory, "data", "variants.bed"),
                            os.path.join(polishing_round_directory, "data", "variants.vcf")
                        )
                        job.name = "smrtanalysis_summarize_polishing." + job_name_suffix
                        job.samples=[sample]
                        jobs.append(job)

        return jobs

    def pacbio_tools_assembly_stats(self):

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    sample_cutoff_mer_size = "_".join([sample.name, cutoff_x, mer_size_text])
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    blast_directory = os.path.join(mer_size_directory, "blast")
                    mummer_file_prefix = os.path.join(mer_size_directory, "mummer", sample.name + ".")
                    report_directory = os.path.join(mer_size_directory, "report")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "consensus.fasta")
                    # Number of unique run-smartcells per sample
                    smartcells = len(set([(readset.run, readset.smartcell) for readset in sample.readsets]))

                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + report_directory, samples=[sample]),
                        # Generate table(s) and figures
                        pacbio_tools.assembly_stats(
                            os.path.join(preassembly_directory, "filtered_shortreads.fa"),
                            os.path.join(preassembly_directory, "filtered_longreads.fa"),
                            os.path.join(preassembly_directory, "corrected.fasta"),
                            os.path.join(sample.name, "filtering", "data", "filtered_summary.csv"),
                            fasta_consensus,
                            sample.name,
                            cutoff_x + "_" + mer_size,
                            sample.readsets[0].estimated_genome_size,
                            smartcells,
                            report_directory
                        )
                    ], name="pacbio_tools_assembly_stats." + sample_cutoff_mer_size))

                    report_file = os.path.join(report_directory, "PacBioAssembly.pacbio_tools_assembly_stats.md")
                    jobs.append(
                        Job(
                            [
                                os.path.join(report_directory, "summaryTableReads.tsv"),
                                os.path.join(report_directory, "summaryTableReads2.tsv"),
                                os.path.join(report_directory, "summaryTableAssembly.tsv"),
                                os.path.join(report_directory, "pacBioGraph_readLengthScore.pdf"),
                                os.path.join(report_directory, "pacBioGraph_readLengthScore.jpeg"),
                                os.path.join(report_directory, "pacBioGraph_histoReadLength.pdf"),
                                os.path.join(report_directory, "pacBioGraph_histoReadLength.jpeg"),
                                fasta_consensus + ".gz"
                            ],
                            [report_file],
                            [['pacbio_tools_assembly_stats', 'module_pandoc']],
                            command="""\
cp {fasta_consensus}.gz {report_directory}/ && \\
total_subreads=`grep -P '^"Total subreads"\t"' {report_directory}/summaryTableReads.tsv | cut -f2 | sed 's/"//g'` && \\
average_subreads_length=`grep -P '^"Average subread length"\t"' {report_directory}/summaryTableReads.tsv | cut -f2 | sed 's/"//g'` && \\
summary_table_reads2=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, sprintf("%\\47d", $2), sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%\\47d", $5), sprintf("%\\47d", $6), sprintf("%\\47d", $7), sprintf("%\\47d", $8), sprintf("%\\47d", $9), sprintf("%\\47d", $10), sprintf("%\\47d", $11)}}}}' {report_directory}/summaryTableReads2.tsv | sed 's/^#//'` && \\
summary_table_assembly=`awk -F"\t" '{{OFS="\t"; if (NR==1) {{print; gsub(/[^\t]/, "-")}} print}}' {report_directory}/summaryTableAssembly.tsv | sed 's/"//g' | sed 's/\t/|/g'` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable total_subreads="$total_subreads" \\
  --variable average_subreads_length="$average_subreads_length" \\
  --variable summary_table_reads2="$summary_table_reads2" \\
  --variable summary_table_assembly="$summary_table_assembly" \\
  --variable smartcells="{smartcells}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                                fasta_consensus=fasta_consensus,
                                report_directory=report_directory,
                                smartcells=smartcells,
                                report_template_dir=self.report_template_dir,
                                basename_report_file=os.path.basename(report_file),
                                report_file=report_file
                            ),
                            report_files=[os.path.relpath(report_file, mer_size_directory)],
                            name="pacbio_tools_assembly_stats_report." + sample_cutoff_mer_size)
                    )

        return jobs

    def blast(self):
        """
        Blast polished assembly against nr using dc-megablast.
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    blast_directory = os.path.join(mer_size_directory, "blast")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))

                    polishing_round_directory = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds))
                    sample_cutoff_mer_size = "_".join([sample.name, cutoff_x, mer_size_text])
                    blast_report = os.path.join(blast_directory, "blast_report.csv")

                    # Blast contigs against nt
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + blast_directory, samples=[sample]),
                        blast.dcmegablast(
                            os.path.join(polishing_round_directory, "data", "consensus.fasta"),
                            "7",
                            blast_report,
                            os.path.join(polishing_round_directory, "data", "coverage.bed"),
                            blast_directory
                        )
                    ], name="blast_dcmegablast." + sample_cutoff_mer_size))

                    # Get fasta file of best hit.
                    job = blast.blastdbcmd(
                        blast_report,
                        "$(grep -v '^#' < " + blast_report + " | head -n 1 | awk -F '\\t' '{print $2}' | sed 's/gi|\([0-9]*\)|.*/\\1/' | tr '\\n' '  ')",
                        os.path.join(blast_directory, "nt_reference.fasta"),
                    )
                    job.name = "blast_blastdbcmd." + sample_cutoff_mer_size
                    job.samples=[sample]
                    jobs.append(job)

                    report_directory = os.path.join(mer_size_directory, "report")
                    report_file = os.path.join(report_directory, "PacBioAssembly.blast.md")
                    jobs.append(
                        Job(
                            [os.path.join(blast_directory, "blastCov.tsv"), os.path.join(blast_directory, "contigsCoverage.tsv")],
                            [report_file],
                            [['blast', 'module_pandoc']],
                            command="""\
cp {blast_directory}/blastCov.tsv {report_directory}/ && \\
cp {blast_directory}/contigsCoverage.tsv {report_directory}/ && \\
blast_table=`awk -F"\t" '{{OFS="\t"; if (NR==1) {{print; gsub(/[^\t]/, "-")}} print}}' {report_directory}/blastCov.tsv | sed 's/|/\\\\\\\\|/g' | sed 's/\t/|/g' | head -21` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable blast_table="$blast_table" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                                blast_directory=blast_directory,
                                report_directory=report_directory,
                                report_template_dir=self.report_template_dir,
                                basename_report_file=os.path.basename(report_file),
                                report_file=report_file
                            ),
                            report_files=[os.path.relpath(report_file, mer_size_directory)],
                            name="blast_report." + sample_cutoff_mer_size)
                    )

        return jobs

    def mummer(self):
        """
        Using MUMmer, align polished assembly against best hit from blast job. Also align polished assembly against itself to detect structure variation such as repeats, etc.
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    fasta_reference = os.path.join(mer_size_directory, "blast", "nt_reference.fasta")
                    mummer_directory = os.path.join(mer_size_directory, "mummer")
                    mummer_file_prefix = os.path.join(mummer_directory, sample.name + ".")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))

                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "consensus.fasta")
                    sample_cutoff_mer_size = "_".join([sample.name, cutoff_x, mer_size_text])
                    sample_cutoff_mer_size_nucmer = sample_cutoff_mer_size + "-nucmer"

                    # Run nucmer
                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + mummer_directory, samples=[sample]),
                        mummer.reference(
                            mummer_file_prefix + "nucmer",
                            fasta_reference,
                            fasta_consensus,
                            sample_cutoff_mer_size_nucmer,
                            mummer_file_prefix + "nucmer.delta",
                            mummer_file_prefix + "nucmer.delta",
                            mummer_file_prefix + "dnadiff",
                            mummer_file_prefix + "dnadiff.delta",
                            mummer_file_prefix + "dnadiff.delta.snpflank"
                        )
                    ], name="mummer_reference." + sample_cutoff_mer_size_nucmer))

                    jobs.append(concat_jobs([
                        Job(command="mkdir -p " + mummer_directory, samples=[sample]),
                        mummer.self(
                            mummer_file_prefix + "nucmer.self",
                            fasta_consensus,
                            sample_cutoff_mer_size_nucmer + "-self",
                            mummer_file_prefix + "nucmer.self.delta",
                            mummer_file_prefix + "nucmer.self.delta"
                        )
                    ], name="mummer_self." + sample_cutoff_mer_size_nucmer))

                    report_directory = os.path.join(mer_size_directory, "report")
                    report_file = os.path.join(report_directory, "PacBioAssembly.mummer.md")
                    jobs.append(
                        Job(
                            [mummer_file_prefix + "nucmer.self.delta.png", mummer_file_prefix + "nucmer.delta.png"],
                            [report_file],
                            [['mummer', 'module_pandoc']],
                            command="""\
cp {mummer_file_prefix}nucmer.self.delta.png {mummer_file_prefix}nucmer.delta.png {report_directory}/ && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable sample="{sample}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                                mummer_file_prefix=mummer_file_prefix,
                                sample=sample.name,
                                report_directory=report_directory,
                                report_template_dir=self.report_template_dir,
                                basename_report_file=os.path.basename(report_file),
                                report_file=report_file
                            ),
                            report_files=[os.path.relpath(report_file, mer_size_directory)],
                            name="mummer_report." + sample_cutoff_mer_size)
                    )

        return jobs

    def compile(self):
        """
        Compile assembly stats of all conditions used in the pipeline (useful when multiple assemblies are performed).
        """

        jobs = []

        for sample in self.samples:
            # Generate table
            job = pacbio_tools.compile(
                sample.name,
                sample.name,
                sample.readsets[0].estimated_genome_size,
                sample.name + ".compiledStats.csv"
            )
            # Job input files (all consensus.fasta) need to be defined here since only sample directory is given to pacbio_tools.compile
            job.input_files = [os.path.join(
                sample.name,
                coverage_cutoff + "X",
                "merSize" + mer_size,
                "polishing" + str(polishing_round),
                "data",
                "consensus.fasta")
                for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list')
                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list')
                for polishing_round in range(1, config.param('DEFAULT', 'polishing_rounds', type='posint') + 1)
            ]

            job.name = "pacbio_tools_compile." + sample.name
            job.samples=[sample]
            jobs.append(job)

        return jobs

    def circlator(self):
        """
        Circularize the assembly contigs if possible.
        User should launch this step after making sure the quality of the assembly is acceptable.
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    circlator_directory = os.path.join(mer_size_directory, "circlator")
                    circlator_file = os.path.join(circlator_directory, "output")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))

                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "consensus.fasta")

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + circlator_directory, samples=[sample]),
                    circlator.circularize(fasta_consensus, os.path.join(preassembly_directory, "corrected.fastq"), circlator_file)
                ], name = "circlator." + sample.name))

        return jobs

    def basemodification(self):
        """
        Run ipdSummary.py for in silico detection of modified bases
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    basemodification_directory = os.path.join(mer_size_directory, "basemodification")
                    basemodification_file_prefix = os.path.join(basemodification_directory, sample.name + ".basemod")

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))
                        
                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds - 1), "data", "consensus.fasta")
                    aligned_reads = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds), "data", "aligned_reads.sorted.cmp.h5")
                    output_gff = basemodification_file_prefix

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + basemodification_directory),
                    smrtanalysis.basemodification(fasta_consensus, aligned_reads, basemodification_file_prefix, mer_size_directory, polishing_rounds)                 
                ], name = "basemodification." + sample.name))

        return jobs


    def motifMaker(self):
        """
        Run motifMaker to generate motif_summary.csv
        """

        jobs = []

        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)
                preassembly_directory = os.path.join(coverage_directory, "preassembly", "data")

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)
                    basemodification_directory = os.path.join(mer_size_directory, "basemodification")
                    basemodification_file_prefix = os.path.join(basemodification_directory, sample.name + ".basemod")
                    motifMaker_directory = os.path.join(mer_size_directory, "motifMaker")
                    motifMaker_file = os.path.join(motifMaker_directory, (sample.name + ".motif_summary.csv"))
                    output_gff = basemodification_file_prefix + ".gff"

                    polishing_rounds = config.param('DEFAULT', 'polishing_rounds', type='posint')
                    if polishing_rounds > 4:
                        _raise(SanitycheckError("Error: polishing_rounds \"" + str(polishing_rounds) + "\" is invalid (should be between 1 and 4)!"))
                        
                    fasta_consensus = os.path.join(mer_size_directory, "polishing" + str(polishing_rounds - 1), "data", "consensus.fasta")
                    output_gff = basemodification_file_prefix + ".gff"

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + motifMaker_directory),
                    smrtanalysis.motifMaker(fasta_consensus, basemodification_file_prefix, mer_size_directory, polishing_rounds, motifMaker_file)
                ], name = "motifMaker." + sample.name))

        return jobs



    def report_jobs(self):
        """
        Overwrite core pipeline report_jobs method to perform it on every sample/coverage_cutoff/mer_size
        """
        for sample in self.samples:
            for coverage_cutoff in config.param('DEFAULT', 'coverage_cutoff', type='list'):
                cutoff_x = coverage_cutoff + "X"
                coverage_directory = os.path.join(sample.name, cutoff_x)

                for mer_size in config.param('DEFAULT', 'mer_sizes', type='list'):
                    mer_size_text = "merSize" + mer_size
                    mer_size_directory = os.path.join(coverage_directory, mer_size_text)

                    super(PacBioAssembly, self).report_jobs(os.path.join(self.output_dir, mer_size_directory))

    @property
    def steps(self):
        return [
            self.smrtanalysis_filtering,
            self.pacbio_tools_get_cutoff,
            self.preassembly,
            self.assembly,
            self.polishing,
            self.pacbio_tools_assembly_stats,
            self.blast,
            self.mummer,
            self.compile,
            self.circlator,
            self.basemodification,
            self.motifMaker
        ]

if __name__ == '__main__':
    PacBioAssembly()
