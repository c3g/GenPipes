#!/usr/bin/env python

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
import argparse
import logging
import os
import re
import socket
import string
import sys
import collections

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

# MUGQIC Modules
from core.config import global_config_parser, _raise, SanitycheckError
from core.job import Job, concat_jobs
from core.pipeline import Pipeline
from bfx.design import parse_design_file
from bfx.readset import parse_illumina_readset_file, parse_nanopore_readset_file
from bfx.sample_tumor_pairs import *

from bfx import metrics
from bfx import bvatools
from bfx import verify_bam_id
from bfx import picard
from bfx import trimmomatic
from bfx import variantBam
from bfx import samtools
from bfx import rmarkdown
from bfx import jsonator
from bfx import bash_cmd as bash

log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all MUGQIC pipelines (readsets, samples, remote log, etc.)
class MUGQICPipeline(Pipeline):

    def __init__(self, *args, readsets_file=None, design_file=None, **kwargs):
        # Add pipeline specific arguments
        self._readsets = None
        self._contrasts = None
        self._readsets_file = readsets_file
        self._design_file = design_file
        self._samples = None
        super(MUGQICPipeline, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-r", "--readsets", dest="readsets_file",
                                    help="readset file", type=argparse.FileType('r'), required=True)
        cls._argparser.add_argument("-d", "--design", dest="design_file", help="design file",
                                    type=argparse.FileType('r'))

        cls._argparser.description = "Version: " + cls.genpipes_version() + \
                                     "\n\nFor more documentation, visit our website: https://bitbucket.org/mugqic/genpipes/"
        cls._argparser.add_argument("-v", "--version", action="version",
                                    version="genpipes " + cls.genpipes_version(), help="show the version information and exit")
        return cls._argparser

    @property
    def contrasts(self):
        if getattr(self, "_contrasts") is None:
                self._contrasts = parse_design_file(self.design_file, self.samples)
        return self._contrasts

    @property
    def design_file(self):
        if self._design_file is None:
            raise MissingInputError("Design file is required for this pipeline/protocol"
                                    "look at help for details")
        return self._design_file

    @property
    def readsets(self):
        raise NotImplementedError

    @property
    def readsets_file(self):
        return self._readsets_file

    @property
    def samples(self):
        if self._samples is None:
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
        return self._samples

    def mugqic_log(self):
        if 'NO_MUGQIC_REPORT' in os.environ:
            return None
        server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi"
        listName = {}
        for readset in self.readsets:
            if readset.sample.name in listName:
                listName[readset.sample.name] +="."+readset.name
            else:
                listName[readset.sample.name] = readset.sample.name+"."+readset.name

        # The unique identifier is computed from:
        # - Pipeline name
        # - Username (at runtime, not at script creation)
        # - Server name
        # - Readset File

        hostName = socket.gethostname()
        serverIP = socket.gethostbyname(hostName)
        pipelineName = self.__class__.__name__
        readsetFiles = ",".join(listName.values())
        uniqueIdentifier = "{serverIP}-{pipelineName}-{readsetFiles}" \
            .format(serverIP=serverIP, pipelineName=pipelineName, readsetFiles=readsetFiles) \
            .replace("'", "''")

        request = '&'.join([
            "hostname=" + hostName,
            "ip=" + serverIP,
            "pipeline=" + pipelineName,
            "steps=" + ",".join([step.name for step in self.step_to_execute]),
            "samples=" + str(len(self.samples))
        ])
        # that is crazy, to have to rely on the bash interface/arguments that deep in the code.
        self.job_scheduler.write("""
{separator_line}
# Call home with pipeline statistics
{separator_line}
LOG_MD5=$(echo $USER-'{uniqueIdentifier}' | md5sum | awk '{{ print $1 }}')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet '{server}?{request}&md5=$LOG_MD5' -O /dev/null || echo "${{bold}}${{yellow}}Warning:${{normal}}${{yellow}} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${{normal}}"
""".format(separator_line = "#" + "-" * 79, server=server, request=request, uniqueIdentifier=uniqueIdentifier))

    def submit_jobs(self):
        super(MUGQICPipeline, self).submit_jobs()
        if self.jobs and self.job_scheduler.name.lower() in ["pbs", "batch", "slurm"]:
            self.mugqic_log()




# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Nanopore(MUGQICPipeline):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def readsets(self):
        if getattr(self, "_readsets") is None:
            self._readsets = parse_nanopore_readset_file(self.readsets_file)
        return self._readsets


# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Illumina(MUGQICPipeline):

    def __init__(self, *args, **kwargs):
        super(Illumina, self).__init__(*args, **kwargs)

    @property
    def readsets(self):
        if getattr(self, "_readsets") is None:
            self._readsets = parse_illumina_readset_file(self.readsets_file)
        return self._readsets


    @property
    def run_type(self):
        run_types = [readset.run_type for readset in self.readsets]
        if len(set(run_types)) == 1 and re.search("^(PAIRED|SINGLE)_END$", run_types[0]):
            return run_types[0]
        else:
            _raise(SanitycheckError("Error: readset run types " + ","
                                    .join(["\"" + run_type + "\"" for run_type in run_types]) +
            " are invalid (should be all PAIRED_END or all SINGLE_END)!"))


    def samtools_bam_sort(self):
        """
        Sorts bam by readname prior to picard_sam_to_fastq step in order to minimize memory consumption.
        If bam file is small and the memory requirements are reasonable, this step can be skipped.
        """

        jobs = []
        for readset in self.readsets:
            # If readset FASTQ files are available, skip this step
            if not readset.fastq1:
                if readset.bam:
                    sortedBamDirectory = os.path.join(
                        self.output_dir,
                        "temporary_bams",
                        readset.sample.name
                    )
                    sortedBamPrefix = os.path.join(
                        sortedBamDirectory,
                        readset.name + ".sorted"
                    )

                    mkdir_job = bash.mkdir(sortedBamDirectory, remove=True)

                    sort_job = samtools.sort(
                        readset.bam,
                        sortedBamPrefix,
                        sort_by_name = True
                    )
                    sort_job.removable_files = [sortedBamPrefix + ".bam"]

                    jobs.append(
                        concat_jobs([
                            mkdir_job,
                            sort_job
                        ], name="samtools_bam_sort."+readset.name, samples=[readset.sample])
                    )
                else:
                    _raise(SanitycheckError("Error: BAM file not available for readset \"" + readset.name + "\"!"))
        return jobs


    def picard_sam_to_fastq(self):
        """
        Convert SAM/BAM files from the input readset file into FASTQ format
        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
        """
        jobs = []
        analyses_dir = os.path.join("analyses")
        for readset in self.readsets:
            # If readset FASTQ files are available, skip this step
            sym_link_job = []
            if not readset.fastq1:
                if readset.bam:
                    ## check if bam file has been sorted:
                    sortedBam = os.path.join(
                        self.output_dir,
                        "temporary_bams",
                        readset.sample.name,
                        readset.name + ".sorted.bam"
                    )
                    candidate_input_files = [
                        [sortedBam],
                        [readset.bam]
                    ]
                    [bam] = self.select_input_files(candidate_input_files)

                    rawReadsDirectory = os.path.join(
                        self.output_dir,
                        "raw_reads",
                        readset.sample.name,
                    )
                    if readset.run_type == "PAIRED_END":
                        fastq1 = os.path.join(rawReadsDirectory, readset.name + ".pair1.fastq.gz")
                        fastq2 = os.path.join(rawReadsDirectory, readset.name + ".pair2.fastq.gz")
                    elif readset.run_type == "SINGLE_END":
                        fastq1 = os.path.join(rawReadsDirectory, readset.name + ".single.fastq.gz")
                        fastq2 = None
                    else:
                        _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

                    mkdir_job = bash.mkdir(rawReadsDirectory)
                    jobs.append(
                        concat_jobs([
                            mkdir_job,
                            picard.sam_to_fastq(
                                bam,
                                fastq1,
                                fastq2
                            )
                        ], name="picard_sam_to_fastq."+readset.name, samples=[readset.sample])
                    )
                else:
                    _raise(SanitycheckError("Error: BAM file not available for readset \"" + readset.name + "\"!"))
        return jobs

    def trimmomatic(self):
        """
        Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
        it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
        reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
        only Adapter1 is used and left unchanged.

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """
        jobs = []
        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = global_config_parser.param('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
            adapter_job = None
            if not adapter_fasta:
                adapter_fasta = trim_file_prefix + "adapters.fa"
                if readset.run_type == "PAIRED_END":
                    if readset.adapter1 and readset.adapter2:
                        # WARNING: Reverse-complement and swap readset adapters for Trimmomatic Palindrome strategy
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Prefix/1
{sequence1}
>Prefix/2
{sequence2}
END
`""".format(adapter_fasta=adapter_fasta, sequence1=readset.adapter2.translate(str.maketrans("ACGTacgt","TGCAtgca"))[::-1], sequence2=readset.adapter1.translate(str.maketrans("ACGTacgt","TGCAtgca"))[::-1]))
                    else:
                        _raise(SanitycheckError("Error: missing adapter1 and/or adapter2 for PAIRED_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Single
{sequence}
END
`""".format(adapter_fasta=adapter_fasta, sequence=readset.adapter1))
                    else:
                        _raise(SanitycheckError("Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))

            trim_stats = trim_file_prefix + "stats.csv"
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_fastq1 = os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair1.fastq.gz")
                    candidate_fastq2 = os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".pair2.fastq.gz")
                    candidate_input_files.append([candidate_fastq1, candidate_fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([os.path.join(self.output_dir, "raw_reads", readset.sample.name, readset.name + ".single.fastq.gz")])
                [fastq1] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            if adapter_job:
                job = concat_jobs([adapter_job, job])

            jobs.append(concat_jobs([
                # Trimmomatic does not create output directory by default
                Job(command="mkdir -p " + trim_directory, samples=[readset.sample]),
                job
            ], name="trimmomatic." + readset.name, samples=[readset.sample]))
        return jobs

    def merge_trimmomatic_stats(self):
        """
        The trim statistics per readset are merged at this step.
        """

        read_type = "Paired" if self.run_type == 'PAIRED_END' else "Single"
        readset_merge_trim_stats = os.path.join("metrics", "trimReadsetTable.tsv")
        job = concat_jobs([Job(command="mkdir -p metrics"), Job(command="echo 'Sample\tReadset\tRaw {read_type} Reads #\tSurviving {read_type} Reads #\tSurviving {read_type} Reads %' > ".format(read_type=read_type) + readset_merge_trim_stats)])
        for readset in self.readsets:
            trim_log = os.path.join("trim", readset.sample.name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                # Retrieve readset raw and surviving reads from trimmomatic log using ugly Perl regexp
                perl_command = "perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/{readset.sample.name}\t{readset.name}\t\\1\t\\2/'".format(readset=readset)
            elif readset.run_type == "SINGLE_END":
                perl_command = "perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/{readset.sample.name}\t{readset.name}\t\\1\t\\2/'".format(readset=readset)

            job = concat_jobs([
                job,
                Job(
                    [trim_log],
                    [readset_merge_trim_stats],
                    module_entries=[['merge_trimmomatic_stats', 'module_perl']],
                    # Create readset trimming stats TSV file with paired or single read count using ugly awk
                    command="""\
grep ^Input {trim_log} | \\
{perl_command} | \\
awk '{{OFS="\t"; print $0, $4 / $3 * 100}}' \\
  >> {readset_merge_trim_stats}""".format(
                        trim_log=trim_log,
                        perl_command=perl_command,
                        readset_merge_trim_stats=readset_merge_trim_stats
                    ),
                    samples=[readset.sample]
                )
            ])

        sample_merge_trim_stats = os.path.join("metrics", "trimSampleTable.tsv")
        report_file = os.path.join("report", "Illumina.merge_trimmomatic_stats.md")
        return [concat_jobs([
            job,
            Job(
                [readset_merge_trim_stats],
                [sample_merge_trim_stats],
                # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                command="""\
cut -f1,3- {readset_merge_trim_stats} | awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_trim_stats}""".format(
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats
                )
            ),
            Job(
                [sample_merge_trim_stats],
                [report_file, os.path.join("report", "trimReadsetTable.tsv"), os.path.join("report", "trimSampleTable.tsv")],
                [['merge_trimmomatic_stats', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} report/ && \\
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_trim_stats}` && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trailing_min_quality={trailing_min_quality} \\
  --variable min_length={min_length} \\
  --variable read_type={read_type} \\
  --variable trim_readset_table="$trim_readset_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                    trailing_min_quality=global_config_parser.param('trimmomatic', 'trailing_min_quality', param_type='int'),
                    min_length=global_config_parser.param('trimmomatic', 'min_length', param_type='posint'),
                    read_type=read_type,
                    report_template_dir=self.report_template_dir,
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file]
            )], name="merge_trimmomatic_stats")]

    def verify_bam_id(self):
        """
        verifyBamID is a software that verifies whether the reads in particular file match previously known
        genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
        as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external
        genotypes are available. When external genotypes are not available, verifyBamID still robustly
        detects sample swaps.
        """

        # Known variants file
        population_AF = global_config_parser.param('verify_bam_id', 'population_AF', required=False)
        candidate_input_files = [[global_config_parser.param('verify_bam_id', 'verifyBamID_variants_file', required=False)]]
        candidate_input_files.append([global_config_parser.param('verify_bam_id', 'verifyBamID_variants_file', required=False) + ".gz"])
        [known_variants_annotated] = self.select_input_files(candidate_input_files)
        verify_bam_id_directory = "verify_bam_id"
        variants_directory = "variants"

        jobs = []

        verify_bam_results = []

        jobs.append(
            Job(
                [known_variants_annotated],
                [variants_directory, verify_bam_id_directory],
                command="mkdir -p " + variants_directory + " " + verify_bam_id_directory,
                name = "verify_bam_id_create_directories"
        ))

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)

            candidate_input_files = [[os.path.join(alignment_directory, sample.name + ".sorted.dup.recal.bam")]]
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")])
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]) # this one is for RnaSeq pipeline
            [input_bam] = self.select_input_files(candidate_input_files)

            output_prefix = os.path.join(verify_bam_id_directory, sample.name)

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            # Run verifyBamID
            job = verify_bam_id.verify(
                input_bam,
                known_variants_annotated,
                output_prefix
            )
            job.name = "verify_bam_id." + sample.name
            job.samples = [sample]

            jobs.append(job)

            verify_bam_results.extend([output_prefix + ".selfSM" ])

        # Coverage bed is null if whole genome experiment
        target_bed=coverage_bed if coverage_bed else ""

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input            = verify_bam_results ,
                job_name             = "verify_bam_id_report",
                input_rmarkdown_file = os.path.join(self.report_template_dir, "Illumina.verify_bam_id.Rmd"),
                samples              = self.samples,
                render_output_dir    = 'report',
                module_section       = 'report',
                prerun_r             = 'source_dir="' + verify_bam_id_directory + '"; report_dir="report" ; params=list(verifyBamID_variants_file="' + known_variants_annotated  + '", dbnsfp_af_field="' + population_AF + '", coverage_bed="' + target_bed + '");'
            )
        )

        return jobs

    def cram_output(self):
        """
        Generate long term storage version of the final alignment files in CRAM format
        Using this function will include the orginal final bam file into the  removable file list
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)

            candidate_input_files = [[os.path.join(alignment_directory, sample.name + ".sorted.dup.recal.bam")]] #DNAseq; TumorPair
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.dedup.bam")]) #DNAseq; ChIPseq; MethylSeq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.mdup.bam")]) # RnaSeq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.dup.bam")]) # TumorPair
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".merded.mdup.bam")]) # ChIPseq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + "matefixed.sorted.bam")]) # DNAseq_high_coverage
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".realigned.sorted.bam")]) # TumorPair; DNAseq_high_coverage
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".merged.bam")]) # ChIPSeq; HiCseq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.filtered.bam")]) # ChIPseq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted_noRG.bam")]) # MethylSeq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.UMI.bam")]) # MethylSeq
            candidate_input_files.append([os.path.join(alignment_directory, sample.name + ".sorted.bam")]) # raw bam
            input_bam = self.select_input_files(candidate_input_files)[0]

            output_cram = re.sub("\.bam$", ".cram", input_bam)

            if "recal" in input_bam:
                job = variantBam.run(
                    input_bam,
                    output_cram,
                )
            else:
                # Run samtools
                job = samtools.view(
                    input_bam,
                    output_cram,
                    options=global_config_parser.param('samtools_cram_output', 'options'),
                    removable=False
                )

            job.name = "samtools_cram_output." + sample.name
            job.samples = [sample]
            job.removable_files = input_bam

            jobs.append(job)

        return jobs


class Error(Exception):
    """
    base error for common pipelines class
    """
    pass


class MissingInputError(Error):
    pass
