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
import re
import socket
import string
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

# MUGQIC Modules
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

from bfx import metrics
from bfx import picard
from bfx import trimmomatic

log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all MUGQIC pipelines (readsets, samples, remote log, etc.)
class MUGQICPipeline(Pipeline):

    def __init__(self):
        self.version = open(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "VERSION"), 'r').read().split('\n')[0]

        # Add pipeline specific arguments
        self.argparser.description = "Version: " + self.version + "\n\nFor more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/"
        self.argparser.add_argument("-v", "--version", action="version", version="mugqic_pipelines " + self.version, help="show the version information and exit")

        super(MUGQICPipeline, self).__init__()

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        if not hasattr(self, "_samples"):
            self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self.readsets]))
        return self._samples

    def mugqic_log(self):
        server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi"
        request = \
            "hostname=" + socket.gethostname() + "&" + \
            "ip=" + socket.gethostbyname(socket.gethostname()) + "&" + \
            "pipeline=" + self.__class__.__name__ + "&" + \
            "steps=" + ",".join([step.name for step in self.step_range]) + "&" + \
            "samples=" + str(len(self.samples))

        print("""
{separator_line}
# Call home with pipeline statistics
{separator_line}
wget "{server}?{request}" --quiet --output-document=/dev/null
""".format(separator_line = "#" + "-" * 79, server=server, request=request))


    def submit_jobs(self):
        super(MUGQICPipeline, self).scheduler.submit(self)
        if self.jobs and self.args.job_scheduler in ["pbs", "batch"]:
            self.mugqic_log()


# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Illumina(MUGQICPipeline):

    def __init__(self):
        self.argparser.add_argument("-r", "--readsets", help="readset file", type=file)
        super(Illumina, self).__init__()

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            if self.args.readsets:
                self._readsets = parse_illumina_readset_file(self.args.readsets.name)
            else:
                self.argparser.error("argument -r/--readsets is required!")
        return self._readsets

    @property
    def run_type(self):
        run_types = [readset.run_type for readset in self.readsets]
        if len(set(run_types)) == 1 and re.search("^(PAIRED|SINGLE)_END$", run_types[0]):
            return run_types[0]
        else:
            raise Exception("Error: readset run types " + ",".join(["\"" + run_type + "\"" for run_type in run_types]) +
            " are invalid (should be all PAIRED_END or all SINGLE_END)!")

    @property
    def contrasts(self):
        if not hasattr(self, "_contrasts"):
            if self.args.design:
                self._contrasts = parse_design_file(self.args.design.name, self.samples)
            else:
                self.argparser.error("argument -d/--design is required!")
        return self._contrasts

    def picard_sam_to_fastq(self):
        """
        Convert SAM/BAM files from the input readset file into FASTQ format
        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
        """
        jobs = []
        for readset in self.readsets:
            # If readset FASTQ files are available, skip this step
            if not readset.fastq1:
                if readset.bam:
                    if readset.run_type == "PAIRED_END":
                        fastq1 = re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)
                        fastq2 = re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                    elif readset.run_type == "SINGLE_END":
                        fastq1 = re.sub("\.bam$", ".single.fastq.gz", readset.bam)
                        fastq2 = None
                    else:
                        raise Exception("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                    job = picard.sam_to_fastq(readset.bam, fastq1, fastq2)
                    job.name = "picard_sam_to_fastq." + readset.name
                    jobs.append(job)
                else:
                    raise Exception("Error: BAM file not available for readset \"" + readset.name + "\"!")
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
            adapter_fasta = config.param('trimmomatic', 'adapter_fasta', required=False, type='filepath')
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
`""".format(adapter_fasta=adapter_fasta, sequence1=readset.adapter2.translate(string.maketrans("ACGTacgt","TGCAtgca"))[::-1], sequence2=readset.adapter1.translate(string.maketrans("ACGTacgt","TGCAtgca"))[::-1]))
                    else:
                        raise Exception("Error: missing adapter1 and/or adapter2 for PAIRED_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!")
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Single
{sequence}
END
`""".format(adapter_fasta=adapter_fasta, sequence=readset.adapter1))
                    else:
                        raise Exception("Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!")


            trim_stats = trim_file_prefix + "stats.csv"
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
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
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
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
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            if adapter_job:
                job = concat_jobs([adapter_job, job])
            jobs.append(concat_jobs([
                # Trimmomatic does not create output directory by default
                Job(command="mkdir -p " + trim_directory),
                job
            ], name="trimmomatic." + readset.name))
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
                    # Create readset trimming stats TSV file with paired or single read count using ugly awk
                    command="""\
grep ^Input {trim_log} | \\
{perl_command} | \\
awk '{{OFS="\t"; print $0, $4 / $3 * 100}}' \\
  >> {readset_merge_trim_stats}""".format(
                        trim_log=trim_log,
                        perl_command=perl_command,
                        readset_merge_trim_stats=readset_merge_trim_stats
                    )
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
                [report_file],
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
                    trailing_min_quality=config.param('trimmomatic', 'trailing_min_quality', type='int'),
                    min_length=config.param('trimmomatic', 'min_length', type='posint'),
                    read_type=read_type,
                    report_template_dir=self.report_template_dir,
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file]
            )], name="merge_trimmomatic_stats")]
