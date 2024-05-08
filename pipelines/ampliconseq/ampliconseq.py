#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import math
import os
import re
import sys
from os.path import basename

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs
import utils.utils

from pipelines import common
from bfx import bash_cmd as bash
from bfx import tools
from bfx import dada2
from bfx import flash
from bfx import vsearch
from bfx import trimmomatic

log = logging.getLogger(__name__)

class AmpliconSeq(common.Illumina):
    """
    Amplicon-Seq Pipeline
    ================

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=argparse.FileType('r'))
        super(AmpliconSeq, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
            'merge_directory': os.path.relpath(os.path.join(self.output_dir, 'merge'), self.output_dir),
            'dada2_analysis_directory': os.path.relpath(os.path.join(self.output_dir, 'dada2_Analysis'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
        }
        return dirs

    def trimmomatic16S(self):
        """
        MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
        it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

        This step takes as input files:
        1. MiSeq paired-End FASTQ files from the readset file
        """

        jobs = []
        #We'll trim the first 5 nucleotides anyway (to account for quality bias)
        headcropValue=5
        for readset in self.readsets:
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = config.param('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
            if not adapter_fasta:
                adapter_fasta = trim_file_prefix + "adapters.fa"
                if readset.primer1 and readset.primer2:
                    #concatenate all adpaters and primers
                    primAdapList = readset.primer1 + ";" + readset.primer2
                    #convert into a list
                    primAdapList = primAdapList.split(";")

                    ambiChar=('N','R','Y','K','M','S','W','B','D','H','V')
                    #check the highest position of any ambiguous nucleotide
                    for item in primAdapList:
                        for char in ambiChar:
                            checkAmbiPos = item.rfind('%s'%char)
                            if (checkAmbiPos + 1) > headcropValue:
                                headcropValue = checkAmbiPos + 1

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append(
                        [
                            re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                            re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)
                        ]
                    )
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic16S(
                    fastq1,
                    fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    trim_log,
                    headcropValue
                )
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic16S(
                    fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    trim_log,
                    headcropValue
                )
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        # Trimmomatic does not create output directory by default
                        bash.mkdir(trim_directory),
                        job
                    ],
                    name="trimmomatic16S." + readset.name,
                    samples=[readset.sample]
                )
            )
        return jobs

    def merge_trimmomatic_stats16S(self):
        """
        The trim statistics per readset are merged at this step.
        """

        read_type = "Paired" if self.run_type == 'PAIRED_END' else "Single"
        readset_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimReadsetTable.tsv")
        job = concat_jobs(
            [
                bash.mkdir(self.output_dirs["metrics_directory"]),
                Job(command="echo 'Sample\tReadset\tRaw {read_type} Reads #\tSurviving {read_type} Reads #\tSurviving {read_type} Reads %' > ".format(read_type=read_type) + readset_merge_trim_stats)
            ]
        )
        for readset in self.readsets:
            trim_log = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                # Retrieve readset raw and surviving reads from trimmomatic log using ugly Perl regexp
                perl_command = "perl -pe 's/^Input Read Pairs: (\d+).*Both Surviving: (\d+).*Forward Only Surviving: (\d+).*$/{readset.sample.name}\t{readset.name}\t\\1\t\\2/'".format(readset=readset)
            elif readset.run_type == "SINGLE_END":
                perl_command = "perl -pe 's/^Input Reads: (\d+).*Surviving: (\d+).*$/{readset.sample.name}\t{readset.name}\t\\1\t\\2/'".format(readset=readset)

            job = concat_jobs(
                [
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
                        ),
                        samples=[readset.sample]
                    )
                ]
            )

        sample_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimSampleTable.tsv")
        report_file = os.path.join(self.output_dirs["report_directory"], "Illumina.merge_trimmomatic_stats.md")
        return [
            concat_jobs(
                [
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
mkdir -p {report_dir} && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} {report_dir}/ && \\
trim_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_trim_stats}` && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable trailing_min_quality="NA" \\
  --variable min_length="NA" \\
  --variable read_type={read_type} \\
  --variable trim_readset_table="$trim_readset_table_md" \\
  --to markdown \\
  > {report_file}""".format(
                            read_type=read_type,
                            report_template_dir=self.report_template_dir,
                            readset_merge_trim_stats=readset_merge_trim_stats,
                            sample_merge_trim_stats=sample_merge_trim_stats,
                            basename_report_file=os.path.basename(report_file),
                            report_dir=self.output_dirs["report_directory"],
                            report_file=report_file
                        ),
                        report_files=[report_file]
                    )
                ],
                name="merge_trimmomatic_stats16S"
            )
        ]

    def flash(self, flash_stats_file=None):
        """
        Merge paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/).
        """
        jobs = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")
            merge_directory = os.path.join(self.output_dirs["merge_directory"], readset.sample.name)
            flash_fastq = os.path.join(merge_directory, readset.name + ".flash_pass2.extendedFrags.fastq") if flash_stats_file else os.path.join(merge_directory, readset.name + ".flash.extendedFrags.fastq")
            flash_log = os.path.join(merge_directory, readset.name + ".flash_pass2.log") if flash_stats_file else os.path.join(merge_directory, readset.name + ".flash.log")
            flash_hist = os.path.join(merge_directory, readset.name + ".flash_pass2.hist") if flash_stats_file else os.path.join(merge_directory, readset.name + ".flash.hist")
            job_name_prefix = "flash_pass2." if flash_stats_file else "flash_pass1."

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!"))

            job = flash.flash(
                fastq1,
                fastq2,
                flash_fastq,
                readset.name,
                flash_log,
                flash_hist,
                flash_stats_file
            )
            job.samples = [readset.sample]

            jobs.append(concat_jobs([
                # FLASh does not create output directory by default
                Job(command="mkdir -p " + merge_directory),
                job
            ], name=job_name_prefix + readset.sample.name))

        return jobs

    def flash_pass1(self):
        jobs = self.flash()
        return jobs

    def flash_pass2(self):
        flash_stats_file = os.path.join(self.output_dirs["metrics_directory"], "FlashLengths.tsv")
        jobs = self.flash(flash_stats_file)
        return jobs

    def merge_flash_stats(self):
        """
        The paired end merge statistics per readset are merged at this step.
        """

        readset_merge_flash_stats = os.path.join(self.output_dirs["metrics_directory"], "mergeReadsetTable.tsv")
        job = concat_jobs(
            [
                bash.mkdir(self.output_dirs["metrics_directory"]),
                Job(command="echo 'Sample\\tReadset\\tTrim Paired Reads #\\tMerged Paired Reads #\\tMerged Paired Reads %' > " + readset_merge_flash_stats)
            ]
        )

        for readset in self.readsets:
            flash_log = os.path.join(self.output_dirs["merge_directory"], readset.sample.name, readset.name + ".flash_pass2.log")

            job = concat_jobs(
                [
                    job,
                    Job(
                        command="""\
printf '{sample}\\t{readset}\\t' \\
  >> {stats}""".format(
                            sample=readset.sample.name,
                            readset=readset.name,
                            stats=readset_merge_flash_stats,
                        ),
                        samples=[readset.sample]
                    )
                ]
            )

            # Retrieve merge statistics using re search in python.
            python_command = """\
python -c 'import re; \\
  import sys; \\
  log_file = open("{flash_log}","r"); \\
  merge_stat=[]; \\
  merge_stat.append([i.split()[3] for i in log_file if re.search("Total pairs",i)][0]); \\
  log_file.seek(0); \\
  merge_stat.append([i.split()[3] for i in log_file if re.search("Combined pairs",i)][0]); \\
  log_file.seek(0); \\
  merge_stat.append([i.split()[3] for i in log_file if re.search("Percent combined",i)][0][:-1]); \\
  log_file.close(); \\
  print "\t".join(merge_stat)'""".format(
                flash_log=flash_log
            )

            job = concat_jobs(
                [
                    job,
                    Job(
                        [flash_log],
                        [readset_merge_flash_stats],
                        [['merge_flash_stats', 'module_python']],
                        # Create readset merging stats TSV file with paired read count using python.
                        command="""\
{python_command} \\
  >> {readset_merge_flash_stats}""".format(
                            python_command=python_command,
                            readset_merge_flash_stats=readset_merge_flash_stats
                        )
                    )
                ]
            )

        sample_merge_flash_stats = os.path.join(self.output_dirs["metrics_directory"], "mergeSampleTable.tsv")
        report_file = os.path.join(self.output_dirs["report_directory"], "Illumina.flash_stats.md")
        return [
            concat_jobs(
                [
                    job,
                    Job(
                        [readset_merge_flash_stats],
                        [sample_merge_flash_stats],
                        # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                        command="""\
cut -f1,3- {readset_merge_flash_stats} | \\
awk -F"\\t" '{{OFS="\\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Trim Reads #", "Merged Reads #", "Merged %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_flash_stats}""".format(
                            readset_merge_flash_stats=readset_merge_flash_stats,
                            sample_merge_flash_stats=sample_merge_flash_stats
                        )
                    ),
                    Job(
                        [sample_merge_flash_stats],
                        [report_file],
                        [['flash', 'module_pandoc']],
                        command="""\
mkdir -p {report_dir} && \\
cp {readset_merge_flash_stats} {sample_merge_flash_stats} {report_dir}/ && \\
merge_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_flash_stats}` && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable min_overlap="{min_overlap}" \\
  --variable max_overlap="{max_overlap}" \\
  --variable read_type="{read_type}" \\
  --variable merge_readset_table="$merge_readset_table_md" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                            min_overlap=config.param('flash', 'min_overlap', param_type='int'),
                            max_overlap=config.param('flash', 'max_overlap', param_type='int'),
                            read_type="Paired",
                            report_template_dir=self.report_template_dir,
                            readset_merge_flash_stats=readset_merge_flash_stats,
                            sample_merge_flash_stats=sample_merge_flash_stats,
                            basename_report_file=os.path.basename(report_file),
                            report_dir=self.output_dirs["report_directory"],
                            report_file=report_file
                        ),
                        report_files=[report_file]
                    )
                ],
                name="merge_flash_stats",
                samples=self.samples
            )
        ]

    def ampliconLengthParser(self):
        """
        look at FLASH output to set amplicon lengths input for dada2. As minimum elligible length, a given length needs to have at least 1% of the total number of amplicons
        """
        jobs = []
        readset_merge_flash_stats = os.path.join(self.output_dirs["metrics_directory"], "FlashLengths.tsv")

        job = concat_jobs(
            [
                bash.mkdir(self.output_dirs["metrics_directory"]),
                Job(command="echo 'Sample\tReadset\tMinimum Amplicon Length\tMaximum Amplicon Length\tMinimum Flash Overlap\tMaximum Flash Overlap' > " + readset_merge_flash_stats)
            ]
        )

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type + "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!"))

            flash_hist = os.path.join(self.output_dirs["merge_directory"], readset.sample.name, readset.name + ".flash.hist")
            job = concat_jobs(
                [
                    job,
                    Job(
                        [fastq1, flash_hist],
                        [readset_merge_flash_stats],
                        command="""\
frag_length=$(zcat {fastq} | head -n2 | tail -n1 | awk '{{print length($0)}}'; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi)
minCount=$(cut -f2 {hist} | sort -n | awk ' {{ sum+=$1;i++ }} END {{ print sum/100; }}' | cut -d"." -f1)
minLen=$(awk -F'\\t' -v var=$minCount '$2>var' {hist} | cut -f1 | sort -g | head -n1)
maxLen=$(awk -F'\\t' -v var=$minCount '$2>var' {hist} | cut -f1 | sort -gr | head -n1)
minFlashOverlap=$(( 2 * frag_length - maxLen ))
maxFlashOverlap=$(( 2 * frag_length - minLen ))
printf "{sample}\\t{readset}\\t${{minLen}}\\t${{maxLen}}\\t${{minFlashOverlap}}\\t${{maxFlashOverlap}}\\n" \\
  >> {stats}""".format(
                            fastq=fastq1,
                            hist=flash_hist,
                            sample=readset.sample.name,
                            readset=readset.name,
                            stats=readset_merge_flash_stats,
                            ),
                        samples=[readset.sample]
                    )
                ]
            )

        job.name = "ampliconLengthParser.run"
        jobs.append(job)

        return jobs

    def asva(self):
        """
        check for design file (required for PCA plots)
        """

        try:
            designFile = os.path.relpath(self.args.design.name, self.args.output_dir)
        except:
            self.argparser.error("argument -d/--design is required!")

        jobs = []

        #Create folders in the output folder
        dada2_directory = self.output_dirs["dada2_analysis_directory"]
        lnkRawReadsFolder = os.path.join(dada2_directory, "trim")
        ampliconLengthFile = os.path.join(self.output_dirs["metrics_directory"], "FlashLengths.tsv")

        #We'll link the readset fastq files into the raw_reads folder just created
        raw_reads_jobs = []
        dada2_inputs = []
        for readset in self.readsets:
            readSetPrefix = os.path.join(lnkRawReadsFolder, readset.name)

            trimmedReadsR1 = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.pair1.fastq.gz")
            trimmedReadsR2 = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.pair2.fastq.gz")

            if readset.run_type == "PAIRED_END":
                left_or_single_reads = readSetPrefix + ".pair1.fastq.gz"
                dada2_inputs.append(left_or_single_reads)

                right_reads = readSetPrefix + ".pair2.fastq.gz"
                dada2_inputs.append(right_reads)

                raw_reads_jobs.append(
                    concat_jobs(
                        [
                            Job(
                                [trimmedReadsR1],
                                [left_or_single_reads],
                                command="ln -nsf " + os.path.abspath(os.path.join(self.args.output_dir, trimmedReadsR1)) + " " + left_or_single_reads,
                                samples=[readset.sample]
                            ),
                            Job(
                                [trimmedReadsR2],
                                [right_reads],
                                command="ln -nsf " + os.path.abspath(os.path.join(self.args.output_dir, trimmedReadsR2)) + " " + right_reads,
                                samples=[readset.sample]
                            )
                        ]
                    )
                )

            #single reads will mainly be for PacBio CCS although I didn't test it yet
            elif readset.run_type == "SINGLE_END":
                left_or_single_reads = readSetPrefix + ".single.fastq.gz"
                raw_reads_jobs.append(
                    Job(
                        [trimmedReadsR1],
                        [left_or_single_reads],
                        command="ln -nsf " + trimmedReadsR1 + " " + left_or_single_reads,
                        samples=[readset.sample]
                    )
                )
                dada2_inputs.append(left_or_single_reads)

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(lnkRawReadsFolder)
                ] + raw_reads_jobs + [
                    dada2.dada2(
                        dada2_inputs,
                        ampliconLengthFile,
                        lnkRawReadsFolder,
                        designFile,
                        dada2_directory
                    )
                ],
                name="dada2.run"
            )
        )
        return jobs

    @property
    def steps(self):
        return [
            self.trimmomatic16S,
            self.merge_trimmomatic_stats16S,
            self.flash_pass1,
            self.ampliconLengthParser,
            self.flash_pass2,
            self.merge_flash_stats,
            self.asva
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        AmpliconSeq()
