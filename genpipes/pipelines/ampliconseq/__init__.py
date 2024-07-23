################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import os
import re

# GenPipes Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs
from .. import common

from ...bfx import (
    bash_cmd as bash,
    dada2,
    flash,
    multiqc,
    trimmomatic
    )

log = logging.getLogger(__name__)

class AmpliconSeq(common.Illumina):
    """
    Amplicon-Seq Pipeline
    ================

    """

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        # Add pipeline specific arguments
        super(AmpliconSeq, self).__init__(*args, **kwargs)

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
    
    @property
    def multiqc_inputs(self):
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    def trimmomatic16S(self):
        """
        MiSeq raw reads adapter & primers trimming and basic QC is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
        it is used first. Else, Adapter1, Adapter2, Primer1 and Primer2 columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. Sequences are reversed-complemented and swapped.

        This step takes as input files:
        1. MiSeq paired-End FASTQ files from the readset file
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        #We'll trim the first 5 nucleotides anyway (to account for quality bias)
        headcropValue=5
        for readset in self.readsets:
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = global_conf.global_get('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
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
                        bash.mkdir(trim_directory),
                        bash.mkdir(link_directory),
                        job,
                        bash.ln(
                            os.path.relpath(trim_log, link_directory),
                            os.path.join(link_directory, readset.name + ".trim.log"),
                            input = trim_log
                        )
                    ], 
                    name="trimmomatic16S." + readset.name
                    )
                )
            self.multiqc_inputs.append(os.path.join(link_directory, readset.name + ".trim.log"))

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
                Job(command="echo -e 'Sample\\tReadset\\tRaw {read_type} Reads #\\tSurviving {read_type} Reads #\\tSurviving {read_type} Reads %' > ".format(read_type=read_type) + readset_merge_trim_stats)
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
            ])

        sample_merge_trim_stats = os.path.join("metrics", "trimSampleTable.tsv")
        return [concat_jobs(
            [
                job,
                Job(
                    [readset_merge_trim_stats],
                    [sample_merge_trim_stats],
                    # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                    command="""\
cut -f1,3- {readset_merge_trim_stats} | awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_trim_stats} && \\
mkdir -p report && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} report/""".format(
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats
                )
            )
            ], name="merge_trimmomatic_stats16S")]

    def flash(self, flash_stats_file=None):
        """
        Merge paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/).
        """
        jobs = []
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

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

            if flash_stats_file:
                link_job = concat_jobs(
                    [
                        bash.mkdir(
                            link_directory
                        ),
                        bash.ln(
                            os.path.relpath(flash_log, link_directory),
                            os.path.join(link_directory, readset.name + ".flash.log"),
                            input = flash_log
                        ),
                        bash.ln(
                            os.path.relpath(flash_hist, link_directory),
                            os.path.join(link_directory, readset.name + ".flash.hist"),
                            input = flash_hist
                        )
                    ]
                )
                self.multiqc_inputs.append(os.path.join(link_directory, readset.name + ".flash.log"))
                self.multiqc_inputs.append(os.path.join(link_directory, readset.name + ".flash.hist"))
            else:
                link_job = None

            job.samples = [readset.sample]

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(merge_directory),
                        job,
                        link_job
                    ], 
                    name=job_name_prefix + readset.sample.name
                )
            )
            
        return jobs

    def flash_pass1(self):
        """
        Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). Overlapping regions between paired-end reads are found and 
        then merged into a continuous strand.
        """
        jobs = self.flash()
        return jobs

    def flash_pass2(self):
        """
        Merges paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/). The second pass uses statistics obtained from the first pass
        to adjust merging.
        """
        flash_stats_file = os.path.join(self.output_dirs["metrics_directory"], "FlashLengths.tsv")
        jobs = self.flash(flash_stats_file)
        return jobs

    def merge_flash_stats(self):
        """
        Merges statistics from both flash passes.
        """

        readset_merge_flash_stats = os.path.join(self.output_dirs["metrics_directory"], "mergeReadsetTable.tsv")
        job = concat_jobs(
            [
                bash.mkdir(self.output_dirs["metrics_directory"]),
                Job(command="echo -e 'Sample\\tReadset\\tTrim Paired Reads #\\tMerged Paired Reads #\\tMerged Paired Reads %' > " + readset_merge_flash_stats)
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
            ])

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
                        [
                            ['merge_flash_stats', 'module_python']
                        ],
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

        sample_merge_flash_stats = os.path.join("metrics", "mergeSampleTable.tsv")

        return [concat_jobs(
            [
                job,
                Job(
                    [readset_merge_flash_stats],
                    [sample_merge_flash_stats],
                    # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                    command="""\
cut -f1,3- {readset_merge_flash_stats} | \\
awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Trim Reads #", "Merged Reads #", "Merged %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_flash_stats} && \\
mkdir -p report && \\
cp {readset_merge_flash_stats} {sample_merge_flash_stats} report/""".format(
                    readset_merge_flash_stats=readset_merge_flash_stats,
                    sample_merge_flash_stats=sample_merge_flash_stats
                )
            )
            ], name="merge_flash_stats", samples=self.samples)]

    def ampliconLengthParser(self):
        """
        Looks at FLASH output statistics to set input amplicon lengths for dada2. Minimum lengths are set by ensuring that they represent 
        at least 1% of the total number of amplicons.
        """
        jobs = []
        readset_merge_flash_stats = os.path.join(self.output_dirs["metrics_directory"], "FlashLengths.tsv")

        job = concat_jobs(
            [
                bash.mkdir(self.output_dirs["metrics_directory"]),
                Job(command="echo -e 'Sample\\tReadset\\tMinimum Amplicon Length\\tMaximum Amplicon Length\\tMinimum Flash Overlap\\tMaximum Flash Overlap' > " + readset_merge_flash_stats)
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
frag_length=$(zcat {fastq} | head -n2 | tail -n1 | awk '{{print length($0)}}')
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
                    )
                ],
                samples=[readset.sample],
                name="ampliconLengthParser.run"
            )

        jobs.append(job)

        return jobs

    def asva(self):
        """
        Checks for design files required for PCA plots, sets up directories, links readset fastq files, and initiates 
        [DADA2](https://benjjneb.github.io/dada2/). 

        DADA2 is used to infer sequence variants of microbial communities.
        """

        designFile = os.path.relpath(self.design_file.name, self.output_dir)

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
                            bash.ln(
                                os.path.relpath(trimmedReadsR1, lnkRawReadsFolder),
                                left_or_single_reads,
                                input=trimmedReadsR1
                            ),
                            bash.ln(
                                os.path.relpath(trimmedReadsR2, lnkRawReadsFolder),
                                right_reads,
                                input=trimmedReadsR2
                            )
                        ],
                        samples=[readset.sample]
                    )
                )

            #single reads will mainly be for PacBio CCS although I didn't test it yet
            elif readset.run_type == "SINGLE_END":
                left_or_single_reads = readSetPrefix + ".single.fastq.gz"
                raw_reads_jobs.append(
                    concat_jobs(
                        [
                            bash.ln(
                                os.path.relpath(trimmedReadsR1, lnkRawReadsFolder),
                                left_or_single_reads,
                                input=trimmedReadsR1
                            )
                        ],
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
    
    def multiqc(self):
        """
        A quality control report for all samples is generated.
        For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)
        """
        jobs = []
        
        input_links = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
            
        output = os.path.join(self.output_dirs['report_directory'], f"AmpliconSeq.multiqc")

        job = multiqc.run(
            input_links,
            output,
            ini_section='multiqc'
            )
        
        job.name = "multiqc"
        jobs.append(job)

        return jobs        


    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return {"default": [
                self.trimmomatic16S,
                self.merge_trimmomatic_stats16S,
                self.flash_pass1,
                self.ampliconLengthParser,
                self.flash_pass2,
                self.merge_flash_stats,
                self.asva,
                self.multiqc
                ]
            }


def main(parsed_args):
    """
    The function that will call this pipeline!
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file

    pipeline = AmpliconSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, container=container)

    pipeline.submit_jobs()
