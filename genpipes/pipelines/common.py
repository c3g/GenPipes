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
import argparse
import logging
import os
import re
import socket
import sys
import collections

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

# GenPipes Modules
from ..core.config import global_conf, _raise, SanitycheckError
from ..core.job import Job, concat_jobs, pipe_jobs
from ..core.pipeline import Pipeline
from ..core.design import parse_design_file
from ..core.readset import parse_illumina_readset_file, parse_nanopore_readset_file
from ..core.sample_tumor_pairs import *

from ..bfx import (
    adapters,
    bash_cmd as bash,
    bvatools,
    bwa2,
    fastp,
    gatk4,
    job2json_project_tracking,
    picard,
    rmarkdown,
    samtools,
    sambamba,
    skewer,
    trimmomatic,
    variantBam,
    verify_bam_id
)

log = logging.getLogger(__name__)

# Abstract pipeline gathering common features of all GinPipes pipelines (readsets, samples, remote log, etc.)
class GenPipesPipeline(Pipeline):

    def __init__(self, *args, readsets_file=None, design_file=None, **kwargs):
        # Add pipeline specific arguments
        self._readsets = None
        self._contrasts = None
        self._readsets_file = readsets_file
        self._design_file = design_file
        self._samples = None
        super(GenPipesPipeline, self).__init__(*args, **kwargs)

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
            raise MissingInputError("Design file is required for this pipeline/protocol. "
                                    "Look at help for details")
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

    def genpipes_log(self):
        if 'NO_GENPIPES_REPORT' in os.environ:
            return None
        server = "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi"
        list_name = {}
        for readset in self.readsets:
            if readset.sample.name in list_name:
                list_name[readset.sample.name]+="."+readset.name
            else:
                list_name[readset.sample.name]=readset.sample.name+"."+readset.name

        # The unique identifier is computed from:
        # - Pipeline name
        # - Username (at runtime, not at script creation)
        # - Server name
        # - Readset File

        host_name = socket.gethostname()
        server_ip = socket.gethostbyname(host_name)
        pipeline_name = self.__class__.__name__
        readset_files = ",".join(list_name.values())
        unique_identifier = f"{server_ip}-{pipeline_name}-{readset_files}".replace("'", "''")

        request = '&'.join([
            "hostname=" + host_name,
            "ip=" + server_ip,
            "pipeline=" + pipeline_name,
            "steps=" + ",".join([step.name for step in self.step_to_execute]),
            "samples=" + str(len(self.samples))
        ])
        # that is crazy, to have to rely on the bash interface/arguments that deep in the code.
        self.job_scheduler.write("""
{separator_line}
# Call home with pipeline statistics
{separator_line}
LOG_MD5=$(echo $USER-'{unique_identifier}' | md5sum | awk '{{ print $1 }}')
if test -t 1; then ncolors=$(tput colors); if test -n "$ncolors" && test $ncolors -ge 8; then bold="$(tput bold)"; normal="$(tput sgr0)"; yellow="$(tput setaf 3)"; fi; fi
wget --quiet '{server}?{request}&md5=$LOG_MD5' -O /dev/null || echo "${{bold}}${{yellow}}Warning:${{normal}}${{yellow}} Genpipes ran successfully but was not send telemetry to mugqic.hpc.mcgill.ca. This error will not affect genpipes jobs you have submitted.${{normal}}"
""".format(separator_line = "#" + "-" * 79, server=server, request=request, unique_identifier=unique_identifier))
        self.job_scheduler.flush()
        log.debug("Pipeline stats call home written")

    def submit_jobs(self):
        super(GenPipesPipeline, self).submit_jobs()
        if self.jobs and self.job_scheduler.name.lower() in ["pbs", "batch", "slurm"]:
            self.genpipes_log()



# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Nanopore(GenPipesPipeline):

    def __init__(self, *args, **kwargs):
        super(Nanopore, self).__init__(*args, **kwargs)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
            'SVariants_directory': os.path.relpath(os.path.join(self.output_dir, 'SVariants'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
        }
        return dirs

    @property
    def readsets(self):
        if getattr(self, "_readsets") is None:
            self._readsets = parse_nanopore_readset_file(self.readsets_file)
        return self._readsets


# Abstract pipeline gathering common features of all Illumina sequencing pipelines (trimming, etc.)
# Specific steps must be defined in Illumina children pipelines.
class Illumina(GenPipesPipeline):

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
                    _raise(SanitycheckError(f"""Error: BAM file not available for readset "{readset.name}"!"""))
        return jobs


    def picard_sam_to_fastq(self):
        """
        Convert SAM/BAM files from the input readset file into FASTQ format
        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
        """
        jobs = []
        
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
                        f"{readset.name}.sorted.bam"
                    )
                    candidate_input_files = [
                        [sortedBam],
                        [readset.bam]
                    ]
                    [bam] = self.select_input_files(candidate_input_files)

                    rawReadsDirectory = os.path.join(
                        self.output_dirs['raw_reads_directory'],
                        readset.sample.name,
                    )
                    if readset.run_type == "PAIRED_END":
                        fastq1 = os.path.join(rawReadsDirectory, f"{readset.name}.pair1.fastq.gz")
                        fastq2 = os.path.join(rawReadsDirectory, f"{readset.name}.pair2.fastq.gz")
                    elif readset.run_type == "SINGLE_END":
                        fastq1 = os.path.join(rawReadsDirectory, f"{readset.name}.single.fastq.gz")
                        fastq2 = None
                    else:
                        _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

                    mkdir_job = bash.mkdir(rawReadsDirectory)
                    jobs.append(
                        concat_jobs([
                            mkdir_job,
                            picard.sam_to_fastq(
                                bam,
                                fastq1,
                                fastq2
                                )
                            ],
                            name=f"picard_sam_to_fastq.{readset.name}",
                            samples=[readset.sample],
                            readsets=[readset]
                            )
                        )
                else:
                    _raise(SanitycheckError(f"Error: BAM file not available for readset {readset.name}!"))
        return jobs

    def gatk_sam_to_fastq(self):
        """
        Converts SAM/BAM files from the input readset file into FASTQ format,
        if FASTQ files are not already specified in the readset file. 
        Do nothing otherwise.
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
                        f"{readset.name}.sorted.bam"
                    )
                    candidate_input_files = [
                        [sortedBam],
                        [readset.bam]
                    ]
                    [bam] = self.select_input_files(candidate_input_files)

                    rawReadsDirectory = os.path.join(
                        self.output_dirs['raw_reads_directory'],
                        readset.sample.name,
                    )
                    if readset.run_type == "PAIRED_END":
                        fastq1 = os.path.join(rawReadsDirectory, f"{readset.name}.pair1.fastq.gz")
                        fastq2 = os.path.join(rawReadsDirectory, f"{readset.name}.pair2.fastq.gz")
                    elif readset.run_type == "SINGLE_END":
                        fastq1 = os.path.join(rawReadsDirectory, f"{readset.name}.single.fastq.gz")
                        fastq2 = None
                    else:
                        _raise(SanitycheckError(f"Error: run type {readset.run_type} is invalid for readset {readset.name} (should be PAIRED_END or SINGLE_END)!"))

                    mkdir_job = bash.mkdir(rawReadsDirectory)
                    jobs.append(
                        concat_jobs([
                            mkdir_job,
                            gatk4.sam_to_fastq(
                                bam,
                                fastq1,
                                fastq2
                            )
                        ],
                            name=f"gatk_sam_to_fastq.{readset.name}",
                            samples=[readset.sample],
                            readsets=[readset]
                        ),
                    )
                else:
                    _raise(SanitycheckError(f"Error: BAM file not available for readset {readset.name} !"))
        return jobs

    def trimmomatic(self):
        """
        Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', parameter 'adapter_fasta'),
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
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"
            link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = global_conf.global_get('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
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
                        _raise(SanitycheckError(f"""Error: missing adapter1 and/or adapter2 for PAIRED_END readset "{readset.name}", or missing adapter_fasta parameter in config file!"""))
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Single
{sequence}
END
`""".format(adapter_fasta=adapter_fasta, sequence=readset.adapter1))
                    else:
                        _raise(SanitycheckError(f"""Error: missing adapter1 for SINGLE_END readset "{readset.name}", or missing adapter_fasta parameter in config file!"""))

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
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            if adapter_job:
                job = concat_jobs([adapter_job, job])

            jobs.append(concat_jobs(
                [
                    # Trimmomatic does not create output directory by default
                    bash.mkdir(trim_directory),
                    bash.mkdir(link_directory),
                    job,
                    bash.ln(
                    os.path.relpath(trim_log, link_directory),
                    os.path.join(link_directory, readset.name + ".trim.log"),
                    trim_log
                    )
                ],
                name="trimmomatic." + readset.name,
                samples=[readset.sample],
                readsets=[readset]
                ))
        return jobs

    def merge_trimmomatic_stats(self):
        """
        The trim statistics per readset are merged at this step.
        """

        read_type = "Paired" if self.run_type == 'PAIRED_END' else "Single"
        readset_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimReadsetTable.tsv")
        job = concat_jobs([
            bash.mkdir(self.output_dirs['metrics_directory']),
            Job(command=f"echo 'Sample\\tReadset\\tRaw {read_type} Reads #\\tSurviving {read_type} Reads #\\tSurviving {read_type} Reads %' > {readset_merge_trim_stats}")
            ])
        for readset in self.readsets:
            trim_log = os.path.join(self.output_dirs["trim_directory"], readset.sample.name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                # Retrieve readset raw and surviving reads from trimmomatic log using ugly Perl regexp
                perl_command = f"perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/{readset.sample.name}\\t{readset.name}\\t\\1\\t\\2/'"
            elif readset.run_type == "SINGLE_END":
                perl_command = f"perl -pe 's/^Input Reads: (\\d+).*Surviving: (\\d+).*$/{readset.sample.name}\\t{readset.name}\\t\\1\\t\\2/'"

            job = concat_jobs(
                [
                    job,
                    Job(
                        [trim_log],
                        [readset_merge_trim_stats],
                        module_entries=[['merge_trimmomatic_stats', 'module_perl']],
                        # Create readset trimming stats TSV file with paired or single read count using ugly awk
                        command="""\
grep ^Input {trim_log} | \\
{perl_command} | \\
awk '{{OFS="\\t"; print $0, $4 / $3 * 100}}' \\
>> {readset_merge_trim_stats}""".format(
                            trim_log=trim_log,
                            perl_command=perl_command,
                            readset_merge_trim_stats=readset_merge_trim_stats
                        ),
                        samples=[readset.sample],
                        readsets=[readset]
                    )
                ]
            )

        sample_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimSampleTable.tsv")
        return [
            concat_jobs(
                [
                    job,
                    Job(
                        [readset_merge_trim_stats],
                        [sample_merge_trim_stats],
                        # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                        command="""\
cut -f1,3- {readset_merge_trim_stats} | awk -F"\\t" '{{OFS="\\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Raw Reads #", "Surviving Reads #", "Surviving %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_trim_stats}""".format(
                    readset_merge_trim_stats=readset_merge_trim_stats,
                    sample_merge_trim_stats=sample_merge_trim_stats
                )
            ),
            Job(
                [sample_merge_trim_stats],
                [os.path.join("report", "trimReadsetTable.tsv"), os.path.join("report", "trimSampleTable.tsv")],
                command="""\
mkdir -p report && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} report/""".format(
                            readset_merge_trim_stats=readset_merge_trim_stats,
                            sample_merge_trim_stats=sample_merge_trim_stats
                        ),
                    )
                ],
                name="merge_trimmomatic_stats"
            )
        ]
        
    def build_adapter_file(self, directory, readset):
        adapter_file = os.path.join(directory, "adapter.tsv")
        if readset.run_type == "SINGLE_END":
            if readset.adapter1:
                adapter_job = Job(
                    command="""\
`cat > {adapter_file} << END
>Adapter\n{sequence}\n
END
`""".format(
                        adapter_file=adapter_file,
                        sequence=readset.adapter1
                    )
                )
            else:
                raise Exception(
                    f"Error: missing adapter1 for SINGLE_END readset {readset.name}, or missing adapter_file parameter in config file!")
        
        elif readset.run_type == "PAIRED_END":
            if readset.adapter1 and readset.adapter2:
                adapter_job = Job(
                    command="""\
`cat > {adapter_file} << END
>Adapter1\n{sequence1}\n
>Adapter2\n{sequence2}\n
END
`""".format(
                        adapter_file=adapter_file,
                        sequence1=readset.adapter1,
                        sequence2=readset.adapter2,
                    )
                )
        
        return adapter_job

    def skewer_trimming(self):
        """
        Trimming using [skewer](https://sourceforge.net/projects/skewer/)
        """

        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join(self.output_dirs['trim_directory'], readset.sample.name)
            link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
            
            trim_file_prefix = os.path.join(output_dir, readset.name)
            trim_log = trim_file_prefix + ".log"

            adapter_file = global_conf.global_get('skewer_trimming', 'adapter_file', required=False, param_type='filepath')
            adapter_job = None

            quality_offset = readset.quality_offset

            if not adapter_file:
                adapter_file = os.path.join(output_dir, "adapter.tsv")
                adapter_job = adapters.create(
                    readset,
                    adapter_file
                )

            fastq1 = ""
            fastq2 = ""
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        readset.name + "."
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + ".single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError(f"Error: run type {readset.run_type} is invalid for readset {readset.name} (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            output_dir,
                            remove=True
                        ),
                        bash.mkdir(link_directory),
                        adapter_job,
                        skewer.trim(
                            fastq1,
                            fastq2,
                            trim_file_prefix,
                            adapter_file,
                            quality_offset
                        ),
                        bash.ln(
                            os.path.relpath(trim_log, link_directory),
                            os.path.join(link_directory, readset.name + ".trim.log"),
                            trim_log
                        )
                    ],
                    name="skewer_trimming." + readset.name,
                    removable_files=[output_dir],
                    samples=[readset.sample],
                    readsets=[readset]
                )
            )

        return jobs

    def trim_fastp(self):
        """
        [Fastp](https://github.com/OpenGene/fastp): A tool designed to provide fast all-in-one preprocessing for FastQ
        files. This tool is developed in C++ with multithreading supported to afford high performance.
        """

        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join(self.output_dirs['trim_directory'], readset.sample.name)
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][readset.sample.name])
            
            trim_json = os.path.join(metrics_directory, f"{readset.name}.trim.json")
            trim_html = os.path.join(metrics_directory, f"{readset.name}.trim.html")
            
            adapter_file = global_conf.global_get('trim_fastp', 'adapter_file', required=False, param_type='filepath')
            adapter_job = None

            if not adapter_file:
                adapter_file = os.path.join(output_dir, "adapter.tsv")
                adapter_job = adapters.create(
                    readset,
                    adapter_file
                )

            fastq1 = ""
            fastq2 = ""
            output1 = ""
            output2 = ""
            trim_file_prefix = os.path.join(output_dir, readset.name)
            
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", "", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([f"{prefix}.pair1.fastq.gz", f"{prefix}.pair2.fastq.gz"])
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        readset.name
                    )
                    candidate_input_files.append([f"{prefix}.pair1.fastq.gz", f"{prefix}.pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                output1 = f"{trim_file_prefix}.trim.pair1.fastq.gz"
                output2 = f"{trim_file_prefix}.trim.pair2.fastq.gz"

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", "", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([f"{prefix}.single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                output1 = f"{trim_file_prefix}.trim.single.fastq.gz"
                output2 = None

            else:
                _raise(SanitycheckError(f"Error: run type {readset.run_type} is invalid for readset {readset.name} (should be PAIRED_END or SINGLE_END)!"))

            job_name = f"trim_fastp.{readset.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        fastp.parse_quality_thirty_metrics_pt(trim_json),
                        job2json_project_tracking.run(
                            input_file=trim_json,
                            pipeline=self,
                            samples=readset.sample.name,
                            readsets=readset.name,
                            job_name=job_name,
                            metrics="bases_over_q30_percent=$bases_over_q30_percent"
                        ),
                        fastp.parse_pre_length_r1_metrics(trim_json),
                        job2json_project_tracking.run(
                            input_file=trim_json,
                            pipeline=self,
                            samples=readset.sample.name,
                            readsets=readset.name,
                            job_name=job_name,
                            metrics="pre_mean_length_r1=$pre_mean_length_r1"
                        ),
                        fastp.parse_post_length_r1_metrics(trim_json),
                        job2json_project_tracking.run(
                            input_file=trim_json,
                            pipeline=self,
                            samples=readset.sample.name,
                            readsets=readset.name,
                            job_name=job_name,
                            metrics="post_mean_length_r1=$post_mean_length_r1"
                        ),
                        fastp.parse_pre_length_r2_metrics(trim_json),
                        job2json_project_tracking.run(
                            input_file=trim_json,
                            pipeline=self,
                            samples=readset.sample.name,
                            readsets=readset.name,
                            job_name=job_name,
                            metrics="pre_mean_length_r2=$pre_mean_length_r2"
                        ),
                        fastp.parse_post_length_r2_metrics(trim_json),
                        job2json_project_tracking.run(
                            input_file=trim_json,
                            pipeline=self,
                            samples=readset.sample.name,
                            readsets=readset.name,
                            job_name=job_name,
                            metrics="post_mean_length_r2=$post_mean_length_r2"
                        ),
                    ])

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            output_dir,
                            remove=True
                        ),
                        bash.mkdir(metrics_directory),
                        adapter_job,
                        fastp.trim(
                            adapter_file,
                            fastq1,
                            fastq2,
                            output1,
                            output2,
                            trim_json,
                            trim_html,
                            ini_section='trim_fastp'
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    removable_files=[output_dir],
                    samples=[readset.sample],
                    readsets=[readset],
                    output_dependency=[output1, output2, trim_json, trim_html]
                )
            )
            self.multiqc_inputs[readset.sample.name].append(
                os.path.join(metrics_directory, f"{readset.name}.trim.json")
            )
            
        return jobs

    def bwa_mem2_samtools_sort(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) with algorithm: bwa mem2.
        BWA output BAM files are then sorted by coordinate using [Samtools](https://www.htslib.org/doc/samtools.html)
        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, f"{readset.name}.trim")
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], readset.sample.name)
            compression_postfix = global_conf.global_get("bwa_mem2_samtools_sort", 'compression')
            readset_prefix = os.path.join(alignment_directory, readset.name, readset.name + ".sorted." + compression_postfix)

            fastq1 = ""
            fastq2 = ""
            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[f"{trim_file_prefix}.pair1.fastq.gz", f"{trim_file_prefix}.pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([f"{prefix}.pair1.fastq.gz", f"{prefix}.pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[f"{trim_file_prefix}.single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub(r"\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([f"{prefix}.single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError(f"Error: run type {readset.run_type} is invalid for readset {readset.name} (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(alignment_directory,readset.name)),
                        pipe_jobs(
                            [
                                bwa2.mem(
                                    fastq1,
                                    fastq2,
                                    read_group="'@RG" + \
                                        "\\tID:" + readset.name + \
                                        "\\tSM:" + readset.sample.name + \
                                        "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                        ("\\tPU:" + readset.sample.name + "." + readset.run + "." + readset.lane if readset.sample.name and readset.run and readset.lane else "") + \
                                        ("\\tCN:" + global_conf.global_get('bwa_mem2_samtools_sort', 'sequencing_center') if global_conf.global_get('bwa_mem2_samtools_sort', 'sequencing_center', required=False) else "") + \
                                        ("\\tPL:" + global_conf.global_get('bwa_mem2_samtools_sort', 'sequencing_technology') if global_conf.global_get('bwa_mem2_samtools_sort', 'sequencing_technology', required=False) else "Illumina") + \
                                        "'",
                                        ini_section='bwa_mem2_samtools_sort'
                                ),
                                samtools.sort(
                                    "/dev/stdin",
                                    readset_prefix,
                                    ini_section='samtools_sort'
                                )
                            ]
                        ),
                        samtools.index(
                            readset_prefix,
                            ini_section='samtools_index'
                        )
                    ],
                    name=f"bwa_mem2_samtools_sort.{readset.name}",
                    samples=[readset.sample],
                    readsets=[readset],
                    removable_files=[]
                )
            )

        return jobs

    def sambamba_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
            # Find input readset BAMs first from previous bwa_mem_sambamba_sort_sam job, then from original BAMs in the readset sheet.
            candidate_readset_bams = [
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.UMI.bam") for readset in sample.readsets],
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam") for readset in sample.readsets],
                [readset.bam for readset in sample.readsets if readset.bam]]

            readset_bams = self.select_input_files(candidate_readset_bams)

            sample_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.bam")
            mkdir_job = bash.mkdir(os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                readset_index = re.sub(r"\.bam$", ".bam.bai", readset_bam)
                sample_index = re.sub(r"\.bam$", ".bam.bai", sample_bam)

                if alignment_directory in readset_bam:
                    bam_link = os.path.relpath(readset_bam, alignment_directory)
                    index_link = os.path.relpath(readset_index, alignment_directory)

                else:
                    bam_link = os.path.relpath(readset_bam, os.path.join(self.output_dir, self.output_dirs['alignment_directory'], sample.name))
                    index_link = os.path.relpath(readset_index, os.path.join(self.output_dir, self.output_dirs['alignment_directory'], sample.name))

                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            bash.ln(
                                bam_link,
                                sample_bam,
                                input=readset_bam
                            ),
                            bash.ln(
                                index_link,
                                sample_index,
                                input=readset_index
                            )
                        ],
                        name=f"symlink_readset_sample_bam.{sample.name}",
                        samples=[sample],
                        readsets=list(sample.readsets)
                    )
                )

            # Sambamba merge fails if a file/symlink with the merged sample name already exists. Remove any existing file before merging.
            elif len(sample.readsets) > 1:
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            bash.rm(sample_bam),
                            bash.rm(re.sub(r"\.bam$", ".bam.bai", sample_bam)),
                            sambamba.merge(
                                readset_bams,
                                sample_bam
                            )
                        ],
                        name="sambamba_merge_sam_files." + sample.name,
                        samples=[sample],
                        readsets=list(sample.readsets),
                        input_dependency=readset_bams
                    )
                )

        return jobs

    def samtools_merge_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Samtools](https://www.htslib.org/).

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        compression_postfix = global_conf.global_get('bwa_mem2_samtools_sort', 'compression')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
            # Find input readset BAMs first from previous bwa_mem_sambamba_sort_sam job, then from original BAMs in the readset sheet.
            candidate_readset_bams = [
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.UMI.{compression_postfix}") for readset in
                 sample.readsets],
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.{compression_postfix}") for readset in
                 sample.readsets],
                [readset.bam for readset in sample.readsets if readset.bam]]

            readset_bams = self.select_input_files(candidate_readset_bams)

            sample_bam = os.path.join(alignment_directory, f"{sample.name}.sorted.{compression_postfix}")
            mkdir_job = bash.mkdir(os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if compression_postfix == 'bam':
                    readset_index = re.sub(r"\.bam$", ".bam.bai", readset_bam)
                    sample_index = re.sub(r"\.bam$", ".bam.bai", sample_bam)

                elif compression_postfix == 'cram':
                    readset_index = re.sub(r"\.cram$", ".cram.crai", readset_bam)
                    sample_index = re.sub(r"\.cram$", ".cram.crai", sample_bam)

                if alignment_directory in readset_bam:
                    bam_link = os.path.relpath(readset_bam, alignment_directory)
                    index_link = os.path.relpath(readset_index, alignment_directory)

                else:
                    bam_link = os.path.relpath(readset_bam, os.path.join(self.output_dir, self.output_dirs['alignment_directory'], sample.name))
                    index_link = os.path.relpath(readset_index, os.path.join(self.output_dir, self.output_dirs['alignment_directory'], sample.name))

                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            bash.ln(
                                bam_link,
                                sample_bam,
                                input=readset_bam
                            ),
                            bash.ln(
                                index_link,
                                sample_index,
                                input=readset_index
                            )
                        ],
                        name=f"symlink_readset_sample_bam.{sample.name}",
                        samples=[sample],
                        readsets=list(sample.readsets)
                    )
                )

            elif len(sample.readsets) > 1:
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            samtools.merge(
                                sample_bam,
                                readset_bams,
                                ini_section='samtools_merge_bams'
                            )
                        ],
                        name=f"samtools_merge_sam_files.{sample.name}",
                        samples=[sample],
                        input_dependency=readset_bams
                    )
                )

        return jobs


    def gatk_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [GATK](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            metrics_directory = os.path.join(self.output_dirs["metrics_directory"][sample.name])

            # Find input readset CRAMs/BAMs first from previous bwa_mem2_sambamba_sort_sam job, then from original BAMs in the readset sheet.
            candidate_readset_bams = [
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.UMI.bam") for readset in sample.readsets],
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.cram") for readset in sample.readsets],
                [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam") for readset in sample.readsets],
                [readset.bam for readset in sample.readsets if readset.bam]
            ]
            input = self.select_input_files(candidate_readset_bams)

            if global_conf.global_get("gatk_mark_duplicates", 'compression') == "cram":
                output = os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")
                output_index = f"{output}.crai"
            else:
                output =os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")
                output_index = f"{output}.bai"

            metrics_file = os.path.join(metrics_directory, f"{sample.name}.sorted.dup.metrics")

            job_name = f"gatk_mark_duplicates.{sample.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                        [
                            gatk4.parse_duplicate_rate_metrics_pt(metrics_file),
                            job2json_project_tracking.run(
                                input_file=metrics_file,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="duplication_percent=$duplication_percent"
                            ),
                        ]
                )
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
                        bash.mkdir(metrics_directory),
                        gatk4.mark_duplicates(
                            input,
                            output,
                            metrics_file,
                            remove_duplicates="false",
                            create_index=False,
                            ini_section='gatk_mark_duplicates'
                        ),
                        samtools.index(
                            output,
                            ini_section='samtools_index'
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency= [output, output_index, metrics_file]
                )
            )
            self.multiqc_inputs[sample.name].append(
                metrics_file
            )

        return jobs

    def verify_bam_id(self):
        """
        verifyBamID is a software that verifies whether the reads in particular file match previously known
        genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
        as a mixture of two samples. verifyBamID can detect sample contamination and swaps when external
        genotypes are available. When external genotypes are not available, verifyBamID still robustly
        detects sample swaps.
        """

        # Known variants file
        population_af = global_conf.global_get('verify_bam_id', 'population_AF', required=False)
        candidate_input_files = [[global_conf.global_get('verify_bam_id', 'verifyBamID_variants_file', required=False)],
                                 [global_conf.global_get('verify_bam_id', 'verifyBamID_variants_file', required=False) + ".gz"]]
        [known_variants_annotated] = self.select_input_files(candidate_input_files)
        verify_bam_id_directory = "verify_bam_id"
        variants_directory = "variants"

        jobs = []

        verify_bam_results = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)

            candidate_input_files = [[os.path.join(alignment_directory, f"{sample.name}.sorted.dup.recal.bam")],
                                     [os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")],
                                     [os.path.join(alignment_directory, f"{sample.name}.sorted.mdup.bam")]]
            [input_bam] = self.select_input_files(candidate_input_files)

            output_prefix = os.path.join(verify_bam_id_directory, sample.name)

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            # Run verifyBamID
            jobs.append(concat_jobs(
                [
                    bash.mkdir(variants_directory),
                    bash.mkdir(verify_bam_id_directory),
                    verify_bam_id.verify(
                        input_bam,
                        output_prefix
                    )
                ],
                name = f"verify_bam_id.{sample.name}",
                samples = [sample]
                )
            )

            verify_bam_results.extend([f"{output_prefix}.selfSM" ])

        # Coverage bed is null if whole genome experiment
        target_bed=coverage_bed if coverage_bed else ""

        # Render Rmarkdown Report
        jobs.append(
            rmarkdown.render(
                job_input=verify_bam_results ,
                job_name="verify_bam_id_report",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "Illumina.verify_bam_id.Rmd"),
                samples=self.samples,
                readsets=self.readsets,
                render_output_dir='report',
                module_section='report',
                prerun_r=f'source_dir="{verify_bam_id_directory}"; report_dir="report" ; params=list(verifyBamID_variants_file="{known_variants_annotated}", dbnsfp_af_field="{population_af}", coverage_bed="{target_bed}");'
            )
        )

        return jobs

    def cram_output(self):
        """
        Generate long term storage version of the final alignment files in CRAM format.
        Using this function will add the orginal final bam file to the removable file list.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)

            candidate_input_files = [
                [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.recal.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.dedup.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.mdup.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.merged.mdup.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.matefixed.sorted.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.realigned.sorted.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.merged.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted_noRG.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.UMI.bam")],
                [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]

            [input_bam] = self.select_input_files(candidate_input_files)

            output_cram = re.sub(r"\.bam$", ".cram", input_bam)

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
                    options=global_conf.global_get('samtools_cram_output', 'options'),
                    removable=False
                )

            job.name = f"samtools_cram_output.{sample.name}"
            job.samples = [sample]
            job.readsets = [*list(sample.readsets)]
            job.removable_files = [input_bam, f"{input_bam}.bai"]

            jobs.append(job)

        return jobs


class Error(Exception):
    """
    base error for common pipelines class
    """
    pass

class MissingInputError(Error):
    """
    Error for missing input files
    """
    pass
