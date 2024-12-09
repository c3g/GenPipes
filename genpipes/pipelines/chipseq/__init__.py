################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPIpes.
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
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

import csv
import logging
import os
import re

# GenPipes Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs, pipe_jobs
from ...core.design import parse_chipseq_design_file
from ...core.readset import parse_illumina_readset_file
from .. import common

from ...bfx.sequence_dictionary import parse_sequence_dictionary_file

from ...bfx import(
    bash_cmd as bash,
    bedtools,
    bwa2,
    differential_binding,
    gatk4,
    homer,
    macs2,
    multiqc,
    picard,
    sambamba,
    samtools,
    tools,
    trimmomatic,
    ucsc
    )

log = logging.getLogger(__name__)

class ChipSeq(common.Illumina):
    """
ChIP-Seq Pipeline
=================

A pipeline to process ChIP-seq data. The pipeline is designed to handle both paired-end and single-end reads and can be used to process data from any Illumina sequencer. The pipeline uses Trimmomatic to trim reads and remove Illumina adapters, BWA to align reads to the reference genome, Sambamba to sort and filter BAM files, and Picard to mark duplicates and collect quality metrics. The pipeline also uses MACS2 to call peaks, HOMER to annotate peaks, and DiffBind to perform differential binding analysis.

The pipeline takes as input a readset file and a design file. The readset file contains the list of samples and readsets, while the design file contains the list of contrasts to be analyzed. The pipeline outputs BAM files, peak calls, and differential binding results. The pipeline also generates quality metrics and reports for each sample.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.

Attributes:
    output_dirs (dict): Output directory paths
    mark_type_conversion (dict): Conversion of mark type codes to mark type names
    ucsc_genome (str): UCSC genome assembly name
    readsets (list): List of readsets, making sure that the mark name and mark type are defined
    multiqc_inputs (list): List of MultiQC input files
    contrasts (list): List of contrasts
    sequence_dictionary_variant (dict): Sequence dictionary variant
    mappable_genome_size (int): Mappable genome size
Methods:
    trimmomatic: Raw reads quality trimming and removing of Illumina adapters
    merge_trimmomatic_stats: Merge trim statistics per readset
    mapping_bwa_mem_sambamba: Align reads to the reference genome
    sambamba_merge_bam_files: Merge BAM readset files into one file per sample
    sambamba_mark_duplicates: Mark duplicates
    sambamba_view_filter: Filter out unmapped reads and low quality reads
    bedtools_blacklist_filter: Remove reads in blacklist regions from BAM
    metrics: Compute the number of raw/filtered and aligned reads per sample
    macs2_callpeak: Call peaks
    homer_annotate_peaks: Annotate peaks
    diffbind: Perform differential binding analysis
    multiqc: Generate a MultiQC report
Parameters:
    readsets_file (str): Readset file path
    design_file (str): Design file path
    output_dir (str): Output directory path
    run_type (str): Run type (default PAIRED_END)
    samples (list): List of samples
    protocol (str): Type of pipeline (default chipseq)
    """

    def __init__(self, *args, protocol="chipseq", **kwargs):
        self.protocol = protocol
        # Add pipeline specific arguments
        super(ChipSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument(
            "-t",
            "--type",
            help="Type of pipeline (default chipseq)",
            choices=["chipseq", "atacseq"],
            default="chipseq",
            dest='protocol'
            )
        return cls._argparser

    @property
    def output_dirs(self):
        """
        Output directory paths.
        Returns:
            dict: Output directory paths.
        """
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'alignment_output_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'report_output_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'homer_output_directory': os.path.relpath(os.path.join(self.output_dir, 'tags'), self.output_dir),
            'graphs_output_directory': os.path.relpath(os.path.join(self.output_dir, 'graphs'), self.output_dir),
            'tracks_output_directory': os.path.relpath(os.path.join(self.output_dir, 'tracks'), self.output_dir),
            'macs_output_directory': os.path.relpath(os.path.join(self.output_dir, 'peak_call'), self.output_dir),
            'anno_output_directory': os.path.relpath(os.path.join(self.output_dir, 'annotation'), self.output_dir),
            'ihecM_output_directory': os.path.relpath(os.path.join(self.output_dir, 'ihec_metrics'), self.output_dir),
            'dba_output_directory': os.path.relpath(os.path.join(self.output_dir, 'differential_binding'), self.output_dir)
        }
        return dirs

    @property
    def mark_type_conversion(self):
        """
        Conversion of mark type codes to mark type names.
        Returns:
            dict: Conversion of mark type codes to mark type names.
        """
        dirs = {
            'N': 'narrow',
            'B': 'broad',
            'I': 'Input'
        }
        return dirs

    @property
    def ucsc_genome(self):
        """
        UCSC genome assembly name.
        Returns:
            str: UCSC genome assembly name.
        """
        genome_source = global_conf.global_get('DEFAULT', 'source')
        if genome_source == "UCSC":
            genome = global_conf.global_get('DEFAULT', 'assembly')
        else:
            genome = global_conf.global_get('DEFAULT', 'assembly_synonyms')
        return genome

    @property
    def readsets(self):
        """
        List of readsets, making sure that the mark name and mark type are defined.
        Returns:
            list: List of readsets.
        """
        if getattr(self, "_readsets") is None:
            self._readsets = parse_illumina_readset_file(self.readsets_file)
            for readset in self.readsets:
                if not readset.mark_name:
                    _raise(SanitycheckError("Error: missing readset MarkName for " + readset.name))
                elif not readset.mark_type:
                    _raise(SanitycheckError("Error: missing readset MarkType for " + readset.name))
        return self._readsets

    @property
    def multiqc_inputs(self):
        """
        List of MultiQC input files.
        Returns:
            list: List of MultiQC input files.
        """
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    @property
    def contrasts(self):
        """
        List of contrasts.
        Returns:
            list: List of contrasts.
        """
        if self.design_file:
            self._contrast = parse_chipseq_design_file(self.design_file.name, self.samples)
        else:
            self.argparser.error("argument -d/--design is required!")

        return self._contrast

    def sequence_dictionary_variant(self):
        """
        Sequence dictionary variant.
        Returns:
            dict: Sequence dictionary variant.
        """
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(global_conf.global_get('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=True)
        return self._sequence_dictionary_variant

    def mappable_genome_size(self):
        """
        Mappable genome size.
        Returns:
            int: Mappable genome size.
        """
        genome_index = csv.reader(open(global_conf.global_get('DEFAULT', 'genome_fasta', param_type='filepath') + ".fai", 'r'), delimiter='\t')
        # 2nd column of genome index contains chromosome length
        # HOMER and MACS2 mappable genome size (without repetitive features) is about 80 % of total size
        return int(sum(int(chromosome[1]) for chromosome in genome_index) * global_conf.global_get('DEFAULT', 'mappable_genome_size', param_type='float', required=True))

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

        Returns:
            list: List of Trimmomatic jobs.
        """
        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for readset in self.readsets:
            # log.info(readset.mark_name)
            trim_directory = os.path.join("trim", readset.sample.name, readset.mark_name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = global_conf.global_get('trimmomatic', 'adapter_fasta', required=False, param_type='filepath')
            adapter_job = None
            if not adapter_fasta:
                adapter_fasta = trim_file_prefix + "adapters.fa"
                if readset.run_type == "PAIRED_END":
                    if readset.adapter1 and readset.adapter2:
                        # WARNING: Reverse-complement and swap readset adapters for Trimmomatic Palindrome strategy
                        adapter_job = Job(command=f"""\
`cat > {adapter_fasta} << END
>Prefix/1
{readset.adapter2.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]}
>Prefix/2
{readset.adapter1.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]}
END
`"""
                        )
                    else:
                        _raise(SanitycheckError("Error: missing adapter1 and/or adapter2 for PAIRED_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command=f"""\
`cat > {adapter_fasta} << END
>Single
{readset.adapter1}
END
`"""
                        )
                    else:
                        _raise(SanitycheckError("Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!"))

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_fastq1 = os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair1.fastq.gz")
                    candidate_fastq2 = os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair2.fastq.gz")
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
                    candidate_input_files.append([os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".single.fastq.gz")])
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
                bash.mkdir(trim_directory),
                bash.mkdir(link_directory),
                job,
                bash.ln(
                    os.path.relpath(trim_log, link_directory),
                    os.path.join(link_directory, readset.name + ".trim.log"),
                    input = trim_log
                )
            ], name="trimmomatic." + readset.name, samples=[readset.sample]))

            self.multiqc_inputs.append(os.path.join(link_directory, readset.name + ".trim.log"))
        return jobs

    def merge_trimmomatic_stats(self):
        """
        The trim statistics per readset are merged at this step.
        Returns:
            list: List of merge_trimmomatic_stats jobs.
        """

        read_type = "Paired" if self.run_type == 'PAIRED_END' else "Single"
        readset_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimReadsetTable.tsv")
        job = concat_jobs([
            bash.mkdir(self.output_dirs['metrics_directory']),
            Job(
                command=f"""
echo -e "Sample\\tReadset\\tMark Name\\tRaw {read_type} Reads #\\tSurviving {read_type} Reads #\\tSurviving {read_type} Reads %" > {readset_merge_trim_stats}"""
            )
        ])

        for readset in self.readsets:
            trim_log = os.path.join("trim", readset.sample.name, readset.mark_name, readset.name + ".trim.log")
            if readset.run_type == "PAIRED_END":
                # Retrieve readset raw and surviving reads from trimmomatic log using ugly Perl regexp
                perl_command = "perl -pe 's/^Input Read Pairs: (\\d+).*Both Surviving: (\\d+).*Forward Only Surviving: (\\d+).*$/{readset.sample.name}\\t{readset.name}\\t{readset.mark_name}\\t\\1\\t\\2/'".format(
                    readset=readset)
            elif readset.run_type == "SINGLE_END":
                perl_command = "perl -pe 's/^Input Reads: (\\d+).*Surviving: (\\d+).*$/{readset.sample.name}\\t{readset.name}\\t{readset.mark_name}\\t\\1\\t\\2/'".format(
                    readset=readset)

            job = concat_jobs([
                job,
                Job(
                    [trim_log],
                    [readset_merge_trim_stats],
                    module_entries=[['merge_trimmomatic_stats', 'module_perl']],
                    # Create readset trimming stats TSV file with paired or single read count using ugly awk
                    command=f"""\
grep ^Input {trim_log} | \\
{perl_command} | \\
awk '{{OFS="\\t"; print $0, $5 / $4 * 100}}' \\
  >> {readset_merge_trim_stats}""",
                    samples=[readset.sample]
                )
            ])

        sample_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimSampleTable.tsv")
        return [concat_jobs([
            job,
            Job(
                [readset_merge_trim_stats],
                [sample_merge_trim_stats],
                # Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
                command=f"""\
cut -f1,3- {readset_merge_trim_stats} | awk -F"\\t" '{{OFS="\\t"; if (NR==1) {{if ($3=="Raw Paired Reads #") {{paired=1}};print "Sample", "Mark Name", "Raw Reads #", "Surviving Reads #", "Surviving %"}} else {{if (paired) {{$3=$3*2; $4=$4*2}}; sample[$1$2]=$1; markname[$1$2]=$2; raw[$1$2]+=$3; surviving[$1$2]+=$4}}}}END{{for (samplemark in raw){{print sample[samplemark], markname[samplemark], raw[samplemark], surviving[samplemark], surviving[samplemark] / raw[samplemark] * 100}}}}' \\
  > {sample_merge_trim_stats} && \\
mkdir -p {self.output_dirs['report_output_directory']} && \\
cp {readset_merge_trim_stats} {sample_merge_trim_stats} {self.output_dirs['report_output_directory']}/"""
            )
            ],
            name="merge_trimmomatic_stats"
            )
        ]

    def mapping_bwa_mem_sambamba(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem2.
        BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        This step takes as input files:
        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

        Returns:
            list: List of mapping_bwa_mem_sambamba jobs.
        """

        jobs = []
        sequencing_center = global_conf.global_get('mapping_bwa_mem_sambamba', 'sequencing_center', required=False)
        sequencing_technology = global_conf.global_get('mapping_bwa_mem_sambamba', 'sequencing_technology') if global_conf.global_get('mapping_bwa_mem_sambamba', 'sequencing_technology', required=False) else "Illumina"
        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.sample.name, readset.mark_name)
            trim_file_prefix = os.path.join(trim_directory, readset.name)
            alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], readset.sample.name, readset.mark_name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam.bai")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.pair1.fastq.gz", trim_file_prefix + ".trim.pair2.fastq.gz"]
                ]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [
                    [trim_file_prefix + ".trim.single.fastq.gz"]
                ]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            jobs.append(
                concat_jobs([
                    bash.mkdir(os.path.dirname(readset_bam)),
                    pipe_jobs([
                        bwa2.mem(
                            fastq1,
                            fastq2,
                            read_group="'@RG" + \
                                       f"\\tID:{readset.name}" + \
                                       f"\\tSM:{readset.sample.name}" + \
                                       "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                                       ("\\tPU:run" + f"{readset.run}_{readset.lane}" if readset.run and readset.lane else "") + \
                                       (f"\\tCN:{sequencing_center}") + \
                                       (f"\\tPL:{sequencing_technology}") + \
                                       "'",
                            ini_section='mapping_bwa_mem_sambamba'
                        ),
                        sambamba.view(
                            "/dev/stdin",
                            None,
                            options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_view_other_options')
                        ),
                        sambamba.sort(
                            "/dev/stdin",
                            readset_bam,
                            tmp_dir=global_conf.global_get('mapping_bwa_mem_sambamba', 'tmp_dir', required=True),
                            other_options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_sort_other_options', required=False)
                        )
                    ]),
                    sambamba.index(
                        readset_bam,
                        index_bam,
                        other_options=global_conf.global_get('mapping_bwa_mem_sambamba', 'sambamba_index_other_options', required=False)
                    )
                ],
                    name=f"mapping_bwa_mem_sambamba.{readset.name}",
                    samples=[readset.sample]
                )
            )

        return jobs


    def sambamba_merge_bam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        This step takes as input files:
        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file

        Returns:
            list: List of sambamba_merge_bam_files jobs.
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name)
                # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
                readset_bams = [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for
                                readset in sample.readsets if readset.mark_name == mark_name]
                sample_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.bam")

                # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
                if len(readset_bams) == 1:
                    readset_bam = readset_bams[0]
                    readset_index = re.sub(r"\.bam$", ".bam.bai", readset_bam)
                    sample_index = re.sub(r"\.bam$", ".bam.bai", sample_bam)

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    os.path.dirname(sample_bam)
                                ),
                                bash.ln(
                                    os.path.relpath(readset_bam, os.path.dirname(sample_bam)),
                                    sample_bam,
                                    input=readset_bam
                                ),
                                bash.ln(
                                    os.path.relpath(readset_index, os.path.dirname(sample_index)),
                                    sample_index,
                                    input=readset_index
                                )
                            ],
                            name=f"symlink_readset_sample_bam.{sample.name}.{mark_name}",
                            samples=[sample]
                        )
                    )

                elif len(sample.readsets) > 1:
                    jobs.append(
                        concat_jobs([
                            bash.mkdir(
                                os.path.dirname(sample_bam)
                            ),
                            sambamba.merge(
                                readset_bams,
                                sample_bam,
                                ini_section="sambamba_merge_bam_files"
                            )
                        ],
                            name=f"sambamba_merge_bam_files.{sample.name}.{mark_name}",
                            samples=[sample]
                        )
                    )
        return jobs


    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        Returns:
            list: List of sambamba_mark_duplicates jobs.
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                                   mark_name)
                input_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.bam")
                output_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.bam")

                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_bam)),
                        sambamba.markdup(
                            input_bam,
                            output_bam,
                            tmp_dir=global_conf.global_get('sambamba_mark_duplicates', 'tmp_dir', required=True),
                            other_options=global_conf.global_get('sambamba_mark_duplicates', 'other_options', required=False)
                        )
                    ],
                        name=f"sambamba_mark_duplicates.{sample.name}.{mark_name}",
                        samples=[sample]
                    )
                )

        return jobs

    def sambamba_view_filter(self):
        """
        Filter out unmapped reads and low quality reads [Sambamba](http://lomereiter.github.io/sambamba/index.html).

        Returns:
            list: List of sambamba_view_filter jobs.
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name)
                input_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.bam")
                output_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                output_bam_index = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam.bai")
                jobs.append(
                    concat_jobs([
                        bash.mkdir(os.path.dirname(output_bam)),
                        sambamba.view(
                            input_bam,
                            output_bam,
                            f"""-t {global_conf.global_get('sambamba_view_filter', 'threads')} -f bam -F \"not unmapped and not failed_quality_control and mapping_quality >= {global_conf.global_get('sambamba_view_filter', 'min_mapq')}\""""
                        ),
                        sambamba.index(
                            output_bam,
                            output_bam_index
                        )
                    ],
                        name=f"sambamba_view_filter.{sample.name}.{mark_name}",
                        samples=[sample]
                    )
                )

        return jobs

    def bedtools_blacklist_filter(self):
        """
        Remove reads in blacklist regions from bam with bedtools intersect if blacklist file is supplied. Do nothing otherwise.
        Returns:
            list: List of bedtools_blacklist_filter jobs.
        """
        jobs = []

        if global_conf.global_get('bedtools_intersect', 'blacklist', required=False, param_type='filepath'):
            blacklist = global_conf.global_get('bedtools_intersect', 'blacklist', param_type='filepath')

            for sample in self.samples:
                for mark_name in sample.marks:
                    alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name)
                    input_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                    output_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")
                    output_bam_index = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam.bai")

                    jobs.append(
                            concat_jobs(
                                [
                                    bedtools.intersect(
                                        input_bam,
                                        output_bam,
                                        blacklist,
                                        include_header=True
                                        ),
                                    sambamba.index(
                                        output_bam,
                                        output_bam_index
                                        )
                                    ],
                                name=f"bedtools_intersect.{sample.name}.{mark_name}",
                                samples=[sample]
                                )
                            )
        return jobs


    def metrics(self):
        """
        The number of raw/filtered and aligned reads per sample are computed at this stage.
        Returns:
            list: List of metrics jobs.
        """

        jobs = []

        samples_associative_array = []
        inputs_report = []

        metrics_output_directory = self.output_dirs['metrics_directory']
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"")
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name)
                raw_bam_file = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.bam")
                # Select input from blacklist filtered (clean) or just sambamba filtered bam
                filtered_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                clean_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")
                candidate_input_files = [[clean_bam], [filtered_bam]]
                [bam_file] = self.select_input_files(candidate_input_files)
                picard_metrics_output = os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"bam$", "all.metrics", os.path.basename(bam_file)))

                picard_multiple_metrics_job = concat_jobs(
                    [
                        bash.mkdir(os.path.join(metrics_output_directory, sample.name, mark_name)),
                        bash.mkdir(link_directory),
                        picard.collect_multiple_metrics(
                            bam_file,
                            picard_metrics_output,
                            library_type=self.run_type
                        )
                    ]
                )

                for outfile in picard_multiple_metrics_job.report_files:
                    self.multiqc_inputs.append(outfile)
                    picard_multiple_metrics_job = concat_jobs(
                        [
                            picard_multiple_metrics_job,
                            bash.ln(
                                os.path.relpath(outfile, link_directory),
                                os.path.join(link_directory, os.path.basename(outfile)),
                                input=outfile
                            )
                        ]
                    )

                picard_multiple_metrics_job.name = f"picard_collect_multiple_metrics.{sample.name}.{mark_name}"
                picard_multiple_metrics_job.samples=[sample]
                jobs.append(picard_multiple_metrics_job)

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.join(metrics_output_directory, sample.name, mark_name)),
                            sambamba.flagstat(
                                raw_bam_file,
                                os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file)))
                            ),
                            sambamba.flagstat(
                                bam_file,
                                os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(bam_file)))
                            ),
                            bash.ln(
                                os.path.relpath(os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file))), link_directory),
                                os.path.join(link_directory, re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file))),
                                input = os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file)))
                            )
                        ],
                        name=f"metrics_flagstat.{sample.name}.{mark_name}",
                        samples=[sample]
                    )
                )
                self.multiqc_inputs.append(os.path.join(link_directory,  re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file))))
                inputs_report.extend(
                    [
                        os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(raw_bam_file))),
                        os.path.join(metrics_output_directory, sample.name, mark_name, re.sub(r"\.bam$", ".flagstat", os.path.basename(bam_file))),
                        bam_file
                    ]
                )

        trim_metrics_file = os.path.join(metrics_output_directory, "trimSampleTable.tsv")
        metrics_file = os.path.join(metrics_output_directory, "SampleMetrics.tsv")
        report_metrics_file = os.path.join(self.output_dirs['report_output_directory'], "SampleMetrics.tsv")
        if global_conf.global_get('bedtools_intersect', 'blacklist', required=False, param_type='filepath'):
            bam_ext = "sorted.dup.filtered.cleaned.bam"
        else:
            bam_ext = "sorted.dup.filtered.bam"
        flagstat_ext = re.sub(r"\.bam", "", bam_ext)
        
        jobs.append(
            Job(
                inputs_report,
                [report_metrics_file],
                [
                    ['metrics', 'module_sambamba'],
                    ['metrics', 'module_samtools']
                ],
                # Retrieve number of aligned and duplicate reads from sample flagstat files
                # Merge trimming stats per sample with aligned and duplicate stats using ugly awk
                command="""\
mkdir -p {metrics_dir}
cp /dev/null {metrics_file} && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    raw_flagstat_file={metrics_dir}/$sample/$mark_name/$sample.$mark_name.sorted.dup.flagstat
    filtered_flagstat_file={metrics_dir}/$sample/$mark_name/$sample.$mark_name.{flagstat_ext}.flagstat
    bam_file={alignment_dir}/$sample/$mark_name/$sample.$mark_name.{bam_ext}
    raw_supplementarysecondary_reads=`bc <<< $(grep "secondary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    mapped_reads=`bc <<< $(grep "mapped (" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$raw_supplementarysecondary_reads`
    filtered_supplementarysecondary_reads=`bc <<< $(grep "secondary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* secondary.*//')+$(grep "supplementary" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* supplementary.*//')`
    filtered_reads=`bc <<< $(grep "in total" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_reads=`bc <<< $(grep "mapped (" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* mapped (.*)//')-$filtered_supplementarysecondary_reads`
    filtered_mapped_rate=`echo "scale=4; 100*$filtered_mapped_reads/$filtered_reads" | bc -l`
    filtered_dup_reads=`grep "duplicates" $filtered_flagstat_file | sed -e 's/ + [[:digit:]]* duplicates$//'`
    filtered_dup_rate=`echo "scale=4; 100*$filtered_dup_reads/$filtered_mapped_reads" | bc -l`
    filtered_dedup_reads=`echo "$filtered_mapped_reads-$filtered_dup_reads" | bc -l`
    if [[ -s {trim_metrics_file} ]]
      then
        raw_reads=$(grep -P "${{sample}}\\t${{mark_name}}" {trim_metrics_file} | cut -f 3)
        raw_trimmed_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_trimmed_reads" | bc -l`
        raw_trimmed_rate=`echo "scale=4; 100*$raw_trimmed_reads/$raw_reads" | bc -l`
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_trimmed_reads" | bc -l`
      else
        raw_reads=`bc <<< $(grep "in total" $raw_flagstat_file | sed -e 's/ + [[:digit:]]* in total .*//')-$raw_supplementarysecondary_reads`
        raw_trimmed_reads="NULL"
        mapped_reads_rate=`echo "scale=4; 100*$mapped_reads/$raw_reads" | bc -l`
        raw_trimmed_rate="NULL"
        filtered_rate=`echo "scale=4; 100*$filtered_reads/$raw_reads" | bc -l`
    fi
    filtered_mito_reads=$(sambamba view -F "not duplicate" -c $bam_file chrM)
    filtered_mito_rate=$(echo "scale=4; 100*$filtered_mito_reads/$filtered_mapped_reads" | bc -l)
    if [[ $mark_name != "Input" ]]
        then
          chip_bed_file={macs_directory}/$sample/$mark_name/$sample.${{mark_name}}_peaks.*Peak.bed
          nmb_peaks=$(wc -l $chip_bed_file | cut -f 1 -d " ")
          reads_under_peaks=$(samtools view -c -L $chip_bed_file $bam_file)
          frip=$(echo "scale=4; $reads_under_peaks/$filtered_mapped_reads" | bc -l)
        else
          nmb_peaks="NA"
          reads_under_peaks="NA"
          frip="NA"
    fi
    echo -e "$sample\\t$mark_name\\t$raw_reads\\t$raw_trimmed_reads\\t$raw_trimmed_rate\\t$mapped_reads\\t$mapped_reads_rate\\t$filtered_reads\\t$filtered_rate\\t$filtered_mapped_reads\\t$filtered_mapped_rate\\t$filtered_dup_reads\\t$filtered_dup_rate\\t$filtered_dedup_reads\\t$filtered_mito_reads\\t$filtered_mito_rate\\t$nmb_peaks\\t$reads_under_peaks\\t$frip" >> {metrics_file}
  done
done && \\
sed -i -e "1 i\\Sample\\tMark Name\\tRaw Reads #\\tRemaining Reads after Trimming #\\tRemaining Reads after Trimming %\\tAligned Trimmed Reads #\\tAligned Trimmed Reads %\\tRemaining Reads after Filtering #\\tRemaining Reads after Filtering %\\tAligned Filtered Reads #\\tAligned Filtered Reads %\\tDuplicate Reads #\\tDuplicate Reads %\\tFinal Aligned Reads # without Duplicates\\tMitochondrial Reads #\\tMitochondrial Reads %\\tNumber of Peaks\\tReads under Peaks\\tFRIP" {metrics_file} && \\
mkdir -p {report_dir} && \\
cp {metrics_file} {report_metrics_file}""".format(
                    sambamba=global_conf.global_get('DEFAULT', 'module_sambamba'),
                    metrics_dir=metrics_output_directory,
                    macs_directory=self.output_dirs['macs_output_directory'],
                    metrics_file=metrics_file,
                    samples_associative_array=" ".join(samples_associative_array),
                    alignment_dir=self.output_dirs['alignment_output_directory'],
                    flagstat_ext=flagstat_ext,
                    bam_ext=bam_ext,
                    report_dir=self.output_dirs['report_output_directory'],
                    trim_metrics_file=trim_metrics_file,
                    report_metrics_file=report_metrics_file
                ),
                name="metrics_report",
                samples=self.samples,
                removable_files=[report_metrics_file]
            )
        )

        return jobs


    def homer_make_tag_directory(self):
        """
        The Homer Tag directories, used to check for quality metrics, are computed at this step.
        Returns:
            list: List of homer_make_tag_directory jobs.
        """

        jobs = []
        for sample in self.samples:
            for mark_name in sample.marks:
                # add input selection that allows for use of blacklist-filtered and unfiltered bams
                filtered_bam = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                cleaned_bam = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")
                candidate_input_files = [[cleaned_bam], [filtered_bam]]
                [alignment_file] = self.select_input_files(candidate_input_files)
                output_dir = os.path.join(self.output_dirs['homer_output_directory'], sample.name, f"{sample.name}.{mark_name}")
                other_options = global_conf.global_get('homer_make_tag_directory', 'other_options', required=False)
                genome = global_conf.global_get('homer_make_tag_directory', 'genome', required=False) if global_conf.global_get('homer_make_tag_directory', 'genome', required=False) else self.ucsc_genome
                link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", f"{sample.name}.{mark_name}")

                job = concat_jobs(
                    [
                        bash.mkdir(link_directory),
                        homer.makeTagDir(
                            output_dir,
                            alignment_file,
                            genome,
                            restriction_site=None,
                            illuminaPE=False,
                            other_options=other_options
                        )
                    ]
                )

                for outfile in job.report_files:
                    self.multiqc_inputs.append(outfile)
                    job = concat_jobs(
                            [
                                job,
                                bash.ln(
                                    os.path.relpath(outfile, link_directory),
                                    os.path.join(link_directory, os.path.basename(outfile)),
                                    input=outfile
                                    )
                                ]
                            )

                job.name = f"homer_make_tag_directory.{sample.name}.{mark_name}"
                job.removable_files = [output_dir]
                jobs.append(job)

        return jobs

    def qc_metrics(self):
        """
        Sequencing quality metrics as tag count, tag autocorrelation, sequence bias and GC bias are generated.
        Returns:
            list: List of qc_metrics jobs.
        """

        readset_file = os.path.relpath(self.readsets_file.name, self.output_dir)

        output_files = [os.path.join(self.output_dirs['graphs_output_directory'], f"{sample.name}.{mark_name}_QC_Metrics.ps") for sample in self.samples for mark_name in sample.marks]

        jobs = []

        jobs.append(
            Job(
                [os.path.join(self.output_dirs['homer_output_directory'], sample.name, f"{sample.name}.{mark_name}", "tagInfo.txt") for sample in self.samples for mark_name in sample.marks],
                output_files,
                [
                    ['qc_plots_R', 'module_mugqic_tools'],
                    ['qc_plots_R', 'module_R']
                ],
                command="""\
mkdir -p {graphs_dir} && \\
Rscript $R_TOOLS/chipSeqGenerateQCMetrics.R \\
  {readset_file} \\
  {output_dir} && \\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    cp --parents {graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.ps {report_dir}/
    convert -rotate 90 {graphs_dir}/${{sample}}.${{mark_name}}_QC_Metrics.ps {report_dir}/graphs/${{sample}}.${{mark_name}}_QC_Metrics.png
  done
done""".format(
                    samples_associative_array=" ".join(
                        ["[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"" for sample in
                         self.samples]),
                    readset_file=readset_file,
                    output_dir=self.output_dir,
                    report_dir=self.output_dirs['report_output_directory'],
                    graphs_dir=self.output_dirs['graphs_output_directory']
                ),
                name="qc_plots_R",
                samples=self.samples,
                removable_files=output_files
            )
        )
        return jobs

    def homer_make_ucsc_file(self):
        """
        Wiggle Track Format files are generated from the aligned reads using [Homer](http://homer.ucsd.edu/homer/index.html).
        The resulting files can be loaded in browsers like IGV or UCSC.
        Returns:
            list: List of homer_make_ucsc_file jobs.
        """

        jobs = []

        for sample in self.samples:
            for mark_name in sample.marks:
                tag_dir = os.path.join(self.output_dirs['homer_output_directory'], sample.name, f"{sample.name}.{mark_name}")
                bedgraph_dir = os.path.join(self.output_dirs['tracks_output_directory'], sample.name, mark_name)
                bedgraph_file = os.path.join(bedgraph_dir, f"{sample.name}.{mark_name}.ucsc.bedGraph")
                big_wig_output = os.path.join(bedgraph_dir, "bigWig", f"{sample.name}.{mark_name}.bw")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(bedgraph_dir),
                            homer.makeUCSCfile(
                                tag_dir,
                                bedgraph_file
                            )
                        ],
                        name=f"homer_make_ucsc_file.{sample.name}.{mark_name}",
                        samples=[sample],
                        removable_files=[bedgraph_dir]
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(os.path.join(bedgraph_dir, "bigWig")),
                            ucsc.bedGraphToBigWig(
                                bedgraph_file,
                                big_wig_output,
                                header=True,
                                ini_section="homer_make_ucsc_file")
                        ],
                        name=f"homer_make_ucsc_file_bigWig.{sample.name}.{mark_name}",
                        samples=[sample]
                    )
                )

        jobs.append(
            Job(
                [
                    os.path.join(self.output_dirs['tracks_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.ucsc.bedGraph.gz")
                        for sample in self.samples
                            for mark_name in sample.marks
                ],
                [os.path.join(self.output_dirs['report_output_directory'], "tracks.zip")],
                command="""\
mkdir -p {report_dir} && \\
zip -r {report_dir}/tracks.zip {tracks_dir}/*/*/*.ucsc.bedGraph.gz""".format(
                    report_dir=self.output_dirs['report_output_directory'],
                    tracks_dir=self.output_dirs['tracks_output_directory']
                ),
                name="homer_make_ucsc_file_zip"
            )
        )
        return jobs

    def macs2_callpeak(self):
        """
        Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
        The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
        The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
        The default mfold parameter of MACS2 is [10,30].
        Returns:
            list: List of macs2_callpeak jobs.
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            mark_list = []
            # if no Input file
            input_file = []
            if global_conf.global_get('bedtools_intersect', 'blacklist', required=False, param_type='filepath'):
                input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam") for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            else:
                input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam") for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            if len(input_file_list) > 0:
                if len(input_file_list) > 1:
                    raise Exception(f"""Error: Sample "{sample.name}" has more than 1 Input!""")
                input_file = [input_file_list[0]]
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    filtered_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")]
                    cleaned_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")]
                    candidate_input_files = [cleaned_file, filtered_file]
                    mark_file = self.select_input_files(candidate_input_files)
                    output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    ## set macs2 variables:

                    options = "--format " + ("BAMPE" if self.run_type == "PAIRED_END" else "BAM")
                    genome_size = global_conf.global_get('macs2_callpeak', 'genome_size', required=False) if global_conf.global_get('macs2_callpeak', 'genome_size', required=False) else self.mappable_genome_size()
                    output_prefix_name = os.path.join(output_dir, f"{sample.name}.{mark_name}")

                    if mark_type == "B":  # Broad region
                        options += " --broad --nomodel"
                    else:  # Narrow region
                        if input_file:
                            options += " --nomodel"
                        else:
                            options += " --fix-bimodal"

                    options += " --shift " + global_conf.global_get('macs2_callpeak', 'shift') if global_conf.global_get('macs2_callpeak', 'shift') else ""
                    options += " --extsize " + global_conf.global_get('macs2_callpeak', 'extsize') if global_conf.global_get('macs2_callpeak', 'extsize') else ""
                    options += " -p " + global_conf.global_get('macs2_callpeak', 'pvalue') if global_conf.global_get('macs2_callpeak', 'pvalue') else ""
                    output = []
                    peak_file = f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak"
                    peak_bed_file = f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak.bed"
                    output.append(os.path.join(output_dir, peak_file))
                    output.append(os.path.join(output_dir, f"{sample.name}.{mark_name}_peaks.xls"))


                    jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            bash.mkdir(link_directory),
                            macs2.callpeak(
                                options,
                                genome_size,
                                mark_file,
                                input_file,
                                output_prefix_name,
                                output
                            ),
                            bash.ln(
                                os.path.relpath(output[1], link_directory),
                                os.path.join(link_directory, f"{sample.name}.{mark_name}_peaks.xls"),
                                input = output[1]
                            )
                        ],
                            name=f"macs2_callpeak.{sample.name}.{mark_name}",
                            removable_files=[output_dir]
                        )
                    )
                    self.multiqc_inputs.append(os.path.join(link_directory, f"{sample.name}.{mark_name}_peaks.xls"))

                    ## For ihec: exchange peak score by log10 q-value and generate bigBed
                    jobs.append(
                        concat_jobs([
                            Job([os.path.join(output_dir, peak_file)],
                                [os.path.join(output_dir, peak_bed_file)],
                                command=f"""\
awk '{{if ($9 > 1000) {{$9 = 1000}}; printf( \"%s\\t%s\\t%s\\t%s\\t%0.f\\n\", $1,$2,$3,$4,$9)}}' {os.path.join(output_dir, peak_file)} > {os.path.join(output_dir, peak_bed_file)}"""
                            ),
                            ucsc.bedToBigBed(
                                os.path.join(output_dir, peak_bed_file),
                                os.path.join(output_dir, f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak.bb")
                            )
                        ],
                            name=f"macs2_callpeak_bigBed.{sample.name}.{mark_name}"
                        )
                    )
                # Else if mark type is Input
                else:
                    log.warning(f"Mark {mark_name} for Sample {sample.name} is an Input... skipping")

        return jobs

    def macs2_atacseq_callpeak(self):
        """
        Peaks are called using the MACS2 software. Different calling strategies are used for narrow and broad peaks.
        The mfold parameter used in the model building step is estimated from a peak enrichment diagnosis run.
        The estimated mfold lower bound is 10 and the estimated upper bound can vary between 15 and 100.
        The default mfold parameter of MACS2 is [10,30].
        Returns:
            list: List of macs2_callpeak jobs for atacseq.
        """

        jobs = []

        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        for sample in self.samples:
            mark_list = []
            # if no Input file
            input_file = []
            if global_conf.global_get('bedtools_intersect', 'blacklist', required=False, param_type='filepath'):
                input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam") for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            else:
                input_file_list = [os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam") for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            if len(input_file_list) > 0:
                if len(input_file_list) > 1:
                    raise Exception(f"""Error: Sample "{sample.name}" has more than 1 Input!""")
                input_file = [input_file_list[0]]
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    filtered_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")]
                    cleaned_file = [
                        os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")]
                    candidate_input_files = [cleaned_file, filtered_file]
                    mark_file = self.select_input_files(candidate_input_files)
                    output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    ## set macs2 variables:
                    options = "--format " + ("BAMPE" if self.run_type == "PAIRED_END" else "BAM")
                    genome_size = global_conf.global_get('macs2_callpeak', 'genome_size', required=False) if global_conf.global_get('macs2_callpeak', 'genome_size', required=False) else self.mappable_genome_size()
                    output_prefix_name = os.path.join(output_dir, f"{sample.name}.{mark_name}")
                    peak_file = f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak"
                    peak_bed_file = f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak.bed"

                    output = []
                    output.append(os.path.join(output_dir, peak_file))
                    output.append(os.path.join(output_dir, f"{sample.name}.{mark_name}_peaks.xls"))
                    # other_options = " --broad --nomodel --bdg --SPMR --keep-dup all"
                    options += " --nomodel --call-summits"
                    options += " --shift " + global_conf.global_get('macs2_callpeak', 'shift') if global_conf.global_get('macs2_callpeak', 'shift') else " --shift -75 "
                    options += " --extsize " + global_conf.global_get('macs2_callpeak', 'extsize') if global_conf.global_get('macs2_callpeak', 'extsize') else " --extsize 150 "
                    options += " -p " + global_conf.global_get('macs2_callpeak', 'pvalue') if global_conf.global_get('macs2_callpeak', 'pvalue') else " -p 0.01 "

                    jobs.append(
                        concat_jobs([
                            bash.mkdir(output_dir),
                            bash.mkdir(link_directory),
                            macs2.callpeak(
                                options,
                                genome_size,
                                mark_file,
                                input_file,
                                output_prefix_name,
                                output
                            ),
                            bash.ln(
                                os.path.relpath(output[1], link_directory),
                                os.path.join(link_directory, f"{sample.name}.{mark_name}_peaks.xls"),
                                input = output[1]
                            )
                        ],
                            name=f"macs2_callpeak.{sample.name}.{mark_name}",
                            removable_files=[output_dir]
                        )
                    )
                    self.multiqc_inputs.append(os.path.join(link_directory, f"{sample.name}.{mark_name}_peaks.xls"))

                    ## For ihec: exchange peak score by log10 q-value and generate bigBed
                    jobs.append(
                        concat_jobs([
                            Job([os.path.join(output_dir, peak_file)],
                                [os.path.join(output_dir, peak_bed_file)],
                                command=f"""\
awk '{{if ($9 > 1000) {{$9 = 1000}}; printf( \"%s\\t%s\\t%s\\t%s\\t%0.f\\n\", $1,$2,$3,$4,$9)}}' {os.path.join(output_dir, peak_file)} > {os.path.join(output_dir, peak_bed_file)}"""
                                ),
                            ucsc.bedToBigBed(
                                os.path.join(output_dir, peak_bed_file),
                                os.path.join(output_dir, f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak.bb")
                            )
                        ],
                            name=f"macs2_callpeak_bigBed.{sample.name}.{mark_name}"
                        )
                    )
                # Else if mark type is Input
                else:
                    log.warning(f"Mark {mark_name} for Sample {sample.name} is an Input... skipping")

        return jobs

    def differential_binding(self):
        """
        Performs differential binding analysis using [DiffBind](http://bioconductor.org/packages/release/bioc/html/DiffBind.html)
        Merge the results of the analysis in a single csv file.
        html report will be generated to QC samples and check how well differential binding analysis was performed.
        Returns:
            list: List of differential_binding jobs.
        """
        jobs = []
        min_overlap = global_conf.global_get('differential_binding', 'minOverlap')
        min_members = global_conf.global_get('differential_binding', 'minMembers')
        method = global_conf.global_get('differential_binding', 'method')
        # If --design <design_file> option is missing, self.contrasts call will raise an Exception
        readset_file = os.path.relpath(self.readsets_file.name, self.output_dir)
        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)

            # If control samples and treatment samples are less than one diff analysis will not be executed
            for contrast in self.contrasts:
                bam_list = []
                controls_count = len(contrast.controls)
                treatments_count = len(contrast.treatments)
                if controls_count < 2 or treatments_count < 2:
                    log.info(f"At leaset two treatments and  controls should be defined. Skipping differential binding analysis for {contrast.name}...")
                else:
                    if global_conf.global_get('bedtools_intersect', 'blacklist', required=False, param_type='filepath'):
                        bam_ext = "sorted.dup.filtered.cleaned.bam"
                    else:
                        bam_ext = "sorted.dup.filtered.bam"

                    for control in contrast.controls:
                        control_sample_name, control_mark_name = control.split("-.-")

                        for sample in self.samples:
                            input_file_list = [
                                os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.{bam_ext}")
                                for mark_name, mark_type in sample.marks.items() if mark_type == "I" and sample.name == control_sample_name
                            ]
                            bam_list.append(input_file_list)

                            input_file_list = [
                                os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.{bam_ext}")
                                for mark_name, mark_type in sample.marks.items() if mark_type != "I" and sample.name == control_sample_name and mark_name == control_mark_name
                            ]
                            bam_list.append(input_file_list)

                            input_file_list = [
                                os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}_peaks.xls")
                                for mark_name, mark_type in sample.marks.items() if mark_type != "I" and sample.name == control_sample_name and mark_name == control_mark_name
                            ]
                            bam_list.append(input_file_list)

                    for treatment in contrast.treatments:
                        treatment_sample_name, treatment_mark_name = treatment.split("-.-")
                        for sample in self.samples:
                            input_file_list = [
                                os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.{bam_ext}")
                                for mark_name, mark_type in sample.marks.items() if mark_type == "I" and sample.name == treatment_sample_name
                            ]
                            bam_list.append(input_file_list)

                            input_file_list = [
                                os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.{bam_ext}")
                                for mark_name, mark_type in sample.marks.items() if mark_type != "I" and sample.name == treatment_sample_name and mark_name == treatment_mark_name
                            ]
                            bam_list.append(input_file_list)

                            input_file_list = [
                                os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}_peaks.xls")
                                for mark_name, mark_type in sample.marks.items() if mark_type != "I" and sample.name == treatment_sample_name and mark_name == treatment_mark_name
                            ]
                            bam_list.append(input_file_list)

                    bam_list = filter(None, bam_list)
                    bam_list = [item for sublist in bam_list for item in sublist]

                    diffbind_job = differential_binding.diffbind(
                        bam_list,
                        contrast.name,
                        design_file,
                        readset_file,
                        self.output_dirs['dba_output_directory'],
                        self.output_dirs['alignment_output_directory'],
                        bam_ext,
                        self.output_dirs['macs_output_directory'],
                        min_overlap,
                        min_members,
                        method
                    )
                    diffbind_job.samples = self.samples
                    diffbind_job.name = "_".join(("differential_binding.diffbind.contrast", contrast.name))
                    jobs.append(diffbind_job)
        else:
            log.info("Comparison column is not defined. Skipping differential binding analysis...")

        return jobs

    def homer_annotate_peaks(self):
        """
        The peaks called previously are annotated with HOMER(http://homer.ucsd.edu/homer/index.html) using RefSeq annotations for the reference genome.
        Gene ontology and genome ontology analysis are also performed at this stage.
        Returns:
            list: List of homer_annotate_peaks jobs.
        """

        jobs = []

        for sample in self.samples:
            mark_list = []
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    peak_file = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak")
                    output_dir = os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)
                    output_prefix = os.path.join(output_dir, f"{sample.name}.{mark_name}")
                    annotation_file = output_prefix + ".annotated.csv"

                    genome = global_conf.global_get('homer_annotate_peaks', 'genome', required=False) if global_conf.global_get('homer_annotate_peaks', 'genome', required=False) else self.ucsc_genome
                    genome_size = global_conf.global_get('homer_annotate_peaks', 'genome_size', required=False) if global_conf.global_get('homer_annotate_peaks', 'genome_size', required=False) else self.mappable_genome_size()

                    annotate_peaks_job = concat_jobs(
                        [
                            bash.mkdir(output_dir),
                            homer.annotatePeaks(
                                peak_file,
                                genome,
                                output_dir,
                                annotation_file,
                                genome_size
                                )
                            ]
                        )

                    jobs.append(
                        concat_jobs(
                            [
                                annotate_peaks_job,
                                Job(
                                    [annotation_file],
                                    [
                                        output_prefix + ".tss.stats.csv",
                                        output_prefix + ".exon.stats.csv",
                                        output_prefix + ".intron.stats.csv",
                                        output_prefix + ".tss.distance.csv"
                                    ],
                                    [
                                        ['homer_annotate_peaks', 'module_perl'],
                                        ['homer_annotate_peaks', 'module_mugqic_tools']
                                    ],
                                    command=f"""\
perl -MReadMetrics -e 'ReadMetrics::parseHomerAnnotations(
  "{annotation_file}",
  "{output_prefix}",
  {global_conf.global_get('homer_annotate_peaks', 'proximal_distance', param_type='int')},
  {global_conf.global_get('homer_annotate_peaks', 'distal_distance', param_type='int')},
  {global_conf.global_get('homer_annotate_peaks', 'distance5d_lower', param_type='int')},
  {global_conf.global_get('homer_annotate_peaks', 'distance5d_upper', param_type='int')},
  {global_conf.global_get('homer_annotate_peaks', 'gene_desert_size', param_type='int')}
)'""",
                                    removable_files=[os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)],
                                )
                            ],
                            name=f"homer_annotate_peaks.{sample.name}.{mark_name}",
                            samples=[sample]
                        )
                    )

                else:
                    log.warning(f"Mark {mark_name} for Sample {sample.name} is an Input... skipping")

        return jobs

    def homer_find_motifs_genome(self):
        """
        De novo and known motif analysis per design are performed using HOMER.
        Returns:
            list: List of homer_find_motifs_genome jobs.
        """

        jobs = []

        for sample in self.samples:
            for mark_name, mark_type in sample.marks.items():
                # Don't find motifs for broad peaks
                if mark_type == "N":

                    peak_file = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak")
                    output_dir = os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)

                    genome = global_conf.global_get('homer_annotate_peaks', 'genome', required=False) if global_conf.global_get('homer_annotate_peaks', 'genome', required=False) else self.ucsc_genome

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(output_dir),
                                homer.findMotifsGenome(
                                    peak_file,
                                    genome,
                                    output_dir
                                )
                            ],
                            name=f"homer_find_motifs_genome.{sample.name}.{mark_name}",
                            removable_files=[os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name)]
                        )
                    )
                else:
                    log.warning(f"Mark {mark_name} for Sample {sample.name}is not Narrow; homer_find_motifs_genome is run on narrow peaks... skipping")

        return jobs

    def annotation_graphs(self):
        """
        The peak location statistics. The following peak location statistics are generated per design:
        proportions of the genomic locations of the peaks. The locations are: Gene (exon or intron), Proximal ([0;2] kb upstream of a transcription start site), Distal ([2;10] kb upstream of a transcription start site), 5d ([10;100] kb upstream of a transcription start site), Gene desert (>= 100 kb upstream or downstream of a transcription start site), Other (anything not included in the above categories); The distribution of peaks found within exons and introns; The distribution of peak distance relative to the transcription start sites (TSS); the Location of peaks per design.
        Returns:
            list: List of annotation_graphs jobs.
        """

        readset_file = os.path.relpath(self.readsets_file.name, self.output_dir)

        input_files = []
        output_files = []
        samples_associative_array = []
        for sample in self.samples:
            mark_list = []
            for mark_name, mark_type in sample.marks.items():
                if mark_type == "N":
                    annotation_prefix = os.path.join(self.output_dirs['anno_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}")
                    input_files.append(annotation_prefix + ".tss.stats.csv")
                    input_files.append(annotation_prefix + ".exon.stats.csv")
                    input_files.append(annotation_prefix + ".intron.stats.csv")
                    input_files.append(annotation_prefix + ".tss.distance.csv")
                    mark_list.append(mark_name)
            samples_associative_array.append("[\"" + sample.name + "\"]=\"" + " ".join(mark_list) + "\"")

            peak_stats_file = os.path.join(self.output_dirs['anno_output_directory'], sample.name, "peak_stats.csv")
            output_files.append(peak_stats_file)

        jobs = []
        report_output_directory = self.output_dirs['report_output_directory']
        graphs_output_directory = self.output_dirs['graphs_output_directory']

        jobs.append(
            Job(
                input_files,
                output_files,
                [
                    ['annotation_graphs', 'module_mugqic_tools'],
                    ['annotation_graphs', 'module_R'],
                    ['annotation_graphs', 'module_pandoc']
                ],
                command=f"""\
cp /dev/null annotation/peak_stats_AllSamples.csv && \\
mkdir -p {graphs_output_directory} && \\
Rscript $R_TOOLS/chipSeqgenerateAnnotationGraphs.R \\
  {readset_file} \\
  {self.output_dir} && \\
declare -A samples_associative_array=({" ".join(samples_associative_array)}) && \\
first_time=true && \\
for sample in ${{!samples_associative_array[@]}}
do
    if $first_time; then
        head -n 1 annotation/$sample/peak_stats.csv > annotation/peak_stats_AllSamples.csv
        first_time=false
    fi
    tail -n+2 annotation/$sample/peak_stats.csv >> annotation/peak_stats_AllSamples.csv
done && \\
mkdir -p {report_output_directory}/annotation/$sample && \\
cp annotation/peak_stats_AllSamples.csv {report_output_directory}/annotation/peak_stats_AllSamples.csv && \\
peak_stats_table=`LC_NUMERIC=en_CA awk -F "," '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:|-----:|-----:|-----:"}} else {{print $1, $2,  sprintf("%\\47d", $3), $4, sprintf("%\\47.1f", $5), sprintf("%\\47.1f", $6), sprintf("%\\47.1f", $7), sprintf("%\\47.1f", $8)}}}}' annotation/peak_stats_AllSamples.csv` && \\
for sample in ${{!samples_associative_array[@]}}
do
  cp annotation/$sample/peak_stats.csv {report_output_directory}/annotation/$sample/peak_stats.csv && \\
  for mark_name in ${{samples_associative_array[$sample]}}
  do
    cp --parents {graphs_output_directory}/${{sample}}.${{mark_name}}_Misc_Graphs.ps {report_output_directory}/
    convert -rotate 90 {graphs_output_directory}/${{sample}}.${{mark_name}}_Misc_Graphs.ps {report_output_directory}/graphs/${{sample}}.${{mark_name}}_Misc_Graphs.png
  done
done""",
                name="annotation_graphs",
                removable_files=output_files
            )
        )

        return jobs


    def run_spp(self):
        """
        runs spp to estimate NSC and RSC ENCODE metrics. For more information: https://github.com/kundajelab/phantompeakqualtools
        Returns:
            list: List of run_spp jobs.
        """
        jobs = []
        ihec_metrics_dir = self.output_dirs['ihecM_output_directory']

        for sample in self.samples:
            for mark_name in sample.marks:
                alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name)
                filtered_file = [
                    os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")]
                cleaned_file = [
                    os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")]
                candidate_input_files = [cleaned_file, filtered_file]
                [sample_merge_mdup_bam] = self.select_input_files(candidate_input_files)
                output_dir = os.path.join(ihec_metrics_dir, sample.name, mark_name)
                output = os.path.join(output_dir, f"{sample.name}.{mark_name}.crosscor")

                jobs.append(
                    concat_jobs([
                        bash.mkdir(output_dir),
                        Job(
                            [sample_merge_mdup_bam],
                            [output],
                            [
                                ['run_spp', 'module_samtools'],
                                ['run_spp', 'module_mugqic_tools'],
                                ['run_spp', 'module_R']
                            ],
                            command=f"""\
cat /dev/null > {output} && \\
Rscript $R_TOOLS/run_spp.R -c={sample_merge_mdup_bam} -savp -out={output} -rf -tmpdir={global_conf.global_get('run_spp', 'tmp_dir')}"""
                        )
                    ],
                        name=f"run_spp.{sample.name}.{mark_name}")
                )
        samples_associative_array = " ".join(["[\"" + sample.name + "\"]=\"" + " ".join(sample.marks.keys()) + "\"" for sample in self.samples])
        jobs.append(
            Job(
                [os.path.join(ihec_metrics_dir, sample.name, mark_name, f"{sample.name}.{mark_name}.crosscor") for sample in self.samples for mark_name, mark_type in sample.marks.items()],
                [os.path.join(ihec_metrics_dir, sample.name, sample.name + ".crosscor") for sample in self.samples],
                [],
                command=f"""\
declare -A samples_associative_array=({samples_associative_array}) && \\
for sample in ${{!samples_associative_array[@]}}
do
echo -e "Filename\\tnumReads\\testFragLen\\tcorr_estFragLen\\tPhantomPeak\\tcorr_phantomPeak\\targmin_corr\\tmin_corr\\tNormalized SCC (NSC)\\tRelative SCC (RSC)\\tQualityTag)" > {ihec_metrics_dir}/${{sample}}/${{sample}}.crosscor
for mark_name in ${{samples_associative_array[$sample]}}
do
cat {ihec_metrics_dir}/${{sample}}/${{mark_name}}/${{sample}}.${{mark_name}}.crosscor >> {ihec_metrics_dir}/${{sample}}/${{sample}}.crosscor
done
done""",
                name="run_spp_report",
                samples=self.samples
            )
        )

        return jobs

    def ihec_metrics(self):
        """
        Generate IHEC's standard metrics.
        """
        jobs = []

        alignment_dir = self.output_dirs['alignment_output_directory']
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        metrics_to_merge = []

        for sample in self.samples:
            mark_list = []
            # if no Input file
            input_file = {}
            input_file_list = [mark_name for mark_name, mark_type in sample.marks.items() if mark_type == "I"]
            if len(input_file_list) > 0:
                if len(input_file_list) > 1:
                    raise Exception(f"""Error: Sample "{sample.name}" has more than 1 Input!""")
                input_file[sample.name] = input_file_list[0]
            for mark_name, mark_type in sample.marks.items():
                if mark_type != "I":
                    mark_list.append(mark_name)

                    chip_bam = os.path.join(alignment_dir, sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.bam")
                    chip_bed = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}_peaks.{self.mark_type_conversion[mark_type]}Peak.bed")
                    output_dir = os.path.join(self.output_dirs['ihecM_output_directory'], sample.name)
                    crosscor_input = os.path.join(self.output_dirs['ihecM_output_directory'], sample.name, sample.name + ".crosscor")
                    genome = global_conf.global_get('IHEC_chipseq_metrics', 'assembly')

                    if not input_file:
                        input_name = "no_input"
                        input_bam = None
                    else:
                        input_name = input_file[sample.name]  # "".join(input_file.keys())
                        input_bam = os.path.join(alignment_dir, sample.name, input_name, f"{sample.name}.{input_name}.sorted.dup.bam")  # input_file[sample.name]

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(output_dir),
                                tools.sh_ihec_chip_metrics(
                                    chip_bam=chip_bam,
                                    input_bam=input_bam,
                                    sample_name=sample.name,
                                    input_name=input_name,
                                    chip_name=mark_name,
                                    chip_type=self.mark_type_conversion[mark_type],
                                    chip_bed=chip_bed,
                                    output_dir=output_dir,
                                    assembly=genome,
                                    crosscor_input=crosscor_input
                                )
                            ],
                            name=f"IHEC_chipseq_metrics.{sample.name}.{mark_name}",
                            removable_files=[output_dir]
                        )
                    )
                    metrics_to_merge.append(os.path.join(output_dir, mark_name, f"IHEC_chipseq_metrics.{sample.name}.{mark_name}.tsv"))

        metrics_merged = "IHEC_chipseq_metrics_AllSamples.tsv"
        metrics_merged_out = os.path.join(self.output_dirs['ihecM_output_directory'], metrics_merged)
        ihec_multiqc_file = os.path.join(link_directory, "ChipSeq.ihecMetrics_mqc.tsv")

        jobs.append(
            Job(
                input_files=metrics_to_merge,
                output_files=[metrics_merged_out],
                name="merge_ihec_metrics",
                command=f"""\
cp /dev/null {metrics_merged_out} && \\
first_time=true && \\
for sample in {" ".join(metrics_to_merge)}
do
    if $first_time; then
        head -n 1 $sample | cut -f -3,5-17,30-33,35,37,39- > {metrics_merged_out}
        first_time=false
    fi
    tail -n 1 $sample | cut -f -3,5-17,30-33,35,37,39- >> {metrics_merged_out}
    sample_name=`tail -n 1 $sample | cut -f 1`
    input_name=`tail -n 1 $sample | cut -f 4`
    input_chip_type="NA"
    genome_assembly=`tail -n 1 $sample | cut -f 5`
    input_core=`tail -n 1 $sample | cut -f 18-29`
    input_nsc=`tail -n 1 $sample | cut -f 34`
    input_rsc=`tail -n 1 $sample | cut -f 36`
    input_quality=`tail -n 1 $sample | cut -f 38`
    if [[ $input_name != "no_input" ]]
    then
        echo -e "${{sample_name}}\\t${{input_name}}\\t${{input_chip_type}}\\t${{genome_assembly}}\\t${{input_core}}\\tNA\\tNA\\tNA\\t${{input_nsc}}\\t${{input_rsc}}\\t${{input_quality}}\\tNA\\tNA" >> {metrics_merged_out}
    fi
done""",
            )
        )

        return jobs

    def multiqc_report(self):
        """
        A quality control report for all samples is generated.
        For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)
        Returns:
            list: List of multiqc_report
        """
        ## set multiQc config file so we can customize one for every pipeline:
        jobs = []
        input_links = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")

        output = os.path.join(self.output_dirs['report_output_directory'], f"ChipSeq.{self.protocol}.multiqc")

        job = multiqc.run(
            input_links,
            output,
            ini_section='multiqc_report'
        )
        job.name = "multiqc_report"

        jobs.append(job)

        return jobs

    def cram_output(self):
        """
        Generate long term storage version of the final alignment files in CRAM format
        Using this function will include the orginal final bam file into the  removable file list
        Returns:
            list: List of cram_output
        """

        jobs = []

        for sample in self.samples:
            for mark_name in sample.marks:
                filtered_bam = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                clean_bam = os.path.join(self.output_dirs['alignment_output_directory'], sample.name, mark_name, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")
                candidate_input_files = [[clean_bam], [filtered_bam]]
                [input_bam] = self.select_input_files(candidate_input_files)
                output_cram = re.sub(r"\.bam$", ".cram", input_bam)

                # Run samtools
                job = samtools.view(
                    input_bam,
                    output_cram,
                    options=global_conf.global_get('samtools_cram_output', 'options'),
                    removable=False
                )
                job.name = f"cram_output.{sample.name}.{mark_name}"
                job.removable_files = input_bam

                jobs.append(job)

        return jobs

    def gatk_haplotype_caller(self):
        """
        GATK haplotype caller for snps and small indels.
        Returns:
            list: List of gatk_haplotype_caller jobs.
        """

        jobs = []


        for sample in self.samples:
            for mark_name, mark_type in sample.marks.items():
                if not mark_type == "I":
                    alignment_directory = os.path.join(self.output_dirs['alignment_output_directory'], sample.name,
                                               mark_name)
                    haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")

                    macs_output_dir = os.path.join(self.output_dirs['macs_output_directory'], sample.name, mark_name)

                    filtered_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.bam")
                    clean_bam = os.path.join(alignment_directory, f"{sample.name}.{mark_name}.sorted.dup.filtered.cleaned.bam")
                    candidate_input_files = [[clean_bam], [filtered_bam]]
                    [input_bam] = self.select_input_files(candidate_input_files)

                    #peak calling bed file from MACS2 is given here to restrict the variant calling to peaks regions
                    interval_list = None
                    if mark_type == "N":
                        interval_list = os.path.join(macs_output_dir, f"{sample.name}.{mark_name}_peaks.narrowPeak.bed")
                    elif mark_type == "B":
                        interval_list = os.path.join(macs_output_dir, f"{sample.name}.{mark_name}_peaks.broadPeak.bed")


                    mkdir_job = bash.mkdir(
                                haplotype_directory,
                                remove=True
                    )
                    jobs.append(
                    concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    mkdir_job,
                    gatk4.haplotype_caller(
                        input_bam,
                        os.path.join(haplotype_directory, f"{sample.name}.{mark_name}.hc.g.vcf.gz"),
                        interval_list=interval_list

                    )
                ],
                    name=f"gatk_haplotype_caller.{sample.name}.{mark_name}",
                    samples=[sample]
                    )
                    )

        return jobs

    def merge_and_call_individual_gvcf(self):
        """
        Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.
        Returns:
            list: List of merge_and_call_individual_gvcf jobs.
        """

        jobs = []

        for sample in self.samples:
            for mark_name, mark_type in sample.marks.items():
                if not mark_type == "I":
                    alignment_directory = os.path.join(self.output_dirs["alignment_output_directory"], sample.name, mark_name)
                    haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
                    haplotype_file_prefix = os.path.join(haplotype_directory , f"{sample.name}.{mark_name}")
                    output_haplotype_file_prefix = os.path.join(self.output_dirs["alignment_output_directory"], sample.name, mark_name, f"{sample.name}.{mark_name}")

                    jobs.append(
                        concat_jobs(
                            [
                                bash.ln(
                                    os.path.relpath(f"{haplotype_file_prefix}.hc.g.vcf.gz", os.path.dirname(f"{output_haplotype_file_prefix}.hc.g.vcf.gz")),
                                    f"{output_haplotype_file_prefix}.hc.g.vcf.gz",
                                    input=f"{haplotype_file_prefix}.hc.g.vcf.gz"
                                ),
                                bash.ln(
                                    os.path.relpath(f"{haplotype_file_prefix}.hc.g.vcf.gz.tbi", os.path.dirname(f"{output_haplotype_file_prefix}.hc.g.vcf.gz.tbi")),
                                    f"{output_haplotype_file_prefix}.hc.g.vcf.gz.tbi",
                                    input=f"{haplotype_file_prefix}.hc.g.vcf.gz.tbi"
                                ),
                                gatk4.genotype_gvcf(
                                    f"{output_haplotype_file_prefix}.hc.g.vcf.gz",
                                    f"{output_haplotype_file_prefix}.hc.vcf.gz",
                                    options=global_conf.global_get('merge_and_call_individual_gvcf', 'options'),
                                    ini_section='merge_and_call_individual_gvcf'
                                )
                            ],
                            name=f"merge_and_call_individual_gvcf.call.{sample.name}.{mark_name}",
                            samples=[sample]
                        )
                    )
        return jobs

    @property
    def step_list(self):
        """
        List of steps for the pipeline
        Returns:
            list: List of steps
        """
        return self.protocols()[self._protocol]

    def protocols(self):
        """
        Returns the protocol for the pipeline.
        Returns:
            dict: A dictionary of protocols for the pipeline.
        """
        return {'chipseq':
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.mapping_bwa_mem_sambamba,
                self.sambamba_merge_bam_files, #5
                self.sambamba_mark_duplicates,
                self.sambamba_view_filter,
                self.bedtools_blacklist_filter,
                self.metrics,
                self.homer_make_tag_directory,
                self.qc_metrics,
                self.homer_make_ucsc_file,  #12
                self.macs2_callpeak,
                self.homer_annotate_peaks,
                self.homer_find_motifs_genome,
                self.annotation_graphs,
                self.run_spp,
                self.differential_binding, #18
                self.ihec_metrics,
                self.multiqc_report,
                self.cram_output,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf #23
            ], 'atacseq':
            [
                self.picard_sam_to_fastq,
                self.trimmomatic,
                self.merge_trimmomatic_stats,
                self.mapping_bwa_mem_sambamba,
                self.sambamba_merge_bam_files,
                self.sambamba_mark_duplicates,
                self.sambamba_view_filter,
                self.bedtools_blacklist_filter,
                self.metrics,
                self.homer_make_tag_directory,
                self.qc_metrics,
                self.homer_make_ucsc_file,
                self.macs2_atacseq_callpeak,
                self.homer_annotate_peaks,
                self.homer_find_motifs_genome,
                self.annotation_graphs,
                self.run_spp,
                self.differential_binding,
                self.ihec_metrics,
                self.multiqc_report,
                self.cram_output,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf
            ]
        }

def main(parsed_args):
    """
    The function that will call this pipeline!
    Parameters:
        parsed_args (argparse.Namespace): The parsed arguments from the command line.
    """

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    no_json = parsed_args.no_json
    json_pt = parsed_args.json_pt
    force = parsed_args.force
    force_mem_per_cpu = parsed_args.force_mem_per_cpu
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file

    # Specific pipeline options
    protocol = parsed_args.protocol

    pipeline = ChipSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container, protocol=protocol)

    pipeline.submit_jobs()
