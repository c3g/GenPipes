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
# along with GenPipes. If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import argparse
import logging
import math
import os
import re
import gzip

# GenPipes Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs, pipe_jobs
from ...core.sample_tumor_pairs import parse_tumor_pair_file
from .. import common
from ...bfx.sequence_dictionary import parse_sequence_dictionary_file, split_by_size

from ...bfx import (
    amber,
    bash_cmd as bash,
    bcbio_variation_recall,
    bcftools,
    breakseq2,
    bvatools,
    bwa,
    cnvkit,
    cobalt,
    conpair,
    cpsr,
    deliverables,
    delly,
    gatk,
    gatk4,
    gemini,
    gridss,
    gripss,
    htslib,
    job2json_project_tracking,
    linx,
    lumpy,
    manta,
    metasv,
    metrics,
    metric_tools,
    multiqc,
    pave,
    pcgr,
    purple,
    mosdepth,
    sambamba,
    samtools,
    sequence_dictionary,
    sequenza,
    snpeff,
    strelka2,
    svtyper,
    tools,
    vardict,
    varscan,
    vawk,
    vcftools,
    verify_bam_id,
    vt,
    wham
    )

log = logging.getLogger(__name__)

class DnaSeqRaw(common.Illumina):
    """
DNA-Seq Pipeline
================

A pipeline to process DNA sequencing data. The pipeline uses Trimmomatic for quality control, BWA for alignment to a reference genome, Picard for marking duplicates, GATK for indel realignment and variant calling, SnpEff for variant annotation, SnpSift for filtering variants and MultiQC for aggregate reports.

The pipeline contains protocols for processing both germline and somatic sequencing data; high-coverage data from targeted sequencing experiments; SNVs and structural variants.

The pipeline is designed to be run on a cluster and is configured using a configuration file. The pipeline can be run in a single step or in multiple steps. The pipeline can also be run in parallel to process multiple samples simultaneously.

Attributes:
    pairs_file (str): A string containing the path to the pairs file.
    profyle (str): A string containing the path to the profyle file.
    pairs (str): The pairs attribute is a property that returns the parsed tumor pairs file.
    output_dirs (dict): The output_dirs attribute is a property that returns a dictionary of output directories.
    multiqc_inputs (dict): The multiqc_inputs attribute is a property that returns a dictionary of multiqc inputs.
    sequence_dictionary (dict): The sequence_dictionary attribute is a property that returns a parsed sequence dictionary file.
    sequence_dictionary_variant (dict): The sequence_dictionary_variant attribute is a property that returns a parsed sequence dictionary file.
    tumor_pairs (dict): The tumor_pairs attribute is a property that returns a parsed tumor pairs file.
Methods:
    generate_approximate_windows(): The generate_approximate_windows method generates approximate windows based on the number of jobs.
    is_gz_file(): The is_gz_file method checks if a file is a gzipped file.
    sym_link_fastq(): The sym_link_fastq method creates a symbolic link of raw reads fastq files.
    sym_link_fastq_pair(): The sym_link_fastq_pair method creates symbolic links and md5 sums for tumor and normal fastq files.
    build_adapter_file(): The build_adapter_file method builds an adapter file.
    bwa_mem_sambamba_sort_sam(): The bwa_mem_sambamba_sort_sam method aligns filtered reads to a reference genome.
    gatk_fixmate(): The gatk_fixmate method verifies mate-pair information between mates and fixes if needed.
    mark_duplicates(): The mark_duplicates method marks duplicates in aligned reads.
    sym_link_final_bam(): The sym_link_final_bam method creates a symbolic link of the final bam for delivery of data to clients.
    sym_link_final_bam_pair(): The sym_link_final_bam_pair method creates a symbolic link of the final bam for delivery of data to clients.
    conpair_concordance_contamination(): The conpair_concordance_contamination method performs concordance verification and contamination level estimation.
    metrics_dna_picard_metrics(): The metrics_dna_picard_metrics method generates metrics with picard.
Parameters:
    pairs_file (str): A string containing the path to the pairs file.
    profyle (str): A string containing the path to the profyle file.
    """

    def __init__(self, *args, pairs_file=None, profyle=None, **kwargs):
        self.pairs = pairs_file
        self.profyle = profyle
        # Add pipeline specific arguments
        super(DnaSeqRaw, self).__init__(*args, **kwargs)

    @property
    def tumor_pairs(self):
        """
        The tumor_pairs attribute is a property that returns a parsed tumor pairs file.
        Returns:
            dict: A dictionary of tumor pairs.
        """
        # Create property only if somatic protocols called
        if 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
            if not hasattr(self, "_tumor_pairs"):
                self._tumor_pairs = parse_tumor_pair_file(
                    self.pairs.name,
                    self.samples,
                    self.profyle
                )
            return self._tumor_pairs

    @property
    def output_dirs(self):
        """
        Output directory paths.
        Returns:
            dict: Output directory paths.
        """
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir,"raw_reads"), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, "trim"), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, "alignment"), self.output_dir),
            'metrics_directory': {}, #os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'variants_directory': os.path.relpath(os.path.join(self.output_dir, "variants"), self.output_dir),
            'paired_variants_directory': os.path.relpath(os.path.join(self.output_dir,"pairedVariants"), self.output_dir),
            'sv_variants_directory': os.path.relpath(os.path.join(self.output_dir, "SVariants"), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, "report"), self.output_dir)
        }
        # Somatic protocol properties only
        if 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
            for tumor_pair in self.tumor_pairs.values():
                dirs['metrics_directory'][tumor_pair.name] = os.path.relpath(os.path.join(self.output_dir, "metrics", tumor_pair.name), self.output_dir)

        # Sample-level properties used for both germline and somatic protocols
        for sample in self.samples:
            dirs['metrics_directory'][sample.name] = os.path.relpath(
                os.path.join(f"{self.output_dir}/metrics/{sample.name}"),
                self.output_dir
            )
        return dirs

    @property
    def multiqc_inputs(self):
        """
        List of MultiQC input files.
        Returns:
            list: List of MultiQC input files.
        """
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = {}
            for sample in self.samples:
                self._multiqc_inputs[sample.name] = []

            if 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
                for tumor_pair in self.tumor_pairs.values():
                    self._multiqc_inputs[tumor_pair.name] = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    @property
    def sequence_dictionary(self):
        """
        The sequence_dictionary attribute is a property that returns a parsed sequence dictionary file.
        Returns:
            dict: A parsed sequence dictionary file.
        """
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(global_conf.global_get('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=True)
        return self._sequence_dictionary

    def sequence_dictionary_variant(self):
        """
        The sequence_dictionary_variant attribute is a property that returns a parsed sequence dictionary file.
        Returns:
            dict: A parsed sequence dictionary file.
        """
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = parse_sequence_dictionary_file(global_conf.global_get('DEFAULT', 'genome_dictionary', param_type='filepath'), variant=True)
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        """
        Generate approximate windows based on the number of jobs.
        Arguments:
            nb_jobs (int): The number of jobs.
        Returns:
            list: A list of approximate windows.
        """
        if nb_jobs <= len(self.sequence_dictionary):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary]
        else:
            total_length = sum(sequence['length'] for sequence in self.sequence_dictionary)
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary))))
            windows = []

            for sequence in self.sequence_dictionary:
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(f"{sequence['name']}:{str(start)}-{str(end)}")
            return windows

    def is_gz_file(self, name):
        """
        Check if a file is a gzipped file.
        Arguments:
            name (str): The name of the file.
        Returns:
            bool: True if the file is a gzipped file, False otherwise.
        """
        # Cf. https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
        gzip_magic_number = b'\x1f\x8b'
        try:
            with open(name, 'rb') as file:
                file_start = file.read(2)
                return file_start == gzip_magic_number
        except (OSError, IOError):
            return False

    def sym_link_fastq(self):
        """
        Create sym link of raw reads fastq files.
        Returns:
            list: A list of sym link fastq jobs.
        """
        jobs = []
        sym_link_job = ""
        for readset in self.readsets:
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        readset.name
                    )
                    candidate_input_files.append([f"{prefix}.pair1.fastq.gz", f"{prefix}.pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

                sym_link_job = concat_jobs(
                    [
                        deliverables.sym_link(
                            fastq1,
                            readset,
                            self.output_dir,
                            type="raw_reads"
                        ),
                        deliverables.sym_link(
                            fastq2,
                            readset,
                            self.output_dir,
                            type="raw_reads"
                        )
                    ]
                )
                sym_link_job.name = f"sym_link_fastq.paired_end.{readset.name}"

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        readset.name
                        )
                    candidate_input_files.append([f"{prefix}.single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)

                sym_link_job = deliverables.sym_link(
                    fastq1,
                    readset,
                    self.output_dir,
                    type="raw_reads"
                )
                sym_link_job.name = f"sym_link_fastq.single_end.{readset.name}"

            else:
                _raise(SanitycheckError(f"Error: run type {readset.run_type} is invalid for readset {readset.name} (should be PAIRED_END or SINGLE_END)!"))

            sym_link_job.samples = [readset.sample]
            sym_link_job.readsets = [readset]
            jobs.append(sym_link_job)

        return jobs

    def sym_link_fastq_pair(self):
        """
        Create sym links and md5 sums for tumor and normal fastq files.
        Returns:
            list: A list of sym link fastq pair jobs.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            [inputs["Normal"]] = [
                self.select_input_files(
                    [
                        [
                            readset.fastq1,
                            readset.fastq2
                        ],
                        [
                            os.path.join(self.output_dirs['raw_reads_directory'], readset.sample.name, f"{readset.name}.pair1.fastq.gz"),
                            os.path.join(self.output_dirs['raw_reads_directory'],readset.sample.name, f"{readset.name}.pair2.fastq.gz")
                        ]
                    ]
                ) for readset in tumor_pair.readsets[tumor_pair.normal.name]
            ]

            [inputs["Tumor"]] = [
                self.select_input_files(
                    [
                        [
                            readset.fastq1,
                            readset.fastq2
                        ],
                        [
                            os.path.join(self.output_dirs['raw_reads_directory'], readset.sample.name, f"{readset.name}.pair1.fastq.gz"),
                            os.path.join(self.output_dirs['raw_reads_directory'], readset.sample.name, f"{readset.name}.pair2.fastq.gz")
                        ]
                    ]
                ) for readset in tumor_pair.readsets[tumor_pair.tumor.name]
            ]

            for key, input_files in inputs.items():
                for read, input_file in enumerate(input_files):
                    symlink_pair_job = deliverables.sym_link_pair(
                        input_file,
                        tumor_pair,
                        self.output_dir,
                        type="raw_reads",
                        sample=key,
                        profyle=self.profyle
                    )
                    dir_name, file_name = os.path.split(symlink_pair_job.output_files[0])
                    # do not compute md5sum in the readset input directory
                    md5sum_job = deliverables.md5sum(
                        symlink_pair_job.output_files[0],
                        f"{file_name}.md5",
                        dir_name
                    )
                    jobs.append(
                        concat_jobs(
                            [
                                symlink_pair_job,
                                md5sum_job
                            ],
                            name=f"sym_link_fastq.pairs.{str(read)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.tumor, tumor_pair.normal],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def build_adapter_file(self, directory, readset):
        """
        Build an adapter file.
        Arguments:
            directory (str): The directory to store the adapter file.
            readset (Readset): The readset object.
        Returns:
            Job: An adapter job.
        """
        adapter_job = None
        adapter_file = os.path.join(directory, "adapter.tsv")
        if readset.run_type == "SINGLE_END":
            if readset.adapter1:
                adapter_job = Job(
                    command=f"""\
`cat > {adapter_file} << END
>Adapter\n{readset.adapter1}\n
END
`"""
                )
            else:
                raise Exception(f"Error: missing adapter1 for SINGLE_END readset {readset.name} or missing adapter_file parameter in config file!")

        elif readset.run_type == "PAIRED_END":
            if readset.adapter1 and readset.adapter2:
                adapter_job = Job(
                    command=f"""\
`cat > {adapter_file} << END
>Adapter1\n{readset.adapter1}\n
>Adapter2\n{readset.adapter2}\n
END
`"""
                )

        return adapter_job

    def bwa_mem_sambamba_sort_sam(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Sambamba](http://lomereiter.github.io/sambamba/index.html)
        This step takes as input files:
        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files

        Returns:
            list: A list of bwa mem sambamba sort sam jobs.
        """

        jobs = []
        sequencing_center = global_conf.global_get('bwa_mem_sambamba_sort_sam', 'sequencing_center', required=False)
        sequencing_technology = global_conf.global_get('bwa_mem_sambamba_sort_sam', 'sequencing_technology') if global_conf.global_get('bwa_mem_sambamba_sort_sam', 'sequencing_technology', required=False) else "Illumina"
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, f"{readset.name}.trim")
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")
            index_bam = os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam.bai")

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
                        bash.mkdir(os.path.dirname(readset_bam)),
                        pipe_jobs(
                            [
                                bwa.mem(
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
                                        ini_section="bwa_mem_sambamba_sort_sam"
                                ),
                                sambamba.view(
                                    "/dev/stdin",
                                    None,
                                    "-S -f bam"
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    readset_bam,
                                    tmp_dir=global_conf.global_get('bwa_mem_sambamba_sort_sam', 'tmp_dir', required=True),
                                    other_options=global_conf.global_get('bwa_mem_sambamba_sort_sam', 'sambamba_sort_options', required=True),
                                )
                            ]
                        ),
                        sambamba.index(
                            readset_bam,
                            index_bam
                        )
                    ],
                    name=f"bwa_mem_sambamba_sort_sam.{readset.name}",
                    samples=[readset.sample],
                    readsets=[readset],
                    removable_files=[]
                )
            )

        return jobs

    def gatk_fixmate(self):
        """
        Verify mate-pair information between mates and fix if needed.
        This ensures that all mate-pair information is in sync between each read and its mate pair.
        Fix is done using [Picard](http://broadinstitute.github.io/picard/).
        Returns:
            list: A list of gatk fix mate information jobs.
        """
        jobs = []

        compression_postfix = global_conf.global_get('bwa_mem2_samtools_sort', 'compression')
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], sample.name)
            input_file = os.path.join(alignment_directory, f"{sample.name}.sorted.{compression_postfix}")
            output = os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.{compression_postfix}")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            alignment_directory,
                            remove=False,
                        ),
                        gatk4.fix_mate_information(
                            input_file,
                            output,
                            create_index=True,
                            ini_section='gatk_fix_mate_information'
                        ),
                    ],
                    name=f"gatk_fix_mate_information.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    removable_files=[]
                )
            )

        return jobs

    def mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        Returns:
            list: A list of mark duplicates jobs.
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            picard_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])
            alignment_file_prefix = os.path.join(alignment_directory,sample.name)
            readset = sample.readsets[0]

            [input_file] = self.select_input_files(
                [
                    [f"{alignment_file_prefix}.sorted.matefixed.bam"],
                    [f"{alignment_file_prefix}.sorted.realigned.bam"],
                    [f"{alignment_file_prefix}.sorted.bam"],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")]
                ]
            )
            output = f"{alignment_file_prefix}.sorted.dup.bam"
            output_index = f"{alignment_file_prefix}.sorted.dup.bam.bai"
            metrics_file = f"{alignment_file_prefix}.sorted.dup.metrics"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            picard_directory,
                            remove=False,
                        ),
                        gatk4.mark_duplicates(
                            input_file,
                            output,
                            metrics_file,
                            create_index=False,
                            ini_section='mark_duplicates'
                        ),
                        sambamba.index(
                            output,
                            output_index
                        ),
                        Job(
                            [metrics_file],
                            [os.path.join(f"{picard_directory}/{sample.name}.sorted.dup.metrics")],
                            command=f"sed -e 's#.realigned##g' {metrics_file} > {os.path.join(picard_directory, sample.name + '.sorted.dup.metrics')}"
                        )
                    ],
                    name=f"mark_duplicates.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    removable_files=[]
                )
            )
        return jobs

    def sym_link_final_bam(self):
        """
        Create sym link of final bam for delivery of data to clients.
        Returns:
            list: A list of sym link final bam jobs.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            alignment_file_prefix = os.path.join(alignment_directory, sample.name)
            readset = sample.readsets[0]

            [input_bam] = self.select_input_files(
                [
                    [f"{alignment_file_prefix}.sorted.dup.cram"],
                    [f"{alignment_file_prefix}.sorted.dup.bam"],
                    [f"{alignment_file_prefix}.sorted.bam"],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")]
                ]
            )
            [input_bai] = self.select_input_files(
                [
                    [f"{alignment_file_prefix}.sorted.dup.cram.crai"],
                    [f"{alignment_file_prefix}.sorted.dup.bam.bai"],
                    [f"{alignment_file_prefix}.sorted.dup.bai"],
                    [f"{alignment_file_prefix}.sorted.bam.bai"],
                    [f"{alignment_file_prefix}.sorted.bai"],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.filtered.bam.bai")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.filtered.bai")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam.bai")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bai")]
                ]
            )
            jobs.append(
                concat_jobs(
                    [
                        deliverables.md5sum(
                            input_bam,
                            f"{input_bam}.md5",
                            self.output_dir
                        ),
                        deliverables.sym_link(
                            input_bam,
                            sample,
                            self.output_dir,
                            type="alignment"
                        ),
                        deliverables.sym_link(
                            input_bai,
                            sample,
                            self.output_dir,
                            type="alignment"
                        ),
                        deliverables.sym_link(
                            f"{input_bam}.md5",
                            sample,
                            self.output_dir,
                            type="alignment"
                        )
                    ],
                    name=f"sym_link_final_bam.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    removable_files=[]
                )
            )

        return jobs


    def sym_link_final_bam_pair(self):
        """
        Create sym link of the final bam for delivery of data to clients.
        Returns:
            list: A list of sym link final bam pair jobs.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            [inputs["Normal"]] = [
                self.select_input_files(
                    [
                        [
                            os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam"),
                            os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bai")
                        ],
                        [
                            os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam"),
                            os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam.bai")
                        ]
                    ]
                )
            ]

            [inputs["Tumor"]] = [
                self.select_input_files(
                    [
                        [
                            os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam"),
                            os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bai")
                        ],
                        [
                            os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam"),
                            os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam.bai")
                        ]
                    ]
                )
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    f"{input_file}.md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="alignment",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    f"{input_file}.md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="alignment",
                                    sample=key,
                                    profyle=self.profyle
                                )
                            ],
                            name=f"sym_link_final_bam.pairs.{str(idx)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            removable_files=[]
                        )
                    )
        return jobs

    def conpair_concordance_contamination(self):
        """
        Conpair is a fast and robust method dedicated to human tumor-normal studies to perform concordance verification
        (= samples coming from the same individual), as well as cross-individual contamination level estimation in
        whole-genome and whole-exome sequencing experiments. Importantly, the method of estimating contamination in
        the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.
        Returns:
            list: A list of conpair concordance contamination jobs.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'],
                                                          tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = (
                    os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                )

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            metrics_directory = self.output_dirs['metrics_directory'][tumor_pair.name]

            input_normal = os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")
            input_tumor = os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")
            pileup_normal = os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.gatkPileup")
            pileup_tumor = os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.gatkPileup")

            concordance_out = os.path.join(metrics_directory, f"{tumor_pair.tumor.name}.concordance.tsv")
            contamination_out = os.path.join(metrics_directory,f"{tumor_pair.tumor.name}.contamination.tsv")

            samples = [tumor_pair.normal, tumor_pair.tumor]
            job_name = f"conpair_concordance_contamination.{tumor_pair.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        conpair.parse_concordance_metrics_pt(concordance_out),
                        job2json_project_tracking.run(
                            input_file=concordance_out,
                            pipeline=self,
                            samples=",".join([sample.name for sample in samples]),
                            readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="concordance=$concordance"
                        ),
                        conpair.parse_contamination_normal_metrics_pt(contamination_out),
                        job2json_project_tracking.run(
                            input_file=contamination_out,
                            pipeline=self,
                            samples=tumor_pair.normal.name,
                            readsets=",".join([readset.name for readset in tumor_pair.normal.readsets]),
                            job_name=job_name,
                            metrics="contamination=$contamination"
                        ),
                        conpair.parse_contamination_tumor_metrics_pt(contamination_out),
                        job2json_project_tracking.run(
                            input_file=contamination_out,
                            pipeline=self,
                            samples=tumor_pair.tumor.name,
                            readsets=",".join([readset.name for readset in tumor_pair.tumor.readsets]),
                            job_name=job_name,
                            metrics="contamination=$contamination"
                        )
                    ])

            jobs.append(
                concat_jobs(
                    [
                        conpair.pileup(
                            input_normal,
                            pileup_normal
                        )
                    ],
                    name=f"conpair_concordance_contamination.pileup.{tumor_pair.name}.{tumor_pair.normal.name}",
                    samples=[tumor_pair.normal],
                    readsets=list(tumor_pair.normal.readsets)
                )
            )
            jobs.append(
                concat_jobs(
                    [
                        conpair.pileup(
                            input_tumor,
                            pileup_tumor
                        )
                    ],
                    name=f"conpair_concordance_contamination.pileup.{tumor_pair.name}.{tumor_pair.tumor.name}",
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            metrics_directory,
                            remove=False
                        ),
                        conpair.concordance(
                            pileup_normal,
                            pileup_tumor,
                            concordance_out
                        ),
                        conpair.contamination(
                            pileup_normal,
                            pileup_tumor,
                            contamination_out
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=samples,
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    output_dependency=[concordance_out,contamination_out],
                    removable_files=[pileup_normal, pileup_tumor]
                )
            )
            self.multiqc_inputs[tumor_pair.name].extend(
                [
                    concordance_out,
                    contamination_out
                ]
            )

        return jobs

    def metrics_dna_picard_metrics(self):
        """
        Generates metrics with picard, including:
            [CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectMultipleMetrics-Picard-)
            [CollectOxoGMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037428231-CollectOxoGMetrics-Picard-)
            [CollectGcBiasMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036481572-CollectGcBiasMetrics-Picard-)
            [CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectWgsMetrics-Picard-)
            [CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectHsMetrics-Picard-)
            [CollectInsertSizeMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectInsertSizeMetrics-Picard-)
            [CollectSequencingArtifactMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectSequencingArtifactMetrics-Picard-)
            [CollectQualityYieldMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectQualityYieldMetrics-Picard-)
            [CollectQualityByCycle](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectQualityByCycle-Picard-)
            [CollectBaseDistributionByCycle](https://gatk.broadinstitute.org/hc/en-us/articles/360037594031-CollectBaseDistributionByCycle-Picard-)

        Returns:
            list: A list of picard metrics jobs.
        """

        # Check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for sample in self.samples:
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])

            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            readset = sample.readsets[0]

            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, readset.name, f"{readset.name}.sorted.bam")]
                ]
            )
            mkdir_job = bash.mkdir(
                metrics_directory
            )
            job_name = f"picard_collect_multiple_metrics.{sample.name}"
            output_prefix = os.path.join(metrics_directory, f"{sample.name}.all.metrics")

            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        gatk4.parse_bases_over_q30_percent_metrics_pt(f"{output_prefix}.quality_distribution_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.quality_distribution_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="bases_over_q30_percent=$bases_over_q30_percent"
                        ),
                        gatk4.parse_mean_insert_metrics(f"{output_prefix}.insert_size_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.insert_size_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="mean_insert_size=$mean_insert_size"
                        ),
                        gatk4.parse_stdev_insert_metrics(f"{output_prefix}.insert_size_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.insert_size_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="stdev_insert_size=$stdev_insert_size"
                        ),
                        gatk4.parse_mode_insert_metrics(f"{output_prefix}.insert_size_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.insert_size_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="mode_insert_size=$mode_insert_size"
                        ),
                        gatk4.parse_total_read_pairs_metrics(f"{output_prefix}.alignment_summary_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.alignment_summary_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="total_read_pairs=$total_read_pairs"
                        ),
                        gatk4.parse_aligned_pairs_metrics_pt(f"{output_prefix}.alignment_summary_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.alignment_summary_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="aligned_pairs_percent=$aligned_pairs_percent"
                        ),
                        gatk4.parse_high_quality_read_pairs_metrics(f"{output_prefix}.alignment_summary_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.alignment_summary_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="hq_read_pairs=$hq_read_pairs"
                        ),
                        gatk4.parse_chimeras_metrics_pt(f"{output_prefix}.alignment_summary_metrics"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.alignment_summary_metrics",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="chimeras_percent=$chimeras_percent"
                        ),
                    ]
                )

            jobs.append(
                concat_jobs(
                    [
                        mkdir_job,
                        gatk4.collect_multiple_metrics(
                            input_file,
                            output_prefix,
                            library_type=library[sample]
                        ),
                        job_project_tracking_metrics,
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[
                        f"{output_prefix}.alignment_summary_metrics",
                        f"{output_prefix}.insert_size_metrics",
                        f"{output_prefix}.quality_by_cycle_metrics",
                        f"{output_prefix}.quality_distribution_metrics"
                    ]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        mkdir_job,
                        #bash.mkdir(link_directory),
                        gatk4.collect_oxog_metrics(
                            input_file,
                            f"{output_prefix}.oxog_metrics.txt"
                        ),
                    ],
                    name=f"picard_collect_oxog_metrics.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency = [
                        f"{output_prefix}.oxog_metrics.txt"
                    ],
                    removable_files=[]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        mkdir_job,
                        gatk4.collect_gcbias_metrics(
                            input_file,
                            output_prefix
                        ),
                    ],
                    name=f"picard_collect_gcbias_metrics.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[
                        f"{output_prefix}.gcbias_metrics.txt"
                    ],
                    removable_files=[]
                )
            )
            self.multiqc_inputs[sample.name].extend(
                [
                    f"{output_prefix}.alignment_summary_metrics",
                    f"{output_prefix}.insert_size_metrics",
                    f"{output_prefix}.quality_by_cycle_metrics",
                    f"{output_prefix}.quality_distribution_metrics",
                    f"{output_prefix}.oxog_metrics.txt",
                    f"{output_prefix}.gcbias_metrics.txt"
                ]
            )

        return jobs

    def metrics_dna_sample_mosdepth(self):
        """
        Calculate depth stats for captured regions with [Mosdepth](https://github.com/brentp/mosdepth)
        Returns:
            list: A list of mosdepth jobs.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)

            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.filtered.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")],
                ]
            )
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])
            output_prefix = os.path.join(metrics_directory, sample.name)
            region = None
            output_dist = f"{output_prefix}.mosdepth.global.dist.txt"
            output_summary = f"{output_prefix}.mosdepth.summary.txt"

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                sample.readsets[0]
            )
            if coverage_bed:
                region = coverage_bed
                output_dist = f"{output_prefix}.mosdepth.region.dist.txt"
                output_summary = f"{output_prefix}.mosdepth.summary.txt"

            job_name = f"mosdepth.{sample.name}"

            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        mosdepth.parse_dedup_coverage_metrics_pt(f"{output_prefix}.mosdepth.summary.txt"),
                        job2json_project_tracking.run(
                            input_file=f"{output_prefix}.mosdepth.summary.txt",
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="dedup_coverage=$dedup_coverage"
                        )
                    ]
                )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(metrics_directory),
                        mosdepth.run(
                            input_file,
                            output_prefix,
                            True,
                            region
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[output_dist, output_summary],
                    removable_files=[]
                )
            )
            self.multiqc_inputs[sample.name].extend(
                [
                    os.path.join(metrics_directory, os.path.basename(output_dist)),
                    os.path.join(metrics_directory, os.path.basename(output_summary))
                ]
            )

        return jobs

    def metrics_dna_samtools_flagstat(self):
        """
        Outputs flag statistics from BAM file.
        https://lomereiter.github.io/sambamba/docs/sambamba-flagstat.html
        Returns:
            list: A list of samtools flagstat jobs.
        """

        jobs = []
        for sample in self.samples:
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)

            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.primerTrim.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )
            output = os.path.join(metrics_directory, f"{sample.name}.flagstat")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            metrics_directory
                        ),
                        samtools.flagstat(
                            input_file,
                            output
                        ),
                    ],
                    name=f"dna_samtools_flagstat.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[
                        os.path.join(metrics_directory, f"{sample.name}.flagstat")
                    ]
                )
            )
            self.multiqc_inputs[sample.name].append(
                os.path.join(metrics_directory, f"{sample.name}.flagstat")
                )

        return jobs


    def run_multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report.
        MultiQC searches a given directory for analysis logs and compiles a HTML report. 
        It's a general use tool, perfect for summarising the output from numerous bioinformatics tools.
        https://multiqc.info/
        Returns:
            list: A list of multiqc jobs.
        """

        jobs = []

        # Generate multiqc report for dnaseq and tumor only protocols
        if 'germline' in self.protocol or 'tumor_only' in self.protocol:
            multiqc_files_paths = [item for sample in self.samples for item in self.multiqc_inputs[sample.name]]

            multiqc_input_path = [os.path.join(self.output_dirs['metrics_directory'][sample.name]) for sample in self.samples]

            output = os.path.join(self.output_dirs['report_directory'], f"DnaSeq.{self.protocol}.multiqc")

            job = multiqc.run(
                multiqc_input_path,
                output
            )
            job.name = "multiqc.all_samples"
            job.samples = self.samples
            job.readsets = self.readsets
            job.input_files = multiqc_files_paths
            job.input_dependency = [multiqc_files_paths]
            job.removable_files=[]
            jobs.append(job)

        # Generate multiqc reports for somatics protocols excluding tumor only protocol
        elif 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
            report_directory = os.path.join(self.output_dirs['report_directory'])
            for tumor_pair in self.tumor_pairs.values():
                patient_folders = [
                    os.path.join(self.output_dirs['metrics_directory'][tumor_pair.name]),
                    os.path.join(self.output_dirs['metrics_directory'][tumor_pair.normal.name]),
                    os.path.join(self.output_dirs['metrics_directory'][tumor_pair.tumor.name]),
                ]
                multiqc_input_paths = []
                for metric in (
                        self.multiqc_inputs[tumor_pair.name] +
                        self.multiqc_inputs[tumor_pair.normal.name] +
                        self.multiqc_inputs[tumor_pair.tumor.name]):
                    multiqc_input_paths.append(metric)

                output = os.path.join(report_directory, f"{tumor_pair.name}.{self.protocol}.multiqc")
                job = multiqc.run(
                    patient_folders,
                    output
                )

                job.name = f"multiqc.{tumor_pair.name}"
                job.samples = [tumor_pair.normal, tumor_pair.tumor]
                job.readsets = [*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                job.input_files = multiqc_input_paths
                job.input_dependency = [multiqc_input_paths]
                job.removable_files=[]
                jobs.append(job)

        return jobs

    def sym_link_report(self):
        """
        Create a sym link of the MultiQC report for delivery to clients.
        Returns:
            list: A list of sym link report jobs.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] = [os.path.join(self.output_dirs['report_directory'], f"{tumor_pair.name}.{self.protocol}.multiqc.html")]

            for key, input_files in inputs.items():
                for idx, report_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.sym_link_pair(
                                    report_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="metrics",
                                    sample=key,
                                    profyle=self.profyle
                                )
                            ],
                            name=f"sym_link_fastq.report.{str(idx)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def metrics_picard_calculate_hs(self):
        """
        Compute on target percent of hybridisation based capture with [Picard CollectHsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard)
        Returns:
            list: A list of picard calculate hs metrics jobs.
        """

        jobs = []

        for sample in self.samples:
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                if os.path.isfile(re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)):
                    interval_list = re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)
                else:
                    interval_list = re.sub(r"\.[^.]+$", ".interval_list", os.path.basename(coverage_bed))
                    job = gatk4.bed2interval_list(
                        None,
                        coverage_bed,
                        interval_list
                    )
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)

                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                [input_file] = self.select_input_files(
                    [
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                        [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                    ]
                )
                job_name = f"gatk_calculate_hs_metrics.{sample.name}"
                output = os.path.join(metrics_directory, f"{sample.name}.onTarget.tsv")

                job_project_tracking_metrics = []
                if self.project_tracking_json:
                    job_project_tracking_metrics = concat_jobs(
                        [
                            gatk4.parse_bed_bait_set_metrics(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="bed_bait_set=$bed_bait_set"
                            ),
                            gatk4.parse_off_target_metrics_pt(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="off_target_percent=$off_target_percent"
                            ),
                            gatk4.parse_total_reads_metrics(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="total_reads=$total_reads"
                            ),
                            gatk4.parse_dedup_reads_metrics(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="dedup_reads=$dedup_reads"
                            ),
                            gatk4.parse_mean_target_coverage_metrics(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="mean_target_coverage=$mean_target_coverage"
                            ),
                            gatk4.parse_median_target_coverage_metrics(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="median_target_coverage=$median_target_coverage"
                            ),
                            gatk4.parse_dup_rate_metrics_pt(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="duplicate_rate_percent=$duplicate_rate_percent"
                            ),
                            gatk4.parse_low_mapping_rate_metrics_pt(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="low_mapping_rate_percent=$low_mapping_rate_percent"
                            ),
                            gatk4.parse_read_overlap_metrics_pt(output),
                            job2json_project_tracking.run(
                                input_file=output,
                                pipeline=self,
                                samples=sample.name,
                                readsets=",".join([readset.name for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="read_overlap_percent=$read_overlap_percent"
                            ),

                        ]
                    )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(metrics_directory),
                            gatk4.calculate_hs_metrics(
                                input_file,
                                os.path.join(metrics_directory, f"{sample.name}.onTarget.tsv"),
                                interval_list
                            ),
                            job_project_tracking_metrics
                        ],
                        name=job_name,
                        samples=[sample],
                        readsets=[*list(sample.readsets)],
                        output_dependency=[os.path.join(metrics_directory, f"{sample.name}.onTarget.tsv")],
                        removable_files=[]
                    )
                )
                self.multiqc_inputs[sample.name].append(
                        os.path.join(metrics_directory, f"{sample.name}.onTarget.tsv")
                )
        return jobs

    def metrics_vcftools_missing_indiv(self):
        """
        vcftools: --missing_indv: Generates a file reporting the missingness on a per-individual basis. The file has the suffix ".imiss".
        input: bgzipped vcf file
        ouput: missingness flat file
        Returns:
            list: A list of vcftools missing indiv jobs.
        """

        jobs = []

        input_prefix = os.path.join(self.output_dirs['variants_directory'], "allSamples")

        job = vcftools.missing_indv(f"{input_prefix}.hc.vcf.gz", input_prefix)
        job.name = "vcftools_missing_indv.allSamples"
        job.samples = self.samples
        job.readsets = self.readsets
        job.removable_files = []
        jobs.append(job)

        return jobs

    def metrics_vcftools_depth_indiv(self):
        """
        vcftools: --depth: Generates a file containing the mean depth per individual. This file has the suffix ".idepth".
        input: bgzipped vcf file
        ouput: idepth flat file
        Returns:
            list: A list of vcftools depth indiv jobs.
    	"""

        jobs = []

        input_prefix = os.path.join(self.output_dirs['variants_directory'], "allSamples")

        job = vcftools.depth(f"{input_prefix}.hc.vcf.gz", input_prefix)
        job.name = "vcftools_depth.allSamples"
        job.samples = self.samples
        job.readsets = self.readsets
        job.removable_files = []
        jobs.append(job)

        return jobs

    def metrics_gatk_sample_fingerprint(self):
        """
        CheckFingerprint (Picard)
        Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) 
        against a set of known genotypes in the supplied genotype file (in VCF format).
        input: sample SAM/BAM or VCF
        output: fingerprint file
        Returns:
            list: A list of gatk sample fingerprint jobs.
		"""

        jobs = []
        inputs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.hc.vcf.gz")]
                ]
            )
            inputs.append(input_file)

        output = os.path.join("metrics", "sample.fingerprint")
        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        os.path.join("metrics"),
                        remove=False
                    ),
                    gatk4.crosscheck_fingerprint(
                        inputs,
                        output,
                    )
                ],
                name="gatk_crosscheck_fingerprint.sample.AllSamples",
                samples=self.samples,
                readsets=self.readsets,
                removable_files = []
            )
        )

        return jobs

    def metrics_gatk_cluster_fingerprint(self):
        """
        CheckFingerprint (Picard). Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) 
        against a set of known genotypes in the supplied genotype file (in VCF format).
        input: sample SAM/BAM or VCF
        output: fingerprint file
        Returns:
            list: A list of gatk cluster fingerprint jobs.
        """

        jobs = []

        input_file = os.path.join("metrics", "sample.fingerprint")

        output = os.path.join("metrics", "cluster.fingerprint")
        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        os.path.join("metrics"),
                        remove=False
                    ),
                    gatk4.cluster_crosscheck_metrics(
                        input_file,
                        output,
                    ),
                ],
                name="gatk_cluster_crosscheck_metrics.allSamples",
                samples=self.samples,
                readsets=self.readsets,
                removable_files=[]
            )
        )

        return jobs

    def metrics_verify_bam_id(self):
        """
        [VerifyBamID](https://github.com/Griffan/VerifyBamID) is a software that verifies whether the reads in particular file match previously known
        genotypes for an individual (or group of individuals), and checks whether the reads are contaminated
        as a mixture of two samples. VerifyBamID can detect sample contamination and swaps when external
        genotypes are available. When external genotypes are not available, verifyBamID still robustly
        detects sample swaps.
        Returns:
            list: A list of verify bam id jobs.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            metrics_directory = os.path.join(self.output_dirs['metrics_directory'][sample.name])

            [input_file] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.fixmate.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            bed_file = None
            if coverage_bed is not None:
                bed_file = coverage_bed

            output = os.path.join(metrics_directory, f"{sample.name}.selfSM")
            job_name = f"verify_bam_id.{sample.name}"

            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                        verify_bam_id.parse_contamination_freemix_metrics(output),
                        job2json_project_tracking.run(
                            input_file=output,
                            pipeline=self,
                            samples=sample.name,
                            readsets=",".join([readset.name for readset in sample.readsets]),
                            job_name=job_name,
                            metrics="contamination_freemix=$contamination_freemix"
                        )
                    ])

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            metrics_directory,
                            remove=False
                        ),
                        verify_bam_id.verify2(
                            input_file,
                            os.path.join(metrics_directory, sample.name),
                            bed_file
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    output_dependency=[output],
                    removable_files=[]
                )
            )
            self.multiqc_inputs[sample.name].append(
                    output
            )
        return jobs

    def gatk_readset_fingerprint(self):
        """
    	CheckFingerprint (Picard)
        Checks the sample identity of the sequence/genotype data in the provided file (SAM/BAM or VCF) against a set of known genotypes in the supplied genotype file (in VCF format).
        input: sample SAM/BAM or VCF
        output: fingerprint file
        Returns:
            list: A list of gatk readset fingerprint jobs.
    	"""

        jobs = []
        inputs = []

        for readset in self.readsets:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], readset.sample.name, readset.name)
            input_file = os.path.join(alignment_directory, readset.sample.name, f"{readset.name}.sorted.bam")
            inputs.append(input_file)

            output = os.path.join(self.output_dirs['metrics_directory'][readset.sample.name], f"{readset.name}.fingerprint")

        job = gatk4.crosscheck_fingerprint(inputs, output)
        job.name = "gatk_crosscheck_fingerprint.readset"
        job.samples = self.samples
        job.readsets = self.readsets
        job.removable_files = []
        jobs.append(job)

        return jobs

    def gatk_haplotype_caller(self):
        """
        [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) haplotype caller for snps and small indels.
        Returns:
            list: A list of gatk haplotype caller jobs.
        """

        jobs = []

        reference = global_conf.global_get('gatk_haplotype_caller', 'genome_fasta', param_type='filepath')
        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            interval_directory = os.path.join(alignment_directory, "intervals")
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")

            [input_bam] = self.select_input_files(
                [
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.cram")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(alignment_directory, f"{sample.name}.sorted.bam")]
                ]
            )

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            interval_list = []
            if coverage_bed:
                interval_list = os.path.join(interval_directory, os.path.basename(coverage_bed).replace('.bed', '.noALT.interval_list'))

            elif scatter_jobs == 1:
                interval_list = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))

            if scatter_jobs == 1 or coverage_bed:
                jobs.append(
                    concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                haplotype_directory,
                                remove=True
                            ),
                            gatk4.haplotype_caller(
                                input_bam,
                                os.path.join(haplotype_directory, f"{sample.name}.hc.g.vcf.gz"),
                                interval_list
                            )
                        ],
                        name=f"gatk_haplotype_caller.{sample.name}",
                        samples=[sample],
                        readsets=[*list(sample.readsets)],
                        removable_files=[]
                    )
                )

            else:
                interval_list = [os.path.join(interval_directory, f"{idx:04d}-scattered.interval_list") for idx in range(scatter_jobs)]

                # Create one separate job for each of the first sequences
                for idx, intervals in enumerate(interval_list):
                    jobs.append(
                        concat_jobs(
                            [
                                # Create output directory since it is not done by default by GATK tools
                                bash.mkdir(
                                    haplotype_directory,
                                    remove=True
                                ),
                                gatk4.haplotype_caller(
                                    input_bam,
                                    os.path.join(haplotype_directory, f"{sample.name}.{str(idx)}.hc.g.vcf.gz"),
                                    intervals
                                )
                            ],
                            name=f"gatk_haplotype_caller.{sample.name}.{str(idx)}",
                            samples=[sample],
                            readsets=[*list(sample.readsets)],
                            removable_files=[]
                        )
                    )

        return jobs

    def merge_and_call_individual_gvcf(self):
        """
        Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.
        Returns:
            list: A list of merge and call individual gvcf jobs.
        """

        jobs = []
        reference = global_conf.global_get('gatk_haplotype_caller', 'genome_fasta', param_type='filepath')
        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            interval_directory = os.path.join(alignment_directory, "intervals")
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
            haplotype_file_prefix = os.path.join(haplotype_directory, sample.name)
            output_haplotype_file_prefix = os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name)

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub(r"\.[^.]+$", ".noALT.interval_list", os.path.join(interval_directory, coverage_bed))

            elif scatter_jobs == 1:
                interval_list = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))

            if scatter_jobs == 1 or interval_list is not None:
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
                                global_conf.global_get('gatk_genotype_gvcf', 'options')
                            )
                        ],
                        name=f"merge_and_call_individual_gvcf.call.{sample.name}",
                        samples=[sample],
                        readsets=[*list(sample.readsets)],
                        removable_files=[
                            f"{output_haplotype_file_prefix}.hc.g.vcf.gz",
                            f"{output_haplotype_file_prefix}.hc.g.vcf.gz.tbi"
                        ]
                    )
                )
            else:
                gvcfs_to_merge = [f"{haplotype_file_prefix}.{str(idx)}.hc.g.vcf.gz" for idx in range(scatter_jobs)]

                job = gatk4.cat_variants(
                    gvcfs_to_merge,
                    f"{output_haplotype_file_prefix}.hc.g.vcf.gz"
                )
                job.name = f"merge_and_call_individual_gvcf.merge.{sample.name}"
                job.samples = [sample]
                jobs.append(job)

                job = gatk4.genotype_gvcf(
                    f"{output_haplotype_file_prefix}.hc.g.vcf.gz",
                    f"{output_haplotype_file_prefix}.hc.vcf.gz",
                    global_conf.global_get('gatk_genotype_gvcf', 'options')
                )
                job.name = f"merge_and_call_individual_gvcf.call.{sample.name}"
                job.samples = [sample]
                job.readsets = [*list(sample.readsets)]
                jobs.append(job)

        return jobs

    def combine_gvcf(self):
        """
        Combine the per sample gvcfs of haplotype caller into one main file for all samples.
        Returns:
            list: A list of combine gvcf jobs.
        """

        jobs = []
        nb_haplotype_jobs = global_conf.global_get('gatk_combine_gvcf', 'nb_haplotype', param_type='posint')
        nb_maxbatches_jobs = global_conf.global_get('gatk_combine_gvcf', 'nb_batch', param_type='posint')

        interval_list = None

        coverage_bed = bvatools.resolve_readset_coverage_bed(self.samples[0].readsets[0])
        if coverage_bed:
            interval_list = re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)

        mkdir_job = bash.mkdir(self.output_dirs['variants_directory'])

        # merge all sample in one shot
        if nb_maxbatches_jobs == 1:
            if nb_haplotype_jobs == 1 or interval_list is not None:
                jobs.append(
                    concat_jobs(
                        [
                            mkdir_job,
                            gatk4.combine_gvcf(
                                [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name) +".hc.g.vcf.gz" for sample in self.samples],
                                os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.g.vcf.gz")
                            )
                        ],
                        name="gatk_combine_gvcf.AllSamples",
                        samples=self.samples,
                        readsets=self.readsets
                    )
                )
            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_haplotype_jobs - 1, variant=True)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    jobs.append(
                        concat_jobs(
                            [
                                mkdir_job,
                                gatk4.combine_gvcf(
                                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name) + ".hc.g.vcf.gz" for sample in self.samples],
                                    os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz",
                                    intervals=sequences
                                )
                            ],
                            name=f"gatk_combine_gvcf.AllSample.{str(idx)}",
                            samples=self.samples,
                            readsets=self.readsets,
                            removable_files=[
                                os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz",
                                os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz.tbi"
                            ]
                        )
                    )

                # Create one last job to process the last remaining sequences and 'others' sequences
                job = gatk4.combine_gvcf(
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name) + ".hc.g.vcf.gz" for sample in self.samples],
                    os.path.join(self.output_dirs['variants_directory'], "allSamples.others.hc.g.vcf.gz"),
                    exclude_intervals=unique_sequences_per_job_others
                )
                job.name = "gatk_combine_gvcf.AllSample" + ".others"
                job.removable_files = [
                    os.path.join(self.output_dirs['variants_directory'], "allSamples.others.hc.g.vcf.gz"),
                    os.path.join(self.output_dirs['variants_directory'], "allSamples.others.hc.g.vcf.gz.tbi")
                ]
                job.samples=self.samples
                job.readsets = self.readsets
                jobs.append(job)
        else:
            #Combine samples by batch (pre-defined batches number in ini)
            sample_per_batch = int(math.ceil(len(self.samples)/float(nb_maxbatches_jobs)))
            batch_of_sample = [ self.samples[i:(i+sample_per_batch)] for i in range(0,len(self.samples),sample_per_batch) ]
            cpt = 0
            batches = []
            for batch in batch_of_sample :
                if nb_haplotype_jobs == 1 or interval_list is not None:
                    jobs.append(
                        concat_jobs(
                            [
                                mkdir_job,
                                gatk4.combine_gvcf(
                                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name) + ".hc.g.vcf.gz" for sample in batch],
                                    os.path.join(self.output_dirs['variants_directory'], "allSamples.batch" + str(cpt)+".hc.g.vcf.gz")
                                )
                            ],
                            name="gatk_combine_gvcf.AllSamples.batch" + str(cpt),
                            samples=batch,
                            removable_files=[
                                os.path.join(self.output_dirs['variants_directory'], "allSamples.batch" + str(cpt) + ".hc.g.vcf.gz"),
                                os.path.join(self.output_dirs['variants_directory'], "allSamples.batch" + str(cpt) + ".hc.g.vcf.gz.tbi")
                            ]
                        )
                    )
                else:
                    unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_haplotype_jobs - 1, variant=True)

                    # Create one separate job for each of the first sequences
                    for idx, sequences in enumerate(unique_sequences_per_job):
                        jobs.append(
                            concat_jobs(
                                [
                                    mkdir_job,
                                    gatk4.combine_gvcf(
                                        [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name)+".hc.g.vcf.gz" for sample in batch],
                                        os.path.join(self.output_dirs['variants_directory'], "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.gz",
                                        intervals=sequences
                                    )
                                ],
                                name="gatk_combine_gvcf.AllSample" + ".batch" + str(cpt) + "." + str(idx),
                                samples=batch,
                                removable_files=[
                                    os.path.join(self.output_dirs['variants_directory'], "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.gz",
                                    os.path.join(self.output_dirs['variants_directory'], "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.gz.tbi"
                                ]
                            )
                        )

                    # Create one last job to process the last remaining sequences and 'others' sequences
                    job = gatk4.combine_gvcf(
                        [os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name)+".hc.g.vcf.gz" for sample in batch],
                        os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.gz"),
                        exclude_intervals=unique_sequences_per_job_others
                    )
                    job.name="gatk_combine_gvcf.AllSample" + ".batch" + str(cpt) + ".others"
                    job.samples = batch
                    job.removable_files = [
                        os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.gz"),
                        os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.gz.tbi")
                    ]

                    jobs.append(job)

                batches.append("batch" + str(cpt))
                cpt = cpt + 1

            #Combine batches altogether
            if nb_haplotype_jobs == 1 or interval_list is not None:
                job = gatk4.combine_gvcf(
                    [os.path.join(self.output_dirs['variants_directory'], "allSamples." + batch_idx + ".hc.g.vcf.gz") for batch_idx in batches],
                    os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.g.vcf.gz")
                )
                job.name = "gatk_combine_gvcf.AllSamples.batches"
                job.samples = self.samples

                jobs.append(job)

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = split_by_size(self.sequence_dictionary_variant(), nb_haplotype_jobs - 1, variant=True)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    job=gatk4.combine_gvcf(
                        [os.path.join(self.output_dirs['variants_directory'], "allSamples." + batch_idx + "." + str(idx) + ".hc.g.vcf.gz") for batch_idx in batches],
                        os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz",
                        intervals=sequences
                    )
                    job.name = "gatk_combine_gvcf.AllSample" + "." + str(idx)
                    job.samples = self.samples
                    job.removable_files=[
                        os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz",
                        os.path.join(self.output_dirs['variants_directory'], "allSamples") + "." + str(idx) + ".hc.g.vcf.gz.tbi"
                    ]

                    jobs.append(job)

                # Create one last job to process the last remaining sequences and 'others' sequences
                job = gatk4.combine_gvcf(
                    [os.path.join(self.output_dirs['variants_directory'], "allSamples." + batch_idx + ".others.hc.g.vcf.gz") for batch_idx in batches],
                    os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".others.hc.g.vcf.gz"),
                    exclude_intervals=unique_sequences_per_job_others
                )
                job.name = "gatk_combine_gvcf.AllSample" + ".others"
                job.samples = self.samples
                job.removable_files = [
                    os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".others.hc.g.vcf.gz"),
                    os.path.join(self.output_dirs['variants_directory'], "allSamples" + ".others.hc.g.vcf.gz.tbi")
                ]
                jobs.append(job)

        return jobs


    def merge_and_call_combined_gvcf(self):
        """
        Merges the combined gvcfs and also generates a general vcf containing genotypes.
        Returns:
            list: A list of merge and call combined gvcf jobs.
        """

        jobs = []
        nb_haplotype_jobs = global_conf.global_get('gatk_combine_gvcf', 'nb_haplotype', param_type='posint')
        haplotype_file_prefix = os.path.join(self.output_dirs['variants_directory'], "allSamples")
        output_haplotype = os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.g.vcf.gz")
        output_haplotype_genotyped = os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vcf.gz")

        interval_list = None

        coverage_bed = bvatools.resolve_readset_coverage_bed(self.samples[0].readsets[0])
        if coverage_bed:
            interval_list = re.sub(r"\.[^.]+$", ".interval_list", coverage_bed)

        if nb_haplotype_jobs > 1 and interval_list is None:
            unique_sequences_per_job, _ = split_by_size(self.sequence_dictionary_variant(), nb_haplotype_jobs - 1, variant=True)
            gvcfs_to_merge = [f"{haplotype_file_prefix}.{str(idx)}.hc.g.vcf.gz" for idx in range(len(unique_sequences_per_job))]

            gvcfs_to_merge.append(f"{haplotype_file_prefix}.others.hc.g.vcf.gz")

            job = gatk4.cat_variants(
                gvcfs_to_merge,
                output_haplotype
            )
            job.name = "merge_and_call_combined_gvcf.merge.AllSample"
            job.samples = self.samples
            jobs.append(job)

        job = gatk4.genotype_gvcf(
            output_haplotype,
            output_haplotype_genotyped,
            global_conf.global_get('gatk_genotype_gvcf', 'options')
        )
        job.name = "merge_and_call_combined_gvcf.call.AllSample"
        job.samples = self.samples
        job.readsets = self.readsets
        jobs.append(job)

        return jobs

    def variant_recalibrator(self):
        """
        GATK VariantRecalibrator.
        The purpose of the variant recalibrator is to assign a well-calibrated probability to each variant call in a call set.
        You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
        The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
        between SNP call annotations (QD, MQ, HaplotypeScore, and ReadPosRankSum, for example) and the probability that a SNP
        is a true genetic variant versus a sequencing or data processing artifact. This model is determined adaptively based
        on "true sites" provided as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP
        chip array. This adaptive error model can then be applied to both known and novel variation discovered in the call set
        of interest to evaluate the probability that each call is real. The score that gets added to the INFO field of each variant
        is called the VQSLOD. It is the log odds ratio of being a true variant versus being false under the trained Gaussian mixture model.
        Using the tranche file generated by the previous step the ApplyRecalibration walker looks at each variant's VQSLOD value
        and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
        have their filter field annotated with its tranche level. This will result in a call set that simultaneously is filtered
        to the desired level but also has the information necessary to pull out more variants for a higher sensitivity but a
        slightly lower quality level.
        Returns:
            list: A list of variant recalibrator jobs.
        """

        jobs = []

        #generate the recalibration tranche files
        output_directory = os.path.join(self.output_dirs['variants_directory'])
        recal_snps_other_options = global_conf.global_get('variant_recalibrator', 'tranch_other_options_snps')
        recal_indels_other_options = global_conf.global_get('variant_recalibrator', 'tranch_other_options_indels')
        variant_recal_snps_prefix = os.path.join(output_directory, "allSamples.hc.snps")
        variant_recal_indels_prefix = os.path.join(output_directory, "allSamples.hc.indels")

        mkdir_job = bash.mkdir(output_directory)
        jobs.append(
            concat_jobs(
                [
                    mkdir_job,
                    gatk4.variant_recalibrator(
                        os.path.join(output_directory, "allSamples.hc.vcf.gz"),
                        recal_snps_other_options,
                        variant_recal_snps_prefix + ".recal",
                        variant_recal_snps_prefix + ".tranches",
                        variant_recal_snps_prefix + ".R",
                        small_sample_check=True
                    ),
                    gatk4.variant_recalibrator(
                        os.path.join(output_directory, "allSamples.hc.vcf.gz"),
                        recal_indels_other_options,
                        variant_recal_indels_prefix + ".recal",
                        variant_recal_indels_prefix + ".tranches",
                        variant_recal_indels_prefix + ".R",
                        small_sample_check=True
                    )
                ],
                name="variant_recalibrator.tranch.allSamples",
                samples=self.samples,
                readsets=self.readsets
            )
        )

        #aply the recalibration
        apply_snps_other_options = global_conf.global_get('variant_recalibrator', 'apply_other_options_snps')
        apply_indels_other_options = global_conf.global_get('variant_recalibrator', 'apply_other_options_indels')
        variant_apply_snps_prefix = os.path.join(output_directory, "allSamples.hc.snps")
        variant_apply_indels_prefix = os.path.join(output_directory, "allSamples.hc.indels")

        jobs.append(
            concat_jobs(
                [
                    mkdir_job,
                    gatk4.apply_recalibration(
                        os.path.join(output_directory, "allSamples.hc.vcf.gz"),
                        variant_apply_snps_prefix + ".recal",
                        variant_apply_snps_prefix + ".tranches",
                        apply_snps_other_options,
                        variant_apply_snps_prefix + "_raw_indels.vqsr.vcf.gz"
                    ),
                    gatk4.apply_recalibration(
                        variant_apply_snps_prefix + "_raw_indels.vqsr.vcf.gz",
                        variant_apply_indels_prefix + ".recal",
                        variant_apply_indels_prefix + ".tranches",
                        apply_indels_other_options,
                        os.path.join(output_directory, "allSamples.hc.vqsr.vcf.gz")
                    )
                ],
                name="variant_recalibrator.apply.allSamples",
                samples=self.samples,
                readsets=self.readsets
            )
        )
        return jobs

    def vt_decompose_and_normalize(
        self,
        input_vcf=os.path.join("variants", "allSamples.merged.flt.vcf"),
        output_vcf=os.path.join("variants", "allSamples.merged.flt.vt.vcf.gz"),
        job_name="decompose_and_normalize"
        ):
        """
        Variants with multiple alternate alleles will not be handled correctly by gemini (or by the tools used to annotate the variants).
        To reduce the number of false negatives, the authors of gemini strongly recommend that gemini users split, left-align, and trim their variants.
        For more info on preprocessing, see the gemini docs: https://gemini.readthedocs.io/en/latest/content/preprocessing.html
        The tool used for decomposing and normalizing VCFs is vt: https://github.com/atks/vt
        Arguments:
            input_vcf (str): The input vcf file to decompose and normalize. Default is allSamples.merged.flt.vcf.
            output_vcf (str): The output vcf file. Default is allSamples.merged.flt.vt.vcf.gz.
            job_name (str): The name of the job. Default is decompose_and_normalize.
        Returns:
            list: A list of vt decompose and normalize jobs.
        """

        jobs = []
        jobs.append(
            pipe_jobs(
                [
                    vt.decompose_and_normalize_mnps(
                        input_vcf,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_vcf
                    )
                ],
                name=job_name,
                samples=self.samples,
                readsets=self.readsets
            )
        )
        return jobs

    def haplotype_caller_decompose_and_normalize(self):
        """
        Variants with multiple alternate alleles will not be handled correctly by gemini (or by the tools used to annotate the variants).
        To reduce the number of false negatives, the authors of gemini strongly recommend that gemini users split, left-align, and trim their variants.
        For more info on preprocessing, see the gemini docs: https://gemini.readthedocs.io/en/latest/content/preprocessing.html
        The tool used for decomposing and normalizing VCFs is vt: https://github.com/atks/vt
        Returns:
            list: A list of haplotype caller decompose and normalize jobs.
        """

        input_vcf = self.select_input_files([[os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vcf.gz")]])

        job = self.vt_decompose_and_normalize(
            input_vcf,
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.vcf.gz")
        )
        return job

    def filter_nstretches(
        self,
        input_vcf=os.path.join("variants", "allSamples.merged.flt.vcf"),
        output_vcf=os.path.join("variants", "allSamples.merged.flt.NFiltered.vcf"),
        job_name="filter_nstretches"
        ):
        """
        The final .vcf files are filtered for long 'N' INDELs which are sometimes introduced and cause excessive
        memory usage by downstream tools.
        Arguments:
            input_vcf (str): The input vcf file to filter for long 'N' INDELs. Default is allSamples.merged.flt.vcf.
            output_vcf (str): The output vcf file. Default is allSamples.merged.flt.NFiltered.vcf.
            job_name (str): The name of the job. Default is filter_nstretches.
        Returns:
            list: A list of filter nstretches jobs.
        """

        job = tools.filter_long_indel(input_vcf, output_vcf)
        job.name = job_name
        job.samples = self.samples
        job.readsets = self.readsets
        return [job]

    def haplotype_caller_filter_nstretches(self):
        """
        The final haplotype caller .vcf files are filtered for long 'N' INDELs which are sometimes introduced and cause excessive
        memory usage by downstream tools.
        Returns:
            list: A list of haplotype caller filter nstretches jobs.
        """

        # Find input vcf first from VSQR, then from non recalibrate hapotype calleroriginal BAMs in the readset sheet.
        hc_vcf = self.select_input_files(
            [
                [os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vcf")],
                [os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vcf.gz")]
            ]
        )
        job = self.filter_nstretches(
            hc_vcf[0],
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.NFiltered.vcf"),
            "haplotype_caller_filter_nstretches"
        )
        return job

    def flag_mappability(
        self,
        input_vcf=os.path.join("variants", "allSamples.merged.flt.vt.vcf"),
        output_vcf=os.path.join("variants", "allSamples.merged.flt.vt.mil.vcf"),
        job_name="flag_mappability"
        ):
        """
        Mappability annotation. An in-house database identifies regions in which reads are confidently mapped
        to the reference genome.
        Arguments:
            input_vcf (str): The input vcf file to annotate for mappability. Default is allSamples.merged.flt.vt.vcf.
            output_vcf (str): The output vcf file. Default is allSamples.merged.flt.vt.mil.vcf.
            job_name (str): The name of the job. Default is flag_mappability.
        Returns:
            list: A list of flag mappability jobs.
        """
        jobs = []
        jobs.append(
            pipe_jobs(
                [
                    vcftools.annotate_mappability(input_vcf, None),
                    htslib.bgzip_tabix(None, output_vcf),
                ],
                name=job_name,
                samples = self.samples,
                readsets = self.readsets
            )
        )
        return jobs

    def haplotype_caller_flag_mappability(self):
        """
        Mappability annotation applied to haplotype caller vcf.
        An in-house database identifies regions in which reads are confidently mapped
        to the reference genome.
        Returns:
            list: A list of haplotype caller flag mappability jobs.
        """

        job = self.flag_mappability(
            os.path.join("variants", "allSamples.hc.vqsr.vt.vcf.gz"),
            os.path.join("variants", "allSamples.hc.vqsr.vt.mil.vcf.gz"),
            "haplotype_caller_flag_mappability"
        )
        return job

    def snp_id_annotation(
        self,
        input_vcf=os.path.join("variants", "allSamples.merged.flt.vt.mil.vcf.gz"),
        output_vcf=os.path.join("variants", "allSamples.merged.flt.vt.mil.snpId.vcf.gz"),
        job_name="snp_id_annotation"
        ):
        """
        dbSNP annotation. The .vcf files are annotated for dbSNP using the software SnpSift (from the [SnpEff suite](http://snpeff.sourceforge.net/)).
        Arguments:
            input_vcf (str): The input vcf file to annotate for dbSNP. Default is allSamples.merged.flt.vt.mil.vcf.gz.
            output_vcf (str): The output vcf file. Default is allSamples.merged.flt.vt.mil.snpId.vcf.gz.
            job_name (str): The name of the job. Default is snp_id_annotation.
        Returns:
            list: A list of snp id annotation jobs.
        """

        jobs = []
        jobs.append(
            pipe_jobs(
                [
                    snpeff.snpsift_annotate(
                        input_vcf,
                        None
                    ),
                    htslib.bgzip_tabix(
                        None,
                        output_vcf
                    )
                ],
                name=job_name
            )
        )
        return jobs

    def haplotype_caller_snp_id_annotation(self):
        """
        dbSNP annotation applied to haplotype caller vcf.
        The .vcf files are annotated for dbSNP using the software SnpSift (from the [SnpEff suite](http://snpeff.sourceforge.net/)).
        Returns:
            list: A list of haplotype caller snp id annotation jobs.
        """

        job = self.snp_id_annotation(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.vcf.gz"),
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.vcf.gz"),
            "haplotype_caller_snp_id_annotation"
        )

        return job

    def snp_effect(
            self,
            input_file=os.path.join("variants", "allSamples.hc.vqsr.vt.mil.snpId.vcf.gz"),
            output=os.path.join("variants", "allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf"),
            job_name="snp_effect.hc"
    ):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Arguments:
            input_file (str): The input vcf file to annotate for variant effects. Default is allSamples.hc.vqsr.vt.mil.snpId.vcf.gz.
            output (str): The output vcf file. Default is allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.
            job_name (str): The name of the job. Default is snp_effect.hc.
        Returns:
            list: A list of snp effect jobs.
        """

        jobs = []

        if "high_cov" in self.protocol:
            [input_file] = [os.path.join("variants", "allSamples.vt.vcf.gz")]
            [output] = [os.path.join("variants", "allSamples.vt.snpeff.vcf")]
            job_name = "snp_effect.high_cov"

        jobs.append(
            concat_jobs(
                [
                    snpeff.compute_effects(
                        input_file,
                        output,
                        options=global_conf.global_get('compute_effects', 'options', required=False)
                    ),
                    htslib.bgzip_tabix(
                        output,
                        f"{output}.gz"
                    )
                ],
                name=job_name,
                samples=self.samples,
                readsets=self.readsets
            )
        )

        return jobs

    def haplotype_caller_snp_effect(self):
        """
        Variant effect annotation applied to haplotype caller vcf.
        The .vcf files are annotated for variant effects using the [SnpEff](http://snpeff.sourceforge.net/) software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Returns:
            list: A list of haplotype caller snp effect jobs.
        """

        jobs = self.snp_effect(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.vcf.gz"),
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf"),
            "haplotype_caller_snp_effect"
        )

        return jobs

    def dbnsfp_annotation(
            self,
            input_vcf=os.path.join("variants", "allSamples.merged.flt.vt.mil.snpId.snpeff.vcf.gz"),
            output_vcf=os.path.join("variants", "allSamples.merged.flt.vt.mil.snpId.snpeff.dbnsfp.vcf"),
            ini_section='dbnsfp_annotation',
            job_name="dbnsfp_annotation"
        ):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        Arguments:
            input_vcf (str): The input vcf file to annotate with dbNSFP. Default is allSamples.merged.flt.vt.mil.snpId.snpeff.vcf.gz.
            output_vcf (str): The output vcf file. Default is allSamples.merged.flt.vt.mil.snpId.snpeff.dbnsfp.vcf.
            ini_section (str): The section of the ini file to use for the dbNSFP annotation. Default is dbnsfp_annotation.
            job_name (str): The name of the job. Default is dbnsfp_annotation.
        Returns:
            list: A list of dbnsfp annotation jobs.
        """

        jobs = []
        jobs.append(
            concat_jobs(
                [
                    snpeff.snpsift_dbnsfp(
                        input_vcf,
                        output_vcf,
                        ini_section,
                    ),
                    htslib.bgzip_tabix(
                        output_vcf,
                        f"{output_vcf}.gz"
                    )
                ],
                name=job_name,
                samples=self.samples,
                readsets=self.readsets
            )
        )

        return jobs

    def haplotype_caller_dbnsfp_annotation(self):
        """
        Additional SVN annotations applied to haplotype caller vcf.
        Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        Returns:
            list: A list of haplotype caller dbnsfp annotation jobs.
        """
        job = self.dbnsfp_annotation(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.snpeff.vcf.gz"),
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf"),
            "dbnsfp_annotation"
        )
        return job

    def gemini_annotations(
        self,
        input_file=os.path.join("variants", "allSamples.merged.flt.vt.mil.snpId.snpeff.dbnsfp.vcf.gz"),
        output=os.path.join("variants", "allSamples.gemini.db"),
        job_name="gemini_annotations"
        ):
        """
        Load functionally annotated vcf file into a mysql lite annotation database :
        http://gemini.readthedocs.org/en/latest/index.html
        Arguments:
            input_file (str): The input vcf file to load into a mysql lite annotation database. Default is allSamples.merged.flt.vt.mil.snpId.snpeff.dbnsfp.vcf.gz.
            output (str): The output database file. Default is allSamples.gemini.db.
            job_name (str): The name of the job. Default is gemini_annotations.
        Returns:
            list: A list of gemini annotations jobs.
        """

        if "high_cov" in self.protocol:
            input_file = os.path.join(self.output_dirs['variants_directory'], "allSamples.vt.snpeff.vcf.gz")

        job = gemini.gemini_annotations(
            input_file,
            output,
        )
        job.name = job_name
        job.samples = self.samples
        job.readsets = self.readsets
        return [job]

    def haplotype_caller_gemini_annotations(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database :
        http://gemini.readthedocs.org/en/latest/index.html
        Returns:
            list: A list of haplotype caller gemini annotations jobs.
        """

        job = self.gemini_annotations(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.snpeff.dbnsfp.vcf.gz"),
            os.path.join(self.output_dirs['variants_directory'], "allSamples.gemini.db"),
            job_name="gemini_annotations"
        )
        return job

    def split_tumor_only(self):
        """
        Splits the merged VCF produced in previous steps to generate a report on a per-patient basis.
        The merged VCF is split using the bcftools +split function with the removal of homozygous reference calls.
        Creates one VCF per patient to be used for downstream reporting.
        Returns:
            list: A list of split tumor only jobs.
        """

        jobs = []

        [input_file] = self.select_input_files([[os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.vcf.gz")]])
        output = os.path.join(self.output_dirs['variants_directory'], "split")
        output_files = [os.path.join(output, f"{sample.name}.vcf.gz") for sample in self.samples]
        options = global_conf.global_get('split_tumor_only', 'options')

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(output),
                    bcftools.split(
                        input_file,
                        output,
                        options
                    )
                ],
                name="split_tumor_only",
                samples=self.samples,
                readsets=self.readsets,
                output_dependency=output_files
            )
        )
        return jobs

    def filter_tumor_only(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        Returns:
            list: A list of filter tumor only jobs.
        """
        jobs = []

        output_directory = os.path.join(self.output_dirs['variants_directory'], "split")

        for sample in self.samples:
            input_file = os.path.join(
                output_directory,
                f"{sample.name}.vcf.gz"
            )
            output = os.path.join(
                output_directory,
                sample.name,
                f"{sample.name}.annot.vcf.gz"
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            os.path.join(output_directory, sample.name)
                        ),
                        tools.format2pcgr(
                            input_file,
                            output,
                            global_conf.global_get('filter_tumor_only', 'call_filter'),
                            "somatic",
                            sample.name,
                            ini_section='filter_tumor_only'
                        ),
                        htslib.tabix(
                            output,
                            options="-pvcf"
                        )
                    ],
                    name=f"filter_tumor_only.{sample.name}",
                    samples=[sample],
                    readsets=[*list(sample.readsets)]
                )
            )

        return jobs

    def report_cpsr(self):
        """
        Creates a cpsr germline report (https://sigven.github.io/cpsr/)
        input: annotated/filter vcf
        output: html report and addtional flat files
        Returns:
            list: A list of cpsr report jobs.
        """
        jobs = []

        # Set directory, ini_section, job and sample name for dnaseq tumor only protocol
        if 'tumor_only' in self.protocol:
            output_directory = os.path.join(self.output_dirs['variants_directory'], "split")
            ini_section = 'report_cpsr_tumor_only'

            for sample in self.samples:
                job_name = f"report_cpsr_tumor_only.{sample.name}"
                sample_name = sample.name
                samples = [sample]

                input_file = os.path.join(
                    output_directory,
                    sample.name,
                    f"{sample.name}.annot.vcf.gz"
                )
                cpsr_directory = os.path.join(
                    output_directory,
                    sample.name,
                    "cpsr",
                )
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cpsr_directory,
                            ),
                            cpsr.report(
                                input_file,
                                cpsr_directory,
                                sample_name,
                                ini_section=ini_section
                            )
                        ],
                        name=job_name,
                        samples=samples,
                        readsets=[*list(sample.readsets)],
                    )
                )

        else:
            for tumor_pair in self.tumor_pairs.values():
                #Set directory, ini_section, job and sample name for tumor pair Fastpass protocol
                if 'fastpass' in self.protocol:
                    panel_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
                    ini_section = 'report_cpsr_fastpass'
                    job_name = f"report_cpsr_fastpass.{tumor_pair.name}"

                    input_file = os.path.join(
                        panel_directory,
                        f"{tumor_pair.name}.varscan2.germline.annot.flt.vcf.gz"
                    )
                    cpsr_directory = os.path.join(
                        panel_directory,
                        "cpsr"
                    )
                #Set directory, ini_section, job and sample name for tumor pair ensemble protocol
                elif 'ensemble' in self.protocol:
                    ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
                    ini_section = 'report_cpsr'
                    job_name = f"report_cpsr.{tumor_pair.name}"

                    input_file = os.path.join(
                        ensemble_directory,
                        tumor_pair.name,
                        f"{tumor_pair.name}.ensemble.germline.vt.annot.2caller.flt.vcf.gz"
                    )
                    cpsr_directory = os.path.join(
                        ensemble_directory,
                        tumor_pair.name,
                        "cpsr"
                    )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cpsr_directory,
                            ),
                            cpsr.report(
                                input_file,
                                cpsr_directory,
                                tumor_pair.name,
                                ini_section=ini_section
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def report_pcgr(self):
        """
        Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
        input: filtered somatic vcf
        output: html report and addtional flat files
        Returns:
            list: A list of pcgr report jobs.
        """
        jobs = []

        # Set directory, ini_section, job and sample name for dnaseq tumor only protocol
        if 'tumor_only' in self.protocol:
            output_directory = os.path.join(self.output_dirs['variants_directory'], "split")
            ini_section = 'report_pcgr_tumor_only'
            assembly = global_conf.global_get(ini_section, 'assembly')

            for sample in self.samples:
                cpsr_directory = os.path.join(
                    output_directory,
                    sample.name,
                    "cpsr"
                )
                input_cpsr = os.path.join(
                    cpsr_directory,
                    f"{sample.name}.cpsr.{assembly}.json.gz"
                )
                input_file = os.path.join(
                    output_directory,
                    sample.name,
                    f"{sample.name}.annot.vcf.gz"
                )

                input_cna = os.path.join(
                    self.output_dirs['sv_variants_directory'],
                    f"{sample.name}.cnvkit.cna.tsv"
                )

                if f"{input_cna}.pass":
                    output_cna = os.path.join(
                        self.output_dirs['sv_variants_directory'],
                        f"{sample.name}.cnvkit.cna.tsv"
                    )
                else:
                    output_cna = None

                pcgr_directory = os.path.join(
                    output_directory,
                    sample.name,
                    "pcgr"
                )

                if global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] < "2":
                    output = os.path.join(
                        pcgr_directory,
                        f"{sample.name}.pcgr_acmg.{assembly}.flexdb.html"
                    )
                    input_dependencies = [input_file, input_cpsr, output_cna]
                # output file name patterns have changed in pcgr versions >2.0.0
                # cpsr input files changed from v2.1.0
                else:
                    input_cpsr = os.path.join(
                        cpsr_directory,
                        f"{sample.name}.cpsr.{assembly}"
                    )
                    output = os.path.join(
                        pcgr_directory,
                        f"{sample.name}.pcgr.{assembly}.html"
                    )
                    input_dependencies = [input_file, input_cpsr + ".classification.tsv.gz", input_cpsr + ".conf.yaml", output_cna]
                job_name = f"report_pcgr_tumor_only.{sample.name}"

                pcgr_job = concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory,
                        ),
                        pcgr.report(
                            input_file,
                            input_cpsr,
                            pcgr_directory,
                            sample.name,
                            input_cna=output_cna,
                            ini_section=ini_section
                        ),
                        bash.ls(output)
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=[*list(sample.readsets)],
                    input_dependency=input_dependencies,
                    output_dependency=[output]
                )
                
                if self.project_tracking_json:
                    samples = [sample]
                    pcgr_output_file = os.path.join(self.output_dir, "job_output", "report_pcgr", f"{job_name}_{self.timestamp}.o")
                    jobs.append(
                        concat_jobs(
                            [
                                pcgr_job,
                                pcgr.parse_pcgr_passed_variants_pt(pcgr_output_file),
                                job2json_project_tracking.run(
                                    input_file=pcgr_output_file,
                                    pipeline=self,
                                    samples=",".join([sample.name for sample in samples]),
                                    readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                                    job_name=job_name,
                                    metrics="pcgr_passed_variants=$pcgr_passed_variants"
                                )
                            ], 
                            name=job_name,
                            samples=[sample],
                            readsets=[*list(sample.readsets)],
                            input_dependency=input_dependencies,
                            output_dependency=[output]
                        )
                    )
                else:
                    jobs.append(pcgr_job)
                
        else:
            for tumor_pair in self.tumor_pairs.values():
                # Set directory, ini_section, job and sample name for tumor pair Fastpass protocol
                if 'fastpass'  in self.protocol:
                    panel_directory = os.path.join(
                        self.output_dirs['paired_variants_directory'],
                        tumor_pair.name,
                        "panel"
                    )
                    ini_section = 'report_pcgr_fastpass'
                    assembly = global_conf.global_get(ini_section, 'assembly')
                    job_name=f"report_pcgr_fastpass.{tumor_pair.name}"

                    cpsr_directory = os.path.join(
                        panel_directory,
                        "cpsr"
                    )
                    input_file = os.path.join(
                        panel_directory,
                        f"{tumor_pair.name}.varscan2.somatic.annot.flt.vcf.gz"
                    )
                    pcgr_directory = os.path.join(
                        panel_directory,
                        "pcgr"
                    )
                # Set directory, ini_section, job and sample name for tumor pair Ensemble protocol
                elif 'ensemble' in self.protocol:
                    ensemble_directory = os.path.join(
                        self.output_dirs['paired_variants_directory'],
                        "ensemble"
                    )
                    ini_section = 'report_pcgr'
                    assembly = global_conf.global_get( ini_section, 'assembly')
                    job_name = f"report_pcgr.{tumor_pair.name}"

                    cpsr_directory = os.path.join(
                        ensemble_directory,
                        tumor_pair.name,
                        "cpsr"
                    )
                    input_file = os.path.join(
                        ensemble_directory,
                        tumor_pair.name,
                        f"{tumor_pair.name}.ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
                    )

                    pcgr_directory = os.path.join(
                        ensemble_directory,
                        tumor_pair.name,
                        "pcgr"
                    )

                input_cpsr = os.path.join(
                    cpsr_directory,
                    f"{tumor_pair.name}.cpsr.{assembly}.json.gz"
                )

                input_cna = os.path.join(
                    self.output_dirs['sv_variants_directory'],
                    f"{tumor_pair.name}.cnvkit.cna.tsv"
                )

                if input_cna + ".pass":
                    output_cna = os.path.join(
                        self.output_dirs['sv_variants_directory'],
                        f"{tumor_pair.name}.cnvkit.cna.tsv"
                    )
                else:
                    output_cna = None

                if global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] < "2":
                    output = os.path.join(
                        pcgr_directory,
                        f"{tumor_pair.name}.pcgr_acmg.{assembly}.flexdb.html"
                    )
                    input_dependencies = [input_file, input_cpsr, output_cna]
                # output file name patterns have changed in pcgr versions >2.0.0
                # cpsr input files changed from v2.1.0
                else:
                    input_cpsr = os.path.join(
                        cpsr_directory,
                        f"{tumor_pair.name}.cpsr.{assembly}"
                    )
                    output = os.path.join(
                        pcgr_directory,
                        f"{tumor_pair.name}.pcgr.{assembly}.html"
                    )
                    input_dependencies = [input_file, input_cpsr + ".classification.tsv.gz", input_cpsr + ".conf.yaml", output_cna]

                pcgr_job = concat_jobs(
                    [
                        bash.mkdir(
                            pcgr_directory,
                        ),
                        pcgr.report(
                            input_file,
                            input_cpsr,
                            pcgr_directory,
                            tumor_pair.name,
                            input_cna=output_cna,
                            ini_section=ini_section
                        ),
                        bash.ls(output)
                    ],
                    name=job_name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=input_dependencies,
                    output_dependency=[output]
                )
                
                samples = [tumor_pair.normal, tumor_pair.tumor]
                
                if self.project_tracking_json:
                    pcgr_output_file = os.path.join(self.output_dir, "job_output", "report_pcgr", f"{job_name}_{self.timestamp}.o")
                    jobs.append(
                        concat_jobs(
                            [
                                pcgr_job,
                                pcgr.parse_pcgr_passed_variants_pt(pcgr_output_file),
                                job2json_project_tracking.run(
                                    input_file=pcgr_output_file,
                                    pipeline=self,
                                    samples=",".join([sample.name for sample in samples]),
                                    readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                                    job_name=job_name,
                                    metrics="pcgr_passed_variants=$pcgr_passed_variants"
                                )
                            ], 
                            name=job_name, 
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            input_dependency=input_dependencies,
                            output_dependency=[output]
                        )
                    )
                else:
                    jobs.append(pcgr_job)

        return jobs

    def metrics_vcf_stats(
        self,
        variants_file_prefix="variants/allSamples.hc.vqsr.vt.mil.snpId.snpeff",
        job_name="metrics_change_rate"
        ):
        """
        Metrics SNV. Multiple metrics associated to annotations and effect prediction are generated at this step:
        change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect,
        counts by genomic region, SNV quality, coverage, InDel lengths, base changes,  transition-transversion rates,
        summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.
        Returns:
            list: A list of metrics vcf stats jobs.
        """

        job = metrics.vcf_stats(
            f"{variants_file_prefix}.dbnsfp.vcf.gz",
            f"{variants_file_prefix}.dbnsfp.vcf.part_changeRate.tsv",
            f"{variants_file_prefix}.vcf.stats.csv"
        )
        job.name = job_name
        job.samples = self.samples

        return [job]

    def haplotype_caller_metrics_vcf_stats(self):
        """
        Metrics SNV applied to haplotype caller vcf.
        Multiple metrics associated to annotations and effect prediction are generated at this step:
        change rate by chromosome, changes by type, effects by impact, effects by functional class, counts by effect,
        counts by genomic region, SNV quality, coverage, InDel lengths, base changes,  transition-transversion rates,
        summary of allele frequencies, codon changes, amino acid changes, changes per chromosome, change rates.
        Returns:
            list: A list of haplotype caller metrics vcf stats jobs.
        """

        job = self.metrics_vcf_stats(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.vt.mil.snpId.snpeff"),
            "haplotype_caller_metrics_change_rate"
        )
        return job

    def metrics_snv_graph_metrics(
        self,
        variants_file_prefix=os.path.join("variants", "allSamples.merged.flt.mil.snpId"),
        snv_metrics_prefix=os.path.join("variants", "allSamples.SNV"),
        job_name="metrics_snv_graph"
        ):
        """
        Returns:
            list: A list of metrics snv graph metrics jobs.
        """

        report_file = f"{self.output_dirs['report_directory']}/DnaSeq.metrics_snv_graph_metrics.md"
        snv_metrics_files = [
            f"{snv_metrics_prefix}.SummaryTable.tsv",
            f"{snv_metrics_prefix}.EffectsFunctionalClass.tsv",
            f"{snv_metrics_prefix}.EffectsImpact.tsv"
        ]

        job = metrics.snv_graph_metrics(
            f"{snv_metrics_prefix}.snpeff.vcf.statsFile.txt",
            snv_metrics_prefix
        )
        job.output_files = snv_metrics_files

        return [
            concat_jobs(
                [
                    job,
                    Job(
                        snv_metrics_files,
                        [report_file],
                        [
                            [job_name + "_report",
                            'module_pandoc']
                        ],
                        command=f"""\
mkdir -p report && \\
paste \\
  <(echo -e "Number of variants before filter\nNumber of not variants\n%\nNumber of variants processed\nNumber of known variants\n%\nTransitions\nTransversions\nTs Tv ratio\nmissense\nnonsense\nsilent\nmissense silent ratio\nhigh impact\nlow impact\nmoderate impact\nmodifier impact") \\
  <(paste \\
    {snv_metrics_prefix}.SummaryTable.tsv \\
    {snv_metrics_prefix}.EffectsFunctionalClass.tsv \\
    <(sed '1d' {snv_metrics_prefix}.EffectsImpact.tsv) \\
  | sed '1d' | sed 's/\t/\\n/g') \\
  > report/SNV.SummaryTable.tsv
snv_summary_table_md=`sed 's/\t/|/g' report/SNV.SummaryTable.tsv`
pandoc \\
  {self.report_template_dir}/{os.path.basename(report_file)} \\
  --template {self.report_template_dir}/{os.path.basename(report_file)} \\
  --variable snv_summary_table="$snv_summary_table_md" \\
  --to markdown \\
  e {report_file}
for file in SNVQuality  IndelLength CountRegions CountEffects BaseChange codonChange AminoAcidChange changeRate TsTv
do
  for ext in jpeg pdf tsv
  do
  cp \\
    {snv_metrics_prefix}.$file.$ext  \\
    report/SNV.$file.$ext
  done
done
cp {snv_metrics_prefix}.chromosomeChange.zip report/SNV.chromosomeChange.zip""",
                        report_files=[report_file]
                    )
                ],
                name=f"{job_name}_report",
                samples=self.samples,
                readsets=self.readsets
            )
        ]

    def haplotype_caller_metrics_snv_graph_metrics(self):
        """
        See general metrics_vcf_stats !  Applied to haplotype caller vcf.
        Returns:
            list: A list of haplotype caller metrics snv graph metrics jobs.
        """

        jobs = self.metrics_snv_graph_metrics(
            os.path.join(self.output_dirs['variants_directory'], "allSamples.hc.vqsr.mil.snpId"),
            f"{self.output_dirs['metrics_directory']}/allSamples.hc.vqsr.SNV",
            "haplotype_caller_metrics_snv_graph"
        )
        return jobs

    def delly_call_filter(self):
        """
        Delly2 is an integrated structural variant prediction method that can
        discover, genotype and visualize deletions, tandem duplications, inversions and translocations
        at single-nucleotide resolution in short-read massively parallel sequencing data. It uses paired-ends
        and split-reads to sensitively and accurately delineate genomic rearrangements throughout the genome.
        Structural variants can be visualized using Delly-maze and Delly-suave.
        Input: normal and tumor final bams
        Output: bcl file
        Returns:
            list: A list of delly call filter jobs.
        """

        jobs = []
        for sample in self.samples:

            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")

            [input_file] = self.select_input_files(
                [
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.recal.bam")],
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.bam")]
                ]
            )

            sv_types = global_conf.global_get('delly_call_filter', 'sv_types_options').split(",")

            for sv_type in sv_types:
                output_bcf = os.path.join(delly_directory, f"{sample.name}.delly.{str(sv_type)}.bcf")
                output_vcf = os.path.join(delly_directory, f"{sample.name}.delly.{str(sv_type)}.germline.flt.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                delly_directory,
                                remove=True
                            ),
                            delly.call(
                                [input_file],
                                output_bcf,
                                sv_type
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_bcf,
                                        None,
                                        global_conf.global_get('delly_call_filter_germline', 'bcftools_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_vcf
                                    )
                                ]
                            )
                        ],
                        name=f"delly_call_filter.{str(sv_type)}.{sample.name}",
                        samples=[sample]
                    )
                )

        return jobs

    def delly_sv_annotation(self):
        """
        Preprocess and annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        Returns:
            list: A list of delly sv annotation jobs.
        """
        jobs = []

        for sample in self.samples:

            directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            final_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name, sample.name)
            delly_directory = os.path.join(directory, "rawDelly")
            output_vcf = os.path.join(delly_directory, f"{sample.name}.delly.merge.sort.vcf.gz")
            germline_vcf = os.path.join(directory, f"{sample.name}.delly.germline.vcf.gz")

            sv_types = global_conf.global_get('delly_call_filter', 'sv_types_options').split(",")

            input_bcf = []
            for sv_type in sv_types:
                input_bcf.append(os.path.join(delly_directory, f"{sample.name}.delly.{str(sv_type)}.bcf"))

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                bcftools.concat(
                                    input_bcf,
                                    None,
                                    "-O v"
                                ),
                                vt.sort(
                                    "-",
                                    "-",
                                    "-m full"
                                ),
                                htslib.bgzip(
                                    None,
                                    output_vcf
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vawk.single_germline(
                                    output_vcf,
                                    sample.name,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    germline_vcf
                                )
                            ]
                        )
                    ],
                    name=f"sv_annotation.delly.merge_sort_filter.{sample.name}",
                    samples=[sample]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            germline_vcf,
                            f"{final_directory}.delly.germline.snpeff.vcf.gz"
                        )
                    ],
                    name=f"sv_annotation.delly.germline.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def germline_manta(self):
        """
        Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
        analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
        Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
        single efficient workflow.
        Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences in VCF 4.1 format.
        Returns:
            list: A list of germline manta jobs.
        """
        jobs = []

        bed_file = None
        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            manta_directory = os.path.join(pair_directory, "rawManta")
            output_prefix = os.path.join(pair_directory, sample.name)

            [input_file] = self.select_input_files(
                [
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.bam")]
                ]
            )

            manta_germline_output = os.path.join(manta_directory, "results", "variants", "diploidSV.vcf.gz")
            manta_germline_output_tbi = os.path.join(manta_directory, "results", "variants", "diploidSV.vcf.gz.tbi")

            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

            if coverage_bed and not bed_file:
                bed_file = f"{coverage_bed}.gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.sort(
                                coverage_bed,
                                f"{coverage_bed}.sort",
                                "-V -k1,1 -k2,2n -k3,3n",
                                extra="; sleep 180"
                            ),
                            htslib.bgzip(
                                f"{coverage_bed}.sort",
                                bed_file
                            ),
                            htslib.tabix(
                                f"{coverage_bed}.gz",
                                options="-f -p bed"
                            )
                        ],
                        name=f"bed_index.{sample.name}",
                        samples=[sample]
                    )
                )

            output_dep = [manta_germline_output, manta_germline_output_tbi]

            jobs.append(
                concat_jobs(
                    [
                        bash.rm(manta_directory),
                        bash.mkdir(manta_directory, remove=True),
                        manta.manta_config(
                            input_file,
                            None,
                            manta_directory,
                            bed_file
                        ),
                        manta.manta_run(
                            manta_directory,
                            output_dep=output_dep
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output, os.path.dirname(f"{output_prefix}.manta.germline.vcf.gz")),
                            f"{output_prefix}.manta.germline.vcf.gz",
                            input=manta_germline_output,
                            remove=False
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output_tbi, os.path.dirname(f"{output_prefix}.manta.germline.vcf.gz.tbi")),
                            f"{output_prefix}.manta.germline.vcf.gz.tbi",
                            input=manta_germline_output_tbi,
                            remove=False
                        )
                    ],
                    input_dependency=[input_file],
                    name=f"germline_manta.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def manta_sv_annotation(self):
        """
        Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        Returns:
            list: A list of manta sv annotation jobs.
        """
        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name, sample.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            f"{pair_directory}.manta.germline.vcf.gz",
                            f"{pair_directory}.manta.germline.snpeff.vcf.gz"
                        )
                    ],
                    name=f"sv_annotation.manta_germline.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def lumpy_paired_sv(self):
        """
        A probabilistic framework for structural variant discovery.
        Lumpy traditional with paired ends and split reads on tumor normal pair.
        Outputs: bams
        Returns:
            list: A list of lumpy paired sv jobs.
        """
        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            lumpy_directory = os.path.join(pair_directory, "rawLumpy")
            input_normal = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")

            discordants_normal = os.path.join(lumpy_directory, f"{sample.name}.discordants.sorted.bam")

            splitters_normal = os.path.join(lumpy_directory, f"{sample.name}.splitters.sorted.bam")

            output_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.vcf")
            gzip_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.vcf.gz")

            genotype_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.genotyped.vcf")
            genotype_gzip = os.path.join(pair_directory, f"{sample.name}.lumpy.genotyped.vcf.gz")
            germline_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.germline.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(lumpy_directory, remove=True),
                        pipe_jobs(
                            [
                                samtools.view(
                                    input_normal,
                                    None,
                                    "-b -F 1294"
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    discordants_normal,
                                    lumpy_directory,
                                    global_conf.global_get('extract_discordant_reads', 'sambamba_options')
                                ),
                            ]
                        ),
                    ],
                    name=f"extract_discordant_reads.{sample.name}",
                    samples=[sample]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(lumpy_directory, remove=True),
                        pipe_jobs(
                            [
                                samtools.view(
                                    input_normal,
                                    None,
                                    "-h"
                                ),
                                Job(
                                    [None],
                                    [None],
                                    [
                                        ['extract_split_reads', 'module_python'],
                                        ['extract_split_reads', 'module_lumpy']
                                    ],
                                    command="$LUMPY_SCRIPTS/extractSplitReads_BwaMem -i stdin"
                                ),
                                samtools.view(
                                    "-",
                                    None,
                                    " -Sb "
                                ),
                                sambamba.sort(
                                    "/dev/stdin",
                                    splitters_normal,
                                    lumpy_directory,
                                    global_conf.global_get('extract_split_reads', 'sambamba_options')
                                ),
                            ]
                        ),
                    ],
                    name=f"extract_split_reads.{sample.name}",
                    samples=[sample]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(lumpy_directory, remove=True),
                        lumpy.lumpyexpress_pair(input_normal, None, output_vcf, spl_normal=splitters_normal, dis_normal=discordants_normal),
                        htslib.bgzip(output_vcf, gzip_vcf),
                    ],
                    name=f"lumpy_paired_sv_calls.{sample.name}",
                    samples=[sample]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                Job(
                                    [gzip_vcf],
                                    [None],
                                    command=f"zcat {gzip_vcf} | grep -v \"^##contig\""
                                ),
                                bcftools.annotate(
                                    None,
                                    None,
                                    global_conf.global_get('lumpy_paired_sv_calls', 'header_options')
                                ),
                                vt.sort(
                                    "-",
                                    os.path.join(pair_directory, f"{sample.name}.lumpy.sorted.vcf"),
                                    "-m full"
                                )
                            ]
                        ),
                        svtyper.genotyper(
                            None,
                            input_normal,
                            os.path.join(pair_directory, f"{sample.name}.lumpy.sorted.vcf"),
                            genotype_vcf,
                            ini_section="lumpy_paired_sv_calls"
                        ),
                        htslib.bgzip(
                            genotype_vcf,
                            genotype_gzip
                        ),
                        pipe_jobs(
                            [
                                vawk.single_germline(genotype_gzip, sample.name, None),
                                vt.sort("-", "-", "-m full"),
                                htslib.bgzip_tabix(None, germline_vcf),
                            ]
                        )
                    ],
                    name=f"lumpy_paired_sv_calls.genotype.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def lumpy_sv_annotation(self):
        """
        Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        Returns:
            list: A list of lumpy sv annotation jobs.
        """

        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            prefix = os.path.join(self.output_dirs["sv_variants_directory"], sample.name, sample.name)
            germline_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.germline.vcf.gz")

            snppeff_job = snpeff.compute_effects(
                germline_vcf,
                f"{prefix}.lumpy.germline.snpeff.vcf.gz"
            )
            snppeff_job.name = f"sv_annotation.lumpy.germline.{sample.name}"
            snppeff_job.samples = [sample]

            jobs.append(snppeff_job)

        return jobs

    def wham_call_sv(self):
        """
        Wham (Whole-genome Alignment Metrics) to provide a single, integrated framework for both structural variant
        calling and association testing, thereby bypassing many of the difficulties that currently frustrate attempts
        to employ SVs in association testing.
        Outputs: vcf
        Returns:
            list: A list of wham call sv jobs.
        """
        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            wham_directory = os.path.join(pair_directory, "rawWham")
            input_normal = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")

            output_vcf = os.path.join(wham_directory, f"{sample.name}.wham.vcf")
            merge_vcf = os.path.join(wham_directory, f"{sample.name}.wham.merged.vcf")
            genotyped_vcf = os.path.join(pair_directory, f"{sample.name}.wham.merged.genotyped.vcf.gz")
            germline_vcf = os.path.join(pair_directory, f"{sample.name}.wham.germline.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(wham_directory, remove=True),
                        wham.call_sv(input_normal, None, output_vcf),
                        pipe_jobs(
                            [
                                wham.merge(output_vcf, None),
                                Job(
                                    [None],
                                    [merge_vcf],
                                    command=f"sed 's/NONE/{sample.name}/g' | sed -e 's#\"\"#\"#g' > {merge_vcf}"
                                )
                            ]
                        )
                    ],
                    name=f"wham_call_sv.call_merge.{sample.name}",
                    samples=[sample]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(wham_directory, remove=True),
                        pipe_jobs(
                            [
                                Job(
                                    [merge_vcf],
                                    [None],
                                    command=f"cat {merge_vcf} | grep -v \"^##contig\""
                                ),
                                bcftools.annotate(
                                    None,
                                    None,
                                    global_conf.global_get('wham_call_sv', 'header_options')
                                ),
                                vt.sort(
                                    "-",
                                    os.path.join(pair_directory, f"{sample.name}.wham.sorted.vcf"),
                                    "-m full"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                svtyper.genotyper(
                                    None,
                                    input_normal,
                                    os.path.join(pair_directory, f"{sample.name}.wham.sorted.vcf"),
                                    None
                                ),
                                Job(
                                    [None],
                                    [None],
                                    command="sed -e 's#\"\"#\"#g'"
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    genotyped_vcf
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vawk.single_germline(genotyped_vcf, sample.name, None),
                                htslib.bgzip_tabix(None, germline_vcf),
                            ]
                        )
                    ],
                    name=f"wham_call_sv.genotype.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def wham_sv_annotation(self):
        """
        Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        Returns:
            list: A list of wham sv annotation jobs.
        """

        jobs = []
        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            germline_vcf = os.path.join(pair_directory, f"{sample.name}.wham.germline.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(germline_vcf, os.path.join(pair_directory, f"{sample.name}.wham.germline.snpeff.vcf.gz")),
                    ],
                    name=f"sv_annotation.wham.germline.{sample.name}",
                    samples=[sample]
                )
            )
        return jobs

    def cnvkit_batch(self):
        """
        [CNVkit](https://cnvkit.readthedocs.io/en/stable/index.html) is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data.
        Returns:
            list: A list of cnvkit batch jobs.
        """
        jobs = []

        if 'germline' in self.protocol or 'tumor_only' in self.protocol:
            for sample in self.samples:
                input_prefix = os.path.join(self.output_dirs['alignment_directory'], sample.name, sample.name)

                input_normal = None
                [input_tumor] = self.select_input_files(
                    [
                        [f"{input_prefix}.sorted.dup.cram"],
                        [f"{input_prefix}.sorted.dup.bam"],
                        [f"{input_prefix}.sorted.bam"],
                    ]
                )
                sample_name = sample.name
                samples = [sample]
                readsets = [*list(sample.readsets)]

                input_vcf = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.hc.vcf.gz")
                flt_vcf = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.hc.flt.vcf.gz")

                # Set sample_id for cnvkit export function
                sample_id =  sample.name
                normal_id = None

                pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
                cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
                tarcov_cnn = os.path.join(cnvkit_dir, f"{sample.name}.sorted.dup.targetcoverage.cnn")
                antitarcov_cnn = os.path.join(cnvkit_dir, f"{sample.name}.sorted.dup.antitargetcoverage.cnn")

                if 'germline_sv' in self.protocol or 'tumor_only' in self.protocol:
                    filter_options = "-i '%QUAL>=50' -m2 -M2 -v snps"
                else:
                    filter_options = "-f PASS -i '%QUAL>=50' -m2 -M2 -v snps"

                # Set coverage bed if using exome
                coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

                bed = None
                if coverage_bed:
                    bed = coverage_bed

                ref_cnn = os.path.join(cnvkit_dir, sample_name + ".reference.cnn")
                metrics_folder = os.path.join(self.output_dirs['sv_variants_directory'], "cnvkit_reference")
                pool_ref = os.path.join(self.output_dir, metrics_folder, "pooledReference.cnn")

                if os.path.isfile(pool_ref):
                    pool_ref_cnn = pool_ref
                    ref_cnn = None

                else:
                    pool_ref_cnn = None

                vcf_gz = os.path.join(pair_directory, f"{sample.name}.cnvkit.vcf.gz")

                jobs.append(
                    pipe_jobs(
                        [
                            bcftools.view(
                                input_vcf,
                                None,
                                filter_options=filter_options
                            ),
                            bash.sed(
                                None,
                                None,
                                r"-e 's/^\#\#INFO=<ID=AF,Number=A,.*\">/##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency of the ALT allele\">/'"
                            ),
                            htslib.bgzip_tabix(
                                None,
                                flt_vcf
                            )
                        ],
                        name=f"cnvkit_batch.vcf_flt.{sample_name}",
                        samples=samples,
                        readsets=readsets
                    )
                )

                call_cns = os.path.join(cnvkit_dir, f"{sample.name}.call.cns")

                input_cna = os.path.join(self.output_dirs['sv_variants_directory'], sample_name, f"{sample.name}.cnvkit.vcf.gz")
                header = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample.name}.header")
                output_cna_body = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample.name}.cnvkit.body.tsv")
                output_cna = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample.name}.cnvkit.cna.tsv")
                output_check = f"{output_cna}.pass"

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cnvkit_dir,
                                remove=True
                            ),
                            cnvkit.batch(
                                input_tumor,
                                input_normal,
                                cnvkit_dir,
                                tar_dep=tarcov_cnn,
                                antitar_dep=antitarcov_cnn,
                                target_bed=bed,
                                reference=pool_ref_cnn,
                                output_cnn=ref_cnn
                            )
                        ],
                        name=f"cnvkit_batch.{sample_name}",
                        samples=samples,
                        readsets=readsets
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cnvkit_dir,
                                remove=True
                            ),
                            cnvkit.fix(
                                tarcov_cnn,
                                antitarcov_cnn,
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                reference=pool_ref_cnn,
                                ref_cnn=ref_cnn
                            ),
                            cnvkit.segment(
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                os.path.join(cnvkit_dir, f"{sample.name}.cns"),
                                vcf=flt_vcf,
                                sample_id=sample_id,
                                normal_id=normal_id
                            ),
                            cnvkit.segmetrics(
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                os.path.join(cnvkit_dir, f"{sample.name}.cns"),
                                os.path.join(cnvkit_dir, f"{sample.name}.seg.cns"),
                            ),
                        ],
                        name=f"cnvkit_batch.correction.{sample_name}",
                        samples=samples,
                        readsets=readsets
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cnvkit_dir,
                                remove=True
                            ),
                            cnvkit.call(
                                os.path.join(cnvkit_dir, f"{sample.name}.seg.cns"),
                                call_cns
                            ),
                            pipe_jobs(
                                [
                                    cnvkit.export(
                                        call_cns,
                                        None,
                                        sample_id=sample_id
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        vcf_gz
                                    )
                                ]
                            ),
                        ],
                        name=f"cnvkit_batch.call.{sample.name}",
                        samples=samples,
                        readsets=readsets
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            pcgr.create_header(
                                header
                            ),
                            bcftools.query(
                                input_cna,
                                output_cna_body,
                                query_options="-f '%CHROM\\t%POS\\t%END\\t%FOLD_CHANGE_LOG\\n'"
                            ),
                            pcgr.create_input_cna(
                                output_cna_body,
                                header,
                                call_cns,
                                output_cna
                            ),
                            cnvkit.file_check(
                                output_cna,
                                output_check
                            )
                        ],
                        name=f"cnvkit_batch.cna.{sample_name}",
                        samples=samples,
                        readsets=readsets,
                        input_dependency=[input_cna, output_cna],
                        output_dependency=[header, output_cna_body, output_cna]
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                cnvkit_dir,
                                remove=False
                            ),
                            cnvkit.metrics(
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                os.path.join(cnvkit_dir, f"{sample.name}.call.cns"),
                                os.path.join(metrics_folder, f"{sample.name}.metrics.tsv")
                            ),
                            cnvkit.scatter(
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                os.path.join(cnvkit_dir, f"{sample.name}.call.cns"),
                                os.path.join(cnvkit_dir, f"{sample.name}.scatter.pdf"),
                                vcf=flt_vcf,
                                normal=normal_id,
                                tumor=sample_id
                            ),
                            cnvkit.diagram(
                                os.path.join(cnvkit_dir, f"{sample.name}.cnr"),
                                os.path.join(cnvkit_dir, f"{sample.name}.call.cns"),
                                os.path.join(cnvkit_dir, f"{sample.name}.diagram.pdf")
                            )
                        ],
                        name=f"cnvkit_batch.metrics.{sample_name}",
                        samples=samples,
                        readsets=readsets,
                        removable_files=[cnvkit_dir]
                    )
                )

        else:
            for tumor_pair in self.tumor_pairs.values():
                if tumor_pair.multiple_normal == 1:
                    normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                else:
                    normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

                tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

                [input_normal] = self.select_input_files(
                    [
                        [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                        [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.cram")]
                    ]
                )
                [input_tumor] = self.select_input_files(
                    [
                        [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                        [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.cram")]
                    ]
                )

                sample_name = tumor_pair.name

                if 'fastpass' in self.protocol:
                    input_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel", f"{tumor_pair.name}.varscan2.germline.vt.vcf.gz")
                    flt_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel", f"{tumor_pair.name}.varscan2.germline.flt.vcf.gz")
                else:
                    input_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, f"{tumor_pair.name}.vardict.germline.vt.vcf.gz")
                    flt_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, f"{tumor_pair.name}.vardict.germline.flt.vcf.gz")

                sample_id = tumor_pair.tumor.name
                normal_id = tumor_pair.normal.name

                pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample_name)
                cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
                tarcov_cnn = os.path.join(cnvkit_dir, f"{tumor_pair.tumor.name}.sorted.dup.targetcoverage.cnn")
                antitarcov_cnn = os.path.join(cnvkit_dir, f"{tumor_pair.tumor.name}.sorted.dup.antitargetcoverage.cnn")

                # Set coverage bed if using exome
                coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

                bed = None
                if coverage_bed:
                    bed = coverage_bed

            ref_cnn = os.path.join(cnvkit_dir, f"{sample_name}.reference.cnn")
            metrics_folder = os.path.join(self.output_dirs['sv_variants_directory'], "cnvkit_reference")
            pool_ref = os.path.join(self.output_dir, metrics_folder, "pooledReference.cnn")

            if os.path.isfile(pool_ref):
                pool_ref_cnn = pool_ref
                ref_cnn = None

            else:
                pool_ref_cnn = None

            vcf_gz = os.path.join(pair_directory, f"{sample_name}.cnvkit.vcf.gz")

            jobs.append(
                pipe_jobs(
                    [
                        bcftools.view(
                            input_vcf,
                            None,
                            filter_options="-f PASS -i '%QUAL>=50' -m2 -M2 -v snps"
                        ),
                        bash.sed(
                            None,
                            None,
                            r"-e 's/^\#\#INFO=<ID=AF,Number=A,.*\">/##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency of the ALT allele\">/'"
                        ),
                        htslib.bgzip_tabix(
                            None,
                            flt_vcf
                        )
                    ],
                    name=f"cnvkit_batch.vcf_flt.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            call_cns = os.path.join(cnvkit_dir, f"{sample_name}.call.cns")

            input_cna = os.path.join(self.output_dirs['sv_variants_directory'], sample_name, f"{sample_name}.cnvkit.vcf.gz")
            header = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample_name}.header")
            output_cna_body = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample_name}.cnvkit.body.tsv")
            output_cna = os.path.join(self.output_dirs['sv_variants_directory'], f"{sample_name}.cnvkit.cna.tsv")
            output_check = f"{output_cna}.pass"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.batch(
                            input_tumor,
                            input_normal,
                            cnvkit_dir,
                            tar_dep=tarcov_cnn,
                            antitar_dep=antitarcov_cnn,
                            target_bed=bed,
                            reference=pool_ref_cnn,
                            output_cnn=ref_cnn
                        )
                    ],
                    name=f"cnvkit_batch.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=False
                        ),
                        cnvkit.fix(
                            tarcov_cnn,
                            antitarcov_cnn,
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            reference=pool_ref_cnn,
                            ref_cnn=ref_cnn
                        ),
                        cnvkit.segment(
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            os.path.join(cnvkit_dir, f"{sample_name}.cns"),
                            vcf=flt_vcf,
                            sample_id=sample_id,
                            normal_id=normal_id
                        ),
                        cnvkit.segmetrics(
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            os.path.join(cnvkit_dir, f"{sample_name}.cns"),
                            os.path.join(cnvkit_dir, f"{sample_name}.seg.cns"),
                        )
                    ],
                    name=f"cnvkit_batch.correction.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    removable_files=[cnvkit_dir]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.call(
                            os.path.join(cnvkit_dir, f"{sample_name}.seg.cns"),
                            call_cns
                        ),
                        pipe_jobs(
                            [
                                cnvkit.export(
                                    os.path.join(cnvkit_dir, f"{sample_name}.call.cns"),
                                    None,
                                    sample_id=sample_id
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    vcf_gz
                                )
                            ]
                        ),
                    ],
                    name=f"cnvkit_batch.call.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pcgr.create_header(
                            header,
                        ),
                        bcftools.query(
                            input_cna,
                            output_cna_body,
                            query_options="-f '%CHROM\\t%POS\\t%END\\t%FOLD_CHANGE_LOG\\n'"
                        ),
                        pcgr.create_input_cna(
                            output_cna_body,
                            header,
                            call_cns,
                            output_cna
                        ),
                        cnvkit.file_check(
                            output_cna,
                            output_check
                            )
                    ],
                    name=f"cnvkit_batch.cna.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_cna, output_cna],
                    output_dependency=[header, output_cna_body, output_cna],
                    removable_files=[header, output_cna_body]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cnvkit_dir,
                            remove=True
                        ),
                        cnvkit.metrics(
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            os.path.join(cnvkit_dir, f"{sample_name}.call.cns"),
                            os.path.join(metrics_folder, f"{sample_name}.metrics.tsv")
                        ),
                        cnvkit.scatter(
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            os.path.join(cnvkit_dir, f"{sample_name}.call.cns"),
                            os.path.join(cnvkit_dir, f"{sample_name}.scatter.pdf"),
                            vcf=flt_vcf,
                            normal=normal_id,
                            tumor=sample_id
                        ),
                        cnvkit.diagram(
                            os.path.join(cnvkit_dir, f"{sample_name}.cnr"),
                            os.path.join(cnvkit_dir, f"{sample_name}.call.cns"),
                            os.path.join(cnvkit_dir, f"{sample_name}.diagram.pdf")
                        )
                    ],
                    name=f"cnvkit_batch.metrics.{sample_name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def cnvkit_sv_annotation(self):
        """
        Annotate VCF with SnpEff.
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        https://pcingola.github.io/SnpEff/se_introduction/
        Returns:
            list: A list of cnvkit sv annotation jobs.
        """
        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name, sample.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            f"{pair_directory}.cnvkit.vcf.gz",
                            f"{pair_directory}.cnvkit.snpeff.vcf"),
                    ],
                    name=f"sv_annotation.cnvkit.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def run_breakseq2(self):
        """
        [BreakSeq2](https://bioinform.github.io/breakseq2/): Ultrafast and accurate nucleotide-resolution analysis of structural variants.
        Returns:
            list: A list of breakseq2 jobs.
        """

        jobs = []
        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            output_dir = os.path.join(pair_directory, "rawBreakseq2")
            output = os.path.join(pair_directory, "rawBreakseq2", "breakseq.vcf.gz")
            final_vcf = os.path.join(pair_directory, f"{sample.name}.breakseq.germline.vcf.gz")

            [input_file] = self.select_input_files(
                [
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.bam")]
                ]
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_dir, remove=True),
                        breakseq2.run(
                            input_file,
                            sample.name,
                            output_dir
                        ),
                        pipe_jobs(
                            [
                                bcftools.view(
                                    output,
                                    None,
                                    global_conf.global_get('run_breakseq2', 'bcftools_options')
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    final_vcf
                                )
                            ]
                        )
                    ],
                    name=f"run_breakseq2.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def ensemble_metasv(self):
        """
        [MetaSV](http://bioinform.github.io/metasv/) is an integrated SV caller which leverages multiple orthogonal SV signals for high accuracy and resolution.
        MetaSV proceeds by merging SVs from multiple tools for all types of SVs.
        Returns:
            list: A list of ensemble metasv jobs.
        """
        jobs = []

        for sample in self.samples:
            pair_directory = os.path.join(self.output_dirs["sv_variants_directory"], sample.name)
            ensemble_directory = os.path.join(self.output_dirs["sv_variants_directory"], "ensemble", sample.name)

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.dup.bam")],
                    [os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.sorted.bam")]
                ]
            )

            isize_file = os.path.join(self.output_dirs['metrics_directory'][sample.name], f"{sample.name}.all.metrics.insert_size_metrics")
            gatk_vcf = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.hc.vcf.gz")
            gatk_pass = os.path.join(self.output_dirs['alignment_directory'], sample.name, f"{sample.name}.hc.flt.vcf.gz")
            lumpy_vcf = os.path.join(pair_directory, f"{sample.name}.lumpy.germline.vcf.gz")
            manta_vcf = os.path.join(pair_directory, f"{sample.name}.manta.germline.vcf.gz")
            # abs_manta = os.path.abspath(manta_vcf)
            wham_vcf = os.path.join(pair_directory, f"{sample.name}.wham.germline.vcf.gz")
            delly_vcf = os.path.join(pair_directory, f"{sample.name}.delly.germline.vcf.gz")
            cnvkit_vcf = os.path.join(pair_directory, f"{sample.name}.cnvkit.germline.vcf.gz")
            breakseq_vcf = os.path.join(pair_directory, f"{sample.name}.breakseq.germline.vcf.gz")

            if os.path.isfile(isize_file):
                isize_mean, isize_sd = metric_tools.extract_isize(isize_file)

            else:
                isize_mean = 325
                isize_sd = 50

            input_cnvkit = None
            if os.path.isfile(cnvkit_vcf):
                input_cnvkit = cnvkit_vcf

            input_gatk = None
            if os.path.isfile(gatk_vcf):
                input_gatk = gatk_pass

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(ensemble_directory, remove=True),
                            bcftools.view(
                                gatk_vcf,
                                gatk_pass,
                                global_conf.global_get('metasv_ensemble', 'filter_pass_options')
                            ),
                            htslib.tabix(
                                gatk_pass,
                                options="-pvcf"
                            )
                        ],
                        name=f"metasv_ensemble.gatk_pass.{sample.name}",
                        samples=[sample]
                    )
                )

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        metasv.ensemble(
                            lumpy_vcf,
                            manta_vcf,
                            input_cnvkit,
                            wham_vcf,
                            delly_vcf,
                            input_gatk,
                            input_tumor,
                            sample.name,
                            os.path.join(ensemble_directory, "rawMetaSV"),
                            ensemble_directory,
                            isize_mean=str(isize_mean),
                            isize_sd=str(isize_sd),
                            output_vcf=os.path.join(ensemble_directory, "variants.vcf.gz"),
                            breakseq=breakseq_vcf
                        )
                    ],
                    name=f"metasv_ensemble.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def metasv_sv_annotation(self):
        """
        Annotate VCF with [SnpEff](https://pcingola.github.io/SnpEff/snpeff/introduction/).
        SnpEff is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes).
        Returns:
            list: A list of metasv sv annotation jobs.
        """

        jobs = []

        for sample in self.samples:
            ensemble_directory = os.path.join(self.output_dirs["sv_variants_directory"], "ensemble", sample.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            os.path.join(ensemble_directory, "variants.vcf.gz"),
                            os.path.join(ensemble_directory, f"{sample.name}.metasv.snpeff.vcf")
                        )
                    ],
                    name=f"sv_annotation.metasv_ensemble.{sample.name}",
                    samples=[sample]
                )
            )

        return jobs

    def set_interval_list(self):
        """
        Create an interval list with ScatterIntervalsByNs from GATK: [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360041416072-ScatterIntervalsByNs-Picard).
        Used for creating a broken-up interval list that can be used for scattering a variant-calling pipeline in a way that will not cause problems at the edges of the intervals. 
        By using large enough N blocks (so that the tools will not be able to anchor on both sides) we can be assured that the results of scattering and gathering 
        the variants with the resulting interval list will be the same as calling with one large region.
        Returns:
            list: A list of set interval list jobs.
        """
        jobs = []

        reference = global_conf.global_get('gatk_scatterIntervalsByNs', 'genome_fasta', param_type='filepath')
        dictionary = global_conf.global_get('gatk_scatterIntervalsByNs', 'genome_dictionary', param_type='filepath')
        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        if 'tumor_only' in self.protocol or 'germline' in self.protocol:
            for sample in self.samples:
                interval_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name, "intervals")
                output = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.interval_list'))
                interval_list_acgt_noalt = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))

                coverage_bed = bvatools.resolve_readset_coverage_bed(
                    sample.readsets[0]
                )
                if coverage_bed:
                    dictionary = global_conf.global_get('gatk_scatterIntervalsByNs', 'genome_dictionary', param_type='filepath')
                    region = coverage_bed
                    interval_list = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.interval_list'))
                    interval_list_noalt = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.noALT.interval_list'))

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.bed2interval_list(
                                    dictionary,
                                    region,
                                    interval_list
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            interval_list,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_noalt,
                                            '-v "EBV"'
                                        )
                                    ]
                                ),
                            ],
                            name=f"gatk_scatterIntervalsByNs.{sample.name}",
                            samples=[sample],
                            readsets=[*list(sample.readsets)]
                        )
                    )
                elif scatter_jobs == 1:
                    bed_file = os.path.join(interval_directory, os.path.basename(reference).replace('.ACGT.noALT.interval_list', '.bed'))
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.scatterIntervalsByNs(
                                    reference,
                                    output
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            output,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_acgt_noalt,
                                            '-v "EBV"'
                                        )
                                    ]
                                ),
                                gatk4.interval_list2bed(
                                    interval_list_acgt_noalt,
                                    bed_file
                                ),
                            ],
                            name=f"gatk_scatterIntervalsByNs.{sample.name}",
                            samples=[sample],
                            readsets=[*list(sample.readsets)]
                        )
                    )
                else:
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.scatterIntervalsByNs(
                                    reference,
                                    output
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            output,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_acgt_noalt,
                                            '-v "EBV"'
                                        )
                                    ]
                                ),
                                gatk4.splitInterval(
                                    interval_list_acgt_noalt,
                                    interval_directory,
                                    global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
                                )
                            ],
                            name=f"gatk_scatterIntervalsByNs.{sample.name}",
                            samples=[sample],
                            readsets=[*list(sample.readsets)]
                        )
                    )

        if 'somatic' in self.protocol and 'tumor_only' not in self.protocol:
            for tumor_pair in self.tumor_pairs.values():
                interval_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "intervals")
                output = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.interval_list'))

                coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

                interval_list_acgt_noalt = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))
                bed_acgt_noalt = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.bed'))

                if self.protocol == "fastpass":
                    coverage_bed = global_conf.global_get('rawmpileup_panel', 'panel')

                if coverage_bed:
                    region = coverage_bed
                    interval_list = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.interval_list'))
                    interval_list_noalt = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.noALT.interval_list'))
                    bed_noalt = os.path.join(interval_directory, os.path.basename(region).replace('.bed', '.noALT.bed'))
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.bed2interval_list(
                                    dictionary,
                                    region,
                                    interval_list
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            interval_list,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_noalt,
                                            '-v "EBV"'
                                        ),
                                    ]
                                ),
                                gatk4.interval_list2bed(
                                    interval_list_noalt,
                                    bed_noalt
                                ),
                            ],
                            name=f"gatk_scatterIntervalsByNs.{tumor_pair.name}",
                            samples=[tumor_pair.tumor, tumor_pair.normal],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
                elif scatter_jobs == 1:
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.scatterIntervalsByNs(
                                    reference,
                                    output
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            output,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_acgt_noalt,
                                            '-v "EBV"'
                                        )
                                    ]
                                ),
                                gatk4.interval_list2bed(
                                    interval_list_acgt_noalt,
                                    bed_acgt_noalt
                                ),
                            ],
                            name=f"gatk_scatterIntervalsByNs.{tumor_pair.name}",
                            samples=[tumor_pair.tumor, tumor_pair.normal],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
                else:
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(interval_directory),
                                gatk4.scatterIntervalsByNs(
                                    reference,
                                    output
                                ),
                                pipe_jobs(
                                    [
                                        bash.grep(
                                            output,
                                            None,
                                            '-Ev "_GL|_K"'
                                        ),
                                        bash.grep(
                                            None,
                                            interval_list_acgt_noalt,
                                            '-v "EBV"'
                                        )
                                    ]
                                ),
                                gatk4.interval_list2bed(
                                    interval_list_acgt_noalt,
                                    bed_acgt_noalt
                                ),
                                gatk4.splitInterval(
                                    interval_list_acgt_noalt,
                                    interval_directory,
                                    global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
                                )
                            ],
                            name=f"gatk_scatterIntervalsByNs.{tumor_pair.name}",
                            samples=[tumor_pair.tumor, tumor_pair.normal],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def sequenza(self):
        """
        Sequenza is a novel set of tools providing a fast Python script to genotype cancer samples,
        and an R package to estimate cancer cellularity, ploidy, genome-wide copy number profile and infer
        for mutated alleles.
        Returns:
            list: A list of sequenza jobs.
        """
        jobs = []
        nb_jobs = global_conf.global_get('sequenza', 'nb_jobs', param_type='posint')
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            sequenza_directory = os.path.join(pair_directory, "sequenza")
            raw_sequenza_directory = os.path.join(sequenza_directory, "rawSequenza")

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            raw_output = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name)
            output = os.path.join(sequenza_directory, tumor_pair.name)

            if nb_jobs == 1:
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                raw_sequenza_directory,
                                remove=False
                            ),
                            sequenza.bam2seqz(
                                input_normal,
                                input_tumor,
                                global_conf.global_get('sequenza', 'gc_file'),
                                f"{raw_output}.all.seqz.gz",
                                None
                            ),
                            sequenza.bin(
                                f"{raw_output}.all.seqz.gz",
                                f"{output}.all.binned.seqz.gz",
                            )
                        ],
                        name=f"sequenza.create_seqz.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[raw_sequenza_directory]
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                raw_sequenza_directory,
                                remove=True
                            ),
                            sequenza.main(
                                f"{output}.all.binned.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name
                            )
                        ],
                        name=f"sequenza.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[raw_sequenza_directory]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        raw_sequenza_directory,
                                        remove=False
                                    ),
                                    sequenza.bam2seqz(
                                        input_normal,
                                        input_tumor,
                                        global_conf.global_get('sequenza', 'gc_file'),
                                        f"{raw_output}.seqz.{sequence['name']}.gz",
                                        sequence['name']
                                    ),
                                    sequenza.bin(
                                        f"{raw_output}.seqz.{sequence['name']}.gz",
                                        f"{raw_output}.binned.seqz.{sequence['name']}.gz",
                                    )
                                ],
                                name=f"sequenza.create_seqz.{sequence['name']}.{tumor_pair.name}",
                                samples=[tumor_pair.normal, tumor_pair.tumor],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                                removable_files=[raw_sequenza_directory]
                            )
                        )

                seqz_outputs = [f"{raw_output}.binned.seqz.{sequence['name']}.gz" for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary']

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                raw_sequenza_directory,
                                remove=True
                            ),
                            pipe_jobs(
                                [
                                    bash.cat(
                                        seqz_outputs,
                                        None,
                                        zip=True
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "'FNR==1 && NR==1{print;}{ if($1!=\"chromosome\" && $1!=\"MT\" && $1!=\"chrMT\" && $1!=\"chrM\") {print $0} }'"
                                    ),
                                    bash.gzip(
                                        None,
                                        f"{output}.binned.merged.seqz.gz",
                                        options="-cf"
                                    )
                                ]
                            )
                        ],
                        name=f"sequenza.merge_binned_seqz.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                raw_sequenza_directory,
                                remove=True
                            ),
                            sequenza.main(
                                f"{output}.binned.merged.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name,
                                ini_section='sequenza_estimate'
                            )
                        ],
                        name=f"sequenza_estimate.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def sym_link_sequenza(self):
        """
        Sym link of sequenza outputs.
        Returns:
            list: A list of sym link sequenza jobs.
        """
        jobs = []

        inputs = {}

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            inputs["Tumor"] = [
                os.path.join(pair_directory, "sequenza", f"{tumor_pair.name}_chromosome_view.pdf"),
                os.path.join(pair_directory, "sequenza", f"{tumor_pair.name}_genome_view.pdf"),
                os.path.join(pair_directory, "sequenza", f"{tumor_pair.name}_CN_bars.pdf"),
                os.path.join(pair_directory, "sequenza", f"{tumor_pair.name}_CP_contours.pdf"),
                os.path.join(pair_directory, "sequenza", f"{tumor_pair.name}_ploidy_celularity.tsv")
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv/cnv",
                                    sample=key,
                                    profyle=self.profyle
                                )
                            ],
                            name=f"sym_link.sequenza.{str(idx)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def preprocess_vcf(self):
        """
        Preprocess vcf for loading into an annotation database - Gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into Gemini.
        Returns:
            list: A list of preprocess vcf jobs.
        """

        jobs = []

        if 'fastpass' in self.protocol:
            for tumor_pair in self.tumor_pairs.values():
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
                prefix = os.path.join(pair_directory, tumor_pair.name)

                jobs.append(
                    pipe_jobs(
                        [
                            tools.preprocess_varscan(
                                f"{prefix}.varscan2.somatic.vt.vcf.gz",
                                None,
                            ),
                            htslib.bgzip_tabix(
                                None,
                                f"{prefix}.varscan2.somatic.vt.prep.vcf.gz"
                            )
                        ],
                        name=f"preprocess_vcf.panel.somatic.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

                jobs.append(
                    pipe_jobs(
                        [
                            tools.preprocess_varscan(
                                f"{prefix}.varscan2.germline.vt.vcf.gz",
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                f"{prefix}.varscan2.germline.vt.prep.vcf.gz"
                            )
                        ],
                        name=f"preprocess_vcf.panel.germline.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        if 'high_cov' in self.protocol:
            prefix = os.path.join(self.output_dirs['variants_directory'], "allSamples")
            output_preprocess = f"{prefix}.tmp.vcf.gz"
            output_fix = f"{prefix}.fix.vcf.gz"

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                tools.preprocess_varscan(
                                    f"{prefix}.vcf.gz",
                                    None,
                                ),
                                htslib.bgzip(
                                    None,
                                    output_preprocess
                                ),
                            ]
                        ),
                        Job(
                            [output_preprocess],
                            [output_fix],
                            command=f"""\
zgrep -v 'ID=AD_O' {output_preprocess} | awk 'BEGIN {{OFS=\"\\t\"; FS=\"\\t\"}} \
{{
    if (NF > 8)
    {{
        for (i=9;i<=NF;i++)
        {{
            x=split($i,na,\":\");
            if (x > 1)
            {{
                tmp=na[1];
                for (j=2;j<x;j++)
                {{
                    if (na[j] == \"AD_O\")
                    {{
                        na[j]=\"AD\"
                    }};
                    if (na[j] != \".\")
                    {{
                        tmp=tmp\":\"na[j]}}
                    }};
                    $i=tmp
                }}
            }}
        }};
    print $0
}}' | bgzip -cf > {output_fix}"""
                        ),
                        pipe_jobs(
                            [
                                tools.preprocess_varscan(
                                    output_fix,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{prefix}.prep.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    f"{prefix}.prep.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{prefix}.vt.vcf.gz"
                                )
                            ]
                        ),
                    ],
                    name="preprocess_vcf.germline.allSamples",
                    samples=self.samples,
                    readsets=self.readsets,
                    input_dependency=[f"{prefix}.vcf.gz"],
                    output_dependency=[output_fix, output_preprocess, f"{prefix}.vt.vcf.gz"],
                    removable_files=[output_fix, output_preprocess]
                )
            )

        return jobs

    def sym_link_panel(self):
        """
        Create sym links of panel variants for deliverables to the clients.
        Returns:
            list: A list of sym link panel jobs.
        """
        jobs = []

        assembly = global_conf.global_get('report_pcgr_fastpass', 'assembly')

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] = [os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    
                    # pcgr output file pattern differs by version
                    if global_conf.global_get('pcgr', 'module_pcgr').split("/")[2] < "2":
                        pcgr_output = os.path.join(sample_prefix, "pcgr", f"{tumor_pair.name}.pcgr_acmg.{assembly}.flexdb.html")
                    else:
                        pcgr_output = os.path.join(sample_prefix, "pcgr", f"{tumor_pair.name}.pcgr.{assembly}.html")                    
                    
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.vcf.gz"),
                                    tumor_pair, self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.vcf.gz.tbi"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.somatic.annot.flt.vcf.gz"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.somatic.annot.flt.vcf.gz.tbi"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.germline.annot.flt.vcf.gz"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, f"{tumor_pair.name}.varscan2.germline.annot.flt.vcf.gz.tbi"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    os.path.join(sample_prefix, "cpsr", f"{tumor_pair.name}.cpsr.{assembly}.html"),
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    pcgr_output,
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.profyle
                                )
                            ],
                            name=f"sym_link_panel.{str(idx)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def manta_sv_calls(self):
        """
        Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. It is optimized for
        the analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
        Manta discovers, assembles and scores large-scale SVs, medium-sized indels and large insertions within a
        single efficient workflow.
        Outputs: Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
        in VCF 4.1 format.
        Returns:
            list: A list of manta sv calls jobs.
        """
        jobs = []

        reference = global_conf.global_get('manta_sv_calls', 'genome_fasta', param_type='filepath')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            interval_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "intervals")

            manta_directory = os.path.join(pair_directory, "rawManta")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )
            manta_somatic_output = os.path.join(manta_directory, "results", "variants", "somaticSV.vcf.gz")
            manta_germline_output = os.path.join(manta_directory, "results", "variants", "diploidSV.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if coverage_bed:
                interval_bed = coverage_bed

            else:
                interval_bed = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.bed'))

            if interval_bed:
                local_coverage_bed = os.path.join(manta_directory, os.path.basename(interval_bed))
                bed_file = f"{local_coverage_bed}.gz"

            output_dep = [
                manta_somatic_output,
                f"{manta_somatic_output}.tbi",
                manta_germline_output,
                f"{manta_germline_output}.tbi"
            ]

            sed_cmd = Job(
                [os.path.join(manta_directory, "runWorkflow.py")],
                [os.path.join(manta_directory, "runWorkflow.py")],
                command=f"""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {os.path.join(manta_directory, "runWorkflow.py")}"""
            )
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(manta_directory),
                        bash.mkdir(
                            manta_directory,
                            remove=True
                        ),
                        bash.sort(
                            interval_bed,
                            f"{local_coverage_bed}.sort",
                            "-k1,1V -k2,2n -k3,3n",
                            extra="; sleep 15"
                        ),
                        htslib.bgzip(
                            f"{local_coverage_bed}.sort",
                            bed_file
                        ),
                        htslib.tabix(
                            bed_file,
                            "-f -p bed"
                        ),
                        manta.manta_config(
                            input_normal,
                            input_tumor,
                            manta_directory,
                            bed_file
                        ),
                        sed_cmd,
                        manta.manta_run(
                            manta_directory,
                            output_dep=output_dep
                        ),
                        bash.ln(
                            os.path.relpath(manta_somatic_output, os.path.dirname(output_prefix)),
                            f"{output_prefix}.manta.somatic.vcf.gz",
                            input=manta_somatic_output,
                        ),
                        bash.ln(
                            os.path.relpath(manta_somatic_output + ".tbi", os.path.dirname(output_prefix)),
                            f"{output_prefix}.manta.somatic.vcf.gz.tbi",
                            input=f"{manta_somatic_output}.tbi"
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output, os.path.dirname(output_prefix)),
                            f"{output_prefix}.manta.germline.vcf.gz",
                            input=manta_germline_output
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output + ".tbi", os.path.dirname(output_prefix)),
                            f"{output_prefix}.manta.germline.vcf.gz.tbi",
                            input=f"{manta_germline_output}.tbi"
                        )
                    ],
                    name=f"manta_sv.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor, interval_bed]
                )
            )

        return jobs

    def strelka2_paired_somatic(self):
        """
        [Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs
        This implementation is optimized for somatic calling.
        Returns:
            list: A list of strelka2 paired somatic jobs.
        """
        jobs = []

        reference = global_conf.global_get('strelka2_paired_somatic', 'genome_fasta', param_type='filepath')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            interval_directory = os.path.join(pair_directory, "intervals")

            somatic_dir = os.path.join(pair_directory, "rawStrelka2_somatic")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            manta_indels = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if coverage_bed:
                interval_bed = coverage_bed

            else:
                interval_bed = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.bed'))

            if interval_bed:
                local_coverage_bed = os.path.join(somatic_dir, os.path.basename(interval_bed))
                bed_file = f"{local_coverage_bed}.gz"

            output_dep = [
                os.path.join(somatic_dir, "results", "variants", "somatic.snvs.vcf.gz"),
                os.path.join(somatic_dir, "results", "variants", "somatic.indels.vcf.gz")
            ]

            sed_cmd = Job(
                [os.path.join(somatic_dir, "runWorkflow.py")],
                [os.path.join(somatic_dir, "runWorkflow.py")],
                command=f"""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {os.path.join(somatic_dir, "runWorkflow.py")}"""
            )

            jobs.append(
                concat_jobs(
                    [
                        bash.rm(somatic_dir),
                        bash.mkdir(somatic_dir),
                        bash.sort(
                            interval_bed,
                            local_coverage_bed + ".sort",
                            "-V -k1,1 -k2,2n -k3,3n",
                            extra="; sleep 15"
                        ),
                        htslib.bgzip(
                            local_coverage_bed + ".sort",
                            bed_file
                        ),
                        htslib.tabix(
                            bed_file,
                            "-f -p bed"
                        ),
                        strelka2.somatic_config(
                            input_normal,
                            input_tumor,
                            somatic_dir,
                            bed_file,
                            manta_indels
                        ),
                        sed_cmd,
                        strelka2.run(
                            somatic_dir,
                            output_dep=output_dep,
                            ini_section='strelka2_paired_somatic'
                        )
                    ],
                    name=f"strelka2_paired_somatic.call.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor, manta_indels, interval_bed],
                    output_dependency=output_dep
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                bcftools.concat(
                                    output_dep,
                                    None
                                ),
                                pipe_jobs(
                                    [
                                        bash.sed(
                                            None,
                                            None,
                                            "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                        ),
                                        bash.sed(
                                            None,
                                            None,
                                            "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                        ),
                                        bash.sed(
                                            None,
                                            None,
                                            "'s/Number=R/Number=./g'"
                                        ),
                                    ]
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{output_prefix}.strelka2.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    f"{output_prefix}.strelka2.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{output_prefix}.strelka2.vt.vcf.gz"
                                )
                            ]
                        ),
                        tools.fix_genotypes_strelka(
                            f"{output_prefix}.strelka2.vt.vcf.gz",
                            f"{output_prefix}.strelka2.somatic.gt.vcf.gz",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name
                        ),
                        bcftools.view(
                            f"{output_prefix}.strelka2.somatic.gt.vcf.gz",
                            f"{output_prefix}.strelka2.somatic.vt.vcf.gz",
                            global_conf.global_get('strelka2_paired_somatic', 'filter_options')
                        )
                    ],
                    name=f"strelka2_paired_somatic.filter.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def strelka2_paired_germline(self):
        """
        [Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs.
        This implementation is optimized for germline calling in cancer pairs.
        Returns:
            list: A list of strelka2 paired germline jobs.
        """
        jobs = []

        reference = global_conf.global_get('strelka2_paired_germline', 'genome_fasta', param_type='filepath')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            interval_directory = os.path.join(pair_directory, "intervals")

            germline_dir = os.path.join(pair_directory, "rawStrelka2_germline")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if coverage_bed:
                interval_bed = coverage_bed

            else:
                interval_bed = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.bed'))
            if interval_bed:
                local_coverage_bed = os.path.join(germline_dir, os.path.basename(interval_bed))
                bed_file = local_coverage_bed + ".gz"

            output_dep = [os.path.join(germline_dir, "results", "variants", "variants.vcf.gz")]

            sed_cmd = Job(
                [os.path.join(germline_dir, "runWorkflow.py")],
                [os.path.join(germline_dir, "runWorkflow.py")],
                command=f"""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {os.path.join(germline_dir, "runWorkflow.py")}"""
            )
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(germline_dir),
                        bash.mkdir(germline_dir),
                        bash.sort(
                            interval_bed,
                            f"{local_coverage_bed}.sort",
                            "-V -k1,1 -k2,2n -k3,3n",
                            extra="; sleep 15"
                        ),
                        htslib.bgzip(
                            f"{local_coverage_bed}.sort",
                            bed_file
                        ),
                        htslib.tabix(
                            bed_file,
                            "-f -p bed"
                        ),
                        strelka2.germline_config(
                            [input_normal, input_tumor],
                            germline_dir,
                            bed_file
                        ),
                        sed_cmd,
                        strelka2.run(
                            germline_dir,
                            output_dep=output_dep,
                            ini_section='strelka2_paired_germline'
                        )
                    ],
                    name=f"strelka2_paired_germline.call.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor, interval_bed],
                    output_dependency=output_dep
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                pipe_jobs(
                                    [
                                        bash.cat(
                                            os.path.join(germline_dir, "results", "variants", "variants.vcf.gz"),
                                            None,
                                            zip=True
                                        ),
                                        bash.sed(
                                            None,
                                            None,
                                            "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                        ),
                                        bash.sed(
                                            None,
                                            None,
                                            "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                        ),
                                        bash.sed(
                                            None,
                                            None,
                                            's/Number=R/Number=./g'""
                                        ),
                                    ]
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{output_prefix}.strelka2.germline.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    f"{output_prefix}.strelka2.germline.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    f"{output_prefix}.strelka2.germline.gt.vcf.gz"
                                )
                            ]
                        ),
                        bcftools.view(
                            f"{output_prefix}.strelka2.germline.gt.vcf.gz",
                            f"{output_prefix}.strelka2.germline.vt.vcf.gz",
                            global_conf.global_get('strelka2_paired_germline', 'filter_options')
                        )
                    ],
                    name=f"strelka2_paired_germline.filter.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def strelka2_paired_snpeff(self):
        """
        [Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs.
        This implementation is optimized for germline calling in cancer pairs.
        Returns:
            list: A list of strelka2 paired snpeff jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        bcftools.split(
                            os.path.join(pair_directory, f"{tumor_pair.name}.strelka2.germline.vt.vcf.gz"),
                            pair_directory,
                            global_conf.global_get('strelka2_paired_snpeff', 'split_options'),
                        )
                    ],
                    name="strelka2_paired_germline_snpeff.split." + tumor_pair.name,
                    input_dependency=[
                        os.path.join(pair_directory, f"{tumor_pair.name}.strelka2.germline.vt.vcf.gz")
                    ],
                    output_dependency=[
                        os.path.join(pair_directory, f"{tumor_pair.normal.name}.vcf.gz"),
                        os.path.join(pair_directory, f"{tumor_pair.tumor.name}.vcf.gz")
                    ],
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            os.path.join(pair_directory, f"{tumor_pair.normal.name}.vcf.gz"),
                            os.path.join(pair_directory, f"{tumor_pair.normal.name}.snpeff.vcf"),
                            options=global_conf.global_get('strelka2_paired_snpeff', 'options')
                        ),
                        htslib.bgzip_tabix(
                            os.path.join(pair_directory, f"{tumor_pair.normal.name}.snpeff.vcf"),
                            os.path.join(pair_directory, f"{tumor_pair.normal.name}.snpeff.vcf.gz")
                        )
                    ],
                    name="strelka2_paired_snpeff.normal." + tumor_pair.name,
                    samples=[tumor_pair.normal],
                    readsets=list(tumor_pair.normal.readsets)
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".vcf.gz"),
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".snpeff.vcf"),
                            options=global_conf.global_get('strelka2_paired_snpeff', 'options')
                        ),
                        htslib.bgzip_tabix(
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".snpeff.vcf"),
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".snpeff.vcf.gz")
                        )
                    ],
                    name="strelka2_paired_snpeff.tumor." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
        return jobs

    def purple(self, sv=False):
        """
        PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

        It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
        somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.
        Returns:
            list: A list of purple jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_dir = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            metrics_dir = self.output_dirs['metrics_directory'][tumor_pair.name]

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            annotate_snv = None
            if os.path.join(pair_dir, f"{tumor_pair.name}.strelka2.somatic.vt.vcf.gz"):
                somatic_snv = os.path.join(pair_dir, f"{tumor_pair.name}.strelka2.somatic.purple.vcf.gz")
                annotate_snv = os.path.join(pair_dir, f"{tumor_pair.name}.strelka2.somatic.pave.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                            purple.strelka2_convert(
                                os.path.join(pair_dir, f"{tumor_pair.name}.strelka2.somatic.vt.vcf.gz"),
                                somatic_snv
                            ),
                            pave.run(
                                somatic_snv,
                                tumor_pair.tumor.name,
                                pair_dir,
                                annotate_snv,
                                ini_section='pave_annotate'
                            )
                        ],
                    name=f"purple.annotate_strelka2.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            gripss_vcf = None
            gripss_filtered_vcf = None
            germline_hotspots = None
            somatic_hotspots = global_conf.global_get('purple', 'somatic_hotspots', param_type='filepath')
            driver_gene_panel = global_conf.global_get('purple', 'driver_gene_panel', param_type='filepath')

            if sv:
                pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
                gridss_directory = os.path.join(pair_dir, "gridss")
                gripss_vcf = os.path.join(gridss_directory, f"{tumor_pair.tumor.name}.gripss.somatic.vcf.gz")
                gripss_filtered_vcf = os.path.join(gridss_directory, f"{tumor_pair.tumor.name}.gripss.filtered.somatic.vcf.gz")
                germline_hotspots = global_conf.global_get('purple', 'germline_hotspots', param_type='filepath')

            purple_dir = os.path.join(pair_dir, "purple")
            amber_dir = os.path.join(purple_dir, "rawAmber")
            cobalt_dir = os.path.join(purple_dir, "rawCobalt")
            ensembl_data_dir = global_conf.global_get('purple', 'ensembl_data_dir', param_type='dirpath')

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            amber_dir,
                            remove=True
                        ),
                        amber.run(
                            input_normal,
                            input_tumor,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            amber_dir,
                        )
                    ],
                    name="purple.amber." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            cobalt_dir,
                            remove=True
                        ),
                        cobalt.run(
                            input_normal,
                            input_tumor,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            cobalt_dir,
                        )
                    ],
                    name="purple.cobalt." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
            purple_purity_output = os.path.join(purple_dir, f"{tumor_pair.tumor.name}.purple.purity.tsv")
            purple_qc_output = os.path.join(purple_dir, f"{tumor_pair.tumor.name}.purple.qc")
            samples = [tumor_pair.normal, tumor_pair.tumor]
            job_name = f"purple.purity.{tumor_pair.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                    purple.parse_purity_metrics_pt(purple_purity_output),
                    job2json_project_tracking.run(
                        input_file=purple_purity_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in samples]),
                        readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="purity=$purity"
                        )
                    ])

            jobs.append(
                concat_jobs(
                    [
                        purple.run(
                            amber_dir,
                            cobalt_dir,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            purple_dir,
                            ensembl_data_dir,
                            annotate_snv,
                            gripss_vcf,
                            gripss_filtered_vcf,
                            somatic_hotspots,
                            germline_hotspots,
                            driver_gene_panel
                        ),
                        bash.mkdir(
                            metrics_dir
                        ),
                        bash.ln(
                            os.path.relpath(purple_purity_output.strip(), metrics_dir),
                            os.path.join(metrics_dir, os.path.basename(purple_purity_output)),
                            input=purple_purity_output
                        ),
                        bash.ln(
                            os.path.relpath(purple_qc_output.strip(), metrics_dir),
                            os.path.join(metrics_dir, os.path.basename(purple_qc_output)),
                            input=purple_qc_output
                        ),
                        job_project_tracking_metrics
                    ],
                    name=job_name,
                    samples=samples,
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
            self.multiqc_inputs[tumor_pair.name].extend(
                [
                    os.path.join(self.output_dirs['metrics_directory'][tumor_pair.name], os.path.basename(purple_purity_output)),
                    os.path.join(self.output_dirs['metrics_directory'][tumor_pair.name], os.path.basename(purple_qc_output))
                ]
            )

        return jobs
    def germline_varscan2(self):
        """
        [VarScan](https://dkoboldt.github.io/varscan/) caller for insertions and deletions.
        Returns:
            list: A list of germline varscan2 jobs.
        """

        jobs = []

        scatter_jobs = global_conf.global_get('germline_varscan2', 'nb_jobs', param_type='posint')
        if scatter_jobs > 50:
            log.warning("Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        compression_postfix = global_conf.global_get('bwa_mem2_samtools_sort', 'compression')
        variants_directory = os.path.join(self.output_dirs["variants_directory"])
        varscan_directory = os.path.join(variants_directory, "rawVarScan")
        output = os.path.join(variants_directory, "allSamples.vcf.gz")

        samples_file = 'sample_list.tsv'
        sample_list = open(samples_file, 'w')
        bam_list = [os.path.join(self.output_dirs["alignment_directory"], sample.name, f"{sample.name}.sorted.fixmate.{compression_postfix}") for sample in self.samples]
        bed_file = ""
        for sample in self.samples:
            sample_list.write(f"{sample.name}\n")
            bed_file = bvatools.resolve_readset_coverage_bed(sample.readsets[0])

        if bed_file or scatter_jobs == 1:
            jobs.append(
                concat_jobs(
                [
                    bash.mkdir(varscan_directory),
                    pipe_jobs(
                        [
                            samtools.mpileup(
                                bam_list,
                                None,
                                region=global_conf.global_get('germline_varscan2', 'regions'),
                                regionFile=bed_file,
                                ini_section='germline_varscan2'
                            ),
                            varscan.mpileupcns(
                                None,
                                None,
                                samples_file
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output
                            )
                        ]
                    ),
                ],
                name= "germline_varscan2.single",
                samples=self.samples,
                readsets=self.readsets
                )
            )
        else:
            input_vcfs = []
            for sequence in self.sequence_dictionary_variant():
                if sequence['type'] == 'primary':
                    output = os.path.join(varscan_directory,
                                               f"allSamples.{sequence['name']}.vcf.gz")
                    input_vcfs.append(output)

                    jobs.append(
                        pipe_jobs(
                        [
                            samtools.mpileup(
                                bam_list,
                                None,
                                region=sequence['name'],
                                regionFile=None,
                                ini_section='germline_varscan2'
                            ),
                            varscan.mpileupcns(
                                None,
                                None,
                                sample_list,
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output
                            )
                        ],
                            name = f"germline_varscan2.{sequence['name']}",
                            samples= self.samples,
                            readsets=self.readsets,
                            input_dependency=bam_list,
                            output_dependency=[output]
                        )
                    )
            output = os.path.join(variants_directory, "allSamples.vcf.gz")
            jobs.append(
                concat_jobs(
                    [
                        gatk4.cat_variants(
                            input_vcfs,
                            output
                        )
                    ],
                    name = "gatk_cat_germline_varscan2",
                    samples = self.samples,
                    readsets = self.readsets,
                    input_dependency = input_vcfs,
                    output_dependency = [output]
                )
            )

        return jobs

    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        Returns:
            list: A list of rawmpileup jobs.
        """

        reference = global_conf.global_get('rawpileup', 'genome_fasta', param_type='filepath')

        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
        if scatter_jobs > 50:
            log.warning("Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            if 'fastpass' in self.protocol:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, 'panel')
                ini_section = 'rawmpileup_panel'
                scatter_jobs = global_conf.global_get(ini_section, 'nb_jobs', param_type='posint')
                bed_file = global_conf.global_get(ini_section, 'panel')
                job_name = f"rawmpileup_panel.{tumor_pair.name}"

            else:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

                interval_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "intervals")
                bed_file = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.bed'))

                ini_section='rawmpileup'
                job_name = f"rawmpileup.{tumor_pair.name}"

                coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

                if coverage_bed:
                    bed_file = coverage_bed

            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if scatter_jobs == 1:
                pair_output = os.path.join(varscan_directory, f"{tumor_pair.name}.mpileup")
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            samtools.mpileup(
                                [
                                    input_normal,
                                    input_tumor
                                ],
                                pair_output,
                                regionFile = bed_file,
                                ini_section = ini_section
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        input_dependency=[input_normal, input_tumor, bed_file]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        pair_output = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.mpileup")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=True
                                    ),
                                    samtools.mpileup(
                                        [
                                            input_normal,
                                            input_tumor
                                        ],
                                        pair_output,
                                        regionFile=bed_file,
                                        region=sequence['name'],
                                        ini_section=ini_section
                                    )
                                ],
                                name=f"{job_name}.{sequence['name']}",
                                samples=[tumor_pair.normal, tumor_pair.tumor],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                                input_dependency=[input_normal, input_tumor, bed_file]
                            )
                        )

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by the author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collision with Samtools output during ensemble.
        Returns:
            list: A list of paired varscan2 jobs.
        """

        jobs = []

        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
        if scatter_jobs > 50:
            log.warning("Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if 'fastpass' in self.protocol:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, 'panel')
                ini_section = 'varscan2_somatic_panel'
                job_name = f"varscan2_somatic_panel.{tumor_pair.name}"

            else:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
                ini_section = 'varscan2_somatic'
                job_name = f"varscan2_somatic.{tumor_pair.name}"

            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
            output = os.path.join(varscan_directory, tumor_pair.name)

            if scatter_jobs == 1:
                input_pair = os.path.join(varscan_directory, f"{tumor_pair.name}.mpileup")
                output_snp = os.path.join(varscan_directory, f"{tumor_pair.name}.snp.vcf")
                output_indel = os.path.join(varscan_directory, f"{tumor_pair.name}.indel.vcf")
                output_vcf = os.path.join(varscan_directory, f"{tumor_pair.name}.varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, f"{tumor_pair.name}.varscan2.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=False
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                output_vcf_dep=output_vcf,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel,
                                ini_section=ini_section
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, f"{tumor_pair.name}.snp.vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, f"{tumor_pair.name}.indel.vcf.gz")
                            ),
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        [
                                            os.path.join(varscan_directory, f"{tumor_pair.name}.snp.vcf.gz"),
                                            os.path.join(varscan_directory, f"{tumor_pair.name}.indel.vcf.gz")
                                        ],
                                        None
                                    ),
                                    pipe_jobs(
                                        [
                                            bash.sed(
                                                None,
                                                None,
                                                "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                            ),
                                            bash.sed(
                                                None,
                                                None,
                                                "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                            ),
                                            bash.grep(
                                                None,
                                                None,
                                                "-v \"INFO=<ID=SSC\""
                                            ),
                                            bash.sed(
                                                None,
                                                output_vcf,
                                                "-E \"s/SSC=(.*);//g\""
                                            )
                                        ]
                                    )
                                ]
                            ),
                            htslib.bgzip_tabix(
                                output_vcf,
                                output_vcf_gz
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.tumor, tumor_pair.normal],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[varscan_directory]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        input_pair = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.mpileup")
                        output = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}")
                        output_snp = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.snp.vcf")
                        output_indel = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.indel.vcf")
                        output_vcf = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.varscan2.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.varscan2.vcf.gz")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=False
                                    ),
                                    varscan.somatic(
                                        input_pair,
                                        output,
                                        output_vcf_dep=output_vcf,
                                        output_snp_dep=output_snp,
                                        output_indel_dep=output_indel,
                                        ini_section=ini_section
                                    ),
                                    htslib.bgzip_tabix(
                                        output_snp,
                                        os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.snp.vcf.gz")
                                    ),
                                    htslib.bgzip_tabix(
                                        output_indel,
                                        os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.indel.vcf.gz")
                                    ),
                                    pipe_jobs(
                                        [
                                            bcftools.concat(
                                                [
                                                    os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.snp.vcf.gz"),
                                                    os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.indel.vcf.gz")
                                                ],
                                                None
                                            ),
                                            pipe_jobs(
                                                [
                                                    bash.sed(
                                                        None,
                                                        None,
                                                        "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                                    ),
                                                    bash.sed(
                                                        None,
                                                        None,
                                                        "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                                    ),
                                                    bash.grep(
                                                        None,
                                                        None,
                                                        "-v \"INFO=<ID=SSC\""
                                                    ),
                                                    bash.sed(
                                                        None,
                                                        output_vcf,
                                                        "-E \"s/SSC=(.*);//g\""
                                                    )
                                                ]
                                            )
                                        ]
                                    ),
                                    htslib.bgzip_tabix(
                                        output_vcf,
                                        output_vcf_gz
                                    )
                                ],
                                name=f"{job_name}.{sequence['name']}",
                                samples=[tumor_pair.tumor, tumor_pair.normal],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                                removable_files=[varscan_directory]
                            )
                        )
        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        Returns:
            list: A list of merge varscan2 jobs.
        """

        jobs = []

        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
        if scatter_jobs > 50:
            log.warning("Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if 'fastpass' in self.protocol:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, 'panel')
                scatter_jobs = global_conf.global_get('rawmpileup_panel', 'nb_jobs', param_type='posint')
                job_name = f"merge_varscan2_panel.{tumor_pair.name}"

            else:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
                job_name = f"merge_varscan2.{tumor_pair.name}"

            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            if scatter_jobs == 1:
                all_inputs = [os.path.join(varscan_directory, f"{tumor_pair.name}.varscan2.vcf.gz")]

            else:
                all_inputs = [
                    os.path.join(varscan_directory, f"{tumor_pair.name}.{sequence['name']}.varscan2.vcf.gz")
                    for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary'
                ]

            for input_vcf in all_inputs:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete varscan2 vcf: {input_vcf}\n")

            all_output = os.path.join(pair_directory, f"{tumor_pair.name}.varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, f"{tumor_pair.name}.varscan2.vt.vcf.gz")

            somtic_output_vt = os.path.join(pair_directory, f"{tumor_pair.name}.varscan2.somatic.vt.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, f"{tumor_pair.name}.varscan2.germline.vt.vcf.gz")

            if scatter_jobs == 1:
                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        all_inputs[0],
                                        None
                                    ),
                                    tools.fix_varscan_output(
                                        None,
                                        None
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        all_output
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    vt.decompose_and_normalize_mnps(
                                        all_output,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        all_output_vt
                                    )
                                ]
                            ),
                            bcftools.view(
                                all_output_vt,
                                somtic_output_vt,
                                global_conf.global_get('merge_varscan2', 'somatic_filter_options')
                            ),
                            htslib.tabix(
                                somtic_output_vt,
                                global_conf.global_get('merge_varscan2', 'tabix_options')
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        all_output_vt,
                                        None,
                                        global_conf.global_get('merge_varscan2', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        global_conf.global_get('merge_varscan2', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        germline_output_vt
                                    )
                                ]
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[all_output, f"{all_output}.tbi"]
                    )
                )

            else:
                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        all_inputs,
                                        None
                                    ),
                                    tools.fix_varscan_output(
                                        None,
                                        None
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        all_output
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    vt.decompose_and_normalize_mnps(
                                        all_output,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        all_output_vt
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        all_output_vt,
                                        None,
                                        global_conf.global_get('merge_varscan2', 'somatic_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        somtic_output_vt
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        all_output_vt,
                                        None,
                                        global_conf.global_get('merge_varscan2', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        global_conf.global_get('merge_varscan2', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        germline_output_vt
                                    )
                                ]
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[all_output, f"{all_output}.tbi", all_output_vt, f"{all_output_vt}.tbi"]
                    )
                )
        return jobs

    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        Returns:
            list: A list of paired mutect2 jobs.
        """

        jobs = []

        reference = global_conf.global_get('paired_mutect2', 'genome_fasta', param_type='filepath')
        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        if scatter_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            interval_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "intervals")

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])

            if coverage_bed:
                interval_list = [os.path.join(interval_directory, os.path.basename(coverage_bed).replace('.bed', '.noALT.interval_list'))]

            else:
                interval_list = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list')),

            if scatter_jobs == 1 and interval_list is not None:
                jobs.append(
                    concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                mutect_directory,
                                remove=True
                            ),
                            gatk4.mutect2(
                                input_normal,
                                tumor_pair.normal.name,
                                input_tumor,
                                tumor_pair.tumor.name,
                                os.path.join(mutect_directory, f"{tumor_pair.name}.mutect2.vcf.gz"),
                                os.path.join(mutect_directory, f"{tumor_pair.name}.f1r2.tar.gz"),
                                interval_list=interval_list[0]
                            )
                        ],
                        name=f"gatk_mutect2.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        input_dependency=[input_normal, input_tumor, interval_list[0]]
                    )
                )

            elif scatter_jobs > 1:
                interval_list = [os.path.join(interval_directory, f"{idx:04d}-scattered.interval_list") for idx in range(scatter_jobs)]

                # Create one separate job for each of the first sequences
                for idx, intervals in enumerate(interval_list):
                    jobs.append(
                        concat_jobs(
                            [
                                # Create output directory since it is not done by default by GATK tools
                                bash.mkdir(
                                    mutect_directory,
                                    remove=True
                                ),
                                gatk4.mutect2(
                                    input_normal,
                                    tumor_pair.normal.name,
                                    input_tumor,
                                    tumor_pair.tumor.name,
                                    os.path.join(mutect_directory, f"{tumor_pair.name}.{str(idx)}.mutect2.vcf.gz"),
                                    os.path.join(mutect_directory, f"{tumor_pair.name}.{str(idx)}.f1r2.tar.gz"),
                                    interval_list=intervals
                                )
                            ],
                            name=f"gatk_mutect2.{tumor_pair.name}.{str(idx)}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            input_dependency=[input_normal, input_tumor, intervals]
                        )
                    )
        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2.
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants.
        Returns:
            list: A list of merge mutect2 jobs.
        """

        jobs = []

        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            output_prefix = os.path.join(mutect_directory, tumor_pair.name)

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, f"{tumor_pair.name}.mutect2.vcf.gz")
            output_flt = os.path.join(pair_directory, f"{tumor_pair.name}.mutect2.flt.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, f"{tumor_pair.name}.mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, f"{tumor_pair.name}.mutect2.somatic.vt.vcf.gz")

            if scatter_jobs == 1:
                if global_conf.global_get('gatk_mutect2', 'module_gatk').split("/")[2] > "4":
                    jobs.append(
                        concat_jobs(
                            [
                                gatk4.learn_read_orientation_model(
                                    [os.path.join(mutect_directory, f"{tumor_pair.name}.f1r2.tar.gz")],
                                    os.path.join(pair_directory, f"{tumor_pair.name}.f1r2.tar.gz")
                                ),
                                gatk4.filter_mutect_calls(
                                    os.path.join(mutect_directory, f"{tumor_pair.name}.mutect2.vcf.gz"),
                                    output_flt,
                                    read_orientation=os.path.join(pair_directory, f"{tumor_pair.name}.f1r2.tar.gz")
                                ),
                                pipe_jobs(
                                    [
                                        vt.decompose_and_normalize_mnps(
                                            output_flt,
                                            None
                                        ),
                                        pipe_jobs(
                                            [
                                                bash.grep(
                                                    None,
                                                    None,
                                                    "-v 'GL00'"
                                                ),
                                                bash.grep(
                                                    None,
                                                    None,
                                                    "-Ev 'chrUn|random'"
                                                ),
                                                bash.grep(
                                                    None,
                                                    None,
                                                    "-vE 'EBV|hs37d5'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    r"-e 's#/\.##g'"
                                                )
                                            ]
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_vt_gz
                                        )
                                    ]
                                ),
                                pipe_jobs(
                                    [
                                        bcftools.view(
                                            output_vt_gz,
                                            None,
                                            global_conf.global_get('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name=f"merge_filter_mutect2.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            removable_files=[output_flt, f"{output_flt}.tbi"]
                        )
                    )

                else:
                    input_vcf = os.path.join(mutect_directory, f"{tumor_pair.name}.mutect2.vcf.gz")
                    jobs.append(
                        concat_jobs(
                            [
                                bash.ln(
                                    os.path.relpath(input_vcf, os.path.dirname(output_gz)),
                                    output_gz,
                                    input=input_vcf
                                ),
                                pipe_jobs(
                                    [
                                        vt.decompose_and_normalize_mnps(
                                            output_gz,
                                            None
                                        ),
                                        pipe_jobs(
                                            [
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/Number=R/Number=./g'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    r"-e 's#/\.##g'"
                                                )
                                            ]
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name=f"symlink_mutect_vcf.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

            elif scatter_jobs > 1:
                merge_list = [f"{output_prefix}.{idx}.mutect2.vcf.gz" for idx in range(scatter_jobs)]

                for input_vcf in merge_list:
                    if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                        log.error(f"Incomplete mutect2 vcf: {input_vcf}\n")

                if global_conf.global_get('gatk_mutect2', 'module_gatk').split("/")[2] > "4":
                    output_stats = os.path.join(pair_directory, f"{tumor_pair.name}.mutect2.vcf.gz.stats")
                    output_models = os.path.join(pair_directory, f"{tumor_pair.name}.read-orientation-model.tar.gz")

                    stats = []
                    models = []
                    for idx, _ in enumerate(merge_list):
                        stats.append(
                            os.path.join(mutect_directory, f"{tumor_pair.name}.{str(idx)}.mutect2.vcf.gz.stats"))
                        models.append(
                            os.path.join(mutect_directory, f"{tumor_pair.name}.{str(idx)}.f1r2.tar.gz"))

                    jobs.append(
                        concat_jobs(
                            [
                                gatk4.learn_read_orientation_model(
                                    models,
                                    output_models
                                ),
                                gatk4.cat_variants(
                                    merge_list,
                                    output_gz
                                ),
                                gatk4.merge_stats(
                                    stats,
                                    output_stats
                                ),
                                gatk4.filter_mutect_calls(
                                    output_gz,
                                    output_flt,
                                    read_orientation=output_models
                                ),
                                pipe_jobs(
                                    [
                                        vt.decompose_and_normalize_mnps(
                                            output_flt,
                                            None
                                        ),
                                        pipe_jobs(
                                            [
                                                bash.sed(
                                                    None,
                                                    None,
                                                    r"-e 's#/\.##g'"
                                                )
                                            ]
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_vt_gz
                                        )
                                    ]
                                ),
                                pipe_jobs(
                                    [
                                        bcftools.view(
                                            output_vt_gz,
                                            None,
                                            global_conf.global_get('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name=f"merge_filter_mutect2.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                            removable_files=[output_flt, f"{output_flt}.tbi"]
                        )
                    )

                else:
                    jobs.append(
                        concat_jobs(
                            [
                                pipe_jobs(
                                    [
                                        bcftools.concat(
                                            merge_list,
                                            None,
                                            global_conf.global_get('merge_filter_mutect2', 'bcftools_options')
                                        ),
                                        pipe_jobs(
                                            [
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/TUMOR/" + tumor_pair.tumor.name + "/g'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/NORMAL/" + tumor_pair.normal.name + "/g'"
                                                ),
                                                bash.sed(
                                                    None,
                                                    None,
                                                    "'s/Number=R/Number=./g'"
                                                ),
                                            ]
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_gz
                                        )
                                    ]
                                ),
                                # gatk4.filter_mutect_calls(output_gz, output_flt),
                                pipe_jobs(
                                    [
                                        vt.decompose_and_normalize_mnps(
                                            output_gz,
                                            None
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_vt_gz
                                        )
                                    ]
                                ),
                                pipe_jobs(
                                    [
                                        bcftools.view(
                                            output_vt_gz,
                                            None,
                                            global_conf.global_get('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name=f"merge_filter_mutect2.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove the instance where REF == ALT and REF is modified to 'N' when REF is
        AUPAC nomenclature.
        Returns:
            list: A list of vardict paired jobs.
        """
        jobs = []

        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')
        if scatter_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        reference = global_conf.global_get('vardict_paired', 'genome_fasta', param_type='filepath')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            interval_directory = os.path.join(pair_directory, "intervals")

            vardict_directory = os.path.join(pair_directory, "rawVardict")

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
            if coverage_bed:
                coverage_bed = os.path.join(interval_directory, os.path.basename(coverage_bed).replace('.bed', '.noALT.bed'))
                output = os.path.join(vardict_directory, f"{tumor_pair.name}.vardict.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            pipe_jobs(
                                [
                                    vardict.paired_java(
                                        input_normal,
                                        input_tumor,
                                        tumor_pair.name,
                                        None,
                                        coverage_bed
                                    ),
                                    vardict.testsomatic(
                                        None,
                                        None
                                    ),
                                    vardict.var2vcf(
                                        None,
                                        tumor_pair.normal.name,
                                        tumor_pair.tumor.name,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output
                                    )
                                ]
                            )
                        ],
                        name=f"vardict_paired.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        input_dependency=[input_normal, input_tumor, coverage_bed]
                    )
                )
            elif scatter_jobs == 1:
                interval_list = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))
                bed_file = os.path.join(vardict_directory, os.path.basename(interval_list).replace('.interval_list', '.padded.bed'))
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            tools.dict2beds(
                                interval_list,
                                bed_file
                            ),
                        ],
                        name=f"vardict.genome.beds.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
                output = os.path.join(vardict_directory, f"{tumor_pair.name}.vardict.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            pipe_jobs(
                                [
                                    vardict.paired_java(
                                        input_normal,
                                        input_tumor,
                                        tumor_pair.name,
                                        None,
                                        bed_file
                                    ),
                                    vardict.testsomatic(
                                        None,
                                        None
                                    ),
                                    vardict.var2vcf(
                                        None,
                                        tumor_pair.normal.name,
                                        tumor_pair.tumor.name,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output
                                    )
                                ]
                            )
                        ],
                        name=f"vardict_paired.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        input_dependency=[input_normal, input_tumor, bed_file]
                    )
                )

            elif scatter_jobs > 1:
                interval_list = os.path.join(interval_directory, os.path.basename(reference).replace('.fa', '.ACGT.noALT.interval_list'))
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            tools.dict2beds(
                                interval_list,
                                os.path.join(vardict_directory, os.path.basename(interval_list).replace('.interval_list', '.padded.bed'))
                            ),
                            tools.chunkBedbyFileNumber(
                                os.path.join(vardict_directory, os.path.basename(interval_list).replace('.interval_list', '.padded.bed')),
                                vardict_directory,
                                scatter_jobs
                            )
                        ],
                        name="vardict.genome.beds." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
                for idx in range(1, scatter_jobs + 1):
                    bed = os.path.join(vardict_directory, f"{idx:04d}-scattered.bed")
                    output = os.path.join(vardict_directory, f"{tumor_pair.name}.{idx}.vardict.vcf.gz")

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    vardict_directory,
                                    remove=True
                                ),
                                pipe_jobs(
                                    [
                                        vardict.paired_java(
                                            input_normal,
                                            input_tumor,
                                            tumor_pair.name,
                                            None,
                                            bed
                                        ),
                                        vardict.testsomatic(
                                            None,
                                            None
                                        ),
                                        vardict.var2vcf(
                                            None,
                                            tumor_pair.normal.name,
                                            tumor_pair.tumor.name,
                                            None
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output
                                        )
                                    ]
                                )
                            ],
                            name=f"vardict_paired.{tumor_pair.name}.{str(idx)}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def merge_filter_paired_vardict(self):
        """
        The fully merged vcf is filtered using following steps:
        1. Retain only variants designated as somatic by VarDict: either StrongSomatic or LikelySomatic
        2. Somatics identified in step 1 must have PASS filter
        Returns:
            list: A list of merge filter paired vardict jobs.
        """

        jobs = []
        scatter_jobs = global_conf.global_get('gatk_splitInterval', 'scatter_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_prefix = os.path.join(vardict_directory, tumor_pair.name)

            output_tmp = os.path.join(pair_directory, f"{tumor_pair.name}.vardict.tmp.vcf.gz")
            output = os.path.join(pair_directory, f"{tumor_pair.name}.vardict.vcf.gz")
            output_vt = os.path.join(pair_directory, f"{tumor_pair.name}.vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, f"{tumor_pair.name}.vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, f"{tumor_pair.name}.vardict.germline.vt.vcf.gz")

            if scatter_jobs == 1:
                input_file = os.path.join(vardict_directory, f"{tumor_pair.name}.vardict.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                            bash.ln(
                                os.path.relpath(input_file, os.path.dirname(output_tmp)),
                                output_tmp,
                                input=input_file
                            ),
                            pipe_jobs(
                                [
                                    bash.cat(
                                        output_tmp,
                                        None,
                                        zip=True
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    vt.decompose_and_normalize_mnps(
                                        output,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_vt
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_vt,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'somatic_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_somatic
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_vt,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_germline_loh
                                    )
                                ]
                            )
                        ],
                        name=f"symlink_vardict_vcf.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[output_tmp, f"{output_tmp}.tbi", output, f"{output}.tbi"]
                    )
                )
            else:
                merge_list = [f"{output_prefix}.{idx}.vardict.vcf.gz" for idx in range(1, scatter_jobs + 1)]

                for input_vcf in merge_list:
                    if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                        log.error(f"Incomplete vardict vcf: {input_vcf}\n")

                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        merge_list,
                                        None
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $4) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, \"N\", $5) } {print}'"
                                    ),
                                    bash.awk(
                                        None,
                                        None,
                                        "-F$'\\t' -v OFS='\\t' '$1!~/^#/ && $4 == $5 {next} {print}'"
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    vt.decompose_and_normalize_mnps(
                                        output,
                                        None
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_vt
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_vt,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'somatic_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_somatic
                                    )
                                ]
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_vt,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        global_conf.global_get('merge_filter_paired_vardict', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_germline_loh
                                    )
                                ]
                            )
                        ],
                        name=f"merge_filter_paired_vardict.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        removable_files=[output_tmp, f"{output_tmp}.tbi", output, f"{output}.tbi"]
                    )
                )

        return jobs

    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls.
        Filter ensemble calls to retain only calls overlapping 2 or more callers.
        Returns:
            list: A list of ensemble somatic jobs.
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, f"{tumor_pair.name}.mutect2.somatic.vt.vcf.gz")
            input_strelka2 = os.path.join(input_directory, f"{tumor_pair.name}.strelka2.somatic.purple.vcf.gz")
            input_vardict = os.path.join(input_directory, f"{tumor_pair.name}.vardict.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, f"{tumor_pair.name}.varscan2.somatic.vt.vcf.gz")
            inputs_somatic = [input_mutect2, input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_somatic:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete ensemble vcf: {input_vcf}\n")

            output_ensemble = os.path.join(paired_ensemble_directory, f"{tumor_pair.name}.ensemble.somatic.vt.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            paired_ensemble_directory,
                            remove=True
                        ),
                        # Remove any existing outputs because they cause silent error
                        bash.rm(
                            output_ensemble
                        ),
                        bash.rm(
                            output_ensemble + ".tbi"
                        ),
                        bash.rm(
                            os.path.join(paired_ensemble_directory, f"{tumor_pair.name}.ensemble.somatic.vt-work")
                        ),
                        bcbio_variation_recall.ensemble(
                            inputs_somatic,
                            output_ensemble,
                            global_conf.global_get('bcbio_ensemble_somatic', 'options')
                        )
                    ],
                    name=f"bcbio_ensemble_somatic.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=inputs_somatic
                )
            )

        return jobs

    def ensemble_germline_loh(self):
        """
        Apply Bcbio.variations ensemble approach for Vardict, Samtools and VarScan2 calls.
        Filter ensemble calls to retain only calls overlapping 2 or more callers.
        Returns:
            list: A list of ensemble germline jobs.
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_strelka2 = os.path.join(input_directory, f"{tumor_pair.name}.strelka2.germline.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, f"{tumor_pair.name}.vardict.germline.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, f"{tumor_pair.name}.varscan2.germline.vt.vcf.gz")

            inputs_germline = [input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_germline:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete ensemble vcf: {input_vcf}\n")

            output_ensemble = os.path.join(paired_ensemble_directory, f"{tumor_pair.name}.ensemble.germline.vt.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            paired_ensemble_directory,
                            remove=True
                        ),
                        # Remove any existing outputs because they cause silent error
                        bash.rm(
                            output_ensemble
                        ),
                        bash.rm(
                            output_ensemble + ".tbi"
                        ),
                        bash.rm(
                            os.path.join(paired_ensemble_directory, f"{tumor_pair.name}.ensemble.germline.vt-work")
                        ),
                        bcbio_variation_recall.ensemble(
                            inputs_germline,
                            output_ensemble,
                            global_conf.global_get('bcbio_ensemble_germline', 'options')
                        )
                    ],
                    name=f"bcbio_ensemble_germline.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=inputs_germline
                )
            )

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations.
        Returns:
            list: A list of gatk variant annotator somatic jobs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = global_conf.global_get('gatk_variant_annotator', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            annot_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, "rawAnnotation")
            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.cram")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.cram")]
                ]
            )
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, f"{tumor_pair.name}.ensemble.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, f"{tumor_pair.name}.ensemble.somatic.vt.annot.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                annot_directory,
                                remove=True
                            ),
                            gatk.variant_annotator(
                                input_normal,
                                input_tumor,
                                input_somatic_variants,
                                output_somatic_variants,
                                global_conf.global_get('gatk_variant_annotator_somatic', 'other_options')
                            )
                        ],
                        name=f"gatk_variant_annotator_somatic.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_somatic_variants = os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.somatic.vt.annot.{str(idx)}.vcf.gz")

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    annot_directory,
                                    remove=True
                                ),
                                gatk.variant_annotator(
                                    input_normal,
                                    input_tumor,
                                    input_somatic_variants,
                                    output_somatic_variants,
                                    global_conf.global_get('gatk_variant_annotator_somatic', 'other_options'),
                                    intervals=sequences
                                )
                            ],
                            name=f"gatk_variant_annotator_somatic.{str(idx)}.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                output_somatic_variants = os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.somatic.vt.annot.others.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                annot_directory,
                                remove=True
                            ),
                            gatk.variant_annotator(
                                input_normal,
                                input_tumor,
                                input_somatic_variants,
                                output_somatic_variants,
                                global_conf.global_get('gatk_variant_annotator_somatic', 'other_options'),
                                exclude_intervals=unique_sequences_per_job_others
                            )
                        ],
                        name=f"gatk_variant_annotator_somatic.others.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field.
        Returns:
            list: A list of gatk variant annotator germline jobs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = global_conf.global_get('gatk_variant_annotator', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            annot_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble", tumor_pair.name, "rawAnnotation")
            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.cram")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.cram")]
                ]
            )
            input_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, f"{tumor_pair.name}.ensemble.germline.vt.vcf.gz")

            if nb_jobs == 1:
                output_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, f"{tumor_pair.name}.ensemble.germline.vt.annot.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                annot_directory,
                                remove=True
                            ),
                            gatk.variant_annotator(
                                input_normal,
                                input_tumor,
                                input_germline_variants,
                                output_germline_variants,
                                global_conf.global_get('gatk_variant_annotator_germline', 'other_options'),
                            )
                        ],
                        name=f"gatk_variant_annotator_germline.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_germline_variants = os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.germline.vt.annot.{str(idx)}.vcf.gz")

                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    annot_directory,
                                    remove=True
                                ),
                                gatk.variant_annotator(
                                    input_normal,
                                    input_tumor,
                                    input_germline_variants,
                                    output_germline_variants,
                                    global_conf.global_get('gatk_variant_annotator_germline', 'other_options'),
                                    intervals=sequences
                                )
                            ],
                            name=f"gatk_variant_annotator_germline.{str(idx)}.{tumor_pair.name}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                output_germline_variants = os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.germline.vt.annot.others.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                annot_directory,
                                remove=True
                            ),
                            gatk.variant_annotator(
                                input_normal,
                                input_tumor,
                                input_germline_variants,
                                output_germline_variants,
                                global_conf.global_get('gatk_variant_annotator_germline', 'other_options'),
                                exclude_intervals=unique_sequences_per_job_others
                            )
                        ],
                        name=f"gatk_variant_annotator_germline.others.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs.
        Returns:
            list: A list of merge gatk variant annotator somatic jobs.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = global_conf.global_get('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_somatic = os.path.join(ensemble_directory, tumor_pair.name,
                                          f"{tumor_pair.name}.ensemble.somatic.vt.annot.vcf.gz")
            if nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.somatic.vt.annot.{str(idx)}.vcf.gz") for idx in range(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.somatic.vt.annot.others.vcf.gz"))

                jobs.append(
                    pipe_jobs(
                        [
                            bcftools.concat(
                                vcfs_to_merge,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_somatic
                            )
                        ],
                        name=f"merge_gatk_variant_annotator.somatic.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs.
        Returns:
            list: A list of merge gatk variant annotator germline jobs.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = global_conf.global_get('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_germline = os.path.join(ensemble_directory, tumor_pair.name, f"{tumor_pair.name}.ensemble.germline.vt.annot.vcf.gz")

            if nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation", f"{tumor_pair.name}.ensemble.germline.vt.annot.{str(idx)}.vcf.gz") for idx in range(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, f"{tumor_pair.name}.ensemble.germline.vt.annot.others.vcf.gz"))

                jobs.append(
                    pipe_jobs(
                        [
                            bcftools.concat(
                                vcfs_to_merge,
                                None
                            ),
                            htslib.bgzip_tabix(
                                None,
                                output_germline
                            )
                        ],
                        name=f"merge_gatk_variant_annotator.germline.{tumor_pair.name}",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def filter_germline(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        Returns:
            list: A list of filter germline jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if 'fastpass' in self.protocol:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, 'panel')

                input_file = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.germline.vt.prep.vcf.gz"
                )
                output = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.germline.annot.vcf.gz"
                )
                output_filter = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.germline.annot.flt.vcf.gz"
                )

                ini_section = 'filter_fastpass'
                job_name = f"filter_fastpass_germline.{tumor_pair.name}"

            else:
                pair_directory = os.path.join(
                    self.output_dirs['paired_variants_directory'],
                    "ensemble",
                    tumor_pair.name
                )

                input_file = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.ensemble.germline.vt.annot.vcf.gz"
                )
                output = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.ensemble.germline.vt.annot.2caller.vcf.gz"
                )
                output_filter = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.ensemble.germline.vt.annot.2caller.flt.vcf.gz"
                )

                ini_section = 'filter_ensemble'
                job_name = f"filter_ensemble.germline.{tumor_pair.name}"

            jobs.append(
                concat_jobs(
                    [
                        tools.format2pcgr(
                            input_file,
                            output,
                            global_conf.global_get(ini_section, 'call_filter'),
                            "germline",
                            tumor_pair.tumor.name,
                            ini_section=ini_section
                        ),
                        pipe_jobs(
                            [
                                bcftools.view(
                                    output,
                                    None,
                                    filter_options=global_conf.global_get(ini_section, 'germline_filter_options'),
                                ),
                                bcftools.view(
                                    None,
                                    None,
                                    filter_options=f"-Oz -s ^{tumor_pair.normal.name}"
                                ),
                                bcftools.sort(
                                    None,
                                    output_filter,
                                    sort_options="-Oz"
                                )
                            ]
                        ),
                        htslib.tabix(
                            output_filter,
                            options="-pvcf"
                        )
                    ],
                    name=job_name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def filter_somatic(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        Returns:
            list: A list of filter somatic jobs.
        """
        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if 'fastpass' in self.protocol:
                pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, 'panel')

                input_file = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.somatic.vt.vcf.gz"
                )
                output = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.somatic.annot.vcf.gz"
                )
                output_filter = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.varscan2.somatic.annot.flt.vcf.gz"
                )

                ini_section = 'filter_fastpass'
                job_name = f"filter_fastpass_somatic.{tumor_pair.name}"

            else:
                pair_directory = os.path.join(
                    self.output_dirs['paired_variants_directory'],
                    "ensemble",
                    tumor_pair.name
                )

                input_file = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.ensemble.somatic.vt.annot.vcf.gz"
                )
                output = os.path.join(
                    pair_directory,
                    f"{tumor_pair.name}.ensemble.somatic.vt.annot.2caller.vcf.gz"
                )
                output_filter = os.path.join(
                    pair_directory,

                    f"{tumor_pair.name}.ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
                )

                ini_section = 'filter_ensemble'
                job_name = f"filter_ensemble.somatic.{tumor_pair.name}"

            jobs.append(
                concat_jobs(
                    [
                        tools.format2pcgr(
                            input_file,
                            output,
                            global_conf.global_get(ini_section, 'call_filter'),
                            "somatic",
                            tumor_pair.tumor.name,
                            ini_section=ini_section
                        ),
                        bcftools.view(
                            output,
                            output_filter,
                            filter_options=global_conf.global_get(ini_section, 'somatic_filter_options'),
                        ),
                        htslib.tabix(
                            output_filter,
                            options="-pvcf"
                        )
                    ],
                    name=job_name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )

        return jobs

    def sym_link_ensemble(self):
        """
        Create sym link of ensemble output for delivery of data to clients.
        Returns:
            list: A list of sym link ensemble jobs.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] = [os.path.join(self.output_dirs["paired_variants_directory"], "ensemble", tumor_pair.name, tumor_pair.name)]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    f"{sample_prefix}.ensemble.somatic.vt.annot.vcf.gz",
                                    f"{sample_prefix}.ensemble.somatic.vt.annot.vcf.gz.md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.somatic.vt.annot.vcf.gz.md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.somatic.vt.annot.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.somatic.vt.annot.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.md5sum(
                                    f"{sample_prefix}.ensemble.germline.vt.annot.vcf.gz",
                                    f"{sample_prefix}.ensemble.germline.vt.annot.vcf.gz.md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.germline.vt.annot.vcf.gz.md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.germline.vt.annot.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                ),
                                deliverables.sym_link_pair(
                                    f"{sample_prefix}.ensemble.germline.vt.annot.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.profyle
                                )
                            ],
                            name=f"sym_link_ensemble.{str(idx)}.{tumor_pair.name}.{key}",
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def gridss_paired_somatic(self):
        """
        Performs joint variant calling on tumor/normal samples using [GRIDSS](https://github.com/PapenfussLab/gridss),
        followed by filtering with [GRIPSS](https://github.com/hartwigmedical/hmftools/tree/master/gripss).
        GRIPSS applies a set of filtering and post processing steps on GRIDSS paired tumor-normal output to produce 
        a high confidence set of somatic SV for a tumor sample. GRIPSS processes the GRIDSS output and produces a somatic vcf.
        Returns:
            list: A list of gridss paired somatic jobs.
        """
        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, f"{tumor_pair.normal.name}.sorted.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, f"{tumor_pair.tumor.name}.sorted.bam")]
                ]
            )

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            gridss_directory = os.path.join(pair_directory, "gridss")
            normal_output_prefix = os.path.join(gridss_directory, tumor_pair.normal.name)
            tumor_output_prefix = os.path.join(gridss_directory, tumor_pair.tumor.name)
            gridss_vcf_output = f"{tumor_output_prefix}.gridss.vcf.gz"
            gridds_bam_output = f"{tumor_output_prefix}.assembly.bam"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            gridss_directory,
                            remove=False
                        ),
                        gridss.paired_somatic(
                            input_normal,
                            input_tumor,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            gridss_vcf_output,
                            gridds_bam_output
                        )
                    ],
                    name=f"gridss_paired_somatic.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor]
                )
            )

            gripss_vcf_output = f"{tumor_output_prefix}.gripss.somatic.vcf.gz"
            gripss_filter_vcf_output = f"{tumor_output_prefix}.gripss.filtered.somatic.vcf.gz"
            jobs.append(
                concat_jobs(
                    [
                        gripss.filter(
                            gridss_vcf_output,
                            gripss_vcf_output,
                            gripss_filter_vcf_output,
                            "somatic",
                            sample=tumor_pair.tumor.name,
                            reference=tumor_pair.normal.name,
                        )
                    ],
                    name=f"gripss_filter.somatic.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                )
            )

            gripss_vcf_output = f"{normal_output_prefix}.gripss.germline.vcf.gz"
            gripss_filter_vcf_output = f"{normal_output_prefix}.gripss.filtered.germline.vcf.gz"
            jobs.append(
                concat_jobs(
                    [
                        gripss.filter(
                            gridss_vcf_output,
                            gripss_vcf_output,
                            gripss_filter_vcf_output,
                            "germline -germline",
                            sample=tumor_pair.normal.name,
                            reference=tumor_pair.tumor.name
                        )
                    ],
                    name=f"gripss_filter.germline.{tumor_pair.name}",
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                )
            )
        return jobs

    def purple_sv(self):
        """
        Runs PURPLE with the optional structural variant input VCFs.
        PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

        It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
        somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.
        Returns:
            list: A list of purple structural variant jobs.
        """
        jobs = self.purple(sv=True)
        return jobs

    def linx_annotations_somatic(self):
        """
        [Linx](https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md) is an annotation, interpretation and visualisation tool for structural variants.
        The primary function of Linx is grouping together individual SV calls into distinct events 
        and properly classify and annotating the event to understand both its mechanism and genomic impact.
        Returns:
            list: A list of linx annotations somatic jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            purple_dir = os.path.join(pair_dir, "purple")
            purple_vcf = os.path.join(purple_dir, f"{tumor_pair.tumor.name}.purple.sv.vcf.gz")
            linx_output_dir = os.path.join(pair_dir, "linx")
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(linx_output_dir),
                        linx.somatic(
                            purple_vcf,
                            purple_dir,
                            tumor_pair.tumor.name,
                            linx_output_dir,
                            ini_section="linx_annotations_somatic"
                        )
                    ],
                    name=f"linx_annotations_somatic.{tumor_pair.name}",
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
        return jobs

    def linx_annotations_germline(self):
        """
        Runs [Linx](https://github.com/hartwigmedical/hmftools/blob/master/linx/README.md) in germline mode.
        Linx is an annotation, interpretation and visualisation tool for structural variants.
        The primary function of Linx is grouping together individual SV calls into distinct events 
        and properly classify and annotating the event to understand both its mechanism and genomic impact.
        Returns:
            list: A list of linx annotations germline jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            purple_dir = os.path.join(pair_dir, "purple")
            purple_vcf = os.path.join(purple_dir, f"{tumor_pair.tumor.name}.purple.sv.vcf.gz")
            linx_output_dir = os.path.join(pair_dir, "linx")
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(linx_output_dir),
                        linx.germline(
                            purple_vcf,
                            tumor_pair.tumor.name,
                            linx_output_dir,
                            ini_section="linx_annotations_germline"
                        )
                    ],
                    name=f"linx_annotations_germline.{tumor_pair.name}",
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
        return jobs

    def linx_plot(self):
        """
        Generate Linx Plot of the tumor pair analysis.
        Returns:
            list: A list of linx plot jobs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            linx_output_dir = os.path.join(pair_dir, "linx")
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(linx_output_dir),
                        linx.plot(
                            tumor_pair.tumor.name,
                            linx_output_dir,
                            ini_section="linx_plot"
                        ),
                        bash.touch(os.path.join(linx_output_dir, f"linx_plot.{tumor_pair.name}.Done"))
                    ],
                    name="linx_plot." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets),
                    output_dependency=[os.path.join(linx_output_dir, f"linx_plot.{tumor_pair.name}.Done")]
                )
            )
        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        """
        Returns the protocol for the pipeline.
        Returns:
            dict: A dictionary of protocols for the pipeline.
        """
        return {'germline_snv':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.set_interval_list,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf,
                self.combine_gvcf,
                self.merge_and_call_combined_gvcf,
                self.variant_recalibrator,
                self.haplotype_caller_decompose_and_normalize,
                self.haplotype_caller_flag_mappability,
                self.haplotype_caller_snp_id_annotation,
                self.haplotype_caller_snp_effect,
                self.haplotype_caller_dbnsfp_annotation,
                self.haplotype_caller_gemini_annotations,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.metrics_picard_calculate_hs,
                self.metrics_verify_bam_id,
                self.run_multiqc,
                self.sym_link_fastq,
                self.sym_link_final_bam,
                self.metrics_vcftools_missing_indiv,
                self.metrics_vcftools_depth_indiv,
                self.metrics_gatk_sample_fingerprint,
                self.metrics_gatk_cluster_fingerprint
            ], 'germline_sv':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.sym_link_final_bam,  # 5
                self.set_interval_list,
                self.gatk_haplotype_caller,
                self.merge_and_call_individual_gvcf,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.metrics_picard_calculate_hs,
                self.run_multiqc,
                self.delly_call_filter,
                self.delly_sv_annotation,
                self.germline_manta,
                self.manta_sv_annotation,
                self.lumpy_paired_sv,
                self.lumpy_sv_annotation,
                self.wham_call_sv,
                self.wham_sv_annotation,
                self.cnvkit_batch,
                self.cnvkit_sv_annotation,
                self.run_breakseq2,
	            self.ensemble_metasv,
                self.metasv_sv_annotation
            ], 'germline_high_cov':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.samtools_merge_files,
                self.gatk_fixmate,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.metrics_picard_calculate_hs,
                self.metrics_verify_bam_id,
                self.germline_varscan2,
                self.preprocess_vcf,
                self.snp_effect,
                self.gemini_annotations,
                self.run_multiqc,
                self.cram_output
            ], 'somatic_tumor_only':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.sym_link_final_bam,  # 5
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.metrics_picard_calculate_hs,
                self.metrics_verify_bam_id,
                self.run_multiqc,
                self.set_interval_list,
                self.gatk_haplotype_caller, 
                self.merge_and_call_individual_gvcf,
                self.combine_gvcf,
                self.merge_and_call_combined_gvcf,
                self.variant_recalibrator,
                self.haplotype_caller_decompose_and_normalize,
                self.cnvkit_batch,
                self.split_tumor_only,
                self.filter_tumor_only,
                self.report_cpsr,
                self.report_pcgr
            ], 'somatic_fastpass':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.set_interval_list,
                self.sequenza,
                self.rawmpileup,
                self.paired_varscan2,
                self.merge_varscan2,
                self.preprocess_vcf,
                self.cnvkit_batch,
                self.filter_germline,
                self.report_cpsr,
                self.filter_somatic,
                self.report_pcgr,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.run_multiqc,
                self.sym_link_report,
                self.sym_link_fastq_pair,
                self.sym_link_panel,
                self.cram_output
            ], 'somatic_ensemble':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.set_interval_list,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_mosdepth,
                self.sequenza,
                self.manta_sv_calls,
                self.strelka2_paired_somatic,
                self.strelka2_paired_germline,
                self.strelka2_paired_snpeff,
                self.purple,
                self.rawmpileup,
                self.paired_varscan2,
                self.merge_varscan2,
                self.paired_mutect2,
                self.merge_mutect2,
                self.vardict_paired,
                self.merge_filter_paired_vardict,
                self.ensemble_somatic,
                self.gatk_variant_annotator_somatic,
                self.merge_gatk_variant_annotator_somatic,
                self.ensemble_germline_loh,
                self.gatk_variant_annotator_germline,
                self.merge_gatk_variant_annotator_germline,
                self.cnvkit_batch,
                self.filter_germline,
                self.report_cpsr,
                self.filter_somatic,
                self.report_pcgr,
                self.run_multiqc,
                self.sym_link_fastq_pair,
                self.sym_link_final_bam,
                self.sym_link_report,
                self.sym_link_ensemble,
                self.cram_output
            ], 'somatic_sv':
            [
                self.gatk_sam_to_fastq,
                self.trim_fastp,
                self.bwa_mem2_samtools_sort,
                self.gatk_mark_duplicates,
                self.set_interval_list,
                self.manta_sv_calls,
                self.strelka2_paired_somatic,
                # self.sv_prep,
                self.gridss_paired_somatic,
                self.purple_sv,
                self.linx_annotations_somatic,
                self.linx_annotations_germline,
                self.linx_plot,
                self.run_multiqc,
                self.cram_output
            ]
        }
class DnaSeq(DnaSeqRaw):
    """
    Arguments:
        config_files (list): A list of configuration files.
        genpipes_file (str): The genpipes configuration file.
        steps (str): The steps to run.
        readsets_file (str): The readsets file.
        clean (bool): Clean the output directory.
        force (bool): Force the pipeline to run.
        job_scheduler (str): The job scheduler to use.
        output_dir (str): The output directory.
        protocol (str): The DNAseq analysis type.
        design_file (str): The design file.
        no_json (bool): Do not create a JSON file.
        container (str): The container to use.
        profyle (bool): Adjust deliverables to PROFYLE folder conventions.
        pairs_file (str): The pairs file.
    """
    __doc__ = DnaSeqRaw.__doc__
    def __init__(self, *args, protocol="germline_snv", **kwargs):
        self._protocol = protocol
        # Add pipeline specific arguments
        super(DnaSeq, self).__init__(*args, **kwargs)

    @classmethod
    def argparser(cls, argparser):
        super().argparser(argparser)
        cls._argparser.add_argument("-p", "--pairs", help="pairs file", type=argparse.FileType('r'))
        cls._argparser.add_argument("--profyle", help=("adjust deliverables to PROFYLE folder conventions (Default: False)"), action="store_true")
        cls._argparser.add_argument("-t", "--type", help="DNAseq analysis type", dest='protocol',
                                    choices=["germline_snv", "germline_sv", "germline_high_cov", "somatic_tumor_only", "somatic_fastpass", "somatic_ensemble", "somatic_sv"], default="germline_snv")
        return cls._argparser

def main(parsed_args):
    """
    The function that will call this pipeline!
    Parameters:
        parsed_args (argparse.Namespace): The parsed arguments from the command line.
    """

    # Pipeline config
    config_files = parsed_args.config
    # Genpipes Config

    # Pipeline options
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
    protocol = parsed_args.protocol
    design_file = parsed_args.design_file
    profyle = parsed_args.profyle
    pairs_file = parsed_args.pairs

    pipeline = DnaSeq(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, protocol=protocol, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container, profyle=profyle, pairs_file=pairs_file)

    pipeline.submit_jobs()
