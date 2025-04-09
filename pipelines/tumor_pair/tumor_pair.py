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
import gzip

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, pipe_jobs
from core.sample_tumor_pairs import parse_tumor_pair_file
import utils.utils

from pipelines.dnaseq import dnaseq

#utilizes
from bfx import (
    adapters,
    amber,
    annotations,
    bash_cmd as bash,
    bcbio_variation_recall,
    bcftools,
    bvatools,
    cnvkit,
    cobalt,
    conpair,
    cpsr,
    deliverables,
    delly,
    djerba,
    fastqc,
    gatk,
    gatk4,
    gemini,
    gridss,
    gripss,
    htslib,
    job2json_project_tracking,
    linx,
    manta,
    multiqc,
    pcgr,
    purple,
    qualimap,
    sambamba,
    samtools,
    sequence_dictionary,
    sequenza,
    snpeff,
    strelka2,
    tools,
    vardict,
    varscan,
    vawk,
    vt
    )

log = logging.getLogger(__name__)

class TumorPair(dnaseq.DnaSeqRaw):
    """
    Tumor Pair Pipeline
    =================

    The Tumor Pair pipeline inherits the initial bam preparation steps of the DNA-Seq pipeline except for the indel realignment (IR) step.
    In the tumor pipeline, the IR step utilizes both the normal and tumor bam to further reduce false positives (FPs) in and around indels.
    The tumor pipeline deviates from the DNA-seq pipeline at the variant calling step. At this point, a paired caller is used to call SNVs
    and Indels from the pairs given as input. Additionally, multiple cancer callers are utilized using an ensemble approach and SNVs and
    Indels seen in at least 2 different callers are retained for further investigation.

    Example command:
    python tumor_pair.py -c a.ini b.base.ini -s x-y,z -r readset.tsv -p pairs.csv

    -c ini files: multiple can be specified e.g WGS or exome, or different clusters e.g. base (abacus) or guillimin

    -r readset: derived from GQ lims or made yourself. See : https://genpipes.readthedocs.io/en/latest/get-started/concepts/readset_file.html#docs-readset-file

    -p pairs : format - patient_name,normal_sample_name,tumor_sample_name
    """

    def __init__(self, protocol=None):
        self._protocol = protocol
        self.argparser.add_argument("-p", "--pairs", help="pairs file", type=argparse.FileType('r'))
        self.argparser.add_argument("--profyle", help="adjust deliverables to PROFYLE folder conventions (Default: False)", action="store_true")
        self.argparser.add_argument("-t", "--type", help="Tumor pair analysis type", choices=["fastpass", "ensemble", "sv"], default="ensemble")
        super(TumorPair, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
            'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'paired_variants_directory': os.path.relpath(os.path.join(self.output_dir, 'pairedVariants'), self.output_dir),
            'sv_variants_directory': os.path.relpath(os.path.join(self.output_dir, 'SVariants'), self.output_dir),
            'report': {}
        }
        for tumor_pair in self.tumor_pairs.values():
            dirs['report'][tumor_pair.name] = os.path.relpath(os.path.join(self.output_dir, 'report', tumor_pair.name), self.output_dir)
        return dirs

    @property
    def multiqc_inputs(self):
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = {}
            for tumor_pair in self.tumor_pairs.values():
                self._multiqc_inputs[tumor_pair.name] = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    @property
    def tumor_pairs(self):
        if not hasattr(self, "_tumor_pairs"):
            self._tumor_pairs = parse_tumor_pair_file(
                self.args.pairs.name,
                self.samples,
                self.args.profyle
            )
        return self._tumor_pairs

    def sequence_dictionary_variant(self):
        if not hasattr(self, "_sequence_dictionary_variant"):
            self._sequence_dictionary_variant = sequence_dictionary.parse_sequence_dictionary_file(
                config.param('DEFAULT', 'genome_dictionary', param_type='filepath'),
                variant=True
            )
        return self._sequence_dictionary_variant

    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary_variant()):
            return [sequence['name'] + ":1-" + str(sequence['length']) for sequence in self.sequence_dictionary_variant()]
        else:
            total_length = sum([sequence['length'] for sequence in self.sequence_dictionary_variant()])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary_variant()))))
            windows = []

            for sequence in self.sequence_dictionary_variant():
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence['length'])] for pos in range(1, sequence['length'] + 1, approximate_window_size)]:
                    windows.append(sequence['name'] + ":" + str(start) + "-" + str(end))

        return windows

    def is_gz_file(self, name):
        if not os.path.isfile(name):
            return True
        #if os.stat(name).st_size == 0:
        #    return False

        with gzip.open(name, 'rb') as f:
            try:
                file_content = f.read(1)
                return len(file_content) > 0
            except:
                return False

    def sym_link_fastq_pair(self):
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
                            os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair1.fastq.gz"),
                            os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair2.fastq.gz")
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
                            os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair1.fastq.gz"),
                            os.path.join(self.output_dirs["raw_reads_directory"], readset.sample.name, readset.name + ".pair2.fastq.gz")
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
                        profyle=self.args.profyle
                    )
                    dir_name, file_name = os.path.split(symlink_pair_job.output_files[0])
                    # do not compute md5sum in the readset input directory
                    md5sum_job = deliverables.md5sum(
                        symlink_pair_job.output_files[0],
                        file_name + ".md5",
                        dir_name
                    )
                    jobs.append(
                        concat_jobs(
                            [
                                symlink_pair_job,
                                md5sum_job
                            ],
                            name="sym_link_fastq.pairs." + str(read) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.tumor, tumor_pair.normal],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.

        Note: modified to use both normal and tumor bams to reduce FPs around indels.
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            quality_offsets = self.readsets[0].quality_offset
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            pair_directory = os.path.join(self.output_dirs['alignment_directory'], "realign", tumor_pair.name)

            input_normal = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.normal.name + ".sorted.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_intervals = os.path.join(pair_directory, "all.intervals")
                bam_postfix = ".realigned.all.bam"

                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.all.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                normal_output_bam = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)

                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.all.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                tumor_output_bam = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                pair_directory,
                                remove=True
                            ),
                            bash.chdir(
                                pair_directory
                            ),
                            gatk.realigner_target_creator(
                                input_normal,
                                realign_intervals,
                                output_dir=self.output_dir,
                                input2=input_tumor,
                                fix_encoding=True if quality_offsets == 64 else ""
                            ),
                            gatk.indel_realigner(
                                input_normal,
                                input2=input_tumor,
                                output_dir=self.output_dir,
                                output_norm_dep=[normal_bam,normal_index],
                                output_tum_dep=[tumor_bam,tumor_index],
                                target_intervals=realign_intervals,
                                optional=bam_postfix,
                                fix_encoding=True if quality_offsets == 64 else ""
                            ),
                            bash.chdir(
                                self.output_dir
                            ),
                            # Move sample realign
                            bash.ln(
                                os.path.relpath(normal_bam, os.path.dirname(normal_output_bam)),
                                normal_output_bam,
                                input=normal_bam
                            ),
                            bash.ln(
                                os.path.relpath(normal_index, os.path.dirname(normal_output_index)),
                                normal_output_index,
                                input=normal_index
                            ),
                            bash.ln(
                                os.path.relpath(tumor_bam, os.path.dirname(tumor_output_bam)),
                                tumor_output_bam,
                                input=tumor_bam
                            ),
                            bash.ln(
                                os.path.relpath(tumor_index, os.path.dirname(tumor_output_index)),
                                tumor_output_index,
                                input=tumor_index
                            )
                        ],
                        name="gatk_indel_realigner." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(
                    self.sequence_dictionary,
                    nb_jobs - 1
                )
                normal_realign_directory = os.path.join(normal_alignment_directory, "realign")
                tumor_realign_directory = os.path.join(tumor_alignment_directory, "realign")

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):
                    realign_prefix = os.path.join(pair_directory, str(idx))
                    realign_intervals = realign_prefix + ".intervals"
                    intervals = sequences
                    if str(idx) == 0:
                        intervals.append("unmapped")
                    bam_postfix = ".realigned." + str(idx) + ".bam"
                    normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_index = re.sub("\.bam$", ".bai", normal_bam)
                    tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                    normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam")
                    normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                    tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam")
                    tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                    jobs.append(
                        concat_jobs(
                            [
                                # Create output directory since it is not done by default by GATK tools
                                bash.mkdir(
                                    pair_directory,
                                    remove=True
                                ),
                                bash.mkdir(
                                    normal_realign_directory,
                                    remove=True
                                ),
                                bash.mkdir(
                                    tumor_realign_directory,
                                    remove=True
                                ),
                                bash.chdir(
                                    pair_directory
                                ),
                                gatk.realigner_target_creator(
                                    input_normal,
                                    realign_intervals,
                                    output_dir=self.output_dir,
                                    input2=input_tumor,
                                    intervals=intervals,
                                    fix_encoding=True if quality_offsets == 64 else ""
                                ),
                                gatk.indel_realigner(
                                    input_normal,
                                    input2=input_tumor,
                                    output_dir=self.output_dir,
                                    output_norm_dep=[normal_bam,normal_index],
                                    output_tum_dep=[tumor_bam,tumor_index],
                                    target_intervals=realign_intervals,
                                    intervals=intervals,
                                    optional=bam_postfix,
                                    fix_encoding=True if quality_offsets == 64 else ""
                                ),
                                bash.chdir(
                                    self.output_dir
                                ),
                                # Move sample realign
                                bash.ln(
                                    os.path.relpath(normal_bam, os.path.dirname(normal_output_bam)),
                                    normal_output_bam,
                                    input=normal_bam
                                ),
                                bash.ln(
                                    os.path.relpath(normal_index, os.path.dirname(normal_output_index)),
                                    normal_output_index,
                                    input=normal_index
                                ),
                                bash.ln(
                                    os.path.relpath(tumor_bam, os.path.dirname(tumor_output_bam)),
                                    tumor_output_bam,
                                    input=tumor_bam
                                ),
                                bash.ln(
                                    os.path.relpath(tumor_index, os.path.dirname(tumor_output_index)),
                                    tumor_output_index,
                                    input=tumor_index
                                )
                            ],
                            name="gatk_indel_realigner." + tumor_pair.name + "." + str(idx),
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_intervals = os.path.join(pair_directory, "others.intervals")
                bam_postfix = ".realigned.others.bam"
                normal_bam = os.path.join(pair_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_index = re.sub("\.bam$", ".bai", normal_bam)
                tumor_bam = os.path.join(pair_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_index = re.sub("\.bam$", ".bai", tumor_bam)
                normal_output_bam = os.path.join(normal_realign_directory, tumor_pair.normal.name + ".sorted.realigned.others.bam")
                normal_output_index = re.sub("\.bam$", ".bai", normal_output_bam)
                tumor_output_bam = os.path.join(tumor_realign_directory, tumor_pair.tumor.name + ".sorted.realigned.others.bam")
                tumor_output_index = re.sub("\.bam$", ".bai", tumor_output_bam)

                jobs.append(
                    concat_jobs(
                        [
                            # Create output directory since it is not done by default by GATK tools
                            bash.mkdir(
                                pair_directory,
                                remove=True
                            ),
                            bash.mkdir(
                                normal_realign_directory,
                                remove=True
                            ),
                            bash.mkdir(
                                tumor_realign_directory,
                                remove=True
                            ),
                            bash.chdir(
                                pair_directory
                            ),
                            gatk.realigner_target_creator(
                                input_normal,
                                realign_intervals,
                                output_dir=self.output_dir,
                                input2=input_tumor,
                                exclude_intervals=unique_sequences_per_job_others,
                                fix_encoding=True if quality_offsets == 64 else ""
                            ),
                            gatk.indel_realigner(
                                input_normal,
                                input2=input_tumor,
                                output_dir=self.output_dir,
                                output_norm_dep=[normal_bam, normal_index],
                                output_tum_dep=[tumor_bam, tumor_index],
                                target_intervals=realign_intervals,
                                exclude_intervals=unique_sequences_per_job_others,
                                optional=bam_postfix,
                                fix_encoding=True if quality_offsets == 64 else ""
                            ),
                            bash.chdir(
                                self.output_dir
                            ),
                            bash.ln(
                                os.path.relpath(normal_bam, os.path.dirname(normal_output_bam)),
                                normal_output_bam,
                                input=normal_bam
                            ),
                            bash.ln(
                                os.path.relpath(normal_index, os.path.dirname(normal_output_index)),
                                normal_output_index,
                                input=normal_index
                            ),
                            bash.ln(
                                os.path.relpath(tumor_bam, os.path.dirname(tumor_output_bam)),
                                tumor_output_bam,
                                input=tumor_bam
                            ),
                            bash.ln(
                                os.path.relpath(tumor_index, os.path.dirname(tumor_output_index)),
                                tumor_output_index,
                                input=tumor_index
                            )
                        ],
                        name="gatk_indel_realigner." + tumor_pair.name + ".others",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def sambamba_merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(self.sequence_dictionary, nb_jobs - 1)

                normal_inputs = []
                for idx, _ in enumerate(unique_sequences_per_job):
                    normal_inputs.append(
                        os.path.join(
                            normal_alignment_directory,
                            "realign",
                            tumor_pair.normal.name + ".sorted.realigned." + str(idx) + ".bam"
                        )
                    )
                normal_inputs.append(
                    os.path.join(
                        normal_alignment_directory,
                        "realign",
                        tumor_pair.normal.name + ".sorted.realigned.others.bam"
                    )
                )

                tumor_inputs = []
                for idx, _ in enumerate(unique_sequences_per_job):
                    tumor_inputs.append(
                        os.path.join(
                            tumor_alignment_directory,
                            "realign",
                            tumor_pair.tumor.name + ".sorted.realigned." + str(idx) + ".bam"
                        )
                    )
                tumor_inputs.append(
                    os.path.join(
                        tumor_alignment_directory,
                        "realign",
                        tumor_pair.tumor.name + ".sorted.realigned.others.bam"
                    )
                )

                job = sambamba.merge(
                    normal_inputs,
                    os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")
                )
                job.name = "sambamba_merge_realigned." + tumor_pair.name + "." + tumor_pair.normal.name
                job.samples = [tumor_pair.normal]
                job.readsets = list(tumor_pair.normal.readsets)
                jobs.append(job)

                job = sambamba.merge(
                    tumor_inputs,
                    os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")
                )
                job.name = "sambamba_merge_realigned." + tumor_pair.name + "." + tumor_pair.tumor.name
                job.samples = [tumor_pair.tumor]
                job.readsets = list(tumor_pair.tumor.readsets)
                jobs.append(job)

        return jobs

    def sambamba_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Sambamba](http://lomereiter.github.io/sambamba/index.html).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            [normal_input] = self.select_input_files([
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")],
            ])
            normal_output = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")

            [tumor_input] = self.select_input_files([
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")],
            ])
            tumor_output = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name  + ".sorted.dup.bam")

            job = sambamba.markdup(
                normal_input,
                normal_output,
                config.param('sambamba_mark_duplicates', 'tmp_dir'),
                other_options=config.param('sambamba_mark_duplicates', 'options')
            )
            job.name = "sambamba_mark_duplicates." + tumor_pair.name + "." + tumor_pair.normal.name
            job.samples = [tumor_pair.normal]
            job.readsets = list(tumor_pair.normal.readsets)
            jobs.append(job)

            job = sambamba.markdup(
                tumor_input,
                tumor_output,
                config.param('sambamba_mark_duplicates', 'tmp_dir'),
                other_options=config.param('sambamba_mark_duplicates', 'options')
            )
            job.name = "sambamba_mark_duplicates." + tumor_pair.name + "." + tumor_pair.tumor.name
            job.samples = [tumor_pair.tumor]
            job.readsets = list(tumor_pair.tumor.readsets)
            jobs.append(job)

        return jobs

    def recalibration(self):
        """
        Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
        the quality scores in the QUAL field in each read in the output BAM are more accurate in that
        the reported quality score is closer to its actual probability of mismatching the reference genome.
        Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
        and sequence context, and by doing so, provides not only more accurate quality scores but also
        more widely dispersed ones.
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            normal_prefix = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.")
            tumor_prefix = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.")

            normal_input = normal_prefix + "bam"
            tumor_input = tumor_prefix + "bam"

            normal_print_reads_output = normal_prefix + "recal.bam"
            tumor_print_reads_output = tumor_prefix + "recal.bam"

            normal_base_recalibrator_output = normal_prefix + "recalibration_report.grp"
            tumor_base_recalibrator_output = tumor_prefix + "recalibration_report.grp"

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )
            if coverage_bed:
                interval_list = os.path.join(tumor_alignment_directory, re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed)))

                if not os.path.isfile(interval_list):
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(tumor_alignment_directory),
                                tools.bed2interval_list(
                                    coverage_bed,
                                    interval_list
                                )
                            ],
                            name="interval_list." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

            job = gatk4.base_recalibrator(
                normal_input,
                normal_base_recalibrator_output,
                intervals=interval_list
            )
            job.name = "gatk_base_recalibrator." + tumor_pair.name + "." + tumor_pair.normal.name
            job.samples = [tumor_pair.normal]
            job.readsets = list(tumor_pair.normal.readsets)
            jobs.append(job)
            if config.param('gatk_print_reads', 'module_gatk').split("/")[2] > "4":
                job = concat_jobs(
                        [
                            gatk4.print_reads(
                                normal_input,
                                normal_print_reads_output,
                                normal_base_recalibrator_output
                            ),
                            deliverables.md5sum(
                                normal_print_reads_output,
                                normal_print_reads_output + ".md5",
                                self.output_dir
                            )
                        ]
                    )
            else:
                job = gatk4.print_reads(
                    normal_input,
                    normal_print_reads_output,
                    normal_base_recalibrator_output
                )
            job.name = "gatk_print_reads." + tumor_pair.name + "." + tumor_pair.normal.name
            job.samples = [tumor_pair.normal]
            job.readsets = list(tumor_pair.normal.readsets)
            jobs.append(job)

            job = gatk4.base_recalibrator(
                tumor_input,
                tumor_base_recalibrator_output,
                intervals=interval_list
            )
            job.name = "gatk_base_recalibrator." + tumor_pair.name + "." + tumor_pair.tumor.name
            job.samples = [tumor_pair.tumor]
            job.readsets = list(tumor_pair.tumor.readsets)
            jobs.append(job)

            if config.param('gatk_print_reads', 'module_gatk').split("/")[2] > "4":
                job = concat_jobs(
                        [
                            gatk4.print_reads(
                                tumor_input,
                                tumor_print_reads_output,
                                tumor_base_recalibrator_output
                            ),
                            deliverables.md5sum(
                                tumor_print_reads_output,
                                tumor_print_reads_output + ".md5",
                                self.output_dir
                            )
                        ]
                    )
            else:
                job = gatk4.print_reads(
                    tumor_input,
                    tumor_print_reads_output,
                    tumor_base_recalibrator_output
                )
            job.name = "gatk_print_reads." + tumor_pair.name + "." + tumor_pair.tumor.name
            job.samples = [tumor_pair.tumor]
            job.readsets = list(tumor_pair.tumor.readsets)
            jobs.append(job)

        return jobs

    def sym_link_final_bam(self):
        """
        Create sym link of the final bam for delivery of data to clients.
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
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam"),
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bai")
                        ],
                        [
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam"),
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam.bai")
                        ],
                        [
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bai")
                        ],
                        [
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam.bai")
                        ]
                    ]
                )
            ]

            [inputs["Tumor"]] = [
                self.select_input_files(
                    [
                        [
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam"),
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bai")
                        ],
                        [
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam"),
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam.bai")
                        ],
                        [
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam"),
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bai")
                        ],
                        [
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam"),
                            os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam.bai")
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
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="alignment",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="alignment",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_final_bam.pairs." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def conpair_concordance_contamination(self):
        """
        Conpair is a fast and robust method dedicated to human tumor-normal studies to perform concordance verification
        (= samples coming from the same individual), as well as cross-individual contamination level estimation in
        whole-genome and whole-exome sequencing experiments. Importantly, the method of estimating contamination in
        the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            metrics_directory = self.output_dirs['metrics_directory']

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")
            pileup_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".gatkPileup")
            pileup_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".gatkPileup")

            concordance_out = os.path.join(metrics_directory, tumor_pair.tumor.name + ".concordance.tsv")
            contamination_out = os.path.join(metrics_directory, tumor_pair.tumor.name + ".contamination.tsv")

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
                    name="conpair_concordance_contamination.pileup." + tumor_pair.name + "." + tumor_pair.normal.name,
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
                    name="conpair_concordance_contamination.pileup." + tumor_pair.name + "." + tumor_pair.tumor.name,
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
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
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
                        bash.ln(
                            os.path.relpath(concordance_out, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(concordance_out)),
                            input=concordance_out
                        ),
                        bash.ln(
                            os.path.relpath(contamination_out, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(contamination_out)),
                            input=contamination_out
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
                    concordance_out,
                    contamination_out
                ]
            )

        return jobs

    def rawmpileup_panel(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')
            bedfile = config.param('rawmpileup_panel', 'panel')

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            samtools.mpileup(
                                [
                                    os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                                    os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")
                                ],
                                input_pair,
                                config.param('rawmpileup_panel', 'mpileup_other_options'),
                                regionFile=bedfile
                            )
                        ],
                        name="rawmpileup_panel." + tumor_pair.name + ".all",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=True
                                    ),
                                    samtools.mpileup(
                                        [
                                            os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam"),
                                            os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name, tumor_pair.tumor.name + ".sorted.dup.bam")
                                        ],
                                        pair_output,
                                        config.param('rawmpileup_panel', 'mpileup_other_options'),
                                        region=sequence['name'],
                                        regionFile=bedfile
                                    )
                                ],
                                name="rawmpileup_panel." + tumor_pair.name + "." + sequence['name'],
                                samples=[tumor_pair.tumor, tumor_pair.normal],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                            )
                        )
        return jobs

    def paired_varscan2_panel(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")

                output = os.path.join(varscan_directory, tumor_pair.name)
                output_snp = os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                config.param('varscan2_somatic_panel', 'other_options'),
                                output_vcf_dep=output_vcf_gz,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                            ),
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        [
                                            os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz"),
                                            os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                                        ],
                                        None
                                    ),
                                    pipe_jobs(
                                        [
                                            bash.sed(
                                                None,
                                                None,
                                                "'s/TUMOR/"+ tumor_pair.tumor.name + "/g'"
                                            ),
                                            bash.sed(
                                                None,
                                                None,
                                                "'s/NORMAL/"+ tumor_pair.normal.name + "/g'"
                                            )
                                        ]
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_vcf_gz
                                    )
                                ]
                            )
                        ],
                        name="varscan2_somatic_panel." + tumor_pair.name + ".all",
                        samples=[tumor_pair.tumor, tumor_pair.normal],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:

                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                        output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                        output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=True
                                    ),
                                    varscan.somatic(
                                        input_pair,
                                        output,
                                        config.param('varscan2_somatic_panel', 'other_options'),
                                        output_vcf_dep=output_vcf_gz,
                                        output_snp_dep=output_snp,
                                        output_indel_dep=output_indel
                                    ),
                                    htslib.bgzip_tabix(
                                        output_snp,
                                        os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz")
                                    ),
                                    htslib.bgzip_tabix(
                                        output_indel,
                                        os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")
                                    ),
                                    pipe_jobs(
                                        [
                                            bcftools.concat(
                                                [
                                                    os.path.join(varscan_directory, tumor_pair.name + ".snp." + sequence['name'] + ".vcf.gz"),
                                                    os.path.join(varscan_directory, tumor_pair.name + ".indel." + sequence['name'] + ".vcf.gz")
                                                ],
                                                None
                                            ),
                                            pipe_jobs(
                                                [
                                                    bash.sed(
                                                        None,
                                                        None,
                                                        "'s/TUMOR/"+ tumor_pair.tumor.name + "/g'"
                                                    ),
                                                    bash.sed(
                                                        None,
                                                        None,
                                                        "'s/NORMAL/"+ tumor_pair.normal.name + "/g'"
                                                    )
                                                ]
                                            ),
                                            htslib.bgzip_tabix(
                                                None,
                                                output_vcf_gz
                                            )
                                        ]
                                    )
                                ],
                                name="varscan2_somatic_panel." + tumor_pair.name + "." + sequence['name'],
                                samples=[tumor_pair.tumor, tumor_pair.normal],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                            )
                        )
        return jobs

    def merge_varscan2_panel(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup_panel', 'nb_jobs', param_type='posint')

            if nb_jobs == 1:
                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bash.cat(
                                        os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                        None,
                                        zip=True
                                    ),
                                    tools.fix_varscan_output(
                                        None,
                                        None,
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
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                    )
                                ]
                            ),
                            bcftools.view(
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                                config.param('merge_varscan2', 'somatic_filter_options')
                            ),
                            htslib.tabix(
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                                config.param('merge_varscan2', 'tabix_options', required=False)
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                        None,
                                        config.param('merge_varscan2', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('merge_varscan2', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                                    )
                                ]
                            )
                        ],
                        name="merge_varscan2." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                all_inputs = [
                    os.path.join(varscan_directory, tumor_pair.name + ".varscan2." + sequence['name'] + ".vcf.gz")
                    for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary'
                ]

                for input_vcf in all_inputs:
                    if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                        log.error(f"Incomplete panel varscan2 vcf: {input_vcf}\n")

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
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                    )
                                ]
                            ),
                            bcftools.view(
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                                config.param('merge_varscan2', 'somatic_filter_options')
                            ),
                            htslib.tabix(
                                os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vcf.gz"),
                                config.param('merge_varscan2', 'tabix_options', required=False)
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz"),
                                        None,
                                        config.param('merge_varscan2', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('merge_varscan2', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vcf.gz"),
                                    )
                                ]
                            )
                        ],
                        name="merge_varscan2." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
        return jobs

    def preprocess_vcf_panel(self):
        """
        Preprocess vcf for loading into an annotation database - Gemini : http://gemini.readthedocs.org/en/latest/index.html
        Processes include normalization and decomposition of MNPs by vt (http://genome.sph.umich.edu/wiki/Vt) and
        vcf FORMAT modification for correct loading into Gemini.
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")

            prefix = os.path.join(pair_directory, tumor_pair.name)
            output_somatic = prefix + ".varscan2.somatic.vt.vcf.gz"

            output_germline = prefix + ".varscan2.germline.vt.vcf.gz"

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    prefix + ".varscan2.somatic.vcf.gz" ,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    prefix + ".prep.vt.vcf.gz"
                                )
                            ]
                        ),
                        tools.preprocess_varscan(
                            prefix + ".prep.vt.vcf.gz",
                            output_somatic
                        )
                    ],
                    name="preprocess_vcf_panel.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    prefix + ".varscan2.germline.vcf.gz" ,
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    prefix + ".germline.prep.vt.vcf.gz"
                                )
                            ]
                        ),
                        tools.preprocess_varscan(
                            prefix + ".germline.prep.vt.vcf.gz",
                            output_germline
                        )
                    ],
                    name="preprocess_vcf_panel.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def snp_effect_panel(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        [SnpEff](https://pcingola.github.io/SnpEff/) annotates and predicts the effects of variants on genes (such as amino acid changes).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            input_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf")
            output_somatic_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.snpeff.vcf.gz")

            input_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")
            output_germline = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf")
            output_germline_gz = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.snpeff.vcf.gz")

            cancer_pair_filename = os.path.join(varscan_directory, tumor_pair.name + '.tsv')

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(varscan_directory),
                        Job(
                            [input_somatic, input_germline],
                            [cancer_pair_filename],
                            command="""\
echo -e "{normal_name}\\t{tumor_name}" \\
  > {cancer_pair_filename}""".format(
                                normal_name=tumor_pair.normal.name,
                                tumor_name=tumor_pair.tumor.name,
                                cancer_pair_filename=cancer_pair_filename
                            )
                        )
                    ],
                    name="compute_cancer_effects_somatic.file." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    removable_files=[varscan_directory]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            input_somatic,
                            output_somatic,
                            cancer_sample_file=cancer_pair_filename,
                            ini_section='compute_cancer_effects_somatic',
                            options=config.param('compute_cancer_effects_somatic', 'options')
                        ),
                        htslib.bgzip_tabix(
                            output_somatic,
                            output_somatic_gz
                        )
                    ],
                    name = "compute_cancer_effects_somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            input_germline,
                            output_germline,
                            cancer_sample_file=cancer_pair_filename,
                            ini_section='compute_cancer_effects_germline',
                            options=config.param('compute_cancer_effects_germline', 'options')
                        ),
                        htslib.bgzip_tabix(
                            output_germline,
                            output_germline_gz
                        )
                    ],
                    name="compute_cancer_effects_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
        return jobs

    def gemini_annotations_panel(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database [Gemini] (http://gemini.readthedocs.org/en/latest/index.html).
        """

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel")

            temp_dir = config.param('DEFAULT', 'tmp_dir')
            gemini_prefix = os.path.join(pair_directory, tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        gemini.gemini_annotations(
                            gemini_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                            gemini_prefix + ".somatic.gemini.db", temp_dir
                        )
                    ],
                    name="gemini_annotations.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        gemini.gemini_annotations(
                            gemini_prefix + ".varscan2.germline.vt.snpeff.vcf.gz",
                            gemini_prefix + ".germline.gemini.db",
                            temp_dir
                        )
                    ],
                    name="gemini_annotations.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sym_link_panel(self):
        """
        Create sym links of panel variants for deliverables to the clients.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] =  [os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "panel", tumor_pair.name)]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.vcf.gz",
                                    tumor_pair, self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.somatic.vt.snpeff.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.germline.vt.snpeff.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".varscan2.germline.vt.snpeff.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".somatic.gemini.db",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".germline.gemini.db",
                                    tumor_pair, self.output_dir,
                                    type="snv/panel",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_panel." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def metrics_dna_picard_metrics(self):
        """
        Runs specific QC metrics on DNA data.
        Functions: collect_multiple_metrics, CollectOxoGMetrics and collect_sequencing_artifacts_metrics
        [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html).
        """

        ffpe = config.param('picard_collect_sequencing_artifacts_metrics', 'FFPE', param_type='boolean')

        ##check the library status
        library = {}
        for readset in self.readsets:
            if not readset.sample in library:
                library[readset.sample] = "SINGLE_END"
            if readset.run_type == "PAIRED_END":
                library[readset.sample] = "PAIRED_END"

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            normal_picard_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics, "picard_metrics")
            tumor_picard_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "picard_metrics")

            [normal_input] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")],
                ]
            )

            [tumor_input] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")],
                ]
            )

            mkdir_job_normal = bash.mkdir(
                normal_picard_directory,
                remove=True
            )

            tumor_pair_jobs = []

            normal_samples = [tumor_pair.normal]
            normal_job_name = f"picard_collect_multiple_metrics.{tumor_pair.name}.{tumor_pair.normal.name}"
            normal_output_prefix = os.path.join(normal_picard_directory, tumor_pair.normal.name + ".all.metrics")
            normal_job_project_tracking_metrics = []
            if self.project_tracking_json:
                normal_job_project_tracking_metrics = concat_jobs(
                    [
                    gatk4.parse_bases_over_q30_percent_metrics_pt(f"{normal_output_prefix}.quality_distribution_metrics"),
                    job2json_project_tracking.run(
                        input_file=f"{normal_output_prefix}.quality_distribution_metrics",
                        pipeline=self,
                        samples=",".join([sample.name for sample in normal_samples]),
                        readsets=",".join([readset.name for sample in normal_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="bases_over_q30_percent=$bases_over_q30_percent"
                        )
                    ])

            collect_multiple_metrics_normal_job = concat_jobs(
                [
                    mkdir_job_normal,
                    gatk4.collect_multiple_metrics(
                        normal_input,
                        normal_output_prefix,
                        library_type=library[tumor_pair.normal]
                        ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                        ),
                    normal_job_project_tracking_metrics
                ]
            )
            for outfile in collect_multiple_metrics_normal_job.report_files:
                self.multiqc_inputs[tumor_pair.name].append(outfile)
                collect_multiple_metrics_normal_job = concat_jobs(
                    [
                        collect_multiple_metrics_normal_job,
                        bash.ln(
                            os.path.relpath(outfile, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(outfile)),
                            input=outfile
                        )
                    ]
                )
            collect_multiple_metrics_normal_job.name = normal_job_name
            collect_multiple_metrics_normal_job.samples = normal_samples
            collect_multiple_metrics_normal_job.readsets = list(tumor_pair.normal.readsets)
            tumor_pair_jobs.append(collect_multiple_metrics_normal_job)

            collect_oxog_metrics_normal_job = concat_jobs(
                [
                    mkdir_job_normal,
                    gatk4.collect_oxog_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name + ".oxog_metrics.txt")
                    ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                    ),
                    bash.ln(
                        os.path.relpath(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".oxog_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                        os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.normal.name + ".oxog_metrics.txt"),
                        input=os.path.join(normal_picard_directory, tumor_pair.normal.name + ".oxog_metrics.txt")
                        )
                ]
            )
            self.multiqc_inputs[tumor_pair.name].append(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".oxog_metrics.txt"))
            collect_oxog_metrics_normal_job.name = "picard_collect_oxog_metrics." + tumor_pair.name + "." + tumor_pair.normal.name
            collect_oxog_metrics_normal_job.samples = [tumor_pair.normal]
            collect_oxog_metrics_normal_job.readsets = list(tumor_pair.normal.readsets)
            tumor_pair_jobs.append(collect_oxog_metrics_normal_job)

            collect_wgs_metrics_normal_job = concat_jobs(
                    [
                        mkdir_job_normal,
                        gatk4.collect_wgs_metrics(
                            normal_input,
                            os.path.join(normal_picard_directory, tumor_pair.normal.name + ".wgs_metrics.txt")
                            ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                            ),
                        bash.ln(
                            os.path.relpath(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".wgs_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.normal.name + ".wgs_metrics.txt"),
                            input=os.path.join(normal_picard_directory, tumor_pair.normal.name + ".wgs_metrics.txt")
                            )
                    ]
            )
            self.multiqc_inputs[tumor_pair.name].append(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".wgs_metrics.txt"))
            collect_wgs_metrics_normal_job.name = "picard_collect_wgs_metrics." + tumor_pair.name + "." + tumor_pair.normal.name
            collect_wgs_metrics_normal_job.samples = [tumor_pair.normal]
            collect_wgs_metrics_normal_job.readsets = list(tumor_pair.normal.readsets)
            tumor_pair_jobs.append(collect_wgs_metrics_normal_job)

            collect_gcbias_metrics_normal_job = concat_jobs(
                [
                    mkdir_job_normal,
                    gatk4.collect_gcbias_metrics(
                        normal_input,
                        os.path.join(normal_picard_directory, tumor_pair.normal.name)
                    ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                    )
                ]
            )
            for outfile in collect_gcbias_metrics_normal_job.report_files:
                self.multiqc_inputs[tumor_pair.name].append(outfile)
                collect_gcbias_metrics_normal_job = concat_jobs(
                    [
                        collect_gcbias_metrics_normal_job,
                        bash.ln(
                            os.path.relpath(outfile, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(outfile)),
                            input=outfile
                        )
                    ]
                )
            collect_gcbias_metrics_normal_job.name = "picard_collect_gcbias_metrics." + tumor_pair.name + "." + tumor_pair.normal.name
            collect_gcbias_metrics_normal_job.samples = [tumor_pair.normal]
            collect_gcbias_metrics_normal_job.readsets = list(tumor_pair.normal.readsets)
            tumor_pair_jobs.append(collect_gcbias_metrics_normal_job)

            mkdir_job_tumor = bash.mkdir(
                tumor_picard_directory,
                remove=True
            )

            tumor_samples = [tumor_pair.tumor]
            tumor_job_name = f"picard_collect_multiple_metrics.{tumor_pair.name}.{tumor_pair.tumor.name}"
            tumor_output_prefix = os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".all.metrics")
            tumor_job_project_tracking_metrics = []
            if self.project_tracking_json:
                tumor_job_project_tracking_metrics = concat_jobs(
                    [
                    gatk4.parse_bases_over_q30_percent_metrics_pt(f"{tumor_output_prefix}.quality_distribution_metrics"),
                    job2json_project_tracking.run(
                        input_file=f"{tumor_output_prefix}.quality_distribution_metrics",
                        pipeline=self,
                        samples=",".join([sample.name for sample in tumor_samples]),
                        readsets=",".join([readset.name for sample in tumor_samples for readset in sample.readsets]),
                        job_name=tumor_job_name,
                        metrics="bases_over_q30_percent=$bases_over_q30_percent"
                        )
                    ])

            collect_multiple_metrics_tumor_job = concat_jobs(
                [
                    mkdir_job_tumor,
                    gatk4.collect_multiple_metrics(
                        tumor_input,
                        tumor_output_prefix,
                        library_type=library[tumor_pair.tumor]
                    ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                    ),
                    tumor_job_project_tracking_metrics
                ]
            )
            for outfile in collect_multiple_metrics_tumor_job.report_files:
                self.multiqc_inputs[tumor_pair.name].append(outfile)
                collect_multiple_metrics_tumor_job = concat_jobs(
                    [
                        collect_multiple_metrics_tumor_job,
                        bash.ln(
                            os.path.relpath(outfile, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(outfile)),
                            input=outfile
                        )
                    ]
                )
            collect_multiple_metrics_tumor_job.name = tumor_job_name
            collect_multiple_metrics_tumor_job.samples = tumor_samples
            collect_multiple_metrics_tumor_job.readsets = list(tumor_pair.tumor.readsets)
            tumor_pair_jobs.append(collect_multiple_metrics_tumor_job)

            collect_oxog_metrics_tumor_job = concat_jobs(
                [
                    mkdir_job_tumor,
                    gatk4.collect_oxog_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".oxog_metrics.txt")
                    ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                    ),
                    bash.ln(
                        os.path.relpath(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".oxog_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                        os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.tumor.name + ".oxog_metrics.txt"),
                        input=os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".oxog_metrics.txt")
                        )
                ]
            )
            self.multiqc_inputs[tumor_pair.name].append(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".oxog_metrics.txt"))
            collect_oxog_metrics_tumor_job.name = "picard_collect_oxog_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name
            collect_oxog_metrics_tumor_job.samples = [tumor_pair.tumor]
            collect_oxog_metrics_tumor_job.readsets = list(tumor_pair.tumor.readsets)
            tumor_pair_jobs.append(collect_oxog_metrics_tumor_job)

            collect_wgs_metrics_tumor_job = concat_jobs(
                    [
                        mkdir_job_tumor,
                        gatk4.collect_wgs_metrics(
                            tumor_input,
                            os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".wgs_metrics.txt")
                            ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                            ),
                        bash.ln(
                            os.path.relpath(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".wgs_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.tumor.name + ".wgs_metrics.txt"),
                            input=os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".wgs_metrics.txt")
                            )
                    ]
            )
            self.multiqc_inputs[tumor_pair.name].append(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".wgs_metrics.txt"))
            collect_wgs_metrics_tumor_job.name = "picard_collect_wgs_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name
            collect_wgs_metrics_tumor_job.samples = [tumor_pair.tumor]
            collect_wgs_metrics_tumor_job.readsets = list(tumor_pair.tumor.readsets)
            tumor_pair_jobs.append(collect_wgs_metrics_tumor_job)

            collect_gcbias_metrics_tumor_job = concat_jobs(
                [
                    mkdir_job_tumor,
                    gatk4.collect_gcbias_metrics(
                        tumor_input,
                        os.path.join(tumor_picard_directory, tumor_pair.tumor.name),
                    ),
                    bash.mkdir(
                        self.output_dirs['report'][tumor_pair.name]
                    )
                ]
            )
            for outfile in collect_gcbias_metrics_tumor_job.report_files:
                self.multiqc_inputs[tumor_pair.name].append(outfile)
                collect_gcbias_metrics_tumor_job = concat_jobs(
                    [
                        collect_gcbias_metrics_tumor_job,
                        bash.ln(
                            os.path.relpath(outfile, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(outfile)),
                            input=outfile
                        )
                    ]
                )
            collect_gcbias_metrics_tumor_job.name = "picard_collect_gcbias_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name
            collect_gcbias_metrics_tumor_job.samples = [tumor_pair.tumor]
            collect_gcbias_metrics_tumor_job.readsets = list(tumor_pair.tumor.readsets)
            tumor_pair_jobs.append(collect_gcbias_metrics_tumor_job)

            if ffpe:
                collect_sequencing_artifacts_metrics_normal_job = concat_jobs(
                    [
                        mkdir_job_normal,
                        gatk4.collect_sequencing_artifacts_metrics(
                            normal_input,
                            os.path.join(normal_picard_directory, tumor_pair.normal.name)
                        ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                        ),
                        bash.ln(
                            os.path.relpath(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".bait_bias_summary_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.normal.name + ".bait_bias_summary_metrics.txt"),
                            input=os.path.join(normal_picard_directory, tumor_pair.normal.name + ".bait_bias_summary_metrics.txt")
                            )
                    ]
                )
                self.multiqc_inputs[tumor_pair.name].append(os.path.join(normal_picard_directory, tumor_pair.normal.name + ".bait_bias_summary_metrics.txt"))
                collect_sequencing_artifacts_metrics_normal_job.name = "picard_collect_sequencing_artifacts_metrics." + tumor_pair.name + "." + tumor_pair.normal.name
                collect_sequencing_artifacts_metrics_normal_job.samples = [tumor_pair.normal]
                collect_sequencing_artifacts_metrics_normal_job.readsets = list(tumor_pair.normal.readsets)
                tumor_pair_jobs.append(collect_sequencing_artifacts_metrics_normal_job)

                collect_sequencing_artifacts_metrics_tumor_job = concat_jobs(
                    [
                        mkdir_job_tumor,
                        gatk4.collect_sequencing_artifacts_metrics(
                            tumor_input,
                            os.path.join(tumor_picard_directory, tumor_pair.tumor.name)
                        ),
                        bash.mkdir(self.output_dirs['report'][tumor_pair.name]),
                        bash.ln(
                            os.path.relpath(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".bait_bias_summary_metrics.txt"), self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], tumor_pair.tumor.name + ".bait_bias_summary_metrics.txt"),
                            input=os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".bait_bias_summary_metrics.txt")
                            )
                    ]
                )
                self.multiqc_inputs[tumor_pair.name].append(os.path.join(tumor_picard_directory, tumor_pair.tumor.name + ".bait_bias_summary_metrics.txt"))
                collect_sequencing_artifacts_metrics_tumor_job.name = "picard_collect_sequencing_artifacts_metrics." + tumor_pair.name + "." + tumor_pair.tumor.name
                collect_sequencing_artifacts_metrics_tumor_job.samples = [tumor_pair.tumor]
                collect_sequencing_artifacts_metrics_tumor_job.readsets = list(tumor_pair.tumor.readsets)
                tumor_pair_jobs.append(collect_sequencing_artifacts_metrics_tumor_job)

            jobs.extend(tumor_pair_jobs)
        return jobs

    def metrics_dna_sample_qualimap(self):
        """
        QC alignment metrics generated by [Qualimap](http://qualimap.conesalab.org/).
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            normal_qualimap_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics, "qualimap", tumor_pair.normal.name)
            tumor_qualimap_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "qualimap", tumor_pair.tumor.name)

            [normal_input] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )
            normal_output = os.path.join(normal_qualimap_directory, "genome_results.txt")
            normal_output_histogram = os.path.join(normal_qualimap_directory, "raw_data_qualimapReport", "insert_size_histogram.txt")

            [tumor_input] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )
            tumor_output = os.path.join(tumor_qualimap_directory, "genome_results.txt")
            tumor_output_histogram = os.path.join(tumor_qualimap_directory, "raw_data_qualimapReport", "insert_size_histogram.txt")

            use_bed = config.param('dna_sample_qualimap', 'use_bed', param_type='boolean', required=True)
            options = None
            if use_bed:
                bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
                options = config.param('dna_sample_qualimap', 'qualimap_options') + " --feature-file " + bed
            else:
                options = config.param('dna_sample_qualimap', 'qualimap_options')

            normal_samples = [tumor_pair.normal]
            normal_job_name = f"dna_sample_qualimap.{tumor_pair.name}.{tumor_pair.normal.name}"
            normal_job_project_tracking_metrics = []
            if self.project_tracking_json:
                normal_job_project_tracking_metrics = concat_jobs(
                    [
                    qualimap.parse_median_insert_size_metrics_pt(normal_output),
                    job2json_project_tracking.run(
                        input_file=normal_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in normal_samples]),
                        readsets=",".join([readset.name for sample in normal_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="median_insert_size=$median_insert_size"
                        ),
                    qualimap.parse_mean_insert_size_metrics_pt(normal_output_histogram),
                    job2json_project_tracking.run(
                        input_file=normal_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in normal_samples]),
                        readsets=",".join([readset.name for sample in normal_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="mean_insert_size=$mean_insert_size"
                        ),
                    qualimap.parse_dedup_coverage_metrics_pt(normal_output),
                    job2json_project_tracking.run(
                        input_file=normal_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in normal_samples]),
                        readsets=",".join([readset.name for sample in normal_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="dedup_coverage=$dedup_coverage"
                        ),
                    qualimap.parse_aligned_reads_count_metrics_pt(normal_output),
                    job2json_project_tracking.run(
                        input_file=normal_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in normal_samples]),
                        readsets=",".join([readset.name for sample in normal_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="aligned_reads_count=$aligned_reads_count"
                        )
                    ])

            tumor_pair_jobs = []
            qualimap_normal_job = concat_jobs(
                [
                    bash.mkdir(
                        normal_qualimap_directory,
                        remove=False
                    ),
                    qualimap.bamqc(
                        normal_input,
                        normal_qualimap_directory,
                        normal_output,
                        options
                    ),
                    bash.mkdir(
                       os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap")
                    ),
                    bash.ln(
                        os.path.relpath(normal_qualimap_directory, os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap")),
                        os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap", tumor_pair.normal.name),
                        input=normal_qualimap_directory
                    ),
                    normal_job_project_tracking_metrics
                ]
            )
            self.multiqc_inputs[tumor_pair.name].append(normal_qualimap_directory)
            qualimap_normal_job.name = normal_job_name
            qualimap_normal_job.samples = normal_samples
            qualimap_normal_job.readsets = list(tumor_pair.normal.readsets)
            tumor_pair_jobs.append(qualimap_normal_job)

            tumor_samples = [tumor_pair.tumor]
            tumor_job_name = f"dna_sample_qualimap.{tumor_pair.name}.{tumor_pair.tumor.name}"
            tumor_job_project_tracking_metrics = []
            if self.project_tracking_json:
                tumor_job_project_tracking_metrics = concat_jobs(
                    [
                    qualimap.parse_median_insert_size_metrics_pt(tumor_output),
                    job2json_project_tracking.run(
                        input_file=tumor_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in tumor_samples]),
                        readsets=",".join([readset.name for sample in tumor_samples for readset in sample.readsets]),
                        job_name=tumor_job_name,
                        metrics="median_insert_size=$median_insert_size"
                        ),
                    qualimap.parse_mean_insert_size_metrics_pt(tumor_output_histogram),
                    job2json_project_tracking.run(
                        input_file=tumor_output_histogram,
                        pipeline=self,
                        samples=",".join([sample.name for sample in tumor_samples]),
                        readsets=",".join([readset.name for sample in tumor_samples for readset in sample.readsets]),
                        job_name=normal_job_name,
                        metrics="mean_insert_size=$mean_insert_size"
                        ),
                    qualimap.parse_dedup_coverage_metrics_pt(tumor_output),
                    job2json_project_tracking.run(
                        input_file=tumor_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in tumor_samples]),
                        readsets=",".join([readset.name for sample in tumor_samples for readset in sample.readsets]),
                        job_name=tumor_job_name,
                        metrics="dedup_coverage=$dedup_coverage"
                        ),
                    qualimap.parse_aligned_reads_count_metrics_pt(tumor_output),
                    job2json_project_tracking.run(
                        input_file=tumor_output,
                        pipeline=self,
                        samples=",".join([sample.name for sample in tumor_samples]),
                        readsets=",".join([readset.name for sample in tumor_samples for readset in sample.readsets]),
                        job_name=tumor_job_name,
                        metrics="aligned_reads_count=$aligned_reads_count"
                        )
                    ])

            qualimap_tumor_job = concat_jobs(
                [
                    bash.mkdir(
                        tumor_qualimap_directory,
                        remove=False
                    ),
                    qualimap.bamqc(
                        tumor_input,
                        tumor_qualimap_directory,
                        tumor_output,
                        options
                    ),
                    bash.mkdir(
                       os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap")
                    ),
                    bash.ln(
                        os.path.relpath(tumor_qualimap_directory, os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap")),
                        os.path.join(self.output_dirs['report'][tumor_pair.name], "qualimap", tumor_pair.tumor.name),
                        input=tumor_qualimap_directory
                    ),
                    tumor_job_project_tracking_metrics
                ]
            )
            self.multiqc_inputs[tumor_pair.name].append(tumor_qualimap_directory)
            qualimap_tumor_job.name = tumor_job_name
            qualimap_tumor_job.samples = tumor_samples
            qualimap_tumor_job.readsets = list(tumor_pair.tumor.readsets)
            tumor_pair_jobs.append(qualimap_tumor_job)

            jobs.extend(tumor_pair_jobs)
        return jobs

    def metrics_dna_fastqc(self):
        """
        QCing metrics are generated on the read level using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
                normal_metrics = os.path.join(tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
                normal_metrics = os.path.join(tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)
            normal_fastqc_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", normal_metrics, "fastqc")
            tumor_fastqc_directory = os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.tumor.name, "fastqc")

            [normal_input] = self.select_input_files(
                [
                    # [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.realigned.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            normal_output_dir = normal_fastqc_directory
            normal_file = re.sub(".bam", "", os.path.basename(normal_input))
            normal_output = os.path.join(normal_fastqc_directory, normal_file + "_fastqc.zip")

            [tumor_input] = self.select_input_files(
                [
                    # [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.realigned.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            tumor_output_dir = tumor_fastqc_directory
            tumor_file = re.sub(".bam", "", os.path.basename(tumor_input))
            tumor_output = os.path.join(tumor_fastqc_directory, tumor_file + "_fastqc.zip")

            adapter_file = config.param('fastqc', 'adapter_file', required=False, param_type='filepath')
            normal_adapter_job = None
            tumor_adapter_job = None

            if not adapter_file:
                normal_adapter_job = adapters.create(
                    tumor_pair.normal.readsets[0],
                    os.path.join(normal_output_dir, "adapter.tsv"),
                    fastqc=True
                )
                tumor_adapter_job = adapters.create(
                    tumor_pair.tumor.readsets[0],
                    os.path.join(tumor_output_dir, "adapter.tsv"),
                    fastqc=True
                )

            tumor_pair_jobs = []
            tumor_pair_jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            normal_output_dir,
                            remove=True
                        ),
                        normal_adapter_job,
                        fastqc.fastqc(
                            normal_input,
                            None,
                            normal_output_dir,
                            normal_output,
                            os.path.join(normal_output_dir, "adapter.tsv")
                        ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                        ),
                        bash.ln(
                            os.path.relpath(normal_output, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(normal_output)),
                            input=normal_output
                        )
                    ],
                    name="fastqc." + tumor_pair.name + "." + tumor_pair.normal.name,
                    samples=[tumor_pair.normal],
                    readsets=list(tumor_pair.normal.readsets)
                )
            )
            tumor_pair_jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            tumor_output_dir,
                            remove=True
                        ),
                        tumor_adapter_job,
                        fastqc.fastqc(
                            tumor_input,
                            None,
                            tumor_output_dir,
                            tumor_output,
                            os.path.join(tumor_output_dir, "adapter.tsv")
                        ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                        ),
                        bash.ln(
                            os.path.relpath(tumor_output, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(tumor_output)),
                            input=tumor_output
                        )
                    ],
                    name="fastqc." + tumor_pair.name + "." + tumor_pair.tumor.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
            for job in tumor_pair_jobs:
                self.multiqc_inputs[tumor_pair.name].extend(job.output_files)
            jobs.extend(tumor_pair_jobs)
        return jobs

    def run_pair_multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report
        MultiQC searches a given directory for analysis logs and compiles an HTML report. It's a general-use tool,
        perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).
        """

        jobs = []

        metrics_directory = os.path.join(self.output_dirs['metrics_directory'], "dna")
        for tumor_pair in self.tumor_pairs.values():
            output = os.path.join(metrics_directory, tumor_pair.name + ".multiqc")
            job = multiqc.run(
                self.output_dirs['report'][tumor_pair.name],
                output
            )
            job.name = "multiqc." + tumor_pair.name
            job.samples = [tumor_pair.normal, tumor_pair.tumor]
            job.readsets = [*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
            job.input_files = self.multiqc_inputs[tumor_pair.name]
            jobs.append(job)
        return jobs

    def sym_link_report(self):
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] = [os.path.join(self.output_dirs['metrics_directory'], "dna", tumor_pair.name + ".multiqc.html")]

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
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_fastq.report." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def rawmpileup(self):
        """
        Full pileup (optional). A raw mpileup file is created using samtools mpileup and compressed in gz format.
        One packaged mpileup file is created per sample/chromosome.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                bed_file = coverage_bed

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            if nb_jobs == 1:
                pair_output = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")
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
                                config.param('rawmpileup', 'mpileup_other_options'),
                                regionFile=bed_file
                            )
                        ],
                        name="rawmpileup." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        pair_output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

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
                                        config.param('rawmpileup', 'mpileup_other_options'),
                                        region=sequence['name'],
                                        regionFile=bed_file
                                    )
                                ],
                                name="rawmpileup." + tumor_pair.name + "." + sequence['name'],
                                samples=[tumor_pair.normal, tumor_pair.tumor],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                            )
                        )

        return jobs

    def paired_varscan2(self):
        """
        Variant calling and somatic mutation/CNV detection for next-generation sequencing data.
        Koboldt et al., 2012. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing
        Varscan2 thresholds based on DREAM3 results generated by the author see: https://github.com/dkoboldt/varscan/releases
        SSC INFO field remove to prevent collision with Samtools output during ensemble.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")
            output = os.path.join(varscan_directory, tumor_pair.name)

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            if nb_jobs == 1:
                input_pair = os.path.join(varscan_directory, tumor_pair.name + ".mpileup")

                output_snp = os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf")
                output_indel = os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf")
                output_vcf = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf")
                output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                varscan_directory,
                                remove=True
                            ),
                            varscan.somatic(
                                input_pair,
                                output,
                                config.param('varscan2_somatic', 'other_options'),
                                output_vcf_dep=output_vcf,
                                output_snp_dep=output_snp,
                                output_indel_dep=output_indel
                            ),
                            htslib.bgzip_tabix(
                                output_snp,
                                os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz")
                            ),
                            htslib.bgzip_tabix(
                                output_indel,
                                os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
                            ),
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        [
                                            os.path.join(varscan_directory, tumor_pair.name + ".snp.vcf.gz"),
                                            os.path.join(varscan_directory, tumor_pair.name + ".indel.vcf.gz")
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
                        name="varscan2_somatic." + tumor_pair.name,
                        samples=[tumor_pair.tumor, tumor_pair.normal],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:

                for sequence in self.sequence_dictionary_variant():
                    if sequence['type'] == 'primary':
                        input_pair = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".mpileup")

                        output = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'])
                        output_snp = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf")
                        output_indel = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf")
                        output_vcf = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf")
                        output_vcf_gz = os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")

                        jobs.append(
                            concat_jobs(
                                [
                                    bash.mkdir(
                                        varscan_directory,
                                        remove=True
                                    ),
                                    varscan.somatic(
                                        input_pair,
                                        output,
                                        config.param('varscan2_somatic', 'other_options'),
                                        output_vcf_dep=output_vcf,
                                        output_snp_dep=output_snp,
                                        output_indel_dep=output_indel
                                    ),
                                    htslib.bgzip_tabix(
                                        output_snp,
                                        os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz")
                                    ),
                                    htslib.bgzip_tabix(
                                        output_indel,
                                        os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")
                                    ),
                                    pipe_jobs(
                                        [
                                            bcftools.concat(
                                                [
                                                    os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".snp.vcf.gz"),
                                                    os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".indel.vcf.gz")
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
                                name="varscan2_somatic." + tumor_pair.name + "." + sequence['name'],
                                samples=[tumor_pair.tumor, tumor_pair.normal],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                            )
                        )
        return jobs

    def merge_varscan2(self):
        """
        Merge mpileup files per sample/chromosome into one compressed gzip file per sample.
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            varscan_directory = os.path.join(pair_directory, "rawVarscan2")

            nb_jobs = config.param('rawmpileup', 'nb_jobs', param_type='posint')
            if nb_jobs > 50:
                log.warning(
                    "Number of mpileup jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

            all_inputs = []
            if nb_jobs == 1:
                all_inputs = [os.path.join(varscan_directory, tumor_pair.name + ".varscan2.vcf.gz")]

            else:
                all_inputs = [
                    os.path.join(varscan_directory, tumor_pair.name + "." + sequence['name'] + ".varscan2.vcf.gz")
                    for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary'
                ]

            for input_vcf in all_inputs:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete varscan2 vcf: {input_vcf}\n")

            all_output = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vcf.gz")
            all_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.vt.vcf.gz")

            somtic_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            germline_output_vt = os.path.join(pair_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")

            if nb_jobs == 1:
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
                                    #vt.sort("-", all_output, "-m full"),
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
                                config.param('merge_varscan2', 'somatic_filter_options')
                            ),
                            htslib.tabix(
                                somtic_output_vt,
                                config.param('merge_varscan2', 'tabix_options')
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        all_output_vt,
                                        None,
                                        config.param('merge_varscan2', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('merge_varscan2', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        germline_output_vt
                                    )
                                ]
                            )
                        ],
                        name="merge_varscan2." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
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
                                    #vt.sort("-", all_output, "-m full"),
                                    htslib.bgzip_tabix(
                                        None,
                                        all_output
                                    )
                                ]
                            ),
                            #htslib.tabix(all_output),
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
                                        config.param('varscan2_readcount_fpfilter', 'somatic_filter_options')
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
                                        config.param('varscan2_readcount_fpfilter', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('varscan2_readcount_fpfilter', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        germline_output_vt
                                    )
                                ]
                            )
                        ],
                        name="merge_varscan2." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
        return jobs

    def paired_mutect2(self):
        """
        GATK MuTect2 caller for SNVs and Indels.
        """

        jobs = []

        created_interval_lists = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of mutect jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            interval_list = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(tumor_pair.normal.readsets[0])
            if coverage_bed:
                interval_list = os.path.join(mutect_directory, re.sub("\.[^.]+$", ".interval_list", os.path.basename(coverage_bed)))

                if not interval_list in created_interval_lists:
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(mutect_directory),
                                tools.bed2interval_list(
                                    coverage_bed,
                                    interval_list
                                )
                            ],
                            name="interval_list." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
                    created_interval_lists.append(interval_list)

            if nb_jobs == 1:

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
                                os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"),
                                os.path.join(mutect_directory, tumor_pair.name + ".f1r2.tar.gz"),
                                interval_list=interval_list
                            )
                        ],
                        name="gatk_mutect2." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)

                # Create one separate job for each of the first sequences
                for idx, sequences in enumerate(unique_sequences_per_job):

                    outprefix = tumor_pair.name + "." + str(idx) + ".mutect2"
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
                                    os.path.join(mutect_directory, outprefix + ".vcf.gz"),
                                    os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".f1r2.tar.gz"),
                                    intervals=sequences,
                                    interval_list=interval_list
                                )
                            ],
                            name="gatk_mutect2." + tumor_pair.name + "." + str(idx),
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                # Create one last job to process the last remaining sequences and 'others' sequences
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
                                os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"),
                                os.path.join(mutect_directory, tumor_pair.name + ".others.f1r2.tar.gz"),
                                exclude_intervals=unique_sequences_per_job_others,
                                interval_list=interval_list
                            )
                        ],
                        name="gatk_mutect2." + tumor_pair.name + ".others",
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def merge_mutect2(self):
        """
        Merge SNVs and indels for mutect2.
        Replace TUMOR and NORMAL sample names in vcf to the exact tumor/normal sample names
        Generate a somatic vcf containing only PASS variants.
        """

        jobs = []

        nb_jobs = config.param('gatk_mutect2', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            mutect_directory = os.path.join(pair_directory, "rawMuTect2")
            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            output_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz")
            output_flt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.flt.vcf.gz")
            output_vt_gz = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vt.vcf.gz")
            output_somatic_vt = os.path.join(pair_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                if config.param('gatk_mutect2', 'module_gatk').split("/")[2] > "4":
                    jobs.append(
                        concat_jobs(
                            [
                                gatk4.learn_read_orientation_model(
                                    [os.path.join(mutect_directory, tumor_pair.name + ".f1r2.tar.gz")],
                                    os.path.join(pair_directory, tumor_pair.name + ".f1r2.tar.gz")
                                ),
                                gatk4.filter_mutect_calls(
                                    os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz"),
                                    output_flt,
                                    read_orientation=os.path.join(pair_directory, tumor_pair.name + ".f1r2.tar.gz")
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
                                                    "-e 's#/\.##g'"
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
                                            config.param('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name="merge_filter_mutect2." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                else:
                    input_vcf = os.path.join(mutect_directory, tumor_pair.name + ".mutect2.vcf.gz")
                    jobs.append(
                        concat_jobs(
                            [
                                bash.ln(
                                    os.path.relpath(input_vcf, os.path.dirname(output_gz)),
                                    output_gz,
                                    input=input_vcf
                                ),
                                #gatk4.filter_mutect_calls(output_gz, output_flt),
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
                                                    "-e 's#/\.##g'"
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
                            name="symlink_mutect_vcf." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

            elif nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(
                    self.sequence_dictionary_variant(), nb_jobs - 1)

                # Create one separate job for each of the first sequences
                inputs = []
                for idx, _ in enumerate(unique_sequences_per_job):
                    inputs.append(os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz"))
                inputs.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz"))

                for input_vcf in inputs:
                    if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                        log.error(f"Incomplete mutect2 vcf: {input_vcf}\n")

                if config.param('gatk_mutect2', 'module_gatk').split("/")[2] > "4":

                    output_stats = os.path.join(pair_directory, tumor_pair.name + ".mutect2.vcf.gz.stats")
                    stats = []
                    for idx, _ in enumerate(unique_sequences_per_job):
                        stats.append(
                            os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".mutect2.vcf.gz.stats"))
                    stats.append(os.path.join(mutect_directory, tumor_pair.name + ".others.mutect2.vcf.gz.stats"))

                    output_models = os.path.join(pair_directory, tumor_pair.name + ".read-orientation-model.tar.gz")
                    models = []
                    for idx, sequences in enumerate(unique_sequences_per_job):
                        models.append(
                            os.path.join(mutect_directory, tumor_pair.name + "." + str(idx) + ".f1r2.tar.gz"))
                    models.append(os.path.join(mutect_directory, tumor_pair.name + ".others.f1r2.tar.gz"))

                    jobs.append(
                        concat_jobs(
                            [
                                gatk4.learn_read_orientation_model(
                                    models,
                                    output_models
                                ),
                                gatk4.cat_variants(
                                    inputs,
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
                                                    "-e 's#/\.##g'"
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
                                            config.param('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name="merge_filter_mutect2." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                else:
                    jobs.append(
                        concat_jobs(
                            [
                                pipe_jobs(
                                    [
                                        bcftools.concat(
                                            inputs,
                                            None,
                                            config.param('merge_filter_mutect2', 'bcftools_options')
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
                                                    "-v 'EBV'"
                                                )
                                            ]
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_gz
                                        )
                                    ]
                                ),
                                #gatk4.filter_mutect_calls(output_gz, output_flt),
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
                                            config.param('merge_filter_mutect2', 'filter_options')
                                        ),
                                        htslib.bgzip_tabix(
                                            None,
                                            output_somatic_vt
                                        )
                                    ]
                                )
                            ],
                            name="merge_filter_mutect2." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def strelka2_paired_somatic(self):
        """
        [Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs
        This implementation is optimized for somatic calling.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            somatic_dir = os.path.join(pair_directory, "rawStrelka2_somatic")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            manta_indels = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, "rawManta", "results", "variants", "candidateSmallIndels.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(pair_directory),
                            bash.sort(
                                coverage_bed,
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
                            )
                        ],
                        name="bed_index." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                bed_file=config.param('strelka2_paired_somatic', 'bed_file')

            output_dep = [
                os.path.join(somatic_dir, "results/variants/somatic.snvs.vcf.gz"),
                os.path.join(somatic_dir, "results/variants/somatic.indels.vcf.gz")
            ]

            sed_cmd = Job(
                    [os.path.join(somatic_dir, "runWorkflow.py")],
                    [os.path.join(somatic_dir, "runWorkflow.py")],
                    command="""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {input}""".format(
    input=os.path.join(somatic_dir, "runWorkflow.py")
    )
)
            
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(somatic_dir),
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
                            output_dep=output_dep
                        )
                    ],
                    name="strelka2_paired_somatic.call."+tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor, manta_indels, bed_file],
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
                                        bash.grep(
                                            None,
                                            None,
                                            " -v 'GL00'"
                                        ),
                                        bash.grep(
                                            None,
                                            None,
                                            "-Ev 'chrUn|random'"
                                        ),
                                        bash.grep(
                                            None,
                                            None,
                                            "-v 'EBV'"
                                        )
                                    ]
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    output_prefix + ".strelka2.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.vt.vcf.gz"
                                )
                            ]
                        ),
                        tools.fix_genotypes_strelka(
                            output_prefix + ".strelka2.vt.vcf.gz",
                            output_prefix + ".strelka2.somatic.gt.vcf.gz",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name
                        ),
                        bcftools.view(
                            output_prefix + ".strelka2.somatic.gt.vcf.gz",
                            output_prefix + ".strelka2.somatic.vt.vcf.gz",
                            config.param('strelka2_paired_somatic', 'filter_options')
                        )
                    ],
                    name="strelka2_paired_somatic.filter." + tumor_pair.name,
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
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            germline_dir = os.path.join(pair_directory, "rawStrelka2_germline")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            input = [input_normal, input_tumor]

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(pair_directory),
                            bash.sort(
                                coverage_bed,
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
                            )
                        ],
                        name="bed_index." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                bed_file = config.param('strelka2_paired_germline', 'bed_file')

            output_dep = [os.path.join(germline_dir, "results", "variants", "variants.vcf.gz")]

            sed_cmd = Job(
                    [os.path.join(germline_dir, "runWorkflow.py")],
                    [os.path.join(germline_dir, "runWorkflow.py")],
                    command="""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {input}""".format(
    input=os.path.join(germline_dir, "runWorkflow.py")
    )
)
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(germline_dir),
                        strelka2.germline_config(
                            input,
                            germline_dir,
                            bed_file,
                        ),
                        sed_cmd,
                        strelka2.run(
                            germline_dir,
                            output_dep=output_dep
                        )
                    ],
                    name="strelka2_paired_germline.call."+tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=input + [bed_file],
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
                                        bash.grep(
                                            None,
                                            None,
                                            "-vE 'GL00|hs37d5'"
                                        ),
                                        bash.grep(
                                            None,
                                            None,
                                            "-Ev 'chrUn|random'"
                                        ),
                                        bash.grep(
                                            None,
                                            None,
                                            "-v 'EBV'"
                                        )
                                    ]
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.germline.vcf.gz"
                                )
                            ]
                        ),
                        pipe_jobs(
                            [
                                vt.decompose_and_normalize_mnps(
                                    output_prefix + ".strelka2.germline.vcf.gz",
                                    None
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    output_prefix + ".strelka2.germline.gt.vcf.gz"
                                )
                            ]
                        ),
                        bcftools.view(
                            output_prefix + ".strelka2.germline.gt.vcf.gz",
                            output_prefix + ".strelka2.germline.vt.vcf.gz",
                            config.param('strelka2_paired_germline', 'filter_options')
                        )
                    ],
                    name="strelka2_paired_germline.filter." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def strelka2_paired_germline_snpeff(self):
        """
        [Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small
        cohorts and somatic variation in tumor/normal sample pairs.
        This implementation is optimized for germline calling in cancer pairs.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        bcftools.split(
                            os.path.join(pair_directory, tumor_pair.name + ".strelka2.germline.vt.vcf.gz"),
                            pair_directory,
                            config.param('strelka2_paired_germline_snpeff', 'split_options'),
                        )
                    ],
                    name="strelka2_paired_germline_snpeff.split." + tumor_pair.name,
                    input_dependency=[
                        os.path.join(pair_directory, tumor_pair.name + ".strelka2.germline.vt.vcf.gz")
                    ],
                    output_dependency=[
                        os.path.join(pair_directory, tumor_pair.normal.name + ".vcf.gz"),
                        os.path.join(pair_directory, tumor_pair.tumor.name + ".vcf.gz")
                    ],
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            os.path.join(pair_directory, tumor_pair.normal.name + ".vcf.gz"),
                            os.path.join(pair_directory, tumor_pair.normal.name + ".snpeff.vcf"),
                            options=config.param('strelka2_paired_germline_snpeff', 'options')
                        ),
                        htslib.bgzip_tabix(
                            os.path.join(pair_directory, tumor_pair.normal.name + ".snpeff.vcf"),
                            os.path.join(pair_directory, tumor_pair.normal.name + ".snpeff.vcf.gz")
                        )
                    ],
                    name="strelka2_paired_germline_snpeff.normal." + tumor_pair.name,
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
                            options=config.param('strelka2_paired_germline_snpeff', 'options')
                        ),
                        htslib.bgzip_tabix(
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".snpeff.vcf"),
                            os.path.join(pair_directory, tumor_pair.tumor.name + ".snpeff.vcf.gz")
                        )
                    ],
                    name="strelka2_paired_germline_snpeff.tumor." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
        return jobs

    def vardict_paired(self):
        """
        vardict caller for SNVs and Indels.
        Note: variants are filtered to remove the instance where REF == ALT and REF is modified to 'N' when REF is
        AUPAC nomenclature.
        """

        ##TO DO - the BED system needs to be revisted !!
        jobs = []

        nb_jobs = config.param('vardict_paired', 'nb_jobs', param_type='posint')
        if nb_jobs > 50:
            log.warning("Number of vardict jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        use_bed = config.param('vardict_paired', 'use_bed', param_type='boolean', required=True)
        genome_dictionary = config.param('DEFAULT', 'genome_dictionary', param_type='filepath')

        interval_list = []

        splitjobs_dir = os.path.join(self.output_dirs['paired_variants_directory'], "splitjobs", "vardict" )
        if use_bed:
            for idx in range(nb_jobs):
                interval_list.append(os.path.join(splitjobs_dir, "exome", "interval_list", str(idx).zfill(4) + "-scattered.interval_list"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(splitjobs_dir, "exome", "interval_list"), remove=True),
                        gatk4.bed2interval_list(
                            genome_dictionary,
                            self.samples[0].readsets[0].beds[0],
                            os.path.join(splitjobs_dir, "exome", "interval_list", config.param('vardict_paired', 'assembly') + ".interval_list")
                        ),
                        gatk4.splitInterval(
                            os.path.join(splitjobs_dir, "exome", "interval_list", config.param('vardict_paired', 'assembly') + ".interval_list"),
                            os.path.join(splitjobs_dir, "exome", "interval_list"),
                            nb_jobs,
                            options="--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
                        )
                    ],
                    name="vardict_paired.create_splitjobs",
                    samples=self.samples,
                    readsets=self.readsets
                )
            )

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            if use_bed:
                idx = 0
                for interval in interval_list:
                    bed = re.sub("interval_list$", "bed", interval)
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx).zfill(4) + ".vardict.vcf.gz")
                    jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(
                                    vardict_directory,
                                    remove=True
                                ),
                                gatk4.interval_list2bed(
                                    interval,
                                    bed
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
                            name="vardict_paired." + tumor_pair.name + "." + str(idx).zfill(4),
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
                    idx += 1
            else:
                beds = []
                for idx in range(nb_jobs):
                    beds.append(os.path.join(vardict_directory, "chr." + str(idx) + ".bed"))

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                vardict_directory,
                                remove=True
                            ),
                            vardict.dict2beds(
                                genome_dictionary,
                                beds
                            )
                        ],
                        name="vardict.genome.beds." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
                for idx in range(nb_jobs):
                    output = os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz")
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
                                            beds[idx]
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
                            name="vardict_paired." + tumor_pair.name + "." + str(idx),
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
        """

        jobs = []
        nb_jobs = config.param('vardict_paired', 'nb_jobs', param_type='posint')
        use_bed = config.param('vardict_paired', 'use_bed', param_type='boolean', required=True)

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            vardict_directory = os.path.join(pair_directory, "rawVardict")
            output_tmp = os.path.join(pair_directory, tumor_pair.name + ".vardict.tmp.vcf.gz")
            output = os.path.join(pair_directory, tumor_pair.name + ".vardict.vcf.gz")
            output_vt = os.path.join(pair_directory, tumor_pair.name + ".vardict.vt.vcf.gz")
            output_somatic = os.path.join(pair_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            output_germline_loh = os.path.join(pair_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            if nb_jobs == 1 and use_bed:
                input = os.path.join(vardict_directory, tumor_pair.name + "." + str(0).zfill(4) + ".vardict.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                            bash.ln(
                                os.path.relpath(input, os.path.dirname(output_tmp)),
                                output_tmp,
                                input=input
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
                                    bash.grep(
                                        None,
                                        None,
                                        "-v 'GL00' | grep -Ev 'chrUn|random' | grep -v 'EBV'"
                                    ),
                                    bash.grep(
                                        None,
                                        None,
                                        "-Ev 'chrUn|random'"
                                    ),
                                    bash.grep(
                                        None,
                                        None,
                                        "-v 'EBV'"
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
                                        config.param('merge_filter_paired_vardict', 'somatic_filter_options')
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
                                        config.param('merge_filter_paired_vardict', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('merge_filter_paired_vardict', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_germline_loh
                                    )
                                ]
                            )
                        ],
                        name="symlink_vardict_vcf." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )
            else:
                input_vcfs = []
                for idx in range(nb_jobs):
                    input_vcfs.append(
                        os.path.join(vardict_directory, tumor_pair.name + "." + str(idx) + ".vardict.vcf.gz"))

                for input_vcf in input_vcfs:
                    if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                        log.error(f"Incomplete vardict vcf: {input_vcf}\n")

                jobs.append(
                    concat_jobs(
                        [
                            pipe_jobs(
                                [
                                    bcftools.concat(
                                        input_vcfs,
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
                                        "-v 'EBV'"
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
                                        config.param('merge_filter_paired_vardict', 'somatic_filter_options')
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
                                        config.param('merge_filter_paired_vardict', 'germline_filter_options')
                                    ),
                                    bcftools.view(
                                        None,
                                        None,
                                        config.param('merge_filter_paired_vardict', 'genotype_filter_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_germline_loh
                                    )
                                ]
                            )
                        ],
                        name="merge_filter_paired_vardict." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def ensemble_somatic(self):
        """
        Apply Bcbio.variations ensemble approach for mutect2, Vardict, Samtools and VarScan2 calls.
        Filter ensemble calls to retain only calls overlapping 2 or more callers.
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_mutect2 = os.path.join(input_directory, tumor_pair.name + ".mutect2.somatic.vt.vcf.gz")
            input_strelka2 = os.path.join(input_directory, tumor_pair.name + ".strelka2.somatic.purple.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.somatic.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.somatic.vt.vcf.gz")
            inputs_somatic = [input_mutect2, input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_somatic:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete ensemble vcf: {input_vcf}\n")

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            jobs.append(
                concat_jobs(
                    [
                        # Create output directory since it is not done by default by GATK tools
                        bash.mkdir(
                            paired_ensemble_directory
                        ),
                        # Remove any existing outputs because they cause silent error
                        bash.rm(
                            output_ensemble
                            ),
                        bash.rm(
                            output_ensemble + ".tbi"
                            ),
                        bash.rm(
                            os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.somatic.vt-work")
                            ),
                        bcbio_variation_recall.ensemble(
                            inputs_somatic,
                            output_ensemble,
                            config.param('bcbio_ensemble_somatic', 'options')
                        )
                    ],
                    name="bcbio_ensemble_somatic." + tumor_pair.name,
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
        """

        jobs = []
        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_ensemble_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            input_strelka2 = os.path.join(input_directory, tumor_pair.name + ".strelka2.germline.vt.vcf.gz")
            input_vardict = os.path.join(input_directory, tumor_pair.name + ".vardict.germline.vt.vcf.gz")
            input_varscan2 = os.path.join(input_directory, tumor_pair.name + ".varscan2.germline.vt.vcf.gz")

            inputs_germline = [input_strelka2, input_vardict, input_varscan2]

            for input_vcf in inputs_germline:
                if not self.is_gz_file(os.path.join(self.output_dir, input_vcf)):
                    log.error(f"Incomplete ensemble vcf: {input_vcf}\n")

            output_ensemble = os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt.vcf.gz")

            # if os.path.isdir(os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt-work")):
            #     rm_job = bash.rm(
            #         os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt-work")
            #     )
            #     jobs.append(rm_job)

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
                            os.path.join(paired_ensemble_directory, tumor_pair.name + ".ensemble.germline.vt-work")
                            ),
                        bcbio_variation_recall.ensemble(
                            inputs_germline,
                            output_ensemble,
                            config.param('bcbio_ensemble_germline', 'options')
                        )
                    ],
                    name="bcbio_ensemble_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=inputs_germline
                )
            )

        return jobs

    def gatk_variant_annotator_somatic(self):
        """
        Add vcf annotations to ensemble vcf: Standard and Somatic annotations.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')
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
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")]
                ]
            )
            input_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.vcf.gz")

            if nb_jobs == 1:
                output_somatic_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")

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
                                config.param('gatk_variant_annotator_somatic', 'other_options')
                            )
                        ],
                        name="gatk_variant_annotator_somatic." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) + ".vcf.gz")

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
                                    config.param('gatk_variant_annotator_somatic', 'other_options'),
                                    intervals=sequences
                                )
                            ],
                            name="gatk_variant_annotator_somatic." + str(idx) + "." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                output_somatic_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz")

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
                                config.param('gatk_variant_annotator_somatic', 'other_options'),
                                exclude_intervals=unique_sequences_per_job_others
                            )
                        ],
                        name="gatk_variant_annotator_somatic.others." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def gatk_variant_annotator_germline(self):
        """
        Add vcf annotations to ensemble vcf: most importantly the AD field.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')
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
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")]
                ]
            )
            input_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.vcf.gz")

            if nb_jobs == 1:
                output_germline_variants = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")

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
                                config.param('gatk_variant_annotator_germline', 'other_options'),
                            )
                        ],
                        name="gatk_variant_annotator_germline." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            else:
                unique_sequences_per_job, unique_sequences_per_job_others = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                for idx, sequences in enumerate(unique_sequences_per_job):
                    output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")

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
                                    config.param('gatk_variant_annotator_germline', 'other_options'),
                                    intervals=sequences
                                )
                            ],
                            name="gatk_variant_annotator_germline." + str(idx) + "." + tumor_pair.name,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

                output_germline_variants = os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz")

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
                                config.param('gatk_variant_annotator_germline', 'other_options'),
                                exclude_intervals=unique_sequences_per_job_others
                            )
                        ],
                        name="gatk_variant_annotator_germline.others." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def merge_gatk_variant_annotator_somatic(self):
        """
        Merge annotated somatic vcfs.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_somatic = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            if nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot." + str(idx) +".vcf.gz")
                                  for idx in range(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.others.vcf.gz"))

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
                        name="merge_gatk_variant_annotator.somatic." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def merge_gatk_variant_annotator_germline(self):
        """
        Merge annotated germline and LOH vcfs.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        nb_jobs = config.param('gatk_variant_annotator', 'nb_jobs', param_type='posint')

        for tumor_pair in self.tumor_pairs.values():
            annot_directory = os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation")
            output_germline = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")

            if nb_jobs > 1:
                unique_sequences_per_job, _ = sequence_dictionary.split_by_size(self.sequence_dictionary_variant(), nb_jobs - 1, variant=True)
                vcfs_to_merge = [os.path.join(ensemble_directory, tumor_pair.name, "rawAnnotation", tumor_pair.name + ".ensemble.germline.vt.annot." + str(idx) + ".vcf.gz")
                                 for idx in range(len(unique_sequences_per_job))]

                vcfs_to_merge.append(os.path.join(annot_directory, tumor_pair.name + ".ensemble.germline.vt.annot.others.vcf.gz"))

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
                        name="merge_gatk_variant_annotator.germline." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def filter_ensemble_germline(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz"
            )
            output_2caller = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.vcf.gz"
            )
            output_filter = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.flt.vcf.gz"
            )

            jobs.append(
                concat_jobs(
                    [
                        tools.format2pcgr(
                            input,
                            output_2caller,
                            config.param('filter_ensemble', 'call_filter'),
                            "germline",
                            tumor_pair.tumor.name,
                            ini_section='filter_ensemble'
                        ),
                        pipe_jobs(
                            [
                                bcftools.view(
                                    output_2caller,
                                    None,
                                    filter_options=config.param('filter_ensemble', 'germline_filter_options'),
                                ),
                                bcftools.view(
                                    None,
                                    None,
                                    filter_options="-Oz -s ^" + tumor_pair.normal.name
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
                    name="filter_ensemble.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def report_cpsr(self):
        """
        Creates a CPSR germline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble germline vcf
        output: html report and addtional flat files
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.germline.vt.annot.2caller.flt.vcf.gz"
            )
            cpsr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "cpsr"
            )

            job_name = f"report_cpsr.{tumor_pair.name}"
            cpsr_job = concat_jobs(
                    [
                        bash.mkdir(
                            cpsr_directory,
                        ),
                        cpsr.report(
                            input,
                            cpsr_directory,
                            tumor_pair.name
                        )
                    ],
                    name=job_name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )

            if self.project_tracking_json:
                samples = [tumor_pair.normal, tumor_pair.tumor]
                cpsr_output_file = os.path.join(self.output_dir, "job_output", "report_cpsr", f"{job_name}_{self.timestamp}.o")
                jobs.append(
                    concat_jobs(
                        [
                            cpsr_job,
                            cpsr.parse_cpsr_passed_variants_pt(cpsr_output_file),
                            job2json_project_tracking.run(
                                input_file=cpsr_output_file,
                                pipeline=self,
                                samples=",".join([sample.name for sample in samples]),
                                readsets=",".join([readset.name for sample in samples for readset in sample.readsets]),
                                job_name=job_name,
                                metrics="cpsr_passed_variants=$cpsr_passed_variants"
                            )
                        ],
                        name=job_name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                        input_dependency=[input]
                    )
                )
            else:
                jobs.append(cpsr_job)

        return jobs

    def filter_ensemble_somatic(self):
        """
        Applies custom script to inject FORMAT information - tumor/normal DP and VAP into the INFO field
        the filter on those generated fields.
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz"
            )
            output_2caller = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.vcf.gz"
            )
            output_filter = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
            )

            jobs.append(
                concat_jobs(
                    [
                        tools.format2pcgr(
                            input,
                            output_2caller,
                            config.param('filter_ensemble', 'call_filter'),
                            "somatic",
                            tumor_pair.tumor.name,
                            ini_section='filter_ensemble'
                        ),
                        bcftools.view(
                            output_2caller,
                            output_filter,
                            filter_options=config.param('filter_ensemble', 'somatic_filter_options'),
                        ),
                        htslib.tabix(
                            output_filter,
                            options="-pvcf"
                        )
                    ],
                    name="filter_ensemble.somatic." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )

        return jobs

    def report_pcgr(self):
        """
        Creates a PCGR somatic + germline report (https://sigven.github.io/cpsr/)
        input: filtered ensemble germline vcf
        output: html report and addtionalflat files
        """
        jobs = []

        ensemble_directory = os.path.join(
            self.output_dirs['paired_variants_directory'],
            "ensemble"
        )
        assembly = config.param('report_pcgr', 'assembly')

        for tumor_pair in self.tumor_pairs.values():
            cpsr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "cpsr"
            )
            input_cpsr = os.path.join(
                cpsr_directory,
                tumor_pair.name + ".cpsr." + assembly + ".json.gz"
            )
            input = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.flt.vcf.gz"
            )
            input_cna = os.path.join(
                self.output_dirs['sv_variants_directory'],
                tumor_pair.name,
                tumor_pair.name + ".cnvkit.vcf.gz"
            )
            header = os.path.join(
                self.output_dirs['sv_variants_directory'],
                tumor_pair.name + ".header"
            )
            output_cna_body = os.path.join(
                self.output_dirs['sv_variants_directory'],
                tumor_pair.name + ".cnvkit.body.tsv"
            )
            output_cna = os.path.join(
                self.output_dirs['sv_variants_directory'],
                tumor_pair.name + ".cnvkit.cna.tsv"
            )
            pcgr_directory = os.path.join(
                ensemble_directory,
                tumor_pair.name,
                "pcgr"
            )

            # PCGR does not accept sample IDs longer than 35 characters and uses the sample ID to name output files.
            # For samples that have longer sample IDs the output files will have non-matching names, so create symlinks with full-length names.
            if tumor_pair.name != tumor_pair.name[:35]:
                output = [
                        os.path.join(pcgr_directory, tumor_pair.name[:35] + ".pcgr_acmg." + assembly + ".flexdb.html"),
                        os.path.join(pcgr_directory, tumor_pair.name[:35] + ".pcgr_acmg." + assembly + ".maf"),
                        os.path.join(pcgr_directory, tumor_pair.name[:35] + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv"),
                        os.path.join(pcgr_directory, tumor_pair.name[:35] + ".pcgr_acmg." + assembly + ".cna_segments.tsv.gz")
                    ]
                final_command = concat_jobs(
                        [
                            bash.ln(
                                os.path.relpath(output[0], pcgr_directory),
                                os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".flexdb.html"),
                                input=output[0]
                                ),
                            bash.ln(
                                os.path.relpath(output[1], pcgr_directory),
                                os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".maf"),
                                input=output[1]
                                ),
                            bash.ln(
                                os.path.relpath(output[2], pcgr_directory),
                                os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv"),
                                input=output[2]
                                ),
                            bash.ln(
                                os.path.relpath(output[3], pcgr_directory),
                                os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".cna_segments.tsv.gz"),
                                input=output[3]
                                )
                            ]
                        )
            else:
                output = [
                        os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".flexdb.html"),
                        os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".maf"),
                        os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".snvs_indels.tiers.tsv"),
                        os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".cna_segments.tsv.gz")
                    ]
                final_command = bash.ls(output[0])
            job_name = f"report_pcgr.{tumor_pair.name}"
            pcgr_job = concat_jobs(
                [
                    bash.mkdir(
                        pcgr_directory,
                    ),
                    pcgr.create_header(
                        header,
                    ),
                    bcftools.query(
                        input_cna,
                        output_cna_body,
                        query_options="-f '%CHROM\\t%POS\\t%END\\t%FOLD_CHANGE_LOG\\n'"
                    ),
                    bash.cat(
                        [
                            header,
                            output_cna_body,
                        ],
                        output_cna
                    ),
                    pcgr.report(
                        input,
                        input_cpsr,
                        pcgr_directory,
                        tumor_pair.name,
                        input_cna=output_cna
                    ),
                    final_command
                ],
                name=job_name,
                samples=[tumor_pair.normal, tumor_pair.tumor],
                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                input_dependency=[header, input, input_cna, input_cpsr, output_cna_body],
                output_dependency=[header, output_cna_body, output_cna] + output
            )

            if self.project_tracking_json:
                samples = [tumor_pair.normal, tumor_pair.tumor]
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
                        input_dependency=[header, input, input_cna, input_cpsr, output_cna_body],
                        output_dependency=[header, output_cna_body, output_cna] + output
                    )
                )
            else:
                jobs.append(pcgr_job)

        return jobs

    def compute_cancer_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        if not os.path.exists(ensemble_directory):
            os.makedirs(ensemble_directory)

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            if not os.path.exists(paired_directory):
                os.makedirs(paired_directory)

            input_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz")
            output_somatic = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            with open(cancer_pair_filename, 'w') as cancer_pair:
                cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            paired_directory,
                            remove=True
                        ),
                        snpeff.compute_effects(
                            input_somatic,
                            output_somatic,
                            cancer_sample_file=cancer_pair_filename,
                            options=config.param('compute_cancer_effects_somatic', 'options')
                        ),
                        htslib.bgzip_tabix(
                            output_somatic,
                            output_somatic + ".gz"
                        )
                    ],
                    name="compute_cancer_effects_somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def compute_cancer_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)

            input_germline = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz")
            output_germline = os.path.join(paired_directory,
                                           tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.vcf")

            cancer_pair_filename = os.path.join(paired_directory, tumor_pair.name + '.tsv')
            with open(cancer_pair_filename, 'w') as cancer_pair:
                cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            paired_directory,
                            remove=True
                        ),
                        snpeff.compute_effects(
                            input_germline,
                            output_germline,
                            options=config.param('compute_cancer_effects_germline', 'options')
                        ),
                        htslib.bgzip_tabix(
                            output_germline,
                            output_germline + ".gz"
                        )
                    ],
                    name="compute_cancer_effects_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def ensemble_somatic_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """

        jobs = []


        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.somatic.vt.annot.snpeff.dbnsfp.vcf")

            jobs.append(
                concat_jobs(
                    [
                        snpeff.snpsift_dbnsfp(
                            input_vcf,
                            output_vcf
                        ),
                        htslib.bgzip_tabix(
                            output_vcf,
                            output_vcf + ".gz"
                        )
                    ],
                    name="dbnsfp_annotation.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
        return jobs

    def ensemble_germline_dbnsfp_annotation(self):
        """
        Additional SVN annotations. Provides extra information about SVN by using numerous published databases.
        Applicable to human samples. Databases available include Biomart (adds GO annotations based on gene information)
        and dbNSFP (an integrated database of functional annotations from multiple sources for the comprehensive
        collection of human non-synonymous SNPs. It compiles prediction scores from four prediction algorithms
        (SIFT, Polyphen2, LRT and MutationTaster), three conservation scores (PhyloP, GERP++ and SiPhy)
        and other function annotations).
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            input_vcf = os.path.join(paired_directory, tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.vcf.gz")
            output_vcf = os.path.join(paired_directory,
                                      tumor_pair.name + ".ensemble.germline.vt.annot.snpeff.dbnsfp.vcf")

            jobs.append(
                concat_jobs(
                    [
                        snpeff.snpsift_dbnsfp(
                            input_vcf,
                            output_vcf
                        ),
                        htslib.bgzip_tabix(
                            output_vcf,
                            output_vcf + ".gz"
                        )
                    ],
                    name="dbnsfp_annotation.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sample_gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database: [Gemini](http://gemini.readthedocs.org/en/latest/index.html)
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)


            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            paired_directory,
                            remove=True
                        ),
                        gemini.gemini_annotations(
                            gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                            gemini_prefix + ".somatic.gemini." + gemini_version + ".db",
                            self.output_dir
                        )
                    ],
                    name="gemini_annotations.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sample_gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database: [Gemini](http://gemini.readthedocs.org/en/latest/index.html)
        """
        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        gemini_module = config.param("DEFAULT", 'module_gemini').split(".")
        gemini_version = ".".join([gemini_module[-2], gemini_module[-1]])

        for tumor_pair in self.tumor_pairs.values():
            paired_directory = os.path.join(ensemble_directory, tumor_pair.name)
            gemini_prefix = os.path.join(paired_directory, tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            paired_directory,
                            remove=True
                        ),
                        gemini.gemini_annotations(
                            gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                            gemini_prefix + ".germline.gemini." + gemini_version + ".db",
                            self.output_dir
                        )
                    ],
                    name="gemini_annotations.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sym_link_ensemble(self):
        """
        Create sym link of ensemble output for delivery of data to clients.
        """
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            inputs["Tumor"] =  [os.path.join(self.output_dirs["paired_variants_directory"], "ensemble", tumor_pair.name, tumor_pair.name)]

            for key, input_files in inputs.items():
                for idx, sample_prefix in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz",
                                    sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.somatic.vt.annot.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.md5sum(
                                    sample_prefix + ".ensemble.germline.vt.annot.vcf.gz",
                                    sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.germline.vt.annot.vcf.gz",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    sample_prefix + ".ensemble.germline.vt.annot.vcf.gz.tbi",
                                    tumor_pair,
                                    self.output_dir,
                                    type="snv/ensemble",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_ensemble." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def combine_tumor_pairs_somatic(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.vcf.gz") for tumor_pair in self.tumor_pairs.values()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        Job(
                            [input_merged_vcfs[0]],
                            [output],
                            command="ln -s -f " + os.path.relpath(input_merged_vcfs[0], os.path.dirname(output)) + " " + output
                        )
                    ],
                    name="gatk_combine_variants.somatic.allPairs",
                    samples=self.samples,
                    readsets=self.readsets
                )
            )

        else:

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        gatk.combine_variants(
                            input_merged_vcfs,
                            output
                        )
                    ],
                    name="gatk_combine_variants.somatic.allPairs",
                    samples=self.samples,
                    readsets=self.readsets
                )
            )

        return jobs

    def combine_tumor_pairs_germline(self):
        """
        Combine numerous ensemble vcfs into one vcf for gemini annotations.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_merged_vcfs = [os.path.join(ensemble_directory, tumor_pair.name,
                                          tumor_pair.name + ".ensemble.germline.vt.annot.vcf.gz") for tumor_pair in
                             self.tumor_pairs.values()]
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        if len(input_merged_vcfs) == 1:
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        bash.ln(
                            os.path.relpath(input_merged_vcfs[0], os.path.dirname(output)),
                            output,
                            intput=input_merged_vcfs[0]
                        )
                    ],
                    name="gatk_combine_variants.germline.allPairs",
                    samples=sample_list,
                    readsets=[sample.readsets for sample in sample_list]
                )
            )

        else:
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            ensemble_directory,
                            remove=True
                        ),
                        gatk.combine_variants(
                            input_merged_vcfs,
                            output
                        )
                    ],
                    name="gatk_combine_variants.germline.allPairs",
                    samples=sample_list,
                    readsets=[sample.readsets for sample in sample_list]
                )
            )

        return jobs

    def decompose_and_normalize_mnps_somatic(self):
        """
        Processes include normalization and decomposition of MNPs by [vt](http://genome.sph.umich.edu/wiki/Vt).
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    vt.decompose_and_normalize_mnps(
                        input,
                        output
                    )
                ],
                name="decompose_and_normalize_mnps.somatic.allPairs",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def decompose_and_normalize_mnps_germline(self):
        """
        Processes include normalization and decomposition of MNPs by [vt](http://genome.sph.umich.edu/wiki/Vt).
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.annot.vcf.gz")
        output_vcf = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")

        job = vt.decompose_and_normalize_mnps(input_vcf, output_vcf)
        job.name = "decompose_and_normalize_mnps.germline.allPairs"

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    vt.decompose_and_normalize_mnps(
                        input_vcf,
                        output_vcf
                    )
                ],
                name="decompose_and_normalize_mnps.somatic.allPairs",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def all_pairs_compute_effects_somatic(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.somatic.vt.annot.snpeff.vcf.gz")

        cancer_pair_filename = os.path.join('cancer_snpeff.tsv')
        with open(cancer_pair_filename, 'w') as cancer_pair:
            sample_list = []
            for tumor_pair in self.tumor_pairs.values():
                cancer_pair.write(tumor_pair.normal.name + "\t" + tumor_pair.tumor.name + "\n")
                sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    snpeff.compute_effects(
                        input,
                        output,
                        cancer_sample_file=cancer_pair_filename,
                        options=config.param('compute_cancer_effects_somatic', 'options')
                    ),
                    htslib.bgzip_tabix(
                        output,
                        output_gz
                    )
                ],
                name="compute_effects.somatic.allPairs",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def all_pairs_compute_effects_germline(self):
        """
        Variant effect annotation. The .vcf files are annotated for variant effects using the SnpEff software.
        SnpEff annotates and predicts the effects of variants on genes (such as amino acid changes).
        Modified arguments to consider paired cancer data.
        Applied to all tumor pairs.
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        input = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.vcf.gz")
        output = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf")
        output_gz = os.path.join(ensemble_directory, "allPairs.ensemble.germline.vt.annot.snpeff.vcf.gz")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    snpeff.compute_effects(
                        input,
                        output,
                        options=config.param('compute_cancer_effects_germline', 'options')
                    ),
                    htslib.bgzip_tabix(
                        output,
                        output_gz
                    )
                ],
                name="compute_effects.germline.allPair",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def gemini_annotations_somatic(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database: [Gemini](http://gemini.readthedocs.org/en/latest/index.html).
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    gemini.gemini_annotations(
                        gemini_prefix + ".ensemble.somatic.vt.annot.snpeff.vcf.gz",
                        gemini_prefix + ".somatic.gemini.db",
                        temp_dir
                    )
                ],
                name="gemini_annotations.somatic.allPairs",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def gemini_annotations_germline(self):
        """
        Load functionally annotated vcf file into a mysql lite annotation database: [Gemini](http://gemini.readthedocs.org/en/latest/index.html).
        """

        jobs = []

        ensemble_directory = os.path.join(self.output_dirs['paired_variants_directory'], "ensemble")
        temp_dir = os.path.join(os.getcwd(), ensemble_directory)
        gemini_prefix = os.path.join(ensemble_directory, "allPairs")

        sample_list = []
        for tumor_pair in self.tumor_pairs.values():
            sample_list.extend([tumor_pair.normal, tumor_pair.tumor])

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(
                        ensemble_directory,
                        remove=True
                    ),
                    gemini.gemini_annotations(
                        gemini_prefix + ".ensemble.germline.vt.annot.snpeff.vcf.gz",
                        gemini_prefix + ".germline.gemini.db",
                        temp_dir
                    )
                ],
                name="gemini_annotations.germline.allPairs",
                samples=sample_list,
                readsets=[sample.readsets for sample in sample_list]
            )
        )

        return jobs

    def sequenza(self):
        """
        Sequenza is a novel set of tools providing a fast Python script to genotype cancer samples,
        and an R package to estimate cancer cellularity, ploidy, genome-wide copy number profile and infer
        for mutated alleles.
        """
        jobs = []
        nb_jobs = config.param('sequenza', 'nb_jobs', param_type='posint')
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
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            raw_output = os.path.join(sequenza_directory, "rawSequenza", tumor_pair.name + ".")
            output = os.path.join(sequenza_directory, tumor_pair.name + ".")

            if nb_jobs == 1:
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                raw_sequenza_directory,
                                remove=True
                            ),
                            sequenza.bam2seqz(
                                input_normal,
                                input_tumor,
                                config.param('sequenza', 'gc_file'),
                                raw_output + "all.seqz.gz",
                                None
                            ),
                            sequenza.bin(
                                raw_output + "all.seqz.gz",
                                output + "all.binned.seqz.gz",
                            )
                        ],
                        name="sequenza.create_seqz." + tumor_pair.name,
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
                                output + "all.binned.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name
                            )
                        ],
                        name="sequenza." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
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
                                        remove=True
                                    ),
                                    sequenza.bam2seqz(
                                        input_normal,
                                        input_tumor,
                                        config.param('sequenza', 'gc_file'),
                                        raw_output + "seqz." + sequence['name'] + ".gz",
                                        sequence['name']
                                    ),
                                    sequenza.bin(
                                        raw_output + "seqz." + sequence['name'] + ".gz",
                                        raw_output + "binned.seqz." + sequence['name'] + ".gz",
                                    )
                                ],
                                name="sequenza.create_seqz." + sequence['name'] + "." + tumor_pair.name,
                                samples=[tumor_pair.normal, tumor_pair.tumor],
                                readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                            )
                        )

                seqz_outputs = [raw_output + "binned.seqz." + sequence['name'] + ".gz" for sequence in self.sequence_dictionary_variant() if sequence['type'] == 'primary']

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
                                        output + "binned.merged.seqz.gz",
                                        options="-cf"
                                    )
                                ]
                            )
                        ],
                        name="sequenza.merge_binned_seqz." + tumor_pair.name,
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
                                output + "binned.merged.seqz.gz",
                                sequenza_directory,
                                tumor_pair.name
                            )
                        ],
                        name="sequenza." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def sym_link_sequenza(self):
        """
        Sym link of sequenza outputs.
        """
        jobs = []

        inputs = {}

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)
            print(pair_directory)
            inputs["Tumor"] = [
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_chromosome_view.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_genome_view.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CN_bars.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_CP_contours.pdf"),
                os.path.join(pair_directory, "sequenza", tumor_pair.name + "_ploidy_celularity.tsv")
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
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link.sequenza." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        return jobs

    def purple(self, sv=False):
        """
        PURPLE is a purity ploidy estimator for whole genome sequenced (WGS) data.

        It combines B-allele frequency (BAF) from AMBER, read depth ratios from COBALT,
        somatic variants and structural variants to estimate the purity and copy number profile of a tumor sample.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_dir = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )

            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            somatic_snv = None
            if os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz"):
                somatic_snv = os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.purple.vcf.gz")
                jobs.append(
                    concat_jobs(
                        [
                        purple.strelka2_convert(
                            os.path.join(pair_dir, tumor_pair.name + ".strelka2.somatic.vt.vcf.gz"),
                            somatic_snv,
                        )
                    ],
                    name="purple.convert_strelka2." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            gripss_vcf = None
            gripss_filtered_vcf = None
            somatic_hotspots = None
            germline_hotspots = None
            driver_gene_panel = None

            if sv:
                pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
                gridss_directory = os.path.join(pair_dir, "gridss")
                gripss_vcf = os.path.join(gridss_directory, tumor_pair.tumor.name + ".gripss.somatic.vcf.gz")
                gripss_filtered_vcf = os.path.join(gridss_directory, tumor_pair.tumor.name + ".gripss.filtered.somatic.vcf.gz")
                somatic_hotspots = config.param('purple', 'somatic_hotspots', param_type='filepath')
                germline_hotspots = config.param('purple', 'germline_hotspots', param_type='filepath')
                driver_gene_panel = config.param('purple', 'driver_gene_panel', param_type='filepath')

            purple_dir = os.path.join(pair_dir, "purple")
            amber_dir = os.path.join(purple_dir, "rawAmber")
            cobalt_dir = os.path.join(purple_dir, "rawCobalt")
            ensembl_data_dir = config.param('purple', 'ensembl_data_dir', param_type='dirpath')

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            amber_dir
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
                            cobalt_dir
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
            purple_purity_output = os.path.join(purple_dir, tumor_pair.tumor.name + ".purple.purity.tsv")
            purple_qc_output = os.path.join(purple_dir, tumor_pair.tumor.name + ".purple.qc")
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
                            somatic_snv,
                            gripss_vcf,
                            gripss_filtered_vcf,
                            somatic_hotspots,
                            germline_hotspots,
                            driver_gene_panel
                        ),
                        bash.mkdir(
                            self.output_dirs['report'][tumor_pair.name]
                        ),
                        bash.ln(
                            os.path.relpath(purple_purity_output, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(purple_purity_output)),
                            input=purple_purity_output
                        ),
                        bash.ln(
                            os.path.relpath(purple_qc_output, self.output_dirs['report'][tumor_pair.name]),
                            os.path.join(self.output_dirs['report'][tumor_pair.name], os.path.basename(purple_qc_output)),
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
                    purple_purity_output,
                    purple_qc_output
                ]
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
        Returns: bcf file
        """

        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")

            filename = os.path.join(delly_directory, tumor_pair.name + '.tsv')
            if not os.path.exists(os.path.dirname(filename)):
                os.makedirs(os.path.dirname(filename))

            with open(filename, 'w') as cancer_pair:
                cancer_pair.write(tumor_pair.tumor.name + "\ttumor\n")
                cancer_pair.write(tumor_pair.normal.name + "\tcontrol\n")

            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")

            inputs = [input_tumor, input_normal]

            sv_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            for sv_type in sv_types:
                output_bcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf")
                output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".somatic.flt.vcf.gz")

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                delly_directory,
                                remove=True
                            ),
                            delly.call(
                                inputs,
                                output_bcf,
                                sv_type
                            ),
                            pipe_jobs(
                                [
                                    bcftools.view(
                                        output_bcf,
                                        None,
                                        config.param('delly_call_filter_somatic', 'bcftools_options')
                                    ),
                                    htslib.bgzip_tabix(
                                        None,
                                        output_vcf
                                    )
                                ]
                            )
                        ],
                        name="delly_call_filter." + str(sv_type) + "." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

        return jobs

    def delly_sv_annotation(self):
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            final_directory = os.path.join(self.output_dirs["sv_variants_directory"], tumor_pair.name, tumor_pair.name)
            delly_directory = os.path.join(pair_directory, "rawDelly")
            output_vcf = os.path.join(delly_directory, tumor_pair.name + ".delly.merge.sort.vcf.gz")
            output_flt_vcf = os.path.join(pair_directory, tumor_pair.name + ".delly.merge.sort.flt.vcf.gz")

            sv_types = config.param('delly_call_filter', 'sv_types_options').split(",")

            input_bcf = []
            for sv_type in sv_types:
                input_bcf.append(os.path.join(delly_directory, tumor_pair.name + ".delly." + str(sv_type) + ".bcf"))

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
                                bcftools.view(
                                    output_vcf,
                                    None,
                                    "-f PASS"
                                ),
                                htslib.bgzip(
                                    None,
                                    output_flt_vcf
                                )
                            ]
                        )
                    ],
                    name="sv_annotation.delly.merge_sort_filter." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        vawk.paired_somatic(
                            output_flt_vcf,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            final_directory + ".delly.somatic.vcf"
                        ),
                        htslib.bgzip(
                            final_directory + ".delly.somatic.vcf",
                            final_directory + ".delly.somatic.vcf.gz"
                        ),
                        snpeff.compute_effects(
                            final_directory + ".delly.somatic.vcf",
                            final_directory + ".delly.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            final_directory + ".delly.somatic.snpeff.vcf",
                            final_directory + ".delly.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            final_directory + ".delly.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "DELLY",
                            final_directory + ".delly.somatic.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.delly.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[output_flt_vcf]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        vawk.paired_germline(
                            output_flt_vcf,
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            final_directory + ".delly.germline.vcf"
                        ),
                        htslib.bgzip(
                            final_directory + ".delly.germline.vcf",
                            final_directory + ".delly.germline.vcf.gz"
                        ),
                        snpeff.compute_effects(
                            final_directory + ".delly.germline.vcf",
                            final_directory + ".delly.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            final_directory + ".delly.germline.snpeff.vcf",
                            final_directory + ".delly.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            final_directory + ".delly.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "DELLY",
                            final_directory + ".delly.germline.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.delly.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sym_link_delly(self):
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".delly.somatic.snpeff.annot.vcf",
                pair_directory + ".delly.somatic.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_delly.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".delly.germline.snpeff.annot.vcf",
                pair_directory + ".delly.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_delly.germline." + str(idx) + "." + tumor_pair.name + "." + key,
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
        Returns: Manta accepts input read mappings from BAM or CRAM files and reports all SV and indel inferences
        in VCF 4.1 format.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            manta_directory = os.path.join(pair_directory, "rawManta")
            output_prefix = os.path.join(pair_directory, tumor_pair.name)

            [input_normal] = self.select_input_files(
                [
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")]
                ]
            )
            manta_somatic_output = os.path.join(manta_directory, "results/variants/somaticSV.vcf.gz")
            manta_germline_output = os.path.join(manta_directory, "results/variants/diploidSV.vcf.gz")

            bed_file = None
            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            if coverage_bed:
                local_coverage_bed = os.path.join(pair_directory, os.path.basename(coverage_bed))
                bed_file = local_coverage_bed + ".gz"
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(manta_directory),
                            bash.sort(
                                coverage_bed,
                                local_coverage_bed + ".sort",
                                "-V -k1,1 -k2,2n -k3,3n"
                            ),
                            htslib.bgzip(
                                local_coverage_bed + ".sort",
                                bed_file
                            ),
                            htslib.tabix(
                                bed_file,
                                "-f -p bed"
                            )
                        ],
                        name="bed_index." + tumor_pair.name,
                        samples=[tumor_pair.normal, tumor_pair.tumor],
                        readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                    )
                )

            output_dep = [
                manta_somatic_output,
                manta_somatic_output + ".tbi",
                manta_germline_output,
                manta_germline_output + ".tbi"
            ]

            sed_cmd = Job(
                    [os.path.join(manta_directory, "runWorkflow.py")],
                    [os.path.join(manta_directory, "runWorkflow.py")],
                    command="""\
sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g {input}""".format(
    input=os.path.join(manta_directory, "runWorkflow.py")
    )
)
            jobs.append(
                concat_jobs(
                    [
                        bash.rm(manta_directory),
                        bash.mkdir(
                            manta_directory,
                            remove=True
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
                            output_prefix + ".manta.somatic.vcf.gz",
                            input=manta_somatic_output,
                        ),
                        bash.ln(
                            os.path.relpath(manta_somatic_output + ".tbi", os.path.dirname(output_prefix)),
                            output_prefix + ".manta.somatic.vcf.gz.tbi",
                            input=manta_somatic_output + ".tbi"
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output, os.path.dirname(output_prefix)),
                            output_prefix + ".manta.germline.vcf.gz",
                            input=manta_germline_output
                        ),
                        bash.ln(
                            os.path.relpath(manta_germline_output + ".tbi", os.path.dirname(output_prefix)),
                            output_prefix + ".manta.germline.vcf.gz.tbi",
                            input=manta_germline_output + ".tbi"
                        )
                    ],
                    name="manta_sv." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor, bed_file]
                )
            )

        return jobs

    def manta_sv_annotation(self):

        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            pair_directory + ".manta.somatic.vcf.gz",
                            pair_directory + ".manta.somatic.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".manta.somatic.snpeff.vcf",
                            pair_directory + ".manta.somatic.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            pair_directory + ".manta.somatic.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "MANTA",
                            pair_directory + ".manta.somatic.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.manta_somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

            jobs.append(
                concat_jobs(
                    [
                        snpeff.compute_effects(
                            pair_directory + ".manta.germline.vcf.gz",
                            pair_directory + ".manta.germline.snpeff.vcf"
                        ),
                        annotations.structural_variants(
                            pair_directory + ".manta.germline.snpeff.vcf",
                            pair_directory + ".manta.germline.snpeff.annot.vcf"
                        ),
                        vawk.sv(
                            pair_directory + ".manta.germline.snpeff.annot.vcf",
                            tumor_pair.normal.name,
                            tumor_pair.tumor.name,
                            "MANTA",
                            pair_directory + ".manta.germline.prioritize.tsv"
                        )
                    ],
                    name="sv_annotation.manta_germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )

        return jobs

    def sym_link_manta(self):
        jobs = []

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".manta.somatic.snpeff.annot.vcf",
                pair_directory + ".manta.somatic.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_manta.somatic." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )

        inputs = {}
        for tumor_pair in self.tumor_pairs.values():
            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name, tumor_pair.name)
            inputs["Tumor"] = [
                pair_directory + ".manta.germline.snpeff.annot.vcf",
                pair_directory + ".manta.germline.prioritize.tsv"
            ]

            for key, input_files in inputs.items():
                for idx, input_file in enumerate(input_files):
                    jobs.append(
                        concat_jobs(
                            [
                                deliverables.md5sum(
                                    input_file,
                                    input_file + ".md5",
                                    self.output_dir
                                ),
                                deliverables.sym_link_pair(
                                    input_file + ".md5",
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                ),
                                deliverables.sym_link_pair(
                                    input_file,
                                    tumor_pair,
                                    self.output_dir,
                                    type="sv",
                                    sample=key,
                                    profyle=self.args.profyle
                                )
                            ],
                            name="sym_link_manta.germline." + str(idx) + "." + tumor_pair.name + "." + key,
                            samples=[tumor_pair.normal, tumor_pair.tumor],
                            readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                        )
                    )
        return jobs

    def cnvkit_batch(self):
        """
        CNVkit is a Python library and command-line software toolkit to infer and visualize copy number from
        high-throughput DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and
        custom target panels, and short-read sequencing platforms such as Illumina and Ion Torrent.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)

            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            cnvkit_dir = os.path.join(pair_directory, "rawCNVkit")
            input_normal = os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")
            input_tumor = os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")
            tarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.targetcoverage.cnn")
            antitarcov_cnn = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".sorted.dup.antitargetcoverage.cnn")
            ref_cnn = os.path.join(cnvkit_dir, tumor_pair.name + ".reference.cnn")
            tumor_cns = os.path.join(cnvkit_dir, tumor_pair.tumor.name + ".cns")
            vcf_gz = os.path.join(pair_directory, tumor_pair.name + ".cnvkit.vcf.gz")

            metrics = os.path.join(self.output_dirs['sv_variants_directory'], "cnvkit_reference")
            pool_ref = os.path.join(self.output_dir, metrics, "pooledReference.cnn")

            if os.path.isfile(pool_ref):
                pool_ref_cnn = pool_ref
                ref_cnn = None

            else:
                pool_ref_cnn = None

            coverage_bed = bvatools.resolve_readset_coverage_bed(
                tumor_pair.normal.readsets[0]
            )

            bed = None

            if coverage_bed:
                bed = coverage_bed

            vardict_vcf = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, tumor_pair.name + ".vardict.germline.vt.vcf.gz")

            input_vcf = None
            normal = None
            tumor = None
            if os.path.isfile(vardict_vcf):
                input_vcf = vardict_vcf
                normal = tumor_pair.normal.name
                tumor = tumor_pair.tumor.name

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
                    name="cnvkit_batch." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
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
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            reference=pool_ref_cnn,
                            ref_cnn=ref_cnn
                        ),
                        cnvkit.segment(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            tumor_cns
                        )
                    ],
                    name="cnvkit_batch.correction." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
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
                            tumor_cns,
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns")
                        ),
                        pipe_jobs(
                            [
                                cnvkit.export(
                                    os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                                    None,
                                    sample_id=tumor_pair.tumor.name
                                ),
                                htslib.bgzip_tabix(
                                    None,
                                    vcf_gz
                                )
                            ]
                        )
                    ],
                    name="cnvkit_batch.call." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
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
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(metrics, tumor_pair.name + ".metrics.tsv")
                        ),
                        cnvkit.scatter(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".scatter.pdf"),
                            input_vcf,
                            normal,
                            tumor
                        ),
                        cnvkit.diagram(
                            os.path.join(cnvkit_dir, tumor_pair.name + ".cnr"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".call.cns"),
                            os.path.join(cnvkit_dir, tumor_pair.name + ".diagram.pdf")
                        )
                    ],
                    name="cnvkit_batch.metrics." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)]
                )
            )
        return jobs


    def gridss_paired_somatic(self):
        """
        """
        jobs = []
        for tumor_pair in self.tumor_pairs.values():
            if tumor_pair.multiple_normal == 1:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name, tumor_pair.name)
            else:
                normal_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.normal.name)
            tumor_alignment_directory = os.path.join(self.output_dirs['alignment_directory'], tumor_pair.tumor.name)

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            #svprep_directory = os.path.join(pair_directory, "gridss", "sv_prep")

            [input_normal] = self.select_input_files(
                [
                   # [os.path.join(svprep_directory, tumor_pair.normal.name + ".sv_prep.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.recal.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.dup.bam")],
                    [os.path.join(normal_alignment_directory, tumor_pair.normal.name + ".sorted.bam")]
                ]
            )
            [input_tumor] = self.select_input_files(
                [
                   # [os.path.join(svprep_directory, tumor_pair.tumor.name + ".sv_prep.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.recal.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.dup.bam")],
                    [os.path.join(tumor_alignment_directory, tumor_pair.tumor.name + ".sorted.bam")]
                ]
            )

            pair_directory = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            gridss_directory = os.path.join(pair_directory, "gridss")
            normal_output_prefix = os.path.join(gridss_directory, tumor_pair.normal.name)
            tumor_output_prefix = os.path.join(gridss_directory, tumor_pair.tumor.name)
            gridss_vcf_output = tumor_output_prefix + ".gridss.vcf.gz"
            gridds_bam_output = tumor_output_prefix + ".assembly.bam"

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(
                            gridss_directory,
                            remove=True
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
                    name="gridss_paired_somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                    input_dependency=[input_normal, input_tumor]
                )
            )

            gripss_vcf_output = tumor_output_prefix + ".gripss.somatic.vcf.gz"
            gripss_filter_vcf_output = tumor_output_prefix + ".gripss.filtered.somatic.vcf.gz"
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
                    name="gripss_filter.somatic." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                )
            )

            gripss_vcf_output = normal_output_prefix + ".gripss.germline.vcf.gz"
            gripss_filter_vcf_output = normal_output_prefix + ".gripss.filtered.germline.vcf.gz"
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
                    name="gripss_filter.germline." + tumor_pair.name,
                    samples=[tumor_pair.normal, tumor_pair.tumor],
                    readsets=[*list(tumor_pair.normal.readsets), *list(tumor_pair.tumor.readsets)],
                 )
            )
        return jobs

    def purple_sv(self):
        """
        """
        jobs = self.purple(sv=True)
        return jobs

    def linx_annotations_somatic(self):
        """
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            purple_dir = os.path.join(pair_dir, "purple")
            purple_vcf = os.path.join(purple_dir, tumor_pair.tumor.name + ".purple.sv.vcf.gz")
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
                    name="linx_annotations_somatic." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                )
            )
        return jobs

    def linx_annotations_germline(self):
        """
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            purple_dir = os.path.join(pair_dir, "purple")
            purple_vcf = os.path.join(purple_dir, tumor_pair.tumor.name + ".purple.sv.vcf.gz")
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
                    name="linx_annotations_germline." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets)
                   )
            )
        return jobs

    def linx_plot(self):
        """
        Generate Linx Plot of the tumor pair analysis.
        """
        jobs = []

        for tumor_pair in self.tumor_pairs.values():
            pair_dir = os.path.join(self.output_dirs['sv_variants_directory'], tumor_pair.name)
            linx_output_dir = os.path.join(pair_dir, "linx")
            linx_plot_dir = os.path.join(linx_output_dir, "plot")
            linx_zip = os.path.join(linx_output_dir, tumor_pair.name + ".linx_plot.zip")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(linx_output_dir),
                        linx.plot(
                            tumor_pair.tumor.name,
                            linx_output_dir,
                            ini_section="linx_plot"
                        ),
                        bash.touch(os.path.join(linx_output_dir, "linx_plot." + tumor_pair.name + ".Done")),
                        bash.zip(
                            linx_plot_dir,
                            linx_zip,
                            recursive=True
                        )
                    ],
                    name="linx_plot." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets),
                    output_dependency=[
                        os.path.join(linx_output_dir, "linx_plot." + tumor_pair.name + ".Done"),
                        linx_zip
                        ]
                )
            )
        return jobs
    
    def report_djerba(self):
        """
        Produce Djerba report.
        """
        jobs = []
        ensemble_directory = os.path.join(
            self.output_dirs['paired_variants_directory'],
            "ensemble"
            )
        assembly = config.param('report_pcgr', 'assembly')
        
        for tumor_pair in self.tumor_pairs.values():
            djerba_dir = os.path.join(self.output_dirs['report'][tumor_pair.name], "djerba")
            purple_dir = os.path.join(self.output_dirs['paired_variants_directory'], tumor_pair.name, "purple") # has to be a zipped directory, create zip file as part of job
            purple_zip = os.path.join(djerba_dir, tumor_pair.tumor.name + ".purple.zip")
            
            cpsr_directory = os.path.join(ensemble_directory, tumor_pair.name, "cpsr")
            input_cpsr = os.path.join(cpsr_directory, tumor_pair.name + ".cpsr." + assembly + ".json.gz")
            input_vcf = os.path.join(ensemble_directory, tumor_pair.name, tumor_pair.name + ".ensemble.somatic.vt.annot.2caller.flt.vcf.gz")
            pcgr_directory = os.path.join(djerba_dir, "pcgr")
            input_maf = os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".maf")
            clean_maf =  os.path.join(pcgr_directory, tumor_pair.name + ".pcgr_acmg." + assembly + ".clean.maf") # MAF from pcgr version 1.4.1 required, remove any empty t_depth lines, needs to be gzipped
            
            provenance_decoy = os.path.join(djerba_dir, "provenance_subset.tsv.gz")
            config_file = os.path.join(djerba_dir, tumor_pair.name + ".djerba.ini")
            djerba_script = os.path.join(djerba_dir, "djerba_report." + tumor_pair.name + ".sh")

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(djerba_dir),
                        bash.mkdir(pcgr_directory),
                        pcgr.report(
                            input_vcf,
                            input_cpsr,
                            pcgr_directory,
                            tumor_pair.name,
                            ini_section='report_djerba'
                            ),# add pcgr job to create MAF in correct format (1.4.1), remove chrM, gzip.
                        djerba.clean_maf(
                            input_maf,
                            clean_maf
                            ),
                        bash.zip(
                            purple_dir,
                            purple_zip,
                            recursive=True
                            ),
                        bash.touch(provenance_decoy),
                        djerba.make_config(
                            config_file,
                            tumor_pair.name,
                            tumor_pair.tumor.name,
                            tumor_pair.normal.name,
                            clean_maf + ".gz",
                            purple_zip
                            ),
                        # djerba report requires internet connection. Script is produced but must be executed locally.
                        djerba.make_script(
                            config_file,
                            djerba_dir,
                            djerba_script
                            )
                    ],
                    name="report_djerba." + tumor_pair.name,
                    samples=[tumor_pair.tumor],
                    readsets=list(tumor_pair.tumor.readsets),
                    input_dependency=[input_vcf, os.path.join(purple_dir, tumor_pair.tumor.name + ".purple.purity.tsv")],
                    output_dependency=[config_file, djerba_script]
                    )
                )


        return jobs

    @property
    def steps(self):
        return [
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba,
                self.sambamba_sort,
                self.sambamba_merge_sam_files, #5
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.manta_sv_calls, #10
                self.rawmpileup_panel,
                self.paired_varscan2_panel,
                self.merge_varscan2_panel,
                self.preprocess_vcf_panel,
                self.snp_effect_panel, #15
                self.gemini_annotations_panel,
                self.conpair_concordance_contamination,
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc, #20
                self.sequenza,
                self.run_pair_multiqc,
                self.sym_link_report,
                self.sym_link_fastq_pair,
                self.sym_link_panel #25
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba,
                self.sambamba_sort,
                self.sambamba_merge_sam_files, #5
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.conpair_concordance_contamination, #10
                self.metrics_dna_picard_metrics,
                self.metrics_dna_sample_qualimap,
                self.metrics_dna_fastqc,
                self.sequenza,
                self.manta_sv_calls, #15
                self.strelka2_paired_somatic,
                self.strelka2_paired_germline,
                self.strelka2_paired_germline_snpeff,
                self.purple,
                self.rawmpileup, #20
                self.paired_varscan2,
                self.merge_varscan2,
                self.paired_mutect2,
                self.merge_mutect2,
                self.vardict_paired, #25
                self.merge_filter_paired_vardict,
                self.ensemble_somatic,
                self.gatk_variant_annotator_somatic,
                self.merge_gatk_variant_annotator_somatic,
                self.ensemble_germline_loh, #30
                self.gatk_variant_annotator_germline,
                self.merge_gatk_variant_annotator_germline,
                self.cnvkit_batch,
                self.filter_ensemble_germline,
                self.filter_ensemble_somatic, #35
                self.report_cpsr,
                self.report_pcgr,
                self.run_pair_multiqc,
                self.sym_link_fastq_pair,
                self.sym_link_final_bam, #40
                self.sym_link_report,
                self.sym_link_ensemble,
                self.report_djerba
            ],
            [
                self.picard_sam_to_fastq,
                self.skewer_trimming,
                self.bwa_mem_sambamba,
                self.sambamba_sort,
                self.sambamba_merge_sam_files, #5
                self.gatk_indel_realigner,
                self.sambamba_merge_realigned,
                self.sambamba_mark_duplicates,
                self.recalibration,
                self.manta_sv_calls, #10
                self.strelka2_paired_somatic,
                #self.sv_prep,
                self.gridss_paired_somatic,
                self.purple_sv,
                self.linx_annotations_somatic,
                self.linx_annotations_germline, #15
                self.linx_plot
            ]
        ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        TumorPair(protocol=['fastpass', 'ensemble', 'sv'])
