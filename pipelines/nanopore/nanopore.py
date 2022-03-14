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
import os
import sys
import logging
import collections
import re

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
import utils.utils
from core.config import global_config_parser, SanitycheckError, _raise
from core.job import Job, concat_jobs, pipe_jobs
from bfx.readset import parse_nanopore_readset_file
from pipelines import common
from bfx import minimap2
from bfx import sambamba
from bfx import svim
from bfx import pycoqc
from bfx import tools
from bfx import gatk4
from bfx import bash_cmd as bash

log = logging.getLogger(__name__)

class Nanopore(common.Nanopore):
    """
    Nanopore Pipeline
    ==============

    The Nanopore is used to analyse long reads produced by the Oxford Nanopore Technologies (ONT) sequencers.
    Currently, the pipeline uses minimap2 to align reads to the reference genome. Additionally, it produces
    a QC report that includes an interactive dashboard with data from the basecalling summary file as well
    as the alignment. A step aligning random reads to the NCBI nt database and reporting the species of the
    highest hits is also done as QC.

    Once the QC and alignments have been produced, Picard is used to merge readsets coming from the same
    sample. Finally, SVIM is used to detect Structural Variants (SV) including deletions, insertions and
    translocations. For a full summary of the types of SVs detected, please consult the following [site](
    https://github.com/eldariont/svim#background-on-structural-variants-and-long-reads).

    The SV calls produced by SVIM are saved as VCFs for each sample, which can then be used in downstream
    analyses. No filtering is performed on the SV calls.

    This pipeline currently does not perform base calling and requires both FASTQ and a sequencing_summary
    file produced by a ONT supported basecaller (we recommend Guppy). Additionally, the testing and
    development of the pipeline were focused on genomics applications, and functionality has not been tested
    for transcriptomics or epigenomics datasets.

    For more information on using ONT data for structural variant detection, as well as an alternative
    approach, please consult [this GitHub repository](https://github.com/nanoporetech/pipeline-structural-variation).

    For information on the structure and contents of the Nanopore readset file, please consult [here](
    https://bitbucket.org/mugqic/genpipes/src/master/#markdown-header-nanopore).
    """

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        super(Nanopore, self).__init__(*args, **kwargs)

    def guppy(self):
        """
        Use the Guppy basecaller to perform basecalling on all raw fast5 files.
        Uses the 'flip-flop' basecalling model by default.
        IN DEVELOPMENT
        """
        jobs = []

        return jobs

    def blastqc(self):
        """
        Uses BLAST to perform a basic QC test by aligning 1000bp of randomly selected
        reads to the NCBI nt database in order to detect potential contamination.
        """
        jobs = []

        for readset in self.readsets:

            blast_directory = os.path.join("blastQC", readset.name)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            job = tools.sh_blastQC_ONT(blast_directory, reads_fastq_dir, readset.name)
            job.samples = [readset.sample]
            jobs.append(job)

        return jobs

    def minimap2_align(self):
        """
        Uses minimap2 to align the Fastq reads that passed the minimum QC threshold to
        the provided reference genome. By default, it aligns to GRCh38.
        """
        jobs = []

        for readset in self.readsets:

            alignment_directory = os.path.join("alignment", readset.sample.name, readset.name)
            out_bam = os.path.join(alignment_directory, readset.name + ".sorted.bam")
            out_bai = re.sub("\.bam$", ".bam.bai", out_bam)

            if readset.fastq_files:
                reads_fastq_dir = readset.fastq_files
            else:
                _raise(SanitycheckError("Error: FASTQ file not available for readset \"" + readset.name + "\"!"))

            read_group = "'@RG" + \
                         "\\tID:" + readset.name + \
                         "\\tSM:" + readset.sample.name + \
                         "\\tLB:" + (readset.library if readset.library else readset.sample.name) + \
                         ("\\tPU:run" + readset.run if readset.run else "") + \
                         "\\tPL:Nanopore" + \
                         "'"
            job = concat_jobs(
                [
                    pipe_jobs(
                        [
                            bash.mkdir(os.path.dirname(out_bam)),
                            minimap2.minimap2_ont(
                                reads_fastq_dir,
                                read_group,
                                ini_section= "minimap2_align"
                            ),
                            sambamba.view(
                                "/dev/stdin",
                                None,
                                options="-S -f bam"
                            ),
                            sambamba.sort(
                                "/dev/stdin",
                                out_bam,
                                tmp_dir=global_config_parser.param('minimap2_align', 'tmp_dir', required=True),
                            )
                        ]
                    ),
                    sambamba.index(
                        out_bam,
                        out_bai,
                    )
                ],
                name="minimap2_align." + readset.name,
                samples=[readset.sample]
            )
            jobs.append(job)

        return jobs

    def pycoqc(self):
        """
        Use pycoQC to produce an interactive quality report based on the summary file and
        alignment outputs.
        """
        jobs = []

        for readset in self.readsets:

            pycoqc_directory = os.path.join("pycoQC", readset.name)

            if readset.summary_file:
                in_summary = readset.summary_file
            else:
                _raise(SanitycheckError("Error: summary file not available for readset \"" + readset.name + "\"!"))

            align_directory = os.path.join("alignment", readset.sample.name, readset.name)
            in_bam = os.path.join(align_directory, readset.name + ".sorted.bam")

            jobs.append(
                concat_jobs([
                    bash.mkdir(pycoqc_directory),
                    pycoqc.pycoqc(
                        readset_name=readset.name,
                        input_summary=in_summary,
                        output_directory=pycoqc_directory,
                        input_barcode=None,
                        input_bam=in_bam
                        )
                ],
                    name="pycoqc." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs

    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample.
        Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:
        Aligned and sorted BAM output files from previous minimap2_align step
        """
        jobs = []

        for sample in self.samples:

            alignment_directory = os.path.join("alignment", sample.name)

            # Find input readset BAMs first from previous minimap2_align job,
            readset_bams = self.select_input_files([
                [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets]
            ])

            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            mkdir_job = bash.mkdir(os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:

                readset_bam = readset_bams[0]

                readset_index = re.sub("\.bam$", ".bam.bai", readset_bam)
                sample_index = re.sub("\.bam$", ".bam.bai", sample_bam)

                job = concat_jobs(
                    [
                        mkdir_job,
                        bash.ln(
                            readset_bam,
                            sample_bam
                        ),
                        bash.ln(
                            readset_index,
                            sample_index
                        )
                    ],
                    name="symlink_readset_sample_bam." + sample.name,
                    samples=[sample],
                )
                job.samples = [sample]

            elif len(sample.readsets) > 1:

                job = concat_jobs(
                    [
                        mkdir_job,
                        gatk4.merge_sam_files(readset_bams, sample_bam)
                    ],
                    samples=[sample],
                    name="picard_merge_sam_files." + sample.name
                )

            jobs.append(job)

        return jobs

    def svim(self):
        """
        Use SVIM to perform SV calling on each sample.
        """
        jobs = []

        for sample in self.samples:

            align_directory = os.path.join("alignment", sample.name)
            in_bam = os.path.join(align_directory, sample.name + ".sorted.bam")

            svim_directory = os.path.join("svim", sample.name)

            job = svim.svim_ont(in_bam, svim_directory)
            job.name = "svim." + sample.name
            job.samples = [sample]
            jobs.append(job)

        return jobs

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return { 'default': [
            self.blastqc,
            self.minimap2_align,
            self.pycoqc,
            self.picard_merge_sam_files,
            self.svim
        ]}


def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Check if Genpipes must be ran inside a container
    utils.container_wrapper_argparse(__file__, argv)
    # Build help
    epilog = Nanopore.process_help(argv)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        conflict_handler='resolve', epilog=epilog)

    # populate the parser
    parser = Nanopore.argparser(parser)

    parsed_args = parser.parse_args(argv)

    sanity_check = parsed_args.sanity_check
    loglevel = parsed_args.log
    utils.set_logger(loglevel, sanity_check=sanity_check)

    # Pipeline config
    config_files = parsed_args.config

    # Common Pipeline options
    genpipes_file = parsed_args.genpipes_file
    container = parsed_args.container
    clean = parsed_args.clean
    report = parsed_args.report
    no_json = parsed_args.no_json
    force = parsed_args.force
    job_scheduler = parsed_args.job_scheduler
    output_dir = parsed_args.output_dir
    steps = parsed_args.steps
    readset_file = parsed_args.readsets_file
    design_file = parsed_args.design_file

    pipeline = Nanopore(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file,
                         clean=clean, report=report, force=force, job_scheduler=job_scheduler, output_dir=output_dir,
                         design_file=design_file, no_json=no_json, container=container)

    pipeline.submit_jobs()

if __name__ == '__main__':
    main()
