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
import math
import os
import re
import sys
import subprocess

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs

from ...bfx import bash_cmd as bash
from ...bfx import gq_seq_utils
from ...bfx import rmarkdown
from ...bfx import differential_expression
from ...bfx import tools

from .. import rnaseq

log = logging.getLogger(__name__)


class RnaSeqLight(rnaseq.RnaSeqRaw):
    __doc__ = rnaseq.RnaSeqRaw.__doc__

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        # Add pipeline specific arguments
        super().__init__(*args, **kwargs)

    def kallisto(self):
        """
            Run Kallisto on fastq files for a fast esimate of abundance.
        """
        transcriptome_file = global_conf.get('kallisto', 'transcriptome_idx', param_type="filepath")
        tx2genes_file = global_conf.get('kallisto', 'transcript2genes', param_type="filepath")
        bootstraps = global_conf.get('kallisto', 'bootstraps')
        other_param = global_conf.get('kallisto', 'other_options', required=False)

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

            # PAIRED
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam),
                                                  re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

                job_name = "kallisto." + readset.name
                output_dir = os.path.join(self.output_dir, "kallisto", readset.sample.name)
                parameters = "--bootstrap-samples=" + str(bootstraps)
                parameters = other_param + parameters if other_param else parameters
                job = tools.rnaseqLight_kallisto(fastq1, fastq2, transcriptome_file, tx2genes_file, output_dir,
                                                 parameters, job_name)
                job.samples = [readset.sample]
                jobs.append(job)

            # SINGLE
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)

                job_name = "kallisto." + readset.name
                output_dir = os.path.join(self.output_dir, "kallisto", readset.sample.name)
                fragment_length = global_conf.get('kallisto', 'fragment_length')
                fragment_length_sd = global_conf.get('kallisto', 'fragment_length_sd')
                # warn user to update parameters in ini file?
                # print("Please make sure to update fragment_length and fragment_length_sd in the ini file!")
                parameters = " --single -l " + fragment_length + " -s " + fragment_length_sd
                parameters = other_param + parameters if other_param else parameters
                job = tools.rnaseqLight_kallisto(fastq1, "", transcriptome_file, tx2genes_file, output_dir, parameters,
                                                 job_name)
                job.samples = [readset.sample]
                jobs.append(job)
            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

        return jobs

    def kallisto_count_matrix(self):
        """
        Use the output from Kallisto to create a transcript count matrix.
        """

        jobs = []
        kallisto_directory = "kallisto"
        all_readset_directory = "All_readsets"
        output_dir = os.path.join(self.output_dir, kallisto_directory, all_readset_directory)

        # per trancripts
        input_abundance_files_transcripts = [
            os.path.join(self.output_dir, kallisto_directory, readset.sample.name, "abundance_transcripts.tsv") for
            readset in self.readsets]
        job_name_transcripts = "kallisto_count_matrix.transcripts"
        data_type_transcripts = "transcripts"
        job = tools.r_create_kallisto_count_matrix(input_abundance_files_transcripts, output_dir, data_type_transcripts,
                                                   job_name_transcripts)
        job.samples = self.samples
        jobs.append(job)

        # per genes
        input_abundance_files_genes = [
            os.path.join(self.output_dir, kallisto_directory, readset.sample.name, "abundance_genes.tsv") for readset in
            self.readsets]
        job_name_genes = "kallisto_count_matrix.genes"
        data_type_genes = "genes"
        job = tools.r_create_kallisto_count_matrix(input_abundance_files_genes, output_dir, data_type_genes,
                                                   job_name_genes)
        job.samples = self.samples
        jobs.append(job)

        report_dir = os.path.join(self.output_dir, "report")
        # copy tx2genes file
        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(report_dir),
                    Job(
                        [os.path.join(self.output_dir, "kallisto", "All_readsets", "all_readsets.abundance_genes.csv"),
                         os.path.join(self.output_dir, "kallisto", "All_readsets",
                                      "all_readsets.abundance_transcripts.csv")],
                        [],
                        command="""\
cp \\
  {tx2genes_file} \\
  {report_dir}""".format(
                            tx2genes_file=global_conf.get('kallisto', 'transcript2genes', param_type="filepath"),
                            report_dir=report_dir
                        )
                    )
                ],
                name="report.copy_tx2genes_file",
                samples=self.samples
            )
        )

        # Create kallisto report
        readset_merge_trim_stats = os.path.join("metrics", "trimReadsetTable.tsv")  # set in merge trimmomatic stats

        jobs.append(
            rmarkdown.render(
                job_input=[os.path.join(self.output_dir, "kallisto"
                                        , "All_readsets", "all_readsets.abundance_genes.csv")
                    , os.path.join(self.output_dir, "kallisto", "All_readsets",
                                   "all_readsets.abundance_transcripts.csv")
                    , readset_merge_trim_stats],
                job_name="report.kallisto_count_matrix",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqLight.kallisto.Rmd"),
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='report_dir="report";'
            )
        )

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq_light(self):
        """
        Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight
        """

        jobs = []
        abundance_file = os.path.join(self.output_dir, "kallisto/All_readsets", "all_readsets.abundance_genes.csv")
        # gqSeqUtils function call
        jobs.append(concat_jobs([
            Job(command="mkdir -p exploratory"),
            gq_seq_utils.exploratory_analysis_rnaseq_light(
                abundance_file,
                global_conf.get('gq_seq_utils_exploratory_analysis_rnaseq_light', 'genes', param_type='filepath'),
                "exploratory"
            )
        ], name="gq_seq_utils_exploratory_analysis_rnaseq_light", samples=self.samples))

        jobs.append(
            rmarkdown.render(
                job_input=os.path.join("exploratory", "index.tsv"),
                job_name="report.gq_seq_utils_exploratory_analysis_rnaseq",
                input_rmarkdown_file=os.path.join(self.report_template_dir,
                                                  "RnaSeqLight.gq_seq_utils_exploratory_analysis_rnaseq_light.Rmd"),
                samples=self.samples,
                render_output_dir='report',
                module_section='report',
                prerun_r='report_dir="report";'
            )
        )

        return jobs

    def sleuth_differential_expression(self):
        """
        Performs differential gene expression analysis using [Sleuth](http://pachterlab.github.io/sleuth/).
        Analysis are performed both at a transcript and gene level, using two different tests: LRT and WT.
        """

        # If --design <design file> option is missing, self.contrasts call will raise an Exception

        if self.contrasts:
            design_file = os.path.relpath(self.design_file.name, self.output_dir)
        output_directory = "sleuth"
        count_matrix = os.path.join(self.output_dir, "kallisto", "All_readsets", "all_readsets.abundance_genes.csv")
        tx2gene = global_conf.get('sleuth_differential_expression', 'tx2gene')

        sleuth_job = differential_expression.sleuth(design_file, count_matrix, tx2gene, output_directory)
        sleuth_job.output_files = [os.path.join(output_directory, contrast.name, "results.wt.gene.csv") for contrast in
                                   self.contrasts]
        sleuth_job.samples = self.samples

        return [concat_jobs([
            Job(command="mkdir -p " + output_directory),
            sleuth_job], name="sleuth_differential_expression")]

    @property
    def step_list(self):
        return self.protocols()[self._protocol]

    def protocols(self):
        return {'default': [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.kallisto,
            self.kallisto_count_matrix,
            self.gq_seq_utils_exploratory_analysis_rnaseq_light,
            self.sleuth_differential_expression
        ]}

