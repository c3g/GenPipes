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
import subprocess

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs
import utils.utils

from bfx import (
    bash_cmd as bash,
    differential_expression,
    gq_seq_utils,
    job2json_project_tracking,
    kallisto,
    rmarkdown,
    tools
    )

from pipelines import common
from pipelines.rnaseq import rnaseq

log = logging.getLogger(__name__)

class RnaSeqLight(rnaseq.RnaSeqRaw):
    def __init__(self,protocol=None):
        self._protocol=protocol
        super(RnaSeqLight, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
            'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
            'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
            'kallisto_directory': os.path.relpath(os.path.join(self.output_dir, 'kallisto'), self.output_dir),
            'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
            'exploratory_directory': os.path.relpath(os.path.join(self.output_dir, 'exploratory'), self.output_dir),
            'sleuth_directory': os.path.relpath(os.path.join(self.output_dir, 'sleuth'), self.output_dir),
            'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir)
        }
        return dirs

    def kallisto(self):
        """
        Run Kallisto on fastq files for a fast esimate of abundance.
        """

        transcriptome_file = config.param('kallisto', 'transcriptome_idx', param_type="filepath")
        tx2genes_file = config.param('kallisto', 'transcript2genes', param_type="filepath")
        bootstraps = config.param('kallisto', 'bootstraps')
        other_param = config.param('kallisto', 'other_options', required=False)

        jobs = []

        for sample in self.samples:
            readset_type = ""
            parameters = ""
            input_fastqs = []
            for readset in sample.readsets:
                trim_file_prefix = os.path.join(self.output_dirs["trim_directory"], sample.name, readset.name + ".trim.")
                if not readset_type:
                    readset_type = readset.run_type
                elif not readset_type == readset.run_type:
                    message = f"Sample {sample.name} has mixed single-end and paired-end readset libraries...\n"
                    message += f"Please use only single-end or only paired-end library for samples with multiple readsets."
                    _raise(SanitycheckError(message))

                #PAIRED
                if readset.run_type == "PAIRED_END":
                    candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                    if readset.fastq1 and readset.fastq2:
                        candidate_input_files.append([readset.fastq1, readset.fastq2])
                    if readset.bam:
                        candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                    [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                    input_fastqs.extend([fastq1, fastq2])

                    parameters = "--bootstrap-samples=" + str(bootstraps)

                #SINGLE
                elif readset.run_type == "SINGLE_END":
                    candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                    if readset.fastq1:
                        candidate_input_files.append([readset.fastq1])
                    if readset.bam:
                        candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                    [fastq1] = self.select_input_files(candidate_input_files)
                    input_fastqs.append(fastq1)

                    fragment_length = config.param('kallisto', 'fragment_length', required=True)
                    fragment_length_sd = config.param('kallisto', 'fragment_length_sd', required=True)
                    parameters = "--single -l "+ fragment_length +" -s " + fragment_length_sd + " --bootstrap-samples=" + str(bootstraps)

                else:
                    _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                    "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            output_dir = os.path.join(self.output_dirs["kallisto_directory"], sample.name)
            parameters = " ".join([other_param, parameters]) if other_param else parameters
            job_name = f"kallisto.{sample.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                    kallisto.parse_mean_insert_size_metrics_pt(os.path.join(output_dir, "abundance.h5")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_dir, "abundance.h5"),
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="mean_insert_size=$mean_insert_size"
                        ),
                    kallisto.parse_median_insert_size_metrics_pt(os.path.join(output_dir, "abundance.h5")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_dir, "abundance.h5"),
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="median_insert_size=$median_insert_size"
                        )
                    ])
            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(output_dir),
                        kallisto.quant(
                            input_fastqs,
                            output_dir,
                            transcriptome_file,
                            parameters
                        ),
                        bash.mv(
                            os.path.join(output_dir, "abundance.tsv"),
                            os.path.join(output_dir, "abundance_transcripts.tsv")
                        ),
                        tools.r_transcript_to_gene(
                            os.path.join(output_dir, "abundance_transcripts.tsv"),
                            os.path.join(output_dir, "abundance_genes.tsv"),
                            tx2genes_file
                        ),
                        job_project_tracking_metrics
                    ],
                    input_dependency=input_fastqs,
                    output_dependency=[
                        os.path.join(output_dir, "abundance_transcripts.tsv"),
                        os.path.join(output_dir, "abundance_genes.tsv"),
                        os.path.join(output_dir, "abundance.h5"),
                        os.path.join(output_dir, "kallisto_quant.log")
                    ],
                    name=job_name,
                    samples=[sample],
                    readsets=sample.readsets
                )
            )

        return jobs

    def kallisto_count_matrix(self):
        """
        Use the output from Kallisto to create a transcript count matrix.
        """

        jobs=[]
        output_dir = os.path.join(self.output_dirs["kallisto_directory"], "All_readsets")

        #per trancripts
        input_abundance_files_transcripts = [os.path.join(self.output_dirs["kallisto_directory"], sample.name, "abundance_transcripts.tsv") for sample in self.samples]
        job_name_transcripts="kallisto_count_matrix.transcripts"
        data_type_transcripts="transcripts"
        job = tools.r_create_kallisto_count_matrix(
            input_abundance_files_transcripts,
            output_dir,
            data_type_transcripts,
            job_name_transcripts
        )
        job.samples = self.samples
        job.readsets = self.readsets
        jobs.append(job)

        #per genes
        input_abundance_files_genes = [os.path.join(self.output_dirs["kallisto_directory"], sample.name, "abundance_genes.tsv") for sample in self.samples]
        job_name_genes="kallisto_count_matrix.genes"
        data_type_genes="genes"
        job = tools.r_create_kallisto_count_matrix(
            input_abundance_files_genes,
            output_dir,
            data_type_genes,
            job_name_genes
        )
        job.samples = self.samples
        job.readsets = self.readsets
        jobs.append(job)

        report_dir = self.output_dirs["report_directory"]
        # Copy tx2genes file
        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(report_dir),
                    bash.cp(
                        config.param('kallisto', 'transcript2genes', param_type="filepath"),
                        report_dir
                    )
                ],
                input_dependency=[
                    os.path.join(output_dir, "all_readsets.abundance_genes.csv"),
                    os.path.join(output_dir, "all_readsets.abundance_transcripts.csv")
                ],
                name="report.copy_tx2genes_file",
                samples=self.samples,
                readsets=self.readsets
            )
        )

        # Create kallisto report
        readset_merge_trim_stats = os.path.join(self.output_dirs["metrics_directory"], "trimReadsetTable.tsv") # set in merge trimmomatic stats
        jobs.append(
            rmarkdown.render(
                job_input=[
                    os.path.join(output_dir, "all_readsets.abundance_genes.csv"),
                    os.path.join(output_dir, "all_readsets.abundance_transcripts.csv"),
                    readset_merge_trim_stats
                ],
                job_name="report.kallisto_count_matrix",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqLight.kallisto.Rmd"),
                samples=self.samples,
                readsets=self.readsets,
                render_output_dir=self.output_dirs['report_directory'],
                module_section='report',
                prerun_r=f'report_dir="{self.output_dirs["report_directory"]}";'
            )
        )

        return jobs

    def gq_seq_utils_exploratory_analysis_rnaseq_light(self):
        """
        Exploratory analysis using the gqSeqUtils R package adapted for RnaSeqLight.
        """

        jobs = []
        abundance_file=os.path.join(self.output_dirs["kallisto_directory"], "All_readsets", "all_readsets.abundance_genes.csv")
        # gqSeqUtils function call
        jobs.append(concat_jobs([
            Job(command=f"mkdir -p {self.output_dirs['exploratory_directory']}"),
            gq_seq_utils.exploratory_analysis_rnaseq_light(
                abundance_file,
                config.param('gq_seq_utils_exploratory_analysis_rnaseq_light', 'genes', param_type='filepath'),
                self.output_dirs['exploratory_directory']
            )],
            name="gq_seq_utils_exploratory_analysis_rnaseq_light",
            samples=self.samples,
            readsets=self.readsets
            )
        )

        jobs.append(
            rmarkdown.render(
                job_input=os.path.join(self.output_dirs['exploratory_directory'], "index.tsv"),
                job_name="report.gq_seq_utils_exploratory_analysis_rnaseq",
                input_rmarkdown_file=os.path.join(self.report_template_dir, "RnaSeqLight.gq_seq_utils_exploratory_analysis_rnaseq_light.Rmd"),
                samples=self.samples,
                readsets=self.readsets,
                render_output_dir=self.output_dirs['report_directory'],
                module_section='report',
                prerun_r=f'report_dir="{self.output_dirs["report_directory"]}";'
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
            design_file = os.path.relpath(self.args.design.name, self.output_dir)
        output_directory = self.output_dirs["sleuth_directory"]
        count_matrix = os.path.join(self.output_dirs["kallisto_directory"], "All_readsets", "all_readsets.abundance_genes.csv")
        tx2gene = config.param('sleuth_differential_expression', 'tx2gene')

        sleuth_job = differential_expression.sleuth(design_file, count_matrix, tx2gene, output_directory)
        sleuth_job.output_files = [os.path.join(output_directory, contrast.name, "results.wt.gene.csv") for contrast in self.contrasts]
        sleuth_job.samples = self.samples
        sleuth_job.readsets = self.readsets

        return [
            concat_jobs(
                [
                    bash.mkdir(output_directory),
                    sleuth_job
                ],
                name="sleuth_differential_expression"
            )
        ]

############

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.kallisto,
            self.kallisto_count_matrix,
            self.gq_seq_utils_exploratory_analysis_rnaseq_light,
            self.sleuth_differential_expression
            ]

if __name__ == '__main__':
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        RnaSeqLight()
