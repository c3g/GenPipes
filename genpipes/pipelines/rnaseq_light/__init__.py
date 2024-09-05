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

# MUGQIC Modules
from ...core.config import global_conf, _raise, SanitycheckError
from ...core.job import Job, concat_jobs

from ...bfx import (
    bash_cmd as bash,
    differential_expression,
    gq_seq_utils,
    job2json_project_tracking,
    kallisto,
    multiqc,
    tools
    )

from .. import rnaseq

log = logging.getLogger(__name__)

class RnaSeqLight(rnaseq.RnaSeqRaw):
    """
RNA-Seq Light Pipeline
================

The RNA-Seq Light Pipeline is a lightweight analysis of gene expression in RNA sequencing data. 
The pipeline is based on [Kallisto](https://pachterlab.github.io/kallisto/about.html) and differential expression analysis is performed by [Sleuth](http://pachterlab.github.io/sleuth/).
It is especially useful for quick Quality Control (QC) in gene sequencing studies.
    """

    def __init__(self, *args, protocol=None, **kwargs):
        if protocol is None:
            self._protocol = 'default'
        # Add pipeline specific arguments
        super().__init__(*args, **kwargs)

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
        Run [Kallisto](https://pachterlab.github.io/kallisto/about.html) on fastq files for a fast esimate of abundance.
        """

        transcriptome_file = global_conf.global_get('kallisto', 'transcriptome_idx', param_type="filepath")
        tx2genes_file = global_conf.global_get('kallisto', 'transcript2genes', param_type="filepath")
        bootstraps = global_conf.global_get('kallisto', 'bootstraps')
        other_param = global_conf.global_get('kallisto', 'other_options', required=False)

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
                        candidate_input_files.append([re.sub(r"\.bam$", ".pair1.fastq.gz", readset.bam), re.sub(r"\.bam$", ".pair2.fastq.gz", readset.bam)])
                    [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                    input_fastqs.extend([fastq1, fastq2])

                    parameters = "--bootstrap-samples=" + str(bootstraps)

                #SINGLE
                elif readset.run_type == "SINGLE_END":
                    candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                    if readset.fastq1:
                        candidate_input_files.append([readset.fastq1])
                    if readset.bam:
                        candidate_input_files.append([re.sub(r"\.bam$", ".single.fastq.gz", readset.bam)])
                    [fastq1] = self.select_input_files(candidate_input_files)
                    input_fastqs.append(fastq1)

                    fragment_length = global_conf.global_get('kallisto', 'fragment_length', required=True)
                    fragment_length_sd = global_conf.global_get('kallisto', 'fragment_length_sd', required=True)
                    parameters = "--single -l "+ fragment_length +" -s " + fragment_length_sd + " --bootstrap-samples=" + str(bootstraps)

                else:
                    _raise(SanitycheckError(f"""Error: run type "{readset.run_type}" is invalid for readset "{readset.name}" (should be PAIRED_END or SINGLE_END)!"""))

            output_dir = os.path.join(self.output_dirs["kallisto_directory"], sample.name)
            link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
            parameters = " ".join([other_param, parameters]) if other_param else parameters
            job_name = f"kallisto.{sample.name}"
            job_project_tracking_metrics = []
            if self.project_tracking_json:
                job_project_tracking_metrics = concat_jobs(
                    [
                    kallisto.parse_mean_insert_size_metrics_pt(os.path.join(output_dir, "abundance.h5")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_dir, "abundance.h5"),
                        pipeline=self,
                        samples=sample.name,
                        readsets=",".join([readset.name for readset in sample.readsets]),
                        job_name=job_name,
                        metrics="mean_insert_size=$mean_insert_size"
                        ),
                    kallisto.parse_median_insert_size_metrics_pt(os.path.join(output_dir, "abundance.h5")),
                    job2json_project_tracking.run(
                        input_file=os.path.join(output_dir, "abundance.h5"),
                        pipeline=self,
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
                        bash.mkdir(link_directory),
                        kallisto.quant(
                            input_fastqs,
                            output_dir,
                            transcriptome_file,
                            parameters
                        ),
                        bash.ln(
                            os.path.relpath(os.path.join(output_dir, "kallisto_quant.log"), link_directory),
                            os.path.join(link_directory, sample.name + ".kallisto_quant.log"),
                            input=os.path.join(output_dir, "kallisto_quant.log")
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
            self.multiqc_inputs.append(os.path.join(link_directory, sample.name + ".kallisto_quant.log"))

        return jobs

    def kallisto_count_matrix(self):
        """
        Use the output from Kallisto to create a transcript count matrix.
        Create a summary table to be included in the multiqc report.
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
                        global_conf.global_get('kallisto', 'transcript2genes', param_type="filepath"),
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
        link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")

        inputs = [
            os.path.join(self.output_dirs["metrics_directory"], "trimReadsetTable.tsv"),
            os.path.join(output_dir, "all_readsets.abundance_genes.csv"),
            os.path.join(output_dir, "all_readsets.abundance_transcripts.csv")
            ]
        kallisto_report_file = os.path.join(output_dir, "all_readsets.kallisto_report.tsv")
        kallisto_multiqc_file = os.path.join(link_directory, "all_readsets.kallisto_mqc.txt")
        
        kallisto_report_job = tools.r_create_kallisto_report(
            report_dir,
            inputs,
            kallisto_report_file
        )

        kallisto_multiqc_format_job = Job(
            [kallisto_report_file],
            [kallisto_multiqc_file],
            command="""\
echo -e "# plot_type: 'table'
# section_name: 'Kallisto' 
# headers:
#   RawReads:
#       title: 'Raw Reads'
#       description: 'total number of reads obtained from the sequencer'
#       format: '{{:,.0f}}'
#       placement: 930
#   SurvivingReads:
#       title: 'Surviving Reads'
#       description: 'number of remaining reads after the trimming step'
#       format: '{{:,.0f}}'
#       placement: 940
#   PercSurviving:
#       title: '% Surviving'
#       description: 'Surviving Reads / Raw Reads * 100'
#       placement: 950
#   Transcriptome:
#       title: 'Transcriptome targets'
#       description: 'number of transcript targets'
#       format: '{{:,.0f}}'
#       placement: 960
#   Transcripts:
#       title: 'Transcripts'
#       description: 'number of transcripts with at least 5 reads'
#       format: '{{:,.0f}}'
#       placement: 970
#   TranscriptsReads:
#       title: 'Transcripts reads'
#       description: 'total number of reads covering the transcripts'
#       format: '{{:,.0f}}'
#       placement: 980
#   PercTranscriptsReads:
#       title: '% Transcripts Reads'
#       description: 'Transcripts Reads # / Surviving reads * 100'
#       placement: 990
#   Genes:
#       title: 'Genes'
#       description: 'number of Genes with at least 5 reads'
#       format: '{{:,.0f}}'
#       placement: 1000
#   GenesReads:
#       title: 'Genes Reads'
#       description: 'total number of reads covering the genes'
#       format: '{{:,.0f}}'
#       placement: 1010
#   PercGenesReads:
#       title: '% Genes Reads'
#       description: 'Genes Reads # / Surviving reads * 100'
#       placement: 1020" > {kallisto_multiqc_file}

cat {kallisto_report_file} >> {kallisto_multiqc_file}""".format(
    kallisto_report_file=kallisto_report_file,
    kallisto_multiqc_file=kallisto_multiqc_file
    )
)

        jobs.append(
            concat_jobs(
                [
                    bash.mkdir(link_directory),
                    kallisto_report_job,
                    kallisto_multiqc_format_job
                ],
                name="report.kallisto_count_matrix"
            )
        )
        self.multiqc_inputs.append(kallisto_multiqc_file)

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
                global_conf.global_get('gq_seq_utils_exploratory_analysis_rnaseq_light', 'genes', param_type='filepath'),
                self.output_dirs['exploratory_directory']
            )],
            name="gq_seq_utils_exploratory_analysis_rnaseq_light",
            samples=self.samples,
            readsets=self.readsets
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
        output_directory = self.output_dirs["sleuth_directory"]
        link_directory = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
        count_matrix = os.path.join(self.output_dirs["kallisto_directory"], "All_readsets", "all_readsets.abundance_genes.csv")
        tx2gene = global_conf.global_get('sleuth_differential_expression', 'tx2gene')

        sleuth_job = differential_expression.sleuth(design_file, count_matrix, tx2gene, output_directory)
        sleuth_job.output_files = [os.path.join(output_directory, contrast.name, "results.wt.gene.csv") for contrast in self.contrasts]
        sleuth_job.samples = self.samples
        sleuth_job.readsets = self.readsets

        multiqc_job = bash.mkdir(link_directory)
        for contrast in self.contrasts:
            heatmap = os.path.join(output_directory, contrast.name, "heatmap.topFCgenes.png")
            pca = os.path.join(output_directory, contrast.name, "pca_plot.trx.png")
            sleuth_job.output_files.extend([heatmap, pca])

            multiqc_job = concat_jobs(
                [
                    multiqc_job,
                    bash.ln(
                        os.path.relpath(heatmap, link_directory),
                        os.path.join(link_directory, f"Heatmap_{contrast.name}_mqc.png"),
                        input = heatmap
                    ),
                    bash.ln(
                        os.path.relpath(pca, link_directory),
                        os.path.join(link_directory, f"PCA_{contrast.name}_mqc.png"),
                        input = pca
                    )
                ]
            )
            self.multiqc_inputs.extend(
                [
                    os.path.join(link_directory, f"Heatmap_{contrast.name}_mqc.png"),
                    os.path.join(link_directory, f"PCA_{contrast.name}_mqc.png")
                ]
            )

        return [
            concat_jobs(
                [
                    bash.mkdir(output_directory),
                    sleuth_job,
                    multiqc_job
                ],
                name="sleuth_differential_expression"
            )
        ]
    
    def multiqc(self):
        """
        A quality control report for all samples is generated.
        For more detailed information about MultiQC visit: [MultiQC](http://multiqc.info/)
        """
        jobs = []
        
        input_links = os.path.join(self.output_dirs["metrics_directory"], "multiqc_inputs")
            
        output = os.path.join(self.output_dirs['report_directory'], "RnaSeqLight.multiqc")

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
        return {'default': [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.kallisto,
            self.kallisto_count_matrix,
            self.gq_seq_utils_exploratory_analysis_rnaseq_light,
            self.sleuth_differential_expression,
            self.multiqc
        ]}

def main(parsed_args):
    """
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

    pipeline = RnaSeqLight(config_files, genpipes_file=genpipes_file, steps=steps, readsets_file=readset_file, clean=clean, force=force, force_mem_per_cpu=force_mem_per_cpu, job_scheduler=job_scheduler, output_dir=output_dir, design_file=design_file, no_json=no_json, json_pt=json_pt, container=container)

    pipeline.submit_jobs()

