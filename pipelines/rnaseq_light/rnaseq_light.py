#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import logging
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.design import *
from bfx.readset import *

# from bfx import bedtools
# from bfx import cufflinks
# from bfx import differential_expression
from bfx import gq_seq_utils
# from bfx import htseq
# from bfx import metrics
from bfx import picard
# from bfx import samtools
# from bfx import star
# from bfx import bvatools
from bfx import rmarkdown
from pipelines import common
import utils

from pipelines.rnaseq import rnaseq

log = logging.getLogger(__name__)

class RNAseqLight(rnaseq.RnaSeq):
	def __init__(self):
		super(RNAseqLight, self).__init__()

	def kallisto(self):
		"""
			Run Kallisto on fastq files for a fast esimate of abundance.
		"""
		transcriptome_file = config.param('kallisto', 'transcriptome_idx', type="filepath")
		gtf_file = config.param('DEFAULT', 'gtf_transcript_id', type="filepath")

		jobs = []
		for readset in self.readsets:
			trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")

			#PAIRED
			candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
			if readset.fastq1 and readset.fastq2:
				candidate_input_files.append([readset.fastq1, readset.fastq2])
			if readset.bam:
				candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
			[fastq1, fastq2] = self.select_input_files(candidate_input_files)

			job_name = "kallisto." + readset.name
			output_dir=self.output_dir+"/kallisto/" + readset.sample.name
			job = tools.rnaseqLight_kallisto(fastq1, fastq2, transcriptome_file, gtf_file, output_dir, job_name)
			jobs.append(job)

			#SINGLE
		return jobs

	def mergeKallistoCounts(self):

		kallisto_directory="kallisto"
		input_abundance_files = [os.path.join(kallisto_directory, readset.sample.name, "abundance_genes.tsv") for readset in self.readsets]

		output_dir=self.output_dir+"/kallisto/"
		job_name = "merge_kallisto"
		data_type="genes"


		job=[tools.r_merge_kallisto_counts(input_abundance_files, output_dir, data_type, job_name)]

		return job

	def gq_seq_utils_exploratory_analysis_rnaseq_light(self):
		"""
		Exploratory analysis using the gqSeqUtils R package adapted for RNAseqLight
		"""

		jobs = []

		# gqSeqUtils function call
		sample_fpkm_readcounts = [[
			sample.name,
			os.path.join("kallisto", sample.name, "abundance_transcripts.tsv"),
			os.path.join("kallisto", sample.name, "abundance_genes.tsv")
			# os.path.join("cufflinks", sample.name, "isoforms.fpkm_tracking"),
			# os.path.join("raw_counts", sample.name + ".readcounts.csv")
		] for sample in self.samples]
		jobs.append(concat_jobs([
			Job(command="mkdir -p exploratory"),
			gq_seq_utils.exploratory_analysis_rnaseq_light(
				os.path.join("kallisto", "all_samples.abundance_genes.csv"),
				config.param('gq_seq_utils_exploratory_analysis_rnaseq_light', 'genes', type='filepath'),
				"exploratory"
			)
		], name="gq_seq_utils_exploratory_analysis_rnaseq_light"))

		# Render Rmarkdown Report
		jobs.append(
			rmarkdown.render(
			 job_input            = os.path.join("exploratory", "index.tsv"),
			 job_name             = "gq_seq_utils_exploratory_analysis_rnaseq_report",
			 input_rmarkdown_file = os.path.join(self.report_template_dir, "RnaSeq.gq_seq_utils_exploratory_analysis_rnaseq_light.Rmd") ,
			 render_output_dir    = 'report',
			 module_section       = 'report', # TODO: this or exploratory?
			 prerun_r             = 'report_dir="report";' # TODO: really necessary or should be hard-coded in exploratory.Rmd?
			 )
		)

		# report_file = os.path.join("report", "RnaSeq.kallisto.md")
		# jobs.append(
		# 	Job(
		# 		[os.path.join("cufflinks", "AllSamples","merged.gtf")],
		# 		[report_file],
		# 		command="""\
		# 		mkdir -p report && \\
		# 		zip -r report/cuffAnalysis.zip cufflinks/ cuffdiff/ kallisto/ && \\
		# 		cp \\
		# 		  {report_template_dir}/{basename_report_file} \\
		# 		  {report_file}""".format(
		# 			report_template_dir=self.report_template_dir,
		# 			basename_report_file=os.path.basename(report_file),
		# 			report_file=report_file
		# 		),
		# 		report_files=[report_file],
		# 		name="kallisto_report")
		# )

		return jobs


############

	@property
	def steps(self):
		return [
			self.picard_sam_to_fastq,
			self.trimmomatic,
			self.merge_trimmomatic_stats,
			self.kallisto,
			#merge readsets to samples
			self.mergeKallistoCounts,
			self.gq_seq_utils_exploratory_analysis_rnaseq_light
			]

if __name__ == '__main__':
	RNAseqLight()

