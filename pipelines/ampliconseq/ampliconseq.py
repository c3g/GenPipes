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
from os.path import basename

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.sequence_dictionary import *

from pipelines import common
from bfx import flash
from bfx import trimmomatic
from bfx import qiime


log = logging.getLogger(__name__)

class AmpliconSeq(common.Illumina):
	"""
	Amplicon-Seq Pipeline
	================

	T
	"""

	def flash(self):
		"""
		Merge paired end reads using [FLASh](http://ccb.jhu.edu/software/FLASH/).
		"""
		jobs = []
	
		for readset in self.readsets:
			trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
			merge_directory = os.path.join("merge", readset.sample.name)
			merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")
			merge_file_prefix_log = os.path.join(merge_directory, readset.name + ".log")

		    # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
			if readset.run_type == "PAIRED_END":
				candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
				if readset.fastq1 and readset.fastq2:
					candidate_input_files.append([readset.fastq1, readset.fastq2])
				[fastq1, fastq2] = self.select_input_files(candidate_input_files)	
				
			else:
				raise Exception("Error: run type \"" + readset.run_type +
				"\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!")

			job = flash.flash(
				fastq1,
				fastq2,
				merge_directory,
				merge_file_prefix,
				readset.name,
				merge_file_prefix_log
			)	        

			jobs.append(concat_jobs([
                # Trimmomatic does not create output directory by default
				Job(command="mkdir -p " + merge_directory),
				job
			], name="flash." + readset.sample.name))
			
		return jobs
		
	def merge_flash_stats(self):
		"""
		The paired end merge statistics per readset are merged at this step.
		"""
		
		readset_merge_flash_stats = os.path.join("metrics", "mergeReadsetTable.tsv")
		job = concat_jobs([Job(command="mkdir -p metrics"), Job(command="echo 'Sample\tReadset\tTrim Paired Reads #\tMerged Paired Reads #\tMerged Paired Reads %' > " + readset_merge_flash_stats)])
		
		for readset in self.readsets:
			flash_log = os.path.join("merge", readset.sample.name, readset.name + ".log")
			
			job = concat_jobs([
				job,
				Job(command="""\
		printf '{}\t{}\t' >> {}""".format(readset.sample.name, readset.name, readset_merge_flash_stats)
				   )
				  ])
				  			
			# Retrieve merge statistics using re search in python.

			python_command = """python -c 'import re; import sys; log_file = open("{}","r"); merge_stat=[]; merge_stat.append([i.split()[3] for i in log_file if re.search("Total pairs",i)][0]); log_file.seek(0); merge_stat.append([i.split()[3] for i in log_file if re.search("Combined pairs",i)][0]); log_file.seek(0); merge_stat.append([i.split()[3] for i in log_file if re.search("Percent combined",i)][0][:-1]); log_file.close(); print "\t".join(merge_stat)'""".format(flash_log)
			
			job = concat_jobs([
				job,
				Job(
					[flash_log],
					[readset_merge_flash_stats],
					# Create readset merging stats TSV file with paired read count using python.
					command="""\
{python_command} \\
  >> {readset_merge_flash_stats}""".format(
						flash_log=flash_log,
						python_command=python_command,
						readset_merge_flash_stats=readset_merge_flash_stats
					)
				)
			])

		sample_merge_flash_stats = os.path.join("metrics", "mergeSampleTable.tsv")
		report_file = os.path.join("report", "Illumina.flash_stats.md")
		return [concat_jobs([
			job,
			Job(
				[readset_merge_flash_stats],
				[sample_merge_flash_stats],
				# Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
				command="""\
cut -f1,3- {readset_merge_flash_stats} | awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Trim Reads #", "Merged Reads #", "Merged %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_flash_stats}""".format(
					readset_merge_flash_stats=readset_merge_flash_stats,
					sample_merge_flash_stats=sample_merge_flash_stats
				)
			),
			Job(
				[sample_merge_flash_stats],
				[report_file],
				[['merge_flash_stats', 'module_pandoc']],
				command="""\
mkdir -p report && \\
cp {readset_merge_flash_stats} {sample_merge_flash_stats} report/ && \\
merge_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_flash_stats}` && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable min_overlap={min_overlap} \\
  --variable max_overlap={max_overlap} \\
  --variable read_type={read_type} \\
  --variable merge_readset_table="$merge_readset_table_md" \\
  --to markdown \\
  > {report_file}""".format(
					min_overlap=config.param('flash', 'min_overlap', type='int'),
					max_overlap=config.param('flash', 'max_overlap', type='int'),
					read_type="Paired",
					report_template_dir=self.report_template_dir,
					readset_merge_flash_stats=readset_merge_flash_stats,
					sample_merge_flash_stats=sample_merge_flash_stats,
					basename_report_file=os.path.basename(report_file),
					report_file=report_file
				),
				report_files=[report_file]
			)], name="merge_flash_stats")]
		
	
	def catenate(self):
	
		"""
		Catenate all the reads in one file for further analysis.
	
		This step takes as input files:
	
		1. Merged FASTQ files from previous step flash. 
		"""
		
		jobs = []
		input_files = []
		sample_name = []
		catenate_fasta = os.path.join("catenate", "seqs.fna")
		
		for readset in self.readsets:
			merge_directory = os.path.join("merge", readset.sample.name)
			merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")
			
			# Find input readset FASTQs first from previous FLASh job,				
			input_files.append(merge_file_prefix)					         
			sample_name.append(str(readset.sample.name).replace("_", "."))
			

		if config.param('qiime_catenate', 'map_file'):
		
		
			job = qiime.catenate(
			input_files,
			sample_name,
			catenate_fasta
			)
			
			job.name = "catenate"			
			jobs.append(job)

		
			return jobs
			
		else:			
		
			job = qiime.catenate(
			input_files,
			sample_name,
			catenate_fasta
			)
								
			jobs.append(concat_jobs([
				Job(command="python $AMP_SEQ_HOME/AmpliconSeq_script.py -m map_build -i " + ','.join(sample_name)),
				job		
			], name="catenate"))	
			
			return jobs
					
					
	def uchime(self):
		"""
		Reference based chimera detection is performed using [vsearch](https://github.com/torognes/vsearch)

		This step takes as input files:

		1. Catenated FASTA file from previous step catenate. 
		"""
		
		jobs = []
		
		cat_sequence_fasta = os.path.join("catenate", "seqs.fna")
		
		filter_directory = "catenate_without_chimeras"
		filter_fasta = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")	
		filter_log = os.path.join(filter_directory, "seqs_chimeras_filtered.log")
		
		job = qiime.uchime(
			cat_sequence_fasta,
			filter_fasta
		)
		
		job_log = Job([filter_fasta], [filter_log])
		job_log.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m catenate_stat -i {} -j {}""".format(filter_fasta,filter_log)

		jobs.append(concat_jobs([
			Job(command="mkdir -p " + filter_directory),
			job,
			job_log
		], name="uchime"))
					
		return jobs
	
		
	def merge_uchime_stats(self):
		"""
		The chimeric sequences filtered out statistics per readset are merged at this step.
		"""
		
		readset_merge_uchime_stats = os.path.join("metrics", "uchimeReadsetTable.tsv")
		job = concat_jobs([Job(command="mkdir -p metrics"), Job(command="echo 'Sample\tReadset\tMerged Paired Reads #\tFiltered Paired Reads #\tFiltered Paired Reads %' > " + readset_merge_uchime_stats)])
		
		filter_directory = "catenate_without_chimeras"
		filter_log = os.path.join(filter_directory, "seqs_chimeras_filtered.log")	
		
		for readset in self.readsets:
			flash_log = os.path.join("merge", readset.sample.name, readset.name + ".log")
				  
			job = concat_jobs([
				job,
				Job(command="""\
		printf '{}\t{}\t' >> {}""".format(readset.sample.name, readset.name, readset_merge_uchime_stats)
				   )
				  ])
				  			
			# Retrieve merge statistics using re search in python.

			python_command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m uchime -i {} -j {} -s {}""".format(filter_log,flash_log,str(readset.sample.name))
			
			job = concat_jobs([
				job,
				Job(
					[flash_log, filter_log],
					[readset_merge_uchime_stats],
					[
						['qiime', 'module_ampliconseq']
					],
					# Create readset merging stats TSV file with paired read count using python.
					command="""\
{python_command} \\
  >> {readset_merge_uchime_stats}""".format(
						python_command=python_command,
						readset_merge_uchime_stats=readset_merge_uchime_stats
					)
				)
			])

		sample_merge_uchime_stats = os.path.join("metrics", "uchimeSampleTable.tsv")
		report_file = os.path.join("report", "AmpliconSeq.uchime.md")
		return [concat_jobs([
			job,
			Job(
				[readset_merge_uchime_stats],
				[sample_merge_uchime_stats],
				# Create sample trimming stats TSV file with total read counts (i.e. paired * 2 if applicable) using ugly awk
				command="""\
cut -f1,3- {readset_merge_uchime_stats} | awk -F"\t" '{{OFS="\t"; if (NR==1) {{if ($2=="Raw Paired Reads #") {{paired=1}};print "Sample", "Merged Reads #", "Chimera filtered out Reads #", "Chimera filtered out %"}} else {{if (paired) {{$2=$2*2; $3=$3*2}}; raw[$1]+=$2; surviving[$1]+=$3}}}}END{{for (sample in raw){{print sample, raw[sample], surviving[sample], surviving[sample] / raw[sample] * 100}}}}' \\
  > {sample_merge_uchime_stats}""".format(
					readset_merge_uchime_stats=readset_merge_uchime_stats,
					sample_merge_uchime_stats=sample_merge_uchime_stats
				)
			),
			Job(
				[sample_merge_uchime_stats],
				[report_file],
				[['merge_uchime_stats', 'module_pandoc']],
				command="""\
mkdir -p report && \\
cp {readset_merge_uchime_stats} {sample_merge_uchime_stats} report/ && \\
uchime_readset_table_md=`LC_NUMERIC=en_CA awk -F "\t" '{{OFS="|"; if (NR == 1) {{$1 = $1; print $0; print "-----|-----|-----:|-----:|-----:"}} else {{print $1, $2, sprintf("%\\47d", $3), sprintf("%\\47d", $4), sprintf("%.1f", $5)}}}}' {readset_merge_uchime_stats}` && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable read_type={read_type} \\
  --variable sequence_max_n="{sequence_max_n}" \\
  --variable uchime_readset_table="$uchime_readset_table_md" \\
  --to markdown \\
  > {report_file}""".format(
					read_type="Paired",
					report_template_dir=self.report_template_dir,
					sequence_max_n=config.param('qiime_catenate', 'sequence_max_n'),
					readset_merge_uchime_stats=readset_merge_uchime_stats,
					sample_merge_uchime_stats=sample_merge_uchime_stats,
					basename_report_file=os.path.basename(report_file),
					report_file=report_file
				),
				report_files=[report_file]
			)], name="merge_uchime_stats")]
		
	def otu_ref_picking(self):
		"""
		The OTU picking step (close_ref) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

		This step takes as input file:

		1. Catenated and filtered FASTA file from previous step.

		"""
		
		jobs = []
		
		
		filter_directory = "catenate_without_chimeras"
		filter_fastq = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")			
			
		output_directory = "otus/pick_otus"
		otu_file = os.path.join(output_directory, "seqs_chimeras_filtered_otus.txt")
	
		job = qiime.otu_ref_picking(
			filter_fastq,
			output_directory,
			otu_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otus/"),
		job
	], name="qiime_otu_picking"))		
		
		return jobs
		
	def otu_picking(self):
		"""
		The OTU picking step (de novo) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [VSEARCH] (https://github.com/torognes/vsearch) and [Qiime] (http://qiime.org).

		This step takes as input file:

		1. Catenated and filtered FASTA file from previous step.

		"""
		
		jobs = []
		
		
		filter_directory = "catenate_without_chimeras"
		filter_fastq = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")			
			
		output_directory = "otus/pick_otus"
		otu_file = os.path.join(output_directory, "seqs_chimeras_filtered_otus.txt")		
	
		job = qiime.otu_picking(
			filter_fastq,
			output_directory,
			otu_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otus/"),
		job
	], name="qiime_otu_picking"))		
		
		return jobs

	def otu_rep_picking(self):
		"""
		After picking OTUs, this step pick a representative sequence for each OTU.

		This step takes as input files:

		1. OTU file from previous step 
		2. Catenated and filtered FASTA file from filter_chimeras step.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_picking_directory = os.path.join(otu_directory, "pick_otus")		
		otu_file = os.path.join(otu_picking_directory, "seqs_chimeras_filtered_otus.txt")
		
		filter_directory = "catenate_without_chimeras"
		filter_fasta = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")	
								
		output_directory = os.path.join(otu_directory, "pick_rep_set")
		otu_rep_file = os.path.join(output_directory, "rep_set.fna")
	
		job = qiime.otu_rep_picking(
			otu_file,
			filter_fasta,
			output_directory,
			otu_rep_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otus/pick_rep_set"),
		job
	], name="qiime_rep_picking"))		
		
		return jobs


	def otu_assigning(self):
		"""
		Given a set of OTUS, this step attempts to assign the taxonomy of each OTU using [Uclust] (http://drive5.com/usearch/manual/uclust_algo.html).

		This step takes as input files:

		1. OTU representative sequence file from previous step.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_rep_picking_directory = os.path.join(otu_directory, "pick_rep_set")		
		otu_rep_picking_fasta = os.path.join(otu_rep_picking_directory, "rep_set.fna")
			
		output_directory = os.path.join(otu_directory, "taxonomy_assignment")
		tax_assign_file = os.path.join(output_directory, "rep_set_tax_assignments.txt")
	
		job = qiime.otu_assigning(
			otu_rep_picking_fasta,
			output_directory,
			tax_assign_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otus/taxonomy_assignment/"),
		job
	], name="qiime_otu_assigning"))		
		
		return jobs

	def otu_table(self):
		"""
		This step make a consensus OTU table in biom format. It tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU. 

		This step takes as input files:

		1. OTU picking file.
		2. Taxonomy assignment for each OTU from the previous step.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_picking_directory = os.path.join(otu_directory, "pick_otus")		
		otu_file = os.path.join(otu_picking_directory, "seqs_chimeras_filtered_otus.txt")
		
		tax_assign_directory = os.path.join(otu_directory, "taxonomy_assignment")
		tax_assign_file = os.path.join(tax_assign_directory, "rep_set_tax_assignments.txt")
			
		otu_table_file = os.path.join(otu_directory, "otu_table_prefiltered.biom")
		otu_table_sample_file = os.path.join(otu_directory, "otu_table_prefiltered_sample.biom")
		otu_table_final = os.path.join(otu_directory, "otu_table.biom")
		otu_table_summary = os.path.join(otu_directory, "otu_table.sum")
		
		otu_sample_directory = os.path.join(otu_directory, "sample")
		sample_name_control = os.path.join(otu_sample_directory, "done.txt")
	
		job = qiime.otu_table(
			otu_file,
			tax_assign_file,
			otu_directory,
			otu_table_file,
			otu_table_summary
		)
		
		#Remove singleton
		job_filter = Job([otu_table_file], [otu_table_sample_file])
		job_filter.command = """$QIIME_HOME/filter_otus_from_otu_table.py -i {} -n 2 -o {}""".format(otu_table_file,otu_table_sample_file)
		
		#Remove samples with fewer than 2 observations.
		job_filter2 = Job([otu_table_sample_file], [otu_table_final])
		job_filter2.command = """$QIIME_HOME/filter_samples_from_otu_table.py -i {} -n 2 -o {}""".format(otu_table_sample_file,otu_table_final)
		
		#Sample remained after filtering.
		job_sample = Job([otu_table_summary], [sample_name_control], [['qiime', 'module_ampliconseq']])
		job_sample.command ="""python $AMP_SEQ_HOME/AmpliconSeq_script.py -m sample_name -i {} -j {}""".format(otu_table_summary,otu_sample_directory)
		
		jobs.append(concat_jobs([
		job,
		job_filter,
		job_filter2,
		Job(command="""$QIIME_HOME/biom summarize-table -i {} > {}""".format(otu_table_final,otu_table_summary)),
		Job(command="mkdir -p " + otu_sample_directory),
		job_sample
	], name="otu_table"))	
	
		
		return jobs

	def otu_alignment(self):
		"""
		Align OTU representative sequences using [PyNAST] (http://biocore.github.io/pynast/).

		This step takes as input file:

		1. OTU representative sequence file.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_rep_picking_directory = os.path.join(otu_directory, "pick_rep_set")		
		otu_rep_picking_fasta = os.path.join(otu_rep_picking_directory, "rep_set.fna")
			
		output_directory = os.path.join(otu_directory, "align_seq")
		align_seq_fasta= os.path.join(output_directory, "rep_set_aligned.fasta")
	
		job = qiime.otu_alignment(
			otu_rep_picking_fasta,
			output_directory,
			align_seq_fasta
		)
			
		job.name = "qiime_otu_alignment"			
		jobs.append(job)
		
		return jobs

	def filter_alignment(self):
		"""
		Filter the alignment by removing positions which are gaps in every sequence.

		This step takes as input file:

		1. Alignment sequence file.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		align_seq_directory = os.path.join(otu_directory, "align_seq")		
		align_seq_fasta= os.path.join(align_seq_directory, "rep_set_aligned.fasta")
			
		output_directory = os.path.join(otu_directory, "filter_alignment")
		filter_align_fasta = os.path.join(output_directory, "rep_set_aligned_pfiltered.fasta")
	
		job = qiime.filter_alignment(
			align_seq_fasta,
			output_directory,
			filter_align_fasta
		)
			
		job.name = "filter_alignment"			
		jobs.append(job)
		
		return jobs

	def phylogeny(self):
		"""
		Build a phylogenetic tree from a multiple sequence alignment using [FastTree] (http://www.microbesonline.org/fasttree/).

		This step takes as input file:

		1. Filtered alignment sequence file from previous step.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		filter_align_directory = os.path.join(otu_directory, "filter_alignment")		
		filter_align_fasta= os.path.join(filter_align_directory, "rep_set_aligned_pfiltered.fasta")
			
		output_directory = os.path.join(otu_directory, "phylogenetic_tree")
		phylo_file = os.path.join(output_directory, "rep_phylo.tre")
	
		job = qiime.phylogeny(
			filter_align_fasta,
			output_directory,
			phylo_file
		)
			
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otus/phylogenetic_tree/"),
		job
	], name="qiime_phylogeny"))		
		
		return jobs

	def qiime_report(self):
		"""
		1st part report for taxonomic affiliation. 
		"""
		
		jobs = []
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		phylo_directory = os.path.join(otu_directory, "phylogenetic_tree")
		phylo_file = os.path.join(phylo_directory, "rep_phylo.tre")
		report_file = os.path.join("report", "AmpliconSeq.qiime.md")
		
		if config.param('qiime', 'amplicon_type') == '16s':
			amp_db = 'Greengenes'
		elif config.param('qiime', 'amplicon_type') == '18s':
			amp_db = 'Silva'
		elif config.param('qiime', 'amplicon_type') == 'ITS':
			amp_db = 'UNITE'
		else:
			amp_db = 'Unknow'
				
		jobs.append(Job(
                [otu_table, phylo_file],
                [report_file],
                [['qiime', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable amplicon_type="{amplicon_type}" \\
  --variable similarity="{similarity}" \\
  --variable amplicon_db="{amplicon_db}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    amplicon_type=config.param('qiime', 'amplicon_type'),
                    similarity=config.param('qiime_otu_picking', 'similarity'),
                    amplicon_db=amp_db,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="qiime_report")
        )		
		
		return jobs
		
	def multiple_rarefaction(self):
		"""
		1st step (/4) for rarefaction plot.
		Rarefies OTU table by random sampling (without replacement) at different depth in order to perform rarefaction analysis. 
		You need to provide the minimum/maximum number of sequences per samples and the size of each steps between the min/max of seqs/sample. 

		This step takes as input files:

		1. OTU non rarefied table in biom format.

		"""
		
		jobs = []
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		otus_input = [otu_table]
		
		alpha_directory = "alpha_diversity"
		rarefied_otu_directory = os.path.join(alpha_directory, "rarefied_otu_tables")

	
		job = qiime.multiple_rarefaction(
			otus_input,
			rarefied_otu_directory
		)
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p alpha_diversity"),
		job
	], name="multiple_rarefaction"))		
		
		return jobs


	def alpha_diversity(self):
		"""
		2nd step (/4) for rarefaction plot.
		Calculate alpha diversity on each sample using a variety of alpha diversity metrics (chao1, shannon, observed otus). 

		This step takes as input files:

		1. Multiple OTU rarefied table in biom format from previous step.

		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"
		rarefied_otu_directory = os.path.join(alpha_directory, "rarefied_otu_tables")
		
		alpha_diversity_directory = os.path.join(alpha_directory, "alpha_diversity_compute")
		
		job = qiime.alpha_diversity(
			rarefied_otu_directory,
			alpha_diversity_directory
		)
		
		job.name = "alpha_diversity"			
		jobs.append(job)
		
		return jobs

	def collate_alpha(self):
		"""
		3rd step (/4) for rarefaction plot.
		Merge all the alpha diversity computed in the previous step. 
		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"
		alpha_diversity_directory = os.path.join(alpha_directory, "alpha_diversity_compute")
		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
		alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
		chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
		observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
		shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")
		
		job = qiime.collate_alpha(
			alpha_diversity_directory,
			alpha_diversity_collated_merge_directory,
			chao1_stat,
			observed_species_stat,
			shannon_stat
		)
				
		jobs.append(concat_jobs([
					Job(command="mkdir -p alpha_diversity/alpha_diversity_collated"),
					job
				], name="collate_alpha"))

					
		return jobs

	def sample_rarefaction_plot(self):
		"""
		Last step for rarefaction plot.
		Plot the rarefaction curve for each sample
		"""
		
		jobs = []
		
		otu_directory = "otus"
		otu_sample_directory = os.path.join(otu_directory, "sample")
		
		curve_sample=[]
		alpha_directory = "alpha_diversity"		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")	
		alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")	
		alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
		
		chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
		observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
		shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")
			
		for readset in self.readsets:
		
			try:
				self.select_input_files([[os.path.join(otu_sample_directory,str(readset.sample.name).replace("_", ".")+'.txt')]])
				
				sample_collated_general_directory = os.path.join(alpha_diversity_collated_directory, readset.sample.name)
				sample_map = os.path.join(sample_collated_general_directory, "map.txt")
				sample_collated_directory = os.path.join(sample_collated_general_directory, "stat")
				sample_rarefaction_directory = os.path.join(alpha_diversity_rarefaction_directory, readset.sample.name)
				chao1_dir = os.path.join(sample_collated_directory, "chao1.txt")
				observed_species_dir = os.path.join(sample_collated_directory, "observed_species.txt")
				shannon_dir = os.path.join(sample_collated_directory, "shannon.txt")
				observed_species_file = """{}/average_plots/observed_species{}.png""".format(readset.sample.name,str(readset.sample.name).replace('_','.'))
				curve_sample.append(os.path.join(alpha_diversity_rarefaction_directory,observed_species_file))
						
				job = qiime.sample_rarefaction_plot(
					chao1_stat,
					observed_species_stat,
					shannon_stat,
					sample_collated_directory,	
					sample_map,
					sample_rarefaction_directory,
					curve_sample,
				)
				
				jobs.append(concat_jobs([
						Job(command="mkdir -p " + sample_collated_directory),
						Job(command="""python $AMP_SEQ_HOME/AmpliconSeq_script.py -m map_per_sample -s {} -j {}""".format(readset.sample.name,sample_map)),
						Job(command="""python $AMP_SEQ_HOME/AmpliconSeq_script.py -m sample_rarefaction -i {} -j {} -s {}""".format(chao1_stat,chao1_dir,readset.sample.name)),
						Job(command="""python $AMP_SEQ_HOME/AmpliconSeq_script.py -m sample_rarefaction -i {} -j {} -s {}""".format(observed_species_stat,observed_species_dir,readset.sample.name)),
						Job(command="""python $AMP_SEQ_HOME/AmpliconSeq_script.py -m sample_rarefaction -i {} -j {} -s {}""".format(shannon_stat,shannon_dir,readset.sample.name)),
						job
					], name="sample_rarefaction_plot"))
			
			except:
				pass 
				
		return jobs

	def qiime_report2(self):
		"""
		2nd part report for taxonomic affiliation. Plot rarefaction curve for each sample.
		"""
		
		jobs = []
		inputs = []
		curve_sample = []
		
		otu_directory = "otus"
		otu_sample_directory = os.path.join(otu_directory, "sample")
		
		alpha_directory = "alpha_diversity"			
		alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
		
		report_file = os.path.join("report", "AmpliconSeq.plot_curve_no_rar.md")
		num_sample = 0	
		
		for readset in self.readsets:
			try:
				self.select_input_files([[os.path.join(otu_sample_directory,str(readset.sample.name).replace("_", ".")+'.txt')]])
		
				inputs.append(str(readset.sample.name).replace('_','.'))
				observed_species_file = """{}/average_plots/observed_species{}.png""".format(readset.sample.name,str(readset.sample.name).replace('_','.'))
				curve_sample.append(os.path.join(alpha_diversity_rarefaction_directory,observed_species_file))
				num_sample+=1

			except:
				pass 
				
		jobs.append(Job(
                curve_sample,
                [report_file],
                [['qiime', 'module_pandoc']],
                command="""\            
mkdir -p report/fig/alpha_diversity/ && \\
montage -mode concatenate -tile {plot_dimension} {curve_sample} report/fig/alpha_diversity/alpha.rarefaction_sample.png  && \\
pandoc --to=markdown \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    plot_dimension="3x"+str((num_sample/3)+1),
                    curve_sample=' '.join(curve_sample),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="qiime_report2")
        )		
		
		return jobs
		
	def single_rarefaction(self):
		"""
		This step is recommended. It subsamples (rarefy) all the samples to an equal number of sequences for further comparaison.
		You have to provide the number of sequences to subsample per sample in the configuration file (single_rarefaction_depth).

		This step takes as input files:

		1. OTU table in biom format.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		normalization_method = os.path.join(otu_normalized_directory,"rarefaction.txt")
		
		alpha_directory = "alpha_diversity"	
		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
		
		alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
		chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
		observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
		shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")
		
		alpha_diversity_collated_merge_rarefied_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples_rarefied")
		chao1_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "chao1.txt")
		observed_species_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "observed_species.txt")
		shannon_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "shannon.txt")

	
		job = qiime.single_rarefaction(
			otu_table,
			chao1_rarefied_stat,
			observed_species_rarefied_stat,
			shannon_rarefied_stat,
			otu_normalized_table,
			normalization_method
		)
		
		job_chao1 = Job([chao1_stat], [chao1_rarefied_stat])
		job_chao1.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(chao1_stat,chao1_rarefied_stat,config.param('qiime_single_rarefaction', 'single_rarefaction_depth'))
		
		job_observed_species = Job([observed_species_stat], [observed_species_rarefied_stat])
		job_observed_species.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(observed_species_stat,observed_species_rarefied_stat,config.param('qiime_single_rarefaction', 'single_rarefaction_depth'))
		
		job_shannon = Job([shannon_stat], [shannon_rarefied_stat])
		job_shannon.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(shannon_stat,shannon_rarefied_stat,config.param('qiime_single_rarefaction', 'single_rarefaction_depth'))
						
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otu_normalized"),
		Job(command="mkdir -p " + alpha_diversity_collated_merge_rarefied_directory),
		Job(command="touch " + normalization_method),
		job_chao1,
		job_observed_species,
		job_shannon,
		job
	], name="single_rarefaction"))	
	
	
		return jobs

	def css_normalization(self):
		"""
		This step is recommended. Alternative method for normalization to rarefaction. 
		Performs the CSS Matrix normalization.

		This step takes as input files:

		1. OTU table in biom format.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		normalization_method = os.path.join(otu_normalized_directory,"css.txt")
		
		alpha_directory = "alpha_diversity"	
		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
		
		alpha_diversity_collated_merge_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples")
		chao1_stat = os.path.join(alpha_diversity_collated_merge_directory, "chao1.txt")
		observed_species_stat = os.path.join(alpha_diversity_collated_merge_directory, "observed_species.txt")
		shannon_stat = os.path.join(alpha_diversity_collated_merge_directory, "shannon.txt")
		
		alpha_diversity_collated_merge_rarefied_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples_rarefied")
		chao1_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "chao1.txt")
		observed_species_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "observed_species.txt")
		shannon_rarefied_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "shannon.txt")

	
		job = qiime.css_normalization(
			otu_table,
			chao1_rarefied_stat,
			observed_species_rarefied_stat,
			shannon_rarefied_stat,
			otu_normalized_table,
			normalization_method
		)
		
		job_chao1 = Job([chao1_stat], [chao1_rarefied_stat])
		job_chao1.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(chao1_stat,chao1_rarefied_stat,config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max'))
		
		job_observed_species = Job([observed_species_stat], [observed_species_rarefied_stat])
		job_observed_species.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(observed_species_stat,observed_species_rarefied_stat,config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max'))
		
		job_shannon = Job([shannon_stat], [shannon_rarefied_stat])
		job_shannon.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m single_rarefaction -i {} -j {} -s {}""".format(shannon_stat,shannon_rarefied_stat,config.param('qiime_multiple_rarefaction', 'multiple_rarefaction_max'))
						
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p otu_normalized"),
		Job(command="mkdir -p " + alpha_diversity_collated_merge_rarefied_directory),
		Job(command="touch " + normalization_method),
		job_chao1,
		job_observed_species,
		job_shannon,
		job
	], name="css_normalization"))	
		
		return jobs
										  		
	def rarefaction_plot(self):
		"""
		Last step for rarefaction plot.
		Rarefaction curve for each sample on the same plot. 
		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
		alpha_diversity_collated_merge_rarefied_directory = os.path.join(alpha_diversity_collated_directory, "merge_samples_rarefied")
		chao1_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "chao1.txt")
		observed_species_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "observed_species.txt")
		shannon_stat = os.path.join(alpha_diversity_collated_merge_rarefied_directory, "shannon.txt")
		
		alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
		alpha_diversity_rarefaction_rarefied_directory = os.path.join(alpha_diversity_rarefaction_directory, "merge_samples_rarefied")
		alpha_diversity_rarefaction_file = os.path.join(alpha_diversity_rarefaction_rarefied_directory, "rarefaction_plots.html")
		
		if config.param('qiime_catenate', 'map_file'):
			map_file = config.param('qiime_catenate', 'map_file')
		else:
			map_file = "map.txt"
			
		job = qiime.rarefaction_plot(
			alpha_diversity_collated_merge_rarefied_directory,
			chao1_stat,
			observed_species_stat,
			shannon_stat,
			map_file,
			alpha_diversity_rarefaction_file,
			alpha_diversity_rarefaction_rarefied_directory
		)
		
		job.name = "rarefaction_plot"			
		jobs.append(job)
		
		return jobs

	def summarize_taxa(self):
		"""
		1st step (/3) for taxonomic affiliation plot.
		Summarize information of taxonomic groups within each sample at different taxonomic level. 

		This step takes as input files:

		1. OTU rarefied table in biom format if available.
		2. Else, OTU non rarefied table in biom format.

		"""
		
		jobs = []
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		
		alpha_directory = "alpha_diversity"
		taxonomic_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		taxonomic_phylum = os.path.join(taxonomic_directory, "otu_normalized_table_L2.txt")
		taxonomic_class = os.path.join(taxonomic_directory, "otu_normalized_table_L3.txt")
		taxonomic_order = os.path.join(taxonomic_directory, "otu_normalized_table_L4.txt")
		taxonomic_family = os.path.join(taxonomic_directory, "otu_normalized_table_L5.txt")
		taxonomic_genus = os.path.join(taxonomic_directory, "otu_normalized_table_L6.txt")


	
		job = qiime.summarize_taxa(
			otu_normalized_table,
			taxonomic_directory,
			taxonomic_phylum,
			taxonomic_class,
			taxonomic_order,
			taxonomic_family,
			taxonomic_genus
		)
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p alpha_diversity"),
		job
	], name="summarize_taxa"))		
		
		return jobs				

	def plot_taxa(self):
		"""
		2nd step (/3) for taxonomic affiliation plot.
		Make taxaonomy summary bar plots based on taxonomy assignment. 

		This step takes as input files:

		1. Summarized information from previous step.

		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"
		taxonomic_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		taxonomic_phylum = os.path.join(taxonomic_directory, "otu_normalized_table_L2.txt")
		taxonomic_class = os.path.join(taxonomic_directory, "otu_normalized_table_L3.txt")
		taxonomic_order = os.path.join(taxonomic_directory, "otu_normalized_table_L4.txt")
		taxonomic_family = os.path.join(taxonomic_directory, "otu_normalized_table_L5.txt")
		taxonomic_genus = os.path.join(taxonomic_directory, "otu_normalized_table_L6.txt")
		
		taxonomic_input = [taxonomic_phylum, taxonomic_class, taxonomic_order, taxonomic_family, taxonomic_genus]
		
		alpha_diversity_taxonomy_bar_plot = os.path.join(taxonomic_directory, "bar_charts.html")		
	
		job = qiime.plot_taxa(
			taxonomic_input,
			alpha_diversity_taxonomy_bar_plot,
			taxonomic_directory
		)
		
		job.name = "plot_taxa"			
		jobs.append(job)	
		
		return jobs	

	def plot_heatmap(self):
		"""
		Last step for taxonomic affiliation plot.
		Make heatmap at phylum level. 

		This step takes as input files:

		1. Summarized information from previous step.

		"""
		
		jobs = []
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		
		alpha_directory = "alpha_diversity"
		taxonomic_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		taxonomic_phylum = os.path.join(taxonomic_directory, "otu_normalized_table_L2.txt")
		
		beta_directory = "beta_diversity"
		heatmap_directory = os.path.join(beta_directory, "heatmap")
		
		heatmap_script = os.path.join(heatmap_directory, "OTU_Phylum_to_R.R")
		heatmap_chart = os.path.join(heatmap_directory, "otu_heatmap.png")
		
		heatmap_otu_data_R = os.path.join(heatmap_directory, "OTU_data.txt")
		heatmap_otu_name_R = os.path.join(heatmap_directory, "OTU_name.txt")
		heatmap_otu_tax_R = os.path.join(heatmap_directory, "OTU_tax_final.txt")
		
		heatmap_otu_table = os.path.join(heatmap_directory, "otumat.tsv")	
		heatmap_tax_table = os.path.join(heatmap_directory, "taxmat.tsv")			
				
		taxonomic_input = [taxonomic_phylum]
		
		job = Job(taxonomic_input, [heatmap_script,heatmap_otu_data_R,heatmap_otu_name_R,heatmap_otu_tax_R], [['qiime', 'module_R'],['qiime', 'module_ampliconseq']])
		job.command = """python $AMP_SEQ_HOME/AmpliconSeq_script.py -m plot_heatmap -i {} -j {} -s {}""".format(taxonomic_phylum,heatmap_directory,1)
		
		jobR = Job([heatmap_script,heatmap_otu_data_R,heatmap_otu_name_R,heatmap_otu_tax_R], [heatmap_chart,heatmap_otu_table,heatmap_tax_table], [['qiime', 'module_R'],['qiime', 'module_ampliconseq']])
		jobR.command = "./beta_diversity/heatmap/OTU_Phylum_to_R.R"
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p beta_diversity/heatmap/"),
		job,
		Job(command="chmod +x beta_diversity/heatmap/OTU_Phylum_to_R.R"), 
		jobR
	], name="plot_heatmap"))	
	
		return jobs									

	def krona(self):
		"""
		Plot Krona chart for taxonomic affiliation
		"""
		
		jobs = []
		sample_name =[]
		
		otu_directory = "otus"
		otu_sample_directory = os.path.join(otu_directory, "sample")
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		
		alpha_directory = "alpha_diversity"		
		alpha_diversity_krona_directory = os.path.join(alpha_directory, "krona_chart")
		alpha_diversity_krona_file = os.path.join(alpha_diversity_krona_directory, "krona_chart.html")
		
		for readset in self.readsets:		
			try:
				self.select_input_files([[os.path.join(otu_sample_directory,str(readset.sample.name).replace("_", ".")+'.txt')]])	
				sample_name.append(alpha_diversity_krona_directory+'/'+str(readset.sample.name).replace("_", ".")+'.txt') 
			except:
				pass 
						
		job = qiime.krona(
			otu_normalized_table,
			sample_name,
			alpha_diversity_krona_file,
		)
		
		jobs.append(concat_jobs([
				# Create an output directory
				Job(command="mkdir -p alpha_diversity/krona_chart"),
				Job(command='$QIIME_HOME/biom convert -i {} -o alpha_diversity/table_tax.txt --table-type="OTU table" --to-tsv --header-key taxonomy'.format(otu_normalized_table)),
				Job(command="python $AMP_SEQ_HOME/AmpliconSeq_script.py -m krona -i alpha_diversity/table_tax.txt"),
				job
			], name="krona"))
		return jobs

	def plot_to_alpha(self):
		"""
		Final report 1st part for the Amplicon-Seq pipeline. Display results (taxonomy, heatmap and alpha diversity).
		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"	
				
		alpha_diversity_taxonomy_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		alpha_diversity_taxonomy_bar_plot = os.path.join(alpha_diversity_taxonomy_directory, "bar_charts.html")
		
		alpha_diversity_krona_directory = os.path.join(alpha_directory, "krona_chart")
		alpha_diversity_krona_directory = os.path.join(alpha_diversity_krona_directory, "krona_chart.html")
		
		alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
		alpha_diversity_rarefaction_file = os.path.join(alpha_diversity_rarefaction_directory, "merge_samples_rarefied/rarefaction_plots.html")
		
		beta_directory = "beta_diversity"	
		
		beta_diversity_heatmap_directory = os.path.join(beta_directory, "heatmap")
		beta_diversity_heatmap_plot = os.path.join(beta_diversity_heatmap_directory, "otu_heatmap.png")
		beta_diversity_heatmap_otumat = os.path.join(beta_diversity_heatmap_directory, "otumat.tsv")
		beta_diversity_heatmap_taxmat = os.path.join(beta_diversity_heatmap_directory, "taxmat.tsv")
		
		otu_normalized_directory = "otu_normalized"
		rarefaction_method = os.path.join(otu_normalized_directory,"rarefaction.txt")
		css_method = os.path.join(otu_normalized_directory,"css.txt")
		
		candidate_input_files = [[rarefaction_method]]
		candidate_input_files.append([css_method])	
		normalization_method = self.select_input_files(candidate_input_files)
						
		inputs = [alpha_diversity_taxonomy_bar_plot,alpha_diversity_krona_directory,alpha_diversity_rarefaction_file,beta_diversity_heatmap_plot,beta_diversity_heatmap_otumat,beta_diversity_heatmap_taxmat,normalization_method[0]]
		
		if normalization_method == ['otu_normalized/rarefaction.txt']:			
			report_file = os.path.join("report", "AmpliconSeq.plot_to_alpha_rar.md")	
		else:
			report_file = os.path.join("report", "AmpliconSeq.plot_to_alpha_css.md")


				
		jobs.append(Job(
                inputs,
                [report_file],
                [['qiime', 'module_pandoc']],
                command="""\            
mkdir -p report/fig/alpha_diversity/ && \\
mkdir -p report/fig/beta_diversity/heatmap/ && \\
cp -r alpha_diversity/taxonomic_affiliation/ report/fig/alpha_diversity/taxonomic_affiliation/ && \\
cp -r alpha_diversity/krona_chart/ report/fig/alpha_diversity/krona_chart/ && \\
cp -r alpha_diversity/alpha_rarefaction/merge_samples_rarefied/ report/fig/alpha_diversity/alpha_rarefaction/ && \\
cp beta_diversity/heatmap/otu_heatmap.png report/fig/beta_diversity/heatmap/otu_heatmap.png && \\
cp beta_diversity/heatmap/otumat.tsv report/fig/beta_diversity/heatmap/otumat.tsv && \\
cp beta_diversity/heatmap/taxmat.tsv report/fig/beta_diversity/heatmap/taxmat.tsv && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable single_rarefaction_depth="{single_rarefaction_depth}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    single_rarefaction_depth=config.param('qiime_single_rarefaction', 'single_rarefaction_depth'),  
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="plot_to_alpha")
        )		
		
		return jobs
		
	def beta_diversity(self):
		"""
		1st step (/3) for 2D PCoA plot.
		Calculate beta diversity (pairwise sample dissimilarity) on OTU table. The OTU table has to be normalized. 
		Only works with >= 4 samples

		This step takes as input files:

		1. OTU rarefied table in biom format.
		2. Tree file.

		"""
		
		jobs = []
		
		otu_normalized_directory = "otu_normalized"
		otu_normalized_table = os.path.join(otu_normalized_directory,"otu_normalized_table.biom")
		
		otu_directory = "otus"
		phylogenetic_tree_directory = os.path.join(otu_directory, "phylogenetic_tree")		
		phylogenetic_tree_file= os.path.join(phylogenetic_tree_directory, "rep_phylo.tre")
		
		beta_diversity_directory = "beta_diversity"
		dm_directory = os.path.join(beta_diversity_directory, "dissimilarity_matrix")
		dm_unweighted_file = os.path.join(dm_directory, "unweighted_unifrac_otu_normalized_table.txt")
		dm_weighted_file = os.path.join(dm_directory, "weighted_unifrac_otu_normalized_table.txt")
	
		job = qiime.beta_diversity(
			otu_normalized_table,
			phylogenetic_tree_file,
			dm_directory,
			dm_unweighted_file,
			dm_weighted_file
		)
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p beta_diversity/dissimilarity_matrix/"),
		job
	], name="beta_diversity"))		
		
		return jobs

	def pcoa(self):
		"""
		2nd step (/3) for 2D PCoA plot.
		Compute coordinates pour PCoA 

		This step takes as input file:

		1. Matrix produced in the previous step.

		"""
		
		jobs = []
		
		beta_diversity_directory = "beta_diversity"
		dm_directory = os.path.join(beta_diversity_directory, "dissimilarity_matrix")
		dm_unweighted_file = os.path.join(dm_directory, "unweighted_unifrac_otu_normalized_table.txt")
		dm_weighted_file = os.path.join(dm_directory, "weighted_unifrac_otu_normalized_table.txt")
		
		pcoa_directory = os.path.join(beta_diversity_directory, "principal_coordinates")
		pcoa_unweighted_file = os.path.join(pcoa_directory, "pcoa_unweighted_unifrac_otu_normalized_table.txt")
		pcoa_weighted_file = os.path.join(pcoa_directory, "pcoa_weighted_unifrac_otu_normalized_table.txt")
			
		job = qiime.pcoa(
			dm_unweighted_file,
			dm_weighted_file,
			dm_directory,
			pcoa_directory,
			pcoa_unweighted_file,
			pcoa_weighted_file
		)
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p beta_diversity/principal_coordinates/"),
		job
	], name="pcoa"))		
		
		return jobs

	def pcoa_plot(self):
		"""
		Last step for 2D PCoA plot.

		This step takes as input file:

		1. PCoA from the previous step.

		"""
		
		jobs = []
		
		beta_diversity_directory = "beta_diversity"
		pcoa_directory = os.path.join(beta_diversity_directory, "principal_coordinates")
		pcoa_unweighted_file = os.path.join(pcoa_directory, "pcoa_unweighted_unifrac_otu_normalized_table.txt")
		pcoa_weighted_file = os.path.join(pcoa_directory, "pcoa_weighted_unifrac_otu_normalized_table.txt")

		pcoa_plot_directory = os.path.join(beta_diversity_directory, "2d_plots")	
		beta_diversity_pcoa_unweighted = os.path.join(pcoa_plot_directory, "pcoa_unweighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
		beta_diversity_pcoa_weighted = os.path.join(pcoa_plot_directory, "pcoa_weighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
		
		if config.param('qiime_catenate', 'map_file'):
			map_file = config.param('qiime_catenate', 'map_file')
		else:
			map_file = "map.txt"
							
		job1 = qiime.pcoa_plot(
			pcoa_unweighted_file,
			pcoa_directory,
			map_file,
			beta_diversity_pcoa_unweighted,
			pcoa_plot_directory
		)

		job2 = qiime.pcoa_plot(
			pcoa_weighted_file,
			pcoa_directory,
			map_file,
			beta_diversity_pcoa_weighted,
			pcoa_plot_directory
		)		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p beta_diversity/2d_plots/"),
		job1,
		job2
	], name="pcoa_plot"))		
		
		return jobs

		
	def plot_to_beta(self):
		"""
		Final report's 2nd part for the Amplicon-Seq pipeline. Display results (beta diversity PCoA plots).
		"""
		
		jobs = []
		
		beta_directory = "beta_diversity"	
			
		beta_diversity_pcoa_directory = os.path.join(beta_directory, "2d_plots")
		beta_diversity_pcoa_unweighted = os.path.join(beta_diversity_pcoa_directory, "pcoa_unweighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
		beta_diversity_pcoa_weighted = os.path.join(beta_diversity_pcoa_directory, "pcoa_weighted_unifrac_otu_normalized_table_2D_PCoA_plots.html")
		
		inputs = [beta_diversity_pcoa_unweighted,beta_diversity_pcoa_weighted]
				
		report_file = os.path.join("report", "AmpliconSeq.plot_to_beta.md")

				
		jobs.append(Job(
                inputs,
                [report_file],
                [['qiime', 'module_pandoc']],
                command="""\            
mkdir -p report/fig/beta_diversity/ && \\
cp -r beta_diversity/2d_plots/ report/fig/beta_diversity/2d_plots/ && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="plot_to_beta")
        )		
		
		return jobs
																
				
	@property
	def steps(self):
		return [
			self.trimmomatic,
			self.merge_trimmomatic_stats,
			self.flash,
			self.merge_flash_stats,
			self.catenate,	#5
			self.uchime,
			self.merge_uchime_stats,
			self.otu_ref_picking,
			self.otu_picking,	
			self.otu_rep_picking,	#10
			self.otu_assigning,	
			self.otu_table,	
			self.otu_alignment,	
			self.filter_alignment,
			self.phylogeny,	#15	
			self.qiime_report,	
			self.multiple_rarefaction,
			self.alpha_diversity,	
			self.collate_alpha,
			self.sample_rarefaction_plot,	#20	
			self.qiime_report2,	
			self.single_rarefaction,
			self.css_normalization,
			self.rarefaction_plot,	
			self.summarize_taxa,	#25	
			self.plot_taxa,	
			self.plot_heatmap,	
			self.krona,
			self.plot_to_alpha,	
			self.beta_diversity,	#30
			self.pcoa,
			self.pcoa_plot,	
			self.plot_to_beta
		]

if __name__ == '__main__': 
	AmpliconSeq()
