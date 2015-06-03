#!/usr/bin/env python

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
from bfx.readset import *
from bfx.sequence_dictionary import *

from pipelines import common
from bfx import flash
from bfx import trimmomatic
from bfx import qiime


log = logging.getLogger(__name__)

class Metagenomics(common.Illumina):
	"""
	Metagenomics Pipeline
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
								    
			elif readset.run_type == "SINGLE_END":
				candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
				if readset.fastq1:
					candidate_input_files.append([readset.fastq1])
				[fastq1] = self.select_input_files(candidate_input_files)
				fastq2 = None
				
			else:
				raise Exception("Error: run type \"" + readset.run_type +
				"\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

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
		output_file = os.path.join("catenate", "seqs.fna")
		
		for readset in self.readsets:
			merge_directory = os.path.join("merge", readset.sample.name)
			merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")
			
			# Find input readset FASTQs first from previous FLASh job,				
			input_files.append(merge_file_prefix)					         
			sample_name.append(str(readset.sample.name).replace("_", "."))
		
		job = qiime.catenate(
			input_files,
			sample_name,
			output_file
		)
		
		job.name = "catenate"			
		jobs.append(job)
		
		return jobs
					
	def uchime(self):
		"""
		Reference based chimera detection is performed using [Usearch/Uchime](http://drive5.com/usearch/) (http://drive5.com/usearch/manual/uchime_algo.html).

		This step takes as input files:

		1. Catenated FASTA file from previous step catenate. 
		"""
		
		jobs = []
		
		cat_file_prefix = os.path.join("catenate", "seqs.fna")
		#candidate_input_files = [[cat_file_prefix]]
		output_file = os.path.join("usearch_checked_chimeras", "chimeras.txt")
				
		job = qiime.uchime(
			#self.select_input_files(candidate_input_files)[0],
			cat_file_prefix,
			output_file
		)
		
		job.name = "uchime"			
		jobs.append(job)
		
		return jobs



	def filter_chimeras(self):
		"""
		Filter the input catenate data file by passing the chimeras file created in the previous step.

		This step takes as input files:

		1. Catenated FASTQ file from 3rd step cat_data.
		2. Chimera file from previous step uchime.
		"""
		
		jobs = []
		
		cat_file_prefix = os.path.join("catenate", "seqs.fna")
		chimeras_file_prefix = os.path.join("usearch_checked_chimeras", "chimeras.txt")
		candidate_input_files = [[cat_file_prefix]]
		candidate_input_files.append([chimeras_file_prefix])
		
		filter_directory = "catenate_without_chimeras"
		filter_fastq = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")			
		
		job = qiime.filter_chimeras(
			candidate_input_files[0][0],
			candidate_input_files[1][0],
			filter_fastq
		)
		
		jobs.append(concat_jobs([
				# Create an output directory
				Job(command="mkdir -p " + filter_directory),
				job
			], name="filter_chimeras"))
		return jobs
		
		

	def otu_picking(self):
		"""
		The OTU picking step (de novo) assigns similar sequences to operational taxonomic units (OTUs) by clustering sequences based on a user-defined similarity threshold. Method per default uses [usearch61] (http://drive5.com/usearch/) program wrapped by [Qiime] (http://qiime.org).

		This step takes as input file:

		1. Catenated and filtered FASTA file from previous step.

		"""
		
		jobs = []
		
		
		filter_directory = "catenate_without_chimeras"
		filter_fastq = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")			
		#candidate_input_files = [[filter_fastq]]
			
		output_directory = "otus/pick_otus"
		otu_file = os.path.join(output_directory, "seqs_chimeras_filtered_otus.txt")
		#print self.select_input_files(candidate_input_files)
		
	
		job = qiime.otu_picking(
			filter_fastq,
			output_directory,
			otu_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir otus/"),
		job
	], name="otu_picking"))		
		
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
					
		#candidate_input_files = [[filter_fastq]]
			
		output_directory = os.path.join(otu_directory, "pick_rep_set")
		otu_rep_file = os.path.join(output_directory, "rep_set.fna")
		#print self.select_input_files(candidate_input_files)
	
		job = qiime.otu_rep_picking(
			otu_file,
			filter_fasta,
			output_directory,
			otu_rep_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir otus/pick_rep_set"),
		job
	], name="otu_rep_picking"))		
		
		return jobs


	def otu_assigning(self):
		"""
		Given a set of OTUS, this step attempts to assign the taxonomy of each OTU using [uclust] (http://drive5.com/usearch/manual/uclust_algo.html).

		This step takes as input files:

		1. OTU representative sequence file from previous step.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_rep_picking_directory = os.path.join(otu_directory, "pick_rep_set")		
		otu_rep_picking_fasta = os.path.join(otu_rep_picking_directory, "rep_set.fna")

		#candidate_input_files = [[filter_fastq]]
			
		output_directory = os.path.join(otu_directory, "taxonomy_assignment")
		tax_assign_file = os.path.join(output_directory, "rep_set_tax_assignments.txt")
		#print self.select_input_files(candidate_input_files)
	
		job = qiime.otu_assigning(
			otu_rep_picking_fasta,
			output_directory,
			tax_assign_file
		)
		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir otus/taxonomy_assignment/"),
		job
	], name="otu_assigning"))		
		
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

		#candidate_input_files = [[filter_fastq]]
			
		otu_table_file = os.path.join(otu_directory, "otu_table.biom")
		#print self.select_input_files(candidate_input_files)
	
		job = qiime.otu_table(
			otu_file,
			tax_assign_file,
			otu_directory,
			otu_table_file
		)
		
		job.name = "otu_table"			
		jobs.append(job)
		
		return jobs


				
				
		
	@property
	def steps(self):
		return [
			self.trimmomatic,
			self.merge_trimmomatic_stats,
			self.flash,
			self.merge_flash_stats,
			self.catenate,	
			self.uchime,
			self.filter_chimeras,
			self.otu_picking,
			self.otu_rep_picking,
			self.otu_assigning,
			self.otu_table
		]

if __name__ == '__main__': 
	Metagenomics()
