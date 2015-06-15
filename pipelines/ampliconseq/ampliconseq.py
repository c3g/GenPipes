#!/usr/bin/env python

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
		output_file = os.path.join("catenate", "seqs.fna")
		
		for readset in self.readsets:
			merge_directory = os.path.join("merge", readset.sample.name)
			merge_file_prefix = os.path.join(merge_directory, readset.name + ".extendedFrags.fastq")
			
			# Find input readset FASTQs first from previous FLASh job,				
			input_files.append(merge_file_prefix)					         
			sample_name.append(str(readset.sample.name).replace("_", "."))

		if config.param('qiime', 'map_file'):
		
		
			job = qiime.catenate(
			input_files,
			sample_name,
			output_file
			)
		
			job.name = "catenate"			
			jobs.append(job)
			
			return jobs
			
		else:			
			job = qiime.catenate(
			input_files,
			sample_name,
			output_file
			)
								
			jobs.append(concat_jobs([
				# Create an output directory
				Job(command="python $MUGQIC_INSTALL_HOME/software/AmpliconSeq/AmpliconSeq_script.py -m map_build -i " + ','.join(sample_name)),
				job
			], name="catenate"))	
			
			return jobs
					
	def uchime(self):
		"""
		Reference based chimera detection is performed using [Usearch/Uchime](http://drive5.com/usearch/) (http://drive5.com/usearch/manual/uchime_algo.html).

		This step takes as input files:

		1. Catenated FASTA file from previous step catenate. 
		"""
		
		jobs = []
		
		cat_sequence_fasta = os.path.join("catenate", "seqs.fna")
		output_directory = "usearch_checked_chimeras"
		chimeras_file = os.path.join(output_directory, "chimeras.txt")
				
		job = qiime.uchime(
			cat_sequence_fasta,
			output_directory,
			chimeras_file
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
		
		cat_sequence_fasta = os.path.join("catenate", "seqs.fna")
		chimeras_file = os.path.join("usearch_checked_chimeras", "chimeras.txt")
		
		filter_directory = "catenate_without_chimeras"
		filter_fasta = os.path.join(filter_directory, "seqs_chimeras_filtered.fna")			
		
		job = qiime.filter_chimeras(
			cat_sequence_fasta,
			chimeras_file,
			filter_fasta
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
			
		output_directory = "otus/pick_otus"
		otu_file = os.path.join(output_directory, "seqs_chimeras_filtered_otus.txt")		
	
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
			
		output_directory = os.path.join(otu_directory, "taxonomy_assignment")
		tax_assign_file = os.path.join(output_directory, "rep_set_tax_assignments.txt")
	
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
			
		otu_table_file = os.path.join(otu_directory, "otu_table.biom")
	
		job = qiime.otu_table(
			otu_file,
			tax_assign_file,
			otu_directory,
			otu_table_file
		)
		
		job.name = "otu_table"			
		jobs.append(job)
		
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
			
		job.name = "otu_alignment"			
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
		Build a phylogenetic tree from a multiple sequence alignment.

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
	], name="phylogeny"))		
		
		return jobs
					

	def single_rarefaction(self):
		"""
		This step is optionnal. It subsamples (rarefy) all the samples to an equal number of sequences for further comparaison.
		You have to provide the number of sequences to subsample per sample in the configuration file (single_rarefaction_depth).

		This step takes as input files:

		1. OTU table in biom format.

		"""
		
		jobs = []
				
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")

	
		job = qiime.single_rarefaction(
			otu_table,
			otu_even_table
		)
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir otus_even"),
		job
	], name="single_rarefaction"))	
	
	
		return jobs


	def multiple_rarefaction(self):
		"""
		1st step (/4) for rarefaction plot.
		Rarefies OTU table by random sampling (without replacement) at different depth in order to perform rarefaction analysis. 
		You need to provide the minimum/maximum number of sequences per samples and the size of each steps between the min/max of seqs/sample. 

		This step takes as input files:

		1. OTU rarefied table in biom format if available.
		2. Else, OTU non rarefied table in biom format.

		"""
		
		jobs = []
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		candidate_input_files = [[otu_even_table]]
		candidate_input_files.append([otu_table])
		
		otus_input = self.select_input_files(candidate_input_files)
		
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
		Calculate alpha diversity on each sample using a variety of alpha diversity metrics (chao1, observed otus). 

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
		
		job = qiime.collate_alpha(
			alpha_diversity_directory,
			alpha_diversity_collated_directory
		)
		
		job.name = "collate_alpha"			
		jobs.append(job)
		
		return jobs

	def rarefaction_plot(self):
		"""
		Last step for rarefaction plot.
		Plot the rarefaction curve. 
		"""
		
		jobs = []
		
		alpha_directory = "alpha_diversity"		
		alpha_diversity_collated_directory = os.path.join(alpha_directory, "alpha_diversity_collated")
		
		alpha_diversity_rarefaction_directory = os.path.join(alpha_directory, "alpha_rarefaction")
		
		if config.param('qiime', 'map_file'):
			map_file = config.param('qiime', 'map_file')
		else:
			map_file = "map.txt"
			
		job = qiime.rarefaction_plot(
			alpha_diversity_collated_directory,
			map_file,
			alpha_diversity_rarefaction_directory
		)
		
		job.name = "rarefaction_plot"			
		jobs.append(job)
		
		return jobs

	def summarize_taxa(self):
		"""
		1st step (/2) for taxonomic affiliation plot.
		Summarize information of taxonomic groups within each sample at different taxonomic level. 

		This step takes as input files:

		1. OTU rarefied table in biom format if available.
		2. Else, OTU non rarefied table in biom format.

		"""
		
		jobs = []
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		candidate_input_files = [[otu_even_table]]
		candidate_input_files.append([otu_table])
		
		otus_input = self.select_input_files(candidate_input_files)
		
		alpha_directory = "alpha_diversity"
		taxonomic_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		taxonomic_phylum = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L2.txt")
		taxonomic_class = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L3.txt")
		taxonomic_order = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L4.txt")
		taxonomic_family = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L5.txt")
		taxonomic_genus = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L6.txt")


	
		job = qiime.summarize_taxa(
			otus_input,
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
		Last step for taxonomic affiliation plot.
		Make taxaonomy summary charts based on taxonomy assignment. 

		This step takes as input files:

		1. Summarized information from previous step.

		"""
		
		jobs = []
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		candidate_input_files = [[otu_even_table]]
		candidate_input_files.append([otu_table])
		
		otus_input = self.select_input_files(candidate_input_files)
		
		alpha_directory = "alpha_diversity"
		taxonomic_directory = os.path.join(alpha_directory, "taxonomic_affiliation")
		taxonomic_phylum = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L2.txt")
		taxonomic_class = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L3.txt")
		taxonomic_order = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L4.txt")
		taxonomic_family = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L5.txt")
		taxonomic_genus = os.path.join(taxonomic_directory, os.path.splitext(basename(otus_input[0]))[0]+"_L6.txt")
		
		taxonomic_input = [taxonomic_phylum, taxonomic_class, taxonomic_order, taxonomic_family, taxonomic_genus]		
	
		job = qiime.plot_taxa(
			taxonomic_input,
			taxonomic_directory
		)
		
		job.name = "plot_taxa"			
		jobs.append(job)	
		
		return jobs	
						

	def krona(self):
		"""
		Plot Krona chart for taxonomic affiliation
		"""
		
		jobs = []
		sample_name =[]
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")
		
		otu_directory = "otus"
		otu_table = os.path.join(otu_directory, "otu_table.biom")
		
		candidate_input_files = [[otu_even_table]]
		candidate_input_files.append([otu_table])
		
		otus_input = self.select_input_files(candidate_input_files)
		
		alpha_directory = "alpha_diversity"		
		alpha_diversity_krona_directory = os.path.join(alpha_directory, "krona_chart")
		alpha_diversity_krona_file = os.path.join(alpha_diversity_krona_directory, "krona_chart.html")
		
		
		for readset in self.readsets:				         
			sample_name.append(alpha_diversity_krona_directory+'/'+str(readset.sample.name).replace("_", ".")+'.txt')
			
		job = qiime.krona(
			otus_input,
			sample_name,
			alpha_diversity_krona_file,
		)
		
		jobs.append(concat_jobs([
				# Create an output directory
				Job(command="mkdir -p alpha_diversity/krona_chart"),
				Job(command='$QIIME_HOME/biom convert -i {} -o alpha_diversity/table_tax.txt --table-type="OTU table" --to-tsv --header-key taxonomy'.format(otus_input[0])),
				Job(command="python $MUGQIC_INSTALL_HOME/software/AmpliconSeq/AmpliconSeq_script.py -m krona -i alpha_diversity/table_tax.txt"),
				job
			], name="krona"))
		return jobs

	def beta_diversity(self):
		"""
		1st step (/3) for 2D PCoA plot.
		Calculate beta diversity (pairwise sample dissimilarity) on OTU table. The OTU table has to be rarefied. 

		This step takes as input files:

		1. OTU rarefied table in biom format.
		2. Tree file.

		"""
		
		jobs = []
		
		otu_even_directory = "otus_even"
		otu_even_table = os.path.join(otu_even_directory,"otu_even_table.biom")
		
		otu_directory = "otus"
		phylogenetic_tree_directory = os.path.join(otu_directory, "phylogenetic_tree")		
		phylogenetic_tree_file= os.path.join(phylogenetic_tree_directory, "rep_phylo.tre")
		
		beta_diversity_directory = "beta_diversity"
		dm_directory = os.path.join(beta_diversity_directory, "dissimilarity_matrix")
		dm_unweighted_file = os.path.join(dm_directory, "unweighted_unifrac_otu_even_table.txt")
		dm_weighted_file = os.path.join(dm_directory, "weighted_unifrac_otu_even_table.txt")
	
		job = qiime.beta_diversity(
			otu_even_table,
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
		dm_unweighted_file = os.path.join(dm_directory, "unweighted_unifrac_otu_even_table.txt")
		dm_weighted_file = os.path.join(dm_directory, "weighted_unifrac_otu_even_table.txt")
		
		pcoa_directory = os.path.join(beta_diversity_directory, "principal_coordinates")
	
		job = qiime.pcoa(
			dm_unweighted_file,
			dm_weighted_file,
			dm_directory,
			pcoa_directory
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
		pcoa_unweighted_file = os.path.join(pcoa_directory, "pcoa_unweighted_unifrac_otu_even_table.txt")
		pcoa_weighted_file = os.path.join(pcoa_directory, "pcoa_weighted_unifrac_otu_even_table.txt")

		pcoa_plot_directory = os.path.join(beta_diversity_directory, "2d_plots")	
		
		if config.param('qiime', 'map_file'):
			map_file = config.param('qiime', 'map_file')
		else:
			map_file = "map.txt"
							
		job1 = qiime.pcoa_plot(
			pcoa_unweighted_file,
			pcoa_directory,
			map_file,
			pcoa_plot_directory
		)

		job2 = qiime.pcoa_plot(
			pcoa_weighted_file,
			pcoa_directory,
			map_file,
			pcoa_plot_directory
		)		
		
		jobs.append(concat_jobs([
		# Create an output directory
		Job(command="mkdir -p beta_diversity/2d_plots/"),
		job1,
		job2
	], name="pcoa_plot"))		
		
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
			self.otu_assigning,	#10
			self.otu_table,
			self.otu_alignment,
			self.filter_alignment,
			self.phylogeny,
			self.single_rarefaction,	#15
			self.multiple_rarefaction,
			self.alpha_diversity,
			self.collate_alpha,
			self.rarefaction_plot,	
			self.summarize_taxa,	#20
			self.plot_taxa,
			self.krona,
			self.beta_diversity,
			self.pcoa,
			self.pcoa_plot
		]

if __name__ == '__main__': 
	AmpliconSeq()
