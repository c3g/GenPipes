#!/usr/bin/env python

# Python Standard Modules
import logging
import os

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)

def catenate(
	input_fastq,
	input_name,
	output_fasta
	):

	inputs = input_fastq
	outputs = [output_fasta]

	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/split_libraries_fastq.py \\
  -i {input_files} \\
  --sample_id {sample_name} \\
  -o {dir_output} \\
  -r 30 \\
  -p 0 \\
  -n 100 \\
  --barcode_type 'not-barcoded'""".format(
		input_files=','.join(input_fastq),
		sample_name=','.join(input_name),
		dir_output="catenate/",
		),
		removable_files=[output_fasta]
	)

def uchime(
	input_fasta,
	output_file
	):

	inputs = [input_fasta]
	outputs = [output_file]

	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime'],
			['qiime', 'module_usearch61']
		],

		command="""\
  $QIIME_HOME/identify_chimeric_seqs.py \\
  -i {input_file} \\
  -m {usearch61} \\
  -o {dir_output} \\
  -r {database}""".format(
		input_file=input_fasta,
		usearch61="usearch61",
		dir_output="usearch_checked_chimeras/",
		database=config.param('qiime', 'chimera_database')
		),
		removable_files=[output_file]
	)
	
def filter_chimeras(
	input_fasta,
	chimera_txt,
	output_fasta
	):

	inputs = [input_fasta, chimera_txt]
	outputs = [output_fasta]

	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/filter_fasta.py \\
  -f {input_file} \\
  -o {output_file} \\
  -s {chimera_file} \\
  -n""".format(
		input_file=input_fasta,
		output_file=output_fasta,
		chimera_file=chimera_txt
		),
		removable_files=[output_fasta]
	)

def otu_picking(
	input_without_chimer,
	output_directory,
	output_otus
	):

	inputs = [input_without_chimer]
	outputs = [output_otus]
	
	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/pick_otus.py \\
  -i {input_without_chimer} \\
  -o {output_directory} \\
  -m {method} \\
  -s {similarity_treshold} \\
  --threads {threads_number}""".format(
		input_without_chimer=input_without_chimer,
		output_directory=output_directory,
		method='usearch61',
		similarity_treshold=config.param('qiime', 'similarity'),
		threads_number=config.param('qiime', 'threads')
		),
		removable_files=[output_otus]
	)

		
def otu_rep_picking(
	otu_file,
	filter_fasta,
	output_directory,
	otu_rep_file
	):

	inputs = [otu_file, filter_fasta]
	outputs = [otu_rep_file]
	
	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/pick_rep_set.py \\
  -i {otu_file} \\
  -f {filter_fasta} \\
  -m {method} \\
  -o {output_directory}""".format(
		otu_file=otu_file,
		filter_fasta=filter_fasta,
		method=config.param('qiime', 'rep_set_picking_method'),
		output_directory=otu_rep_file
		),
		removable_files=[otu_rep_file]
	)
			
def otu_assigning(
	otu_rep_picking_fasta,
	output_directory,
	tax_assign_file
	):

	inputs = [otu_rep_picking_fasta]
	outputs = [tax_assign_file]
	
	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/assign_taxonomy.py \\
  -i {otu_rep_picking_fasta} \\
  -r {database_otus} \\
  -t {taxonomy_otus} \\
  -o {output_directory}""".format(
		otu_rep_picking_fasta=otu_rep_picking_fasta,
		database_otus=config.param('qiime', 'reference_seqs_fp'),
		taxonomy_otus=config.param('qiime', 'id_to_taxonomy_fp'),
		output_directory=output_directory
		),
		removable_files=[tax_assign_file]
	)
			
				
def otu_table(
	otu_file,
	tax_assign_file,
	otu_directory,
	otu_table_file
	):

	inputs = [otu_file, tax_assign_file]
	outputs = [otu_table_file]
	
	return Job(
		inputs,
		outputs,
		[
			['qiime', 'module_qiime']
		],

		command="""\
  $QIIME_HOME/make_otu_table.py \\
  -i {otu_file} \\
  -t {tax_assign_file} \\
  -o {otu_table_file}""".format(
		otu_file=otu_file,
		tax_assign_file=tax_assign_file,
		otu_table_file=otu_table_file
		),
		removable_files=[otu_table_file]
	)
		
			