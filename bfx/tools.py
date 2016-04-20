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
import os

# MUGQIC Modules
from core.config import *
from core.job import *

## functions for awk tools ##

## functions for python tools ## 
def py_addLengthRay (file_scaffolds_fasta, length_file, output):
    return Job(
        [file_scaffolds_fasta, length_file],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/addLengthRay.py \\
  -s {scaFile} \\
  -l {lenFile}""".format(
        scaFile=file_scaffolds_fasta,
        lenFile=length_file
        )
    )

def py_blastMatchSca (prefix_scaffolds_fasta, blast_file, output):
    return Job(
        [prefix_scaffolds_fasta + ".fasta", blast_file],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/blastMatchSca.py \\
  -f {scaFile} \\
  -b {blastFile}""".format(
        scaFile=prefix_scaffolds_fasta,
        blastFile=blast_file
        )
    )

def py_equalFastqFile (fastq_ref, fastq, output):
    return Job(
        [fastq_ref, fastq],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/equalFastqFile.py \\
  -r {ref} \\
  -f {fastq}""".format(
        ref=fastq_ref,
        fastq=fastq
        )
    )

def py_rrnaBAMcount (bam, gtf, output, typ="transcript"):
    return Job(
        [bam],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/rrnaBAMcounter.py \\
  -i {bam} \\
  -g {gtf} \\
  -o {output} \\
  -t {typ}""".format(
        bam=bam,
        gtf=gtf,
        output=output,
        typ=typ
        )
    )
# Parse Trinotate output, create best blast annotated file, GO terms and a list of filtered configs
def py_parseTrinotateOutput(trinotate_annotation_report, trinotate_report_genes_prefix, trinotate_report_transcripts_prefix, gene_id_column, transcript_id_column, isoforms_lengths_file, job_name, filters=None):  
    return Job(
        [trinotate_annotation_report, isoforms_lengths_file], 
        [trinotate_report_genes_prefix + '_blast.tsv', trinotate_report_transcripts_prefix + '_blast.tsv' , 
        trinotate_report_genes_prefix + '_go.tsv', trinotate_report_transcripts_prefix + '_go.tsv',
        trinotate_report_transcripts_prefix + '_filtered.tsv'],
        [['DEFAULT', 'module_mugqic_tools'],
         ['DEFAULT', 'module_python']
         ],
        name=job_name,
        command="""\
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_genes_prefix} -i \"{gene_id_column}\" -l {isoforms_lengths} &&
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_transcripts_prefix} -i \"{transcript_id_column}\"{filters}""".format(
        trinotate_annotation_report=trinotate_annotation_report,
        trinotate_report_genes_prefix=trinotate_report_genes_prefix,
        trinotate_report_transcripts_prefix=trinotate_report_transcripts_prefix,
        gene_id_column = gene_id_column,
        isoforms_lengths=isoforms_lengths_file, 
        transcript_id_column=transcript_id_column,
        filters="" if not filters else " -f " + ' and '.join(filters)
        )
    )

def py_parseMergeCsv(input_files, delimiter, output , common, subset=None, exclude=None, left_join=None, sort_by=None, make_names=None, filters=None):  
    return Job(
        input_files, 
        [output],
        [['DEFAULT', 'module_mugqic_tools'],
         ['DEFAULT', 'module_python']
         ], 
        command="""\
$PYTHON_TOOLS/parseMergeCsv.py -i {input_files} \\
      -o {output} \\
      -c {common_columns} \\
      -d {delimiter} {subset}{toexclude}{left_outer_join}{sort_by_field}{make_names}{filters}""".format(
        input_files=" ".join(input_files),
        output=output,
        common_columns=common,        
        delimiter=delimiter, 
        subset=" -s " + subset if subset else "", 
        toexclude=" -x " + exclude if exclude else "", 
        left_outer_join=" -l " if  left_join else "",
        sort_by_field=" -t " + sort_by if sort_by else "",
        make_names=" -n " if make_names else "",
        filters="" if not filters else " -f " + ' and '.join(filters)
        )
    )

def py_filterAssemblyToFastaToTsv(fasta_file, filter_file, fasta_id_column, output):  
    return Job(
        [ fasta_file , filter_file], 
        [ output + "." + ext for ext in ["fasta", "tsv"] ],
        [['DEFAULT', 'module_mugqic_tools'],
         ['DEFAULT', 'module_python']
         ], 
        command="""\
$PYTHON_TOOLS/filterAssemblyToFastaToXls.py -f {fasta_file} \\
-o {output} \\
-l {filter_file} \\
-c {fasta_id_column} """.format(
        fasta_file = fasta_file,
        output=output,
        filter_file=filter_file,        
        fasta_id_column=fasta_id_column
        )
    )


## functions for perl tools ##
def bed2interval_list(dictionary, bed, output):
    return Job(
        [dictionary, bed],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT' , 'module_perl']
        ],
        command="""\
bed2IntervalList.pl \\
  --dict {dictionary} \\
  --bed {bed} \\
  > {output}""".format(
        dictionary=dictionary if dictionary else config.param('DEFAULT', 'genome_dictionary', type='filepath'),
        bed=bed,
        output=output
        )
    )

def dict2beds(dictionary,beds):
    return Job(
        [dictionary],
        beds,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
dict2BEDs.py \\
  --dict {dictionary} \\
  --beds {beds}""".format(
        dictionary=dictionary if dictionary else config.param('DEFAULT', 'genome_dictionary', type='filepath'),
        beds=' '.join(beds)
        )
    )

def preprocess_varscan(input,output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/preprocess.py \\
  --ref-depth RD --alt-depth AD \\
  {input} \\
  | bgzip -cf > {output}""".format(
        input=input,
        output=output
        )
    )

def filter_long_indel(input, output):
    pre_gzip_command = ""
    post_rm_command = ""
    input_filename, input_file_extension = os.path.splitext(input)
    if input_file_extension == ".bgz" :
        pre_gzip_command="""\
zcat {input} > {input_filename} && """.format(
            input=input,
            input_filename=input_filename
        )
        post_rm_command=""" && \
rm {input_filename} """.format(
            input_filename=input_filename
        )
        input_next=input_filename
    else :
        input_next=input
        
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_perl']
        ],
        command="""\
{pre_gzip_command}filterLongIndel.pl \\
  {input} \\
  > {output}{post_rm_command}""".format(
        pre_gzip_command=pre_gzip_command,
        input=input_next,
        output=output,
        post_rm_command=post_rm_command
        ),
        removable_files=[output]
    )


## functions for R tools ##

def r_select_scaffolds(input, output, folder_sca, kmer, name_sample, type_insert, min_insert_size=200):
    return Job(
        input,
        output,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R']
        ],
        command="""\
R --no-save --args \\
  {folder_sca} \\
  {kmer} \\
  {name_sample} \\
  {type_insert} \\
  {min_insert_size} \\
  < $R_TOOLS/puureAnalyseSelectSca.r""".format(
        folder_sca=folder_sca,
        kmer=kmer,
        name_sample=name_sample,
        type_insert=type_insert,
        min_insert_size=min_insert_size
        )
    )

def r_find_cluster(input, output, folder_sca, kmer, unmap_type, name_sample, type_insert, max_insert_size=200, min_mapping_quality=10):
    return Job(
        input,
        output,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R']
        ],
        command="""\
R --no-save --args \\
  {folder_sca} \\
  {kmer} \\
  {name_sample} \\
  {type_insert} \\
  {min_mapping_quality} \\
  {max_insert_size} \\
  < $R_TOOLS/puureAnalyseFindCluster{unmap_type}.r""".format(
        folder_sca=folder_sca,
        kmer=kmer,
        name_sample=name_sample,
        type_insert=type_insert,
        min_mapping_quality=min_mapping_quality,
        max_insert_size=max_insert_size,
        unmap_type=unmap_type
        )
    )

def r_find_insert(input, output, folder_sca, kmer, name_sample, type_insert, mean_coverage=20, max_insert_size=200, min_overlap=2, exclu_file="None"):
    return Job(
        input,
        output,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R']
        ],
        command="""\
R --no-save --args \\
  {folder_sca} \\
  {kmer} \\
  {name_sample} \\
  {type_insert} \\
  {mean_coverage} \\
  {max_insert_size} \\
  {min_overlap} \\
  {exclu_file} \\
  < $R_TOOLS/puureAnalyseFindInsert.r""".format(
        folder_sca=folder_sca,
        kmer=kmer,
        name_sample=name_sample,
        type_insert=type_insert,
        mean_coverage=mean_coverage,
        max_insert_size=max_insert_size,
        min_overlap=min_overlap,
        exclu_file=exclu_file
        )
    )

def r_filter_insert(input, output, folder_sca, kmer, name_sample, type_insert, mean_coverage=20, max_insert_size=200, strand=1, min_num_read=1, mean_read_length=100):
    return Job(
        input,
        output,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R']
        ],
        command="""\
R --no-save --args \\
  {folder_sca} \\
  {kmer} \\
  {name_sample} \\
  {type_insert} \\
  {mean_coverage} \\
  {max_insert_size} \\
  {strand} \\
  {min_num_read} \\
  {mean_read_length} \\
  < $R_TOOLS/puureAnalyseFindInsert.r""".format(
        folder_sca=folder_sca,
        kmer=kmer,
        name_sample=name_sample,
        type_insert=type_insert,
        mean_coverage=mean_coverage,
        max_insert_size=max_insert_size,
        strand=strand,
        min_num_read=min_num_read,
        mean_read_length=mean_read_length
        )
    )
