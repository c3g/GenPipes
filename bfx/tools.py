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
        [
            trinotate_report_genes_prefix + '_blast.tsv', trinotate_report_transcripts_prefix + '_blast.tsv' ,
            trinotate_report_genes_prefix + '_go.tsv', trinotate_report_transcripts_prefix + '_go.tsv',
            trinotate_report_transcripts_prefix + '_filtered.tsv'
        ],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        name=job_name,
        command="""\
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_genes_prefix} -i \"{gene_id_column}\" -l {isoforms_lengths} &&
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_transcripts_prefix} -i \"{transcript_id_column}\"{filters}""".format(
            trinotate_annotation_report=trinotate_annotation_report,
            trinotate_report_genes_prefix=trinotate_report_genes_prefix,
            trinotate_report_transcripts_prefix=trinotate_report_transcripts_prefix,
            gene_id_column=gene_id_column,
            isoforms_lengths=isoforms_lengths_file,
            transcript_id_column=transcript_id_column,
            filters="" if not filters else " -f " + ' and '.join(filters)
        )
    )

def py_parseMergeCsv(input_files, delimiter, output , common, subset=None, exclude=None, left_join=None, sort_by=None, make_names=None, filters=None):
    return Job(
        input_files,
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
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

def py_ampliconSeq(input_files, output_files, function, supplemental_parameters):
    return Job(
        input_files,
        output_files,
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/AmpliconSeq_script.py \\
  -m {function} \\
  {supplemental_parameters}""".format(
            function=function,
            supplemental_parameters=supplemental_parameters
        )
    )

def py_filterAssemblyToFastaToTsv(fasta_file, filter_file, fasta_id_column, output):
    return Job(
        [fasta_file , filter_file],
        [output + "." + ext for ext in ["fasta", "tsv"]],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
$PYTHON_TOOLS/filterAssemblyToFastaToXls.py -f {fasta_file} \\
-o {output} \\
-l {filter_file} \\
-c {fasta_id_column} """.format(
            fasta_file=fasta_file,
            output=output,
            filter_file=filter_file,
            fasta_id_column=fasta_id_column
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

def fix_varscan_output(input, output, options=None):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/fixVS2VCF.py {options} {input} \\
    {output}""".format(
            options=options if options else "",
            input=input if input else "",
            output=output if input else "",
        )
    )


def cpg_cov_stats(input, output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/CpG_coverageStats.py \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
         )
    )

## functions for perl tools ##

def bed2interval_list(dictionary, bed, output):
    return Job(
        [dictionary, bed],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_perl']
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

def vcf2bed(input, output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_perl']
        ],
        command="""\
cat {input} | perl $PERL_TOOLS/vcf2bed.pl - \\
  > {output}""".format(
            input=input,
            output=output
        )
    )

def rnaseqLight_kallisto(fastq_file1, fastq_file2, transcriptome_file, tx2genes_file, output_dir, parameters, job_name):
    return Job(
        [
            fastq_file1,
            fastq_file2,
            transcriptome_file,
            tx2genes_file
        ],
        [
            os.path.join(output_dir, "abundance_transcripts.tsv"),
            os.path.join(output_dir, "abundance_genes.tsv")
        ],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R'],
            ['kallisto', 'module_kallisto']
        ],
        name=job_name,
        command="""\
bash rnaseq_light_kallisto.sh \\
  {fastq_file1} \\
  {fastq_file2} \\
  {transcriptome_file} \\
  {tx2genes_file} \\
  {output_dir} \\
  {parameters}""".format(
            fastq_file1=fastq_file1,
            fastq_file2=fastq_file2,
            transcriptome_file=transcriptome_file,
            tx2genes_file=tx2genes_file,
            output_dir=output_dir,
            parameters=parameters
        )
     )

## functions for R tools ##
def r_create_kallisto_count_matrix(input_abundance_files, output_dir, data_type, job_name):
    return Job(
        input_abundance_files,
        [os.path.join(output_dir, "all_readsets.abundance_" + data_type + ".csv")],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R'],
            ['DEFAULT', 'module_mugqic_R_packages']
        ],
        name=job_name,
        command="""\
R --no-save --args \\
  {input_abundance_files} \\
  {output_dir} \\
  {data_type} \\
  < $R_TOOLS/mergeKallistoCounts.R""".format(
            input_abundance_files=",".join(input_abundance_files),
            output_dir=output_dir,
            data_type=data_type #transcripts or genes
        )
    )

def r_create_kallisto_count_matrix(input_abundance_files, output_dir, data_type, job_name):
    return Job(
        input_abundance_files,
        [os.path.join(output_dir, "all_readsets.abundance_" + data_type + ".csv")],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_R'],
            ['DEFAULT', 'module_mugqic_R_packages']
        ],
        name=job_name,
        command="""\
R --no-save --args \\
  {input_abundance_files} \\
  {output_dir} \\
  {data_type} \\
  < $R_TOOLS/mergeKallistoCounts.R""".format(
            input_abundance_files=",".join(input_abundance_files),
            output_dir=output_dir,
            data_type=data_type #transcripts or genes
        )
    )

def r_select_scaffolds(inputs, outputs, folder_sca, kmer, name_sample, type_insert, min_insert_size=200):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
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

def r_find_cluster(inputs, outputs, folder_sca, kmer, unmap_type, name_sample, type_insert, max_insert_size=200, min_mapping_quality=10):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
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

def r_find_insert(inputs, outputs, folder_sca, kmer, name_sample, type_insert, mean_coverage=20, max_insert_size=200, min_overlap=2, exclu_file="None"):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
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

def r_filter_insert(inputs, outputs, folder_sca, kmer, name_sample, type_insert, mean_coverage=20, max_insert_size=200, strand=1, min_num_read=1, mean_read_length=100):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
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


## functions for bash tools ##

def sh_ihec_rna_metrics(input_bam, input_name, input_picard_dup, output_dir):
    output_metrics=os.path.join(output_dir, input_name+".read_stats.txt")
    output_duplicates=os.path.join(output_dir, input_name+".duplicated.txt")

    return Job(
        [input_bam, input_picard_dup],
        [output_metrics, output_duplicates],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_samtools']
        ],
        command="""\
IHEC_rnaseq_metrics.sh \\
    {input_bam} \\
    {input_name} \\
    {input_picard_dup} \\
    {intergenic_bed} \\
    {rrna_bed} \\
    {output_dir}""".format(
            input_bam=input_bam,
            input_name=input_name,
            input_picard_dup=input_picard_dup,
            intergenic_bed=config.param('IHEC_rnaseq_metrics', 'intergenic_bed', type='filepath',required=True),
            rrna_bed=config.param('IHEC_rnaseq_metrics', 'ribo_rna_bed', type='filepath',required=True),
            output_dir=output_dir
        )
    )

def sh_ihec_chip_metrics(chip_bam, input_bam, sample_name, input_name, chip_type, chip_bed, output_dir, assembly):
    output_metrics=os.path.join(output_dir, "IHEC_metrics_chipseq_"+ sample_name + ".txt")
    output_fingerprints=os.path.join(output_dir, sample_name+".fingerprint.txt")
    output_fingerprints_png=os.path.join(output_dir, sample_name+".fingerprint.png")
    output_dedup_chip_bam=os.path.join(output_dir, sample_name+".dedup.bam")
    output_dedup_chip_bai=os.path.join(output_dir, sample_name+".dedup.bam.bai")
    output_dedup_input_bam=os.path.join(output_dir, sample_name+"_IMPUT.dedup.bam")
    output_dedup_input_bai=os.path.join(output_dir, sample_name+"_IMPUT.dedup.bam.bai")
    output_flagstats=os.path.join(output_dir, sample_name+".markDup_flagstat.txt")
    crosscor_input =os.path.join(output_dir, sample_name + ".crosscor")
    return Job(
        [input_bam, chip_bam, chip_bed, crosscor_input],
        [output_metrics, output_dedup_chip_bam, output_dedup_chip_bai, output_flagstats],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_samtools'],
            ['DEFAULT', 'module_deeptools']
        ],
        command="""\
IHEC_chipseq_metrics_max.sh \\
    -d {chip_bam} \\
    -i {input_bam} \\
    -s {sample_name} \\
    -j {input_name} \\
    -t {chip_type} \\
    -n {threads} \\
    -p {chip_bed} \\
    -o {output_dir} \\
    -a {assembly}""".format(
            input_bam=input_bam,
            input_name=input_name,
            sample_name=sample_name,
            chip_bam=chip_bam,
            chip_type=chip_type,
            threads=config.param('IHEC_chipseq_metrics', 'thread', type='int') if config.param('IHEC_chipseq_metrics', 'thread', type='int',required=False) else 1,
            chip_bed=chip_bed,
            output_dir=output_dir,
            assembly=assembly
        ),
        removable_files=[output_fingerprints,output_fingerprints_png,output_dedup_chip_bam,output_dedup_chip_bam,output_dedup_chip_bai,output_dedup_input_bam,output_dedup_input_bai,output_flagstats]
    )

def sh_fastq_readname_edit(fastq, job_name):
    return Job(
        [fastq],
        [fastq + ".edited.gz"],
        [
            ['DEFAULT', 'module_mugqic_tools'],
        ],
        command="""\
bash FastqReadNameEdit.sh \\
  -i {input_fastq} \\
  -o {output_fastq} \\
  -p {fastq_abs_path}""".format(
            input_fastq=fastq,
            output_fastq=fastq + ".edited.gz",
            fastq_abs_path=os.path.abspath(fastq)
        ),
        name=job_name,
        removable_files = [fastq + ".edited.gz"]
    )

def sh_create_rmap(genome_digest_input, rmap_output, job_name):
    return Job(
        [genome_digest_input],
        [rmap_output],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""
bash createRmapFile.sh \\
  {infile} \\
  {outfile}""".format(
            infile=genome_digest_input,
            outfile=rmap_output
        ),
        name=job_name
    )


def sh_create_baitmap(bait, sorted_bait, annotation, output):
    return Job(
        [output + ".tmp", bait],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""
bash createBaitMapFile.sh \\
  {input_file}
  {bait_file} \\
  {sorted_bait_file} \\
  {annotation} \\
  {tmp_file} \\
  {output_file}""".format(
            input_file=output + ".tmp",
            bait_file=bait,
            sorted_bait_file=sorted_bait,
            annotation=annotation,
            tmp_file=output + ".tmp2",
            output_file=output
        ),
        removable_files=[input],
        name="create_baitmap_file.addAnno." + os.path.basename(output)
    )

def sh_extract_bait_bed(ibed_file, sample_name):
    return Job(
        [ibed_file],
        [ibed_file + ".bait"],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command = """
bash extractBaitBed.sh \\
  {input_file} \\
  {output_tmp} \\
  {output_file}""".format(
            input_file=ibed_file,
            output_tmp=ibed_file + ".bait.tmp",
            output_file=ibed_file + ".bait"
        ),
        name="extract_bait_bed." + sample_name,
        removable_files=[ibed_file + ".bait"]
    )

def sh_extract_capture_bed(ibed_file, sample_name):
    return Job(
        [ibed_file],
        [ibed_file + ".capture"],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command = """
bash extractCaptureBed.sh \\
  {input_file} \\
  {output_tmp} \\
  {output_file}""".format(
            input_file=ibed_file,
            output_tmp=ibed_file + ".capture.tmp",
            output_file=ibed_file + ".capture"
        ),
        name="extract_capture_bed." + sample_name,
        removable_files=[ibed_file + ".capture"]
    )

def clean_otu(otu_table):
    """
    Used by ampliconseq pipeline
    Cleans the OTU table : removes all the lines containing characters (e.g. division, OP3, WS6...)
    """
    bkp_otu_table = re.sub('OTU_data', 'OTU_data_BACKUP', otu_table)
    return Job(
        [otu_table],
        [bkp_otu_table],
        [
            ['DEFAULT', 'module_mugqic_tools'],
        ],
        command="""\
cleanOTUtable.sh \\
 {otu}""".format(
            otu=otu_table,
         )
    )

## methylseq tools

def bismark_combine(input, output):
    return Job(
        [input],
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_perl']
        ],
        command="""\
methylProfile.bismark.pl \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
         )
    )

def cpg_stats(input, cg_stats, lambda_stats, puc19_stats):
    return Job(
        [input],
        [cg_stats, lambda_stats, puc19_stats],
        [
            ['DEFAULT', 'module_mugqic_tools']
        ],
        command="""\
bash cpgStats.sh \\
  {input} \\
  {cg_stats} \\
  {lambda_stats} \\
  {puc19_stats}""".format(
            input=input,
            cg_stats=cg_stats,
            lambda_stats=lambda_stats,
            puc19_stats=puc19_stats
        )
    )

def methylseq_metrics_report(sample_list, inputs, output, target_bed):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        [output],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_samtools']
        ],
        command="""\
bash methylseq_metrics.sh \\
  {sample_list} \\
  {output_file} \\
  {targeted_flag}""".format(
            sample_list=",".join(sample_list),
            output_file=output,
            targeted_flag=1 if target_bed else 0
        )
    )

def methylseq_ihec_metrics_report(sample_name, inputs, output, output_all, target_bed, count):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        [output, output_all],
        [
            ['DEFAULT', 'module_mugqic_tools'],
            ['DEFAULT', 'module_samtools']
        ],
        command="""\
bash IHEC_methylseq_metrics.sh \\
  {sample_name} \\
  {output_file} \\
  {output_all_file} \\
  {targeted_flag} \\
  {counter}""".format(
            sample_name=sample_name,
            output_file=output,
            output_all_file=output_all,
            counter=count,
            targeted_flag=1 if target_bed else 0
        )
    )
