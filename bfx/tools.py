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
import os
import re

# MUGQIC Modules
from core.config import config
from core.job import Job

## functions for awk tools ##

## functions for python tools ##
def py_addLengthRay (file_scaffolds_fasta, length_file, output):
    return Job(
        [file_scaffolds_fasta, length_file],
        [output],
        [
            ['add_length_ray', 'module_mugqic_tools'],
            ['add_length_ray', 'module_python']
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
            ['blast_match_sca', 'module_mugqic_tools'],
            ['blast_match_sca', 'module_python']
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
            ['equal_fastq_file', 'module_mugqic_tools'],
            ['equal_fastq_file', 'module_python']
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
            ['rrna_bam_count', 'module_mugqic_tools'],
            ['rrna_bam_count', 'module_python']
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
            ['parse_trinotate_output', 'module_mugqic_tools'],
            ['parse_trinotate_output', 'module_python']
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
            ['parse_merg_csv', 'module_mugqic_tools'],
            ['parse_merg_csv', 'module_python']
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
            ['py_ampliconSeq', 'module_mugqic_tools'],
            ['py_ampliconSeq', 'module_python']
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
            ['py_filterAssemblyToFastaToTsv', 'module_mugqic_tools'],
            ['py_filterAssemblyToFastaToTsv', 'module_python']
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
            ['dict2beds', 'module_mugqic_tools'],
            ['dict2beds', 'module_python']
        ],
        command="""\
dict2BEDs.py \\
  --dict {dictionary} \\
  --beds {beds}""".format(
            dictionary=dictionary if dictionary else config.param('dict2beds', 'genome_dictionary', param_type='filepath'),
            beds=' '.join(beds)
        )
    )

def preprocess_varscan(input,output):
    return Job(
        [input],
        [output],
        [
            ['preprocess_varscan', 'module_mugqic_tools'],
            ['preprocess_varscan', 'module_python']
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

def fix_varscan_output(input, output=None, options=None):
    return Job(
        [input],
        [output],
        [
            ['fix_varscan_output', 'module_mugqic_tools'],
            ['fix_varscan_output', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/fixVS2VCF.py {options} {input} \\
    {output}""".format(
            options=options if options else "",
            input=input if input else "",
            output=output if input else "",
        )
    )

def fix_genotypes_strelka(input, output, normal, tumor):
        return Job(
            [input],
            [output],
            [
                ['DEFAULT', 'module_mugqic_tools'],
                ['DEFAULT', 'module_python']
            ],
            command="""\
	python $PYTHON_TOOLS/update_genotypes_strelka.py \\
	    -i {input} \\
	    -o {output} \\
	    -n {normal} \\
	    -t {tumor}""".format(
                input=input if input else "",
                output=output if input else "",
                normal=normal,
                tumor=tumor,
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

def bed2interval_list(
    bed,
    output
    ):
    
    return Job(
        [bed],
        [output],
        [
            ['bed2interval_list', 'module_mugqic_tools'],
            ['bed2interval_list', 'module_perl']
        ],
        command="""\
bed2IntervalList.pl \\
  --dict {dictionary} \\
  --bed {bed} \\
  > {output}""".format(
            dictionary=config.param('bed2interval_list', 'genome_dictionary', param_type='filepath'),
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
            ['filter_long_indel', 'module_mugqic_tools'],
            ['filter_long_indel', 'module_perl']
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
            ['vcf2bed', 'module_mugqic_tools'],
            ['vcf2bed', 'module_perl']
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
            os.path.join(output_dir, "abundance_genes.tsv"),
            os.path.join(output_dir, "kallisto_quant.log")
        ],
        [
            ['kallisto', 'module_mugqic_tools'],
            ['kallisto', 'module_R'],
            ['kallisto', 'module_kallisto']
        ],
        name=job_name,
        command="""\
bash rnaseq_light_kallisto.sh \\
  {output_dir} \\
  {parameters} \\
  {transcriptome_file} \\
  {tx2genes_file} \\
  {fastq_file1} \\
  {fastq_file2}""".format(
            output_dir=output_dir,
            parameters=parameters,
            transcriptome_file=transcriptome_file,
            tx2genes_file=tx2genes_file,
            fastq_file1=fastq_file1,
            fastq_file2=fastq_file2
        ),
        multiqc_files=[os.path.join(output_dir, "kallisto_quant.log")]
     )

## functions for R tools ##
def r_create_kallisto_count_matrix(input_abundance_files, output_dir, data_type, job_name):
    return Job(
        input_abundance_files,
        [os.path.join(output_dir, "all_readsets.abundance_" + data_type + ".csv")],
        [
            ['kallisto_count_matrix', 'module_mugqic_tools'],
            ['kallisto_count_matrix', 'module_R'],
            ['kallisto_count_matrix', 'module_mugqic_R_packages']
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
            ['puure_select_scaffolds', 'module_mugqic_tools'],
            ['puure_select_scaffolds', 'module_R']
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
            ['puure_find_cluster', 'module_mugqic_tools'],
            ['puure_find_cluster', 'module_R']
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
            ['puure_find_insert', 'module_mugqic_tools'],
            ['puure_find_insert', 'module_R']
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
            ['puure_filter_insert', 'module_mugqic_tools'],
            ['puure_filter_insert', 'module_R']
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

def methylkit_differential_analysis(design_file, input_files, outputfiles, output_dir):

    suffix = re.sub(".*\.readset_", ".readset_", input_files[0])

    return Job(
        input_files,
        outputfiles,
        [
            ['methylkit_differential_analysis', 'module_R'],
            ['methylkit_differential_analysis', 'module_mugqic_tools']
        ],
        command="""\
R --no-save '--args \\
  -design {design_file} \\
  -outdir {output_folder} \\
  -build {genome} \\
  -suff {input_suffix} \\
  {other_options}' \\
  < $R_TOOLS/methylKit.R""".format(
            design_file=design_file,
            genome=config.param('methylkit_differential_analysis', 'assembly'),
            output_folder=output_dir,
            input_suffix=suffix,
            other_options=config.param('methylkit_differential_analysis', 'other_options')
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
            ['IHEC_rnaseq_metrics', 'module_mugqic_tools'],
            ['IHEC_rnaseq_metrics', 'module_samtools']
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
            intergenic_bed=config.param('IHEC_rnaseq_metrics', 'intergenic_bed', param_type='filepath', required=True),
            rrna_bed=config.param('IHEC_rnaseq_metrics', 'ribo_rna_bed', param_type='filepath', required=True),
            output_dir=output_dir
        )
    )

def sh_ihec_chip_metrics(chip_bam, input_bam, sample_name, input_name, chip_name, chip_type, chip_bed, output_dir, assembly, crosscor_input):
    output_dedup_input_bam = None
    output_dedup_input_bai = None
    output_fingerprints = None
    output_fingerprints_png = None
    # crosscor_input = None
    if input_bam:
        output_dedup_input_bam = os.path.join(output_dir, input_name, sample_name + "." + "Input.dedup.bam")
        output_dedup_input_bai = os.path.join(output_dir, input_name, sample_name + "." + "Input.dedup.bam.bai")
        output_fingerprints = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".fingerprint.tsv")
        output_fingerprints_png = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".fingerprint.png")
        # crosscor_input = os.path.join(output_dir, sample_name + "." + chip_name + ".crosscor")
    output_metrics = os.path.join(output_dir, chip_name, "IHEC_chipseq_metrics." + sample_name + "." + chip_name + ".tsv")
    # output_fingerprints = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".fingerprint.tsv")
    # output_fingerprints_png = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".fingerprint.png")
    output_dedup_chip_bam = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".dedup.bam")
    output_dedup_chip_bai = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".dedup.bam.bai")
    output_flagstats = os.path.join(output_dir, chip_name, sample_name + "." + chip_name + ".markDup_flagstat.txt")
    # crosscor = os.path.join(output_dir, sample_name + ".crosscor")

    return Job(
        [chip_bam, input_bam, chip_bed, crosscor_input],
        [output_metrics, output_fingerprints, output_fingerprints_png, output_dedup_chip_bam, output_dedup_chip_bai, output_flagstats, output_dedup_input_bam, output_dedup_input_bai],
        [
            ['IHEC_chipseq_metrics', 'module_mugqic_tools'],
            ['IHEC_chipseq_metrics', 'module_samtools'],
            ['IHEC_chipseq_metrics', 'module_sambamba'],
            ['IHEC_chipseq_metrics', 'module_deeptools']
        ],
        command="""\
IHEC_chipseq_metrics_max.sh \\
    -d {chip_bam} \\
    -i {input_bam} \\
    -s {sample_name} \\
    -j {input_name} \\
    -t {chip_type} \\
    -c {chip_name} \\
    -n {threads} \\
    -p {chip_bed} \\
    -o {output_dir} \\
    -a {assembly}""".format(
            input_bam=input_bam,
            input_name=input_name,
            sample_name=sample_name,
            chip_bam=chip_bam,
            chip_type=chip_type,
            chip_name=chip_name,
            threads=config.param('IHEC_chipseq_metrics', 'thread', param_type='int'),
            chip_bed=chip_bed,
            output_dir=output_dir,
            assembly=assembly
        ),
        removable_files=[output_fingerprints, output_fingerprints_png, output_dedup_chip_bam, output_dedup_chip_bam, output_dedup_chip_bai, output_dedup_input_bam, output_dedup_input_bai, output_flagstats]
    )

def sh_fastq_readname_edit(
    fastq,
    working_dir,
    job_name
    ):
    return Job(
        [fastq],
        [fastq + ".edited.gz"],
        [
            ['fastq_readname_edit', 'module_mugqic_tools'],
        ],
        command="""\
bash FastqReadNameEdit.sh \\
  -i {input_fastq} \\
  -o {output_fastq} \\
  -p {fastq_abs_path}""".format(
            input_fastq=fastq if os.path.isabs(fastq) else os.path.join(working_dir, fastq),
            output_fastq=fastq + ".edited.gz" if os.path.isabs(fastq) else os.path.join(working_dir, fastq + ".edited.gz"),
            fastq_abs_path=fastq if os.path.isabs(fastq) else os.path.join(working_dir, fastq)
        ),
        name=job_name,
        removable_files = [fastq + ".edited.gz"]
    )

def sh_create_rmap(genome_digest_input, rmap_output, job_name):
    return Job(
        [genome_digest_input],
        [rmap_output],
        [
            ['create_rmap', 'module_mugqic_tools']
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
            ['create_baitmap', 'module_mugqic_tools']
        ],
        command="""
bash createBaitMapFile.sh \\
  {input_file} \\
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
            ['extract_bait_bed', 'module_mugqic_tools']
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
            ['extract_capture_bed', 'module_mugqic_tools']
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


def sh_blastQC_ONT(output_directory, reads_fastq_dir, readset_name):
    return Job(
        [reads_fastq_dir],
        [os.path.join(output_directory, readset_name + "blastHit_20MF_species.txt")],
        [
            ["blastqc", "module_mugqic_tools"],
            ["blastqc", "module_blast"],
            ["blastqc", "module_python"]
        ],
        command="""
bash runBlastQC_ONT.sh \\
  {output_directory} \\
  {reads_fastq_dir} \\
  {readset_name}""".format(
            output_directory=output_directory,
            reads_fastq_dir=reads_fastq_dir,
            readset_name=readset_name
        ),
        name="blastQC." + readset_name,
        removable_files=[os.path.join(output_directory, "subsample_input.trim.blastres"),
                         os.path.join(output_directory, "subsample_input.trim.fasta"),
                         os.path.join(output_directory, "subsample_input.trim.fastq"),
                         os.path.join(output_directory, "subsample_input.trim.qual"),
                         os.path.join(output_directory, "subsample_input.fastq")]
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
            ['clean_otu', 'module_mugqic_tools'],
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
            ['bismark_combine', 'module_mugqic_tools'],
            ['bismark_combine', 'module_perl']
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
            ['cpg_stats', 'module_mugqic_tools']
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

def cpg_cov_stats(input, output):
    return Job(
        [input],
        [output],
        [
            ['cpg_cov_stats', 'module_mugqic_tools'],
            ['cpg_cov_stats', 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/CpG_coverageStats.py \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
        )
    )

def filter_snp_cpg(input, output):
    return Job(
        [input],
        [output],
        [
            ['filter_snp_cpg', 'module_bedops']
        ],
        command="""\
awk '$11>=5' {input} | sort-bed - > {sorted_bed} && \\
zcat {filter_file} | vcf2bed - > {filter_sorted_bed} && \\
bedops --not-element-of \\
  {sorted_bed} \\
  {filter_sorted_bed} \\
  > {output}""". format(
            input=input,
            sorted_bed=os.path.join(config.param('filter_snp_cpg', 'tmp_dir'), os.path.basename(input)+".tmp.sorted.bed"),
            filter_file=config.param('filter_snp_cpg', 'known_variants'),
            filter_sorted_bed=os.path.join(config.param('filter_snp_cpg', 'tmp_dir'), os.path.basename(config.param('filter_snp_cpg', 'known_variants'))+".tmp.sorted.bed"),
            output=output
        )
    )

def prepare_methylkit(input, output):
    cutoff=config.param('prepare_methylkit', 'min_CpG', required=True)
    return Job(
        [input],
        [output],
        [
            ['prepare_methylkit', 'module_mugqic_tools']
        ],
        command="""\
echo -e \"chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\" \\
  > {output} && \\
cat {input} | \\
awk -F"\t" '$11>={cutoff} {{print $1"."$2"\t"$1"\t"$2"\tF\t"$11"\t"$12"\t"(100-$12)}}' \\
  >> {output}""".format(
            input=input,
            output=output,
            cutoff=int(cutoff)
        )
    )

def methylseq_metrics_report(sample_list, inputs, output, target_bed):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        [output],
        [
            ['methylseq_metrics_report', 'module_mugqic_tools'],
            ['methylseq_metrics_report', 'module_samtools']
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
            ['ihec_sample_metrics', 'module_mugqic_tools'],
            ['ihec_sample_metrics', 'module_samtools']
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

def r_somatic_signature(input, output, remove_file=None):
    return Job(
        [input],
        [os.path.join(output,"Alexandrov_weigth.tsv")],
        [
            ['somatic_signature', 'module_mugqic_tools'],
            ['somatic_signature', 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/somaticSignatureAlexandrov.R \\
  {input} \\
  {output} """.format(
        input=input,
        output=output
        ),removable_files=remove_file
    )
