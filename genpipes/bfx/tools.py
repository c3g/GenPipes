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
import os
import re

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

## functions for awk tools ##

## functions for python tools ##
def py_addLengthRay(
        file_scaffolds_fasta, 
        length_file, 
        output,
        ini_section='add_length_ray'
        ):
    return Job(
        [file_scaffolds_fasta, length_file],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/addLengthRay.py \\
  -s {scaFile} \\
  -l {lenFile}""".format(
            scaFile=file_scaffolds_fasta,
            lenFile=length_file
        )
    )

def py_blastMatchSca(
        prefix_scaffolds_fasta, 
        blast_file, 
        output,
        ini_section='blast_match_sca'
        ):
    return Job(
        [prefix_scaffolds_fasta + ".fasta", blast_file],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/blastMatchSca.py \\
  -f {scaFile} \\
  -b {blastFile}""".format(
            scaFile=prefix_scaffolds_fasta,
            blastFile=blast_file
        )
    )

def py_equalFastqFile(
        fastq_ref, 
        fastq, 
        output,
        ini_section='equal_fastq_file'):
    return Job(
        [fastq_ref, fastq],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/equalFastqFile.py \\
  -r {ref} \\
  -f {fastq}""".format(
            ref=fastq_ref,
            fastq=fastq
        )
    )

def py_rrnaBAMcount(
        bam, 
        gtf, 
        output, 
        typ="transcript",
        ini_section='rrna_bam_count'):
    return Job(
        [bam],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
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
def py_parseTrinotateOutput(
        trinotate_annotation_report, 
        trinotate_report_genes_prefix, 
        trinotate_report_transcripts_prefix, 
        gene_id_column, 
        transcript_id_column, 
        isoforms_lengths_file, 
        job_name, 
        filters=None,
        ini_section='parse_trinotate_output'
        ):
    return Job(
        [trinotate_annotation_report, isoforms_lengths_file],
        [
            trinotate_report_genes_prefix + '_blast.tsv', trinotate_report_transcripts_prefix + '_blast.tsv' ,
            trinotate_report_genes_prefix + '_go.tsv', trinotate_report_transcripts_prefix + '_go.tsv',
            trinotate_report_transcripts_prefix + '_filtered.tsv'
        ],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        name=job_name,
        command="""\
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_genes_prefix} -i \"{gene_id_column}\" -l {isoforms_lengths} {other_options} &&
$PYTHON_TOOLS/parseTrinotateOutput.py -r {trinotate_annotation_report} -o {trinotate_report_transcripts_prefix} -i \"{transcript_id_column}\"{filters} {other_options}""".format(
            trinotate_annotation_report=trinotate_annotation_report,
            trinotate_report_genes_prefix=trinotate_report_genes_prefix,
            trinotate_report_transcripts_prefix=trinotate_report_transcripts_prefix,
            gene_id_column=gene_id_column,
            isoforms_lengths=isoforms_lengths_file,
            transcript_id_column=transcript_id_column,
            filters="" if not filters else " -f " + ' and '.join(filters),
            other_options=global_conf.global_get(ini_section, 'other_options')
        )
    )

def py_parseMergeCsv(
    input_files,
    delimiter,
    output,
    common,
    subset=None,
    exclude=None,
    left_join=None,
    sort_by=None,
    make_names=None,
    filters=None,
    ini_section='parse_merge_csv'
    ):
    return Job(
        input_files,
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
$PYTHON_TOOLS/parseMergeCsv.py \\
  -i {input_files} \\
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

def py_ampliconSeq(
    input_files,
    output_files,
    function,
    supplemental_parameters,
    ini_section='py_ampliconSeq'
    ):
    return Job(
        input_files,
        output_files,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/AmpliconSeq_script.py \\
  -m {function} \\
  {supplemental_parameters}""".format(
            function=function,
            supplemental_parameters=supplemental_parameters
        )
    )

def py_filterAssemblyToFastaToTsv(
        fasta_file, 
        filter_file, 
        fasta_id_column, 
        output,
        ini_section='py_filterAssemblyToFastaToTsv'
        ):
    return Job(
        [fasta_file , filter_file],
        [output + "." + ext for ext in ["fasta", "tsv"]],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
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

def dict2beds(
        dictionary,
        bed,
        ini_section='py_processIntervals'
        ):
        
    return Job(
        [dictionary],
        [bed],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python3 $PYTHON_TOOLS/dict2BEDs.py \\
  --dict {dictionary} \\
  --beds {bed}{chunk}{overlap}""".format(
            dictionary=dictionary if dictionary else global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath'),
            chunk=" --chunk " + global_conf.global_get(ini_section, 'chunk') if global_conf.global_get(ini_section, 'chunk') else "",
            overlap=" --overlap " + global_conf.global_get(ini_section, 'overlap') if global_conf.global_get(ini_section, 'overlap') else "",
            bed=bed
        )
    )

def preprocess_varscan(
        input,
        output=None,
        ini_section='preprocess_varscan'
):
    
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/preprocess.py \\
  --ref-depth RD --alt-depth AD \\
  {input} \\
  {output}""".format(
            input=input,
            output=output if output else "",
        )
    )

def fix_varscan_output(
        input,
        output=None,
        options=None,
        ini_section='fix_varscan_output'):
    
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python3 $PYTHON_TOOLS/fixVS2VCF.py {options} {input} \\
    {output}""".format(
            options=options if options else "",
            input=input if input else "",
            output=output if input else "",
        )
    )

def fix_genotypes_strelka(
        input, 
        output, 
        normal, 
        tumor,
        ini_section='DEFAULT'
        ):
        return Job(
            [input],
            [output],
            [
                ['DEFAULT', 'module_mugqic_tools'],
                ['DEFAULT', 'module_python']
            ],
            command="""\
	python3 $PYTHON_TOOLS/update_genotypes_strelka.py \\
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

def format2pcgr(input,
                output,
                filter,
                variant_type,
                tumor,
                ini_section='filter_ensemble'):
    
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python3 $PYTHON_TOOLS/format2pcgr.py \\
	-i {input} \\
	-o {output} \\
	-f {filter} \\
	-v {variant_type} \\
	-t {tumor}""".format(
            input=input if input else "",
            output=output if output else "",
            filter=filter,
            variant_type=variant_type,
            tumor=tumor,
        )
    )

def chunkBedbyFileNumber(
        input,
        output,
        chunk,
        ini_section='vardict_scatter_jobs'
):
    output_list = [os.path.join(output,
                                f"{idx:04d}-scattered.bed") for idx in range(1, chunk + 1)]
    return Job(
        [input],
        output_list,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python3 $PYTHON_TOOLS/chunkBedbyFileNumber.py \\
	--input {input} \\
	--output_dir {output} \\
	--chunk {chunk}""".format(
            input=input,
            output=output,
            chunk=chunk,
        )
    )


## functions for perl tools ##

def bed2interval_list(
    bed,
    output,
    ini_section='bed2interval_list'
    ):

    return Job(
        [bed],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_perl']
        ],
        command="""\
bed2IntervalList.pl \\
  --dict {dictionary} \\
  --bed {bed} \\
  > {output}""".format(
            dictionary=global_conf.global_get(ini_section, 'genome_dictionary', param_type='filepath'),
            bed=bed,
            output=output
        )
    )

def filter_long_indel(
        input, 
        output,
        ini_section='filter_long_indel'
        ):
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
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_perl']
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

def vcf2bed(
        input, 
        output,
        ini_section='vcf2bed'
        ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_perl']
        ],
        command="""\
cat {input} | perl $PERL_TOOLS/vcf2bed.pl - \\
  > {output}""".format(
            input=input,
            output=output
        )
    )

def rnaseqLight_kallisto(
        fastq_file1, 
        fastq_file2, 
        transcriptome_file, 
        tx2genes_file, 
        output_dir, 
        parameters, 
        job_name,
        ini_section='kallisto'
        ):
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
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R'],
            [ini_section, 'module_kallisto']
        ],
        name=job_name,
        command="""\
rm -rf {output_dir} && \\
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
def r_transcript_to_gene(
    abundance_tsv,
    output_tsv,
    tx2genes_file,
    ini_section='kallisto_count_matrix'
    ):
    return Job(
        [abundance_tsv],
        [output_tsv],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R'],
            [ini_section, 'module_mugqic_R_packages']
        ],
        command="""\
Rscript --vanilla $R_TOOLS/abundanceTranscript2geneLevel.R \\
  {abundance_tsv} \\
  {tx2genes}""".format(
            abundance_tsv=abundance_tsv,
            tx2genes=tx2genes_file
        )
    )


def r_create_kallisto_count_matrix(
        input_abundance_files, 
        output_dir, 
        data_type, 
        job_name,
        ini_section='kallisto_count_matrix'
        ):
    return Job(
        input_abundance_files,
        [os.path.join(output_dir, "all_readsets.abundance_" + data_type + ".csv")],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R'],
            [ini_section, 'module_mugqic_R_packages']
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

def r_create_kallisto_report(
        report_dir,
        inputs,
        output_file,
        ini_section='report'
        ):
    return Job(
        inputs,
        [output_file],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
        ],
        command="""\
R --no-save --args \\
    {report_dir} {output_file} \\
    < $R_TOOLS/kallistoReport.R""".format(
        report_dir=report_dir,
        output_file=output_file
    )
)

def r_select_scaffolds(
        inputs, 
        outputs, 
        folder_sca,
        kmer, 
        name_sample, 
        type_insert, 
        min_insert_size=200,
        ini_section='puure_select_scaffolds'
        ):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
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

def r_find_cluster(
        inputs, 
        outputs, 
        folder_sca, 
        kmer, 
        unmap_type, 
        name_sample, 
        type_insert, 
        max_insert_size=200, 
        min_mapping_quality=10,
        ini_section='puure_find_cluster'
        ):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
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

def r_find_insert(
        inputs, 
        outputs, 
        folder_sca, 
        kmer, 
        name_sample, 
        type_insert, 
        mean_coverage=20,
        max_insert_size=200, 
        min_overlap=2, 
        exclu_file="None",
        ini_section='puure_find_insert'):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
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

def r_filter_insert(
        inputs, 
        outputs, 
        folder_sca, 
        kmer, 
        name_sample,
        type_insert, 
        mean_coverage=20, 
        max_insert_size=200, 
        strand=1, 
        min_num_read=1, 
        mean_read_length=100,
        ini_section='puure_filter_insert'
        ):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
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

def methylkit_differential_analysis(
        design_file, 
        input_files, 
        outputfiles, 
        output_dir,
        ini_section='methylkit_differential_analysis'
        ):

    suffix = re.sub(r".*\.readset_", ".readset_", input_files[0])

    return Job(
        input_files,
        outputfiles,
        [
            [ini_section, 'module_R'],
            [ini_section, 'module_mugqic_tools']
        ],
        command="""\
Rscript $R_TOOLS/methylKit.R \\
  -design {design_file} \\
  -outdir {output_folder} \\
  -build {genome} \\
  -suff {input_suffix} \\
  {other_options}""".format(
            design_file=design_file,
            genome=global_conf.global_get(ini_section, 'assembly'),
            output_folder=output_dir,
            input_suffix=suffix,
            other_options=global_conf.global_get(ini_section, 'other_options')
        )
    )


## functions for bash tools ##

def sh_ihec_rna_metrics(
        input_bam, 
        input_name, 
        input_picard_dup, 
        output_dir,
        ini_section='IHEC_rnaseq_metrics'
        ):
    output_metrics=os.path.join(output_dir, input_name+".read_stats.txt")
    output_duplicates=os.path.join(output_dir, input_name+".duplicated.txt")

    return Job(
        [input_bam, input_picard_dup],
        [output_metrics, output_duplicates],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_samtools']
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
            intergenic_bed=global_conf.global_get(ini_section, 'intergenic_bed', param_type='filepath', required=True),
            rrna_bed=global_conf.global_get(ini_section, 'ribo_rna_bed', param_type='filepath', required=True),
            output_dir=output_dir
        )
    )

def sh_ihec_chip_metrics(
        chip_bam, 
        input_bam, 
        sample_name, 
        input_name, 
        chip_name, 
        chip_type, 
        chip_bed, 
        output_dir, 
        assembly, 
        crosscor_input,
        ini_section='IHEC_chipseq_metrics'
        ):
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
        [output_metrics, output_fingerprints, output_fingerprints_png, output_dedup_chip_bam, output_dedup_chip_bai, output_flagstats],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_samtools'],
            [ini_section, 'module_sambamba'],
            [ini_section, 'module_deeptools']
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
            threads=global_conf.global_get(ini_section, 'thread', param_type='int'),
            chip_bed=chip_bed,
            output_dir=output_dir,
            assembly=assembly
        ),
        removable_files=[output_fingerprints, output_fingerprints_png, output_dedup_chip_bam, output_dedup_chip_bam, output_dedup_chip_bai, output_dedup_input_bam, output_dedup_input_bai, output_flagstats]
    )

def sh_fastq_readname_edit(
    input_fastq,
    output_fastq,
    ini_section='fastq_readname_edit'
    ):
    return Job(
        [input_fastq],
        [output_fastq],
        [
            [ini_section, 'module_mugqic_tools'],
        ],
        command="""\
bash FastqReadNameEdit.sh \\
  -i {input_fastq} \\
  -o {output_fastq} \\
  -p {fastq_rel_path}""".format(
            input_fastq=input_fastq,
            output_fastq=output_fastq,
            fastq_rel_path=os.path.relpath(input_fastq, os.path.dirname(output_fastq))
        ),
        removable_files=[output_fastq]
    )

def sh_create_rmap(
        genome_digest_input, 
        rmap_output, 
        job_name,
        ini_section='create_rmap'
        ):
    return Job(
        [genome_digest_input],
        [rmap_output],
        [
            [ini_section, 'module_mugqic_tools']
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

def sh_create_baitmap(
        bait, 
        sorted_bait, 
        annotation, 
        output,
        ini_section='create_baitmap'
        ):
    return Job(
        [output + ".tmp", bait],
        [output],
        [
            [ini_section, 'module_mugqic_tools']
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

def sh_extract_bait_bed(
        ibed_file, 
        sample_name,
        ini_section='extract_bait_bed'
        ):
    return Job(
        [ibed_file],
        [ibed_file + ".bait"],
        [
            [ini_section, 'module_mugqic_tools']
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

def sh_extract_capture_bed(
        ibed_file, 
        sample_name,
        ini_section='extract_capture_bed'
        ):
    return Job(
        [ibed_file],
        [ibed_file + ".capture"],
        [
            [ini_section, 'module_mugqic_tools']
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


def sh_blastQC_ONT(
        output_directory, 
        reads_fastq_dir, 
        readset_name,
        ini_section='blastqc'
        ):
    return Job(
        [reads_fastq_dir],
        [os.path.join(output_directory, readset_name + ".blastHit_20MF_species.txt")],
        [
            [ini_section, "module_mugqic_tools"],
            [ini_section, "module_blast"],
            [ini_section, "module_python"]
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
        removable_files=[
            os.path.join(output_directory, "subsample_input.trim.blastres"),
            os.path.join(output_directory, "subsample_input.trim.fasta"),
            os.path.join(output_directory, "subsample_input.trim.fastq"),
            os.path.join(output_directory, "subsample_input.trim.qual"),
            os.path.join(output_directory, "subsample_input.fastq")
        ]
    )

def clean_otu(
        otu_table,
        ini_section='clean_otu'
        ):
    """
    Used by ampliconseq pipeline
    Cleans the OTU table : removes all the lines containing characters (e.g. division, OP3, WS6...)
    """
    bkp_otu_table = re.sub('OTU_data', 'OTU_data_BACKUP', otu_table)
    return Job(
        [otu_table],
        [bkp_otu_table],
        [
            [ini_section, 'module_mugqic_tools'],
        ],
        command="""\
cleanOTUtable.sh \\
 {otu}""".format(
            otu=otu_table,
         )
    )

## methylseq tools

def bismark_combine(
        input, 
        output,
        ini_section='bismark_combine'
        ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_perl']
        ],
        command="""\
methylProfile.bismark.pl \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
         )
    )

def cpg_stats(
        input, 
        cg_stats, 
        lambda_stats, 
        puc19_stats,
        ini_section='cpg_stats'
        ):
    return Job(
        [input],
        [cg_stats, lambda_stats, puc19_stats],
        [
            [ini_section, 'module_mugqic_tools']
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

def cpg_cov_stats(
        input, 
        output,
        ini_section='cpg_cov_stats'
        ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
python $PYTHON_TOOLS/CpG_coverageStats.py \\
 -i {input} \\
 -o {output}""".format(
            input=input,
            output=output
        )
    )

def filter_snp_cpg(
        input, 
        output,
        ini_section='filter_snp_cpg'
        ):
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bedops']
        ],
        command="""\
awk '$11>={cpg_threshold}' {input} | sort-bed - > {sorted_bed} && \\
zcat {filter_file} | vcf2bed --sort-tmpdir={tmp_dir} --max-mem={ram} - > {filter_sorted_bed} && \\
bedops --not-element-of \\
  {sorted_bed} \\
  {filter_sorted_bed} \\
  > {output}""". format(
            cpg_threshold=global_conf.global_get(ini_section, 'cpg_threshold', required=False) if global_conf.global_get(ini_section, 'cpg_threshold', required=False) else "5",
            input=input,
            sorted_bed=os.path.join(global_conf.global_get(ini_section, 'tmp_dir'), os.path.basename(input)+".tmp.sorted.bed"),
            filter_file=global_conf.global_get(ini_section, 'known_variants'),
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            ram=global_conf.global_get(ini_section, 'ram'),
            filter_sorted_bed=os.path.join(global_conf.global_get(ini_section, 'tmp_dir'), os.path.basename(global_conf.global_get(ini_section, 'known_variants'))+".tmp.sorted.bed"),
            output=output
        )
    )

def gembs_bcf_to_vcf(
    input,
    output,
    ini_section='gembs_bcf_to_vcf'
    ):

    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bcftools'],
            [ini_section, 'module_bedtools']
        ],
        command="""\
mkdir -p {tmp_dir}/tmp && \\
bcftools sort \\
  -T {tmp_dir}/tmp \\
  -m {ram} \\
  -Ov -o {filter_sorted_bed} \\
  {filter_file} && \\
bcftools view {bcftools_options} \\
  {input} \\
  -o {tmp_output} && \\
bcftools sort \\
  -T {tmp_dir}/tmp \\
  -m {ram} \\
  -Ov -o {sorted_tmp_output} \\
  {tmp_output} && \\
bedtools intersect {bedtools_options} \\
  -a {sorted_tmp_output} \\
  -b {filter_sorted_bed} \\
  > {output}""".format(
        filter_file=global_conf.global_get(ini_section, 'known_variants'),
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        ram=global_conf.global_get(ini_section, 'ram'),
        bcftools_options=global_conf.global_get(ini_section, 'bcftools_options'),
        input=input,
        tmp_output=os.path.join(global_conf.global_get(ini_section, 'tmp_dir'), os.path.basename(input) + ".tmp.vcf"),
        sorted_tmp_output=os.path.join(global_conf.global_get(ini_section, 'tmp_dir'), os.path.basename(input) + ".sorted.tmp.vcf"),
        bedtools_options=global_conf.global_get(ini_section, 'bedtools_options'),
        filter_sorted_bed=os.path.join(global_conf.global_get(ini_section, 'tmp_dir'), os.path.basename(global_conf.global_get(ini_section, 'known_variants'))+".tmp.sorted.vcf"),
        output=output
      )
    )

def prepare_methylkit(
        input, 
        output,
        ini_section='prepare_methylkit'
        ):
    cutoff=global_conf.global_get(ini_section, 'min_CpG', required=True)
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_mugqic_tools']
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

def gembs_format_cpg_report(
        input, 
        output,
        ini_section='methylseq_gembs_format_cpg_report'
        ):
    return Job(
            [input],
            [output],
            [
                [ini_section, 'module_mugqic_tools']
            ],
            command="""\
zcat {input} | tail -n +2 | \\
awk -F"\t" -v OFS="\t" '{{print $1, $3, $6, int($10*($11/100)+0.5), int($10*((100-$11)/100)+0.5), $12, $13}}' | \\
gzip -c > {output}""".format(
    input=input,
    output=output
    )
  )

def methylseq_metrics_report(
        sample_list, 
        inputs, 
        output, 
        target_bed,
        ini_section='methylseq_metrics_report'
        ):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        [output],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_samtools']
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

def methylseq_ihec_metrics_report(
        sample_name, 
        inputs, output, 
        output_all, 
        target_bed, 
        count,
        ini_section='ihec_sample_metrics'
        ):
    if not isinstance(inputs, list):
        inputs=[inputs]

    return Job(
        inputs,
        [output, output_all],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_samtools']
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

def r_somatic_signature(
        input, 
        output, 
        remove_file=None,
        ini_section='somatic_signature'):
    return Job(
        [input],
        [os.path.join(output,"Alexandrov_weigth.tsv")],
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R']
        ],
        command="""\
Rscript $R_TOOLS/somaticSignatureAlexandrov.R \\
  {input} \\
  {output} """.format(
        input=input,
        output=output
        ),removable_files=remove_file
    )

def design_report(
        input,
        output):
    return Job(
        [input],
        [output],
        command="""\
echo -e "# plot_type: 'table'
# section_name: 'Differential Analysis Design'
# description: 'The designs used in differential analysis are presented in the following table, which contains the sample names as well as the sample group membership per design. For each experimental design (column name), three conditions/groups are possible: 0, 1 and 2. If a sample is assigned 0, it is not included in a particular analysis. If a sample is assigned 1, the sample is considered as a member of the control group. If a sample is assigned 2, the sample is considered as a member of the test/case group.'" > {output}
cat {design_file} >> {output}""".format(
    design_file=input,
    output=output
)
    )

def diff_exp_report(
        input,
        output,
        section_name,
        contrast):
    return Job(
        [input],
        [output],
        command="""\
echo -e "# plot_type: 'table'
# section_name: {section_name}
# description: 'Differential Expression for contrast {contrast}'
# headers:
#   gene_symbol:
#       title: 'Gene Symbol'
#       description: 'Gene symbol if available, component ID otherwise'
#       placement: 970
#   log_FC:
#       title: 'log2FC'
#       description: 'log2 Fold Change of gene level expression'
#       placement: 980
#   log_CPM:
#       title: 'log2CPM'
#       description: 'log2 Counts Per Million of gene level expression'
#       placement: 990
#   deseq2.p.value:
#       title: 'DESeq p-value'
#       description: 'DESeq nominal p-value'
#       format: '{{:,s}}'
#       placement: 1000
#   deseq2.adj.pvalue:
#       title: 'DESeq adj p-value'
#       description: 'DESeq False Discovery Rate (FDR) adjusted p-value'
#       format: '{{:,s}}'
#       placement: 1010
#   edger.p.value:
#       title: 'edgeR p-value'
#       description: 'edgeR nominal p-value'
#       format: '{{:,s}}'
#       placement: 1020
#   edger.adj.p.value:
#       title: 'edgeR adj p-value'
#       description: 'edgeR False Discovery Rate (FDR) adjusted p-value'
#       format: '{{:,s}}'
#       placement: 1030" > {output}

cut -f1-8 {input} | head -n 11 >> {output}""".format(
        section_name=section_name,
        contrast=contrast,
        output=output,
        input=input
    )
)