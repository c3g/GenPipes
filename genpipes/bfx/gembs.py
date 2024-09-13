################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
from ..core.config import global_conf
from ..core.job import Job

def make_config(genpipes_dir, output):

    config_content = f"""\
# gemBS config file 
# base = {genpipes_dir}

reference = {global_conf.global_get('gembs_prepare', 'genome_fasta')}

index_dir = {global_conf.global_get('gembs_prepare', 'gem3_index_dir')}
sequence_dir = {genpipes_dir}/trim/@BARCODE
bam_dir = {genpipes_dir}/alignment/@BARCODE
bcf_dir = {genpipes_dir}/methylation_call/@BARCODE
extract_dir = {genpipes_dir}/variants/@BARCODE
report_dir = {genpipes_dir}/report

project = {global_conf.global_get('gembs_report', 'project_name')}
species = {global_conf.global_get('default', 'scientific_name')}_{global_conf.global_get('default', 'assembly')}

# Default parameters
threads = 1
cores = 1
jobs = 1

[dbsnp]
dbsnp_files = {global_conf.global_get('default', 'known_variants')}
dbsnp_index = {global_conf.global_get('default', 'gem3_index_dir')}/dbSNP_gemBS.idx
{"dbsnp_selected = " + global_conf.global_get('default', 'selected_snps') if global_conf.global_get('default', 'selected_snps', required = False) else ""}
{"dbsnp_chrom_alias = " + global_conf.global_get('default', 'alias_file') if global_conf.global_get('default', 'alias_file', required = False) else ""}

[mapping]
memory = {global_conf.global_get('gembs_map', 'ram')}
threads = {global_conf.global_get('gembs_map', 'threads')}
cores = {global_conf.global_get('gembs_map', 'cores')}
merge_cores = {global_conf.global_get('gembs_map', 'merge_cores')}
merge_memory = {global_conf.global_get('gembs_map', 'merge_ram')}

# set names of spiked in conversion controls
{"underconversion_sequence = " + global_conf.global_get('gembs_map', 'underconversion_sequence') if global_conf.global_get('gembs_map', 'underconversion_sequence') else ""}
{"overconversion_sequence = " + global_conf.global_get('gembs_map', 'overconversion_sequence') if global_conf.global_get('gembs_map', 'overconversion_sequence') else ""}

# Include a standard configuration file with parameters
# defined for the standard IHEC WGBS pipeline
{"include " + global_conf.global_get('gembs_map', 'standard_IHEC') if global_conf.global_get('gembs_map', 'standard_IHEC', required = False) else ""}
 
[calling]
contig_pool_limit = {global_conf.global_get('gembs_call', 'contig_pool_limit')}
cores = {global_conf.global_get('gembs_call', 'cores')}
call_threads = {global_conf.global_get('gembs_call', 'threads')}
memory = {global_conf.global_get('gembs_call', 'ram')}
left_trim = {global_conf.global_get('gembs_call', 'left_trim')}
right_trim = {global_conf.global_get('gembs_call', 'right_trim')}
contig_list = {global_conf.global_get('gembs_call', 'contig_list')}

[extract]
cores = {global_conf.global_get('gembs_extract', 'cores')}
extract_threads = {global_conf.global_get('gembs_extract', 'extract_threads')}
threads = {global_conf.global_get('gembs_extract', 'threads')}
memory = {global_conf.global_get('gembs_extract', 'ram')}

make_cpg = {global_conf.global_get('gembs_extract', 'make_cpg')}
make_non_cpg = {global_conf.global_get('gembs_extract', 'make_non_cpg')}
make_bedmethyl = {global_conf.global_get('gembs_extract', 'make_bedmethyl')}
make_snps = {global_conf.global_get('gembs_extract', 'make_snps')}"""

    return Job(
            output_files = [output],
            command="""\
echo \"{config_content}\" > {config_file}""".format(
        config_content=config_content,
        config_file=output
            )
        )

def make_metadata(metadata_list, metadata_file):
    
    return Job(
            output_files = [metadata_file],
            command="""\
echo -e 'sampleID,dataset,library,sample,file1,file2 {metadata_list}' > {metadata_file}""".format(
    metadata_list = " ".join(["\\n" + metadata for metadata in metadata_list]),
    metadata_file=metadata_file
                )
            )

def prepare(
        metadata, 
        config_file, 
        output_dir,
        ini_section="gembs_prepare"
        ):

    output = [
            output_dir + "/.gemBS/gemBS.mp",
            ]
    
    return Job(
        [metadata,config_file],
        output,
        [
            [ini_section, 'module_gembs']
        ],
        command="""\
rm -rf {hidden_dir}
gemBS {gembs_flags} {gembs_options} \\
  prepare {flags} {options} \\
  --config {config_file} \\
  --text-metadata {metadata}""".format(
      hidden_dir=output_dir + "/.gemBS",
      gembs_flags=global_conf.global_get(ini_section, 'gembs_flags', required=False),
      gembs_options=global_conf.global_get(ini_section, 'gembs_options', required=False),
      flags=global_conf.global_get(ini_section, 'flags', required=False),
      options=global_conf.global_get(ini_section, 'options', required=False),
      config_file=config_file,
      metadata=metadata
      )
    )

def index(input, output):

    output_info = re.sub(".gem", ".info", output)
    return Job(
        [input],
        [output, output_info],
        [
            ['gembs_index', 'module_gembs']
        ],
        command="""\
gemBS {gembs_flags} {gembs_options} \\
  index {flags} {options}""".format(
    gembs_flags=global_conf.global_get('gembs_index', 'gembs_flags', required=False),
    gembs_options=global_conf.global_get('gembs_index', 'gembs_options', required=False),
    flags=global_conf.global_get('gembs_index', 'flags', required=False),
    options=global_conf.global_get('gembs_index', 'options', required=False)
            )
        )

def map(
        sample, 
        output_dir,
        ini_section="gembs_map"
        ):
    outputs = [
            output_dir + "/" + sample + ".bam",
            output_dir + "/" + sample + ".bam.csi",
            output_dir + "/" + sample + ".bam.md5",
            output_dir + "/" + sample + ".json"
            ]

    return Job(
        [os.path.join(output_dir, ".gemBS", "gemBS.mp")],
        outputs,
        [
            [ini_section, 'module_gembs']
        ],
        command="""\
gemBS {gembs_flags} {gembs_options} \\
  --dir {working_dir} \\
  map {flags} {options} \\
  --barcode {sample} \\
  --tmp-dir {tmp_dir}""".format(
      gembs_flags=global_conf.global_get(ini_section, 'gembs_flags', required=False),
      gembs_options=global_conf.global_get(ini_section, 'gembs_options', required=False),
      working_dir=output_dir,
      flags=global_conf.global_get(ini_section, 'flags', required=False),
      options=global_conf.global_get(ini_section, 'options', required=False),
      sample=sample,
      tmp_dir=global_conf.global_get(ini_section, 'tmp_dir')
      )
    )

def call(
        sample, 
        input, 
        output_prefix,
        ini_section="gembs_call"
        ):
    outputs = [
            output_prefix + ".json",
            output_prefix + ".bcf",
            output_prefix + ".bcf.csi",
            output_prefix + ".bcf.md5"
            ]

    return Job(
        [input],
        outputs,
        [
            [ini_section, 'module_gembs']
        ],
        command = """\
gemBS {gembs_flags} {gembs_options} \\
  --dir {working_dir} \\
  call {flags} {options} \\
  --barcode {sample} \\
  --tmp-dir {tmp_dir} {dbSNP}""".format(
      gembs_flags=global_conf.global_get(ini_section, 'gembs_flags', required=False),
      gembs_options=global_conf.global_get(ini_section, 'gembs_options', required=False),
      working_dir=os.path.dirname(output_prefix),
      flags=global_conf.global_get(ini_section, 'flags', required=False),
      options=global_conf.global_get(ini_section, 'options', required=False),
      sample=sample,
      tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
      dbSNP="-D " + global_conf.global_get(ini_section, 'dbSNP_index') if global_conf.global_get(ini_section, 'dbSNP_index', required=False) else ""
      )
    )

def extract(
        input, 
        sample, 
        output_dir,
        ini_section="gembs_extract"
        ):
    output = [
            output_dir + "/" + sample + "_cpg.bb",
            output_dir + "/" + sample + "_cpg.bed.gz"
            ]

    return Job(
        [input],
        output,
        [
            [ini_section, 'module_gembs']
        ],
        command="""\
gemBS {gembs_flags} {gembs_options} \\
  --dir {working_dir} \\
  extract {flags} {options} \\
  --barcode {sample}""".format(
    gembs_flags=global_conf.global_get(ini_section, 'gembs_flags', required=False),
    gembs_options=global_conf.global_get(ini_section, 'gembs_options', required=False),
    working_dir=output_dir,
    flags=global_conf.global_get(ini_section, 'flags', required=False),
    options=global_conf.global_get(ini_section, 'options', required=False),
    sample=sample
    )
)

def report(
        inputs, 
        output,
        ini_section="gembs_report"
        ):

    output_dir = os.path.dirname(output)
    outputs = [
            output,
            re.sub(".html", ".tex", output),
            output_dir + "/mapping/map_report.tex",
            output_dir + "/calling/call_report.tex"
            ]
    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_gembs']
        ],
        command="""\
gemBS {gembs_flags} {gembs_options} \\
  report {flags} {options}""".format(
    gembs_flags=global_conf.global_get(ini_section, 'gembs_flags', required=False),
    gembs_options=global_conf.global_get(ini_section, 'gembs_options', required=False),
    flags=global_conf.global_get(ini_section, 'flags', required=False),
    options=global_conf.global_get(ini_section, 'options', required=False),
            )
        )
