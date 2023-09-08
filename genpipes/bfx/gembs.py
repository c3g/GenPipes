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

index_dir = {genpipes_dir}/alignment/index
sequence_dir = {genpipes_dir}
bam_dir = {genpipes_dir}/alignment/@BARCODE
bcf_dir = {genpipes_dir}/methylation_call/@BARCODE
extract_dir = {genpipes_dir}/variants/@BARCODE
report_dir = {genpipes_dir}/report

project = GenPipes
species = {global_conf.global_get('default', 'scientific_name')}

# Default parameters
threads = 1
cores = 1
jobs = 1

[index]
threads = {global_conf.global_get('gembs_index', 'threads')}
memory = {global_conf.global_get('gembs_index', 'ram')}

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

[extract]
cores = {global_conf.global_get('gembs_extract', 'cores')}
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
   # with open(output, "w") as config_file:
   #     config_file.write(config_content)

def prepare(metadata, config_file, output_dir):

    prefix = os.path.join(output_dir, "alignment/index", global_conf.global_get('default', 'scientific_name') + ".gemBS.")
    output = [
            output_dir + "/.gemBS/gemBS.mp",
            prefix + "ref",
            prefix + "ref.fai",
            prefix + "ref.gzi",
            prefix + "contig_sizes",
            prefix + "contig_md5"
            ]
    
    return Job(
        [metadata,config_file],
        output,
        [
            ['gembs_prepare', 'module_htslib'],
            ['gembs_prepare', 'module_gembs']
        ],
        command="""\
gemBS prepare {flags} {options} \\
  --config {config_file} \\
  --text-metadata {metadata}""".format(
      flags=global_conf.global_get('gembs_prepare', 'flags', required=False),
      options=global_conf.global_get('gembs_prepare', 'options', required=False),
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
            ['gembs_index', 'module_gembs'],
            ['gembs_index', 'module_htslib']
        ],
        command="""\
gemBS index {flags} {options}""".format(
    flags=global_conf.global_get('gembs_index', 'flags', required=False),
    options=global_conf.global_get('gembs_index', 'options', required=False)
            )
        )

def map(sample, gembs_config, index, output_prefix):
    outputs = [
            output_prefix + ".bam",
            output_prefix + ".bam.csi",
            output_prefix + ".bam.md5",
            output_prefix + ".json"
            ]

    return Job(
        [gembs_config, index],
        outputs,
        [
            ['gembs_map', 'module_gembs'],
            ['gembs_map', 'module_samtools'],
            ['gembs_map', 'module_htslib']
        ],
        command="""\
gemBS map {flags} {options} \\
  --barcode {sample} \\
  --tmp-dir {tmp_dir}""".format(
      flags=global_conf.global_get('gembs_map', 'flags', required=False),
      options=global_conf.global_get('gembs_map', 'options', required=False),
      sample=sample,
      tmp_dir=global_conf.global_get('gembs_map', 'tmp_dir')
      )
    )

def call(sample, input, output_prefix):
    outputs = [
            output_prefix + "_call.json",
            output_prefix + ".bcf",
            output_prefix + ".bcf.csi",
            output_prefix + ".bcf.md5"
            ]

    return Job(
        [input],
        outputs,
        [
            ['gembs_call', 'module_gembs'],
            ['gembs_call', 'module_htslib']
        ],
        command = """\
gemBS call {flags} {options} \\
  --barcode {sample} \\
  --tmp-dir {tmp_dir}""".format(
      flags=global_conf.global_get('gembs_call', 'flags', required=False),
      options=global_conf.global_get('gembs_call', 'options', required=False),
      sample=sample,
      tmp_dir=global_conf.global_get('gembs_call', 'tmp_dir')
      )
    )

def extract(input, sample, output_dir):
    output = [
            output_dir + "/" + sample + "_cpg.bb",
            output_dir + "/" + sample + "_cpg.bed.gz"
            ]

    return Job(
        [input],
        output,
        [
            ['gembs_extract', 'module_gembs']
        ],
        command="""\
gemBS {gembs_flags} {gembs_options} \\
  --dir {working_dir} \\
  extract {flags} {options} \\
  --barcode {sample}""".format(
    gembs_flags=global_conf.global_get('gembs_extract', 'gembs_flags', required=False),
    gembs_options=global_conf.global_get('gembs_extract', 'gembs_options', required=False),
    working_dir=output_dir,
    flags=global_conf.global_get('gembs_extract', 'flags', required=False),
    options=global_conf.global_get('gembs_extract', 'options', required=False),
    sample=sample
    )
)

def report(inputs, output):

    return Job(
        inputs,
        [output],
        [
            ['gembs_report', 'module_gembs'],
            ['gembs_report', 'module_htslib']
        ],
        command="""\
gemBS report {flags} {options}""".format(
    flags=global_conf.global_get('gembs_report', 'flags', required=False),
    options=global_conf.global_get('gembs_report', 'options', required=False),
            )
        )
