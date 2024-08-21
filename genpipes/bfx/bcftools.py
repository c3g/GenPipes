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

# GenPipes Modules
from ..core.config import global_conf
from ..core.job import Job

def add_reject(input, output):
    """
    Adds REJECT to the filter. Used with scalpel when merging common indels
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_add_reject', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v --soft-filter 'REJECT' -e '%TYPE="indel"' \\
  {input}{output}""".format(
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def add_chi2Filter(input, output):
    """
    Adds CHI2FILTER to the filter. Used with scalpel on somatic indels
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_add_chi2Filter', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v --soft-filter 'CHI2FILTER' -e 'INFO/CHI2 > 20.0' \\
  {input}{output}""".format(
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def mpileup(inputs, output, options=None, regions=None, regionFile=None):
    """
    New bcftools mpileup function
    """
    if not isinstance(inputs, list):
        inputs = [inputs]
        
    return Job(
        inputs,
        [output],
        [
            ['bcftools_mpileup', 'module_bcftools']
        ],
        command="""\
bcftools \\
  mpileup {options} \\
  -f {reference_fasta} \\
  {inputs} \\
  {regions} \\
  {regionFile} \\
  {output}""".format(
        options=options if options else "",
        reference_fasta=global_conf.global_get('samtools_mpileup', 'genome_fasta', param_type='filepath'),
        inputs="".join(" \\\n  " + input for input in inputs),
        regions="-r " + regions if regions else "",
        regionFile="-R " + regionFile if regionFile else "",
        output=" \\\n > " + output if output else ""
        )
    )


def call(inputs, output, options=None):
    """
    New bcftools call function
    """
    return Job(
        inputs,
        [output],
        [
            ['bcftools_call', 'module_bcftools']
        ],
        command="""\
bcftools \\
  call {options} \\
  {output} \\
  {inputs}""".format(
        options=options if options else "",
        inputs="".join(" \\\n  " + input for input in inputs),
        output=" \\\n -o " + output if output else ""
        )
    )

def index(inputs, options=None):
    """
    New bcftools index function
    """
    output = inputs + ".csi"

    return Job(
        [inputs],
        [output],
        [
            ['bcftools_index', 'module_bcftools']
        ],
        command="""\
bcftools \\
  index -f {options} \\
  {inputs}""".format(
        options=options if options else "",
        inputs=inputs,
        )
    )

def concat(inputs, output, options=None):
    """
    Concatenate or combine VCF/BCF files
    """
    return Job(
        inputs,
        [output],
        [
            ['bcftools_concat', 'module_bcftools'],
        ],
        command="""\
bcftools \\
  concat -a {options} \\
  {output} \\
  {inputs}""".format(
        options=options if options else "",
        inputs="".join(" \\\n  " + input for input in inputs),
        output=" \\\n -o " + output if output else ""
        )
    )

def view(input, output, filter_options=None):
    """
    Generalized view 
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_view', 'module_bcftools']
        ],
        command="""\
bcftools \\
  view {filter_options} \\
  {output}{input}""".format(
        filter_options=filter_options if filter_options else "",
        input=" \\\n " + input if input else "",
        output=" \\\n -o " + output if output else ""
        )
    )

def sort(input, output, sort_options=None):
    """
    Generalized sort
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_sort', 'module_bcftools']
        ],
        command="""\
bcftools \\
  sort {sort_options} \\
  {output}{input}""".format(
        sort_options=sort_options if sort_options else "",
        input=" \\\n " + input if input else "",
        output=" \\\n -o " + output if output else ""
        )
    )

def query(input, output, query_options=None):
    """
    Generalized sort
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_sort', 'module_bcftools']
        ],
        command="""\
bcftools \\
  query {query_options} \\
  {output}{input}""".format(
        query_options=query_options if query_options else "",
        input=" \\\n " + input if input else "",
        output=" \\\n -o " + output if output else ""
        )
    )

def filter(input, output, filter_options):
    """
    Generalized filter function
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_filter', 'module_bcftools']
        ],
        command="""\
bcftools \\
  filter -m '+' -O v {filter_options} \\
  {input}{output}""".format(
        input=" \\\n " + input if input else "",
        filter_options=filter_options,
        output=" \\\n  > " + output if output else ""
        )
    )

def annotate(input, output, options):
    """
    Generalized merge
    """
    return Job(
        [input],
        [output],
        [
            ['bcftools_annotate', 'module_bcftools']
        ],
        command="""\
bcftools \\
  annotate {options} \\
  {input}{output}""".format(
        options=options if options else "",
        input=" \\\n " + input if input else "",
        output=" \\\n -o " + output if output else ""
        )
    )

def consensus(input, output, options, ini_section='bcftools_consensus'):
    """
    Consensus creation
    """
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools \\
  consensus {options} \\
  {input}\\
  {output}""".format(
        options=options if options else "",
        input=input,
        output="> " + output if output else ""
        )
    )

def norm(input, output, options, ini_section='bcftools_norm'):
    """
    VCF normalisation
    """
    return Job(
        [input],
        [output],
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools \\
  norm {options} \\
  {input} \\
  {output}""".format(
        options=options if options else "",
        input=input,
        output="> " + output if output else ""
        )
    )

def split(input, output, options, ini_section='bcftools_split'):
    """
    split vcf by sample
    """
    
    return Job(
        [input],
        output,
        [
            [ini_section, 'module_bcftools']
        ],
        command="""\
bcftools \\
  +split {options} \\
  {input} \\
  -o {output}""".format(
        options=options if options else "",
        input=input,
        output=output
        )
    )
