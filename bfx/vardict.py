################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
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

# MUGQIC Modules
from core.config import *
from core.job import *

def paired(
        input_normal,
        input_tumor,
        tumor_name,
        output=None,
        region=[],
        ini_section='vardict_paired'
):
    return Job(
        [input_normal, input_tumor, region],
        [output],
        [
        [ini_section, 'module_vardict'],
        [ini_section, 'module_samtools'],
        [ini_section, 'module_perl'],
        [ini_section, 'module_R']
        ],
        command="""\
vardict \\
  -G {reference} \\
  -N {tumor_name} \\
  -b "{paired_samples}" \\
  {vardict_options}{region}{output}""".format(
        reference=config.param(ini_section, 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        vardict_options=config.param(ini_section, 'vardict_options'),
        region=" \\\n  " + region if region else "",
        output=" \\\n  > " + output if output else ""
        )
    )

def paired_java(
        input_normal,
        input_tumor,
        tumor_name,
        output=None,
        region=[],
        ini_section='vardict_paired'
):
    return Job(
        [input_normal, input_tumor, region],
        [output],
        [
        [ini_section, 'module_java'],
        [ini_section, 'module_vardict_java'],
        [ini_section, 'module_samtools'],
        [ini_section, 'module_perl'],
        [ini_section, 'module_R']
        ],
        command="""\
java {java_other_options} -Djava.io.tmpdir={tmp_dir} -Xms768m -Xmx{ram} -classpath {classpath} \\
  -G {reference} \\
  -N {tumor_name} \\
  -b "{paired_samples}" \\
  {vardict_options}{region}{output}""".format(
        tmp_dir=config.param(ini_section, 'tmp_dir'),
        reference=config.param(ini_section, 'genome_fasta', param_type='filepath'),
        tumor_name=tumor_name,
        paired_samples=input_tumor + "|" + input_normal,
        java_other_options=config.param(ini_section, 'java_other_options'),
        ram=config.param(ini_section, 'ram'),
        classpath=config.param(ini_section, 'classpath'),
        vardict_options=config.param(ini_section, 'vardict_options'),
        region=" \\\n  " + region if region else "",
        output=" \\\n  > " + output if output else ""
        )
    )

def testsomatic(
        input=None,
        output=None,
        ini_section='vardict_paired'
):
    return Job(
        [input],
        [output],
        [
        [ini_section, 'module_vardict_java'],
        [ini_section, 'module_R']
        ],
        command="""\
$VARDICT_BIN/testsomatic.R {input} {output}""".format(
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
        )
    )


def var2vcf(
        output,
        normal_name,
        tumor_name,
        input=None,
        ini_section='vardict_paired'
):
    return Job(
        [input],
        [output],
        [
        [ini_section, 'module_vardict_java'],
        [ini_section, 'module_perl']
        ],
        command="""\
perl $VARDICT_BIN/var2vcf_paired.pl \\
    -N "{pairNames}" \\
    {var2vcf_options}{input}{output}""".format(
        pairNames=tumor_name + "|" + normal_name,
        var2vcf_options=config.param(ini_section, 'var2vcf_options'),
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
        )
    )

def dict2beds(
        dictionary,
        beds,
        ini_section='vardict_paired'
):
    return Job(
        [dictionary],
        beds,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_python']
        ],
        command="""\
dict2BEDs.py \\
  --dict {dictionary} \\
  --beds {beds} {dict2bed_options}""".format(
        dictionary=dictionary if dictionary else config.param(ini_section, 'genome_dictionary', param_type='filepath'),
        beds=' '.join(beds),
        dict2bed_options=config.param(ini_section, 'dict2bed_options')
        )
    )

