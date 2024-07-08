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
from ..core.job import Job, concat_jobs

def compute_effects(
        input,
        output,
        split=False,
        cancer_sample_file=None,
        options=[],
        ini_section='compute_effects'
):
    
    output_stats = f"{output}.stats.csv"
    output_stats_html = f"{output}.stats.html"
    
    outputs = []
    if not isinstance(output, list):
        outputs = [output]
    
    outputs += [output_stats, output_stats_html]
    
    job = Job(
        [input],
        outputs,
        [
            ['compute_effects', 'module_java'],
            ['compute_effects', 'module_snpeff']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $SNPEFF_HOME/snpEff.jar eff {options} \\
  {cancer_sample_file} \\
  -c $SNPEFF_HOME/snpEff.config \\
  -i vcf \\
  -o vcf \\
  -csvStats {output_stats} \\
  -stats {output_stats_html} \\
  {reference_snpeff_genome} \\
  {input}{output}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        options=options if options else "",
        cancer_sample_file="-cancerSamples " + cancer_sample_file if cancer_sample_file else "",
        output_stats=output_stats,
        output_stats_html=output_stats_html,
        reference_snpeff_genome=global_conf.global_get(ini_section, 'snpeff_genome'),
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

    if split:
        split_output_stats = output + ".statsFile.txt"
        split_job = Job(
            [output_stats],
            [split_output_stats],
            [['compute_effects', 'module_mugqic_tools']],
            command="""\
splitSnpEffStat.awk \\
  {output_stats} \\
  {output_part} \\
  {split_output_stats}""".format(
            output_stats=output_stats,
            output_part=output + ".part",
            split_output_stats=split_output_stats
            )
        )

        job = concat_jobs([job, split_job])

    return job

def snpsift_annotate(input, output):
    return Job(
        [input],
        [output],
        [
            ['snpsift_annotate', 'module_java'],
            ['snpsift_annotate', 'module_snpeff']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $SNPEFF_HOME/SnpSift.jar annotate \\
  {db_snp} \\
  {input}{output}""".format(
        tmp_dir=global_conf.global_get('snpsift_annotate', 'tmp_dir'),
        java_other_options=global_conf.global_get('snpsift_annotate', 'java_other_options'),
        ram=global_conf.global_get('snpsift_annotate', 'ram'),
        db_snp=global_conf.global_get('snpsift_annotate', 'known_variants', param_type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
        ),
        removable_files=[output]
    )

def snpsift_dbnsfp(input, output, ini_section='dbnsfp_annotation'):
    return Job(
        [input],
        [output],
        [
            [ ini_section, 'module_java'],
            [ ini_section, 'module_snpeff']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $SNPEFF_HOME/SnpSift.jar dbnsfp \\
  -v -db {db_nsfp} \\
  {input}{output}""".format(
        tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
        java_other_options=global_conf.global_get(ini_section, 'java_other_options'),
        ram=global_conf.global_get(ini_section, 'ram'),
        db_nsfp=global_conf.global_get(ini_section, 'dbnsfp', param_type='filepath'),
        input=input,
        output=" \\\n  > " + output if output else ""
        )
    )

def snpsift_intervals_index(input, intervals_file, output=None, job_name="filter_vcf_snpsift"):
    return Job(
        [input, intervals_file],
        [output] if output else [],
        [
            ['snpsift_dbnsfp', 'module_java'],
            ['snpsift_dbnsfp', 'module_snpeff']
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $SNPEFF_HOME/SnpSift.jar intidx \\
    {input} {intervals_file}{output}""".format(
            tmp_dir=global_conf.global_get('snpsift_dbnsfp', 'tmp_dir', required=False),
            java_other_options=global_conf.global_get('snpsift_dbnsfp', 'java_other_options', required=False),
            ram=global_conf.global_get('snpsift_dbnsfp', 'ram', required=False),
            input=input,
            intervals_file=intervals_file,
            output=" \\\n  > " + output if output else ""
        ),
        name=job_name
    )

def snpeff_annotate(input_vcf, output_vcf, metrics_prefix):
    inputs = [input_vcf]
    output = [output_vcf, metrics_prefix + ".html", metrics_prefix + ".genes.txt"]

    return Job(
        input_files=inputs,
        output_files=output,
        module_entries=[
            ['snpeff_annotate', 'module_snpeff'],
            ['snpeff_annotate', 'module_java']
        ],

        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} \\
  -jar $SNPEFF_HOME/snpEff.jar ann \\
  -v -c $SNPEFF_HOME/snpEff.config \\
  -s {metrics_prefix} \\
  {reference_snpeff_genome} \\
  {input_vcf} > \\
  {output_vcf}""".format(
            tmp_dir=global_conf.global_get('snpeff_annotate', 'tmp_dir', required=False),
            java_other_options=global_conf.global_get('snpeff_annotate', 'java_other_options', required=False),
            ram=global_conf.global_get('snpeff_annotate', 'ram', required=False),
            metrics_prefix=metrics_prefix + ".html",
            reference_snpeff_genome=global_conf.global_get('snpeff_annotate', 'snpeff_genome'),
            input_vcf=input_vcf,
            output_vcf=output_vcf
            ),
        )
