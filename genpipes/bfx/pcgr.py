################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job


def report(input_vcf,
           cpsr_report,
           output_dir,
           tumor_id,
           input_cna,
           ini_section='report_pcgr'
           ):

    if global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] >= "1":
        call = 'pcgr'
    else:
        call = 'pcgr.py'
    if global_conf.global_get('report_pcgr', 'module_pcgr').split("/")[2] >= "2":
        return report2(input_vcf, cpsr_report, output_dir, tumor_id, input_cna, ini_section=ini_section)
    else:
        tumor_id = tumor_id[:35]

        return Job(
            [
                input_vcf,
                input_cna
            ],
            output_dir,
            [
                [ini_section, 'module_pcgr'],
            ],
            command="""\
if [ -e {input_cna}.pass ]; then
    export input_cna="--input_cna {input_cna}"
 else
    export input_cna=""
fi && \\
{call} {options} \\
    {tumor_type} \\
    {assay} \\
    {tumor_options} \\
    {normal_options} \\
    {mutsig_options} \\
    {tmb_options} \\
    {msi_options} \\
    --input_vcf {input_vcf} \\
    --cpsr_report {cpsr_report} \\
    $input_cna \\
    --pcgr_dir $PCGR_DATA \\
    --output_dir {output_dir} \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id}""".format(
            call=call,
            options=global_conf.global_get(ini_section, 'options'),
            tumor_type=global_conf.global_get(ini_section, 'tumor_type'),
            assay=global_conf.global_get(ini_section, 'assay'),
            tumor_options=global_conf.global_get(ini_section, 'tumor_options'),
            normal_options=global_conf.global_get(ini_section, 'normal_options'),
            mutsig_options=global_conf.global_get(ini_section, 'mutsig_options'),
            tmb_options=global_conf.global_get(ini_section, 'tmb_options'),
            msi_options=global_conf.global_get(ini_section, 'msi_options'),
            input_vcf=input_vcf,
            cpsr_report=cpsr_report,
            input_cna=input_cna,
            output_dir=output_dir,
            assembly=global_conf.global_get(ini_section, 'assembly'),
            tumor_id=tumor_id
        )
    )

def report2(input_vcf,
           cpsr_prefix,
           output_dir,
           tumor_id,
           input_cna,
           ini_section='report_pcgr'
           ):
    
    if cpsr_prefix:
        cpsr_input = f"{cpsr_prefix}.classification.tsv.gz"
        cpsr_yaml = f"{cpsr_prefix}.conf.yaml"
    # use tmp dir for pcgr to avoid disk quota issues caused by bcftools tmp dir settings
    return Job(
        [
            input_vcf,
            input_cna
        ],
        output_dir,
        [
            [ini_section, 'module_pcgr'],
        ],
            command="""\
if [ -e {input_cna}.pass ]; then
    export input_cna="--input_cna {input_cna}"
 else
    export input_cna=""
fi && \\
mkdir -p {tmp_dir}/pcgr && \\
pcgr {options} \\
    {tumor_type} \\
    {assay} \\
    {tumor_options} \\
    {normal_options} \\
    {mutsig_options} \\
    {tmb_options} \\
    {msi_options} \\
    --input_vcf {input_vcf} \\
    {cpsr_report} {cpsr_yaml} \\
    $input_cna \\
    --refdata_dir $PCGR_DATA \\
    --vep_dir $PCGR_VEP_CACHE \\
    --output_dir {tmp_dir}/pcgr \\
    --genome_assembly {assembly} \\
    --sample_id {tumor_id} && \\
cp -r {tmp_dir}/pcgr {output_dir}""".format(
            options=global_conf.global_get(ini_section, 'options_v2'),
            tumor_type=global_conf.global_get(ini_section, 'tumor_type'),
            assay=global_conf.global_get(ini_section, 'assay'),
            tumor_options=global_conf.global_get(ini_section, 'tumor_options'),
            normal_options=global_conf.global_get(ini_section, 'normal_options'),
            mutsig_options=global_conf.global_get(ini_section, 'mutsig_options'),
            tmb_options=global_conf.global_get(ini_section, 'tmb_options'),
            msi_options=global_conf.global_get(ini_section, 'msi_options'),
            input_vcf=input_vcf,
            cpsr_report="--input_cpsr " + cpsr_input if cpsr_input else "",
            cpsr_yaml="--input_cpsr_yaml " + cpsr_yaml if cpsr_yaml else "", 
            input_cna=input_cna,
            tmp_dir=global_conf.global_get(ini_section, 'tmp_dir'),
            output_dir=os.path.dirname(output_dir),
            assembly=global_conf.global_get(ini_section, 'assembly'),
            tumor_id=tumor_id
        )
    )

def create_header(output):
    return Job(
        command=f"""\
`cat > {output} << END
Chromosome\tStart\tEnd\tSegment_Mean\tnMajor\tnMinor
END`"""
        )

def create_input_cna(
        cna_body,
        header,
        cnvkit_calls, 
        output
        ):
    return Job(
        [cna_body, cnvkit_calls],
        [output],
        [],
        command=f"""\
 cat {header} > {output}
 while read line; do 
    LOCUS=$(echo $line | awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3}}')
    grep "$LOCUS" {cnvkit_calls} | awk 'BEGIN {{OFS="\\t"}} {{print $1, $2, $3, $5, $8, $9}}' >> {output}; done < {cna_body}"""
)

def parse_pcgr_passed_variants_pt(input_file):
    """
    Parse PCGR passed variants.
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export pcgr_passed_variants=`grep "pcgr-gene-annotate - INFO - Number of PASSed variant calls:" {input_file} | awk -F': ' '{{print $2}}'`"""
        )
