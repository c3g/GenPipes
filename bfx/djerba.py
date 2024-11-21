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
from core.job import Job

def make_config(
        output,
        tumor_pair_name,
        tumor_id,
        normal_id,
        maf_input,
        purple_input,
        assay="WGTS",
        ini_section = 'report_djerba'
        ):

    config_content = f"""\
[core]
archive_name = djerba
archive_url = http://$username:$password@$address:$port
attributes = research 
author = C3G Author
configure_priority = 100
depends_configure = 
depends_extract = 
document_config = document_config.json
extract_priority = 100
render_priority = 100
report_id = {tumor_pair_name}
report_version = 1

[provenance_helper]
attributes = research
assay = {assay}
tumour_id = {tumor_id}
normal_id = {normal_id}
sample_name_tumour = None
sample_name_normal = None
sample_name_aux = None
project = {config.param(ini_section, 'project_name')}
donor = {tumor_pair_name}
provenance_input_path = provenance_input.tsv.gz

[case_overview]
assay = testing
attributes = research
assay_description = Whole Genome Sequencing - Tumor Pair Pipeline
primary_cancer = {config.param(ini_section, 'cancer_type') if config.param(ini_section, 'cancer_type') else "Unknown"}
site_of_biopsy = Unknown
study = {config.param(ini_section, 'project_name', required = False) if config.param(ini_section, 'project_name') else "MoH"}
patient_study_id = {tumor_pair_name}
donor = {tumor_pair_name}
tumour_id = {tumor_id}
normal_id = {normal_id}
requisition_approved = 2185-07-18

[treatment_options_merger]
attributes = research,supplementary
configure_priority = 300
depends_configure = 
render_priority = 50

[wgts.snv_indel]
apply cache = {config.param(ini_section, 'apply_cache') if config.param(ini_section, 'apply_cache') else "False"}
update cache = {config.param(ini_section, 'update_cache') if config.param(ini_section, 'update_cache') else "False"}
oncokb cache = {config.param(ini_section, 'oncokb_cache') if config.param(ini_section, 'oncokb_cache') else ""}
attributes = research
configure_priority = 700
depends_configure = 
depends_extract = 
extract_priority = 800
render_priority = 700
maf_path = {maf_input}
oncotree_code = {config.param(ini_section, 'cancer_type') if config.param(ini_section, 'cancer_type') else ""}
tumour_id = {tumor_id}
normal_id = {normal_id}
whizbam_project = COL

[wgts.cnv_purple]
apply cache = {config.param(ini_section, 'apply_cache') if config.param(ini_section, 'apply_cache') else "False"}
update cache = {config.param(ini_section, 'update_cache') if config.param(ini_section, 'update_cache') else "False"}
oncokb cache = {config.param(ini_section, 'oncokb_cache') if config.param(ini_section, 'oncokb_cache') else ""}
attributes = research
configure_priority = 900
tumour_id = {tumor_id}
oncotree_code = {config.param(ini_section, 'cancer_type') if config.param(ini_section, 'cancer_type') else ""}
purple_zip = {purple_input}
whizbam_project=OCTCAP
assay = {assay}

[gene_information_merger]
attributes = research,supplementary
configure_priority = 2000
depends_configure = 
render_priority = 2000

[supplement.body]
assay = {assay}
attributes = research
depends_configure = 
depends_extract = 
configure_priority = 1200
extract_priority = 1200
render_priority = 1200
{"custom_template_dir = " + config.param(ini_section, 'custom_template_dir') if config.param(ini_section, 'custom_template_dir') else ""}
"""
    return Job(
        output_files = [output],
        command="""\
echo \"{config_content}\" > {config_file}""".format(
    config_content=config_content,
    config_file=output
        )
    )

def clean_maf(
        input_maf,
        output_maf
        ):
    
    output = [output_maf + ".gz"]
    return Job(
        [input_maf],
        [output],
        [],
        command="""\
col=\\$( awk -v RS='\t' '/t_depth/{{print NR; exit}}' {input_maf} ) && \\
awk -F'\t' -v col=\\$col '! ( \\$col=="" )' {input_maf} > {output_maf} && \\
gzip {output_maf}""".format(
        input_maf=input_maf,
        output_maf=output_maf
        )
    )

def make_script(
        config_file,
        output_dir,
        output_file,
        ini_section = 'report_djerba'
        ):
    
    return Job(
        [config_file],
        [output_file],
        [],
        command="""\
echo "module purge && module load {module_djerba} {module_wkhtmltopdf} && \\
export ONCOKB_TOKEN={oncokb_token} \\

\\${DJERBA_BIN}/djerba.py {djerba_options} report \\
    -i {config_file} \\
    -o {output_dir} \\
    -p --no-archive" > {output_file}""".format(
        module_djerba=config.param(ini_section, 'module_djerba'),
        module_wkhtmltopdf=config.param(ini_section, 'module_wkhtmltopdf'),
        oncokb_token=config.param(ini_section, 'oncokb_token'),
        djerba_options=config.param(ini_section, 'djerba_options', required=False),
        config_file=config_file,
        output_dir=output_dir,
        output_file=output_file
        )
    )