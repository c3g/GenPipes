################################################################################
# Copyright (C) 2025 C3G, The Victor Phillip Dahdaleh Institute of Genomic Medicine at McGill University
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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

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
project = {global_conf.global_get(ini_section, 'project_name', required = False)}
donor = {tumor_pair_name}
provenance_input_path = provenance_input.tsv.gz

[case_overview]
assay = {assay}
attributes = research
assay_description = WGS - DnaSeq Pipeline - Somatic Ensemble Protocol
primary_cancer = {global_conf.global_get(ini_section, 'cancer_type', required = False) if global_conf.global_get(ini_section, 'cancer_type', required = False) else "Unknown"}
site_of_biopsy = Unknown
study = {global_conf.global_get(ini_section, 'project_name', required = False) if global_conf.global_get(ini_section, 'project_name', required = False) else "Unknown"}
patient_study_id = {tumor_pair_name}
donor = {tumor_pair_name}
tumour_id = {tumor_id}
normal_id = {normal_id}
requisition_approved = 0000-00-00

[treatment_options_merger]
attributes = research,supplementary
configure_priority = 300
depends_configure = 
render_priority = 50

[wgts.snv_indel]
apply cache = {global_conf.global_get(ini_section, 'apply_cache', required = False) if global_conf.global_get(ini_section, 'apply_cache', required = False) else "False"}
update cache = {global_conf.global_get(ini_section, 'update_cache', required = False) if global_conf.global_get(ini_section, 'update_cache', required = False) else "False"}
oncokb cache = {global_conf.global_get(ini_section, 'oncokb_cache', required = False) if global_conf.global_get(ini_section, 'oncokb_cache', required = False) else ""}
attributes = research
configure_priority = 700
depends_configure = 
depends_extract = 
extract_priority = 800
render_priority = 700
maf_path = {maf_input}
oncotree_code = {global_conf.global_get(ini_section, 'cancer_type', required = False) if global_conf.global_get(ini_section, 'cancer_type', required = False) else ""}
tumour_id = {tumor_id}
normal_id = {normal_id}
whizbam_project = COL

[wgts.cnv_purple]
apply cache = {global_conf.global_get(ini_section, 'apply_cache', required = False) if global_conf.global_get(ini_section, 'apply_cache', required = False) else "False"}
update cache = {global_conf.global_get(ini_section, 'update_cache', required = False) if global_conf.global_get(ini_section, 'update_cache', required = False) else "False"}
oncokb cache = {global_conf.global_get(ini_section, 'oncokb_cache', required = False) if global_conf.global_get(ini_section, 'oncokb_cache', required = False) else ""}
attributes = research
configure_priority = 900
tumour_id = {tumor_id}
oncotree_code = {global_conf.global_get(ini_section, 'cancer_type', required = False) if global_conf.global_get(ini_section, 'cancer_type', required = False) else ""}
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
{"template_dir = " + global_conf.global_get(ini_section, 'template_dir', required = False) if global_conf.global_get(ini_section, 'template_dir', required = False) else ""}
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
col=$( awk -v RS='\\t' '/t_depth/{{print NR; exit}}' {input_maf} ) && \\
awk -F'\\t' -v col=$col '! ( $col=="" )' {input_maf} > {output_maf} && \\
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
echo "module purge && \\
module load {module_djerba} {module_wkhtmltopdf} && \\
export ONCOKB_TOKEN={oncokb_token} && \\
export DJERBA_CORE_HTML_DIR={html_directory} && \\

djerba.py {djerba_options} report \\
    -i {config_file} \\
    -o {output_dir} \\
    -p --no-archive" > {output_file}""".format(
        module_djerba=global_conf.global_get(ini_section, 'module_djerba'),
        module_wkhtmltopdf=global_conf.global_get(ini_section, 'module_wkhtmltopdf'),
        oncokb_token=global_conf.global_get(ini_section, 'oncokb_token'),
        html_directory=global_conf.global_get(ini_section, 'custom_html_directory'),
        djerba_options=global_conf.global_get(ini_section, 'djerba_options', required=False),
        config_file=config_file,
        output_dir=output_dir,
        output_file=output_file
        )
    )