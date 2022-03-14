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
import logging

# MUGQIC Modules
from core.config import *
from core.job import *

log = logging.getLogger(__name__)


def run_ncovtools(output_filtered_bam,
    output_primer_trimmed_bam,
    output_consensus,
    output_variants,
    readset_file,
    metadata,
    ncovtools_directory,
    ncovtools_config,
    output_dir,
    ini_section='prepare_report'):

    return Job(
                input_files=[output_filtered_bam, output_primer_trimmed_bam, output_consensus, output_variants],
                output_files=[],
                module_entries=[
                    [ini_section, 'module_ncovtools']
                ],
                command="""\\
module purge && \\
module load {ncovtools} && \\
echo "Preparing to run ncov_tools..." && \\
ln -sf {reference_genome} {ncovtools_directory}/{reference_genome_file} && \\
NEG_CTRL=$(grep -Ei "((negctrl|ext)|ntc)|ctrl_neg|neg" {readset_file} | awk '{{pwet=pwet", ""\\""$1"\\""}} END {{print substr(pwet,2)}}') && \\
echo "data_root: data
platform: \\"{platform}\\"
run_name: \\"{run_name}\\"
reference_genome: {reference_genome_file}
amplicon_bed: {amplicon_bed}
primer_bed: {primer_bed}
offset: 0
completeness_threshold: 0.9
bam_pattern: \\"{{data_root}}/{{sample}}{bam_pattern_extension}\\"
primer_trimmed_bam_pattern: \\"{{data_root}}/{{sample}}{primer_trimmed_bam_pattern_extension}\\"
consensus_pattern: \\"{{data_root}}/{{sample}}{consensus_pattern_extension}\\"
variants_pattern: \\"{{data_root}}/{{sample}}{variants_pattern_extension}\\"
metadata: \\"{metadata}\\"
negative_control_samples: [$NEG_CTRL]
primer_prefix: \\"{primer_prefix}\\" 
assign_lineages: true" > {ncovtools_config} && \\
echo "Running ncov_tools..." && \\
cd {ncovtools_directory} && \\
snakemake --unlock --configfile {ncovtools_config_local} --cores {nb_threads} -s $NCOVTOOLS_SNAKEFILE
snakemake --rerun-incomplete --configfile {ncovtools_config_local} --cores {nb_threads} -s $NCOVTOOLS_SNAKEFILE all
snakemake --rerun-incomplete --configfile {ncovtools_config_local} --cores {nb_threads} -s $NCOVTOOLS_SNAKEFILE all_qc_summary
snakemake --rerun-incomplete --configfile {ncovtools_config_local} --cores {nb_threads} -s $NCOVTOOLS_SNAKEFILE all_qc_analysis
cd {output_dir} && \\
module purge""".format(
                    ncovtools=global_config_parser.param(ini_section, 'module_ncovtools'),
                    reference_genome=global_config_parser.param(ini_section, 'reference_genome', required=True),
                    ncovtools_directory=ncovtools_directory,
                    reference_genome_file=global_config_parser.param(ini_section, 'reference_genome', required=True).split("/")[-1],
                    readset_file=readset_file,
                    platform=global_config_parser.param(ini_section, 'platform', required=True),
                    run_name=global_config_parser.param(ini_section, 'run_name', required=True),
                    amplicon_bed=global_config_parser.param(ini_section, 'amplicon_bed', required=True),
                    primer_bed=global_config_parser.param(ini_section, 'primer_bed', required=True),
                    bam_pattern_extension=re.sub(r"^.*?\.", ".", output_filtered_bam),
                    primer_trimmed_bam_pattern_extension=re.sub(r"^.*?\.", ".", output_primer_trimmed_bam),
                    consensus_pattern_extension=re.sub(r"^.*?\.", ".", output_consensus),
                    variants_pattern_extension=re.sub(r"^.*?\.", ".", output_variants),
                    metadata=os.path.basename(metadata),
                    primer_prefix=global_config_parser.param(ini_section, 'primer_prefix'),
                    ncovtools_config=ncovtools_config,
                    ncovtools_config_local=os.path.basename(ncovtools_config),
                    nb_threads=global_config_parser.param(ini_section, 'nb_threads'),
                    output_dir=output_dir)
                )
