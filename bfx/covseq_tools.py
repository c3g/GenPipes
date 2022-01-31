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


def covid_collect_metrics(readset_file, covid_collect_metrics_inputs, ini_section='prepare_table'):

    return Job(
                input_files=covid_collect_metrics_inputs,
                output_files=[os.path.join("metrics", "metrics.csv"), os.path.join("metrics", "host_contamination_metrics.tsv"), os.path.join("metrics", "host_removed_metrics.tsv"), os.path.join("metrics", "kraken2_metrics.tsv")],
                module_entries=[
                    [ini_section, 'module_R'],
                    [ini_section, 'module_CoVSeQ_tools']
                ],
                command="""\\
echo "Collecting metrics..." && \\
covid_collect_metrics.sh -r {readset_file}""".format(
    readset_file=readset_file
    )
                )

def covid_collect_nanopore_metrics(sample_name, consensus, fq_stats, pickle, output, ini_section='covseq_metrics'):

    return Job(
                input_files=[consensus, fq_stats, pickle],
                output_files=[output],
                module_entries=[
                    [ini_section, 'module_python'],
                    [ini_section, 'module_CoVSeQ_tools']
                ],
                command="""\\
covid_collect_nanopore_metrics.py -s {sample_name} -c {consensus} -fqs {fq_stats} -pk {pickle} -o {output}""".format(
    sample_name=sample_name,
    consensus=consensus,
    fq_stats=fq_stats,
    pickle=pickle,
    output=output
    )
                )

def generate_report_tables(readset_file_report, output_name_pattern, ini_section='prepare_report'):
    metrics=os.path.join("metrics", "metrics.csv")
    host_contamination_metrics=os.path.join("metrics", "host_contamination_metrics.tsv")

    return Job(
                input_files=[metrics, host_contamination_metrics],
                output_files=[output_name_pattern + ".csv", output_name_pattern + ".tsv"],
                module_entries=[
                    [ini_section, 'module_R'],
                    [ini_section, 'module_CoVSeQ_tools']
                ],
                command="""\\
echo "Generating report tables..." && \\
generate_report_tables.R --report_readset={readset_file_report} --metrics={metrics} --host_contamination_metrics={host_contamination_metrics} --output_name_pattern={output_name_pattern}""".format(
    readset_file_report=readset_file_report,
    metrics=metrics,
    host_contamination_metrics=host_contamination_metrics,
    output_name_pattern=output_name_pattern
    )
                )

def render_report(software_version, run_metadata, output_name_pattern, caller, ini_section='prepare_report'):
    if caller == "ivar":
        report_template="$RUN_REPORT"
    elif caller == "freebayes":
        report_template="$RUN_REPORT_FREEBAYES"
    else:
        report_template=""

    return Job(
                input_files=[software_version, run_metadata],
                output_files=[output_name_pattern + ".pdf"],
                module_entries=[
                    [ini_section, 'module_R'],
                    [ini_section, 'module_CoVSeQ_tools']
                ],
                command="""\\
echo "Rendering report..." && \\
Rscript -e "report_path <- tempfile(fileext = '.Rmd'); file.copy('{report_template}', report_path, overwrite = TRUE); rmarkdown::render(report_path, output_file='{output_name_pattern}.pdf', output_format = 'all', output_dir='$(pwd)/report', knit_root_dir='$(pwd)')" """.format(
    report_template=report_template,
    output_name_pattern=output_name_pattern
    )
                )
