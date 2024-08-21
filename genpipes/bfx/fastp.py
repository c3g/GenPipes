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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def basic_qc(
    input1,
    input2,
    output_json_path,
    output_html_path=None,
    overrepresentation_analysis=True,
    ini_section='fastp'
    ):

    if input2:
        inputs = [input1, input2]
    else:
        inputs = [input1]

    num_threads = global_conf.global_get(ini_section, 'threads', required=False, param_type='posint')
    output_files = filter(None, [output_json_path, output_html_path])

    return Job(
        inputs,
        output_files,
        module_entries=[
            [ini_section, 'module_fastp']
        ],
        command="""\
fastp -V \\
  {rds1} \\
  {rds2} \\
  {cpus} \\
  {json} \\
  {html}""".format(
            rds1 = "\\\n  --in1 "    + input1,
            rds2 = "\\\n  --in2 "    + input2 if input2 else "",
            orep = "\\\n  --overrepresentation_analysis" if overrepresentation_analysis else "",
            cpus = "\\\n  --thread " + str(num_threads) if num_threads else "",
            json = "\\\n  --json "   + output_json_path if output_json_path else "",
            html = "\\\n  --html "   + output_html_path if output_html_path else "",
        )
    )

def trim(
        adapter,
        input1,
        input2,
        output1,
        output2,
        output_json_path,
        output_html_path=None,
        ini_section = None
):

    if input2:
        inputs = [input1, input2]
    else:
        inputs = [input1]
        
    if output2:
        outputs = [output1, output2, output_json_path, output_html_path]
    else:
        outputs = [output1, output_json_path, output_html_path]

    num_threads = global_conf.global_get(ini_section, 'threads', required=False, param_type='posint')
    parameters = global_conf.global_get(ini_section, 'options', required=False)

    return Job(
        inputs,
        outputs,
        module_entries=[
            [ini_section, 'module_fastp']
        ],
        command="""\
fastp -V \\
  {options} \\
  {adapter} \\
  {rds1} \\
  {rds2} \\
  {out1} \\
  {out2} \\
  {cpus} \\
  {json} \\
  {html}""".format(
            options = "  "                  + parameters,
            adapter = "  --adapter_fasta "  + adapter,
            rds1 = "  --in1 "               + input1,
            rds2 = "  --in2 "               + input2 if input2 else "",
            out1 = "  --out1 "              + output1,
            out2 = "  --out2 "              + output2 if output2 else "",
            cpus = "  --thread "            + str(num_threads) if num_threads else "",
            json = "  --json "              + output_json_path if output_json_path else "",
            html = "  --html "              + output_html_path if output_html_path else "",
        )
    )

def parse_quality_thirty_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export bases_over_q30_percent=`jq -r '.summary.after_filtering.q30_rate' {input_file}`"""
        )

def parse_pre_length_r1_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export pre_mean_length_r1=`jq -r '.summary.before_filtering.read1_mean_length' {input_file}`"""
        )

def parse_post_length_r1_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export post_mean_length_r1=`jq -r '.summary.after_filtering.read1_mean_length' {input_file}`"""
        )

def parse_pre_length_r2_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export pre_mean_length_r2=`jq -r '.summary.before_filtering.read2_mean_length' {input_file}`"""
        )

def parse_post_length_r2_metrics(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export post_mean_length_r2=`jq -r '.summary.after_filtering.read2_mean_length' {input_file}`"""
        )