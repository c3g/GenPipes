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

# Python Standard Modules
import logging
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job 

def bamqc(input_bam, output_directory, output, options, ini_section='qualimap'):

    inputs = [input_bam]
    raw_data_qualimap_dir = os.path.join(os.path.dirname(output), "raw_data_qualimapReport")
    outputs = [
        output,
        os.path.join(raw_data_qualimap_dir, "coverage_histogram.txt"),
        os.path.join(raw_data_qualimap_dir, "insert_size_histogram.txt"),
        os.path.join(raw_data_qualimap_dir, "genome_fraction_coverage.txt"),
        os.path.join(raw_data_qualimap_dir, "mapped_reads_gc-content_distribution.txt")
    ]

    return Job(
        inputs,
        outputs,
        [
            [ini_section, 'module_java'],
            [ini_section, 'module_qualimap'],
        ],
        command="""\
qualimap bamqc {other_options} \\
  -bam {input_bam} -outdir {output_directory} \\
  --java-mem-size={ram}""".format(
            input_bam=input_bam,
            output_directory=output_directory,
            other_options=options,
#            bed="\\\n  --feature-file " + bed if bed else "",
            ram=global_conf.global_get(ini_section, 'ram'),
        ),
        removable_files=[],
        report_files=outputs
    )

def rnaseq(input_bam, output_directory, output):

    inputs = [input_bam]
    outputs = [output]

    return Job(
        inputs,
        outputs,
        [
            ['DEFAULT', 'module_java'],
            ['qualimap', 'module_qualimap'],
        ],
        command="""\
qualimap rnaseq \\
  -bam {input_bam} \\
  -gtf {gtf} \\
  -outdir {output_directory} \\
  -oc {output} \\
  --java-mem-size={ram} \\
  {other_options}""".format(
            input_bam=input_bam,
            gtf=global_conf.global_get('qualimap', 'gtf', param_type='filepath'),
            output_directory=output_directory,
            output=output,
            ram=global_conf.global_get('qualimap', 'ram'),
            other_options=global_conf.global_get('qualimap', 'other_options')
        ),
        removable_files=[]
    )

def multibamqc(inputs, output_directory):

    outputs = [os.path.join(output_directory, "report.html")]

    job = Job(
        inputs,
        outputs,
        [
            ['qualimap_multibamqc', 'module_qualimap'],
            ['qualimap_multibamqc', 'module_R']
        ],
        command="""\
qualimap multi-bamqc \\
  -d {output_directory}/multi-bamqc_list.txt \\
  -outdir {output_directory} \\
  -outfile {outfile}""".format(
            output_directory=output_directory,
            outfile=os.path.join(output_directory, "report.html")
        ),
        removable_files=[]
    )

    job = concat_jobs([
        Job(command="""\
for i in {input_files}; do \\
  path1=$(dirname $i); \\
  path2=$(dirname $path1); \\
  echo -e \"$(basename $path2)\t$path1\"; \\
done > {output_directory}/multi-bamqc_list.txt""".format(
                input_files=" ".join(inputs),
                output_directory=output_directory
            )),
        job
    ])

    return job

def parse_median_insert_size_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export median_insert_size=`awk '{{if ($0 ~ /median insert size/) {{match($0,/[0-9]+$/,value); print value[0]}}}}' {input_file}`"""
        )

def parse_mean_insert_size_metrics_pt(input_file):
    """
    Calculus coming from https://stackoverflow.com/questions/46086663/how-to-get-mean-and-standard-deviation-from-a-frequency-distribution-table
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
module load mugqic/python/3.10.4 &&
export mean_insert_size=`python <<EOF
import numpy as np
dataset = np.genfromtxt(fname="{input_file}", delimiter="\t", skip_header=1)
frequencies = dataset[:, 1]
values = dataset[:, 0]
print(np.around(np.average(values, weights=frequencies), decimals=1))
EOF`"""
        )

def parse_dedup_coverage_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export dedup_coverage=`awk '{{if ($0 ~ /mean coverage/) {{match($0,/[0-9]+.[0-9]+/,value); printf "%.2f", value[0]}}}}' {input_file}`"""
        )

def parse_aligned_reads_count_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export aligned_reads_count=`awk '{{if ($0 ~ /number of mapped reads/) {{match($0,/[0-9]+.* /,value); print value[0]}}}}' {input_file} | head -n 1 | sed -e 's@,@@g' | sed -e 's@ @@g'`"""
        )
