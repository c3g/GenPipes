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

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def pileup(input_bam, output):
    """
    """
    return Job(
        [input_bam],
        [output],
        [
            ['conpair_concordance_contamination', 'module_java'],
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_gatk'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
python3 $CONPAIR_SCRIPTS/run_gatk_pileup_for_sample.py -t {tmp_dir} \\
  -m {ram} \\
  -G $GATK_JAR \\
  -D $CONPAIR_DIR \\
  -R {reference_sequence} \\
  -M {markers} \\
  -B {input} \\
  -O {output} {other_options}""".format(
        tmp_dir=global_conf.global_get('conpair_concordance_contamination', 'tmp_dir'),
        ram=global_conf.global_get('conpair_concordance_contamination', 'ram'),
        reference_sequence=global_conf.global_get('conpair_concordance_contamination', 'genome_fasta', param_type='filepath'),
        markers=global_conf.global_get('conpair_concordance_contamination', 'markers_bed'),
        input=input_bam,
        other_options=global_conf.global_get('conpair_concordance_contamination', 'other_options', required= False),
        output=output
        )
    )

def concordance(input_normal, input_tumor, output):
    """
    """
    return Job(
        [input_normal, input_tumor],
        [output],
        [
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
rm -f {output} && \\
python $CONPAIR_SCRIPTS/verify_concordance.py {options} \\
  --markers {markers} \\
  --normal_pileup {input_normal} \\
  --tumor_pileup {input_tumor} \\
  --outfile {output}""".format(
        options=global_conf.global_get('conpair_concordance_contamination', 'concord_options', required=False),
        markers=global_conf.global_get('conpair_concordance_contamination', 'markers_txt'),
        input_normal=input_normal,
        input_tumor=input_tumor,
        output=output
        )
    )

def contamination(input_normal, input_tumor, output):
    """
    """
    return Job(
        [input_normal, input_tumor],
        [output],
        [
            ['conpair_concordance_contamination', 'module_python'],
            ['conpair_concordance_contamination', 'module_conpair']
        ],
        command="""\
rm -f {output} && \\
python3 $CONPAIR_SCRIPTS/estimate_tumor_normal_contamination.py {options} \\
  --markers {markers} \\
  --normal_pileup {input_normal} \\
  --tumor_pileup {input_tumor} \\
  --outfile {output}""".format(
        options=global_conf.global_get('conpair_concordance_contamination', 'contam_options'),
        markers=global_conf.global_get('conpair_concordance_contamination', 'markers_txt'),
        input_normal=input_normal,
        input_tumor=input_tumor,
        output=output
        )
    )

def parse_concordance_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export concordance=`awk '{{if ($0 ~ /^Concordance:/) {{match($0,/[0-9]+.[0-9]+/,value); print value[0]}}}}' {input_file}`"""
        )

def parse_contamination_normal_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export contamination=`awk '{{if ($0 ~ /^Normal/) {{match($0,/[0-9]+.[0-9]+/,value); print value[0]}}}}' {input_file}`"""
    )

def parse_contamination_tumor_metrics_pt(input_file):
    """
    """
    return Job(
        [input_file],
        [],
        [],
        command=f"""\
export contamination=`awk '{{if ($0 ~ /^Tumor/) {{match($0,/[0-9]+.[0-9]+/,value); print value[0]}}}}' {input_file}`"""
    )
