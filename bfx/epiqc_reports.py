#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os
# MUGQIC Modules

from core.job import *
from core.config import *

def bigwiginfo_report(bigwiginfo_file, output_file):
    return Job(
        [bigwiginfo_file],
        [output_file],
        [['epiqc_report', 'module_python']],
        command="""\
python /home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc/epiqc_2021/genpipes/bfx/epiqc_report.py \\
-b {bigwiginfo_file} \\
-cb {chromCount} \\
-bc1 {low_alert_bases_covered} \\
-bc2 {medium_alert_bases_covered} \\
-o {output_file}""".format(
            bigwiginfo_file=bigwiginfo_file,
            chromCount=config.param('epiqc_report', 'chromcount_threshold'),
            low_alert_bases_covered=config.param('epiqc_report', 'low_alert_bases_covered'),
            medium_alert_bases_covered=config.param('epiqc_report', 'medium_alert_bases_covered'),
            output_file=output_file)
    )
def signal_to_noise_report(signal_noise_file, report_file):
    return Job(
                        [signal_noise_file],
                        [report_file],
                        [['epiqc_report', 'module_python']],
                        command="""\
python /home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc/epiqc_2021/genpipes/bfx/epiqc_report.py \\
-s {signal_noise_file} \\
-s1 {percent1} \\
-s2 {percent2} \\
-st1 {signal_noise_threshold_M} \\
-st2 {signal_noise_threshold_L} \\
-o {output_file}""".format(
                            signal_noise_file=signal_noise_file,
                            percent1=config.param('signal_noise', 'percent1'),
                            percent2=config.param('signal_noise', 'percent2'),
                            signal_noise_threshold_M=config.param('epiqc_report', 'signal_noise_threshold_M'),
                            signal_noise_threshold_L=config.param('epiqc_report', 'signal_noise_threshold_L'),
                            output_file=report_file))

def chromimpute_report(eval_file, report_file):
    return Job(
                        [eval_file],
                        [report_file],
                        [['epiqc_report', 'module_python']],
                        command="""\
python /home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc/epiqc_2021/genpipes/bfx/epiqc_report.py \\
-c {eval_file} \\
-p1 {percent1} \\
-p2 {percent2} \\
-ct1 {chromimpute_threshold_M} \\
-ct2 {chromimpute_threshold_L} \\
-o {output_file}""".format(
                            eval_file=eval_file,
                            percent1=config.param('chromimpute_eval', 'percent1'),
                            percent2=config.param('chromimpute_eval', 'percent2'),
                            chromimpute_threshold_M=config.param('epiqc_report', 'chromimpute_threshold_M'),
                            chromimpute_threshold_L=config.param('epiqc_report', 'chromimpute_threshold_L'),
                            output_file=report_file))


def epigeec_report(input_matrix_file, output_heatmap_file, report_dir):
    return Job(
            [input_matrix_file],
            [output_heatmap_file],
            [['epiqc_report', 'module_python']],
            name="epigeec_report",
            command="""\
python /home/pubudu/projects/rrg-bourqueg-ad/pubudu/epiqc/epiqc_2021/genpipes/bfx/epiqc_report.py \\
-e {correlation_matrix} \\
-o {output_dir}""".format(
                correlation_matrix=input_matrix_file,
                output_dir=report_dir))