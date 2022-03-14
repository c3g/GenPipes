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
import os
# MUGQIC Modules

from core.job import *
from core.config import *

def bigwiginfo_report(bigwiginfo_file, output_file):
    return Job(
        [bigwiginfo_file],
        [output_file],
        [['epiqc_report', 'module_python'],
         ['epiqc_report', 'module_mugqic_tools']],
        command="""\
python $PYTHON_TOOLS/epiqc_report.py \\
-b {bigwiginfo_file} \\
-cb {chromCount} \\
-bc1 {low_alert_bases_covered} \\
-bc2 {medium_alert_bases_covered} \\
-o {output_file}""".format(
            bigwiginfo_file=bigwiginfo_file,
            chromCount=global_config_parser.param('epiqc_report', 'chromcount_threshold'),
            low_alert_bases_covered=global_config_parser.param('epiqc_report', 'low_alert_bases_covered'),
            medium_alert_bases_covered=global_config_parser.param('epiqc_report', 'medium_alert_bases_covered'),
            output_file=output_file)
    )
def signal_to_noise_report(signal_noise_file, report_file, chr_name):
    return Job(
                        [signal_noise_file],
                        [report_file],
                        [['epiqc_report', 'module_python'],
                         ['epiqc_report', 'module_mugqic_tools']],
                        command="""\
python $PYTHON_TOOLS/epiqc_report.py \\
-s {signal_noise_file} \\
-s1 {percent1} \\
-s2 {percent2} \\
-chr {chromosome} \\
-st1 {signal_noise_threshold_M} \\
-st2 {signal_noise_threshold_L} \\
-o {output_file}""".format(
                            signal_noise_file=signal_noise_file,
                            percent1=global_config_parser.param('signal_noise', 'percent1'),
                            percent2=global_config_parser.param('signal_noise', 'percent2'),
                            chromosome =chr_name,
                            signal_noise_threshold_M=global_config_parser.param('epiqc_report', 'signal_noise_threshold_M'),
                            signal_noise_threshold_L=global_config_parser.param('epiqc_report', 'signal_noise_threshold_L'),
                            output_file=report_file))

def chromimpute_report(eval_file, report_file):
    return Job(
                        [eval_file],
                        [report_file],
                        [['epiqc_report', 'module_python'],
                         ['epiqc_report', 'module_mugqic_tools']],
                        command="""\
python $PYTHON_TOOLS/epiqc_report.py \\
-c {eval_file} \\
-p1 {percent1} \\
-p2 {percent2} \\
-ct1 {chromimpute_threshold_M} \\
-ct2 {chromimpute_threshold_L} \\
-o {output_file}""".format(
                            eval_file=eval_file,
                            percent1=global_config_parser.param('chromimpute_eval', 'percent1'),
                            percent2=global_config_parser.param('chromimpute_eval', 'percent2'),
                            chromimpute_threshold_M=global_config_parser.param('epiqc_report', 'chromimpute_threshold_M'),
                            chromimpute_threshold_L=global_config_parser.param('epiqc_report', 'chromimpute_threshold_L'),
                            output_file=report_file))


def epigeec_report(input_matrix_file, output_heatmap_file, report_dir):

    return Job(
            [input_matrix_file],
            [output_heatmap_file],
            [['epiqc_report', 'module_python'],
             ['epiqc_report', 'module_mugqic_tools']],
            name="epigeec_report",
            command="""\
python $PYTHON_TOOLS/epiqc_report.py \\
-e {correlation_matrix} \\
-o {output_dir}""".format(
                correlation_matrix=input_matrix_file,
                output_dir=report_dir))

def final_report(bigwiginfo_report, chromimpute_report, signalnoise_report, final_report, sample, histone):

    #generate the final report TSV file by integrating all report files from each step.
    #Do not change anything in report files, or report generating scripts since the final report is
    #completely depend on the string matches in individual reports

    #first extract first and second columns of the report and check whether they match with predefined strings
    #and then generate the report accordingly
    return Job(
            [bigwiginfo_report, chromimpute_report, signalnoise_report],
            [final_report],
            [['epiqc_report', 'module_python']],
            name="epiqc_report",
            command="""\
bigwig=$(cut -d ' ' -f1 {bigwiginfo_report} | head -n 1)
bigwig_title=$(cut -d ' ' -f2 {bigwiginfo_report} | head -n 1)

signalnoise=$(cut -d ' ' -f1 {signalnoise_report} | head -n 1)

chromimpute=$(cut -d ' ' -f1 {chromimpute_report} | head -n 1)

if [ "$bigwig" == "BigWigInfo:" ] ; then
 if [ "$bigwig_title" == "Chromosome_count:" ] ; then
  bigwig1=$(cut -d ' ' -f3 {bigwiginfo_report} | head -n 1) &&
  echo $bigwig1
fi
if [ "$bigwig_title" == "Whole_genome_bases_covered:" ] ; then
 bigwig2=$(cut -d ' ' -f3 {bigwiginfo_report} | head -n 2) &&
 echo $bigwig2
fi
else
 bigwig1="FILE_ERROR" &&
 bigwig2="FILE_ERROR" &&
 echo "Please check whether the BigWigInfo report file is properly structured or not empty"
fi &&

if [ "$chromimpute" == "ChromImpute:" ] ; then
 impute=$(cut -d ' ' -f3 {chromimpute_report} | head -n 1) &&
 echo $impute
else
 impute="FILE_ERROR" &&
 impute="FILE_ERROR" &&
 echo "Please check whether the ChromImpute report file is properly structured or not empty"
fi &&

if [ "$signalnoise" == "Signal_to_noise:" ] ; then
 signal=$(cut -d ' ' -f3 {signalnoise_report} | head -n 1) &&
 echo $signal
else
 signal="FILE_ERROR" &&
 signal="FILE_ERROR" &&
 echo "Please check whether the Signal_to_noise report file is properly structured or not empty"
fi &&

if [ "$bigwig1" == "HIGH_LEVEL_ALERT:" ] ; then
 echo -e "{sample}\t{histone}\tHIGH_LEVEL_ALERT" >> {final_report}
else
 if [[ "$bigwig2" == "MEDIUM_LEVEL_ALERT" || "$impute" == "MEDIUM_LEVEL_ALERT" || "$signal" == "MEDIUM_LEVEL_ALERT" ]] ; then
  echo -e "{sample}\t{histone}\tMEDIUM_LEVEL_ALERT" >> {final_report}
 elif [[ "$bigwig2" == "LOW_LEVEL_ALERT" || "$impute" == "LOW_LEVEL_ALERT" || "$signal" == "LOW_LEVEL_ALERT" ]] ; then
  echo -e "{sample}\t{histone}\tLOW_LEVEL_ALERT" >> {final_report}
 elif [[ "$bigwig2" == "FILE_ERROR" || "$impute" == "FILE_ERROR" || "$signal" == "FILE_ERROR" ]] ; then
  echo -e "{sample}\t{histone}\tFILE_ERROR" >> {final_report}
 else
  echo -e "{sample}\t{histone}\tPASSED" >> {final_report}
 fi
fi""".format(
            bigwiginfo_report=bigwiginfo_report,
            chromimpute_report=chromimpute_report,
            signalnoise_report = signalnoise_report,
            final_report=final_report,
            sample=sample,
            histone=histone))