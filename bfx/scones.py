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

from core.config import *
from core.job import *

def scones_pair(bined_file, output_basename, window):
    
    models= [0,2,3,4,5,6,7]
    scones_outputs=[output_basename + "_Model_" + str(i) + "_GenomicRatios.pdf" for i in models ] + [output_basename + "_Model_" + str(i) + "_CNVcalls.txt" for i in models ]
    
    return Job(
        [bined_file],
        scones_outputs,
        [
            ['scones_pair', 'module_R'],
            ['scones_pair', 'module_scones']
        ],
        command="""\
        Rscript $SCONES  {options} \\
        -f {bined_file} \\
        -o {output_basename} \\
        -c {GC_map_bed} \\
        -b {window} """.format(
            options=global_config_parser.param('scones_pair', 'other_options') if global_config_parser.param('scones_pair', 'other_options') else "",
            bined_file=bined_file,
            output_basename=output_basename,
            GC_map_bed=global_config_parser.param('scones_pair', 'gc_map_bedfile', param_type='filepath', required=True),
            window=window
        )
    )


def scones_filter(scones_calls,pair_name, output):
    
    return Job(
        [scones_calls],
        [output],
        [
            ['scones_filter', 'module_scones']
        ],
        command="""\
        filterOut.sh \\
        {scones_calls} \\
        {output} \\
        {pair_name} """.format(
            scones_calls=scones_calls,
            output=output,
            pair_name=pair_name
        )
    )

def scones_annotate(scones_calls_filtered, output_basename, tmp_basename):
    
    scones_outputs=[output_basename + ".counts.filteredSV.annotate.txt", output_basename + ".other.filteredSV.annotate.txt", output_basename + ".TumS.filteredSV.annotate.txt"]
    
    return Job(
        [scones_calls_filtered],
        scones_outputs,
        [
            ['scones_annotate', 'module_scones']
        ],
        command="""\
        filterAnnotCNV.sh \\
        {scones_calls_filtered} \\
        {excluded_regions} \\
        {genes} \\
        {DGV} \\
        {microsat} \\
        {repeatMasker} \\
        {AutosomeSize} \\
        {output_basename} \\
        {tmp_basename} """.format(
            scones_calls_filtered=scones_calls_filtered,
            excluded_regions=global_config_parser.param('scones_annotate', 'excluded_regions_bed', param_type='filepath', required=True),
            genes=global_config_parser.param('scones_annotate', 'genes_bed', param_type='filepath', required=True),
            DGV=global_config_parser.param('scones_annotate', 'dgv_bed', param_type='filepath', required=True),
            microsat=global_config_parser.param('scones_annotate', 'microsat_bed', param_type='filepath', required=True),
            repeatMasker=global_config_parser.param('scones_annotate', 'repeat_masker_bed', param_type='filepath', required=True),
            AutosomeSize=global_config_parser.param('scones_annotate', 'autosome_size_file', param_type='filepath', required=True),
            output_basename=output_basename,
            tmp_basename=tmp_basename
        )
    )
