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
from core.config import *
from core.job import *

def bedGraphToBigWig(input_bed_graph, output_wiggle, header=True, ini_section='ucsc'):

    # If using GRCh37 assembly then create temporary dict file
    # only using autosomal + X + Y chromosomes and adding chr (e.g. chr1, instead of 1)
    if (global_config_parser.param('DEFAULT', 'assembly') == 'GRCh37') :
        chromosome_size_file = re.sub(".fa.fai", ".withchr.fa.fai",
                                      global_config_parser.param(ini_section,
                                                                 'chromosome_size', param_type='filepath'))
        chromosome_prefix = "chr"
        chromosome_sed = "| sed 's/MT/chrM/'"
    else :
        chromosome_size_file = global_config_parser.param(ini_section, 'chromosome_size', param_type='filepath')
        chromosome_prefix = ""
        chromosome_sed = ""

    # Check it the input is a real bedGrah (i.e. contains the bedGraph header : track type=bedGraph)
    # or if it is just a regular bed file (i.e. no bedGraph header)
    if header :
        remove_head_command="""\
({open} {input_bed_graph} | head -n 1 > {input_bed_graph}.head.tmp ; ec=$?; if [ "$ec" -eq 141 ]; then exit 0; else exit "$ec"; fi) && \\
{open} {input_bed_graph} | awk ' NR > 1 ' | sort  --temporary-directory={temp_dir} -k1,1 -k2,2n | \\
awk '{{if($0 !~ /^[A-W]/) print "{chrom}"$0; else print $0}}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random" {chromosome_sed} | \\
awk '{{printf "%s\\t%d\\t%d\\t%4.4g\\n", $1,$2,$3,$4}}' > {input_bed_graph}.body.tmp && \\
cat {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp > {input_bed_graph}.sorted && \\
rm {input_bed_graph}.head.tmp {input_bed_graph}.body.tmp""".format(
            open="zcat" if (os.path.splitext(input_bed_graph)[1] == ".gz") else "cat",
            input_bed_graph=input_bed_graph,
            temp_dir=global_config_parser.param(ini_section, 'tmp_dir', required=True),
            chrom=chromosome_prefix,
            chromosome_sed=chromosome_sed
        )
    else:
        remove_head_command="""\
{open} {input_bed_graph} | sort --temporary-directory={temp_dir} -k1,1 -k2,2n | \\
awk '{{if($0 !~ /^[A-W]/) print "{chrom}"$0; else print $0}}' | grep -v "GL\|lambda\|pUC19\|KI\|\KN\|random" {chromosome_sed} | \\
awk '{{printf "%s\\t%d\\t%d\\t%4.4g\\n", $1,$2,$3,$4}}' > {input_bed_graph}.sorted""".format(
            open="zcat" if (os.path.splitext(input_bed_graph)[1] == ".gz") else "cat",
            input_bed_graph=input_bed_graph,
            temp_dir=global_config_parser.param(ini_section, 'tmp_dir', required=True),
            chrom=chromosome_prefix,
            chromosome_sed=chromosome_sed
        )

    return Job(
        [input_bed_graph],
        [input_bed_graph+".sorted", output_wiggle],
        [
            [ini_section, 'module_ucsc']
        ],
        command="""\
{remove_head_command} && \\
bedGraphToBigWig \\
  {input_bed_graph}.sorted \\
  {chromosome_size} \\
  {output_wiggle}""".format(
            remove_head_command=remove_head_command,
            chromosome_size=chromosome_size_file,
            input_bed_graph=input_bed_graph,
            output_wiggle=output_wiggle
        ),
        removable_files=[input_bed_graph+".sorted"]
    )

def bedToBigBed(bed_file, bigBed_file, ini_section='ucsc'):

    return Job(
        [bed_file],
        [bigBed_file],
        [
            [ini_section, 'module_ucsc']
        ],
        command="""\
bedToBigBed \\
  {bed_file} \\
  {chromosome_size} \\
  {bigBed_file}""".format(
            chromosome_size=global_config_parser.param(ini_section, 'chromosome_size', param_type='filepath'),
            bed_file=bed_file,
            bigBed_file=bigBed_file
        )
    )
