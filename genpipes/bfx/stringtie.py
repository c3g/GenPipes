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
import os

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

def stringtie(input_bam, output_directory, gtf=None, abund=False):
    ## Get param from config file
    stranded = global_conf.global_get('stringtie', 'strand_info')

    ## Check library type for strand information
    if stranded.lower() == "fr-firststrand":
        strd_cmd = "--rf"
    elif stranded.lower() == "fr-secondstrand":
        strd_cmd = "--fr"
    elif stranded.lower == "unstranded":
        strd_cmd = ""
    else:
        raise Exception("Strand info\"" + stranded + "\" unrecognized")

    ## Define output files depending on whether or not abundances will be calculated
    if abund:
        out_files = [os.path.join(output_directory , "transcripts.gtf"),
        os.path.join(output_directory , "abundance.tab")]
    else:
        out_files = [os.path.join(output_directory ,"transcripts.gtf")]

    return Job(
        [input_bam, gtf],
        out_files,
        [["stringtie", "module_stringtie"]],
        command="""\
mkdir -p {output_directory} && \\
stringtie -v {other_options} {strd_cmd} {gtf} {abund_cmd} \\
  -p {num_threads} \\
  -m {min_length} \\
  -o {outgtf} \\
  {input_bam}""".format(
      output_directory=output_directory,
      other_options=global_conf.global_get('stringtie', 'other_options', required=False),
      strd_cmd="\\\n  " + strd_cmd if strd_cmd else " ",
      gtf="\\\n  -G " + gtf if gtf else " ",
      abund_cmd="\\\n  -eB -A " + os.path.join(output_directory, "abundance.tab") if abund else " ",
      num_threads=global_conf.global_get('stringtie', 'threads', param_type='posint'),
      min_length=global_conf.global_get('stringtie', 'min_length', param_type='posint'),
      outgtf=os.path.join(output_directory, "transcripts.gtf"),
      input_bam=input_bam
        )
    )

def stringtie_merge(gtf_list, output_prefix, gtf=None):

    return Job(
        [gtf_list],
        [os.path.join(output_prefix, "merged.gtf")],
        [["stringtie", "module_stringtie"]],
        command="""\
mkdir -p {output_directory} &&\\
stringtie --merge {other_options}{gtf} \\
  -o {outfile} \\
  {gtf_list}""".format(
      output_directory=os.path.dirname(output_prefix),
      other_options=global_conf.global_get('stringtie_merge', 'other_options', required=False),
      gtf="\\\n -G " + gtf if gtf else "",
      outfile=os.path.join(output_prefix, "merged.gtf"),
      gtf_list=gtf_list
        )
    )
