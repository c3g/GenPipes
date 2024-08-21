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

import os

from ..core.config import global_conf
from ..core.job import Job

#Assuming that the module file has already been written and the scripts have been added to path

def genome_gz(output, ini_section='sequenza'):
    return Job(
        [None],
        [output],
        [
            [ini_section, 'module_sequenza_utils'],
            [ini_section, 'module_R'],
        ],
        command="""\\
sequenza-utils \\
    GC-windows -w {window} \\
    {input} \\
    | gzip > \\
    {out}""".format(
        input=global_conf.global_get('samtools_mpileup', 'genome_fasta', param_type='filepath'),
        out=output,
        window=global_conf.global_get(ini_section, 'window_length')
    ),
  )


def bam2seqz_mpileup(normal_gz, tumor_gz, genome, output, ini_section='sequenza'):
    return Job(
        [normal_gz, tumor_gz],
        [output],
        [
            [ini_section, 'module_python'],
            [ini_section, 'module_R'],
        ],
        command="""\\
sequenza-utils \\
    bam2seqz -p {pileup_options} \\
    -gc  {gen}   \\
    -n {normal}  \\
    -t {tumor}{out}""".format(
        gen=genome,
        normal=normal_gz,
        tumor=tumor_gz,
        pileup_options=global_conf.global_get(ini_section, 'pileup_options'),
        out=" \\\n > " + output if output else ""
        )
    )

def bam2seqz(normal, tumor, genome, output, chr=None, ini_section='sequenza'):
    return Job(
        [normal, tumor],
        [output],
        [
            [ini_section, 'module_sequenza_utils'],
            [ini_section, 'module_samtools'],
            [ini_section, 'module_htslib'],
        ],
        command="""\\
sequenza-utils \\
    bam2seqz {pileup_options} --samtools samtools --tabix tabix \\
    {chr} \\
    -gc {gen} \\
    --fasta {reference_sequence} \\
    --normal {normal} \\
    --tumor {tumor} \\
    --output {out}""".format(
        chr="\\\n    --chromosome " + chr if chr else "",
        gen=genome,
        reference_sequence=global_conf.global_get(ini_section, 'genome_fasta', param_type='filepath'),
        normal=normal,
        tumor=tumor,
        pileup_options=global_conf.global_get(ini_section, 'pileup_options'),
        out=output
        )
    )

def bin(seqz_gz, output, ini_section='sequenza'):
    return Job(
        [seqz_gz],
        [output],
        [
            [ini_section, 'module_sequenza_utils'],
            [ini_section, 'module_R'],
        ],
        command="""\\
sequenza-utils  \\
    seqz_binning  \\
    -w {window}  \\
    -s {seqz_gz} \\
    -o {output}""".format(
        window=global_conf.global_get(ini_section, 'bin_window_size'),
        seqz_gz=seqz_gz,
        output=output,
        )
    )

def main(
        seqz,
        output_folder,
        sample_name,
        ini_section='sequenza'
):
    output_dep = [os.path.join(output_folder, sample_name + "_chromosome_view.pdf"),
                  os.path.join(output_folder, sample_name + "_genome_view.pdf"),
                  os.path.join(output_folder, sample_name + "_CN_bars.pdf"),
                  os.path.join(output_folder, sample_name + "_CP_contours.pdf"),
                  os.path.join(output_folder, sample_name + "_segments.txt"),
                  os.path.join(output_folder, sample_name + "_ploidy_celularity.tsv")]
    return Job(
        [seqz],
        output_dep,
        [
            [ini_section, 'module_mugqic_tools'],
            [ini_section, 'module_R'],
        ],
        command="""\\
Rscript $R_TOOLS/RunSequenza_analysis.R \\
    {input}   \\
    {OUTPUT_FOLDER}   \\
    {OUTPUT_BASE_NAME}""".format(
            input=seqz,
            OUTPUT_FOLDER=output_folder,
            OUTPUT_BASE_NAME=sample_name
        )
    )

def filter(calls, pair_name, output, ini_section='sequenza'):
     return Job(
         [calls],
         [output],
         [
             [ini_section, 'module_mugqic_tools']
         ],
         command="""\
sequenza_filterOut.sh \\
    {scones_calls} \\
    {output} \\
    {pair_name} """.format(
             scones_calls=calls,
             output=output,
             pair_name=pair_name
         )
     )

def annotate(calls_filtered, output_basename, tmp_basename, ini_section='scones_annotate'):
     scones_outputs = [output_basename + ".counts.filteredSV.annotate.txt",
                       output_basename + ".other.filteredSV.annotate.txt",
                       output_basename + ".TumS.filteredSV.annotate.txt"]
    
     return Job(
         [calls_filtered],
         scones_outputs,
         [
             ['sequenza', 'module_mugqic_tools']
         ],
         command="""\
sequenza_filterAnnotCNV.sh \\
    {scones_calls_filtered} \\
    {excluded_regions} \\
    {genes} \\
    {DGV} \\
    {microsat} \\
    {repeatMasker} \\
    {AutosomeSize} \\
    {output_basename} \\
    {tmp_basename} """.format(
             scones_calls_filtered=calls_filtered,
             excluded_regions=global_conf.global_get(ini_section, 'excluded_regions_bed', param_type='filepath', required=True),
             genes=global_conf.global_get(ini_section, 'genes_bed', param_type='filepath', required=True),
             DGV=global_conf.global_get(ini_section, 'dgv_bed', param_type='filepath', required=True),
             microsat=global_conf.global_get(ini_section, 'microsat_bed', param_type='filepath', required=True),
             repeatMasker=global_conf.global_get(ini_section, 'repeat_masker_bed', param_type='filepath', required=True),
             AutosomeSize=global_conf.global_get(ini_section, 'autosome_size_file', param_type='filepath', required=True),
             output_basename=output_basename,
             tmp_basename=tmp_basename
         )
     )
