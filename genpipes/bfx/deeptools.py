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
import logging

# MUGQIC Modules
from ..core.config import global_conf
from ..core.job import Job

log = logging.getLogger(__name__)

### Start from here ###

def bamcoverage(input_bam, output_file, strand=None):
    return Job(
        input_files=[input_bam],
        output_files=[output_file],
        module_entries=[['default', 'module_deeptools']
        ],
        # need to add option for clusters,
        # each other_option individually &
        # no need to have forward and reverse commands, all in one, dependent on the ini.
        command="""\
bamCoverage --verbose \\
    --outFileFormat bigwig \\
    --numberOfProcessors {cpu} \\
    --binSize {bs} \\
    --normalizeUsing {nu} \\
    --bam {input_bam} \\
    --outFileName {output_file} {strand}""".format(
            output_file=output_file,
            input_bam=input_bam,
            cpu=global_conf.global_get('wiggle', 'cluster_cpu', required=True),
            bs=global_conf.global_get('wiggle', 'bin_size', required=True),
            nu=global_conf.global_get('wiggle', 'norm_using', required=True),
            strand="--filterRNAstrand " + strand if strand else "",
        )
  )


def multibamsummary(all_bam_files, summ_matrix, ini_section='deeptools_QC'):
    return Job(
        input_files=all_bam_files,
        output_files=[summ_matrix],
        module_entries=[['default', 'module_deeptools']
        ],
        command="""\
multiBamSummary bins --verbose \\
    --numberOfProcessors {cpu} \\
    --bamfiles {all_bam_files} \\
    --outFileName {summ_matrix}""".format(
            cpu=global_conf.global_get('deeptools_QC', 'cluster_cpu', required=True),
            all_bam_files=" ".join(all_bam_files),
            summ_matrix=summ_matrix,
        )
  )


def plotcorrelation(input_matrix, corr_plot, corr_table, ini_section='deeptools_QC'):
    return Job(
        input_files=[input_matrix],
        output_files=[corr_plot, corr_table],
        module_entries=[['default', 'module_deeptools']
        ],
        command="""\
plotCorrelation --verbose \\
    --corData {input_matrix} \\
    --corMethod spearman \\
    --whatToPlot heatmap \\
    --plotFileFormat pdf \\
    --plotFile {corr_plot} \\
    --outFileCorMatrix {corr_table}""".format(
            input_matrix=input_matrix,
            corr_plot=corr_plot,
            corr_table=corr_table,
        )
  )


def plotfingerplot(any_bam_file, fingerprint_plot, fingerprint_matrix, ini_section='deeptools_QC'):
    return Job(
        input_files=any_bam_file,
        output_files=[fingerprint_plot, fingerprint_matrix],
        module_entries=[['default', 'module_deeptools']],
        command="""\
plotFingerprint --verbose \\
    {options} \\
    --plotFileFormat pdf \\
    --centerReads \\
    --numberOfProcessors {cpu} \\
    --bamfiles {any_bam_file} \\
    --plotFile {fingerprint_plot} \\
    --outRawCounts {fingerprint_matrix}""".format(
            options=global_conf.global_get(ini_section, 'options'),
            cpu=global_conf.global_get('deeptools_QC', 'cluster_cpu', required=True),
            any_bam_file=" ".join(any_bam_file),
            fingerprint_plot=fingerprint_plot,
            fingerprint_matrix=fingerprint_matrix,
        )
  )
