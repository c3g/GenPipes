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
from core.config import config
from core.job import Job

def plotfingerprint(chip_bam, input_bam, sample_name, input_name, chip_name, chip_type, output_dir):
    if input_bam:
        inputs = [chip_bam, input_bam]
        labels = [sample_name + "." + chip_name, sample_name + "." + input_name]
    else:
        inputs = [chip_bam]
        labels = [sample_name + "." + chip_name]

    if chip_type == "narrow":
        bin_size = 200
    elif chip_type == "broad":
        bin_size = 1000
    bin_size = config.param('plotfingerprints', 'bin_size', required=False, type='int') if config.param('plotfingerprints', 'bin_size', required=False, type='int') else bin_size

    plot = os.path.join(output_dir, sample_name, sample_name + "." + chip_name + ".fingerprint.png")
    table = os.path.join(output_dir, sample_name, sample_name + "." + chip_name + ".fingerprint.tsv")

    return Job(
        inputs,
        [plot, table],
        module_entries=[['deeptools', 'module_deeptools']],
        command="""\
plotFingerprint --bamfiles {inputs}{jsd_sample} \\
    --binSize {bin_size} \\
    --labels {labels} \\
    --outQualityMetrics {table} \\
    --plotFile {plot}""".format(
        inputs=" ".join(inputs),
        jsd_sample=" --JSDsample " + input_bam if input_bam else "",
        bin_size=bin_size,
        labels=" ".join(labels),
        table=table,
        plot=plot
        )
    )
