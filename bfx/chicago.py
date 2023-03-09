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

# MUGQIC Modules
from core.job import *


def makeDesignFiles (rmapfile, baitmapfile, file_prefix, design_dir=".", other_options=""):

    command = """makeDesignFiles.py {other_options} \\
    --rmapfile {rmapfile} \\
    --baitmapfile {baitmapfile} \\
    --designDir {design_dir} \\
    --outfilePrefix {file_prefix}""".format(
        other_options=other_options,
        rmapfile=rmapfile,
        baitmapfile=baitmapfile,
        design_dir=design_dir,
        file_prefix=file_prefix)

    return Job(
                input_files=[rmapfile, baitmapfile],
                output_files=[file_prefix + ".npb", file_prefix + ".poe", file_prefix + ".nbpb"],
                module_entries=[['create_design_files', 'module_chicago'], ['create_design_files', 'module_python']],
                name="create_design_files." + os.path.basename(file_prefix),
                command=command
                )


def bam2chicago (bam, baitmap, rmap, sample, other_options=""):

    command = """bam2chicago.sh {other_options} \\
    {bam} \\
    {baitmap} \\
    {rmap} \\
    {sample}""".format(
        other_options=other_options,
        bam=bam,
        baitmap=baitmap,
        rmap=rmap,
        sample=sample
        )

    return Job(
                input_files=[bam, baitmap, rmap],
                output_files=[os.path.join(sample, os.path.basename(sample) + ".chinput")],
                module_entries=[['create_input_files', 'module_chicago'], ['create_input_files', 'module_bedtools']],
                name="create_input_files." + os.path.basename(sample),
                command=command
                )


def runChicago(design_dir, sample, output_dir, design_file_prefix, other_options=""):

    command = """runChicago.R {other_options} \\
    -e seqMonk,interBed,washU_text,washU_track \\
    --export-order score \\
    --design-dir {design_dir} \\
    -o {output} \\
    {input} \\
    {output_prefix}""".format(
        other_options=other_options,
        design_dir=design_dir,
        output=os.path.join(output_dir, sample),
        input=os.path.join(design_dir, sample, sample + ".chinput"),
        output_prefix=sample
        )

    return Job(
                input_files=[os.path.join(design_dir, design_file_prefix + ".baitmap"),
                            os.path.join(design_dir, design_file_prefix + ".npb"),
                            os.path.join(design_dir, design_file_prefix + ".poe"),
                            os.path.join(design_dir, design_file_prefix + ".nbpb"),
                            os.path.join(design_dir, sample, sample + ".chinput")],
                output_files=[os.path.join(output_dir, sample, "data", sample + ".Rds"),
                            os.path.join(output_dir, sample, "data", sample + ".ibed"),
                            os.path.join(output_dir, sample, "data", sample + "_washU_track.txt"),
                            os.path.join(output_dir, sample, "data", sample + "_washU_text.txt"),
                            os.path.join(output_dir, sample, "data", sample + "_washU_track.txt.gz.tbi"),
                            os.path.join(output_dir, sample, "data", sample + "_seqmonk.txt")],
                module_entries=[['runChicago', 'module_chicago'], ['runChicago', 'module_R']],
                name="runChicago." + sample,
                command=command
                )

def runChicago_featureOverlap(featuresBed, sample, output_dir, design_file_prefix, other_options=""):

    command = """runChicago.R {other_options} \\
    --features-only \\
    --en-feat-list {featuresBed} \\
    --en-full-cis-range \\
    --en-trans \\
    -o {output_dir} \\
    {input} \\
    {output_prefix}""".format(
        other_options=other_options,
        featuresBed=featuresBed,
        output_dir=os.path.join(output_dir, sample),
        input=os.path.join(output_dir, sample, "data", sample + ".Rds"),
        output_prefix=sample
        )

    return Job(
                input_files=[os.path.join(output_dir, sample, "data", sample + ".Rds"), featuresBed],
                output_files=[os.path.join(output_dir, sample, "enrichment_data", sample + "_feature_overlaps.pdf"),
                            os.path.join(output_dir, sample, "enrichment_data", sample + "_feature_overlaps.txt")],
                module_entries=[['runChicago', 'module_chicago'], ['runChicago', 'module_R']],
                name="runChicago_featureOverlap." + sample,
                command=command
                )
