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

# MUGQIC Modules
from core.job import *


def makeDesignFiles (rmapfile, baitmapfile, outfilePrefix, designDir = ".", other_options=""):

    command = """makeDesignFiles.py {other_options}\\
    --rmapfile {rmapfile}\\
    --baitmapfile {baitmapfile}\\
    --designDir {designDir}\\
    --outfilePrefix {outfilePrefix}""" .format(other_options = other_options, 
        rmapfile = rmapfile, 
        baitmapfile = baitmapfile,
        designDir = designDir, 
        outfilePrefix = outfilePrefix)

    return Job(input_files = [rmapfile, baitmapfile],
            output_files = [outfilePrefix + ".nbp", outfilePrefix + ".nbpb", outfilePrefix + ".poe"],
            module_entries = [['create_design_files', 'module_chicago']],
            name = "create_design_files." + os.path.basename(outfilePrefix),
            command = command,
            )


def bam2chicago (bam, baitmap, rmap, sample, other_options=""):

    command = """bam2chicago.sh {other_options}\\
    {bam}\\
    {baitmap}\\
    {rmap}\\
    {sample}""" .format(other_options = other_options,
        bam = bam, 
        baitmap = baitmap,
        rmap = rmap, 
        sample = sample
        )

    return Job(input_files = [bam, baitmap, rmap],
            output_files = [os.path.join(sample, os.path.basename(sample) + ".chinput")],
            module_entries = [['create_input_files', 'module_chicago'], ['create_input_files', 'module_bedtools']],
            name = "create_input_files." + os.path.basename(sample),
            command = command,
            )


def runChicago(design_dir, sample, output_prefix, design_file_prefix, other_options=""):


#runChicago.R -e seqMonk,interBed,washU_text,washU_track --export-order score --design-dir input_files input_files/BT416P4_HiC/BT416P4_HiC.chinput chicago_output

## to assess featues:
#--features-only --en-full-cis-range --en-trans
 
    command = """runChicago.R {other_options}\\
    -e seqMonk,interBed,washU_text,washU_track \\
    --export-order score  \\
    --design-dir {design_dir}\\
    {input} \\
    {output_prefix}\\""" .format(
        other_options = other_options,
        design_dir = design_dir,
        input = os.path.join(design_dir, sample, sample + ".chinput"), 
        output_prefix = output_prefix
        )

    return Job(input_files = [os.path.join(design_dir, design_file_prefix + ".baitmap"), os.path.join(design_dir, design_file_prefix + ".npb"), os.path.join(design_dir, design_file_prefix + ".poe"), os.path.join(design_dir, design_file_prefix + ".nbpb"), os.path.join(design_dir, sample, sample + ".chinput")],
            output_files = [output_prefix],
            module_entries = [['runChicago', 'module_chicago'], ['runChicago', 'module_R']],
            name = "runChicago." + sample,
            command = command,
            )

