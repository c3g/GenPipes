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

from core.job import *

def wigToBigWig(wigFile, chromSizes):
    output_bigWig = wigFile + ".bigWig"
    
    return Job(
        [wigFile],
        [output_bigWig],
        [['ucsc', 'module_ucsc']],
        name="wigToBigWig",
        command = """\
wigToBigWig \\
  {wigFile} \\
  {chromSizes} \\
  {bigWigFile}""".format(
        wigFile=wigFile,
        chromSizes=chromSizes,
        bigWigFile=output_bigWig
        )
    )

def bigWigInfo(bigWigFile):
    output = "bigwiginfo_"+ bigWigFile + ".txt"

    return Job(
        [bigWigFile],
        [output],
        [['ucsc', 'module_ucsc']],
        name="bigWigInfo",
        command = """\
bigWigInfo \\
  {bigWigFile} > {output}""".format(
        bigWigFile=bigWigFile,
        output=output
        )
    )