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

## functions for awk tools ##

## functions for python tools ##
def extract_isize(file):
    with open(file, "r") as fd:
        lines = fd.readlines()
        for i in range(0, len(lines)):
            line = lines[i]
            if line.startswith('MEDIAN_INSERT_SIZE'):
                ne = lines[i + 1]
                isize_mean = ne.split("\t")[4]
                isize_sd = ne.split("\t")[5]
                break

    return isize_mean, isize_sd
