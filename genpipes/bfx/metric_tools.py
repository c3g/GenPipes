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
