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
import csv
import logging
import os
import re

# MUGQIC Modules
from .sample import Sample

log = logging.getLogger(__name__)

class SampleTumorPair(object):

    def __init__(self, name, normal, tumor, readsets, multiple_normal, pair_profyle, normal_profyle, tumor_profyle):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
            self._normal = normal
            self._tumor = tumor
            self._readsets = readsets
            self._multiple_normal = multiple_normal
            self._pair_profyle = pair_profyle
            self._normal_profyle = normal_profyle
            self._tumor_profyle = tumor_profyle
        else:
            raise Exception("Error: tumor pair  name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

    @property
    def name(self):
        return self._name

    @property
    def normal(self):
        return self._normal

    @property
    def tumor(self):
        return self._tumor

    @property
    def readsets(self):
        return self._readsets

    @property
    def multiple_normal(self):
        return self._multiple_normal

    @property
    def pair_profyle(self):
        return self._pair_profyle

    @property
    def normal_profyle(self):
        return self._normal_profyle

    @property
    def tumor_profyle(self):
        return self._tumor_profyle



def parse_tumor_pair_file(tumor_pair_file, samples, profyle=False):
    samples_dict = dict((sample.name, sample) for sample in samples)
    readsets_dict = dict((sample.name, sample.readsets) for sample in samples)

    tumor_pairs = dict()
    seen = {}
    dup_normal_check = {}

    pair_csv = csv.reader(open(tumor_pair_file, 'r'), delimiter=',')
    for line in pair_csv:
        normal = samples_dict[line[1]]

        if normal.name not in seen:
            seen[normal.name] = 1
            dup_normal_check[normal.name] = 0
        else:
            if seen[normal.name] == 1:
                dup_normal_check[normal.name] = 1
            seen[normal.name] += 1
    
    log.info("Parse Tumor Pair file " + tumor_pair_file + " ...")
    pair_csv = csv.reader(open(tumor_pair_file, 'r'), delimiter=',')
    for line in pair_csv:
        sample_name = line[0]
        normal = samples_dict[line[1]]
        tumor = samples_dict[line[2]]

        if profyle == True:
            profyle_normal = normal.name.split("_")
            profyle_tumor = tumor.name.split("_")
            profyle_pair = profyle_normal[1] + "_" + profyle_tumor[1]

            sample_tumor_pair = SampleTumorPair(
                sample_name,
                normal,
                tumor,
                readsets_dict,
                dup_normal_check[normal.name],
                profyle_pair,
                profyle_normal[0] + "_" + profyle_normal[1],
                profyle_tumor[0] + "_" + profyle_tumor[1]
            )
        else:
            sample_tumor_pair = SampleTumorPair(
                sample_name,
                normal,
                tumor,
                readsets_dict,
                dup_normal_check[normal.name],
                None,
                None,
                None
            )

        tumor_pairs[sample_name] = sample_tumor_pair

    log.info(str(len(tumor_pairs)) + " tumor pair" + ("s" if len(tumor_pairs) > 1 else "") + " parsed")
    return tumor_pairs
