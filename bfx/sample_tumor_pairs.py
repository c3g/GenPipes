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
import csv
import logging
import os
import re
import ConfigParser

# MUGQIC Modules
from sample import *

log = logging.getLogger(__name__)

class SampleTumorPair(object):

    def __init__(self, name, normal, tumor):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
            self._normal = normal
            self._tumor = tumor
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

def parse_tumor_pair_file(tumor_pair_file, samples):
    samples_dict = dict((sample.name, sample) for sample in samples)
    tumor_pairs = dict()

    log.info("Parse Tumor Pair file " + tumor_pair_file + " ...")
    pair_csv = csv.reader(open(tumor_pair_file, 'rb'), delimiter=',')
    for line in pair_csv:
        sample_name = line[0]
        sample_tumor_pair = SampleTumorPair(sample_name, samples_dict[line[1]], samples_dict[line[2]])
        tumor_pairs[sample_name] = sample_tumor_pair

    log.info(str(len(tumor_pairs)) + " tumor pair" + ("s" if len(tumor_pairs) > 1 else "") + " parsed")
    return tumor_pairs
