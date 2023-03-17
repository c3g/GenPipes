################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import re

# MUGQIC Modules

log = logging.getLogger(__name__)

class SampleDOvEEPair(object):

    def __init__(self, name, saliva, brush, readsets, multiple_saliva):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
            self._saliva = saliva
            self._brush = brush
            self._readsets = readsets
            self._multiple_saliva = multiple_saliva
        else:
            raise Exception("Error: DOvEE pair  name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

    @property
    def name(self):
        return self._name

    @property
    def saliva(self):
        return self._saliva

    @property
    def brush(self):
        return self._brush

    @property
    def readsets(self):
        return self._readsets

    @property
    def multiple_saliva(self):
        return self._multiple_saliva


def parse_dovee_pair_file(dovee_pair_file, samples):
    samples_dict = dict((sample.name, sample) for sample in samples)
    readsets_dict = dict((sample.name, sample.readsets) for sample in samples)

    dovee_pairs = dict()
    seen = {}
    dup_saliva_check = {}

    pair_csv = csv.reader(open(dovee_pair_file, 'r'), delimiter=',')
    for line in pair_csv:
        saliva = samples_dict[line[1]]

        if saliva.name not in seen:
            seen[saliva.name] = 1
            dup_saliva_check[saliva.name] = 0
        else:
            if seen[saliva.name] == 1:
                dup_saliva_check[saliva.name] = 1
            seen[saliva.name] += 1
    
    log.info("Parse DOvEE Pair file " + dovee_pair_file + " ...")
    pair_csv = csv.reader(open(dovee_pair_file, 'r'), delimiter=',')
    for line in pair_csv:
        sample_name = line[0]
        saliva = samples_dict[line[1]]
        brush = samples_dict[line[2]]

        sample_dovee_pair = SampleDOvEEPair(
            sample_name,
            saliva,
            brush,
            readsets_dict,
            dup_saliva_check[saliva.name]
            )

        dovee_pairs[sample_name] = sample_dovee_pair

    log.info(str(len(dovee_pairs)) + " dovee pair" + ("s" if len(dovee_pairs) > 1 else "") + " parsed")
    return dovee_pairs
