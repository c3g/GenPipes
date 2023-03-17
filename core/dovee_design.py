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

# MUGQIC Modules
from core.config import _raise
from core.sample import *

log = logging.getLogger(__name__)


class Contrast(object):

    def __init__(self, name):
        self._name = name
        self._salivas = []
        self._brushes = []

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def salivas(self):
        return self._salivas

    @salivas.setter
    def salivas(self, value):
        self._salivas = value

    @property
    def brushes(self):
        return self._brushes

    @brushes.setter
    def brushes(self, value):
        self._brushes = value


def parse_dovee_design_file(design_file, samples):

    design_csv = csv.DictReader(open(design_file, 'r'), delimiter='\t')

    # Skip first column which is Sample
    contrasts = Contrast(design_csv.fieldnames[1]) 

    for line in design_csv:

        sample_name = line['Sample']
        contrasts_name = line['Source']
        matching_samples = [sample for sample in samples if sample.name == sample_name]
        if matching_samples:
            # There should be only one matching sample
            sample = matching_samples[0]
        else:
            _raise(SanitycheckError("Error: sample " + sample_name + " in design file " + design_file + " not found in pipeline samples!"))

        sample_contrast_type = line['Source']
            # Empty or '0' types are ignored
        if not sample_contrast_type or sample_contrast_type == "0":
            pass
        elif sample_contrast_type == "saliva":
            contrasts.salivas.append(sample)
        elif sample_contrast_type == "brush":
            contrasts.brushes.append(sample)
        else:
            _raise(SanitycheckError("Error: invalid value for sample " + sample_name + " and contrast " + contrasts_name + " in design file " + design_file + " (should be 'saliva' for saliva, 'brush' for brush, '0' or '' to be ignored)!"))

    log.info("Contrast " + contrasts.name + " (salivas: " + str(len(contrasts.salivas)) + ", brushes: " + str(len(contrasts.brushes)) + ") created")

    return contrasts
