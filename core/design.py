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

# MUGQIC Modules
from core.config import _raise
from core.sample import *

log = logging.getLogger(__name__)


class Contrast(object):

    def __init__(self, name):
        self._name = name
        self._controls = []
        self._treatments = []

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def controls(self):
        return self._controls

    @controls.setter
    def controls(self, value):
        self._controls = value

    @property
    def treatments(self):
        return self._treatments

    @treatments.setter
    def treatments(self, value):
        self._treatments = value


def parse_chipseq_design_file(design_file, samples):
    design_csv = csv.DictReader(open(design_file, 'r'), delimiter='\t')
    # Skip first column which is Sample
    contrasts = [Contrast(name) for name in design_csv.fieldnames[2:]]

    for line in design_csv:
        sample_name = line['Sample']
        markname = line['MarkName']
        matching_samples = [sample.name + "-.-" + mark_name for sample in samples for mark_name in sample.marks
                            if (sample.name == sample_name and mark_name == markname)]
        if matching_samples:

            # There should be only one matching sample and mark name
           sample = matching_samples[0]

        else:
            _raise(SanitycheckError("Error: Sample " + sample_name + " and MarkName " + mark_name + " in design file " +
                                    design_file + " not found in pipeline samples!"))
        for contrast in contrasts:
            sample_contrast_type = line[contrast.name]
            # Empty or '0' types are ignored
            if not sample_contrast_type or sample_contrast_type == "0":
                pass
            elif sample_contrast_type == "1":
                contrast.controls.append(sample)
            elif sample_contrast_type == "2":
                contrast.treatments.append(sample)
            else:
                _raise(SanitycheckError("Error: invalid value for sample " + sample_name + " and MarkName " + mark_name
                                        + " and contrast " + contrast.name + " in design file " + design_file +
                                        " (should be '1' for control, '2' for treatment, '0' or '' to be ignored)!"))
    for contrast in contrasts:
        log.info("Contrast " + contrast.name + " (controls: " + str(len(contrast.controls)) + ", treatments: "
                 + str(len(contrast.treatments)) + ") created")
    log.info(str(len(contrasts)) + " contrast" + ("s" if len(contrasts) > 1 else "") + " parsed\n")
    return contrasts



def parse_design_file(design_file, samples):

    design_csv = csv.DictReader(open(design_file, 'r'), delimiter='\t')

    # Skip first column which is Sample
    contrasts = [Contrast(name) for name in design_csv.fieldnames[1:]]

    for line in design_csv:

        sample_name = line['Sample']
        matching_samples = [sample for sample in samples if sample.name == sample_name]
        if matching_samples:
            # There should be only one matching sample
            sample = matching_samples[0]
        else:
            _raise(SanitycheckError("Error: sample " + sample_name + " in design file " + design_file + " not found in pipeline samples!"))

        for contrast in contrasts:
            sample_contrast_type = line[contrast.name]
            # Empty or '0' types are ignored
            if not sample_contrast_type or sample_contrast_type == "0":
                pass
            elif sample_contrast_type == "1":
                contrast.controls.append(sample)
            elif sample_contrast_type == "2":
                contrast.treatments.append(sample)
            else:
                _raise(SanitycheckError("Error: invalid value for sample " + sample_name + " and contrast " + contrast.name + " in design file " + design_file + " (should be '1' for control, '2' for treatment, '0' or '' to be ignored)!"))

    for contrast in contrasts:
        log.info("Contrast " + contrast.name + " (controls: " + str(len(contrast.controls)) + ", treatments: " + str(len(contrast.treatments)) + ") created")
    log.info(str(len(contrasts)) + " contrast" + ("s" if len(contrasts) > 1 else "") + " parsed\n")

    return contrasts
