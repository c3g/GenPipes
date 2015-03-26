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

import re

class Sample:

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: sample name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

        self._readsets = []

    def show(self):
        print("Sample -- name: " + self._name + ", readsets: " +
            ", ".join([readset.name for readset in self._readsets]))

    @property
    def name(self):
        return self._name

    @property
    def readsets(self):
        return self._readsets

    def readsets_by_name(self, name):
        return [readset for readset in self.readsets if readset.name == name]

    def add_readset(self, readset):
        if self.readsets_by_name(readset.name):
            raise Exception("Error: readset name \"" + readset.name +
                "\" already exists for sample \"" + self.name + "\"!")
        else:
            self.readsets.append(readset)
            readset._sample = self
