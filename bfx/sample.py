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
import re

# MUGQIC Modules
from core.config import _raise, SanitycheckError

class UniqueName(type):
    """
    This ensure that only one copy if each object
    with the same name exist
    """
    def __call__(cls, name, *args, **kwargs):
        if name not in cls._registered_samples:
            self = cls.__new__(cls, name, *args, **kwargs)
            cls.__init__(self, name, *args, **kwargs)
            cls._registered_samples[name] = self
        return cls._registered_samples[name]

    def __init__(cls, name, *args, **kwargs):
        super().__init__(name, *args, **kwargs)
        cls._registered_samples = {}

class Sample(metaclass=UniqueName):

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            _raise(SanitycheckError("Sample Error: sample name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!"))

        self._readsets = []

        self._json_file = name + ".json"

        self._marks = {}

    def show(self):
        print("Sample -- name: " + self._name + ", readsets: " +
            ", ".join([readset.name for readset in self._readsets]))

    @property
    def name(self):
        return self._name

    @property
    def marks(self):
        return self._marks

    @property
    def analysis_name(self):
        return self._analysis_name

    @property
    def readsets(self):
        return self._readsets

    @property
    def json_file(self):
        return self._json_file

    def readsets_by_name(self, name):
        return [readset for readset in self.readsets if readset.name == name]

    def add_readset(self, readset):
        if self.readsets_by_name(readset.name):
            _raise(SanitycheckError("Sample Error: readset name \"" + readset.name +
                "\" already exists for sample \"" + self.name + "\"!"))
        else:
            self.readsets.append(readset)
            readset._sample = self

    def add_mark(self, mark_name, mark_type):
        if mark_name not in self.marks:
            self.marks[mark_name] = mark_type


class NanoporeSample(Sample):
    """docstring for NanoporeSample"""
    @property
    def run(self):
        return self._run

    @property
    def flowcell(self):
        return self._flowcell

    @property
    def library(self):
        return self._library

    @property
    def summary_file(self):
        return self._summary_file

    @property
    def fastq_files(self):
        return self._fastq_files

    @property
    def fast5_files(self):
        return self._fast5_files

    @property
    def analysis_name(self):
        return self._analysis_name

    @property
    def barcode(self):
        return self._barcode
