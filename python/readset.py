#!/usr/bin/env python

# Python Standard Modules
import csv
import re

class Readset:

    def __init__(self, name, run_type):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: readset name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

        if run_type in ("PAIRED_END", "SINGLE_END"):
            self._run_type = run_type
        else:
            raise Exception("Error: readset run_type \"" + run_type +
                "\" is invalid (should be \"PAIRED_END\" or \"SINGLE_END\")!")

    def show(self):
        print 'Basic -- name: ' + self._name + ', run_type: ' + self._run_type

    @property
    def name(self):
        return self._name

    @property
    def run_type(self):
        return self._run_type

    @property
    def sample(self):
        return self._sample

    @property
    def bam(self):
        return self._bam

    @property
    def fastq1(self):
        return self._fastq1

    @property
    def fastq2(self):
        return self._fastq2

    @property
    def library(self):
        return self._library

    @property
    def run(self):
        return self._run

    @property
    def lane(self):
        return self._lane

    @property
    def adaptor1(self):
        return self._adaptor1

    @property
    def adaptor2(self):
        return self._adaptor2

    @property
    def quality_offset(self):
        return self._quality_offset

    @property
    def beds(self):
        return self._beds

readset = Readset('readset1', 'SINGLE_END')
readset.show()
readset.bam = "toto.bam"
print readset.bam

#Readset.parse_readset_file("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/samples.tsv")
