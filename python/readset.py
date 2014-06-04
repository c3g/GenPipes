#!/usr/bin/env python

# Python Standard Modules
import csv
import os
import re

# MUGQIC Modules
from sample import *

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
        print 'Readset -- name: ' + self._name + ', run_type: ' + self._run_type

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

    @staticmethod
    def parse_readset_file(readset_file):
        readsets = []
        samples = []

        readset_csv = csv.DictReader(open(readset_file, 'rb'), delimiter='\t')
        for line in readset_csv:
            sample_name = line['Sample']
            sample_names = [sample.name for sample in samples]
            if sample_name in sample_names:
                # Sample already exists
                sample = samples[sample_names.index(sample_name)]
            else:
                # Create new sample
                sample = Sample(sample_name)
                samples.append(sample)

            # Create readset and add it to sample
            readset = Readset(line['ReadSet'], line['RunType'])

            # Readset file paths are either absolute or relative to the readset file
            # Convert them to absolute paths and check if files exist
            for format in ("BAM", "FASTQ1", "FASTQ2"):
                if line[format]:
                    if not os.path.isabs(line[format]):
                        line[format] = os.path.dirname(os.path.abspath(readset_file)) + os.sep + line[format]

                    if not os.path.isfile(line[format]):
                        raise Exception("Error in parse_readset_file: \"" + line[format] +
                            "\" does not exist or is not a valid plain file!")

            readset.bam = line['BAM']
            readset.fastq1 = line['FASTQ1']
            readset.fastq2 = line['FASTQ2']
            readset.library = line['Library']
            readset.run = line['Run']
            readset.lane = line['Lane']
            readset.adaptor1 = line['Adaptor1']
            readset.adaptor2 = line['Adaptor2']
            readset.quality_offset = line['QualityOffset']
            if line['BED'] != None and line['BED'] != "":
                readset.beds = line['BED'].split(";")
            else:
                readset.beds = []

            readsets.append(readset)
            sample.add_readset(readset)

        return readsets

#readset = Readset("readset1", "SINGLE_END")
#readset.show()
#readset.bam = "toto.bam"
#print readset.bam
#
#sample = Sample("sample1")
#sample.show()
#sample.add_readset(readset)
#sample.show()
#
#Readset.parse_readset_file("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/samples.tsv")
