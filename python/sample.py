#!/usr/bin/env python

import os
import re

from readset import *

class Sample:

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: sample name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

        self._readsets = []

    def show(self):
        print 'Basic -- name: ' + self._name + ", readsets: " + \
            ", ".join([readset.name for readset in self._readsets])

    @property
    def name(self):
        return self._name

    @property
    def readsets(self):
        return self._readsets

    def get_readset_by_name(self, name):
        readset_names = [readset.name for readset in self._readsets]
        # Should find 1 readset at most since readset names are unique
        if name in readset_names:
            return self._readsets[readset_names.index(name)]
        else:
            return None

    def add_readset(self, readset):
        if self.get_readset_by_name(readset.name):
            raise Exception("Error: readset name \"" + readset.name +
                "\" already exists for sample \"" + self.name + "!")
        else:
            self._readsets.append(readset)
            readset.sample = self

    @staticmethod
    def parse_readset_file(readset_file):
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
                if line[format] != None and line[format] != "":
                    if not os.path.isabs(line[format]):
                        print line[format]
                        line[format] = os.path.abspath(os.path.dirname(readset_file) + os.sep + line[format])
                        print line[format]

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

            sample.add_readset(readset)

        return samples


sample1 = Sample('sample1')
sample1.show()
sample1.bam = "toto.bam"
print sample1.bam
readset1 = Readset("readset1", "SINGLE_END")
sample1.add_readset(readset1)
sample1.show()

samples = Sample.parse_readset_file("/lb/project/mugqic/projects/jfillon_pipelines/dnaseq/bam2fastq/samples.tsv")
print ", ".join([sample.name for sample in samples])
for sample in samples:
    print ", ".join([readset.bam for readset in sample.readsets])
