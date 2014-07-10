#!/usr/bin/env python

# Python Standard Modules
import csv
import logging
import os
import re

# MUGQIC Modules
from sample import *

log = logging.getLogger(__name__)

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
        print('Readset -- name: ' + self._name + ', run_type: ' + self._run_type)

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

def parse_readset_file(readset_file):
    readsets = []
    samples = []

    log.info("Parse readset file " + readset_file + " ...")
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
        # Convert them to absolute paths
        for format in ("BAM", "FASTQ1", "FASTQ2"):
            if line[format] and not os.path.isabs(line[format]):
                line[format] = os.path.dirname(os.path.abspath(readset_file)) + os.sep + line[format]

        readset.bam = line['BAM']
        readset.fastq1 = line['FASTQ1']
        readset.fastq2 = line['FASTQ2']
        readset.library = line['Library']
        readset.run = line['Run']
        readset.lane = line['Lane']
        readset.adaptor1 = line['Adaptor1']
        readset.adaptor2 = line['Adaptor2']
        readset.quality_offset = int(line['QualityOffset'])
        if line['BED']:
            readset.beds = line['BED'].split(";")
        else:
            readset.beds = []

        readsets.append(readset)
        sample.add_readset(readset)

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(readsets)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

def parse_nanuq_readset_file(readset_file):
    readsets = []
    samples = []

    log.info("Parse Nanuq readset file " + readset_file + " ...")
    readset_csv = csv.DictReader(open(readset_file, 'rb'), delimiter=',', quotechar='"')
    for line in readset_csv:
        if line['Status'] and line['Status'] == "Data is valid":
            sample_name = line['Name']
            sample_names = [sample.name for sample in samples]
            if sample_name in sample_names:
                # Sample already exists
                sample = samples[sample_names.index(sample_name)]
            else:
                # Create new sample
                sample = Sample(sample_name)
                samples.append(sample)
    
            # Create readset and add it to sample
            readset = Readset(line['Filename Prefix'], line['Run Type'])
    
            readset.library = line['Library Barcode']
            readset.run = line['Run']
            readset.lane = line['Region']
    
            readset.adaptor1 = line['Adaptor Read 1 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)']
            readset.adaptor2 = line['Adaptor Read 2 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)']
            readset.quality_offset = int(line['Quality Offset'])
            if line['BED Files']:
                readset.beds = line['BED Files'].split(";")
            else:
                readset.beds = []
    
            file_prefix = "raw_reads/{sample_name}/run{readset.run}_{readset.lane}/{sample_name}.{readset.library}.{readset.quality_offset}.".format(sample_name=sample_name, readset=readset)
    
            if line['BAM']:
                line['BAM'] = file_prefix + "bam"
                line['FASTQ1'] = ""
                line['FASTQ2'] = ""
            elif line['FASTQ1']:
                if line['FASTQ2']:
                    line['FASTQ1'] = file_prefix + "pair1.fastq.gz"
                    line['FASTQ2'] = file_prefix + "pair2.fastq.gz"
                else:
                    line['FASTQ1'] = file_prefix + "single.fastq.gz"

            # Readset file paths are either absolute or relative to the readset file
            # Convert them to absolute paths
            for format in ['BAM', 'FASTQ1', 'FASTQ2']:
                if line[format] and not os.path.isabs(line[format]):
                    line[format] = os.path.abspath(line[format])

            readset.bam = line['BAM']
            readset.fastq1 = line['FASTQ1']
            readset.fastq2 = line['FASTQ2']

            readsets.append(readset)
            sample.add_readset(readset)

        else:
            log.warning("Sample Name " + line['Name'] + ", Run ID " + line['Run'] + ", Lane " + line['Region'] + " data is not valid... skipping")

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(readsets)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets
