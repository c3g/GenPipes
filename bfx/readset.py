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
from collections import namedtuple
import csv
import logging
import os
import re

# MUGQIC Modules
from .run_processing_aligner import BwaRunProcessingAligner, StarRunProcessingAligner 
from .sample import Sample, NanoporeSample
from core.config import global_config_parser, _raise, SanitycheckError

log = logging.getLogger(__name__)

class Readset(object):

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: readset name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

    @property
    def name(self):
        return self._name

    @property
    def sample(self):
        return self._sample


class IlluminaReadset(Readset):

    def __init__(self, name, run_type):
        super(IlluminaReadset, self).__init__(name)

        if run_type in ("PAIRED_END", "SINGLE_END"):
            self._run_type = run_type
        else:
            raise Exception("Error: readset run_type \"" + run_type +
                "\" is invalid (should be \"PAIRED_END\" or \"SINGLE_END\")!")

        self.fastq1 = None
        self.fastq2 = None

    @property
    def run_type(self):
        return self._run_type

    @property
    def bam(self):
        if not hasattr(self, "_bam"):
            return None
        else:
            return self._bam

    @property
    def umi(self):
        if not hasattr(self, "_umi"):
            return None
        else:
            return self._umi

    @property
    def bigwig(self):
        if not hasattr(self, "_bigwig"):
            return None
        else:
            return self._bigwig

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
    def adapter1(self):
        return self._adapter1

    @property
    def adapter2(self):
        return self._adapter2

    @property
    def primer1(self):
        if not hasattr(self, "_primer1"):
            return None
        else:
            return self._primer1

    @property
    def primer2(self):
        if not hasattr(self, "_primer2"):
            return None
        else:
            return self._primer2

    @property
    def quality_offset(self):
        return self._quality_offset

    @property
    def beds(self):
        return self._beds

    # For ChIP-Seq only
    @property
    def mark_name(self):
        if not hasattr(self, "_mark_name"):
            return None
        else:
            return self._mark_name

    @property
    def mark_type(self):
        if not hasattr(self, "_mark_type"):
            return None
        else:
            return self._mark_type

    # @property
    # def input_sample(self):
    #     if not hasattr(self, "_input_sample"):
    #         return None
    #     else:
    #         return self._input_sample

def parse_illumina_readset_file(illumina_readset_file):
    readsets = []
    samples = []

    if isinstance(illumina_readset_file, str):
       readset_fp = open(illumina_readset_file, 'r')
    else:
        readset_fp = illumina_readset_file
        illumina_readset_file = readset_fp.name

    log.info("Parse Illumina readset file " + illumina_readset_file + " ...")

    # Check for duplicate readsets in file
    dup_found_message = checkDuplicateReadsets(illumina_readset_file)
    if dup_found_message:
        # Display message with the path to the corrected readset file and end the pipeline
        log.error(dup_found_message)
        exit(18)

    # If no duplicate readset was found, then parse the readset file
    readset_csv = csv.DictReader(open(illumina_readset_file, 'r'), delimiter='\t')
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
        readset = IlluminaReadset(line['Readset'], line['RunType'])

        # Readset file paths are either absolute or relative to the readset file
        # Convert them to absolute paths
        for format in ("BAM", "FASTQ1", "FASTQ2"):
            if line.get(format, None):
                line[format] = os.path.expandvars(line[format])
                if not os.path.isabs(line[format]):
                    line[format] = os.path.dirname(os.path.abspath(os.path.expandvars(illumina_readset_file))) + os.sep + line[format]
                line[format] = os.path.normpath(line[format])

        readset._bam = line.get('BAM', None)
        readset._umi = line.get('UMI', None)
        readset._bigwig = line.get('BIGWIG', None)
        readset.fastq1 = line.get('FASTQ1', None)
        readset.fastq2 = line.get('FASTQ2', None)
        readset._library = line.get('Library', None)
        readset._run = line.get('Run', None)
        readset._lane = line.get('Lane', None)
        readset._adapter1 = line.get('Adapter1', None)
        readset._adapter2 = line.get('Adapter2', None)
        #ASVA add-on
        readset._primer1 = line.get('primer1', None)
        readset._primer2 = line.get('primer2', None)
        #remove the adapter from the primer sequences
        if readset._primer1 :
            readset._primer1 = readset._primer1.replace(readset._adapter1,"")
        if readset._primer2 :
            readset._primer2 = readset._primer2.replace(readset._adapter2,"")

        readset._quality_offset = int(line['QualityOffset']) if line.get('QualityOffset', None) else None
        readset._beds = line['BED'].split(";") if line.get('BED', None) else []

        # For ChIP-Seq only
        if line.get('MarkType', None):
            readset._mark_name = line.get('MarkName', None)
            if re.search("^[NBI]$", line.get('MarkType', None)):
                readset._mark_type = line.get('MarkType', None)
            else:
                raise Exception("Mark Error: MarkType \"" + line.get('MarkName', None) +
                    "\" is invalid (should be either N, B or I)!")
        # readset._input_sample = line.get('InputSample', None)

        readsets.append(readset)
        sample.add_readset(readset)
        sample.add_mark(readset.mark_name, readset.mark_type)

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

class IlluminaRawReadset(IlluminaReadset):

    def __init__(self, name, run_type):
        super(IlluminaRawReadset, self).__init__(name, run_type)

    @property
    def index(self):
        return self._index

    @property
    def sample_number(self):
        return self._sample_number

    @property
    def aligner(self):
        return self._aligner

    @property
    def aligner_reference_index(self):
        return self._aligner_reference_index

    @property
    def reference_file(self):
        return self._reference_file

    @property
    def is_rna(self):
        return self._is_rna

    @property
    def annotation_files(self):
        if not hasattr(self, "_annotation_files"):
            return None
        else:
            return self._annotation_files

    @property
    def genomic_database(self):
        return self._genomic_database

    @property
    def project(self):
        return self._project

    @property
    def library_source(self):
        return self._library_source

    @property
    def library_type(self):
        return self._library_type

    @property
    def operator(self):
        return self._operator

    @property
    def recipe(self):
        return self._recipe

    @property
    def control(self):
        return self._control

    @property
    def description(self):
        return self._description

    @property
    def flow_cell(self):
        return self._flow_cell


def parse_illumina_raw_readset_files(output_dir, run_type, nanuq_readset_file, casava_sheet_file, lane, genome_root, nb_cycles):
    readsets = []
    samples = []
    GenomeBuild = namedtuple('GenomeBuild', 'species assembly')

    # Parsing Nanuq readset sheet
    log.info("Parse Nanuq Illumina readset file " + nanuq_readset_file + " ...")
    readset_csv = csv.DictReader(open(nanuq_readset_file, 'rb'), delimiter=',', quotechar='"')
    genome_build = None
    for line in readset_csv:
        current_lane = line['Region']

        if int(current_lane) != lane:
            continue

        sample_name = line['Name']

        # Always create a new sample
        sample = Sample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        readset = IlluminaRawReadset(line['ProcessingSheetId'], run_type)
        readset._quality_offset = 33
        readset._library = line['Library Barcode']
        readset._library_source = line['Library Source']
        readset._library_type = line['Library Type']
        readset._genomic_database = line['Genomic Database']

        readset._run = line['Run']
        readset._lane = current_lane
        readset._sample_number = str(len(readsets) + 1)

        readset._is_rna = re.search("RNA|cDNA", readset.library_source) or (readset.library_source == "Library"
                                                                            and re.search("RNA", readset.library_type))

        if line['BED Files']:
            readset._beds = line['BED Files'].split(";")
        else:
            readset._beds = []

        readsets.append(readset)
        sample.add_readset(readset)

    # Parsing Casava sheet
    log.info("Parsing Casava sample sheet " + casava_sheet_file + " ...")
    casava_csv = csv.DictReader(open(casava_sheet_file, 'rb'), delimiter=',')
    for line in casava_csv:
        if int(line['Lane']) != lane:
            continue
        processing_sheet_id = line['SampleID']
        readset = [x for x in readsets if x.name == processing_sheet_id][0]
        readset._flow_cell = line['FCID']
        readset._index = line['Index']
        readset._description = line['Description']
        readset._control = line['Control']
        readset._recipe = line['Recipe']
        readset._operator = line['Operator']
        readset._project = line['SampleProject']

    # Searching for a matching reference for the specified species
    for readset in readsets:
        m = re.search("(?P<build>\w+):(?P<assembly>[\w\.]+)", readset.genomic_database)
        genome_build = None
        if m:
            genome_build = GenomeBuild(m.group('build'), m.group('assembly'))

        if genome_build is not None:
            folder_name = os.path.join(genome_build.species + "." + genome_build.assembly)
            current_genome_folder = genome_root + os.sep + folder_name

            if readset.is_rna:
                readset._aligner = StarRunProcessingAligner(output_dir, current_genome_folder, nb_cycles)
            else:
                readset._aligner = BwaRunProcessingAligner(output_dir, current_genome_folder)

            aligner_reference_index = readset.aligner.get_reference_index()
            annotation_files = readset.aligner.get_annotation_files()
            reference_file = os.path.join(current_genome_folder,
                                          "genome",
                                          folder_name + ".fa")
            if reference_file and os.path.isfile(reference_file):
                if aligner_reference_index and (os.path.isfile(aligner_reference_index) or os.path.isdir(aligner_reference_index)):
                    readset._aligner_reference_index = aligner_reference_index
                    readset._annotation_files = annotation_files
                    readset._reference_file = reference_file
                    readset._bam = os.path.join(output_dir,
                                                "Aligned." + readset.lane,
                                                'alignment',
                                                readset.sample.name,
                                                'run' + readset.run + "_" + readset.lane,
                                                readset.sample.name + "." + readset.library + ".sorted")
                else:
                    log.warning("Unable to access the aligner reference file: '" + aligner_reference_index +
                                "' for aligner: '" + readset.aligner.__class__.__name__ + "'")
            else:
                log.warning("Unable to access the reference file: '" + reference_file + "'")

        if readset.bam is None and len(readset.genomic_database) > 0:
            log.info("Skipping alignment for the genomic database: '" + readset.genomic_database + "'")

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets


class PacBioReadset(Readset):

    @property
    def run(self):
        return self._run

    @property
    def smartcell(self):
        return self._smartcell

    @property
    def protocol(self):
        return self._protocol

    @property
    def nb_base_pairs(self):
        return self._nb_base_pairs

    @property
    def estimated_genome_size(self):
        if self._estimated_genome_size:
            return self._estimated_genome_size
        else:
            raise Exception("Error: readset \"" + self.name + "\" estimated_genome_size is not defined!")

    @property
    def bas_files(self):
        return self._bas_files

    @property
    def bax_files(self):
        return self._bax_files

def parse_pacbio_readset_file(pacbio_readset_file):
    readsets = []
    samples = []

    log.info("Parse PacBio readset file " + pacbio_readset_file + " ...")
    readset_csv = csv.DictReader(open(pacbio_readset_file, 'rb'), delimiter='\t')
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
        readset = PacBioReadset(line['Readset'])

        # Readset file paths are either absolute or relative to the readset file
        # Convert them to absolute paths
        for format in ("BAS", "BAX"):
            if line.get(format, None):
                abs_files = []
                for file in line[format].split(","):
                    file = os.path.expandvars(file)
                    if not os.path.isabs(file):
                        file = os.path.dirname(os.path.abspath(os.path.expandvars(pacbio_readset_file))) + os.sep + file
                    abs_files.append(os.path.normpath(file))
                line[format] = ",".join(abs_files)

        readset._run = line.get('Run', None)
        readset._smartcell = line.get('Smartcell', None)
        readset._protocol = line.get('Protocol', None)
        readset._nb_base_pairs = int(line['NbBasePairs']) if line.get('NbBasePairs', None) else None
        readset._estimated_genome_size = int(line['EstimatedGenomeSize']) if line.get('EstimatedGenomeSize', None) else None
        readset._bas_files = line['BAS'].split(",") if line.get('BAS', None) else []
        readset._bax_files = line['BAX'].split(",") if line.get('BAX', None) else []

        readsets.append(readset)
        sample.add_readset(readset)

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

class NanoporeReadset(Readset):

    @property
    def sample(self):
        return self._sample

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


def parse_nanopore_readset_file(nanopore_readset_file):
    readsets = []
    samples = []

    if isinstance(nanopore_readset_file, str):
       readset_fp = open(nanopore_readset_file, 'r')
    else:
        readset_fp = nanopore_readset_file
        nanopore_readset_file = readset_fp.name

    log.info("Parse Nanopore readset file " + nanopore_readset_file + " ...")
    # Check for duplicate readsets in file
    dup_found_message = checkDuplicateReadsets(nanopore_readset_file)
    if dup_found_message:
        # Display message with the path to the corrected readset file and end the pipeline
        log.error(dup_found_message)
        exit(18)

    readset_csv = csv.DictReader(open(nanopore_readset_file, 'r'), delimiter='\t')
    for line in readset_csv:
        sample_name = line['Sample']
        sample_names = [sample.name for sample in samples]
        # Create new sample
        sample = NanoporeSample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        readset = NanoporeReadset(line['Readset'])

        # Readset file paths are either absolute or relative to the readset file
        # Convert them to absolute paths
        for format in ("Summary", "FASTQ", "FAST5"):
            if line.get(format, None):
                abs_files = []
                for file in line[format].split(","):
                    file = os.path.expandvars(file)
                    if not os.path.isabs(file):
                        file = os.path.dirname(os.path.abspath(os.path.expandvars(nanopore_readset_file))) + os.sep + file
                    abs_files.append(os.path.normpath(file))
                line[format] = ",".join(abs_files)

        readset._sample = Sample(sample_name)
        readset._run = line.get('Run', None)
        sample._run = readset._run
        readset._flowcell = line.get('Flowcell', None)
        sample._flowcell = readset._flowcell
        readset._library = line.get('Library', None)
        sample._library = readset._library
        readset._summary_file = line['Summary'] if line.get('Summary', None) else []
        sample._summary_file = readset._summary_file
        readset._fastq_files = line['FASTQ'] if line.get('FASTQ', None) else []
        sample._fastq_files = readset._fastq_files
        readset._fast5_files = line['FAST5'] if line.get('FAST5', None) else []
        sample._fast5_files = readset._fast5_files
        readset._analysis_name = line.get('AnalysisName', None)
        sample._analysis_name = readset._analysis_name
        sample._barcode = line.get('Barcode', None)

        readsets.append(readset)
        sample.add_readset(readset)

    # log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

def checkDuplicateReadsets(readset_file):
    readset_csv = csv.DictReader(open(readset_file, 'r'), delimiter='\t')

    readset_dict = {}
    for readset_key in [line['Readset'] for line in readset_csv]:
        if readset_key in readset_dict:
            readset_dict[readset_key] += 1
        else:
            readset_dict[readset_key] = 1
    duplicate_readsets = [readset_name for readset_name, readset_count in readset_dict.items() if readset_count > 1]

    # If duplicate readsets are found
    execption_message = ""
    if len(duplicate_readsets) > 0:
        # Rebuild a readset file with unique readset IDs
        genpipes_proposed_readset_file = os.path.join(
            os.path.splitext(os.path.basename(readset_file))[0] + ".genpipes.txt"
        )
        # Set the header
        csv_headers = readset_csv.fieldnames
        writer = csv.DictWriter(
            open(genpipes_proposed_readset_file, 'w'),
            delimiter=str('\t'),
            fieldnames=csv_headers
        )
        writer.writeheader()
        # Set the counter for already written duplicates
        dup_written = {}
        readset_csv = csv.DictReader(open(readset_file, 'r'), delimiter='\t')
        for line in readset_csv:
            # If current redset has no duplicate, just write the line as is
            if readset_dict[line['Readset']] == 1:
                writer.writerow(line)
            # If current readset has duplicates
            else:
                current_readset = line['Readset']
                # Get the number of replicates
                rep_total = readset_dict[current_readset]
                # Get how much was already written
                if not current_readset in dup_written:
                    dup_written[current_readset] = 0
                # Use counter to ensure uniqueness of readset ID
                line['Readset'] += "_" + str(dup_written[current_readset] + 1)
                # Write the corrected line
                writer.writerow(line)
                dup_written[current_readset] += 1
        # Prepare the error message before ending the pipeline
        exception_message = "Error: Readsets should be unique !!\n\tDuplicates found in the readset file \"" + readset_file + "\": \n"
        exception_message += "\t\t\"" + "\", \"".join(duplicate_readsets)+ "\"\n"
        exception_message += "  You should either upadate your readset file with unique readset names,\n"
        exception_message += "  or use \"" + os.path.realpath(genpipes_proposed_readset_file) + "\" (automatically built by the pipeline upon \"" + readset_file + "\"" + " with unique readset IDs)"

    return execption_message

