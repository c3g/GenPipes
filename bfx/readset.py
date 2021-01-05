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
import xml.etree.ElementTree as Xml
import sys
import subprocess
from Bio.Seq import Seq

# MUGQIC Modules
from .run_processing_aligner import BwaRunProcessingAligner, StarRunProcessingAligner 
from .sample import Sample, NanoporeSample
from core.config import config, _raise, SanitycheckError

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
    def sample_tag(self):
        return self._sample_tag

    @property
    def gender(self):
        return self._gender

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
        if not hasattr(self, "_adapter1"):
            return None
        else:
            return self._adapter1

    @property
    def adapter2(self):
        if not hasattr(self, "_adapter2"):
            return None
        else:
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
    def index_name(self):
        return self._index_name

    @property
    def index(self):
        return self._index

    @property
    def indexes(self):
        return self._indexes

    @property
    def index_type(self):
        return self._index_type

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
    def project_id(self):
        return self._project_id

    @property
    def project(self):
        return self._project

    @property
    def library_source(self):
        return self._library_source

    @property
    def protocol(self):
        return self._protocol

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


def parse_illumina_raw_readset_files(
    output_dir,
    run_dir,
    run_type,
    readset_file,
    lane,
    genome_root,
    nb_cycles,
    index1cycles,
    index2cycles,
    seqtype
    ):

    readsets = []
    samples = []
    skipped_db = []
    GenomeBuild = namedtuple('GenomeBuild', 'species assembly')

    # Parsing Clarity event file
    log.info("Parse Clarity Illumina event file " + readset_file + " ...")

    readset_csv = csv.DictReader(open(readset_file, 'rb'), delimiter='\t', quotechar='"')

    for line in readset_csv:
        protocol_file = config.param('DEFAULT', 'library_protocol_file', type='filepath', required=False)
        if not (protocol_file and os.path.isfile(protocol_file)):
            protocol_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'library_protocol_list.csv')
        protocol_csv = csv.DictReader(open(protocol_file, 'rb'), delimiter=',', quotechar='"')

        adapter_file = config.param('DEFAULT', 'adapter_type_file', type='filepath', required=False)
        if not (adapter_file and os.path.isfile(adapter_file)):
            adapter_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_types.txt')
        adapter_csv = csv.reader(open(adapter_file, 'rb'), delimiter=',', quotechar='"')

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')
        index_csv = csv.reader(open(index_file, 'rb'), delimiter=',', quotechar='"')

        current_lane = line['Position'].split(':')[0]

        if int(current_lane) != lane:
            continue

        sample_name = line['SampleName']

        # Always create a new sample
        sample = Sample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        readset = IlluminaRawReadset(line['SampleName']+"_"+line['LibraryLUID'], run_type)
        readset._quality_offset = 33
        readset._description = line['Index'].split(' ')[0]
        readset._index_name = line['Index']
        readset._library = line['LibraryLUID']
        readset._sample_tag = line['Sample Tag']
        readset._gender = line['Gender']

        for protocol_line in protocol_csv:
            if protocol_line['Clarity Step Name'] == line['LibraryProcess']:
                readset._protocol = line['LibraryProcess']
                readset._library_source = protocol_line['Processing Protocol Name']
                readset._library_type = protocol_line['Library Structure']

                if readset.index_name:
                    if re.search("SI-*", readset._description):
                        key = readset._description
                        for index_line in index_csv:
                            if index_line and index_line[0] == key:
                                readset._index = "-".join(index_line[1:])
                                break
                        else:
                            _raise(SanitycheckError("Could not find index " + key + " in index file " + index_file + " Aborting..."))
                    else:
                        key = readset._description.split("-")[0]
                        index_to_use = ""
                        for idx, index in enumerate(readset._description.split("-")):
                            for index_line in index_csv:
                                if index_line and index_line[0] == index:
                                    if idx > 0 and len(index_line[1]) > 0:
                                        index_to_use += "-"
                                    index_to_use += index_line[1]
                                    break
                            else:
                                _raise(SanitycheckError("Could not find index " + index + " in index file " + index_file + " Aborting..."))
                        readset._index = index_to_use
                    for adapter_line in adapter_csv:
                        if adapter_line:
                            if adapter_line[0] == key:
                                readset._library_type = adapter_line[1]  # TruSeq, Nextera, TenX...
                                readset._index_type = adapter_line[2]    # SINGLEINDEX or DUALINDEX
                                break
                    else:
                        _raise(SanitycheckError("Could not find adapter " + key + " in adapter file " + adapter_file + " Aborting..."))
                    # At this point, inedx_file, adapter_file  protocol_file were succesfully parsed, then exit the loop !
                    break
            else:
                _raise(SanitycheckError("Could not find protocol " + line['LibraryProcess'] + " (from event file " + readset_file + ") in protocol library file " + protocol_file + " for readset " + readset.name + " Aborting..."))

        readset._genomic_database = line['Reference']

        readset._run = Xml.parse(os.path.join(run_dir, "RunInfo.xml")).getroot().find('Run').get('Number')
        readset._lane = current_lane
        readset._sample_number = str(len(readsets) + 1)

        readset._flow_cell = Xml.parse(os.path.join(run_dir, "RunInfo.xml")).getroot().find('Run').find('Flowcell').text
                            #Xml.parse(os.path.join(run_dir, "RunParameters.xml")).getroot().find('RfidsInfo').find('FlowCellSerialBarcode').text
        readset._control = "N"
        readset._recipe = None
        readset._operator = None
        readset._project = line['ProjectName']
        readset._project_id = line['ProjectLUID']
        readset._is_rna = re.search("RNA|cDNA", readset.library_source) or (readset.library_source == "Library" and re.search("RNA", readset.library_type))

        readset._beds = line['Capture REF_BED'].split(";")[1:] if line['Capture REF_BED'] and line['Capture REF_BED'] != "N/A" else []
        
        fastq_file_pattern = os.path.join(
            output_dir,
            "Unaligned." + readset.lane,
            "Project_" + readset.project_id,
            "Sample_" + readset.name,
            readset.name + '_S' + readset.sample_number + "_L00" + readset.lane + "_R{read_number}_001.fastq.gz"
        )
        readset.fastq1 = fastq_file_pattern.format(read_number=1)
        readset.fastq2 = fastq_file_pattern.format(read_number=2) if run_type == "PAIRED_END" else None
        readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1) if index1cycles else None
        readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2) if index2cycles else None

        readset._indexes = get_index(readset, index1cycles, index2cycles, seqtype) if index1cycles else None

        if any(s in readset.protocol for s in ["10X_scRNA", "Single Cell RNA"]):
            readset._is_scrna = True
        else:
            readset._is_scrna = False

        m = re.search("(?P<build>\w+):(?P<assembly>[\w\.]+)", readset.genomic_database)
        genome_build = None
        if m:
            genome_build = GenomeBuild(m.group('build'), m.group('assembly'))

        if genome_build is not None:
            folder_name = os.path.join(genome_build.species + "." + genome_build.assembly)
            current_genome_folder = os.path.join(genome_root, folder_name)

            if readset.is_rna:
                readset._aligner = run_processing_aligner.StarRunProcessingAligner(output_dir, current_genome_folder, nb_cycles)
            else:
                readset._aligner = run_processing_aligner.BwaRunProcessingAligner(output_dir, current_genome_folder)

            aligner_reference_index = readset.aligner.get_reference_index()
            annotation_files = readset.aligner.get_annotation_files()
            reference_file = os.path.join(
                current_genome_folder,
                "genome",
                folder_name + ".fa"
            )
            if reference_file and os.path.isfile(reference_file):
                if aligner_reference_index and (os.path.isfile(aligner_reference_index) or os.path.isdir(aligner_reference_index)):
                    readset._aligner_reference_index = aligner_reference_index
                    readset._annotation_files = annotation_files
                    readset._reference_file = reference_file
                    readset._bam = os.path.join(
                        output_dir,
                        "Aligned." + readset.lane,
                        'alignment',
                        sample_name,
                        'run' + readset.run + "_" + readset.lane,
                        sample_name + "." + readset.library + ".sorted"
                    )

                else:
                    log.warning("Unable to access the aligner reference file: '" + aligner_reference_index +
                                "' for aligner: '" + readset.aligner.__class__.__name__ + "'")
            else:
                log.warning("Unable to access the reference file: '" + reference_file + "'")

        if readset.bam is None and len(readset.genomic_database) > 0 and readset.genomic_database not in skipped_db:
            skipped_db.append(readset.genomic_database)

        readsets.append(readset)
        sample.add_readset(readset)

    if len(skipped_db) > 0:
        log.info("Skipping alignment for the genomic database: '" + "', '".join(skipped_db) + "'")
    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

class MGIReadset(Readset):
    
    def __init__(self, name, run_type):
        super(MGIReadset, self).__init__(name)

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
    def sample_tag(self):
        return self._sample_tag

    @property
    def gender(self):
        if not hasattr(self, "_gender"):
            return None
        else:
            return self._gender

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

def parse_mgi_readset_file(
    mgi_readset_file
    ):

    readsets = []
    samples = []

    log.info("Parsing MGI readset file " + mgi_readset_file + " ...")
    readset_csv = csv.DictReader(open(mgi_readset_file, 'rb'), delimiter='\t')
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
        readset = MGIReadset(line['Readset'], line['RunType'])

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

        readsets.append(readset)
        sample.add_readset(readset)

class MGIRawReadset(MGIReadset):

    def __init__(self, name, run_type):
        super(MGIRawReadset, self).__init__(name, run_type)

    @property
    def index_name(self):
        return self._index_name

    @property
    def index(self):
        return self._index

    @property
    def indexes(self):
        return self._indexes

    @property
    def index_number(self):
        return self._index_number

    @property
    def index_type(self):
        return self._index_type

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
    def is_scrna(self):
        return self._is_scrna

    @property
    def library_source(self):
        return self._library_source

    @property
    def is_mgi_index(self):
       return self._is_mgi_index

    @property
    def library_type(self):
        return self._library_type

    @property
    def protocol(self):
        return self._protocol

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
    def flow_cell(self):
        return self._flow_cell

def parse_mgi_raw_readset_files(
    readset_file,
    bioinfo_file,
    run_type,
    seqtype,
    flowcell,
    lane,
    nb_cycles,
    output_dir
    ):

    readsets = []
    samples = []
    GenomeBuild = namedtuple('GenomeBuild', 'species assembly')

    # Parsing MGI readset sheet
    log.info("Parse MGI readset file " + readset_file + " ...")
    readset_csv = csv.DictReader(open(readset_file, 'rb'), delimiter=',', quotechar='"')
    for line in readset_csv:

        adapter_file = config.param('DEFAULT', 'adapter_type_file', type='filepath', required=False)
        if not (adapter_file and os.path.isfile(adapter_file)):
            adapter_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_types.txt')
        adapter_csv = csv.reader(open(adapter_file, 'rb'), delimiter=',', quotechar='"')

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')
        index_csv = csv.reader(open(index_file, 'rb'), delimiter=',', quotechar='"')

        current_pool = line['Pool_ID']
        current_lane = line['Lane']
        current_flowcell = line['Flowcell_ID']

        if current_pool == "FAIL" or int(current_lane) != int(lane) or current_flowcell != flowcell:
            continue

        sample_name = line['Sample_Name']

        # Always create a new sample
        sample = Sample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        readset = MGIRawReadset(line['Sample_Name'] + "_" + line['Library'], run_type)
        readset._library = line['Library']
        readset._index_name = line['Index']

        # Dual Index
        if re.search("-", readset.index_name) and not re.search("SI-", readset.index_name):
            key = readset.index_name.split('-')[0]
            for idx, index in enumerate(readset.index_name.split("-")):
                for index_line in index_csv:
                    if index_line and index_line[0] == index:
                        if idx > 0 and len(index_line[1]) > 0:
                            index2 = index_line[1]
                        else:
                            index1 = index_line[1]
                        break
                else:
                     _raise(SanitycheckError("Could not find index " + index + " in index file " + index_file + " Aborting..."))
            index_from_lims = [{'INDEX1':index1, 'INDEX2':index2}] 
        # Single-index (pooled or not)
        else:
            key = readset.index_name
            for index_line in index_csv:
                if index_line and index_line[0] == key:
                    index_from_lims = [{'INDEX1':index, 'INDEX2':""} for index in index_line[1:] if len(index) > 0]
                    break
            else:
                _raise(SanitycheckError("Could not find index " + key + " in index file file " + index_file + " Aborting..."))
        readset._index = index_from_lims

        for adapter_line in adapter_csv:
            if adapter_line :
                if adapter_line[0] == key:
                    readset._library_type = adapter_line[1]  # TruSeq, Nextera, TenX...
                    readset._index_type = adapter_line[2]    # SINGLEINDEX or DUALINDEX
                    break
        else:
            _raise(SanitycheckError("Could not find adapter " + key + " in adapter file " + adapter_file + " Aborting..."))

        readset._project = line['Project']
        readset._project_id = line['Project_ID']
        readset._protocol = line['Protocol']
        readset._library_source = line['Library_Source']
        readset._pool_id = line['Pool_ID']
        readset._run = line['RUN_ID']
        readset._sequencer_name = line['Sequencer']
        readset._sequencer_id = line['SequencerID']
        readset._flow_cell = flowcell
        readset._lane = current_lane
        readset._sample_number = str(len(readsets) + 1)
        readset._quality_offset = 33

        readset._genomic_database = config.param('DEFAULT', 'scientific_name', required=True) + "." + config.param('DEFAULT', 'assembly', required=True)

        readset._is_rna = re.search("RNA|cDNA", readset.library_source) or (readset.library_source == "Library" and re.search("RNA", readset.library_type))

        if line.get('BED Files', None):
            readset._beds = line['BED Files'].split(";")
        else:
            readset._beds = []

        fastq_file_pattern = os.path.join(
            output_dir,
            "Unaligned." + readset.lane,
            "Project_" + readset.project,
            "Sample_" + readset.name,
            readset.name + '_S' + readset.sample_number + "_L00" + readset.lane + "_R{read_number}_001.fastq.gz"
        )
        readset.fastq1 = fastq_file_pattern.format(read_number=1)
        readset.fastq2 = fastq_file_pattern.format(read_number=2) if run_type == "PAIRED_END" else None
        readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1)
        readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2)

#        readsets.append(readset)
#        sample.add_readset(readset)

        # Searching for a matching reference for the specified species
#    for readset in readsets:
        readset._indexes = get_index(readset, bioinfo_file, seqtype)

        genome_root = config.param('DEFAULT', 'genome_root', type="dirpath")

        if readset.protocol == "10X_scRNA":
            readset._is_scrna = True
        else:
            readset._is_scrna = False

        if len(readset.genomic_database) > 0 and not readset.is_scrna:
            folder_name = readset.genomic_database
            current_genome_folder = os.path.join(genome_root, folder_name)

            if readset.is_scrna:
                readset._aligner = run_processing_aligner.CellrangerRunProcessingAligner(output_dir, current_genome_folder)
            elif readset.is_rna:
                bioinfo_csv = csv.reader(open(self.bioinfo_files[lane], 'rb'))
                for row in bioinfo_csv:
                    if row[0] == "Read1 Cycles":
                        return row[1]
                _raise(SanitycheckError("Could not get Read 1 Cycles from " + self.bioinfo_files[lane]))
                readset._aligner = run_processing_aligner.StarRunProcessingAligner(output_dir, current_genome_folder, int(nb_cycles))
            else:
                readset._aligner = run_processing_aligner.BwaRunProcessingAligner(output_dir, current_genome_folder)

            aligner_reference_index = readset.aligner.get_reference_index(readset)
            annotation_files = readset.aligner.get_annotation_files()
            reference_file = os.path.join(
                current_genome_folder,
                "genome",
                folder_name + ".fa"
            )
            if reference_file and os.path.isfile(reference_file):
                if aligner_reference_index and (os.path.isfile(aligner_reference_index) or os.path.isdir(aligner_reference_index)):
                    readset._aligner_reference_index = aligner_reference_index
                    readset._annotation_files = annotation_files
                    readset._reference_file = reference_file
                    readset._bam = os.path.join(
                        output_dir,
                        "Aligned." + readset.lane,
                        'alignment',
                        sample_name,
                        'run' + readset.run + "_" + readset.lane,
                        sample_name + "." + readset.library + ".sorted"
                    )

                else:
                    log.warning("Unable to access the aligner reference file: '" + aligner_reference_index +
                                "' for aligner: '" + readset.aligner.__class__.__name__ + "'")
            else:
                log.warning("Unable to access the reference file: '" + reference_file + "'")

        elif readset.is_scrna:
            log.info("Skipping alignment for scRNA assay on " + readset.name)
        elif readset.bam is None and len(readset.genomic_database) > 0:
            log.info("Skipping alignment for the genomic database: '" + readset.genomic_database + "'")

        readsets.append(readset)
        sample.add_readset(readset)

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

def get_index(
        readset,
        index1cycles,
        index2cycles,
        seqtype
        ):
        """
        """

        indexes = []
        index1 = ""
        index2 = ""
        index1seq = ""
        index2seq = ""

        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        index_str_pattern = "grep '%s,' %s | head -n1"
#        if len(readset.index.split('-')) > 2:    # pool of single-index barcodes
#            indexes = readset.index.split('-')
#        else:                                    # either regular single-index or dual-index
#            indexes = [readset.index.replace('-', '')]
        if re.search("-", readset.index_name) and not re.search("SI-", readset.index_name):
            index1 = readset.index_name.split("-")[0]
            index2 = readset.index_name.split("-")[1]
            index1_str = subprocess.check_output(index_str_pattern % (index1, index_file), shell=True).strip()
            index2_str = subprocess.check_output(index_str_pattern % (index2, index_file), shell=True).strip()
            index1_seq = index1_str.split(",")[1:]
            index2_seq = index2_str.split(",")[1:]
            if len(index1_seq) == len(index2_seq):
                for (index1seq, index2seq) in zip(index1_seq, index2_seq):
                    [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = sub_get_index(readset, index1seq, index2seq, index1cycles, index2cycles, seqtype)
                    indexes.append({
                        'SAMPLESHEET_NAME': readset.name,
                        'LIBRARY': readset.library,
                        'PROJECT': readset.project_id,
                        'INDEX_NAME': readset.index_name,
                        'INDEX1': actual_index1seq,
                        'INDEX2': actual_index2seq,
                        'ADAPTERi7' : adapteri7,
                        'ADAPTERi5' : adapteri5
                    })
            else:
                error_msg = "Error: BAD INDEX DEFINITION for " + readset.name + "...\nDUALINDEX barcodes do not contain the same number of index sequences :\n"
                error_msg += index1 + " : " + ",".join(index1_seq) + "\n"
                error_msg += index2 + " : " + ",".join(index2_seq)
                _raise(SanitycheckError(error_msg))

        else:
            index = readset.index_name
            index_str = subprocess.check_output(index_str_pattern % (index, index_file), shell=True).strip()
            index_seq = index_str.split(",")[1:]
            char = ord("A")
            for seq in index_seq:
                if readset.library_type == "tenX_sc_RNA_v1":
                    index2seq = seq
                else:
                    index1seq = seq
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = sub_get_index(readset, index1seq, index2seq, index1cycles, index2cycles, seqtype)
                indexes.append({
                    'SAMPLESHEET_NAME': readset.name if len(index_seq) == 1 else readset.name + "_" + chr(char),
                    'LIBRARY': readset.library,
                    'PROJECT': readset.project_id,
                    'INDEX_NAME': readset.index_name,
                    'INDEX1': actual_index1seq,
                    'INDEX2': actual_index2seq,
                    'ADAPTERi7' : adapteri7,
                    'ADAPTERi5' : adapteri5
                })
                char += 1

        return indexes

def sub_get_index(
        readset,
        index1seq,
        index2seq,
        index1cycles,
        index2cycles,
        seqtype
        ):
        """
        """
        index_file = config.param('DEFAULT', 'index_settings_file', type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

        actual_index1seq='';
        actual_index2seq='';

        primer_seq_pattern = "grep -A8 '%s' %s | grep '_IDX_' | awk -F':' '{print $2}' | tr -d \"35\'\- \" | awk -F',' '{print $%d}'"

        index1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 1), shell=True).strip()
        index1_primeroffset = int(subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 1") % (seqtype, index_file, 2), shell=True).strip())
        index2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 1), shell=True).strip()
        index2_primeroffset = int(subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index 2") % (seqtype, index_file, 2), shell=True).strip())

        main_seq_pattern = "grep -A4 '^%s:' %s | grep -E \"^3'|^5'\" | head -n 1 | sed \"s/5'//g\"  | sed \"s/3'//g\" | tr -d \" '\-\" | grep %s"
        actual_seq_pattern = "echo %s | awk -F\"%s\" '{print $%d}' | sed \"%s\" | cut -c %s"

        present = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, "-c " + index1_primer), shell=True).strip()
        if present == "1" :

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, index1_primer), shell=True).strip()

            if seqtype in ["dnbseqg400", "dnbseqt7"]:
                if len(index1seq) < int(index1cycles):
                    index1_primer_seq = index1_primer[:len(index1seq)-int(index1cycles)]
                else:
                    index1_primer_seq = index1_primer
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, index1_primer_seq, 2, "s/\[i7\]/"+index1seq+"/g", str(index1_primeroffset+1)+"-"+str(index1_primeroffset+int(index1cycles))+" | tr 'ATGC' 'TACG' | rev"), shell=True).strip()
            else:
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, index1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(index1_primeroffset+1)+"-"+str(index1_primeroffset+int(index1cycles))), shell=True).strip()

            if index2cycles:
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, index2_primer), shell=True).strip()

                if seqtype in ["hiseqx", "hiseq4000", "iSeq", "dnbseqg400", "dnbseqt7"]:
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, index2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", str(index2_primeroffset+1)+"-"+str(index2_primeroffset+int(index2cycles))), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, index2_primer, 2, "s/\[i5\]/"+index2seq+"/g", str(index2_primeroffset+1)+"-"+str(index2_primeroffset+int(index2cycles))), shell=True).strip()

        else :
            indexn1_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "1"), shell=True).strip()
            indexn1_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N1").replace("_DIGIT_", "2"), shell=True).strip()
            indexn2_primer = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "1"), shell=True).strip()
            indexn2_primeroffset = subprocess.check_output(primer_seq_pattern.replace("_IDX_", "Index N2").replace("_DIGIT_", "2"), shell=True).strip()

            main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, indexn1_primer), shell=True).strip()

            if seqtype in ["dnbseqg400", "dnbseqt7"]:
                if len(index1seq) < int(index1cycles):
                    indexn1_primer_seq = indexn1_primer[:len(index1seq)-int(index1cycles)]
                else:
                    indexn1_primer_seq = indexn1_primer
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn1_primer_seq, 2, "s/\[i7\]/"+index1seq+"/g", str(indexn1_primeroffset+1)+"-"+str(indexn1_primeroffset+int(index1cycles)))+" | tr 'ATGC' 'TACG' | rev", shell=True).strip()
            else:
                actual_index1seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn1_primer, 2, "s/\[i7\]/"+index1seq+"/g", str(indexn1_primeroffset+1)+"-"+str(indexn1_primeroffset+int(index1cycles))), shell=True).strip()

            if index2cycles:
                main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, indexn2_primer), shell=True).strip()

                if seqtype in ["hiseqx", "hiseq4000", "iSeq", "dnbseqg400", "dnbseqt7"]:
                    actual_index2seq = subprocess.check_output(actual_seq_pattern.replace("| cut -c", "| rev | cut -c") % (main_seq, indexn2_primer, 1, "s/\[i5c\]/$(echo "+index2seq+" | tr 'ATGC' 'TACG' )/g", str(indexn2_primeroffset+1)+"-"+str(indexn2_primeroffset+index2cycles)), shell=True).strip()
                else :
                    actual_index2seq = subprocess.check_output(actual_seq_pattern % (main_seq, indexn2_primer, 2, "s/\[i5\]/"+index2seq+"/g", str(indexn2_primeroffset+1)+"-"+str(indexn2_primeroffset+int(index2cycles))), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern % (readset.library_type, index_file, "\[i7\]"), shell=True).strip()
        adapteri7 = subprocess.check_output("echo \"%s\" | awk -F\'\\\[i7\\\]\' '{print $1}' | awk -F\'\\\]\' '{print $NF}'" % (main_seq), shell=True).strip()

        main_seq = subprocess.check_output(main_seq_pattern.replace("| head -n 1 |", "|") % (readset.library_type, index_file, "\[i5c\]"), shell=True).strip()
        adapteri5 = subprocess.check_output("echo \"%s\" | awk -F\'\\\[i5c\\\]\' '{print $2}' | awk -F\'\\\[\' '{print $1}' | rev" % (main_seq), shell=True).strip()

        return [actual_index1seq, actual_index2seq, adapteri7, adapteri5]

