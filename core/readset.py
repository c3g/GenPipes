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
from collections import namedtuple, UserDict
import csv
import logging
import os
import re
import xml.etree.ElementTree as Xml
import sys
import subprocess
import json
from itertools import zip_longest

# MUGQIC Modules
from .run_processing_aligner import BwaRunProcessingAligner, StarRunProcessingAligner, CellrangerRunProcessingAligner, AtacRunProcessingAligner
from .sample import Sample, RunProcessingSample, NanoporeSample
from .config import config, _raise, SanitycheckError

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))

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
        if not hasattr(self, "_reference_file"):
            return None
        else:
            return self._reference_file
    
    @property
    def dictionary_file(self):
        return self._dictionary_file

    @property
    def is_rna(self):
        return self._is_rna
    
    @property
    def is_10x(self):
        return self._is_10x

    @property
    def is_atac(self):
        return self._is_atac

    @property
    def is_scrna(self):
        return self._is_scrna

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
    def species(self):
        return self._species

    @property
    def project_id(self):
        return self._project_id

    @property
    def project(self):
        return self._project
    
    @property
    def hercules_project_id(self):
        return self._hercules_project_id
    
    @property
    def hercules_project_name(self):
        return self._hercules_project_name

    @property
    def external_project_id(self):
        return self._external_project_id

    @property
    def external_project_name(self):
        return self._external_project_name
    
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
    def pool_fraction(self):
        return self._pool_fraction

    @property
    def control(self):
        return self._control

    @property
    def description(self):
        return self._description

    @property
    def species(self):
        return self._species

    @property
    def flow_cell(self):
        return self._flow_cell

    @property
    def report_files(self):
        if not hasattr(self, "_report_files"):
            self._report_files = {}
        return self._report_files

    @report_files.setter
    def report_files(self, value):
        self._report_files = value

class MGIRawReadset(IlluminaRawReadset):
    def __init__(self, name, run_type):
        super(MGIRawReadset, self).__init__(name, run_type)

def parse_freezeman_readset_file(
    readset_file,
    run,
    run_type,
    lane,
    seqtype,
    read1cycles,
    read2cycles,
    index1cycles,
    index2cycles,
    output_dirs,
    platform
    ):

    readsets = []
    samples = []
    skipped_db = []
    GenomeBuild = namedtuple('GenomeBuild', 'species assembly')

    if platform == "illumina":
        platform = "Illumina"
    elif 'mgi' in platform:
        platform = "MGI"
    else:
        _raise(SanitycheckError("Unknow sequencing platform : " + platform + ". Aborting..."))


    # Parsing Freezeman event file
    log.info("Parsing Freezeman info file " + readset_file + " for readset in lane " + lane + "...")

    with open(readset_file, "r") as jff:
        readset_json = json.load(jff)

    for rdst in readset_json['samples']:
        adapter_file = config.param('DEFAULT', 'adapter_type_file', param_type='filepath', required=False)
        if not (adapter_file and os.path.isfile(adapter_file)):
            adapter_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_types.txt')
        adapter_csv = csv.reader(open(adapter_file, 'r'), delimiter=',', quotechar='"')

        protocol_file = config.param('DEFAULT', 'library_protocol_file', param_type='filepath', required=False)
        if not (protocol_file and os.path.isfile(protocol_file)):
            protocol_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'library_protocol_list.csv')
        protocol_csv = csv.DictReader(open(protocol_file, 'r'), delimiter=',', quotechar='"')

        current_lane = str(rdst['lane'])

        if int(current_lane) != int(lane):
            continue

        sample_name = rdst['sample_name']

        # Always create a new sample
        sample = RunProcessingSample(sample_name)
        samples.append(sample)

        library = str(rdst['derived_sample_obj_id'])
        # Create readset and add it to sample
        if readset_json['platform'] == 'ILLUMINA':
            if not platform == "Illumina":
                _raise(SanitycheckError(f"Platform conflict ({platform}) ! Are you sur you provided an Illumina info file ?"))
            readset = IlluminaRawReadset(sample_name+"_"+library, run_type)
        elif readset_json['platform'] == 'DNBSEQ':
            if not platform == "MGI":
                _raise(SanitycheckError(f"Platform conflict ({platform}) ! Are you sur you provided an MGI info file ?"))
            readset = MGIRawReadset(sample_name+"_"+library, run_type)
        else:
            _raise(SanitycheckError("Unknow sequencing platform in info file : " + readset_json['platform'] + ". Aborting..."))
        readset._quality_offset = 33
        readset._index_name = rdst['index_name']
        readset._description = rdst['index_set_name']
        readset._index = [{'INDEX1':index1, 'INDEX2':index2} for index1, index2 in zip_longest(rdst['index_sequence_3_prime'], rdst['index_sequence_3_prime'])]
        readset._library = library
        readset._gender = rdst['expected_sex']
        readset._pool_fraction = rdst['pool_volume_ratio']

        readset._protocol = rdst['library_type']
        for protocol_line in protocol_csv:
            if protocol_line['Clarity Step Name'] == readset.protocol:
                readset._library_source = protocol_line['Processing Protocol Name']
                break
        else:
            _raise(SanitycheckError("Could not find protocol '" + readset.protocol + "' (from event file " + readset_file + ") in protocol library file " + protocol_file + " for readset " + readset.name + " Aborting..."))

        key = readset.index_name.split('-')[0] if re.search("-", readset.index_name) and not re.search("SI-", readset.index_name) else readset.index_name
        if readset.index_name:
            for adapter_line in adapter_csv:
                if adapter_line:
                    if adapter_line[0] == key:
                        readset._library_type = adapter_line[1]  # TruSeq, Nextera, TenX...
                        readset._index_type = adapter_line[2]    # SINGLEINDEX or DUALINDEX
                        break
            else:
                _raise(SanitycheckError("Could not find adapter " + key + " in adapter file " + adapter_file + ". Aborting..."))

        readset._genomic_database = rdst['reference'] if 'reference' in rdst else None
        readset._species = rdst['taxon_name']

        readset._run = run
        readset._lane = current_lane
        readset._sample_number = str(len(readsets) + 1)

        readset._flow_cell = readset_json['container_barcode']
        readset._control = "N"
        readset._recipe = None
        readset._operator = None
        readset._project = rdst['project_name']
        readset._project_id = str(rdst['project_obj_id'])
        readset._external_project = rdst['external_project_name']
        readset._external_project_id = rdst['external_project_id']
        readset._is_rna = re.search("RNA|cDNA", readset.library_source) or (readset.library_source == "Library" and re.search("RNA", readset.library_type))
        readset._is_10x = False
        readset._is_atac = False
        # if "10x" in readset.protocol:
        #     readset._is_10x = True
        # if "ATAC" in readset.protocol:
        #     readset._is_atac = True
        if any(s in readset.protocol for s in ["10X_scRNA", "Single Cell RNA"]):
            readset._is_scrna = True
        else:
            readset._is_scrna = False

        if 'capture_bed' in rdst and rdst['capture_bed'] and rdst['capture_bed'] != "N/A":
            readset._beds = rdst['capture_bed'].split(";")[1:]
            if not rdst['capture_bed'].split(";")[0] == readset.genomic_database:
                readset._genomic_database = rdst['capture_bed'].split(";")[0]
        else:
            readset._beds = []

        fastq_file_pattern = os.path.join(
            output_dirs[readset.lane][f"Unaligned.{readset.lane}_directory"],
            "Project_" + readset.project_id,
            "Sample_" + readset.name,
            readset.name + '_S' + readset.sample_number + "_L00" + readset.lane + "_R{read_number}_001.fastq.gz"
        )
        readset.fastq1 = fastq_file_pattern.format(read_number=1)
        readset.fastq2 = fastq_file_pattern.format(read_number=2) if run_type == "PAIRED_END" else None
        readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1) if index1cycles else None
        readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2) if index2cycles else None

        readset._indexes = get_index(readset, index1cycles, index2cycles, seqtype) if index1cycles else None

        # Searching for a matching reference for the specified species
        genome_root = config.param('DEFAULT', 'genome_root', param_type="dirpath")

        m = re.search("(?P<build>\w+):(?P<assembly>[\w\.]+)", readset.genomic_database) if readset.genomic_database else None
        genome_build = None
        if m:
            genome_build = GenomeBuild(m.group('build'), m.group('assembly'))
        # Setting default human ref if needed
        elif "Homo sapiens" in readset.species:
            genome_build = GenomeBuild("Homo_sapiens", "GRCh38")
        # Setting default mouse ref if needed
        elif "Mus musculus" in readset.species:
            genome_build = GenomeBuild("Mus_musculus", "GRCm38")

        if genome_build is not None:
            folder_name = os.path.join(genome_build.species + "." + genome_build.assembly)
            current_genome_folder = os.path.join(genome_root, folder_name)

            if readset.is_10x:
                if readset.is_rna:
                    readset._aligner = CellrangerRunProcessingAligner(output_dirs, current_genome_folder, platform)
                elif readset.is_atac:
                    readset._aligner = AtacRunProcessingAligner(output_dirs, current_genome_folder, platform)
            elif readset.is_rna:
                if readset.is_scrna:
                    readset._aligner = StarRunProcessingAligner(output_dirs, current_genome_folder, int(read2cycles), platform)
                else:
                    readset._aligner = StarRunProcessingAligner(output_dirs, current_genome_folder, int(read1cycles), platform)                
            else:
                readset._aligner = BwaRunProcessingAligner(output_dirs, current_genome_folder, platform)

            aligner_reference_index = readset.aligner.get_reference_index()
            annotation_files = readset.aligner.get_annotation_files()
            reference_file = readset.aligner.get_reference_file()
            dictionary_file = readset.aligner.get_dictionary_file()
            if reference_file and os.path.isfile(reference_file):
                if aligner_reference_index and (os.path.isfile(aligner_reference_index) or os.path.isdir(aligner_reference_index)):
                    readset._aligner_reference_index = aligner_reference_index
                    readset._annotation_files = annotation_files
                    readset._reference_file = reference_file
                    readset._dictionary_file = dictionary_file
                    readset._bam = os.path.join(
                        output_dirs[readset.lane][f"Aligned.{readset.lane}_directory"],
                        'alignment',
                        sample_name,
                        'run' + readset.run + "_" + readset.lane,
                        sample_name + "_" + readset.library + ".sorted"
                    )

                else:
                    log.warning("Unable to access the aligner reference file: '" + aligner_reference_index + "' for aligner: '" + readset.aligner.__class__.__name__ + "'")
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

def parse_clarity_readset_file(
    readset_file,
    run,
    run_type,
    lane,
    seqtype,
    read1cycles,
    read2cycles,
    index1cycles,
    index2cycles,
    output_dirs,
    platform
    ):

    readsets = []
    samples = []
    skipped_db = []
    GenomeBuild = namedtuple('GenomeBuild', 'species assembly')

    # Parsing Clarity event file
    log.info(f"Parsing Clarity event file {readset_file} for readset in lane {lane}...")

    readset_csv = csv.DictReader(open(readset_file, 'r'), delimiter='\t', quotechar='"')

    for line in readset_csv:
        protocol_file = config.param('DEFAULT', 'library_protocol_file', param_type='filepath', required=False)
        if not (protocol_file and os.path.isfile(protocol_file)):
            protocol_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'library_protocol_list.csv')
        protocol_csv = csv.DictReader(open(protocol_file, 'r'), delimiter=',', quotechar='"')

        adapter_file = config.param('DEFAULT', 'adapter_type_file', param_type='filepath', required=False)
        if not (adapter_file and os.path.isfile(adapter_file)):
            adapter_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_types.txt')
        adapter_csv = csv.reader(open(adapter_file, 'r'), delimiter=',', quotechar='"')

        index_file = config.param('DEFAULT', 'index_settings_file', param_type='filepath', required=False)
        if not (index_file and os.path.isfile(index_file)):
            index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')
        index_csv = csv.reader(open(index_file, 'r'), delimiter=',', quotechar='"')

        current_lane = line['Position'].split(':')[0]

        if int(current_lane) != int(lane):
            continue

        sample_name = line['SampleName']

        # Always create a new sample
        sample = RunProcessingSample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        if platform == 'illumina':
            readset = IlluminaRawReadset(line['SampleName']+"_"+line['LibraryLUID'], run_type)
        elif 'mgi' in platform:
            readset = MGIRawReadset(line['SampleName']+"_"+line['LibraryLUID'], run_type)
        readset._quality_offset = 33
        readset._description = line['Index'].split(' ')[0]
        readset._index_name = line['Index']
        readset._library = line['LibraryLUID']
        readset._gender = line['Gender'] if line['Gender'] else 'N/A'

        for protocol_line in protocol_csv:
            if protocol_line['Clarity Step Name'] == line['LibraryProcess']:
                readset._protocol = line['LibraryProcess']
                readset._library_source = protocol_line['Processing Protocol Name']
                readset._library_type = protocol_line['Library Structure']

                if readset.index_name:
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
            _raise(SanitycheckError("Could not find protocol '" + line['LibraryProcess'] + "' (from event file " + readset_file + ") in protocol library file " + protocol_file + " for readset " + readset.name + " Aborting..."))

        readset._genomic_database = line['Reference']
        readset._species = line['Species']

        readset._run = run
        readset._lane = current_lane
        readset._sample_number = str(len(readsets) + 1)
        readset._flow_cell = line['ContainerName']
        readset._control = "N"
        readset._recipe = None
        readset._operator = None
        readset._project = line['ProjectName']
        readset._project_id = line['ProjectLUID']
        readset._pool_fraction = float(line['Pool Fraction'])

        readset._is_rna = re.search("RNA|cDNA", readset.library_source) or (readset.library_source == "Library" and re.search("RNA", readset.library_type))
        readset._is_10x = False
        readset._is_atac = False
        # if "10x" in readset.protocol:
        #     readset._is_10x = True
        # if "ATAC" in readset.protocol:
        #     readset._is_atac = True
        if any(s in readset.protocol for s in ["10X_scRNA", "Single Cell RNA"]):
            readset._is_scrna = True
        else:
            readset._is_scrna = False

        if line['Capture REF_BED'] and line['Capture REF_BED'] != "N/A":
            readset._beds = line['Capture REF_BED'].split(";")[1:]
            if not line['Capture REF_BED'].split(";")[0] == readset.genomic_database:
                readset._genomic_database = line['Capture REF_BED'].split(";")[0]
        else:
            readset._beds = []

        fastq_file_pattern = os.path.join(
            output_dirs[readset.lane][f"Unaligned.{readset.lane}_directory"],
            "Project_" + readset.project_id,
            "Sample_" + readset.name,
            readset.name + '_S' + readset.sample_number + "_L00" + readset.lane + "_R{read_number}_001.fastq.gz"
        )
        readset.fastq1 = fastq_file_pattern.format(read_number=1)
        readset.fastq2 = fastq_file_pattern.format(read_number=2) if run_type == "PAIRED_END" else None
        readset.index_fastq1 = re.sub("_R1_", "_I1_", readset.fastq1) if index1cycles else None
        readset.index_fastq2 = re.sub("_R2_", "_I2_", readset.fastq2) if index2cycles else None

        readset._indexes = get_index(readset, index1cycles, index2cycles, seqtype) if index1cycles else None

        # Searching for a matching reference for the specified species
        genome_root = config.param('DEFAULT', 'genome_root', param_type="dirpath")

        m = re.search("(?P<build>\w+):(?P<assembly>[\w\.]+)", readset.genomic_database)
        genome_build = None
        if m:
            genome_build = GenomeBuild(m.group('build'), m.group('assembly'))
        # Setting default human ref if needed
        elif "Homo sapiens" in readset.species:
            genome_build = GenomeBuild("Homo_sapiens", "GRCh38")
        # Setting default mouse ref if needed
        elif "Mus musculus" in readset.species:
            genome_build = GenomeBuild("Mus_musculus", "GRCm38")

        if genome_build is not None:
            folder_name = os.path.join(genome_build.species + "." + genome_build.assembly)
            current_genome_folder = os.path.join(genome_root, folder_name)

            if platform == "illumina":
                common_platform = "Illumina"
            elif "mgi" in platform:
                common_platform = "MGI"

            if readset.is_10x:
                if readset.is_rna:
                    readset._aligner = CellrangerRunProcessingAligner(output_dirs, current_genome_folder, common_platform)
                elif readset.is_atac:
                    readset._aligner = AtacRunProcessingAligner(output_dirs, current_genome_folder, common_platform)
            elif readset.is_rna:
                if readset.is_scrna:
                    readset._aligner = StarRunProcessingAligner(output_dirs, current_genome_folder, int(read2cycles), common_platform)
                else:
                    readset._aligner = StarRunProcessingAligner(output_dirs, current_genome_folder, int(read1cycles), common_platform)                
            else:
                readset._aligner = BwaRunProcessingAligner(output_dirs, current_genome_folder, common_platform)

            aligner_reference_index = readset.aligner.get_reference_index()
            annotation_files = readset.aligner.get_annotation_files()
            reference_file = readset.aligner.get_reference_file()
            dictionary_file = readset.aligner.get_dictionary_file()
            if reference_file and os.path.isfile(reference_file):
                if aligner_reference_index and (os.path.isfile(aligner_reference_index) or os.path.isdir(aligner_reference_index)):
                    readset._aligner_reference_index = aligner_reference_index
                    readset._annotation_files = annotation_files
                    readset._reference_file = reference_file
                    readset._dictionary_file = dictionary_file
                    readset._bam = os.path.join(
                        output_dirs[readset.lane][f"Aligned.{readset.lane}_directory"],
                        'alignment',
                        sample_name,
                        'run' + readset.run + "_" + readset.lane,
                        sample_name + "_" + readset.library + ".sorted"
                    )

                else:
                    log.warning("Unable to access the aligner reference file: '" + aligner_reference_index + "' for aligner: '" + readset.aligner.__class__.__name__ + "'")
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
    readset_csv = csv.DictReader(open(pacbio_readset_file, 'r'), delimiter='\t')
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
    Builds an array of dict defining all the indexes of the readset
    """

    indexes = []
    index1 = ""
    index2 = ""
    index1seq = ""
    index2seq = ""

    index_file = config.param('DEFAULT', 'index_settings_file', param_type='filepath', required=False)
    if not (index_file and os.path.isfile(index_file)):
        index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')

    index_str_pattern = "grep '%s,' %s | head -n1"
    if re.search("-", readset.index_name) and not re.search("SI-", readset.index_name):
        index1 = readset.index_name.split("-")[0]
        index2 = readset.index_name.split("-")[1]
        index1_str = subprocess.check_output(index_str_pattern % (index1, index_file), shell=True, text=True).strip()
        index2_str = subprocess.check_output(index_str_pattern % (index2, index_file), shell=True, text=True).strip()
        index1_seq = index1_str.split(",")[1:]
        index2_seq = index2_str.split(",")[1:]
        if len(index1_seq) == len(index2_seq):
            for (index1seq, index2seq) in zip(index1_seq, index2_seq):
                [actual_index1seq, actual_index2seq, adapteri7, adapteri5] = sub_get_index(readset, index1seq, index2seq, index1cycles, index2cycles, seqtype)
                indexes.append(
                    {
                        'SAMPLESHEET_NAME': readset.name,
                        'LIBRARY': readset.library,
                        'PROJECT': readset.project_id,
                        'INDEX_NAME': readset.index_name,
                        'INDEX1': actual_index1seq,
                        'INDEX2': actual_index2seq,
                        'ADAPTERi7' : adapteri7,
                        'ADAPTERi5' : adapteri5
                    }
                )
        else:
            error_msg = "Error: BAD INDEX DEFINITION for " + readset.name + "...\nDUALINDEX barcodes do not contain the same number of index sequences :\n"
            error_msg += index1 + " : " + ",".join(index1_seq) + "\n"
            error_msg += index2 + " : " + ",".join(index2_seq)
            _raise(SanitycheckError(error_msg))

    else:
        index = readset.index_name
        index_str = subprocess.check_output(index_str_pattern % (index, index_file), shell=True, text=True).strip()
        index_seq = index_str.split(",")[1:]
        char = ord("A")
        for seq in index_seq:
            if readset.library_type == "tenX_sc_RNA_v1" or readset.library_type == "TELL-Seq" or readset.library_type == "SHARE-Seq_ATAC" or readset.library_type == "SHARE-Seq_RNA":
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
    Constructs the actual sequence of the indexes
    """
    index_file = config.param('DEFAULT', 'index_settings_file', param_type='filepath', required=False)
    if not (index_file and os.path.isfile(index_file)):
        index_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'adapter_settings_format.txt')
    index_fh = open(index_file, 'r')
    index_line = index_fh.readline()
    while index_line:
        if (len(index_line.strip().split(', ')) > 1) and (seqtype in index_line.strip().split(', ')):
            readt1_def = index_fh.readline()
            readt2_def = index_fh.readline()
            indext1_def = index_fh.readline().split(' - ')[1]
            indext2_def = index_fh.readline().split(' - ')[1]
            readn1_def = index_fh.readline()
            readn2_def = index_fh.readline()
            indexn1_def = index_fh.readline().split(' - ')[1]
            indexn2_def = index_fh.readline().split(' - ')[1]
            [indext1_primer, indext1_primeroffset] = indext1_def.split(',')
            [indext2_primer, indext2_primeroffset] = indext2_def.split(',')
            [indexn1_primer, indexn1_primeroffset] = indexn1_def.split(',')
            [indexn2_primer, indexn2_primeroffset] = indexn2_def.split(',')
            break
        else:
            index_line = index_fh.readline()

    actual_index1seq=''
    actual_index2seq=''

    # Get the main sequence patterns
    while index_line:
        if index_line.startswith("##"):
            library_type_def = index_fh.readline()
            if library_type_def.split(':')[0] == readset.library_type:
                empty_line = index_fh.readline()
                fwd_line_def = index_fh.readline()
                index1_main_seq = re.sub("[\s|\-|']", '', re.sub("3'", "", re.sub("5'", "", fwd_line_def)))
                adapteri7 = index1_main_seq.split("[i7]")[0].split("]")[-1]
                if indext1_primer in fwd_line_def:
                    index1_primer = indext1_primer
                    index1_primeroffset = int(indext1_primeroffset)
                else:
                    index1_primer = indexn1_primer
                    index1_primeroffset = int(indexn1_primeroffset)

                rev_line_def = index_fh.readline()
                index2_main_seq = re.sub("[\s|\-|']", '', re.sub("3'", "", re.sub("5'", "", rev_line_def)))
                adapteri5 = index2_main_seq.split("[i5c]")[1].split("[")[0][::-1]
                if (indext2_primer in fwd_line_def) or (indext2_primer in rev_line_def):
                    index2_primer = indext2_primer
                    index2_primeroffset = int(indext2_primeroffset)
                    if indext2_primer in fwd_line_def:
                        index2_main_seq = index1_main_seq
                        adapteri5 = index2_main_seq.split("[i5]")[1].split("[")[0].replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]
                else:
                    index2_primer = indexn2_primer
                    index2_primeroffset = int(indexn2_primeroffset)
                    
                break
        else:
            index_line = index_fh.readline()
    index_fh.close()

    if index1cycles:
        if len(index1seq) < int(index1cycles):
            index1_primer_seq = index1_primer[:len(index1seq)-int(index1cycles)]
        else:
            index1_primer_seq = index1_primer
        index1_primer_seq = index1_primer
        if index1_primer_seq:
            if readset.library_type == 'tenX_sc_RNA_v1' or readset.library_type == 'TELL-Seq' or readset.library_type == "SHARE-Seq_ATAC" or readset.library_type == "SHARE-Seq_RNA":
                actual_index1seq = ""
            elif seqtype in ["dnbseqg400", "dnbseqt7"] and readset.run_type == "PAIRED_END":
                actual_index1seq = re.sub("\[", "", re.sub("\]", "", re.sub("i7", index1seq, index1_main_seq.split(index1_primer_seq)[1])))[index1_primeroffset:index1_primeroffset+len(index1seq)].replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1][:int(index1cycles)]
            else:
                actual_index1seq = re.sub("\[", "", re.sub("\]", "", re.sub("i7", index1seq, index1_main_seq.split(index1_primer_seq)[1])))[index1_primeroffset:index1_primeroffset+int(index1cycles)]

    if index2cycles:
        if len(index2seq) < int(index2cycles):
            index2_primer_seq = index2_primer[:len(index2seq)-int(index2cycles)]
        else:
            index2_primer_seq = index2_primer
        index2_primer_seq = index2_primer
        if index2_primer_seq:
            if seqtype in ["hiseqx", "hiseq4000", "iSeq"] or (seqtype in ["dnbseqg400", "dnbseqt7"] and readset.run_type == "PAIRED_END"):
                actual_index2seq = re.sub("\[", "", re.sub("\]", "", re.sub("i5c", index2seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper(), index2_main_seq.split(index2_primer_seq)[0])))[::-1][index2_primeroffset:index2_primeroffset+int(index2cycles)]
            else:
                actual_index2seq = re.sub("\[", "", re.sub("\]", "", re.sub("i5", index2seq, index2_main_seq.split(index2_primer_seq)[1])))[index2_primeroffset:index2_primeroffset+int(index2cycles)]

    return [actual_index1seq, actual_index2seq, adapteri7, adapteri5]
