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

import os
import json

# MUGQIC Modules
from core.config import *
from core.job import *

# Start creating the json dump for the passed sample
def create(pipeline, sample):
    # First check if json file has already been created (during a previous pipeline execution for instance)
    # If it does :
    if os.path.exists(os.path.join(pipeline.output_dir, "json", sample.json_file)):
        with open(os.path.join(pipeline.output_dir, "json", sample.json_file), 'r') as json_file:
            current_json = json.load(json_file)

    # If the json file has not been created yet :
    else :
        # Retrieve all unique module version values in config files
        # assuming that all module key names start with "module_"
        modules = []
        for section in config.sections():
            for name, value in config.items(section):
                if re.search("^module_", name) and value not in modules:
                    modules.append(value)

        # From the modules, build the list of softwares with names and versions
        softwares = []
        for module in modules:
            softwares.append({
                'name' : module.split("/")[1],
                'version' : module.split("/")[2]
            })

        # Prepare the general information hash
        general_info = {}
        if pipeline.__class__.__name__ == "AmpliconSeq":
            general_info = {
                'amplicon_type' : config.param("DEFAULT", 'amplicon_type'),
                'db_name' : config.param("DEFAULT", 'db_name'),
                'db_version' : config.param("DEFAULT", 'db_version'),
                'similarity_threshold' : config.param("DEFAULT", 'similarity_threshold')
            }
        elif pipeline.__class__.__name__ == "PacBioAssembly":
            general_info = {
                'library_type' : config.param("DEFAULT", 'library_type'),
                'blast_db' : config.param("DEFAULT", 'blast_db')
            }
        else :
            general_info = {
                'analysed_species' : config.param("DEFAULT", 'scientific_name'),
                'assembly_used' : config.param("DEFAULT", 'assembly')
            }
            if config.param("DEFAULT", 'dbsnp_version') : general_info['dbSNP_version'] = config.param("DEFAULT", 'dbsnp_version')

        current_json = json.dumps(
            {'sample' : {
                'name' : sample.name,
                'readset' : [{
                    "name" : readset.name,
                    "library" : readset.library,
                    "runType" : readset.run_type,
                    "run" : readset.run,
                    "lane" : readset.lane,
                    "adapter1" : readset.adapter1,
                    "adapter2" : readset.adapter2,
                    "qualityoffset" : readset.quality_offset,
                    "bed" : [bed for bed in readset.beds],
                    "fastq1" : os.path.realpath(readset.fastq1),
                    "fastq2" : os.path.realpath(readset.fastq2),
                    "bam" : os.path.realpath(readset.bam),
                } for readset in sample.readsets],
                'pipeline' : {
                    'name' : pipeline.__class__.__name__,
                    'general_information': general_info,
                    'software' : [{
                        'name' : software['name'],
                        'version' : software['version']
                    } for software in softwares],
                    'step' : []
                }
            }}, indent=4)

        if not os.path.exists(os.path.join(pipeline.output_dir, "json")):
            os.makedirs(os.path.join(pipeline.output_dir, "json"))

        # Print to file
        with open(os.path.join(pipeline.output_dir, "json", sample.json_file), 'w') as out_json:
            out_json.write(current_json)

    return current_json