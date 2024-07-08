################################################################################
# Copyright (C) 2014, 2024 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of GenPipes.
#
# GenPipes is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GenPipes is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with GenPipes.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import os
import re
import json

# MUGQIC Modules
from ..core.config import global_conf

log = logging.getLogger(__name__)

# Start creating the json dump for the passed sample
def create(pipeline, sample):
    jsonator_version = "1.0.2"

    # Retrieve the project name fome the config file, if not specified then use the parent folder name of where the pipeline has been launched
    if global_conf.global_get("DEFAULT", 'project_name', required=False):
        project_name = global_conf.global_get("DEFAULT", 'project_name', required=False)
    else:
        project_name = os.path.basename(pipeline.output_dir)

    # Prepare the general information hash
    general_info = {}
    if pipeline.__class__.__name__ == "AmpliconSeq":
        general_info = {
            'amplicon_type' : global_conf.global_get("DEFAULT", 'amplicon_type'),
            'db_name' : global_conf.global_get("DEFAULT", 'db_name'),
            'db_version' : global_conf.global_get("DEFAULT", 'db_version'),
            'similarity_threshold' : global_conf.global_get("DEFAULT", 'similarity_threshold')
        }
    elif pipeline.__class__.__name__ == "PacBioAssembly":
        general_info = {
            'library_type' : global_conf.global_get("DEFAULT", 'library_type'),
            'blast_db' : global_conf.global_get("DEFAULT", 'blast_db')
        }
    elif pipeline.__class__.__name__ == "Nanopore":
        general_info = {
            'instrument': global_conf.global_get("DEFAULT", 'instrument_type'),
            'blast_db': global_conf.global_get("DEFAULT", 'blast_db')
        }
    elif pipeline.__class__.__name__ == "RnaSeqDeNovoAssembly":
        general_info = {
            'swissprot_db' : global_conf.global_get("DEFAULT", 'swissprot_db'),
            'uniref_db' : global_conf.global_get("DEFAULT", 'uniref_db'),
            'pfam_db' : global_conf.global_get("DEFAULT", 'pfam_db')
        }
    elif pipeline.__class__.__name__ == "IlluminaRunProcessing":
        general_info = {
            'analysed_species' : global_conf.global_get("DEFAULT", 'scientific_name'),
            'assembly_used' : global_conf.global_get("DEFAULT", 'assembly'),
            'assembly_source' : global_conf.global_get("DEFAULT", 'source')
        }
    else :
        general_info = {
            'analysed_species' : global_conf.global_get("DEFAULT", 'scientific_name'),
            'assembly_used' : global_conf.global_get("DEFAULT", 'assembly'),
            'assembly_source' : global_conf.global_get("DEFAULT", 'source')
        }
    if global_conf.global_get("DEFAULT", 'dbsnp_version', required=False) : general_info['dbsnp_version'] = global_conf.global_get("DEFAULT", 'dbsnp_version', required=False)
    general_info['server'] = global_conf.global_get("DEFAULT", 'cluster_server', required=True)
    general_info['analysis_folder'] = pipeline.output_dir + "/"

    # Prepare the software hash by first retrieving all unique module version values in config files
    # assuming that all module key names start with "module_"
    modules = []
    for section in global_conf.sections():
        for name, value in global_conf.items(section):
            if re.search("^module_", name) and value not in modules:
                modules.append(value)

    # Then from the modules, build the list of softwares with names and versions
    softwares = []
    for module in modules:
        softwares.append(
            {
                'name' : module.split("/")[-2],
                'version' : module.split("/")[-1]
            }
        )

    general_info['pipeline_version'] = pipeline.version

    # Check if 'force_jobs' is 'True'
    # Or    if the json file has not been created yet :
    if pipeline.force_jobs or not os.path.exists(os.path.join(pipeline.output_dir, "json", sample.json_file)) or os.stat(os.path.join(pipeline.output_dir, "json", sample.json_file)).st_size == 0:
        # Then (re-)create it !!
        log.debug(f"Creating analysis JSON for sample {sample.name}")
        json_hash = create_json(
            sample,
            pipeline,
            project_name,
            general_info,
            softwares,
            jsonator_version
        )

    # If the json file has already been created (during a previous pipeline execution for instance) :
    else :
        log.debug(f"Opening for update existing analysis JSON for sample '{sample.name}'.")
        try:
            with open(os.path.join(pipeline.output_dir, "json", sample.json_file), 'r') as json_file:
                current_json_hash = json.load(json_file)
            if current_json_hash:
                json_hash = update_json(
                    current_json_hash,
                    sample,
                    pipeline,
                    general_info,
                    softwares
                )
            else:
                log.debug(f"Existing analysis JSON for sample '{sample.name}' is empty... re-creating it.")
                json_hash = create_json(
                    sample,
                    pipeline,
                    project_name,
                    general_info,
                    softwares,
                    jsonator_version
                )

        except json.decoder.JSONDecodeError:
            log.debug(f"Update failed... re-creating analysis JSON for sample '{sample.name}'.")
            json_hash = create_json(
                sample,
                pipeline,
                project_name,
                general_info,
                softwares,
                jsonator_version
            )
    
    current_json = json.dumps(json_hash, indent=4, default=str)

    if not os.path.exists(os.path.join(pipeline.output_dir, "json")):
        os.makedirs(os.path.join(pipeline.output_dir, "json"))

    # Print to file
    filepath = os.path.join(pipeline.output_dir, "json", sample.json_file)
    with open(filepath, 'w') as out_json:
        out_json.write(current_json)
    out_json.close()

    return filepath


def create_json(
    sample,
    pipeline,
    project,
    general_info,
    softwares,
    jsonator_version
    ):
    if pipeline.__class__.__name__ == "PacBioAssembly":
        json_hash = {
            'version': jsonator_version,
            'project': project,
            'submission_date': "",      # Create a submission time entry and let it empty : will be updated as the bash script is launched
            'sample_name' : sample.name,
            'readset' : [
                {
                    "name" : readset.name,
                    "run" : readset.run,
                    "smartcell" : readset.smartcell,
                    "protocol" : readset.protocol,
                    "nb_base_pairs" : readset.nb_base_pairs,
                    "estimated_genome_size" : readset.estimated_genome_size,
                    "bas" : [os.path.realpath(bas) for bas in readset.bas_files],
                    "bax" : [os.path.realpath(bax) for bax in readset.bax_files]
                } for readset in sample.readsets
            ],
            'pipeline' : {
                'name' : pipeline.__class__.__name__,
                'general_information': general_info,
                'software' : [
                    {
                        'name' : software['name'],
                        'version' : software['version']
                    } for software in softwares
                ],
                'step': []
            }
        }
    elif pipeline.__class__.__name__ == "Nanopore" or pipeline.__class__.__name__ == "NanoporeCoVSeq":
        json_hash = {
            'version': jsonator_version,
            'project': project,
            'submission_date': "",
            # Create a submission time entry and let it empty : will be updated as the bash script is launched
            'sample_name': sample.name,
            'readset': [
                {
                    "name": readset.name,
                    "run": readset.run,
                    "flowcell": readset.flowcell,
                    "library": readset.library,
                    "summary": readset.summary_file,
                } for readset in sample.readsets
            ],
            'pipeline': {
                'name': pipeline.__class__.__name__,
                'general_information': general_info,
                'software': [
                    {
                        'name': software['name'],
                        'version': software['version']
                    } for software in softwares
                ],
                'step': []
            }
        }
    else :
        json_hash = {
            'version': jsonator_version,
            'project': project,
            'submission_date': "",      # Create a submission time entry and let it empty : will be updated as the bash script is launched
            'sample_name' : sample.name,
            'readset' : [
                {
                    "name" : readset.name,
                    "library" : readset.library,
                    "runType" : readset.run_type,
                    "run" : readset.run,
                    "lane" : readset.lane,
                    "adapter1" : readset.adapter1,
                    "adapter2" : readset.adapter2,
                    "qualityoffset" : readset.quality_offset,
                    "bed" : [bed for bed in readset.beds],
                    "fastq1" : os.path.realpath(readset.fastq1) if readset.fastq1 else "",
                    "fastq2" : os.path.realpath(readset.fastq2) if readset.fastq2 else "",
                    "bam" : os.path.realpath(readset.bam) if readset.bam else ""
                } for readset in sample.readsets
            ],
            'pipeline' : {
                'name' : pipeline.__class__.__name__,
                'general_information': general_info,
                'software' : [
                    {
                        'name' : software['name'],
                        'version' : software['version']
                    } for software in softwares
                ],
                'step': []
            }
        }
    for step in pipeline.step_to_execute:
        # First verify if the step is meant to be "jsonified"
        jsonify_step = False
        for job in step.jobs:
            if sample in job.samples:
                jsonify_step = True
        if jsonify_step:
            json_hash['pipeline']['step'].append(
                {
                    'name': step.name,
                    'job': [
                        {
                            "name": job.name,
                            "id": job.id,
                            "command": re.sub("\\\\\n", "", job.command_with_modules),
                            "input_file": job.input_files,
                            "output_file": job.output_files,
                            "dependency": [dependency_job.id for dependency_job in job.dependency_jobs]
                        } for job in step.jobs if sample in job.samples
                    ]
                }
            )
    return json_hash

def update_json(
        current_json_hash,
        sample,
        pipeline,
        general_info,
        softwares
    ):

    # Then check if information is up-to-date by comparing it with the previously retrieved informations
    for info_key in general_info.keys():
        if not info_key in current_json_hash['pipeline']['general_information'] or current_json_hash['pipeline']['general_information'][info_key] != general_info[info_key] :
            current_json_hash['pipeline']['general_information'][info_key] = general_info[info_key]

    # And do the same checking with the list of softwares
    for soft in softwares:
        soft_found = False
        for jsoft in current_json_hash['pipeline']['software']:
            if soft['name'] == jsoft['name']:
                soft_found = True
                if soft['version'] != jsoft['version']:
                    jsoft['version'] = soft['version']
        if not soft_found:
            current_json_hash['pipeline']['software'].append(
                {
                    'name' : soft['name'],
                    'version' : soft['version']
                }
            )

    # Finally check if the requested steps/jobs are already in the JSON :
    #   if so  : update them with the current information
    #   if not : add them to the json object
    for step in pipeline.step_to_execute:
        # First make sure the step is meant to be "jsonified"
        jsonify_step = False
        if step.jobs:
            for job in step.jobs:
                if sample in job.samples:
                    jsonify_step = True

        if jsonify_step:
            # Then check if the step is found in the current json
            step_found = False
            for jstep in current_json_hash['pipeline']['step']:
                if step.name == jstep['name']:
                    step_found = True

            # If step is found, then remove it from the json object (so that it can be replaced by the new one if needed)
            if step_found:
                for i in range(len(current_json_hash['pipeline']['step'])):
                    if current_json_hash['pipeline']['step'][i]['name'] == step.name:
                        del current_json_hash['pipeline']['step'][i]
                        break

            # Now it is time to add the current step record (with its jobs) to the json object
            current_json_hash['pipeline']['step'].append(
                {
                    'name': step.name,
                    'job': [
                        {
                            "name": job.name,
                            "id": job.id,
                            "command": re.sub("\\\\\n", "", job.command_with_modules),
                            "input_file": job.input_files,
                            "output_file": job.output_files,
                            "dependency": [dependency_job.id for dependency_job in job.dependency_jobs]
                        } for job in step.jobs if sample in job.samples
                    ]
                }
            )
    return current_json_hash
