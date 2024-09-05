################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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

def create(pipeline, sample):
    """
    Starts creating the json dump for the passed sample
    """

    json_folder = os.path.join(pipeline.output_dir, "json")
    try:
        json_file = os.path.join(json_folder, f"{pipeline.__class__.__name__}.{pipeline.protocol}_{pipeline.timestamp}.json")
    except AttributeError:
        json_file = os.path.join(json_folder, f"{pipeline.__class__.__name__}_{pipeline.timestamp}.json")

    if not os.path.isabs(pipeline.output_dir):
        pipeline_output_dir = os.path.abspath(pipeline.output_dir)
    else:
        pipeline_output_dir = pipeline.output_dir
    # /!\ Can't os.path.join otherwise the 'server://' disappears
    path_prefix = global_conf.global_get("DEFAULT", 'cluster_server', required=True) + "://" + pipeline_output_dir

    log.debug(f"Updating project_tracking JSON {json_file} for sample '{sample.name}'")

    with open(json_file, 'r') as j_file:
        current_json_hash = json.load(j_file)

    # Existing sample in json
    sample_hash_position = [pos for pos, val in enumerate(current_json_hash['sample']) if val['sample_name'] == sample.name]
    if not sample_hash_position:
        sample_json = {
            'sample_name': sample.name,
            'readset': []
            }
        current_json_hash['sample'].append(sample_json)
    else:
        sample_json = current_json_hash['sample'][sample_hash_position[0]]

    for readset in sample.readsets:
        readset_hash_position = [pos for pos, val in enumerate(sample_json['readset']) if val['readset_name'] == readset.name]
        if not readset_hash_position:
            readset_json = {
                'readset_name': readset.name,
                'job': []
                }
            sample_json['readset'].append(readset_json)
        else:
            readset_json = sample_json['readset'][readset_hash_position[0]]
        for step in pipeline.step_to_execute:
            for job in step.jobs:
                if readset in job.readsets:
                    job_hash_position = [pos for pos, val in enumerate(readset_json['job']) if val['job_name'] == job.name]
                    if not job_hash_position:
                        job_json = {
                            'job_name': job.name,
                            'job_start': None,
                            'job_stop': None,
                            'job_status': None,
                            'file': []
                            }
                        readset_json['job'].append(job_json)
                    else:
                        job_json = readset_json['job'][job_hash_position[0]]
                    for output_file in job.output_files:
                        file_hash_position = [pos for pos, val in enumerate(job_json['file']) if val['file_name'] == os.path.basename(output_file)]
                        if not file_hash_position:
                            # /!\ Can't os.path.join for location_uri otherwise the 'server://' disappears
                            _, extension = os.path.splitext(os.path.basename(output_file))
                            if extension:
                                file_json = {
                                    'location_uri': f'{path_prefix}/{output_file}',
                                    'file_name': os.path.basename(output_file)
                                    }
                                job_json['file'].append(file_json)
                        else:
                            file_json = job_json['file'][file_hash_position[0]]

    with open(json_file, 'w', encoding='utf-8') as j_file:
        json.dump(current_json_hash, j_file, ensure_ascii=False, indent=4)

    return json_file

def init(
    operation_name,
    operation_config_version,
    operation_cmd_line,
    operation_config_md5sum,
    operation_config_data,
    pipeline_output_dir,
    timestamp
    ):
    """
    Initializes the json for project_tracking database
    """

    json_output = {
        'project_name': global_conf.global_get("DEFAULT", 'project_name', required=False) if global_conf.global_get("DEFAULT", 'project_name', required=False) else None,
        'operation_config_name': 'genpipes_ini',
        'operation_config_version': operation_config_version.strip(),
        'operation_config_md5sum': operation_config_md5sum,
        'operation_config_data': ''.join(operation_config_data),
        'operation_platform': global_conf.global_get("DEFAULT", 'cluster_server', required=True),
        'operation_cmd_line': operation_cmd_line,
        'operation_name': f'GenPipes_{operation_name}',
        'sample': []
    }

    json_folder = os.path.join(pipeline_output_dir, "json")

    if not os.path.exists(json_folder):
        os.makedirs(json_folder)

    json_file = os.path.join(json_folder, f"{operation_name}_{timestamp}.json")
    log.debug(f"Initializing project_tracking JSON: {json_file}")
    with open(json_file, 'w', encoding='utf-8') as j_file:
        json.dump(json_output, j_file, ensure_ascii=False, indent=4)

    return json_file
