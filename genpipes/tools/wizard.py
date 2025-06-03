#!/usr/bin/env python3

import json
import os
import sys

import questionary 
from jinja2 import Environment 

def load_guide (file_path):
    """
    Load wizard JSON files by filename: {general,deployment,pipeline,protocol,command,step}_guide.json
    """
    full_file_path = os.path.join(os.path.dirname(__file__), "wizard_json", file_path)
    with open(full_file_path) as file:
        return json.load(file)
    
class Wizard:

    def __init__ (self, start_file):
        #store variables: pipeline_name, protocol_name, {r,j,c,o,d,p,s,g}_command, scheduler_server_name, 
        #server_in, path_custom_ini, final_command, step_range
        self.variables = {}

        #for rendering
        self.env = Environment()

        #load initial tree
        self.current_guide = load_guide(start_file)

        #for json file switching
        self.current_file = start_file

        #determines starting node based on node id
        self.current_node_id = self.current_guide["_meta"]["entry_point"]