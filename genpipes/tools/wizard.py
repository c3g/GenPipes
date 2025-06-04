#!/usr/bin/env python3

import json
import os
import sys

import questionary 
from jinja2 import Environment 

def load_guide (file_path):
    """
    Load wizard JSON files by filename: {general, deployment, pipeline, protocol, command, step}_guide.json
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

    def apply_variables (self, message):
        """
        Fill the placeholder {{...}} in the message with the data from the current variable
        """
        return self.env.from_string(message).render(**self.variables)
    
    def find_node (self, node_id):
        """
        Find and return the node dictionary with the given id of the node
        """
        for n in self.current_guide["nodes"]:
            if n["id"] == node_id:
                return n
        raise RuntimeError(f"Node '{node_id}' not found in {self.current_file}")
