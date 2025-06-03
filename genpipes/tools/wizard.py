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
        #hold variables: pipeline_name, protocol_name, 
        self.variables = {}
        self.env = Environment()
        self.current_guide = load_guide(start_file)
        self.current_file = start_file
        self.current_node_id = self.current_guide["_meta"]["entry_point"]


class WizardEngine:
    def __init__(self, start_file):
        # state will hold variables like pipeline_name, c_command, etc.
        self.state = {}
        self.env = Environment()
        self.current_guide = load_guide(start_file)
        self.current_file = start_file
        # entry_point is a node ID (string)
        self.current_node_id = self.current_guide["_meta"]["entry_point"]