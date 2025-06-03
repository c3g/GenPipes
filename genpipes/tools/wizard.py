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
    