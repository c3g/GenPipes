#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def fastareformat (input, output):
    cmd = "fastareformat " + input + " > " + output 
    job = Job(output_files=[output], command=cmd, module_entries=[['DEFAULT' , 'module_exonerate']])
    
    return job




