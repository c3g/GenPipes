#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def fastareformat (input, output):
    job = Job(
        input_files=[input],
        output_files=[output], 
        command="fastareformat " + input + " > " + output, 
        module_entries=[['DEFAULT' , 'module_exonerate']]
    )
    
    return job




