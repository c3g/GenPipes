#!/usr/bin/env python

# Python Standard Modules
import logging
import re

# MUGQIC Modules

log = logging.getLogger(__name__)

def parse_sequence_dictionary_file(sequence_dictionary_file):
    sequence_dictionary = []

    log.info("Parse sequence dictionary " + sequence_dictionary_file + " ...")

    with open(sequence_dictionary_file) as sdf:
        for line in sdf:
            parsed_line = re.search("^@SQ\tSN:([^\t]+)\tLN:(\d+)", line)
            if parsed_line:
                sequence_dictionary.append({'name': parsed_line.group(1), 'length': int(parsed_line.group(2))})

    log.info(str(len(sequence_dictionary)) + " sequences parsed\n")

    return sequence_dictionary
