################################################################################
# Copyright (C) 2014, 2022 GenAP, McGill University and Genome Quebec Innovation Centre
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

# Python Standard Modules
import logging
import re

# MUGQIC Modules

log = logging.getLogger(__name__)

def parse_sequence_dictionary_file(sequence_dictionary_file, variant=False):
    sequence_dictionary = []

    log.info("Parse sequence dictionary " + sequence_dictionary_file + " ...")

    with open(sequence_dictionary_file) as sdf:
        for line in sdf:
            parsed_line = re.search("^@SQ\tSN:([^\t]+)\tLN:(\d+)", line)
            if parsed_line:
                if variant:
                    if "_" in parsed_line.group(1) or "." in parsed_line.group(1) or "EBV" in parsed_line.group(1) or "hs37d5" in parsed_line.group(1):
                        sequence_dictionary.append({'name': parsed_line.group(1), 'length': int(parsed_line.group(2)), 'type': str('alt')})
                        continue
                    else:
                        sequence_dictionary.append({'name': parsed_line.group(1), 'length': int(parsed_line.group(2)),'type': str('primary')})
                else:
                    sequence_dictionary.append({'name': parsed_line.group(1), 'length': int(parsed_line.group(2)), 'type': str('primary')})

    log.info(str(len(sequence_dictionary)) + " sequences parsed\n")

    return sequence_dictionary

def split_by_size(sequence_dictionary, nbSplits, variant=False):
    split_list = []
    altsToExclude = []

    total_genome = 0
    for sequence in sequence_dictionary:
        if variant and sequence['type'] == 'alt':
            altsToExclude.append(sequence['name'])
            
        if variant and sequence['type'] != 'alt':
            total_genome += sequence['length']
            
        else:
            total_genome += sequence['length']

    total = 0
    used_size = 0
    toExcludeChr = []
    currentChrs = []
    for sequence in sequence_dictionary:
        blockSize = int((total_genome - used_size)/nbSplits)
        # Stop if we already reached our limit.
        # This can gappen since the size of chromosomes vary
        if len(split_list) == nbSplits:
            break

        currentChrs.append(sequence['name'])
        total += sequence['length']
        if total > blockSize:
            split_list.append(currentChrs)
            toExcludeChr.extend(currentChrs)
            currentChrs = []
            used_size += total
            total = 0

    # If the split gave a round number remove the last block and set it in other
    if len(split_list) == len(sequence_dictionary):
        toRemove = split_list.pop()
        for sequence in toRemove:
          toExcludeChr.pop()
        
    if variant and altsToExclude:
        toExcludeChr = toExcludeChr + altsToExclude

    return split_list,toExcludeChr
