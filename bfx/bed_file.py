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

def parse_bed_file(bed_file):
    bed_intervals = []

    log.info("Parse bed file" + bed_file + " ...")
    
    total=0
    with open(bed_file) as bf:
        for line in bf:
            parsed_line = line.split()
            if parsed_line:
                bed_intervals.append({'chr': parsed_line[0], 'start':  parsed_line[1], 'end': parsed_line[2] ,'length': int(int(parsed_line[2])-int(parsed_line[1]))})
                total+=int(parsed_line[2])-int(parsed_line[1])

    log.info(str(len(bed_intervals)) + " Intervals parsed\n")

    return bed_intervals, total

def split_by_size(bed_intervals, interval_size, nbSplits, output="tmp"):
    
    blockSize = int(interval_size/nbSplits)
    
    total = 0
    splitNum = 0
    bed_file_list=[]
    out=open(output+"."+str(splitNum)+".bed",'w')
    bed_file_list.append(output+"."+str(splitNum)+".bed")
    for intervals in bed_intervals:
        # Stop if we already reached our limit.
        # This can gappen since the size of chromosomes vary
        out.write(intervals['chr']+"\t"+intervals['start']+"\t"+intervals['end']+"\n")
        total+=intervals['length']
        if total > blockSize:
            out.close()
            splitNum += 1
            out=open(output+"."+str(splitNum)+".bed",'w')
            bed_file_list.append(output+"."+str(splitNum)+".bed")
            total = 0
    out.close()
        

    return bed_file_list
