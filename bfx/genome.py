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
import re

# MUGQIC Modules
from core.config import *
from core.job import *



def chr_names(genome_dict_file):
    """
    extracts chromosome ids from the genome_dictionary file.
    Returns array of chr names in the order they are in genome_dictionary.
    """

    chrs = []
    genome_dict = open(genome_dict_file, 'r')

    for line in genome_dict:
        mychr = line.split()[1]
        if mychr.startswith("SN:"):
            mychr = re.sub("^SN:", "", mychr)
            chrs.append(mychr)
    
    genome_dict.close()

    return chrs


def chr_sizes(genome_dict_file):
    """
    extracts chromosome ids with their sizes from the genome_dictionary file.
    Returns dictionary of chr:size in the order they are in genome_dictionary.
    """

    chrs = {}
    genome_dict = open(genome_dict_file, 'r')

    for line in genome_dict:
        mychr = line.split()[1]
        if mychr.startswith("SN:"):
            mychr = re.sub("^SN:", "", mychr)
        else:
            continue
        size = line.split()[2]
        if size.startswith("LN:"):
            size = int(re.sub("^LN:", "", size))

        chrs[mychr] = size
    
    genome_dict.close()

    return chrs

def chr_names_conv(genome_dict_file, chrX=True, chrY=False, chrM=False, otherPatterns=None):
    """
    extracts chromosome ids from the genome_dictionary file.
    Removes all non conventional chromosomes containing "random", "hap", "chrUn", "GL", "NT_".
    Can also remove chrX, chrM, chrY if params set to false.
    To remove other chr patterns list them as an array in otherChr.
    Returns array of chr names in the order they are in genome_dictionary.
    """

    ## extract chrs using chr_names then filter unwanted ones:
    chrs = chr_names(genome_dict_file)
    chrs_conv = []
    remove = ["hap", "random", "chrUn", "EBV", "GL", "NT_"]

    if not chrY:
        remove.append("chrY")
        remove.append("Y")

    if not chrM:
        remove.append("chrM")
        remove.append("MT")

    if not chrX:
        remove.append("chrX")
        remove.append("X")

    if otherPatterns is not None:
        remove.extend(otherPatterns)

    for mychr in chrs:
        if not any(x in mychr for x in remove):
            chrs_conv.append(mychr)

    return chrs_conv


def genome_size(genome_dict_file):
    """
    Estimates genome size by adding all chrs and contig sizes in genome_dict_file.
    """
    chrs = chr_sizes(genome_dict_file)
    return sum(chrs.values())


def genome_size_conv(genome_dict_file):
    """
    Estimates genome size by adding conventional chr sizes in genome_dict_file. It ignores contigs.
    """
    chrs = chr_sizes(genome_dict_file)
    chrs_conv = chr_names_conv(genome_dict_file)

    for key, value in list(chrs.items()):
        if key not in chrs_conv:
            del chrs[key]

    return sum(chrs.values())
