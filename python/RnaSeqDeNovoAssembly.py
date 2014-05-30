#!/usr/bin/env python

import argparse

from Pipeline import *
from Trimmomatic import *

samples = ["sampleA", "sampleB"]
readsets = ["readsetA_1", "readsetA_2", "readsetB_1", "readsetB_2"]

def trim(readset):
    return trimmomatic(readset, readset + ".trim")

def normalization(readset):
    return normalize(readset + ".trim", readset + ".trim.normalized")

def align():
    return trinity([readset + ".trim.normalized" for readset in readsets])

def blast():
    return blastx("Trinity.fasta")

def abundance(sample):
    return rsem("Trinity.fasta", sample)

def annotate():
    return trinotate("Trinity.fasta")

def deliverable():
    return nozzle(["Trinity.fasta", "Trinity_stats.csv", "trinotate.tsv"] + ["rsem_" + sample + ".fpkm" for sample in samples], "report")

step_dict_map = [
    {"name": trim, "loop": readsets},
    {"name": normalization, "loop": readsets},
    {"name": align, "loop": GLOBAL},
    {"name": blast, "loop": GLOBAL},
    {"name": abundance, "loop": samples},
    {"name": annotate, "loop": GLOBAL},
    {"name": deliverable, "loop": GLOBAL}
]

pipeline = Pipeline([], step_dict_map, "1-7")
pipeline.parser.add_argument("-d", "--design", help="design file", type=file)
args = pipeline.parser.parse_args()
pipeline.show()
