#!/usr/bin/env python

from Pipeline import *
from Trimmomatic import *

samples = ["sampleA", "sampleB"]
readsets = ["readsetA_1", "readsetA_2", "readsetB_1", "readsetB_2"]

def trim(readset):
    return trimmomatic(readset, readset + ".output", "")

def align(readset):
    return trinity(readset + ".output")

def blast(sample):
    return blastx("Trinity.fasta")

def deliverable(sample):
    return nozzle("stats.csv", "report")

steps = [
    {"name": trim, "loop": readsets},
    {"name": align, "loop": readsets, "parents": [trim]},
    {"name": blast, "loop": samples},
    {"name": deliverable, "loop": GLOBAL, "parents": [align, blast]}
]

#obj1 = Pipeline('rnaSeqDenovoAssembly', [], ["trim", "align", "blast", "abundance", "dge", "trinotate", "deliverable"], "1,3,5-7")
#obj1.show()
obj2 = Pipeline('rnaSeqDenovoAssembly2', [], steps, "1,2,3-4")
obj2.show()
print steps[2]
