#!/usr/bin/env python

from Job import *

def trimmomatic(input1, output1, output2):
    job = Job(None, [input1], [output1, output2])
    job.command = "java -jar Trimmomatic.jar " + input1 + " > " + output1
    return job

def trinity(input1):
    job = Job(None, [input1], ["Trinity.fasta"])
    job.command = "trinity " + input1 + " > Trinity.fasta"
    return job

def blastx(input1):
    job = Job(None, [input1], ["blastx_nr.tsv"])
    job.command = "blastx " + input1 + " > blastx_nr.tsv"
    return job

def nozzle(input1, output1):
    job = Job(None, [input1], [output1])
    job.command = "nozzle " + input1 + " > " + output1
    return job
