#!/usr/bin/env python

from Job import *

def trimmomatic(input1, output1):
    job = Job(None, [input1], [output1])
    job.command = "java -jar Trimmomatic.jar " + input1 + " > " + output1
    return job

def normalize(input1, output1):
    job = Job(None, [input1], [output1])
    job.command = "normalize " + input1 + " > " + output1
    return job

def trinity(normalized_readsets):
    job = Job(None, normalized_readsets, ["Trinity.fasta", "Trinity_stats.csv"])
    job.command = "trinity " + " ".join(normalized_readsets) + " > Trinity.fasta"
    return job

def blastx(input1):
    job = Job(None, [input1], ["blastx_nr.tsv"])
    job.command = "blastx " + input1 + " > blastx_nr.tsv"
    return job

def rsem(reference, fasta):
    job = Job(None, [reference, fasta], ["rsem_" + fasta + ".fpkm"])
    job.command = "rsem -db " + reference + " -input " + fasta
    return job

def trinotate(input1):
    job = Job(None, [input1], ["trinotate.tsv"])
    job.command = "trinotate " + input1 + " > trinotate.tsv"
    return job

def nozzle(input_files, output1):
    job = Job(None, input_files, [output1])
    job.command = "nozzle " + " ".join(input_files) + " > " + output1
    return job
