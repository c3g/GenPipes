#!/usr/bin/env python

# Python Standard Modules
import logging
import math
import os
import re
import sys

# Append mugqic_pipeline directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bio.readset import *
from bio.sequence_dictionary import *

from bio import bvatools
from bio import bwa
from bio import gatk
from bio import gq_seq_utils
from bio import igvtools
from bio import metrics
from bio import picard
from bio import samtools
from bio import snpeff
from bio import tools
from bio import vcftools
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class Puure(dnaseq.DnaSeq):
  
    @property  
    def extract_sclip(self):
        jobs = []
        
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            sclip_directory = os.path.join("sclip", sample.name)
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            
            Job(command="insertSize=\$(awk 'NR==8 {print (\$1+(\$2*3))}' "+ alignment_file_prefix + ".sorted.dup.all.metrics.insert_size_metrics)")
            
            job = bvatools.extract_sclip(alignment_file_prefix + "sorted.dup.bam", sclip_file_prefix)
            job.name = "extract_sclip." + sample.name
            jobs.append(job)
            
            job = concat_jobs([Job(command="mkdir -p " + sclip_directory),
                job
            ], name="extract_sclip." + sample.name)
            
            job = concat_jobs([job,
                samtools.index(sclip_file_prefix + "sc.bam", sclip_file_prefix + "sc.bai"),
            ], name="extract_sclip." + sample.name)
            
            job = concat_jobs([job,
                samtools.index(sclip_file_prefix + "scOthers.bam", sclip_file_prefix + "scOthers.bai"),
            ], name="extract_sclip." + sample.name)
            
            job = concat_jobs([job,
                igvtools.compute_tdf(sclip_file_prefix + "sc.bam", sclip_file_prefix + "sc.tdf"),
            ], name="extract_sclip." + sample.name)
            
            jobs.append(job)
        return jobs
    
    def extract_bam_unmap(self):
        jobs = []

        for sample in self.samples:
            sclip_directory = os.path.join("sclip", sample.name)
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="mkdir -p " + extract_directory)

            ## extract Orphan
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "scOthers.bam", extract_file_prefix + "ORPHAN.bam", "-b -h -f 12 -F 256"),
                 samtools.sort(extract_file_prefix + "ORPHAN.bam", extract_file_prefix + "ORPHAN.sName", True)
                ])    
            ], name="extract_ORPHAN." + sample.name)
            
            jobs.append(job)
            
            ## extract OEA close to sclip
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "scOEAUNMAP.1.bam", "-b -h -f 68 -F 264"),
                 samtools.sort(extract_file_prefix + "scOEAUNMAP.1.bam", extract_file_prefix + "scOEAUNMAP.1.sName", True)
                ])    
            ], name="extract_scOEAUNMAP.1." + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "scOEAUNMAP.2.bam", "-b -h -f 132 -F 264"),
                 samtools.sort(extract_file_prefix + "scOEAUNMAP.2.bam", extract_file_prefix + "scOEAUNMAP.2.sName", True)
                ])    
            ], name="extract_scOEAUNMAP.2." + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "scOEAMAP.bam", "-b -h -f 8 -F 1284"),
                 samtools.sort(extract_file_prefix + "scOEAMAP.bam", extract_file_prefix + "scOEAMAP.sName", True)
                ])    
            ], name="extract_scOEAMAP." + sample.name)
            
            jobs.append(job)
            
        return jobs
    
    def extract_fastq_unmap(self):
        jobs = []

        for sample in self.samples:
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="mkdir -p " + extract_directory)

            ## create fastq of ORPHAN
            job = picard.sam_to_fastq(extract_file_prefix + "ORPHAN.sName.bam", extract_file_prefix + "ORPHAN.1.fastq.gz", extract_file_prefix + "ORPHAN.2.fastq.gz")
            job.name = "extract_fastq_ORPHAN." + sample.name
            jobs.append(job)
            
            ## create fastq of OEA close to sclip
            file_fastq=extract_file_prefix + "scOEAUNMAP.1.fastq"
            job1 = pipe_jobs([
                 samtools.view(extract_file_prefix + "scOEAUNMAP.1.sName.bam"),
                 tools.awk_samToFastq(output=file_fastq)    
            ], name="extract_fastq_SCOEA1" + sample.name)
            job2 = concat_jobs([
                 Job(input_files=[file_fastq], output_files=[file_fastq+".gz"], command="gzip -c " + file_fastq + " > " + file_fastq + ".gz"),
                 Job(command="rm " + file_fastq)
            ], name="extract_fastq_SCOEA1" + sample.name)
            job = concat_jobs([
                 job1,
                 job2
            ], name="extract_fastq_SCOEA1" + sample.name)
            jobs.append(job)
            
            file_fastq=extract_file_prefix + "scOEAUNMAP.2.fastq"
            job1 = pipe_jobs([
                 samtools.view(extract_file_prefix + "scOEAUNMAP.2.sName.bam"),
                 tools.awk_samToFastq(output=file_fastq)    
            ], name="extract_fastq_SCOEA2" + sample.name)
            job2 = concat_jobs([
                 Job(input_files=[file_fastq], output_files=[file_fastq+".gz"], command="gzip -c " + file_fastq + " > " + file_fastq + ".gz"),
                 Job(command="rm " + file_fastq)
            ], name="extract_fastq_SCOEA2" + sample.name)
            job = concat_jobs([
                 job1,
                 job2
            ], name="extract_fastq_SCOEA2" + sample.name)
            jobs.append(job)
            
            file_fastq=extract_file_prefix + "scOEAMAP.1.fastq"
            job1 = pipe_jobs([
                 samtools.view(extract_file_prefix + "scOEAMAP.bam", "-f 64 -q " + config.param('DEFAULT', 'minMappingQuality')),
                 tools.awk_samToFastq(output=file_fastq)    
            ], name="extract_fastq_SCOEAMAP1" + sample.name)
            job2 = concat_jobs([
                 Job(input_files=[file_fastq], output_files=[file_fastq+".gz"], command="gzip -c " + file_fastq + " > " + file_fastq + ".gz"),
                 Job(command="rm " + file_fastq)
            ], name="extract_fastq_SCOEAMAP1" + sample.name)
            job = concat_jobs([
                 job1,
                 job2
            ], name="extract_fastq_SCOEAMAP1" + sample.name)
            jobs.append(job) 
            
            file_fastq=extract_file_prefix + "scOEAMAP.2.fastq"
            job1 = pipe_jobs([
                 samtools.view(extract_file_prefix + "scOEAMAP.bam", "-f 128 -q " + config.param('DEFAULT', 'minMappingQuality')),
                 tools.awk_samToFastq(output=file_fastq)    
            ], name="extract_fastq_SCOEAMAP2" + sample.name)
            job2 = concat_jobs([
                 Job(input_files=[file_fastq], output_files=[file_fastq+".gz"], command="gzip -c " + file_fastq + " > " + file_fastq + ".gz"),
                 Job(command="rm " + file_fastq)
            ], name="extract_fastq_SCOEAMAP2" + sample.name)
            job = concat_jobs([
                 job1,
                 job2
            ], name="extract_fastq_SCOEAMAP2" + sample.name)
            jobs.append(job)              
          
            ## equal fastq file beetween OEAMAP and OEAUNMAP (due to the filter of mapq)
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "scOEAMAP.2.fastq.gz", extract_file_prefix + "scOEAUNMAP.1.fastq.gz", extract_file_prefix + "scOEAUNMAP.1.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "scOEAUNMAP.1.fastq.gz")
            ], name="equal_fastq_SCOEA1" + sample.name)
            jobs.append(job)              
          
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "scOEAMAP.1.fastq.gz", extract_file_prefix + "scOEAUNMAP.2.fastq.gz", extract_file_prefix + "scOEAUNMAP.2.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "scOEAUNMAP.2.fastq.gz")
            ], name="equal_fastq_SCOEA2" + sample.name)
            jobs.append(job)    
            
        return jobs        
     
    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.bwa_mem_picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.fix_mate_by_coordinate,
            self.picard_mark_duplicates,
            self.metrics,
            self.extract_sclip,
            self.extract_bam_unmap,
            self.extract_fastq_unmap,
        ]

Puure().submit_jobs()
