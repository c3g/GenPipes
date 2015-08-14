#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
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
import math
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.sequence_dictionary import *

from bfx import blast
from bfx import blat
from bfx import bvatools
from bfx import bwa
from bfx import exonerate
from bfx import gatk
from bfx import igvtools
from bfx import metrics
from bfx import picard
from bfx import ray
from bfx import samtools
from bfx import tools
from pipelines import common


log = logging.getLogger(__name__)

class Puure(common.Illumina):

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#some fonctions used by "insert steps"
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    def get_job_sam_to_fastq(self, file_fastq):
        job=Job(
          output_files=[file_fastq],
          command="awk ' BEGIN { FS=\"\t\" } { print \"@\"$1; print $10; print \"+\";  print $11}' | gzip -c > " + file_fastq
        ) 
        return job

    def get_job_max_insert_size(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
        job=Job(
           input_files=[alignment_file_prefix + "all.metrics.insert_size_metrics"],
           command="maxInsertSize=$(awk 'NR==8 {print ($1+($2*3))}' "+ alignment_file_prefix + "all.metrics.insert_size_metrics)"
        )
        return job

    def get_job_min_insert_size(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
        job=Job(
           input_files=[alignment_file_prefix + "all.metrics.insert_size_metrics"],
           command="minInsertSize=$(awk 'NR==8 {print ($1-($2*3))}' "+ alignment_file_prefix + "all.metrics.insert_size_metrics)"
        )
        return job

    def get_job_mean_cov(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
        job=Job(
           input_files=[alignment_file_prefix + "all.coverage.sample_summary"],
           command="meanCov=$( grep " + sample.name + " " + alignment_file_prefix + "all.coverage.sample_summary | cut -f3)"
        )
        return job

    def get_job_mean_read_length(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
        job=Job(
           input_files=[alignment_file_prefix + "all.metrics.alignment_summary_metrics"],
           command="meanReadLen=$(awk 'NR==10 {print $16}' " + alignment_file_prefix + "all.metrics.alignment_summary_metrics)"
        )
        return job

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
#"alignment steps" if puure is running with fastq  (it's a copy/paste of DnaSeq without fixemate, recal and a subset of metrics)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
        return self._sequence_dictionary

    def bwa_mem_picard_sort_sam(self):
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name + ".sorted.bam")
            
            if readset.run_type == "PAIRED_END":
                fastq1 = trim_file_prefix + "pair1.fastq.gz"
                fastq2 = trim_file_prefix + "pair2.fastq.gz"
            elif readset.run_type == "SINGLE_END":
                fastq1 = trim_file_prefix + "single.fastq.gz"
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")
            
            job = concat_jobs([
                Job(command="mkdir -p " + alignment_directory),
                pipe_jobs([
                    bwa.mem(
                        fastq1,
                        fastq2,
                        read_group="'@RG" + \
                            "\tID:" + readset.name + \
                            "\tSM:" + readset.sample.name + \
                            "\tLB:" + readset.library + \
                            "\tPU:run" + readset.run + "_" + readset.lane + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        readset_bam,
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam." + readset.name)
            
            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM, along with its index.
            if len(readset.sample.readsets) == 1:
                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                sample_bam = os.path.join(alignment_directory, readset.sample.name + ".sorted.bam")
                sample_index = re.sub("\.bam$", ".bai", sample_bam)
                job = concat_jobs([
                    job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + os.path.relpath(readset_bam, os.path.dirname(sample_bam)) + " " + sample_bam),
                    Job([readset_bam], [sample_index], command="ln -s -f " + os.path.relpath(readset_index, os.path.dirname(sample_index)) + " " + sample_index)
                ], name=job.name)
            
            jobs.append(job)
        return jobs

    def picard_merge_sam_files(self):
        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name + ".sorted.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")
                
                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def gatk_indel_realigner(self):
        jobs = []
        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")
        
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            
            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + realign_directory),
                    gatk.realigner_target_creator(input, realign_intervals),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals),
                    # Create sample realign symlink since no merging is required
                    Job([output_bam], [sample_output_bam], command="ln -s -f " + os.path.relpath(output_bam, os.path.dirname(sample_output_bam)) + " " + sample_output_bam)
                ], name="gatk_indel_realigner." + sample.name))
            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job = [sequence['name'] for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                # Create one separate job for each of the first sequences
                for sequence in unique_sequences_per_job:
                    realign_prefix = os.path.join(realign_directory, sequence)
                    realign_intervals = realign_prefix + ".intervals"
                    output_bam = realign_prefix + ".bam"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + realign_directory),
                        gatk.realigner_target_creator(input, realign_intervals, intervals=[sequence]),
                        gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=[sequence])
                    ], name="gatk_indel_realigner." + sample.name + "." + sequence))
                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + realign_directory),
                    gatk.realigner_target_creator(input, realign_intervals, exclude_intervals=unique_sequences_per_job),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=["unmapped"], exclude_intervals=unique_sequences_per_job)
                ], name="gatk_indel_realigner." + sample.name + ".others"))
        return jobs

    def merge_realigned(self):
        jobs = []
        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")
            
            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                realigned_bams = [os.path.join(realign_directory, sequence['name'] + ".bam") for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                realigned_bams.append(os.path.join(realign_directory, "others.bam"))
                
                job = picard.merge_sam_files(realigned_bams, merged_realigned_bam)
                job.name = "merge_realigned." + sample.name
                jobs.append(job)
        return jobs

    def picard_mark_duplicates(self):
        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "realigned.qsorted.bam"
            output = alignment_file_prefix + "sorted.dup.bam"
            metrics_file = alignment_file_prefix + "sorted.dup.metrics"
            
            job = picard.mark_duplicates([input], output, metrics_file)
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def metrics(self):
        jobs = []
        for sample in self.samples:
            file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
            input = file_prefix + "bam"

            job = picard.collect_multiple_metrics(input, file_prefix + "all.metrics")
            job.name = "picard_collect_multiple_metrics." + sample.name
            jobs.append(job)

            # Compute genome coverage
            job = gatk.depth_of_coverage(input, file_prefix + "all.coverage")
            job.name = "gatk_depth_of_coverage.genome." + sample.name
            jobs.append(job)

            job = igvtools.compute_tdf(input, input + ".tdf")
            job.name = "igvtools_compute_tdf." + sample.name
            jobs.append(job)
        return jobs
  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
#"insert steps" when bam files and metrics are done
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    def extract_sclip(self):
        jobs = []
        
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            sclip_directory = os.path.join("sclip", sample.name)
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name)
            
            job = concat_jobs([
                Job(command="if [ ! -d " + sclip_directory + " ]; then mkdir -p " + sclip_directory + "; fi"),
                self.get_job_max_insert_size(sample),
                bvatools.extract_sclip(alignment_file_prefix + "sorted.dup.bam", sclip_file_prefix, "$maxInsertSize"),
                samtools.index(sclip_file_prefix + ".sc.bam"),
                samtools.index(sclip_file_prefix + ".scOthers.bam"),
                igvtools.compute_tdf(sclip_file_prefix + ".sc.bam", sclip_file_prefix + ".sc.tdf")
            ], name="extract_sclip_" + sample.name)
            
            jobs.append(job)
        return jobs

    def extract_bam_unmap(self):
        jobs = []
        
        for sample in self.samples:
            sclip_directory = os.path.join("sclip", sample.name)
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="if [ ! -d " + extract_directory + " ]; then mkdir -p " + extract_directory + "; fi")
            ## extract Orphan
            job = concat_jobs([
            jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "scOthers.bam", extract_file_prefix + "ORPHAN.bam", "-b -h -f 12 -F 256"),
                 samtools.sort(extract_file_prefix + "ORPHAN.bam", extract_file_prefix + "ORPHAN.sName", True)
                ])    
            ], name="extract_bam_ORPHAN_" + sample.name)
            
            jobs.append(job)
            
            ## extract OEA close to sclip
            job = concat_jobs([
            jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAUNMAP.1.bam", "-b -h -f 68 -F 264"),
                 samtools.sort(extract_file_prefix + "OEAUNMAP.1.bam", extract_file_prefix + "OEAUNMAP.1.sName", True)
                ])    
            ], name="extract_bam_OEAUNMAP1_" + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
            jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAUNMAP.2.bam", "-b -h -f 132 -F 264"),
                 samtools.sort(extract_file_prefix + "OEAUNMAP.2.bam", extract_file_prefix + "OEAUNMAP.2.sName", True)
                ])    
            ], name="extract_bam_OEAUNMAP2_" + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
            jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAMAP.bam", "-b -h -f 8 -F 1284"),
                 samtools.sort(extract_file_prefix + "OEAMAP.bam", extract_file_prefix + "OEAMAP.sName", True)
                ])    
            ], name="extract_bam_OEAMAP_" + sample.name)
            
            jobs.append(job)
            
        return jobs

    def extract_fastq_orphan(self):
        jobs = []

        for sample in self.samples:
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="if [ ! -d " + extract_directory + " ]; then mkdir -p " + extract_directory + "; fi")
            
            ## create fastq of ORPHAN
            job = picard.sam_to_fastq(extract_file_prefix + "ORPHAN.sName.bam", extract_file_prefix + "ORPHAN.1.fastq.gz", extract_file_prefix + "ORPHAN.2.fastq.gz")
            job.name = "extract_fastq_ORPHAN_" + sample.name
            jobs.append(job)
            
        return jobs        

    def extract_fastq_oea_sclip(self):
        jobs = []

        for sample in self.samples:
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="if [ ! -d " + extract_directory + " ]; then mkdir -p " + extract_directory + "; fi")
            
            ## create fastq of OEA close to sclip
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAUNMAP.1.sName.bam"),
                 self.get_job_sam_to_fastq(extract_file_prefix + "OEAUNMAP.1.fastq.gz")
            ], name="extract_fastq_OEA1_" + sample.name)
            jobs.append(job)
            
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAUNMAP.2.sName.bam"),
                 self.get_job_sam_to_fastq(extract_file_prefix + "OEAUNMAP.2.fastq.gz")
            ], name="extract_fastq_OEA2_" + sample.name)
            jobs.append(job)
            
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAMAP.bam", options="-f 64 -q " + config.param('DEFAULT', 'min_mapping_quality')),
                 self.get_job_sam_to_fastq(extract_file_prefix + "OEAMAP.1.fastq.gz")
            ], name="extract_fastq_OEAMAP1_" + sample.name)
            jobs.append(job)
            
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAMAP.bam", options="-f 128 -q " + config.param('DEFAULT', 'min_mapping_quality')),
                 self.get_job_sam_to_fastq(extract_file_prefix + "OEAMAP.2.fastq.gz")
            ], name="extract_fastq_OEAMAP2_" + sample.name)
            jobs.append(job)
            
            ## equal fastq file beetween OEAMAP and OEAUNMAP (due to the filter of mapq)
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "OEAMAP.2.fastq.gz", extract_file_prefix + "OEAUNMAP.1.fastq.gz", extract_file_prefix + "OEAUNMAP.1.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "OEAUNMAP.1.fastq.gz")
            ], name="equal_fastq_OEA1_" + sample.name)
            jobs.append(job)              
            
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "OEAMAP.1.fastq.gz", extract_file_prefix + "OEAUNMAP.2.fastq.gz", extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "OEAUNMAP.2.fastq.gz")
            ], name="equal_fastq_OEA2_" + sample.name)
            jobs.append(job)
            
            ## create fastq sclip
            jobMkdir = Job(command="if [ ! -d " + extract_directory + " ]; then mkdir -p " + extract_directory + "; fi")
            jobFastq = Job(
                input_files=[sclip_file_prefix+"scSequences.txt"],
                output_files=[extract_file_prefix+"sclip.1.fastq.gz"],
                command="awk 'NR>1 {if ($3==\"+\") { print \"@\"$4; print $5 ;print \"+\";  print $6}}' " + 
                        sclip_file_prefix + "scSequences.txt " +
                        "| gzip -c > " + extract_file_prefix+"sclip.1.fastq.gz",
                name="fastq_sclip1_" + sample.name
            )
            job = concat_jobs([
                jobMkdir,
                jobFastq
            ], name="fastq_sclip1_" + sample.name)
            jobs.append(job)
            
            jobFastq = Job(
                input_files=[sclip_file_prefix+"scSequences.txt"],
                output_files=[extract_file_prefix+"sclip.2.fastq.gz"],
                command="awk 'NR>1 {if ($3==\"-\") { print \"@\"$4; print $5 ;print \"+\";  print $6}}' " + 
                        sclip_file_prefix + "scSequences.txt " +
                        "| gzip -c > " + extract_file_prefix+"sclip.2.fastq.gz",
                name="fastq_sclip2_" + sample.name
            )
            job = concat_jobs([
                jobMkdir,
                jobFastq
            ], name="fastq_sclip2_" + sample.name)
            jobs.append(job)
            
        return jobs        

    def assembly_of_unmap(self):
        jobs = []
        
        for sample in self.samples:
            ray_directory1 = os.path.join("scaffolds", sample.name, "ray")
            ray_directory2 = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobRM = Job(command="if [ -d " + ray_directory2 + " ]; then rm -r " + ray_directory2 + "; fi")
            jobMkdir = Job(command="if [ ! -d " + ray_directory1 + " ]; then mkdir -p " + ray_directory1 + "; fi")
            
            jobRay = ray.ray(
                    ray_directory2, 
                    [extract_file_prefix + "ORPHAN.1.fastq.gz"], 
                    [extract_file_prefix + "ORPHAN.2.fastq.gz"], 
                    [
                     extract_file_prefix + "OEAUNMAP.1.equal.fastq.gz", 
                     extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz", 
                     extract_file_prefix + "sclip.1.fastq.gz", 
                     extract_file_prefix + "sclip.2.fastq.gz"
                    ]
            )
            jobRay.name="ray_" + sample.name
                        
            jobFormat = concat_jobs([
                Job(
                     input_files=[os.path.join(ray_directory2, "Scaffolds.fasta")],
                     command="sed -i '/^$/d' " + os.path.join(ray_directory2, "Scaffolds.fasta") + " && sed -i 's/scaffold-//g' " + os.path.join(ray_directory2, "Scaffolds.fasta")
                ),
                tools.py_addLengthRay(
                     os.path.join(ray_directory2, "Scaffolds.fasta"), 
                     os.path.join(ray_directory2, "ScaffoldLengths.txt"), 
                     os.path.join(ray_directory2, "Scaffolds.fasta.length")
                ),
                Job(
                     input_files=[
                        os.path.join(ray_directory2, "Scaffolds.fasta.length"), 
                        os.path.join(ray_directory2, "Scaffolds.fasta")
                     ],
                     command="rm -r " + os.path.join(ray_directory2, "Scaffolds.fasta") + 
                             " && mv " + os.path.join(ray_directory2, "Scaffolds.fasta.length") + " " + os.path.join(ray_directory2, "Scaffolds.fasta")
                ),
                exonerate.fastareformat(
                  os.path.join(ray_directory2, "Scaffolds.fasta"), 
                  os.path.join(ray_directory2, "Scaffolds2.fasta")
                ),
                Job(
                     input_files=[os.path.join(ray_directory2, "Scaffolds2.fasta")],
                     output_files=[os.path.join(ray_directory2, "Scaffolds.fasta")],
                     command="rm -r " + os.path.join(ray_directory2, "Scaffolds.fasta") + 
                             " && mv " + os.path.join(ray_directory2, "Scaffolds2.fasta")+ " " + os.path.join(ray_directory2, "Scaffolds.fasta")
                )
            ], name="format_" + sample.name)
            
            jobIndex = concat_jobs([
                samtools.faidx(os.path.join(ray_directory2, "Scaffolds.fasta")),
                bwa.index(os.path.join(ray_directory2, "Scaffolds.fasta"))
            ], name="index_" + sample.name)
            
            job = concat_jobs([
               jobRM,
               jobMkdir,
               jobRay,
               jobFormat,
               jobIndex
            ], name="rayAssembly_formatOutput_and_indexFiles_" + sample.name)
            jobs.append(job)
        
        return jobs

    def map_on_scaffolds(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            scaffolds_file = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "Scaffolds.fasta")
            
            #map Orphan read
            job = concat_jobs([
                Job(command="if [ ! -d " + cov_directory + " ]; then mkdir -p " + cov_directory + "; fi"),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "ORPHAN.1.fastq.gz",
                        extract_file_prefix + "ORPHAN.2.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_orphan" \
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:orphan" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=scaffolds_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        os.path.join(cov_directory, "ORPHAN.bam"),
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_ORPHAN_" + sample.name)
            jobs.append(job)
            
            #map OEA read
            job = concat_jobs([
                Job(command="if [ ! -d " + cov_directory + " ]; then mkdir -p " + cov_directory + "; fi"),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "OEAUNMAP.1.equal.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_scoea1"\
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:scoea1" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=scaffolds_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        os.path.join(cov_directory, "OEAUNMAP.1.bam"),
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_OEA1_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                Job(command="if [ ! -d " + cov_directory + " ]; then mkdir -p " + cov_directory + "; fi"),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_scoea2" \
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:scoea2" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=scaffolds_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        os.path.join(cov_directory, "OEAUNMAP.2.bam"),
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_OEA2_" + sample.name)
            jobs.append(job)
            
            #map sclip read
            job = concat_jobs([
                Job(command="if [ ! -d " + cov_directory + " ]; then mkdir -p " + cov_directory + "; fi"),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "sclip.1.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_sclip1" \
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:sclip1" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=scaffolds_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        os.path.join(cov_directory, "sclip.1.bam"),
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_sclip1_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                Job(command="if [ ! -d " + cov_directory + " ]; then mkdir -p " + cov_directory + "; fi"),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "sclip.2.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_sclip2" \
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:sclip2" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'",
                        ref=scaffolds_file
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        os.path.join(cov_directory, "sclip.2.bam"),
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_sclip2_" + sample.name)
            jobs.append(job)
        
        return jobs

    def merge_and_cov_scaffolds(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            
            job = picard.merge_sam_files(
                [
                 os.path.join(cov_directory, "sclip.1.bam"),
                 os.path.join(cov_directory, "sclip.2.bam"),
                 os.path.join(cov_directory, "OEAUNMAP.1.bam"),
                 os.path.join(cov_directory, "OEAUNMAP.2.bam"),
                 os.path.join(cov_directory, "ORPHAN.bam"),
                ],
                os.path.join(cov_directory, "readunmap.bam")
            )
            job.name = "covSca_merge_" + sample.name
            jobs.append(job)
            
            job = bvatools.depth_of_coverage(
                os.path.join(cov_directory, "readunmap.bam"), 
                os.path.join(cov_directory, "readunmap.cov.txt"), 
                [], 
                os.path.join(ray_directory, "Scaffolds.fasta"),
                "--gc --ommitN --minMappingQuality " + config.param('DEFAULT', 'min_mapping_quality') + " --threads " + config.param('merge_and_cov_scaffolds', 'threads')
            )
            job.name = "covSca_" + sample.name
            jobs.append(job)
        
        return jobs

    def blast_scaffolds_on_nt(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            
            job = blast.blastn_on_db(
                "nt", 
                os.path.join(ray_directory, "Scaffolds.fasta"), 
                os.path.join(ray_directory, "Scaffolds.fasta.blastn.xml"),
                "-evalue 1e-10 -max_target_seqs 20 -outfmt 5"
            )
            job.name = "blastn_sca_on_nt_" + sample.name
            jobs.append(job)
            
            job = tools.py_blastMatchSca(
                os.path.join(ray_directory, "Scaffolds"), 
                os.path.join(ray_directory, "Scaffolds.fasta.blastn.xml"), 
                os.path.join(ray_directory, "Scaffolds.fasta.blast")
            )
            job.name = "blastMatch_sca_" + sample.name
            jobs.append(job)
        
        return jobs

    def blat_scaffolds_on_ref(self):
        jobs = []
        
        for sample in self.samples:
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            
            job = concat_jobs([
            blat.blat_dna_vs_dna(
                   config.param('DEFAULT', 'genome_fasta'), 
                   os.path.join(ray_directory, "Scaffolds.fasta"), 
                   os.path.join(ray_directory, "Scaffolds.fasta.refGenome.psl")
                ),
                Job(
                 input_files=[os.path.join(ray_directory, "Scaffolds.fasta.refGenome.psl")],
                 output_files=[os.path.join(ray_directory, "Scaffolds.fasta.refGenome.noHeader.psl")],
                 command="awk 'NR>5 {print $0} ' " + os.path.join(ray_directory, "Scaffolds.fasta.refGenome.psl") + " > " + os.path.join(ray_directory, "Scaffolds.fasta.refGenome.noHeader.psl")
                )
            ], name = "blat_sca_on_ref_" + sample.name)
            jobs.append(job)
        
        return jobs

    def find_insertion(self, type_insert):
        jobs = []
        
        for sample in self.samples:
            main_directory = "./"
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            insert_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "insert" + type_insert)
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="if [ ! -d " + insert_directory + " ]; then mkdir -p " + insert_directory + "; fi")
            jobMaxInsert = self.get_job_max_insert_size(sample)
            jobMinInsert = self.get_job_min_insert_size(sample)
            jobMeanCov = self.get_job_mean_cov(sample)
            jobMeanReadLength = self.get_job_mean_read_length(sample)
            
            job = concat_jobs([
                jobMkdir,
                jobMinInsert,
                tools.r_select_scaffolds(
                  [
                   os.path.join(cov_directory, "readunmap.cov.txt"), 
                   os.path.join(ray_directory, "Scaffolds.fasta.refGenome.noHeader.psl"), 
                   os.path.join(ray_directory, "Scaffolds.fasta.blast")
                  ],
                  [
                   os.path.join(insert_directory, "scaffolds.tab"), 
                   os.path.join(insert_directory, "scaffolds.toDelete.tab")
                  ],
                  main_directory,
                  config.param('ray', 'kmer'),
                  sample.name,
                  type_insert,
                  "$minInsertSize"
                )
            ], name="analyse_scaffolds_" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                tools.r_find_cluster(
                  [
                   os.path.join(insert_directory, "scaffolds.tab"), 
                   extract_file_prefix + "OEAMAP.bam", 
                   os.path.join(cov_directory, "OEAUNMAP.1.bam"), 
                   os.path.join(cov_directory, "OEAUNMAP.2.bam")
                  ],
                  [
                   os.path.join(insert_directory, "cluster.OEA.tab"), 
                   os.path.join(insert_directory, "cluster.OEA.fusion.tab")
                  ],
                  main_directory,
                  config.param('ray', 'kmer'),
                  "OEA",
                  sample.name,
                  type_insert,
                  "$maxInsertSize",
                  config.param('DEFAULT', 'min_mapping_quality')
                )
            ], name="analyse_cluster_OEA_" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                tools.r_find_cluster(
                  [
                   os.path.join(insert_directory, "scaffolds.tab"), 
                   extract_file_prefix + "OEAMAP.bam", 
                   os.path.join(cov_directory, "OEAUNMAP.1.bam"), 
                   os.path.join(cov_directory, "OEAUNMAP.2.bam")
                  ],
                  [
                   os.path.join(insert_directory, "cluster.OEA.tab"), 
                   os.path.join(insert_directory, "cluster.OEA.fusion.tab")
                  ],
                  main_directory,
                  config.param('ray', 'kmer'),
                  "Sclip",
                  sample.name,
                  type_insert,
                  "$maxInsertSize",
                  config.param('DEFAULT', 'min_mapping_quality')
                )
            ], name="analyse_cluster_sclip_" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                jobMeanCov,
                tools.r_find_insert(
                    [
                     os.path.join(insert_directory, "scaffolds.tab"),
                     extract_file_prefix + "OEAMAP.bam", 
                     os.path.join(cov_directory, "OEAUNMAP.1.bam"),
                     os.path.join(cov_directory, "OEAUNMAP.2.bam")
                    ],
                    [
                     os.path.join(insert_directory, "cluster.OEA.tab"), 
                     os.path.join(insert_directory, "cluster.OEA.fusion.tab")
                    ],
                    main_directory,
                    config.param('ray', 'kmer'),
                    sample.name,
                    type_insert,
                    "$meanCov",
                    "$maxInsertSize",
                    config.param('DEFAULT', 'min_overlap_for_cluster'),
                    config.param('DEFAULT', 'genome_mappability_bed_indexed')
                )
            ], name="analyse_find_insert_" + type_insert + "_" + sample.name)
            jobs.append(job)
            
        return jobs

    def find_virus(self):
        return self.find_insertion("Virus")

    def find_human(self):
        return self.find_insertion("Human")

    def find_other(self):
        return self.find_insertion("Other")

    def find_none(self):
        return self.find_insertion("None")

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.bwa_mem_picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.picard_mark_duplicates,
            self.metrics,
            self.extract_sclip,
            self.extract_bam_unmap,
            self.extract_fastq_orphan,
            self.extract_fastq_oea_sclip,
            self.assembly_of_unmap,
            self.map_on_scaffolds,
            self.merge_and_cov_scaffolds,
            self.blast_scaffolds_on_nt,
            self.blat_scaffolds_on_ref,
            self.find_virus,
            self.find_human,
            self.find_other,
            self.find_none,
        ]

if __name__ == "__main__": 
    # PUUre has not been completely tested
    raise NotImplementedError
    Puure()
