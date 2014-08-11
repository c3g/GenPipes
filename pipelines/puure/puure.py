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

from bio import blast
from bio import blat
from bio import bvatools
from bio import bwa
from bio import exonerate
from bio import gatk
from bio import gq_seq_utils
from bio import igvtools
from bio import metrics
from bio import picard
from bio import ray
from bio import samtools
from bio import snpeff
from bio import tools
from bio import vcftools
from pipelines.dnaseq import dnaseq

log = logging.getLogger(__name__)

class Puure(dnaseq.DnaSeq):
  
    def get_job_max_insert_size(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
        job=Job(command="maxInsertSize=\$(awk 'NR==8 {print (\$1+(\$2*3))}' "+ alignment_file_prefix + ".sorted.dup.recall.all.metrics.insert_size_metrics)")
        return job

    def get_job_min_insert_size(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
        job=Job(command="minInsertSize=\$(awk 'NR==8 {print (\$1-(\$2*3))}' "+ alignment_file_prefix + ".sorted.dup.recall.all.metrics.insert_size_metrics)")
        return job

    def get_job_mean_cov(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
        job=Job(command="meanCov=\$( grep " + sample.name + " " + alignment_file_prefix + ".sorted.dup.recall.all.coverage.sample_summary | cut -f3")
        return job

    def get_job_mean_read_length(self, sample):
        alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
        job=Job(command="meanReadLen=\$(awk 'NR==10 {print \$16}' " + alignment_file_prefix + "sorted.dup.recall.all.metrics.alignment_summary_metrics")
        return job

    def extract_sclip(self):
        jobs = []
        
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            sclip_directory = os.path.join("sclip", sample.name)
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            
            job = concat_jobs([Job(command="mkdir -p " + sclip_directory),
                self.get_job_max_insert_size(sample)
            ], name="extract_sclip_" + sample.name)
            
            job = concat_jobs([job,
                bvatools.extract_sclip(alignment_file_prefix + "sorted.dup.recall.bam", sclip_file_prefix, "\$maxInsertSize")
            ], name="extract_sclip_" + sample.name)
            
            job = concat_jobs([job,
                samtools.index(sclip_file_prefix + "sc.bam", sclip_file_prefix + "sc.bai"),
            ], name="extract_sclip_" + sample.name)
            
            job = concat_jobs([job,
                samtools.index(sclip_file_prefix + "scOthers.bam", sclip_file_prefix + "scOthers.bai"),
            ], name="extract_sclip_" + sample.name)
            
            job = concat_jobs([job,
                igvtools.compute_tdf(sclip_file_prefix + "sc.bam", sclip_file_prefix + "sc.tdf"),
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
            
            jobMkdir = Job(command="mkdir -p " + extract_directory)
            ## extract Orphan
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "scOthers.bam", extract_file_prefix + "ORPHAN.bam", "-b -h -f 12 -F 256"),
                 samtools.sort(extract_file_prefix + "ORPHAN.bam", extract_file_prefix + "ORPHAN.sName", True)
                ])    
            ], name="extract_ORPHAN_" + sample.name)
            
            jobs.append(job)
            
            ## extract OEA close to sclip
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAUNMAP.1.bam", "-b -h -f 68 -F 264"),
                 samtools.sort(extract_file_prefix + "OEAUNMAP.1.bam", extract_file_prefix + "OEAUNMAP.1.sName", True)
                ])    
            ], name="extract_OEAUNMAP1_" + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAUNMAP.2.bam", "-b -h -f 132 -F 264"),
                 samtools.sort(extract_file_prefix + "OEAUNMAP.2.bam", extract_file_prefix + "OEAUNMAP.2.sName", True)
                ])    
            ], name="extract_OEAUNMAP2_" + sample.name)
            
            jobs.append(job)
            
            job = concat_jobs([
	        jobMkdir,
                concat_jobs([
                 samtools.view(sclip_file_prefix + "sc.bam", extract_file_prefix + "OEAMAP.bam", "-b -h -f 8 -F 1284"),
                 samtools.sort(extract_file_prefix + "OEAMAP.bam", extract_file_prefix + "OEAMAP.sName", True)
                ])    
            ], name="extract_OEAMAP_" + sample.name)
            
            jobs.append(job)
            
        return jobs

    def extract_fastq_unmap(self):
        jobs = []

        for sample in self.samples:
            extract_directory = os.path.join("extract", sample.name)
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            sclip_file_prefix = os.path.join("sclip", sample.name, sample.name + ".")
            
            jobMkdir = Job(command="mkdir -p " + extract_directory)
            
            ## create fastq of ORPHAN
            job = picard.sam_to_fastq(extract_file_prefix + "ORPHAN.sName.bam", extract_file_prefix + "ORPHAN.1.fastq.gz", extract_file_prefix + "ORPHAN.2.fastq.gz")
            job.name = "extract_fastq_ORPHAN_" + sample.name
            jobs.append(job)
            
            
            
            ## create fastq of OEA close to sclip
            file_fastq=extract_file_prefix + "OEAUNMAP.1.fastq.gz"
            jobAwk = Job([],[file_fastq])
            jobAwk.command("awk ' BEGIN { FS=\"\t\" } { print \"@\"\$1; print \$10; print \"+\";  print \$11}' | gzip -c > " + file_fastq)
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAUNMAP.1.sName.bam"),
                 jobAwk    
            ], name="extract_fastq_SCOEA1_" + sample.name)
            jobs.append(job)
            
            file_fastq=extract_file_prefix + "OEAUNMAP.2.fastq.gz"
            jobAwk = Job([],[file_fastq])
            jobAwk.command("awk ' BEGIN { FS=\"\t\" } { print \"@\"\$1; print \$10; print \"+\";  print \$11}' | gzip -c > " + file_fastq)
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAUNMAP.2.sName.bam"),
                 jobAwk    
            ], name="extract_fastq_SCOEA2_" + sample.name)
            jobs.append(job)
            
            file_fastq=extract_file_prefix + "OEAMAP.1.fastq.gz"
            jobAwk = Job([],[file_fastq])
            jobAwk.command("awk ' BEGIN { FS=\"\t\" } { print \"@\"\$1; print \$10; print \"+\";  print \$11}' | gzip -c > " + file_fastq)
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAMAP.bam", "-f 64 -q " + config.param('DEFAULT', 'min_mapping_quality')),
                 jobAwk    
            ], name="extract_fastq_SCOEAMAP1_" + sample.name)
            jobs.append(job)
            
            file_fastq=extract_file_prefix + "OEAMAP.2.fastq.gz"
            jobAwk = Job([],[file_fastq])
            jobAwk.command("awk ' BEGIN { FS=\"\t\" } { print \"@\"\$1; print \$10; print \"+\";  print \$11}' | gzip -c > " + file_fastq)
            job = pipe_jobs([
                 samtools.view(extract_file_prefix + "OEAMAP.bam", "-f 128 -q " + config.param('DEFAULT', 'min_mapping_quality')),
                 jobAwk    
            ], name="extract_fastq_SCOEAMAP2_" + sample.name)
            jobs.append(job)
            
            ## equal fastq file beetween OEAMAP and OEAUNMAP (due to the filter of mapq)
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "OEAMAP.2.fastq.gz", extract_file_prefix + "OEAUNMAP.1.fastq.gz", extract_file_prefix + "OEAUNMAP.1.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "OEAUNMAP.1.fastq.gz")
            ], name="equal_fastq_SCOEA1_" + sample.name)
            jobs.append(job)              
            
            job = concat_jobs([
                 tools.py_equalFastqFile(extract_file_prefix + "OEAMAP.1.fastq.gz", extract_file_prefix + "OEAUNMAP.2.fastq.gz", extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz"),
                 Job(command="rm " + extract_file_prefix + "OEAUNMAP.2.fastq.gz")
            ], name="equal_fastq_SCOEA2_" + sample.name)
            jobs.append(job)
            
            ## create fastq sclip
            job = Job([sclip_file_prefix+"scSequences.txt"],[extract_file_prefix+"sclip.1.fastq.gz"])
            job.name = "sclip1_fastq_" + sample.name
            job.command("awk 'NR>1 {if (\$3==\"+\") { print \"@\"\$4; print \$5 ;print \"+\";  print \$6}} '" + 
                       sclip_file_prefix + "scSequences.txt " + 
                       "| gzip -c > " + extract_file_prefix+"sclip.1.fastq.gz"
                       )
            jobs.append(job)
            
            job = Job([sclip_file_prefix+"scSequences.txt"],[extract_file_prefix+"sclip.2.fastq.gz"])
            job.name = "sclip2_fastq_" + sample.name
            job.command("awk 'NR>1 {if (\$3==\"-\") { print \"@\"\$4; print \$5 ;print \"+\";  print \$6}} '" + 
                       sclip_file_prefix + "scSequences.txt " + 
                       "| gzip -c > " + extract_file_prefix+"sclip.2.fastq.gz"
                       )
            jobs.append(job)
            
        return jobs        

    def assembly_of_unmap(self):
        jobs = []
        
        for sample in self.samples:
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobRM = Job(command="rm -r " + ray_directory)
            jobRay = ray.ray(
                ray_directory, 
                [extract_file_prefix + "ORPHAN.1.fastq.gz"], 
                [extract_file_prefix + "ORPHAN.2.fastq.gz"], 
                [extract_file_prefix + "OEAUNMAP.1.equal.fastq.gz", 
                 extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz", 
                 extract_file_prefix + "sclip.1.fastq.gz", 
                 extract_file_prefix + "sclip.2.fastq.gz"
                ]
            )         
            job = concat_jobs([jobRM,
                jobRay
            ], name="ray_" + sample.name)
            
            jobFormat = Job(command="sed -i '/^$/d' " + ray_directory + "Scaffolds.fasta " + "&& sed -i 's/scaffold-//g' " + ray_directory + "Scaffolds.fasta")
            jobFormat = concat_jobs([jobFormat,
                tools.py_addLengthRay(ray_directory + "Scaffolds.fasta", ray_directory + "ScaffoldLengths.txt", ray_directory + "Scaffolds.fasta.length")
            ], name="format_" + sample.name)
            jobFormat = concat_jobs([jobFormat,
                Job(command="rm -r " + ray_directory + "Scaffolds.fasta" + " && mv " + ray_directory + "Scaffolds.fasta.length " + ray_directory + "Scaffolds.fasta")
            ], name="format_" + sample.name)
            jobFormat = concat_jobs([jobFormat,
                exonerate.fastareformat(ray_directory + "Scaffolds.fasta", ray_directory + "Scaffolds2.fasta")
            ], name="format_" + sample.name)
            jobFormat = concat_jobs([jobFormat,
                Job(command="rm -r " + ray_directory + "Scaffolds.fasta" + " && mv " + ray_directory + "Scaffolds2.length " + ray_directory + "Scaffolds.fasta")
            ], name="format_" + sample.name)
            
            jobIndex = concat_jobs([
                samtools.faidx(ray_directory + "Scaffolds.fasta"),
                bwa.index(ray_directory + "Scaffolds.fasta")
            ], name="index_" + sample.name)
            
            job = concat_jobs([
                jobRay,
                concat_jobs([
                  jobFormat,
                  jobIndex
                ], name="formatOutput_and_index_" + sample.name)
            ], name="ray_formatOutput_and_index_" + sample.name)
            jobs.append(job)
        
        return jobs

    def map_on_scaffolds(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            #map Orphan read
            job = concat_jobs([
                Job(command="mkdir " + cov_directory),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "ORPHAN.1.fastq.gz",
                        extract_file_prefix + "ORPHAN.2.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_orphan"\
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:orphan" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        cov_directory + "ORPHAN.bam",
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_ORPHAN_" + sample.name)
            jobs.append(job)
            
            #map OEA read
            job = concat_jobs([
                Job(command="mkdir " + cov_directory),
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
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        cov_directory + "OEAUNMAP.1.bam",
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_OEA1_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                Job(command="mkdir " + cov_directory),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "OEAUNMAP.2.equal.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_scoea2"\
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:scoea2" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        cov_directory + "OEAUNMAP.2.bam",
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_OEA2_" + sample.name)
            jobs.append(job)
            
            #map sclip read
            job = concat_jobs([
                Job(command="mkdir " + cov_directory),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "sclip.1.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_sclip1"\
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:sclip1" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        cov_directory + "sclip.1.bam",
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam_sclip1_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                Job(command="mkdir " + cov_directory),
                pipe_jobs([
                    bwa.mem(
                        extract_file_prefix + "sclip.2.fastq.gz",
                        read_group="'@RG" + \
                            "\tID:" + sample.name + "_ray_sclip2"\
                            "\tSM:" + sample.name + \
                            "\tLB:" + sample.name + \
                            "\tPU:sclip2" + \
                            "\tCN:" + config.param('bwa_mem', 'sequencing_center') + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        cov_directory + "sclip.2.bam",
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
            
            job = picard.merge_sam_files([
              cov_directory + "sclip.1.bam",
              cov_directory + "sclip.2.bam",
              cov_directory + "OEAUNMAP.1.bam",
              cov_directory + "OEAUNMAP.2.bam",
              cov_directory + "ORPHAN.bam",
            ],cov_directory + "readunmap.bam")
            job.name = "covSca_merge_" + sample.name
            jobs.append(job)
            
            job = bvatools.depth_of_coverage(
                cov_directory + "readunmap.bam", 
                cov_directory + "readunmap.cov.txt", 
                [], 
                ray_directory + "Scaffolds.fasta",
                "--gc --ommitN --minMappingQuality " + config.param('DEFAULT', 'min_mapping_quality') + "--threads 5")
            job.name = "covSca_" + sample.name
            jobs.append(job)
        
        return jobs

    def blast_scaffolds_on_nt(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            
            job = blast.blast_on_db(
                "nt", 
                ray_directory + "Scaffolds.fasta", 
                ray_directory + "Scaffolds.fasta.blastn.xml",
                "-evalue 1e-10 -max_target_seqs 20 -outfmt 5"
            )
            job.name = "blastn_sca_on_nt_" + sample.name
            jobs.append(job)
            
            job = tools.py_blastMatchSca(
                ray_directory + "Scaffolds", 
                ray_directory + "Scaffolds.fasta.blastn.xml", 
                ray_directory + "Scaffolds.fasta.blast"
            )
            job.name = "blastMatch_sca_" + sample.name
            jobs.append(job)
        
        return jobs

    def blat_scaffolds_on_ref(self):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            
            job = blat.blat_dna_vs_dna(
                config.param('DEFAULT', 'genome_fasta'), 
                ray_directory + "Scaffolds.fasta", 
                ray_directory + "Scaffolds.fasta.refGenome.psl"
            )
            job.name = "blat_sca_on_ref_" + sample.name
            jobs.append(job)
        
        return jobs

    def find_insertion(self, type_insert):
        jobs = []
        
        for sample in self.samples:
            cov_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "cov")
            insert_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'), "insert" + type_insert)
            ray_directory = os.path.join("scaffolds", sample.name, "ray", "ray" + config.param('ray', 'kmer'))
            extract_file_prefix = os.path.join("extract", sample.name, sample.name + ".")
            
            jobMaxInsert = self.get_job_max_insert_size(sample)
            jobMinInsert = self.get_job_min_insert_size(sample)
            jobMeanCov = self.get_job_mean_cov(sample)
            jobMeanReadLength = self.get_job_mean_read_length(sample)
            
            job = concat_jobs([
                jobMinInsert,
                tools.r_select_scaffolds(
                  [cov_directory + "readunmap.cov.txt", ray_directory + "Scaffolds.fasta.refGenome.noHeader.psl", ray_directory + "Scaffolds.fasta.blast"],
                  [insert_directory + "scaffolds.tab", insert_directory + "scaffolds.toDelete.tab"],
                  sample.name,
                  type_insert,
                  "\$minInsertSize"
                )
            ], name="analyse_scaffolds_" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                tools.r_find_cluster(
                  [insert_directory + "scaffolds.tab", extract_file_prefix + "OEAMAP.bam", cov_directory + "OEAUNMAP.1.bam", cov_directory + "OEAUNMAP.2.bam"],
                  [insert_directory + "cluster.OEA.tab", insert_directory + "cluster.OEA.fusion.tab"],
                  "OEA",
                  sample.name,
                  type_insert,
                  "\$maxInsertSize",
                  config.param('DEFAULT', 'min_mapping_quality')
                )
            ], name="analyse_cluster_OEA" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                tools.r_find_cluster(
                  [insert_directory + "scaffolds.tab", extract_file_prefix + "OEAMAP.bam", cov_directory + "OEAUNMAP.1.bam", cov_directory + "OEAUNMAP.2.bam"],
                  [insert_directory + "cluster.OEA.tab", insert_directory + "cluster.OEA.fusion.tab"],
                  "sclip",
                  sample.name,
                  type_insert,
                  "\$maxInsertSize",
                  config.param('DEFAULT', 'min_mapping_quality')
                )
            ], name="analyse_cluster_sclip" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                concat_job([jobMeanCov,
                  tools.r_find_insert(
                    [insert_directory + "scaffolds.tab", extract_file_prefix + "OEAMAP.bam", cov_directory + "OEAUNMAP.1.bam", cov_directory + "OEAUNMAP.2.bam"],
                    [insert_directory + "cluster.OEA.tab", insert_directory + "cluster.OEA.fusion.tab"],
                    sample.name,
                    type_insert,
                    "\$meanCov",
                    "\$maxInsertSize",
                    config.param('DEFAULT', 'min_overlap_for_cluster'),
                    config.param('DEFAULT', 'genome_mappability_bed_indexed')
                  )
                ], name="analyse_find_insert_tmp" + type_insert + "_" + sample.name)
            ], name="analyse_find_insert" + type_insert + "_" + sample.name)
            jobs.append(job)
            
            job = concat_jobs([
                jobMaxInsert,
                concat_job([jobMeanCov,
                  concat_job([jobMeanReadLength,
                    tools.r_filter_insert(
                      [insert_directory + "scaffolds.tab", extract_file_prefix + "OEAMAP.bam", cov_directory + "OEAUNMAP.1.bam", cov_directory + "OEAUNMAP.2.bam"],
                      [insert_directory + "cluster.OEA.tab", insert_directory + "cluster.OEA.fusion.tab"],
                      sample.name,
                      type_insert,
                      "\$meanCov",
                      "\$maxInsertSize",
                      config.param('DEFAULT', 'num_starnd_for_insertion'),
                      config.param('DEFAULT', 'min_read_for_insertion'),
                      "\$meanReadLen"
                    )
                  ], name="analyse_find_insert_tmp1" + type_insert + "_" + sample.name)
                ], name="analyse_find_insert_tmp2" + type_insert + "_" + sample.name)
            ], name="analyse_find_insert" + type_insert + "_" + sample.name)
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
            self.fix_mate_by_coordinate,
            self.picard_mark_duplicates,
            self.recalibration,
            self.metrics,
            self.extract_sclip,
            self.extract_bam_unmap,
            self.extract_fastq_unmap,
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
   Puure().submit_jobs()
