#!/usr/bin/env python

######
# Copyright etc
#####

# Python Standard Modules
import logging
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
import utils.utils

from bfx import trimmer
from bfx import bwa
from bfx import samtools
from bfx import locatit
from bfx import vardict

from bfx import bash_cmd as bash
from core.config import config

log = logging.getLogger(__name__)

class DOvEE_gene():
    """
    DOvEE Gene Pipeline
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments? 
        super(DOvEE_gene, self).__init__(protocol)

    @property
    def output_dirs(self):
        dirs = {
                'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
                'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
                'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
                'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
                'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
                'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
                }
        return dirs

    def trimmer():
        """
        Read trimming with AGeNT Trimmer.
        """

        jobs = []

        for readset in self.readsets:
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim."
            fastq1 = readset.fastq1
            fastq2 = readset.fastq2
            job = trimmer.trimmer(
                fastq1,
                fastq2,
                trim_file_prefix
                )

             jobs.append(
                 concat.jobs(
                     [
                         command="mkdir =p " + trim_directory,
                         job
                        ], name="trimmer." + readset.name
                     )                               
                 )
            return jobs

    def bwa_mem_samtools_sort(self):
        """
        The trimmed reads are aligned to the reference genome with bwa mem, followed by sorting with samtools.
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            
            fastq1 = ""
            fastq2 = ""
            # Find input readset FASTQs first from previous trimmer job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + "pair1.fastq.gz", prefix + "pair2.fastq.gz"])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    prefix = os.path.join(
                        self.output_dirs["raw_reads_directory"],
                        readset.sample.name,
                        re.sub("\.bam$", ".", os.path.basename(readset.bam))
                    )
                    candidate_input_files.append([prefix + ".single.fastq.gz"])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.dirname(readset_bam)),
                        pipe_jobs(
                            [
                                bwa.mem(
                                    fastq1,
                                    fastq2,
                                    read_group="'@RG" + \
                                        "\\tID:" + readset.name + \
                                        "\\tSM:" + readset.sample.name + \
                                        "\\tLB:" + readset.library + \
                                        "\\tPU:" + readset.library + \
                                        ("\\tCN:" + config.param('bwa_mem', 'sequencing_center')) + \
                                        ("\\tPL:" + config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem_sambamba_sort_sam', 'sequencing_technology', required=False) else "Illumina") + \
                                        "'",
                                        ini_section="bwa_mem"
                                ),
                                samtools.view(
                                    "/dev/stdin",
                                    None,
                                    ""
                                ),
                                samtools.sort(
                                    "/dev/stdin",
                                    readset_bam,
                                    other_options=config.param('samtools_sort', 'sort_options', required=False),
                                )
                            ]
                        )
                    ],
                    name="bwa_mem_samtools_sort." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs

    def locatit_dedup_bam(self):
        """
        Deduplicated bam files with AGeNT locatit in hybrid or duplex mode.
        """

        jobs =[]

        for sample in self.samples:
           alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
           input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

           if 'saliva' in sample.name:
               output_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.bam")
               output_hybrid = os.path.join(alignment_directory, sample.name + "dedup.hybrid.bam")
               covered_bed = config.param('locatit', 'covered_bedv8', paramtype='filepath')

               jobs.append(
                       concat_jobs(
                           [
                               locatit.dedup(
                                   input_bam,
                                   output_duplex,
                                   covered_bed,
                                   duplex
                                   ),
                               locatit.dedup(
                                   input_bam,
                                   output_hybrid,
                                   covered_bed,
                                   hybrid
                                   )
                            ],
                           name='locatit_dedup' + sample.name
                           )
                       )
            elif 'brush' in sample.name:
                output_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.bam")
                covered_bed = config.param('locatit', 'covered_bedv7', paramtype='filepath')

                jobs.append(
                        locatit.dedup(
                            input_bam,
                            output_duplex,
                            covered_bed,
                            duplex
                            )
                        )
    

    def samtools_sort(self):
        """
        Sort deduplicated bams by coordinate with samtools.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            
            if 'saliva' in sample.name: #temporary solution
                input_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.bam")
                input_hybrid = os.path.join(alignment_directory, sample.name + "dedup.hybrid.bam")

                jobs.append(
                        concat_jobs(
                            [
                                samtools.sort(
                                    input_duplex,
                                    resub(".duplex.", ".duplex.sorted.", input_duplex)
                                    ),
                                samtools.sort(
                                    input_hybrid,
                                    resub(".hybrid.", ".hybrid.sorted.", input_hybrid)
                                    )
                            ],
                            name = "samtools.sort" + sample.name
                        )
                    )
            elif 'brush' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.bam")

                jobs.append(
                        samtools.sort(
                            input_duplex,
                            resub(".duplex.", ".duplex.sorted.", input_duplex)
                            )
                        )
        return jobs

    def samtools_index(self):
        """
        Index deduplicated and sorted bams with samtools.
        """

        jobs=[]

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            if 'saliva' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.sorted.bam")
                input_hybrid = os.path.join(alignment_directory, sample.name + "dedup.hybrid.sorted.bam")

                jobs.append(
                        concat.jobs(
                            [
                                samtools.index(
                                    input_duplex
                                    ),
                                samtools.index(
                                    input_hybrid
                                    )
                            ],
                            name = "samtools.index" + sample.name
                            )
                        )
            elif 'brush' in sample.name:
                input_duplex = os.path.join(alignment_directory, sample.name + "dedup.duplex.sorted.bam")

                jobs.append(
                        samtools.index(
                            input_duplex
                            )
                        )
            return jobs

    def vardict_single(self):
        """
        Variant calling with vardict.
        """

        jobs = []
        
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            variant_directory = os.path.join(self.output_dirs['variant_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.sorted.bam")
            
            if sample == 'brush':
               freq=0.001
               region=  #target bed file

               jobs.append(
                       vardict.single_java(
                           input_bam,
                           sample.name,
                           None,
                           sv=True,
                           freq,
                           region
                           )
                       )
            elif sample == 'saliva':
                freq=0.1
                region= #targetv8 bed file

                jobs.append(
                        vardict.single_java(
                            input_bam,
                            sample.name,
                            None,
                            sv=False,
                            freq,
                            region
                            )
                        )

            return Jobs

    @property
    def steps(self):
        return [
                self.trimmer,
                self.bwa_mem_samtools_sort,
                self.locatit_dedup_bam,
                self.samtools_sort,
                self.samtools_index,
                self.vardict_single,
                ]
