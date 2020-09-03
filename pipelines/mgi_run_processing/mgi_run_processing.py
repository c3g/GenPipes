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
from __future__ import print_function, division, unicode_literals, absolute_import
import argparse
import logging
import os
import re
import sys
import itertools
import subprocess
import xml.etree.ElementTree as Xml
import math
import subprocess
import csv
import collections
import json

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs

from bfx.readset import parse_mgi_raw_readset_files
from bfx import bvatools
from bfx import picard
from bfx import fastqc
from bfx import tools
from bfx import run_processing_tools
from bfx import bash_cmd as bash

from pipelines import common

log = logging.getLogger(__name__)

class MGIRunProcessing(common.MUGQICPipeline):
    """
    MGI Run Processing Pipeline
    ================================

    The standard MUGQIC MGI Run Processing pipeline uses fastq files produced
    by the sequencer, then does demultiplexing. Finally, the
    pipeline runs some QCs on the raw data, on the fastq and on the alignment.

    Sample Sheets
    -------------

    The pipeline uses one input sample sheet.
    CURRENTLY BASED ON MGI RUN PROCESSING GOOGLE SHEET
    A csv file having the following columns :

    - Sample
    - Readset
    - Library
    - Project
    - Project ID
    - Protocol
    - Index
    - Pool ID
    - Run ID
    - Flowcell ID
    - Lane
    - Sequencer
    - Sequencer ID

    Example:
        Sample,Readset,Library,Project,Project ID,Protocol,Index,PoolID,RunID,FlowcellID,Lane,Sequencer,SequencerID
        LSPQ_Viral_Culture_dil_10-1_10cycles,LSPQ_Viral_Culture_dil_10-1_10cycles_PROD_000034-A01,PROD_000034-A01,LSPQ,,CleanPlex_MGI,1,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-2_10cycles,LSPQ_Viral_Culture_dil_10-2_10cycles_PROD_000034-B01,PROD_000034-B01,LSPQ,,CleanPlex_MGI,2,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-3_10cycles,LSPQ_Viral_Culture_dil_10-3_10cycles_PROD_000034-C01,PROD_000034-C01,LSPQ,,CleanPlex_MGI,3,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Nasal_Swab_Neg_ctl_10cycles,LSPQ_Nasal_Swab_Neg_ctl_10cycles_PROD_000034-A02,PROD_000034-A02,LSPQ,,CleanPlex_MGI,25,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-4_13cycles,LSPQ_Viral_Culture_dil_10-4_13cycles_PROD_000034-A03,PROD_000034-A03,LSPQ,,CleanPlex_MGI,28,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        LSPQ_Viral_Culture_dil_10-5_13cycles,LSPQ_Viral_Culture_dil_10-5_13cycles_PROD_000034-B03,PROD_000034-B03,LSPQ,,CleanPlex_MGI,29,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        L00241026_dil_10-1_13cycles,L00241026_dil_10-1_13cycles_PROD_000034-E03,PROD_000034-E03,LSPQ,,CleanPlex_MGI,33,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        L00241026_dil_10-2_13cycles,L00241026_dil_10-2_13cycles_PROD_000034-F03,PROD_000034-F03,LSPQ,,CleanPlex_MGI,34,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01
        Arctic_RT_reaction_13cycles,Arctic_RT_reaction_13cycles_PROD_000034-B04,PROD_000034-B04,LSPQ,,CleanPlex_MGI,4,LSPQ_Pool_01,1004MG01B,V300035341,2,Marie Curie,01

    """

    def __init__(self, protocol=None):
        self._protocol=protocol
        self.copy_job_inputs = []
        self.argparser.add_argument("-d", "--run", help="Run directory (mandatory)", required=True, dest="run_dir")
        self.argparser.add_argument("--lane", help="Lane number (mandatory)", type=int, required=True, dest="lane_number")
        self.argparser.add_argument("-r", "--readsets", help="Sample sheet for the MGI run to process (mandatory)", type=file, required=True)

        super(MGIRunProcessing, self).__init__(protocol)

    @property
    def readsets(self):
        if not hasattr(self, "_readsets"):
            self._readsets = self.load_readsets()
        return self._readsets

    @property
    def is_paired_end(self):
        if not hasattr(self, "_is_paired_end"):
            if self.readsets[0].run_type == "PAIRED_END":
                self._is_paired_end = True
            else:
                self._is_paired_end = False
        return self._is_paired_end

    @property
    def instrument(self):
        """
        The instrument id from the run folder
        """
        if not hasattr(self, "_instrument"):
#            if "_orig" in os.path.basename(os.path.dirname(self.run_dir)):
#                self._instrument = os.path.basename(os.path.dirname(self.run_dir)).split("_")[0]
#            else:
#                self._instrument = os.path.basename(os.path.dirname(self.run_dir))
            self._instrument = self.get_instrument()
        return self._instrument

    @property
    def seqtype(self):
        if not hasattr(self, "_seqtype"):
            self._seqtype = self.get_seqtype()
        return self._seqtype

    @property
    def year(self):
        """
        """
        if not hasattr(self, "_year"):
            self._year = "2020"
        return self._year 

    @property
    def run_id(self):
        """
        The run id from the readset objects.
        """
        if not hasattr(self, "_run_id"):
            #for readset in self.readsets:
            #    log.error(readset.run)
            runs = set([readset.run for readset in self.readsets])
            #log.error(runs)
            if len(list(runs)) > 1:
                _raise(SanitycheckError("Error: more than one run were parsed in the sample sheet... " + runs))
            else:
                self._run_id = list(runs)[0]
        return self._run_id

    @property
    def run_dir(self):
        if self.args.run_dir:
            return self.args.run_dir
        else:
            _raise(SanitycheckError("Error: missing '-d/--run_dir' option!"))

    @property
    def flowcell_id(self):
        """
        The flow cell ID from the run folder
        """
        if not hasattr(self, "_flowcell_id"):
            self._flowcell_id = os.path.basename(self.run_dir.rstrip('/'))
        return self._flowcell_id


    @property
    def lane_number(self):
        if self.args.lane_number:
            return self.args.lane_number
        else:
            _raise(SanitycheckError("Error: missing '--lane' option!"))

    @property
    def readset_file(self):
        if self.args.readsets:
            return self.args.readsets.name
        else:
            _raise(SanitycheckError("Error: missing '-r/--readsets' argument !"))

    @property
    def report_hash(self):
        if not hasattr(self, "_report_hash"):
            self._report_hash = {
                "run" : self.run_id,
                "instrument" : self.instrument,
                "flowcell" : self.flowcell_id,
                "lane" : self.lane_number,
                "seqtype" : self.seqtype,
                "sequencing_method" : "PAIRED_END" if self.is_paired_end else "SINGLE_END",
                "steps" : [],
                "run_validation" : []
            }
        return self._report_hash
    
    @property
    def report_inputs(self):
        if not hasattr(self, "_report_inputs"):
            self._report_inputs = {
                'qc' : {},
                'blast' : {},
                'mark_dup' : {},
                'align' : {}
            }
        return self._report_inputs

    def fastq(self):
        """
        *** In the future, may generate the fastq files from the raw CAL files. ***
        When fastq files are already generated by the sequencer, then copy all the fastq files in processing folder.
        The fastq (and metrics) files which match the readsets of the sample sheet are put in the Unaligned folder and renamed accordingly to their record in the sheet.
        All the fastq (and metrics) files which don't match any of the readsets in the sample sheet are put in the Unexpected folder without renaming.
        """

        jobs = []

        # First we copy all the run folder into the Unexpected folder
        unaligned_dir = os.path.join(self.output_dir, "Unaligned." + str(self.lane_number))
        unexpected_dir = os.path.join(unaligned_dir, 'Unexpected')
        copy_job = concat_jobs([
            bash.mkdir(unexpected_dir),
            bash.cp(
                os.path.join(self.run_dir, "L0" + str(self.lane_number), "."),
                unexpected_dir,
                recursive=True
            )
        ])

        # Then do the necessary moves and renamings, from the Unexpected folder to the Unaligned folders 
        fastq_job = Job()
        for readset in self.readsets:
            output_dir = os.path.join(self.output_dir, "Unaligned." + readset.lane, 'Project_' + readset.project, 'Sample_' + readset.name)

            copy_job.output_files.extend([
                os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.gz"),
                os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.fqStat.txt")
            ])
            fastq_job = concat_jobs([
                fastq_job,
                bash.mkdir(output_dir),
                bash.mv(
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.gz"),
                    readset.fastq1,
                ),
                bash.ln(
                    readset.fastq1,
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.gz")
                ),
                bash.touch(readset.fastq1),
                bash.mv(
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.fqStat.txt"),
                    re.sub("gz", "fqStat.txt", readset.fastq1),
                ),
                bash.ln(
                    re.sub("gz", "fqStat.txt", readset.fastq1),
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_1.fq.fqStat.txt")
                ),
                bash.touch(re.sub("gz", "fqStat.txt", readset.fastq1))
            ])

            if readset.run_type == "PAIRED_END":
                copy_job.output_files.extend([
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.gz"),
                    os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.fqStat.txt")
                ])
                fastq_job = concat_jobs([
                    fastq_job,
                    bash.mv(
                        os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.gz"),
                        readset.fastq2,
                    ),
                    bash.ln(
                        readset.fastq2,
                        os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.gz")
                    ),
                    bash.touch(readset.fastq2),
                    bash.mv(
                        os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.fqStat.txt"),
                        re.sub("gz", "fqStat.txt", readset.fastq2),
                    ),
                    bash.ln(
                        re.sub("gz", "fqStat.txt", readset.fastq2),
                        os.path.join(unexpected_dir, readset.flow_cell + "_L0" + readset.lane + "_" + readset.index + "_2.fq.fqStat.txt")
                    ),
                    bash.touch(re.sub("gz", "fqStat.txt", readset.fastq2))
            ])

#        copy_job.name = "fastq_copy." + self.run_id + "." + str(self.lane_number)
#        copy_job.samples = self.samples
##        copy_job.output_dependency = fastq_job.output_files

#        fastq_job.name = "fastq_rename." + self.run_id + "." + str(self.lane_number)
#        fastq_job.samples = self.samples

#        jobs.append(copy_job)
#        jobs.append(fastq_job)
        jobs.append(
            concat_jobs(
                [
                    copy_job,
                    fastq_job
                ],
                name="fastq." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
        )

        self.add_copy_job_inputs(jobs)
        return jobs

    def qc_graphs(self):
        """ 
        Generate some QC Graphics and a summary XML file for each sample using 
        [BVATools](https://bitbucket.org/mugqic/bvatools/).

        Files are created in a 'qc' subfolder of the fastq directory. Examples of
        output graphic:

        - Per cycle qualities, sequence content and sequence length;
        - Known sequences (adaptors);
        - Abundant Duplicates;
        """
        jobs = []

        for readset in self.readsets:
            output_dir = os.path.join(os.path.dirname(readset.fastq1), "qc")
            region_name = readset.name + "_" + readset.flow_cell + "_L0" + readset.lane + "_" + readset.index

            file1 = readset.fastq1
            file2 = readset.fastq2
            type = "FASTQ"

            jobs.append(
                concat_jobs([
                    bash.mkdir(output_dir),
                    bvatools.readsqc(
                        file1,
                        file2,
                        type,
                        region_name,
                        output_dir
                )],
                name="qc." + readset.name + ".qc." + self.run_id + "." + str(self.lane_number),
                samples=[readset.sample]
            ))
            self.report_inputs['qc'][readset.name] = os.path.join(output_dir, "mpsQC_" + region_name + "_stats.xml")

        self.add_to_report_hash("qc_graphs", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def fastqc(self):
        """
        """
        jobs = []
        for readset in self.readsets:
            input_dict = {
                readset.fastq1 : [
                    os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq1))),
                    os.path.join(os.path.dirname(readset.fastq1), "fastqc.R1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq1)))
#                ],
#                readset.index_fastq1 : [
#                    os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq1))),
#                    os.path.join(os.path.dirname(readset.index_fastq1), "fastqc.I1", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq1)))
                ]
            }
            if readset.run_type == "PAIRED_END":
                input_dict[readset.fastq2] = [
                    os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.fastq2))),
                    os.path.join(os.path.dirname(readset.fastq2), "fastqc.R2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.fastq2))),
                ]
#                if readset.index_type == "DUALINDEX":
#                    input_dict[readset.index_fastq2] = [
#                        os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.zip", os.path.basename(readset.index_fastq2))),
#                        os.path.join(os.path.dirname(readset.index_fastq2), "fastqc.I2", re.sub(".fastq.gz", "_fastqc.html", os.path.basename(readset.index_fastq2))),
#                    ]
            for input, outputs in input_dict.items():
                job_suffix = re.search('.*fastqc\.([IR][12])', os.path.dirname(outputs[0])).group(1)
                jobs.append(
                    concat_jobs([
                        bash.mkdir(
                            os.path.dirname(outputs[0]),
                            remove=True
                        ),
                        fastqc.fastqc(
                            input,
                            None,
                            outputs,
                            adapter_file=None,
                            use_tmp=True
                    )],
                    name="fastqc." + readset.name + "_" + job_suffix,
                    samples=[readset.sample]
                ))

        self.add_to_report_hash("fastqc", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def blast(self):
        """ 
        Run blast on a subsample of the reads of each sample to find the 20 most
        frequent hits.

        The `runBlast.sh` tool from MUGQIC Tools is used. The number of reads to
        subsample can be configured by sample or for the whole lane. The output will be
        in the `Blast_sample` folder, under the Unaligned folder.
        """
        jobs = []

        nb_blast_to_do = config.param('blast', 'nb_blast_to_do', type="posint")
        is_nb_blast_per_lane = config.param('blast', 'is_nb_for_whole_lane', type="boolean")

        if is_nb_blast_per_lane:
            nb_blast_to_do = int(nb_blast_to_do) // len(self.readsets)

        nb_blast_to_do = max(1, nb_blast_to_do)

        for readset in self.readsets:
            output_prefix = os.path.join(
                self.output_dir,
                "Unaligned." + readset.lane,
                "Blast_sample",
                readset.name + "_" + readset.flow_cell + "_L0" + readset.lane + "_" + readset.index
            )
            output = output_prefix + '.R1.RDP.blastHit_20MF_species.txt'
            self.report_inputs['blast'][readset.name] = os.path.join(output)
            current_jobs = [
                bash.mkdir(os.path.dirname(output))
            ]

            fasta_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.fasta".format(
                nb_blast_to_do=nb_blast_to_do
            )
            result_file = output_prefix + ".R1.subSampled_{nb_blast_to_do}.blastres".format(
                nb_blast_to_do=nb_blast_to_do
            )

            inputs = [readset.fastq1, readset.fastq2]
            command = "runBlast.sh " + str(nb_blast_to_do) + " " + output_prefix + " " + readset.fastq1 + " "
            if readset.fastq2:
                command += readset.fastq2
            current_jobs.append(
                Job(
                    inputs,
                    [output],
                    [
                        ["blast", "module_mugqic_tools"],
                        ["blast", "module_blast"]
                    ],
                    command=command
                )
            )

            # rRNA estimate using silva blast db, using the same subset of reads as the "normal" blast
            rrna_db = config.param('blast', 'rrna_db', required=False)
            if readset.is_rna and rrna_db:
                rrna_result_file = result_file + "Rrna"
                rrna_output = output_prefix + ".R1.subSampled_{nb_blast_to_do}.rrna".format(
                    nb_blast_to_do=nb_blast_to_do)
                command = """\
blastn \\
  -query {fasta_file} \\
  -db {db} \\
  -out {result_file} \\
  -perc_identity 80 \\
  -num_descriptions 1 \\
  -num_alignments 1""".format(
                    fasta_file=fasta_file,
                    result_file=rrna_result_file,
                    db=rrna_db
                )
                current_jobs.append(
                    Job(
                        [],
                        [],
                        [["blast", "module_blast"]],
                        command=command
                    )
                )

                command = """\
echo '{db}' \\
  > {output}""".format(
                    db=rrna_db,
                    output=rrna_output
                )
                current_jobs.append(
                    Job(
                        [],
                        [output],
                        [],
                        command=command
                    )
                    
                )

                command = """\
grep ">" {result_file} | \\
wc -l >> {output}""".format(
                    result_file=rrna_result_file,
                    output=rrna_output
                )
                current_jobs.append(
                    Job(
                        [],
                        [output],
                        [],
                        command=command
                    )
                )

                command = """\
grep ">" {fasta_file} | \\
wc -l >> {output}""".format(
                    fasta_file=fasta_file,
                    output=rrna_output
                )
                current_jobs.append(
                    Job(
                        [],
                        [output],
                        [],
                        command=command
                    )
                )

            # merge all blast steps of the readset into one job
            job = concat_jobs(
                current_jobs,
                name="blast." + readset.name + ".blast." + self.run_id + "." + str(self.lane_number),
                samples = [readset.sample]
            )
            jobs.append(job)

        self.add_to_report_hash("blast", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def align(self):
        """
        Align the reads from the fastq file, sort the resulting .bam and create an index
        of that .bam.

        An basic aligment is performed on a sample when the `SampleRef` field of the
        MGI sample sheet match one of the regexp in the configuration file and the
        corresponding genome (and indexes) are installed.

        `STAR` is used as a splice-junctions aware aligner when the sample
        `library_source` is `cDNA` or contains `RNA`; otherwise `BWA_mem` is used to
        align the reads.
        """
        jobs = []
        for readset in [readset for readset in self.readsets if readset.bam]:
            job = readset.aligner.get_alignment_job(readset)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_to_report_hash("align", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def picard_mark_duplicates(self):
        """
        Runs Picard mark duplicates on the sorted bam file.
        """
        jobs = []

        for readset in [readset for readset in self.readsets if readset.bam]:
            input_file_prefix = readset.bam + '.'
            input = input_file_prefix + "bam"
            output = input_file_prefix + "dup.bam"
            metrics_file = input_file_prefix + "dup.metrics"
            self.report_inputs['mark_dup'][readset.name] = metrics_file

            job = picard.mark_duplicates(
                [input],
                output,
                metrics_file
            )
            job.name = "picard_mark_duplicates." + readset.name + ".dup." + self.run_id + "." + str(self.lane_number)
            job.samples = [readset.sample]
            jobs.append(job)

        self.add_to_report_hash("picard_mark_duplicates", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def metrics(self):
        """
        This step runs a series of multiple metrics collection jobs and the output bam
        from mark duplicates.

        - Picard CollectMultipleMetrics: A collection of picard metrics that runs at the
        same time to save on I/O.
            - CollectAlignmentSummaryMetrics,
            - CollectInsertSizeMetrics,
            - QualityScoreDistribution,
            - MeanQualityByCycle,
            - CollectBaseDistributionByCycle
        - BVATools DepthOfCoverage: Using the specified `BED Files` in the sample sheet,
        calculate the coverage of each target region.
        - Picard CalculateHsMetrics: Calculates a set of Hybrid Selection specific
        metrics from the BAM file. The bait and interval list is automatically created
        from the specicied `BED Files`.
        """
        jobs = []

        for readset in [readset for readset in self.readsets if readset.bam]:
            job_list = readset.aligner.get_metrics_jobs(readset)
            for job in job_list:
                job.samples = [readset.sample]
            jobs.extend(job_list)

            if readset.is_rna:
                self.report_inputs['align'][readset.name] = [
                    readset.bam + '.metrics.alignment_summary_metrics',
                    readset.bam + '.metrics.insert_size_metrics',
                    os.path.join(os.path.dirname(readset.bam), readset.sample.name + "." + readset.library + '.rnaseqc.sorted.dup.metrics.tsv'),
                    readset.bam + '.metrics.rRNA.tsv'
                ]
            else:
                self.report_inputs['align'][readset.name] = [
                    readset.bam + '.metrics.alignment_summary_metrics',
                    readset.bam + '.metrics.insert_size_metrics',
                    readset.bam + ".dup.metrics",
                    readset.bam + ".metrics.verifyBamId.tsv",
                    readset.bam + ".metrics.targetCoverage.txt"
                ]

        self.add_to_report_hash("metrics", jobs)
        self.add_copy_job_inputs(jobs)
        return self.throttle_jobs(jobs)

    def md5(self):
        """
        Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
        util.

        One checksum file is created for each file.
        """
        jobs = []

        for readset in self.readsets:
            current_jobs = [
                Job(
                    [readset.fastq1],
                    [readset.fastq1 + ".md5"],
                    command="md5sum -b " + readset.fastq1 + " > " + readset.fastq1 + ".md5"
                )
            ]

            # Second read in paired-end run
            if readset.fastq2:
                current_jobs.append(
                    Job(
                        [readset.fastq2],
                        [readset.fastq2 + ".md5"],
                        command="md5sum -b " + readset.fastq2 + " > " + readset.fastq2 + ".md5"
                    )
                )

            # Alignment files
            if readset.bam:
                current_jobs.append(
                    Job(
                        [readset.bam + ".bam"],
                        [readset.bam + ".bam.md5"],
                        command="md5sum -b " + readset.bam + ".bam" + " > " + readset.bam + ".bam.md5"
                    )
                )
                current_jobs.append(
                    Job(
                        [],
                        [readset.bam + ".bai.md5"],
                        command="md5sum -b " + readset.bam + ".bai" + " > " + readset.bam + ".bai.md5"
                    )
                )

            job = concat_jobs(
                current_jobs,
                name="md5." + readset.name + ".md5." + self.run_id + "." + str(self.lane_number),
                samples=[readset.sample]
            )

            jobs.append(job)

        if config.param('md5', 'one_job', required=False, type="boolean"):
            job = concat_jobs(
                jobs,
                name="md5." + self.run_id + "." + str(self.lane_number),
                samples=self.readsets
            )
            self.add_copy_job_inputs([job])
            return [job]
        else:
            self.add_copy_job_inputs(jobs)
            return jobs

    def copy(self):
        """
        Copy the whole processing foler to where they can be serve or loaded into a LIMS
        """
        copy_job = concat_jobs(
            [
                bash.mkdir(
                    os.path.join(
                        config.param("copy", "mgi_runs_root", type="dirpath"),
                        self.seqtype,
                        self.year
                    )
                ),
                bash.cp(
                    os.path.join(self.output_dir, "."),
                    os.path.join(
                        config.param("copy", "mgi_runs_root", type="dirpath"),
                        self.seqtype,
                        self.year,
                        self.instrument + "_" + self.flowcell_id + "_" + self.run_id
                    ),
                    recursive=True,
                    input_dependency=self.copy_job_inputs
                )
            ],
            name="copy." + self.run_id + "." + str(self.lane_number),
            samples=self.samples
        )
        return [copy_job]

    def copy_legacy(self):
        """
        Copy processed files to another place where they can be served or loaded into a
        LIMS.

        The destination folder and the command used can be set in the configuration
        file.

        An optional notification can be sent before the copy. The command used is in the configuration file.
        """
        inputs = self.copy_job_inputs
        jobs_to_concat = []

        # Notification
        output1 = os.path.join(self.output_dir, "notificationProcessingComplete." + str(self.lane_number) + ".out")
        output2 = os.path.join(self.output_dir, "notificationCopyStart." + str(self.lane_number) + ".out")

        notification_command = config.param('copy', 'notification_command', required=False)
        if notification_command:
            job = Job(
                inputs,
                [output1, output2],
                name="start_copy_notification." + self.run_id + "." + str(self.lane_number),
                samples=self.samples
            )
            job.command = notification_command.format(
                technology=config.param('copy', 'technology'),
                output_dir=self.output_dir,
                run_id=self.run_id,
                output1=output1,
                output2=output2,
                lane_number=self.lane_number
            )
            jobs_to_concat.append(job)

        # Actual copy
        full_destination_folder = config.param('copy', 'destination_folder', type="dirpath") + os.path.basename(self.run_dir)
        output = os.path.join(full_destination_folder, "copyCompleted." + str(self.lane_number) + ".out")

        exclude_bam = config.param('copy', 'exclude_bam', required=False, type='boolean')
        exclude_fastq_with_bam = config.param('copy', 'exclude_fastq_with_bam', required=False, type='boolean')
        if exclude_bam and exclude_fastq_with_bam:
            log.warn("Excluding both BAM and fastq files")

        excluded_files = []

        if exclude_bam or exclude_fastq_with_bam:
            for readset in [readset for readset in self.readsets if readset.bam]:
                if exclude_bam:
                    excluded_files.append(readset.bam + ".bam*")
                    excluded_files.append(readset.bam + ".bai*")
                if exclude_fastq_with_bam and not exclude_bam:
                    excluded_files.append(readset.fastq1)
                    if readset.fastq2:
                        excluded_files.append(readset.fastq2)

        if self.run_dir != self.output_dir:
            copy_command_run_folder = config.param('copy', 'copy_command', required=False).format(
                exclusion_clauses="",
                lane_number=self.lane_number,
                run_id=self.run_id,
                source=self.run_dir,
                run_name=os.path.basename(self.run_dir)
            )
            jobs_to_concat.append(
                Job(
                    inputs,
                    [output],
                    command=copy_command_run_folder,
                    samples=self.samples
                )
            )

        copy_command_output_folder = config.param('copy', 'copy_command', required=False).format(
            exclusion_clauses="\\\n".join(
                [" --exclude '" + excludedfile.replace(self.output_dir + os.sep, "") + "'" for excludedfile in excluded_files]),
            lane_number=self.lane_number,
            run_id=self.run_id,
            source=self.output_dir,
            run_name=os.path.basename(self.run_dir)
        )
        jobs_to_concat.append(
            Job(
                inputs,
                [output],
                command=copy_command_output_folder,
                samples=self.samples
            )
        )
        jobs_to_concat.append(
            Job(
                command="touch " + output,
                samples=self.samples
            )
        )

        job = concat_jobs(
            jobs_to_concat,
            name="copy." + self.run_id + "." + str(self.lane_number)
        )

        return [job]

    def report(self):
        """
        Generate en JSON file reporting the whole pipeline
        """
        jobs = []

        report_dir = os.path.join(self.output_dir, "report")
        sample_report_dir  = os.path.join(report_dir, "sample_json")

        run_validation_inputs = []

        for readset in self.readsets:
            output_file = os.path.join(sample_report_dir, readset.name + ".report.json")
            jobs.append(
                concat_jobs([
                    bash.mkdir(sample_report_dir),
                    tools.run_validation_sample_report(
                        readset,
                        self.report_inputs,
                        output_file
                    )],
                    name="sample_report." + readset.name + "." + self.run_id + "." + str(readset.lane)
                )
            )
            run_validation_inputs.append(output_file)

        run_validation_report_json = os.path.join(report_dir, self.run_id + "." + str(self.lane_number) + ".run_validation_report.json")
        general_information_file = os.path.join(self.output_dir, self.run_id + "." + str(self.lane_number) + ".general_information.json")
        with open(general_information_file, 'w') as out_json:
            json.dump(self.report_hash, out_json, indent=4)

        log.error("\n".join(run_validation_inputs))
        jobs.append(
            concat_jobs([
                bash.mkdir(report_dir),
                tools.run_validation_aggregate_report(
                    general_information_file,
                    run_validation_inputs,
                    run_validation_report_json
                )],
                name="report." + self.run_id + "." + str(self.lane_number)
            )
        )

        return jobs

    #
    # Utility methods
    #

    def add_copy_job_inputs(self, jobs):
        for job in jobs:
            # we first remove dependencies of the current job, since we will have a dependency on that job
            self.copy_job_inputs = [item for item in self.copy_job_inputs if item not in job.input_files]
            self.copy_job_inputs.extend(job.output_files)

    def add_to_report_hash(self, step_name, jobs=[]):
        report_hash = {
            "step_name" : step_name,
            "jobs" : [{
                "job_name" : job.name,
                "input_files" : [os.path.relpath(input_file, os.path.join(self.output_dir, "report")) for input_file in job.input_files],
                "output_files" : [os.path.relpath(output_file, os.path.join(self.output_dir, "report")) for output_file in job.output_files],
                "command" : job.command,
                "modules" : job.modules
            } for job in jobs]
        }
        self.report_hash["steps"].append(report_hash)

    def get_read2cycles(self):
        """
        Parse the BioInfo.csv file for the number of read cycles in the run
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read1 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not Read1 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_read2cycles(self):
        """
        Parse the BioInfo.csv file for the number of read cycles in the run
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Read2 Cycles":
                return row[1]
        _raise(SanitycheckError("Could not Read2 Cycles from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_instrument(self):
        """
        Parse the BioInfo.csv file for the instrument name the run has been running on
        """
        bioinfo_csv = csv.reader(open(os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv"), 'rb'))
        for row in bioinfo_csv:
            if row[0] == "Machine ID":
                return row[1]
        _raise(SanitycheckError("Could not find intrument from " + os.path.join(self.run_dir, "L0" + str(self.lane_number), "BioInfo.csv")))

    def get_seqtype(self):
        """
        Determine which kind of sequencing (iseq, miseq, novaseq, hiseqx, hiseq4000 or hiseq2500) was performed,
        depending on the instrument used for the run
        """

        instrument = self.instrument
        instrument_file = config.param('DEFAULT', 'instrument_list_file', type='filepath', required=False)
        if not (instrument_file and os.path.isfile(instrument_file)):
            instrument_file = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "resources", 'instrument_list.csv')

        return subprocess.check_output("grep -m1 '"+instrument+"' %s | awk -F',' '{print $3}'" % instrument_file, shell=True).strip()

    def load_readsets(self):
        """
        Parse the sample sheet and return a list of readsets.
        """

        return parse_mgi_raw_readset_files(
            self.readset_file,
            self.run_dir,
            self.flowcell_id,
            int(self.lane_number),
            self.output_dir
        )

    def submit_jobs(self):
        super(MGIRunProcessing, self).submit_jobs()

    def throttle_jobs(self, jobs):
        """
        Group jobs of the same task (same name prefix) if they exceed the configured threshold number.
        """

        max_jobs_per_step = config.param('default', 'max_jobs_per_step', required=False, type="int")
        jobs_by_name = collections.OrderedDict()
        reply = []

        # group jobs by task (name)
        for job in jobs:
            jobs_by_name.setdefault(job.name.split(".", 1)[0], []).append(job)

        # loop on all task
        for job_name in jobs_by_name:
            current_jobs = jobs_by_name[job_name]
            if max_jobs_per_step and 0 < max_jobs_per_step < len(current_jobs):
                # we exceed the threshold, we group using 'number_task_by_job' jobs per group
                number_task_by_job = int(math.ceil(len(current_jobs) / max_jobs_per_step))
                merged_jobs = []
                for x in range(max_jobs_per_step):
                    if x * number_task_by_job < len(current_jobs):
                        merged_jobs.append(concat_jobs(
                            current_jobs[x * number_task_by_job:min((x + 1) * number_task_by_job, len(current_jobs))],
                            job_name + "." + str(x + 1) + "." + self.run_id + "." + str(self.lane_number)))
                reply.extend(merged_jobs)
            else:
                reply.extend(current_jobs)
        return reply

    @property
    def steps(self):
        return [
            self.fastq,
            self.qc_graphs,
            self.fastqc,
            self.blast,
            self.align,
            self.picard_mark_duplicates,
            self.metrics,
        #    self.md5,
            self.copy,
            self.report
        ]

if __name__ == '__main__':

    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        MGIRunProcessing()
