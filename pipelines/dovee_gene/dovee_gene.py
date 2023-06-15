#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2023 GenAP, McGill University and Genome Quebec Innovation Centre
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
import os
import re
import sys
import argparse

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config, _raise, SanitycheckError
from core.job import Job, concat_jobs, pipe_jobs
import utils.utils
from pipelines import common

from bfx import agent
from bfx import bwa
from bfx import samtools
from bfx import vardict
from bfx import hmm
from bfx import ichorCNA
from bfx import fastp
from bfx import mosdepth
from bfx import picard2
from bfx import multiqc
from bfx import bcftools
from bfx import conpair

from bfx import bash_cmd as bash
from core.sample_dovee_pairs import parse_dovee_pair_file
from core.dovee_design import parse_dovee_design_file 

log = logging.getLogger(__name__)

class DOvEE_gene(common.Illumina):
    """
    DOvEE Gene Pipeline
    """
    def __init__(self, protocol=None):
        self._protocol = protocol
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="File indicating whether sample is saliva or brush", type=argparse.FileType('r')) # only needed for vardict protocol
        self.argparser.add_argument("-p", "--pairs", help="File with sample pairing information", type=argparse.FileType('r')) # only needed for copy number protocol and for conpair compare bams step in vardict protocol
        self.argparser.add_argument("-t", "--type", help="Type of pipeline (default vardict)", choices=["vardict", "copy-number"], default="vardict") 
        super(DOvEE_gene, self).__init__(protocol)

    @property
    def multiqc_inputs(self):
        if not hasattr(self, "_multiqc_inputs"):
            self._multiqc_inputs = []
        return self._multiqc_inputs

    @multiqc_inputs.setter
    def multiqc_inputs(self, value):
        self._multiqc_inputs = value

    @property
    def output_dirs(self):
        dirs = {
                'raw_reads_directory': os.path.relpath(os.path.join(self.output_dir, 'raw_reads'), self.output_dir),
                'trim_directory': os.path.relpath(os.path.join(self.output_dir, 'trim'), self.output_dir),
                'alignment_directory': os.path.relpath(os.path.join(self.output_dir, 'alignment'), self.output_dir),
                'metrics_directory': os.path.relpath(os.path.join(self.output_dir, 'metrics'), self.output_dir),
                'variants_directory': os.path.relpath(os.path.join(self.output_dir, 'variants'), self.output_dir),
                'report_directory': os.path.relpath(os.path.join(self.output_dir, 'report'), self.output_dir),
                'wig_directory' : os.path.relpath(os.path.join(self.output_dir, 'wig'), self.output_dir),
                'cna_directory' : os.path.relpath(os.path.join(self.output_dir, 'cna'), self.output_dir)
                }
        return dirs

    @property
    def dovee_pairs(self):
        if not hasattr(self, "_dovee_pairs"):
            self._dovee_pairs = parse_dovee_pair_file(
                    self.args.pairs.name,
                    self.samples,
                    )
        return self._dovee_pairs

    @property
    def contrasts(self):
        if not hasattr(self, "_contrasts"):
            if self.args.design:
                self._contrasts = parse_dovee_design_file(self.args.design.name, self.samples)
            else:
                self.argparser.error("argument -d/--design is required for contrast")
        return self._contrasts

    def trimmer(self):
        """
        Read trimming with AGeNT Trimmer from Agilent Genomics NextGen tool kit. 
        """

        jobs = []

        for readset in self.readsets:
            trim_directory = os.path.join(self.output_dirs["trim_directory"], readset.sample.name)
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim")
            fastq1 = readset.fastq1
            fastq2 = readset.fastq2

            # Trimmer does not create STATS file if one already exists. Remove any existing STATS file.
            job_rm = Job(
                    command="""\
    if [ -f {STATS} ]; then
        rm {STATS}
    fi""".format(
                    STATS=trim_file_prefix + "_STATS.properties"
                    )
                )
            
            job = agent.trimmer(
                fastq1,
                fastq2,
                trim_file_prefix
                )

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(trim_directory),
                            job_rm,
                            job,
                            bash.ln(os.path.relpath(trim_file_prefix + "_R1.fastq.gz", os.path.dirname(trim_file_prefix + ".pair1.fastq.gz")), 
                                trim_file_prefix + ".pair1.fastq.gz",
                                trim_file_prefix + "_R1.fastq.gz"), # trimmer only allows specifying outfile prefix to which it adds _R1/2.fastq.gz
                            bash.ln(os.path.relpath(trim_file_prefix + "_R2.fastq.gz", os.path.dirname(trim_file_prefix + ".pair2.fastq.gz")), 
                                trim_file_prefix + ".pair2.fastq.gz",
                                trim_file_prefix + "_R2.fastq.gz")
                        ], name="agent_trimmer." + readset.name,
                        samples = [readset.sample]
                        )
                 )
        return jobs
    
    def fastp(self):
        """
        Generate basic QC metrics for trimmed reads with fastp.
        """

        jobs = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim.")
            trim1 = trim_file_prefix + "pair1.fastq.gz"
            trim2 = trim_file_prefix + "pair2.fastq.gz"

            output_json_path = os.path.join(os.path.dirname(trim1), readset.name + ".trim.fastp.json")
            output_html_path = os.path.join(os.path.dirname(trim1), readset.name + ".trim.fastp.html")

            link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(link_directory),
                            fastp.basic_qc(
                                trim1,
                                trim2,
                                output_json_path,
                                output_html_path
                                ),
                            bash.ln(
                                os.path.relpath(output_json_path, link_directory),
                                os.path.join(link_directory, readset.name + ".trim.fastp.json"),
                                output_json_path
                                )
                        ],
                        name = "fastp." + readset.name,
                        samples = [readset.sample]
                        )
                    )
            self.multiqc_inputs.append(output_json_path)

        return jobs

    def bwa_mem_samtools_sort(self):
        """
        The trimmed reads are aligned to the reference genome with bwa mem, followed by sorting with samtools.
        The bwa mem output is piped directory into samtools to avoid saving the intermediate file.
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join(self.output_dirs['trim_directory'], readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join(self.output_dirs["alignment_directory"], readset.sample.name, readset.name)
            readset_bam = os.path.join(alignment_directory, readset.name + ".sorted")
            
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

            else:
                _raise(SanitycheckError("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!"))

            jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(alignment_directory),
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
                                        ("\\tPL:" + config.param('bwa_mem', 'sequencing_technology') if config.param('bwa_mem', 'sequencing_technology', required=False) else "\\tPL:Illumina") + \
                                        "'",
                                    ini_section="bwa_mem"
                                ),
                                samtools.view(
                                    "-",
                                    None,
                                    "-b -@ 7"
                                ),
                                samtools.sort(
                                    "-",
                                    readset_bam,
                                )
                            ]
                        )
                    ],
                    name="bwa_mem." + readset.name,
                    samples=[readset.sample]
                )
            )

        return jobs

    def samtools_merge(self):
        """
        Merges sorted bam files for each sample, if there are multiple readsets for the sample. If only a single readset file is found, a symlink to the bam is created. 
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            readset_bams = [os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets]
            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                
                # remove any existing symlink
                job_rm = Job(
                    command="""\
if [ -L {LINK} ]; then
    rm {LINK}
fi""".format(
                    LINK=sample_bam
                    )
                )

                jobs.append(
                        concat_jobs(
                            [
                                job_rm,
                                bash.ln(
                                    os.path.relpath(readset_bam, os.path.dirname(sample_bam)),
                                    sample_bam,
                                    input=readset_bam
                                    )
                                ], name = "symlink_readset_sample_bam." + sample.name,
                            samples = [sample]
                            )
                        )

            # If multiple readsets exist for the sample, merge readset bams into one sample bam
            elif len(sample.readsets) > 1:
                job = samtools.merge(
                        sample_bam,
                        readset_bams
                        )
                job.name = "samtools_merge." + sample.name
                job.samples = [sample]

                jobs.append(job)

        return jobs
                        
    def locatit_dedup_bam(self):
        """
        Deduplicate bam files with AGeNT locatit in hybrid (saliva only) or duplex (both saliva and brush) mode.
        Used for SureSelect samples in vardict protocol.
        Saliva and brush samples are identified via the design file. 
        """

        jobs = []
        metrics_inputs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            sorted_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")

            if sample in self.contrasts.salivas:
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                index_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bai")
                output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
                covered_bed = config.param('agent_locatit', 'covered_bedv8', param_type='filepath')
                source = "_saliva"

                jobs.append(
                        concat_jobs(
                            [
                                agent.locatit(
                                    input_bam,
                                    output_duplex,
                                    covered_bed,
                                    "v2Duplex"
                                    ),
                                agent.locatit(
                                    input_bam,
                                    output_hybrid,
                                    covered_bed,
                                    "v2Hybrid"
                                    ),
                                bash.ln(
                                    os.path.relpath(output_duplex, os.path.dirname(sorted_bam)),
                                    sorted_bam,
                                    output_duplex
                                    ),
                                bash.ln(
                                    os.path.relpath(index_duplex, os.path.dirname(sorted_bam)),
                                    sorted_bam + ".bai",
                                    index_duplex
                                    ),
                                bash.ln(
                                    os.path.relpath(re.sub(".bam", ".properties", output_duplex), os.path.dirname(sorted_bam)),
                                    os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.properties"),
                                    re.sub(".bam", ".properties", output_duplex)
                                    )
                            ],
                            name='agent_locatit.' + sample.name,
                            samples=[sample]
                            )
                        )
            elif sample in self.contrasts.brushes:
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                index_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bai")
                covered_bed = config.param('agent_locatit', 'covered_bedv7', param_type='filepath')
                source = "_brush"

                jobs.append(
                        concat_jobs(
                            [
                                agent.locatit(
                                    input_bam,
                                    output_duplex,
                                    covered_bed,
                                    "v2Duplex"
                                    ),
                              #  bash.ln(
                              #      os.path.relpath(output_duplex, os.path.dirname(sorted_bam)),
                              #      sorted_bam,
                              #      output_duplex
                              #      ),
                              #  bash.ln(
                              #      os.path.relpath(index_duplex, os.path.dirname(sorted_bam)),
                              #      sorted_bam + ".bai",
                              #      index_duplex
                              #      ),
                                bash.ln(
                                    os.path.relpath(re.sub(".bam", ".properties", output_duplex), os.path.dirname(sorted_bam)),
                                    os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.properties"),
                                    re.sub(".bam", ".properties", output_duplex)
                                    )
                                ],
                            name='agent_locatit.' + sample.name,
                            samples=[sample]
                            )
                        )
                        
            else:
                _raise(SanitycheckError("Error: sample \"" + sample.name +
                "\" is not included in the design file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))

            metrics_inputs.append(os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.properties"))

        locatit_metrics = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", "locatit_dedup.duplex.stats_mqc.txt")
        
        parse_job=Job(
                metrics_inputs,
                [locatit_metrics],
                command="""\
echo -e "# plot_type: 'table'
# section_name: 'LocatIt'
# description: 'stats from deduplication with AGeNT LocatIt'
# headers:
#   DuplexConsensusReads:
#       title: 'Duplex Consensus Reads'
#       description: 'reads called as duplex consensus'
#       format: '{{:,.0f}}'
#       placement: 900
#   SingleConsensusReads:
#       title: 'Single Consensus Reads'
#       description: 'reads called as single consensus'
#       format: '{{:,.0f}}'
#   DupRate:
#       title: 'Duplication Rate'
#       description: 'duplication rate calculated by locatit' 
#       format: '{{:,.3f}}'
#       placement: 1020
Sample\\tDuplexConsensusReads\\tSingleConsensusReads\\tDupRate" > {output}
for f in {input}; do
    sample=$(basename $f .dedup.duplex.properties)
    duplex_consensus=$(grep "DUPLEX_VS_SINGLE_CONSENSUS_READS" $f | awk -F'=' '{{ print $2 }}' | awk -F',' '{{ print $1 }}' )
    single_consensus=$(grep "DUPLEX_VS_SINGLE_CONSENSUS_READS" $f | awk -F',' '{{ print $2 }}')
    dup_rate=$(grep "DUPLICATE_RATE" $f | awk -F'=' '{{ print $2 }}')

    echo -e "$sample\\t$duplex_consensus\\t$single_consensus\\t$dup_rate" >> {output}
done""".format(
    input=" ".join(["" + input for input in metrics_inputs]),
    output=locatit_metrics
            )
        )

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")),
                        parse_job
                        ],
                    name = "parse_locatit_metrics"
                    )
                )

        self.multiqc_inputs.append(locatit_metrics)

        return jobs

    def locatit_hybrid_dedup(self):
        """
        Deduplicate bam files with locatit in hybrid mode only. Used for low pass copy-number protocol.
        """

        jobs = []
        metrics_inputs = []

        for sample in self.samples:
           alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
           input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
           output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
           index_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bai")
           sorted_bam = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.sorted.bam")
           
           if sample.name in [pair.brush.name for pair in self.dovee_pairs.values()]:
               source = "_brush"
           elif sample.name in [pair.saliva.name for pair in self.dovee_pairs.values()]:
               source = "_saliva"
           else:
               _raise(SanitycheckError("Error: sample \"" + sample.name + "\" is not included in the pairs file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))

           job = agent.locatit(
                   input_bam,
                   output_hybrid,
                   None, # covered bed not needed
                   "v2Hybrid"
                   )
            
           jobs.append(
                   concat_jobs(
                       [
                           job,
                           bash.ln(
                               os.path.relpath(output_hybrid, os.path.dirname(sorted_bam)),
                               sorted_bam,
                               output_hybrid
                               ),
                           bash.ln(
                               os.path.relpath(index_hybrid, os.path.dirname(sorted_bam)),
                               sorted_bam + ".bai",
                               index_hybrid
                               ),
                           bash.ln(
                               os.path.relpath(re.sub(".bam", ".properties", output_hybrid), os.path.dirname(sorted_bam)),
                               os.path.join(alignment_directory, sample.name + source + ".dedup.hybrid.properties"),
                               re.sub(".bam", ".properties", output_hybrid),
                               )
                        ],
                        name='agent_locatit.' + sample.name,
                        samples=[sample]
                        )
                    )
           
           metrics_inputs.append(os.path.join(alignment_directory, sample.name + source + ".dedup.hybrid.properties"))
                
        locatit_metrics = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", "locatit_dedup.hybrid.stats_mqc.txt")
        
        parse_job=Job(
                metrics_inputs,
                [locatit_metrics],
                command="""\
echo -e "# plot_type: 'table'
# section_name: 'LocatIt'
# description: 'stats from deduplication with AGeNT LocatIt'
# headers:
#   DuplexConsensusReads:
#       title: 'Duplex Consensus Reads'
#       description: 'reads called as duplex consensus'
#       format: '{{:,.0f}}'
#       placement: 900
#   SingleConsensusReads:
#       title: 'Single Consensus Reads'
#       description: 'reads called as single consensus'
#       format: '{{:,.0f}}'
#   DupRate:
#       title: 'Duplication Rate'
#       description: 'duplication rate calculated by locatit' 
#       format: '{{:,.3f}}'
#       placement: 1020
Sample\\tDuplexConsensusReads\\tSingleConsensusReads\\tDupRate" > {output}
for f in {input}; do
    sample=$(basename $f .dedup.hybrid.properties)
    duplex_consensus=$(grep "DUPLEX_VS_SINGLE_CONSENSUS_READS" $f | awk -F'=' '{{ print $2 }}' | awk -F',' '{{ print $1 }}' )
    single_consensus=$(grep "DUPLEX_VS_SINGLE_CONSENSUS_READS" $f | awk -F',' '{{ print $2 }}')
    dup_rate=$(grep "DUPLICATE_RATE" $f | awk -F'=' '{{ print $2 }}')

    echo -e "$sample\\t$duplex_consensus\\t$single_consensus\\t$dup_rate" >> {output}
done""".format(
    input=" ".join(["" + input for input in metrics_inputs]),
    output=locatit_metrics
            )
        )

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")),
                        parse_job
                        ],
                    name = "parse_locatit_metrics"
                    )
                )

        self.multiqc_inputs.append(locatit_metrics)

        return jobs

    def creak_hybrid_dedup(self):
        """
        Deduplicate bam files with creak in hybrid mode only. Used for low pass copy-number protocol.
        """

        jobs = []
        metrics_inputs = []
        
        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
            index_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bai")
            sorted_bam = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.sorted.bam")

            if sample.name in [pair.brush.name for pair in self.dovee_pairs.values()]:
                source = "_brush"
            elif sample.name in [pair.saliva.name for pair in self.dovee_pairs.values()]:
                source = "_saliva"
            else:
                _raise(SanitycheckError("Error: sample \"" + sample.name + "\" is not included in the pairs file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))
            
            job = agent.creak(
                    input_bam,
                    output_hybrid,
                    "HYBRID",
                    None
                    )

            jobs.append(
                    concat_jobs(
                        [
                            job,
                            bash.ln(
                                os.path.relpath(output_hybrid, os.path.dirname(sorted_bam)),
                                sorted_bam,
                                output_hybrid
                                ),
                            bash.ln(
                                os.path.relpath(index_hybrid, os.path.dirname(sorted_bam)),
                                sorted_bam + ".bai",
                                index_hybrid
                                ),
                            bash.ln(
                                os.path.relpath(re.sub(".bam", ".stats", output_hybrid), os.path.dirname(sorted_bam)),
                                os.path.join(alignment_directory, sample.name + source + ".dedup.hybrid.stats"),
                                re.sub(".bam", ".stats", output_hybrid),
                                )
                        ],
                        name='agent_creak.' + sample.name,
                        samples=[sample]
                        )
                    )

            metrics_inputs.append(os.path.join(alignment_directory, sample.name + source + ".dedup.hybrid.stats"))

        creak_metrics = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", "creak_dedup.hybrid.stats_mqc.txt")
        
        parse_job=Job(
                metrics_inputs,
                [creak_metrics],
                command="""\
echo -e "# plot_type: 'table'
# section_name: 'CReaK'
# description: 'stats from deduplication with AGeNT CReaK'
# headers:
#   ReadPairs:
#       title: 'Read Pairs'
#       description: 'correctly-paired read pairs for MBC Consensus calling'
#       format: '{{:,.0f}}'
#       placement: 900
#   DuplexConsensusPairs:
#       title: 'Duplex Consensus Pairs'
#       description: 'read pairs called as duplex consensus'
#       format: '{{:,.0f}}'
#       placement: 920
#   SingleConsensusPairs:
#       title: 'Single Consensus Pairs'
#       description: 'read pairs called as single consensus'
#       format: '{{:,.0f}}'
#   DupReadsPairs:
#       title: 'Duplicated Read Pairs'
#       description: 'read pairs marked as dups during consensus calling'
#       format: '{{:,.0f}}'
#       placement: 1010
#   DupRate:
#       title: 'Duplication Rate'
#       description: 'duplication rate calculated by dividing Dup Reads Pairs by Read Pairs' 
#       format: '{{:,.3f}}'
#       placement: 1020
Sample\\tReadPairs\\tDuplexConsensusPairs\\tSingleConsensusPairs\\tDupReadsPairs\\tDupRate" > {output}
for f in {input}; do
    sample=$(basename $f .dedup.hybrid.stats)
    read_pairs=$(grep "correctly-paired read pairs for MBC Consensus calling:" $f | awk '{{ print $9 }}')
    duplex_consensus=$(grep "read pairs called as duplex consensus:" $f | awk '{{ print $8 }}')
    single_consensus=$(grep "read pairs called as single consensus:" $f | awk '{{ print $8 }}')
    dup_read_pairs=$(grep "read pairs marked as dups during consensus calling:" $f | awk '{{ print $10 }}')
    dup_rate=$(awk -v r=$read_pairs -v dups=$dup_read_pairs BEGIN'{{ print dups / r }}')

    echo -e "$sample\\t$read_pairs\\t$duplex_consensus\\t$single_consensus\\t$dup_read_pairs\\t$dup_rate" >> {output}
done""".format(
    input=" ".join(["" + input for input in metrics_inputs]),
    output=creak_metrics
            )
        )

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")),
                        parse_job
                    ],
                    name = "parse_creak_metrics"
                    )
                )

        self.multiqc_inputs.append(creak_metrics)
        return jobs

    def creak_dedup_bam(self):
        """
        Deduplicate bam files with AGeNT creak in hybrid (saliva only) or duplex (both saliva and brush) mode.
        Testing for SureSelect samples in vardict protocol.
        Saliva and brush samples are identified via the design file. 
        """

        jobs = []

        metrics_inputs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            sorted_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")

            if sample in self.contrasts.salivas:
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                index_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bai")
                output_hybrid = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.bam")
                covered_bed = config.param('agent_creak', 'covered_bedv8', param_type='filepath')
                source = "_saliva"

                jobs.append(
                        concat_jobs(
                            [
                                agent.creak(
                                    input_bam,
                                    output_duplex,
                                    "DUPLEX",
                                    covered_bed
                                    ),
                                agent.creak(
                                    input_bam,
                                    output_hybrid,
                                    "HYBRID",
                                    covered_bed
                                    ),
                                bash.ln(
                                    os.path.relpath(output_duplex, os.path.dirname(sorted_bam)),
                                    sorted_bam,
                                    output_duplex
                                    ),
                                bash.ln(
                                    os.path.relpath(index_duplex, os.path.dirname(sorted_bam)),
                                    sorted_bam + ".bai",
                                    index_duplex
                                    ),
                                bash.ln(
                                    os.path.relpath(re.sub(".bam", ".stats", output_duplex), os.path.dirname(sorted_bam)),
                                    os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.stats"),
                                    re.sub(".bam", ".stats", output_duplex),
                                    )
                            ],
                            name='agent_creak.' + sample.name,
                            samples=[sample]
                            )
                        )
            elif sample in self.contrasts.brushes:
                output_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                index_duplex = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bai")
                covered_bed = config.param('agent_creak', 'covered_bedv7', param_type='filepath')
                source = "_brush"
                
                jobs.append(
                        concat_jobs(
                            [
                                agent.creak(
                                    input_bam,
                                    output_duplex,
                                    "DUPLEX",
                                    covered_bed
                                    ),
                                bash.ln(
                                    os.path.relpath(re.sub(".bam", ".stats", output_duplex), os.path.dirname(index_duplex)),
                                    os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.stats"),
                                    re.sub(".bam", ".stats", output_duplex),
                                    )
                                #,     # skip this symlink creation here and create symlink instead after subsample to keep names between samples consistent.
                                #bash.ln(
                                #    os.path.relpath(output_duplex, os.path.dirname(sorted_bam)),
                                #    sorted_bam,
                                #    output_duplex
                                #    ),
                                #bash.ln(
                                #    os.path.relpath(index_duplex, os.path.dirname(sorted_bam)),
                                #    sorted_bam + ".bai",
                                #    index_duplex
                                #    )
                                ],
                            name='agent_creak.' + sample.name,
                            samples=[sample]
                            )
                        )
                        
            else:
                _raise(SanitycheckError("Error: sample \"" + sample.name +
                "\" is not included in the design file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))
            
            metrics_inputs.append(os.path.join(alignment_directory, sample.name + source + ".dedup.duplex.stats"))

        creak_metrics = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs", "creak_dedup.duplex.stats_mqc.txt")
        
        parse_job=Job(
                metrics_inputs,
                [creak_metrics],
                command="""\
echo -e "# plot_type: 'table'
# section_name: 'CReaK'
# description: 'stats from deduplication with AGeNT CReaK'
# headers:
#   ReadPairs:
#       title: 'Read Pairs'
#       description: 'correctly-paired read pairs for MBC Consensus calling'
#       format: '{{:,.0f}}'
#       placement: 900
#   DuplexConsensusPairs:
#       title: 'Duplex Consensus Pairs'
#       description: 'read pairs called as duplex consensus'
#       format: '{{:,.0f}}'
#       placement: 920
#   SingleConsensusPairs:
#       title: 'Single Consensus Pairs'
#       description: 'read pairs called as single consensus'
#       format: '{{:,.0f}}'
#   DupReadsPairs:
#       title: 'Duplicated Read Pairs'
#       description: 'read pairs marked as dups during consensus calling'
#       format: '{{:,.0f}}'
#       placement: 1010
#   DupRate:
#       title: 'Duplication Rate'
#       description: 'duplication rate calculated by dividing Dup Reads Pairs by Read Pairs' 
#       format: '{{:,.3f}}'
#       placement: 1020
Sample\\tReadPairs\\tDuplexConsensusPairs\\tSingleConsensusPairs\\tDupReadsPairs\\tDupRate" > {output}
for f in {input}; do
    sample=$(basename $f .dedup.duplex.stats)
    read_pairs=$(grep "correctly-paired read pairs for MBC Consensus calling:" $f | awk '{{ print $9 }}')
    duplex_consensus=$(grep "read pairs called as duplex consensus:" $f | awk '{{ print $8 }}')
    single_consensus=$(grep "read pairs called as single consensus:" $f | awk '{{ print $8 }}')
    dup_read_pairs=$(grep "read pairs marked as dups during consensus calling:" $f | awk '{{ print $10 }}')
    dup_rate=$(awk -v r=$read_pairs -v dups=$dup_read_pairs BEGIN'{{ print dups / r }}')

    echo -e "$sample\\t$read_pairs\\t$duplex_consensus\\t$single_consensus\\t$dup_read_pairs\\t$dup_rate" >> {output}
done""".format(
    input=" ".join(["" + input for input in metrics_inputs]),
    output=creak_metrics
            )
        )

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")),
                        parse_job
                    ],
                    name = "parse_creak_metrics"
                    )
                )

        self.multiqc_inputs.append(creak_metrics)

        return jobs

    def samtools_subsample(self):
        """
        Vardict protocol only.
        Subsample the deduplicated bams to 1.5M reads if the number of overlapping reads exceeds this number (brush samples only).
        Number of overlapping reads is assessed with samtools view -c, from which a fraction is calculated for subsampling.
        If fraction >= 1, the bam is not subsampled.
        """
        jobs = []

        for sample in self.samples:
            if sample in self.contrasts.brushes:
                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.bam")
                sorted_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")

                covered_bed = config.param('agent_creak', 'covered_bedv7', param_type='filepath')
                seed = config.param('samtools_subsample', 'seed')

                jobs.append(
                        concat_jobs(
                            [
                            Job(
                                command="""\
READS=`samtools view -c -L {covered_bed} {input_bam}` &&
if [[ $READS -ge 1500000 ]]; then 
    FRACTION=`awk -v nr=$READS \'BEGIN{{ print 1500000 / nr }}\'`
else
    FRACTION=1; fi""".format(
                                sorted_bam=sorted_bam,
                                covered_bed=covered_bed,
                                input_bam=input_bam
                                )
                            ),
                            samtools.view(
                                input_bam,
                                sorted_bam,
                                "-b --subsample-seed " + seed + " --subsample $FRACTION"
                                ),
                            samtools.index(
                                sorted_bam
                                        ) 
                            ],
                        name = "samtools_subsample." + sample.name,
                        samples = [sample]
                        )
                    )

        return jobs

    def samtools_sort(self):
        """
        Sort deduplicated bams by coordinate and index with samtools. Note: Not needed when using -S flag with locatit or using creak for deduplication.
        """

        jobs = []

        dovee_protocol = self.args.type

        if dovee_protocol == "vardict":

            for sample in self.samples:
                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                bam_file_prefix = os.path.join(alignment_directory, sample.name + ".dedup.")

                if sample in self.contrasts.brushes:
                    input_dedup = bam_file_prefix + "duplex.bam" 
                    output_dedup = bam_file_prefix + "duplex.sorted"
    
                    jobs.append(
                            concat_jobs(
                                [
                                    samtools.sort(
                                        input_dedup,
                                        output_dedup
                                        ),
                                    samtools.index(
                                        output_dedup + ".bam"
                                        ) 
                                ], name = "samtools_sort." + sample.name,
                                samples = [sample]
                                )
                            )
    
                elif sample in self.contrasts.salivas:
                    input_duplex = bam_file_prefix + "duplex.bam"
                    input_hybrid = bam_file_prefix + "hybrid.bam"
                    output_duplex = bam_file_prefix + "duplex.sorted"
                    output_hybrid = bam_file_prefix + "hybrid.sorted"
    
                    jobs.append(
                            concat_jobs(
                                [
                                    samtools.sort(
                                        input_duplex,
                                        output_duplex
                                        ),
                                    samtools.index(
                                        output_duplex + ".bam"
                                        ),
                                    samtools.sort(
                                        input_hybrid,
                                        output_hybrid
                                        ),
                                    samtools.index(
                                        output_hybrid + ".bam"
                                        )
                                ], name = "samtools_sort." + sample.name,
                                samples = [sample]
                                )
                            )

                else:
                    _raise(SanitycheckError("Error: sample \"" + sample.name +
                    "\" is not included in the design file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))

        elif dovee_protocol == "copy-number":
            for sample in self.samples:
                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                bam_file_prefix = os.path.join(alignment_directory, sample.name + ".dedup.")
                input_dedup = bam_file_prefix + "hybrid.bam" 
                output_dedup = bam_file_prefix + "hybrid.sorted"

                jobs.append(
                        concat_jobs(
                            [
                                samtools.sort(
                                    input_dedup,
                                    output_dedup
                                    ),
                                samtools.index(
                                    output_dedup + ".bam"
                                    ) 
                            ], name = "samtools_sort." + sample.name,
                            samples = [sample]
                            )
                        )
        return jobs

    def mosdepth(self):
        """
        Calculate depth stats for captured regions with mosdepth.
        """

        jobs = []

        dovee_protocol = self.args.type

        mosdepth_directory = os.path.join(self.output_dirs['metrics_directory'], "mosdepth")
        link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
        
        if dovee_protocol == "vardict":
            for sample in self.samples:
                alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
                output_prefix = os.path.join(mosdepth_directory, sample.name)
                input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")
                if sample in self.contrasts.salivas:
                    region=config.param('mosdepth', 'region_bed', param_type='filepath')
                    source="_saliva"
                elif sample in self.contrasts.brushes:
                    region=config.param('mosdepth', 'region_bed', param_type='filepath')
                    source="_brush"

                jobs.append(
                        concat_jobs(
                            [
                                bash.mkdir(mosdepth_directory),
                                bash.mkdir(link_directory),
                                mosdepth.mosdepth(
                                    input_bam,
                                    output_prefix,
                                    True,
                                    region
                                    ),
                                bash.ln(
                                    os.path.relpath(output_prefix + ".mosdepth.region.dist.txt", link_directory),
                                    os.path.join(link_directory, sample.name + source + ".mosdepth.region.dist.txt"),
                                    output_prefix + ".mosdepth.region.dist.txt"
                                    )
                            ],
                        name = "mosdepth." + sample.name,
                        samples = [sample]
                        )
                    )
                self.multiqc_inputs.append(output_prefix + ".mosdepth.region.dist.txt")
        
        elif dovee_protocol == "copy-number":
            for sample_pair in self.dovee_pairs.values():
                
                input_brush = os.path.join(self.output_dirs['alignment_directory'], sample_pair.brush.name, sample_pair.brush.name + ".dedup.hybrid.sorted.bam")
                input_saliva = os.path.join(self.output_dirs['alignment_directory'], sample_pair.saliva.name, sample_pair.saliva.name + ".dedup.hybrid.sorted.bam")
                output_prefix_brush = os.path.join(mosdepth_directory, sample_pair.brush.name)
                output_prefix_saliva = os.path.join(mosdepth_directory, sample_pair.saliva.name)
                region="1000" # window size

                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(mosdepth_directory),
                            bash.mkdir(link_directory),
                            mosdepth.mosdepth(
                                input_brush,
                                output_prefix_brush,
                                True,
                                region
                                ),
                            bash.ln(
                                os.path.relpath(output_prefix_brush + ".mosdepth.region.dist.txt", link_directory),
                                os.path.join(link_directory, sample_pair.brush.name + "_brush" + ".mosdepth.region.dist.txt"),
                                output_prefix_brush + ".mosdepth.region.dist.txt"
                                )
                            ],
                        name = "mosdepth." + sample_pair.brush.name,
                        samples = [sample_pair.brush]
                        )
                    )
                self.multiqc_inputs.append(output_prefix_brush + ".mosdepth.region.dist.txt")
                
                jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(mosdepth_directory),
                            bash.mkdir(link_directory),
                            mosdepth.mosdepth(
                                input_saliva,
                                output_prefix_saliva,
                                True,
                                region
                                ),
                            bash.ln(
                                os.path.relpath(output_prefix_saliva + ".mosdepth.region.dist.txt", link_directory),
                                os.path.join(link_directory, sample_pair.saliva.name + "_saliva" + ".mosdepth.region.dist.txt"),
                                output_prefix_saliva + ".mosdepth.region.dist.txt"
                                )
                            ],
                        name = "mosdepth." + sample_pair.saliva.name,
                        samples = [sample_pair.saliva]
                        )
                    )
                self.multiqc_inputs.append(output_prefix_saliva + ".mosdepth.region.dist.txt")
        return jobs

    def picard_metrics(self):
        """
        Collect on and off target metrics for SureSelect samples with Picard.
        Collect gc bias metrics for SureSelect samples with Picard.
        """
        jobs = []

        picard_directory = os.path.join(self.output_dirs['metrics_directory'], "picard")
        link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
        
        targetv7 = config.param('vardict_single', 'target_filev7', param_type='filepath')  #target bed file 
        targetv8 = config.param('vardict_single', 'target_filev8', param_type='filepath')  #target bed file

        outputv7 = os.path.join(picard_directory, "targetv7.interval_list")
        outputv8 = os.path.join(picard_directory, "targetv8.interval_list")

        jobs.append(
                concat_jobs(
                    [
                        bash.mkdir(picard_directory),
                        picard2.bed2interval_list(
                            None,
                            targetv7,
                            outputv7
                            ),
                        picard2.bed2interval_list(
                            None,
                            targetv8,
                            outputv8
                            )
                        ],
                    name = "picard.bed2interval_list"
                    )
                )

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")

            output = os.path.join(picard_directory, sample.name + ".picard_HS_metrics.txt")

            if sample in self.contrasts.brushes:
                intervals = outputv7
                source = "_brush"

            elif sample in self.contrasts.salivas:
                intervals = outputv8
                source = "_saliva"

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(link_directory),
                            picard2.calculate_hs_metrics(
                                input_bam,
                                output,
                                intervals
                                ),
                            bash.ln(
                                os.path.relpath(output, link_directory),
                                os.path.join(link_directory, sample.name + source + ".picard_HS_metrics.txt"),
                                input=output
                                )
                        ], name = "picard_calculate_hs_metrics." + sample.name,
                        samples = [sample]
                        )
                    )
            self.multiqc_inputs.append(picard_directory)
            
            output_prefix = os.path.join(picard_directory, sample.name)

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(link_directory),
                            picard2.collect_gcbias_metrics(
                                input_bam,
                                output_prefix
                                ),
                            bash.ln(
                                os.path.relpath(output_prefix + ".gcbias_metrics.txt", link_directory),
                                os.path.join(link_directory, sample.name + source + ".gcbias_metrics.txt"),
                                input=output_prefix + ".gcbias_metrics.txt"
                                ),
                            bash.ln(
                                os.path.relpath(output_prefix + ".gcbias_summary_metrics.txt", link_directory),
                                os.path.join(link_directory, sample.name + source + ".gcbias_summary_metrics.txt"),
                                input=output_prefix + ".gcbias_summary_metrics.txt"
                                )
                            ],
                        name = "picard_collect_gcbias_metrics." + sample.name,
                        samples = [sample]
                        )
                    )
                    
        return jobs

    def vardict_single(self):
        """
        Variant calling with vardict in single mode.
        """

        jobs = []

        if self.contrasts:
            design_file = os.path.relpath(self.args.design.name, self.output_dir)

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            variants_directory = os.path.join(self.output_dirs['variants_directory'], sample.name)
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.duplex.sorted.bam")
            output = os.path.join(variants_directory, sample.name + ".vcf")

            if sample in self.contrasts.brushes:
                freq=0.001
                region=config.param('vardict_single', 'target_filev7', param_type='filepath')  #target bed file
                nosv=True

            elif sample in self.contrasts.salivas:
                freq=0.1
                region=config.param('vardict_single', 'target_filev8', param_type='filepath')  #target bed file
                nosv=False

            else:
                _raise(SanitycheckError("Error: sample \"" + sample.name +
                "\" is not included in the design file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))

            vardict_job = vardict.single_java(
                    input_bam,
                    sample.name,
                    None,
                    nosv=nosv,
                    freq=freq,
                    region=region
                    )

            teststrandbias_job = vardict.teststrandbias(
                    None,
                    None
                    )

            var2vcf_job = vardict.var2vcf_valid(
                    output,
                    sample.name,
                    freq,
                    None
                    )

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(variants_directory),
                            pipe_jobs(
                                [
                                    vardict_job,
                                    teststrandbias_job,
                                    var2vcf_job
                                    ]
                                )
                            ], 
                        name = "vardict_single." + sample.name,
                        samples = [sample]
                        )
                    )

        return jobs

    def hmm_readCounter(self):
        """
        Counting number of reads in non-overlapping windows of fixed width directly from BAM files with HMM Copy Utils readCounter:
        https://github.com/shahcompbio/hmmcopy_utils
        Run on both saliva and brush, as both wigs needed as input for IchorCNA.
        """

        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join(self.output_dirs['alignment_directory'], sample.name)
            wig_directory = os.path.join(self.output_dirs['wig_directory'], sample.name) 
            input_bam = os.path.join(alignment_directory, sample.name + ".dedup.hybrid.sorted.bam") 
            output = os.path.join(wig_directory, sample.name + ".out.wig")

            jobs.append(
                    concat_jobs(
                    [
                        bash.mkdir(wig_directory),
                        hmm.readCounter(
                            input_bam,
                            output
                            )
                    ], 
                    name = "hmm_readCounter." + sample.name,
                    samples = [sample]
                )
            )

        return jobs

    def run_ichorCNA(self):
        """
        https://github.com/broadinstitute/ichorCNA
        'Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.'
        Run once per patient using paired saliva and brush sample wigs, requires knowing which saliva and brush came from same patient.
        """

        jobs = []
        
        for sample_pair in self.dovee_pairs.values(): # pairs in pair file must be complete
            wig_directory = self.output_dirs['wig_directory']
            input_brush = os.path.join(wig_directory, sample_pair.brush.name, sample_pair.brush.name + ".out.wig")
            input_saliva = os.path.join(wig_directory, sample_pair.saliva.name, sample_pair.saliva.name + ".out.wig")
            output_dir = os.path.join(self.output_dirs['cna_directory'], sample_pair.name)
    
            jobs.append(
                    concat_jobs(
                        [
                        bash.mkdir(output_dir),
                        ichorCNA.run_ichorCNA(
                            input_brush,
                            input_saliva,
                            sample_pair.name,
                            output_dir
                            )
                        ], 
                        name = "run_ichorCNA." + sample_pair.name,
                        samples = [sample_pair.saliva, sample_pair.brush]
                    )
                )
    
        return jobs

    def bcftools_stats(self):
        """
        Collect stats number and types of variants in vcf files produced by vardict_single.
        """
        
        jobs = []

        for sample in self.samples:
            variants_directory = os.path.join(self.output_dirs['variants_directory'], sample.name)
            input = os.path.join(variants_directory, sample.name + ".vcf")
            output_directory = os.path.join(self.output_dirs['metrics_directory'], "bcftools")
            link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
            output = os.path.join(output_directory, sample.name + ".bcftools.stats")
            
            if sample in self.contrasts.salivas:
                source="_saliva"
            elif sample in self.contrasts.brushes:
                source="_brush"
            else:
                _raise(SanitycheckError("Error: sample \"" + sample.name +
                "\" is not included in the design file! Cannot determine whether \"" + sample.name + "\" is a saliva or brush sample."))
            
            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(output_directory),
                            bash.mkdir(link_directory),
                            bcftools.stats(
                                input,
                                output                                
                                ),
                            bash.ln(
                                os.path.relpath(output, link_directory),
                                os.path.join(link_directory, sample.name + source + ".bcftools.stats"),
                                input=output
                                )
                        ], name = "bcftools_stats." + sample.name,
                        samples = [sample]
                    )
                )
        
            self.multiqc_inputs.append(output)
        
        return jobs

    def conpair_concordance(self):
        """
        Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification
        (= samples coming from the same individual).
        Run once per patient to ensure that saliva and brush samples are correctly assigned.
        Requires pair file.
        """

        jobs = []

        conpair_directory = os.path.join(self.output_dirs['metrics_directory'], 'concordance')
        link_directory = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")

        for sample_pair in self.dovee_pairs.values():
            input_brush = os.path.join(self.output_dirs['alignment_directory'], sample_pair.brush.name, sample_pair.brush.name + ".dedup.duplex.sorted.bam")
            input_saliva = os.path.join(self.output_dirs['alignment_directory'], sample_pair.saliva.name, sample_pair.saliva.name + ".dedup.duplex.sorted.bam") 
            output_dir = os.path.join(conpair_directory, sample_pair.name)
            pileup_brush = os.path.join(output_dir, sample_pair.brush.name + ".pileup")
            pileup_saliva = os.path.join(output_dir, sample_pair.saliva.name + ".pileup")
            concordance_out = os.path.join(output_dir, sample_pair.name + ".concordance.out")

            jobs.append(
                    concat_jobs(
                        [
                            bash.mkdir(
                                output_dir
                                ),
                            bash.mkdir(link_directory),
                            conpair.pileup(
                                input_brush,
                                pileup_brush
                                ),
                            conpair.pileup(
                                input_saliva,
                                pileup_saliva
                                ),
                            conpair.concordance(
                                pileup_saliva,
                                pileup_brush,
                                concordance_out
                                ),
                            bash.ln(
                                os.path.relpath(concordance_out, link_directory),
                                os.path.join(link_directory, sample_pair.name + ".concordance.out"),
                                input=concordance_out
                                )
                        ],
                        name = "conpair_concordance_contamination." + sample_pair.name,
                        samples = [sample_pair.brush, sample_pair.saliva]
                        )
                    )
            self.multiqc_inputs.append(concordance_out)
        
        return jobs

    def multiqc(self):
        """
        Aggregate results from bioinformatics analyses across many samples into a single report
        MultiQC searches a given directory for analysis logs and compiles a HTML report. It's a general use tool,
        perfect for summarising the output from numerous bioinformatics tools [MultiQC](https://multiqc.info/).
        """

        jobs = []

        input_links = os.path.join(self.output_dirs['metrics_directory'], "multiqc_inputs")
        output = os.path.join(self.output_dirs['metrics_directory'], "multiqc")

        job = multiqc.run(
            [input_links],
            output
            )
        job.name = "multiqc"
        job.input_dependency = self.multiqc_inputs
        jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
                [
                    self.trimmer,
                    self.fastp,
                    self.bwa_mem_samtools_sort,
                    self.samtools_merge,
                    self.locatit_dedup_bam,
                   # self.creak_dedup_bam,  # protocol can be switched to use agent creak deduplcation, instead of locatit
                    self.samtools_subsample,
                    self.mosdepth,
                    self.picard_metrics,
                    self.vardict_single,
                    self.bcftools_stats,
                    self.conpair_concordance,
                    self.multiqc
                ],
                [
                    self.trimmer, 
                    self.fastp, 
                    self.bwa_mem_samtools_sort,
                    self.samtools_merge,
                    self.locatit_hybrid_dedup,
                    #self.creak_hybrid_dedup,
                    self.mosdepth,
                    self.hmm_readCounter,
                    self.run_ichorCNA,
                    self.multiqc
                 ]
            ]


if __name__ == "__main__":
    argv = sys.argv
    if '--wrap' in argv:
        utils.utils.container_wrapper_argparse(argv)
    else:
        DOvEE_gene(protocol=['vardict', 'copy-number'])
